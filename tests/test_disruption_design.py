import json
import unittest
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import webapp
from crispr4p.annotations import GeneAnnotation, GenomeAnnotations
from crispr4p.disruption import (
    StopCassette,
    build_donor,
    has_junction_pam,
    load_cassettes,
    recut_sites,
    target_gene,
    target_strand,
)
from crispr4p.models import DesignResult
from crispr4p.primers import (
    InsertionChecks,
    InsertionPrimerPair,
    JunctionPrimerPair,
    PrimerNotFoundError,
    insertion_checks,
    insertion_primers,
    overlap_oligos,
)
from crispr4p.resources import read_fasta
from crispr4p.service import Crispr4pService
from crispr4p.spedit import reverse_complement
from crispr4p.web_views import (
    annotation_rows,
    cassette_data,
    donor_data,
    render_design,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA = PROJECT_ROOT / "data"
CANDIDATES = DATA / "stop_cassettes.json"
FASTA = DATA / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa"
FIRST = "TGATGACGTGATGACCTAGTAGG"
OLD_FIRST = "TAGTAGCCTAGTGACCTAGTAGG"
ADE6_GUIDE = "ACATTGGCTTACGACGGTCG"
ADE6_REVERSE_GUIDE = "GTGGCGACAGGGACACCTCG"
BUB1_GUIDE = "TATCAGATTGCTCGGCCACA"


def gene(gene_id, name, strand, gene_type="protein_coding_gene"):
    return GeneAnnotation(
        gene_id=gene_id,
        name=name,
        gene_type=gene_type,
        start=100,
        end=300,
        strand=strand,
        viability="viable",
    )


def unsafe_cassette():
    return SimpleNamespace(
        orient=lambda strand: (
            reverse_complement(OLD_FIRST) if strand == "-" else OLD_FIRST
        )
    )


class DisruptionDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cassettes = load_cassettes(CANDIDATES)
        records = read_fasta(FASTA)
        cls.reference = records["III"].sequence
        cls.references = tuple(record.sequence for record in records.values())

    def test_packaged_candidates(self):
        self.assertEqual(tuple(range(1, 11)), tuple(c.id for c in self.cassettes))
        self.assertFalse(any(has_junction_pam(c.sequence) for c in self.cassettes))
        self.assertEqual(FIRST, self.cassettes[0].sequence)
        self.assertEqual(
            "GTGATGAGTGATGACTAGTGAGG",
            self.cassettes[-1].sequence,
        )

        data = json.loads(CANDIDATES.read_text(encoding="utf-8"))
        self.assertEqual(1, data["version"])
        self.assertEqual("ASM294v2.26", data["assembly"])
        for cassette in self.cassettes:
            for reference in self.references:
                self.assertNotIn(cassette.sequence, reference)
                self.assertNotIn(
                    reverse_complement(cassette.sequence),
                    reference,
                )

    def test_candidate_details(self):
        cassette = self.cassettes[0]

        self.assertEqual("TGATGACGTGATGACCTAGT", cassette.guide)
        self.assertEqual("AGG", cassette.pam)
        self.assertEqual(45.0, cassette.gc_percent)
        self.assertEqual(3, len(cassette.frames))
        self.assertEqual(
            "CCTACTAGGTCATCACGTCATCA",
            cassette.orient("-"),
        )
        self.assertEqual(FIRST, cassette.orient("+"))

    def test_invalid_candidate(self):
        with self.assertRaisesRegex(ValueError, "23 nt"):
            StopCassette(1, "TAG")
        with self.assertRaisesRegex(ValueError, "tandem stops"):
            StopCassette(1, "ACGTACGTACGTACGTACGTAGG")
        with self.assertRaisesRegex(ValueError, "beside the junction"):
            StopCassette(1, OLD_FIRST)

    def test_junction_recut(self):
        sites = recut_sites(
            self.reference,
            (1316791, 1316792),
            ADE6_GUIDE,
            unsafe_cassette(),
            "+",
        )

        self.assertTrue(has_junction_pam(OLD_FIRST))
        self.assertEqual(1, len(sites))
        self.assertEqual("ACATTGGCTTACGACGGTAG", sites[0].target)
        self.assertEqual("TAG", sites[0].pam)
        self.assertEqual(1, sites[0].mismatches)
        reverse_sites = recut_sites(
            self.reference,
            (1317806, 1317807),
            ADE6_REVERSE_GUIDE,
            unsafe_cassette(),
            "-",
        )
        self.assertEqual("-", reverse_sites[0].target_strand)
        self.assertEqual(1, reverse_sites[0].mismatches)
        self.assertFalse(
            recut_sites(
                self.reference,
                (1316791, 1316792),
                ADE6_GUIDE,
                self.cassettes[0],
            )
        )

    def test_disruption_donor(self):
        donor = build_donor(
            self.reference,
            (1316791, 1316792),
            self.cassettes[0],
            "+",
            80,
        )

        self.assertEqual(self.reference[1316711:1316791], donor.left_arm)
        self.assertEqual(FIRST, donor.insert)
        self.assertEqual(self.reference[1316791:1316871], donor.right_arm)
        self.assertEqual(183, donor.total_length)
        self.assertEqual(reverse_complement(donor.sequence), donor.reverse)

        oligos = donor.oligos
        self.assertEqual(donor.left_arm + donor.insert, oligos.forward)
        self.assertEqual(
            reverse_complement(donor.insert + donor.right_arm),
            oligos.reverse,
        )
        self.assertEqual(donor.insert, oligos.overlap)
        self.assertEqual(
            reverse_complement(oligos.overlap),
            oligos.reverse[-len(oligos.overlap):],
        )
        self.assertEqual((103, 103), (len(oligos.forward), len(oligos.reverse)))
        self.assertEqual(donor.sequence, oligos.product)
        self.assertEqual(183, oligos.product_length)

        minus = build_donor(
            self.reference,
            (1316791, 1316792),
            self.cassettes[0],
            "-",
            80,
        )
        self.assertEqual(reverse_complement(FIRST), minus.insert)
        self.assertEqual(minus.insert, minus.oligos.overlap)
        self.assertEqual(minus.sequence, minus.oligos.product)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            build_donor(
                self.reference,
                (1316791, 1316792),
                self.cassettes[0],
                "+",
                0,
            )
        with self.assertRaisesRegex(ValueError, "complete homology arms"):
            build_donor(
                self.reference,
                (40, 41),
                self.cassettes[0],
                "+",
                80,
            )

    def test_construction_oligos_cover_every_cassette_orientation(self):
        for cassette in self.cassettes:
            for strand in ("+", "-"):
                donor = build_donor(
                    self.reference,
                    (1316791, 1316792),
                    cassette,
                    strand,
                    80,
                )
                oligos = donor.oligos
                self.assertEqual(
                    (103, 103, 23),
                    (
                        len(oligos.forward),
                        len(oligos.reverse),
                        len(oligos.overlap),
                    ),
                )
                self.assertEqual(donor.insert, oligos.overlap)
                self.assertEqual(donor.sequence, oligos.product)

    def test_construction_oligo_validation(self):
        with self.assertRaisesRegex(ValueError, "invalid nucleotide"):
            overlap_oligos("AACN", 1, 3)
        with self.assertRaisesRegex(ValueError, "within the donor"):
            overlap_oligos("AACG", 2, 2)

    def test_insertion_primer_inputs(self):
        answer = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 250,
            "PRIMER_LEFT_0": [130, 20],
            "PRIMER_RIGHT_0": [379, 20],
        }
        designer = Mock(return_value=answer)
        pair = insertion_primers(
            "A" * 300 + "C" * 300,
            (300, 301),
            primer_designer=designer,
        )
        sequence_args, global_args = designer.call_args.args

        self.assertEqual(
            "A" * 300 + "C" * 300,
            sequence_args["SEQUENCE_TEMPLATE"],
        )
        self.assertEqual(
            [[220, 160]],
            sequence_args["SEQUENCE_EXCLUDED_REGION"],
        )
        self.assertEqual(
            [[0, 220, 380, 220]],
            sequence_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"],
        )
        self.assertEqual(
            [[200, 300]],
            global_args["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(0, global_args["PRIMER_PICK_INTERNAL_OLIGO"])
        self.assertEqual(1, global_args["PRIMER_NUM_RETURN"])
        self.assertEqual("FORWARD", pair.forward)
        self.assertEqual("REVERSE", pair.reverse)
        self.assertEqual(250, pair.wt_product_size)
        self.assertEqual(273, pair.disrupted_product_size)
        self.assertEqual(130, pair.forward_start)
        self.assertEqual(379, pair.reverse_end)

        no_pair = Mock(return_value={"PRIMER_PAIR_NUM_RETURNED": 0})
        with self.assertRaises(PrimerNotFoundError):
            insertion_primers(
                "A" * 300 + "C" * 300,
                (300, 301),
                primer_designer=no_pair,
            )

    def test_junction_primer_inputs(self):
        spanning = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 298,
            "PRIMER_LEFT_0": [130, 20],
            "PRIMER_RIGHT_0": [427, 20],
        }
        left = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": reverse_complement(FIRST),
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 59.1,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 193,
        }
        right = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": FIRST,
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.1,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 151,
        }
        designer = Mock(side_effect=(spanning, left, right))

        checks = insertion_checks(
            "A" * 300 + "C" * 300,
            (300, 301),
            FIRST,
            primer_designer=designer,
        )

        self.assertEqual(3, designer.call_count)
        left_args, left_settings = designer.call_args_list[1].args
        right_args, right_settings = designer.call_args_list[2].args
        edited = "A" * 300 + FIRST + "C" * 300
        self.assertEqual(edited, left_args["SEQUENCE_TEMPLATE"])
        self.assertEqual("FORWARD", left_args["SEQUENCE_PRIMER"])
        self.assertEqual(
            reverse_complement(FIRST),
            left_args["SEQUENCE_PRIMER_REVCOMP"],
        )
        self.assertEqual(FIRST, right_args["SEQUENCE_PRIMER"])
        self.assertEqual("REVERSE", right_args["SEQUENCE_PRIMER_REVCOMP"])
        self.assertEqual("check_primers", left_settings["PRIMER_TASK"])
        self.assertEqual(57.0, left_settings["PRIMER_MIN_TM"])
        self.assertEqual(63.0, left_settings["PRIMER_MAX_TM"])
        self.assertEqual(
            [[193, 193]],
            left_settings["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(
            [[151, 151]],
            right_settings["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(193, checks.left.product_size)
        self.assertEqual(151, checks.right.product_size)

        failed = Mock(side_effect=(spanning, {"PRIMER_PAIR_NUM_RETURNED": 0}))
        with self.assertRaises(PrimerNotFoundError):
            insertion_checks(
                "A" * 300 + "C" * 300,
                (300, 301),
                FIRST,
                primer_designer=failed,
            )

        problem = dict(left)
        problem["PRIMER_RIGHT_0_PROBLEMS"] = "high hairpin stability"
        rejected = Mock(side_effect=(spanning, problem))
        with self.assertRaisesRegex(PrimerNotFoundError, "high hairpin"):
            insertion_checks(
                "A" * 300 + "C" * 300,
                (300, 301),
                FIRST,
                primer_designer=rejected,
            )

    def test_real_insertion_primers(self):
        service = Crispr4pService.from_project_data()
        guide = service.design_gene("ade6").guides[0]
        pair = service.insertion_primers(
            guide.chromosome,
            guide.cut_coordinates,
        )
        cut = guide.cut_coordinates[0]

        self.assertEqual("GCAACTCTGCGATGCATTCA", pair.forward)
        self.assertEqual("TGCGTACTACCATCACTGCA", pair.reverse)
        self.assertAlmostEqual(59.551009575121896, pair.forward_tm)
        self.assertAlmostEqual(59.10793666168439, pair.reverse_tm)
        self.assertEqual(298, pair.wt_product_size)
        self.assertEqual(321, pair.disrupted_product_size)
        self.assertNotEqual(
            -1,
            self.reference.find(pair.forward, cut - 300, cut - 80),
        )
        self.assertNotEqual(
            -1,
            self.reference.find(
                reverse_complement(pair.reverse),
                cut + 80,
                cut + 300,
            ),
        )

    def test_real_junction_primers(self):
        service = Crispr4pService.from_project_data()
        guide = service.design_gene("ade6").guides[0]

        for cassette in self.cassettes:
            checks = service.insertion_checks(
                guide.chromosome,
                guide.cut_coordinates,
                cassette.id,
                "+",
            )
            self.assertEqual(checks.spanning.forward, checks.left.forward)
            self.assertEqual(
                reverse_complement(cassette.sequence),
                checks.left.reverse,
            )
            self.assertEqual(cassette.sequence, checks.right.forward)
            self.assertEqual(checks.spanning.reverse, checks.right.reverse)
            self.assertEqual(193, checks.left.product_size)
            self.assertEqual(151, checks.right.product_size)

        reverse_guide = service.design_gene("bub1").guides[0]
        for cassette in self.cassettes:
            reverse_checks = service.insertion_checks(
                reverse_guide.chromosome,
                reverse_guide.cut_coordinates,
                cassette.id,
                "-",
            )
            self.assertEqual(cassette.sequence, reverse_checks.left.reverse)
            self.assertEqual(
                reverse_complement(cassette.sequence),
                reverse_checks.right.forward,
            )
            self.assertEqual(197, reverse_checks.left.product_size)
            self.assertEqual(123, reverse_checks.right.product_size)
        with self.assertRaisesRegex(ValueError, "coding strand"):
            service.insertion_checks(
                guide.chromosome,
                guide.cut_coordinates,
                1,
                None,
            )

    def test_insertion_primer_endpoint(self):
        service = Mock()
        service.insertion_primers.return_value = InsertionPrimerPair(
            forward="FORWARD",
            reverse="REVERSE",
            forward_tm=59.4,
            reverse_tm=60.2,
            wt_product_size=250,
            insert_length=23,
        )
        handler = webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)
        handler.send_response = Mock()
        handler.send_header = Mock()
        handler.end_headers = Mock()
        handler.wfile = BytesIO()

        with patch.object(webapp, "create_service", return_value=service):
            handler.serve_insertion_primers(
                {
                    "chromosome": ["III"],
                    "cut_left": ["1316791"],
                    "cut_right": ["1316792"],
                }
            )

        service.insertion_primers.assert_called_once_with(
            "III",
            (1316791, 1316792),
            arm_length=80,
            insert_length=23,
            window=300,
        )
        handler.send_response.assert_called_once_with(200)
        self.assertIn(
            ("Cache-Control", "no-store"),
            [call.args for call in handler.send_header.call_args_list],
        )
        self.assertEqual(
            {
                "forward": "FORWARD",
                "reverse": "REVERSE",
                "forward_tm": 59.4,
                "reverse_tm": 60.2,
                "wt_product_size": 250,
                "disrupted_product_size": 273,
            },
            json.loads(handler.wfile.getvalue()),
        )

    def test_junction_primer_endpoint(self):
        spanning = InsertionPrimerPair(
            forward="FORWARD",
            reverse="REVERSE",
            forward_tm=59.4,
            reverse_tm=60.2,
            wt_product_size=298,
            insert_length=23,
        )
        checks = InsertionChecks(
            spanning=spanning,
            left=JunctionPrimerPair(
                "FORWARD",
                "CASSETTE_REVERSE",
                59.4,
                59.1,
                193,
            ),
            right=JunctionPrimerPair(
                "CASSETTE_FORWARD",
                "REVERSE",
                59.1,
                60.2,
                151,
            ),
        )
        service = Mock()
        service.insertion_checks.return_value = checks
        handler = webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)
        handler.send_response = Mock()
        handler.send_header = Mock()
        handler.end_headers = Mock()
        handler.wfile = BytesIO()

        with patch.object(webapp, "create_service", return_value=service):
            handler.serve_insertion_primers(
                {
                    "chromosome": ["III"],
                    "cut_left": ["1316791"],
                    "cut_right": ["1316792"],
                    "cassette_id": ["1"],
                    "coding_strand": ["+"],
                }
            )

        service.insertion_checks.assert_called_once_with(
            "III",
            (1316791, 1316792),
            1,
            "+",
            arm_length=80,
            window=300,
        )
        self.assertEqual(
            {
                "forward": "FORWARD",
                "reverse": "REVERSE",
                "forward_tm": 59.4,
                "reverse_tm": 60.2,
                "wt_product_size": 298,
                "disrupted_product_size": 321,
                "left_junction": {
                    "forward": "FORWARD",
                    "reverse": "CASSETTE_REVERSE",
                    "forward_tm": 59.4,
                    "reverse_tm": 59.1,
                    "product_size": 193,
                },
                "right_junction": {
                    "forward": "CASSETTE_FORWARD",
                    "reverse": "REVERSE",
                    "forward_tm": 59.1,
                    "reverse_tm": 60.2,
                    "product_size": 151,
                },
            },
            json.loads(handler.wfile.getvalue()),
        )

    def test_minus_strand_donor(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("bub1")
        guide = result.guides[0]
        annotation = service.annotate_guide(guide)
        choices = service.cassette_choices(
            (guide,),
            (annotation,),
            result.name,
        )
        donors = service.disruption_donors(
            (guide,),
            (annotation,),
            choices,
            80,
            result.name,
        )
        donor = donors[0][0]
        cut = guide.cut_coordinates[0]

        self.assertEqual(BUB1_GUIDE, guide.seed)
        self.assertEqual(1, guide.strand)
        self.assertEqual((1314935, 1314937), guide.pam_coordinates)
        self.assertEqual((1314931, 1314932), guide.cut_coordinates)
        self.assertEqual("-", donor.coding_strand)
        self.assertEqual(self.reference[cut - 80:cut], donor.left_arm)
        self.assertEqual(self.reference[cut:cut + 80], donor.right_arm)
        self.assertEqual(reverse_complement(FIRST), donor.insert)
        self.assertEqual(
            reverse_complement(donor.right_arm)
            + FIRST
            + reverse_complement(donor.left_arm),
            donor.reverse,
        )

    def test_target_strand(self):
        plus = gene("plus", "target", "+")
        minus = gene("minus", "other", "-")
        non_coding = gene("rna", None, "-", "lncRNA_gene")

        self.assertEqual(
            "+",
            target_strand(SimpleNamespace(genes=(plus, minus)), "target"),
        )
        self.assertEqual(
            "-",
            target_strand(SimpleNamespace(genes=(minus,))),
        )
        self.assertIsNone(
            target_strand(SimpleNamespace(genes=(plus, minus)))
        )
        annotation = SimpleNamespace(genes=(plus, non_coding))
        self.assertEqual(non_coding, target_gene(annotation, "rna"))
        self.assertIsNone(target_strand(annotation, "rna"))

    def test_service_uses_packaged_catalog(self):
        service = Crispr4pService.from_project_data()

        self.assertEqual(self.cassettes, service.cassettes)
        self.assertIs(service.cassettes, service.cassettes)

    def test_service_filters_recut_sites(self):
        service = Crispr4pService.from_project_data(
            cassettes=(unsafe_cassette(),) + self.cassettes,
        )
        result = service.design_gene("ade6")
        guide = result.guides[0]
        annotation = service.annotate_guide(guide)

        choices = service.cassette_choices(
            (guide,),
            (annotation,),
            "ade6",
        )
        donors = service.disruption_donors(
            (guide,),
            (annotation,),
            choices,
            80,
            "ade6",
        )

        self.assertEqual((self.cassettes,), choices)
        self.assertEqual(10, len(donors[0]))
        self.assertEqual(183, donors[0][0].total_length)

    def test_noncoding_target(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("SPNCRNA.1438")
        annotations = service.annotate_guides(result.guides)
        choices = service.cassette_choices(
            result.guides,
            annotations,
            result.name,
        )
        donors = service.disruption_donors(
            result.guides,
            annotations,
            choices,
            80,
            result.name,
        )
        rows = annotation_rows(result.guides, annotations, result.name)
        template = (PROJECT_ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )
        page = render_design(
            result,
            annotations,
            template,
            cassette_choices=choices,
            disruption_donors=donors,
        )

        self.assertTrue(all(not group for group in choices))
        self.assertTrue(all(not group for group in donors))
        self.assertTrue(all(row["coding_target"] is False for row in rows))
        self.assertTrue(all(row["coding_strand"] is None for row in rows))
        self.assertIn("Selected target is non-coding.", page)
        self.assertIn('annotation.coding_target === false', page)

    def test_compact_web_data(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("ade6")
        annotations = service.annotate_guides(result.guides)
        choices = service.cassette_choices(
            result.guides,
            annotations,
            result.name,
        )
        donors = service.disruption_donors(
            result.guides,
            annotations,
            choices,
            80,
            result.name,
        )
        template = (PROJECT_ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )
        page = render_design(
            result,
            annotations,
            template,
            cassette_choices=choices,
            disruption_donors=donors,
        )

        def browser_data(name):
            start = page.index(f"var {name} = ") + len(f"var {name} = ")
            end = page.index(";", start)
            return page[start:end], json.loads(page[start:end])

        cassette_json, cassettes = browser_data("cassette_data")
        donor_json, arms = browser_data("donor_arms")
        cassette = cassettes["catalog"]["1"]
        donor = donors[0][0]
        insert = cassette["sequence"]
        if arms[0]["coding_strand"] == "-":
            insert = reverse_complement(insert)
        sequence = arms[0]["left_arm"] + insert + arms[0]["right_arm"]
        overlap_start = len(sequence) - arms[0]["hrrv_length"]
        browser_forward = sequence[:arms[0]["hrfw_length"]]
        browser_reverse = reverse_complement(sequence[overlap_start:])
        browser_overlap = sequence[
            overlap_start:arms[0]["hrfw_length"]
        ]

        self.assertEqual(donor.sequence, sequence)
        self.assertEqual(donor.reverse, reverse_complement(sequence))
        self.assertEqual(donor.oligos.forward, browser_forward)
        self.assertEqual(donor.oligos.reverse, browser_reverse)
        self.assertEqual(donor.oligos.overlap, browser_overlap)
        self.assertEqual(103, arms[0]["hrfw_length"])
        self.assertEqual(103, arms[0]["hrrv_length"])
        self.assertEqual(23, arms[0]["overlap_length"])
        self.assertEqual(183, arms[0]["hr_product_length"])
        self.assertEqual(10, len(cassettes["catalog"]))
        self.assertEqual(144, len(cassettes["choices"]))
        self.assertEqual(144, len(arms))
        self.assertLess(len(cassette_json.encode("utf-8")), 10000)
        self.assertLess(len(donor_json.encode("utf-8")), 50000)

    def test_web_section(self):
        annotations = GenomeAnnotations.from_files(
            DATA / "Schizosaccharomyces_pombe_all_chromosomes.gff3",
            DATA / "gene_viability.tsv",
        )
        annotation = annotations.annotate_cut("III", (1316791, 1316792))
        guide = ADE6_GUIDE
        primer = (
            guide,
            "FORWARD",
            "REVERSE",
            (1316795, 1316797),
            1,
            "TGG",
        )
        result = DesignResult.from_legacy(
            (
                [[guide, primer, 5, 1, 1, 1, 1, 1, 1]],
                ("HR_FORWARD", "HR_REVERSE", "DELETED_DNA"),
                [],
                "ade6",
                "III",
                "1316337",
                "1317995",
            )
        )
        template = (PROJECT_ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )

        cassettes = cassette_data((self.cassettes,))
        donors = tuple(
            build_donor(
                self.reference,
                (1316791, 1316792),
                cassette,
                "+",
                80,
            )
            for cassette in self.cassettes
        )
        arms = donor_data((donors,))
        guide_rows = annotation_rows(result.guides, (annotation,), "ade6")
        page = render_design(
            result,
            (annotation,),
            template,
            cassette_choices=(self.cassettes,),
            disruption_donors=(donors,),
        )

        self.assertEqual(10, len(cassettes["catalog"]))
        self.assertEqual(list(range(1, 11)), cassettes["choices"][0])
        self.assertEqual(
            "TGA* TGA* CGT GAT GAC CTA GTA",
            cassettes["catalog"]["1"]["frames"][0],
        )
        self.assertEqual(80, arms[0]["arm_length"])
        self.assertEqual(donors[0].left_arm, arms[0]["left_arm"])
        self.assertEqual(donors[0].right_arm, arms[0]["right_arm"])
        self.assertNotIn("sequence", arms[0])
        self.assertEqual("+", guide_rows[0]["coding_strand"])
        self.assertIn("Stop-Cassette Disruption Design", page)
        self.assertEqual(
            3,
            page.count('class="toggle_header project_workflow_header"'),
        )
        self.assertEqual(
            3,
            page.count('class="workflow_section project_workflow"'),
        )
        self.assertEqual(6, page.count('class="toggle_button"'))
        self.assertEqual(1, page.count('aria-expanded="true"'))
        self.assertEqual(5, page.count('aria-expanded="false"'))
        self.assertEqual(5, page.count('style="display: none;"'))
        self.assertIn(
            'aria-expanded="true" aria-controls="primers_table"',
            page,
        )
        for section_id in (
            "selected_guide_and_cut",
            "deletion_design",
            "spedit_oligos",
            "stop_cassette_design",
            "wild_type_restoration",
        ):
            self.assertIn(
                'aria-expanded="false" aria-controls="{}"'.format(
                    section_id
                ),
                page,
            )
        self.assertNotIn(
            'id="selected_guide_and_cut" '
            'class="workflow_section project_workflow"',
            page,
        )
        self.assertNotIn(
            'id="deletion_design" class="workflow_section project_workflow"',
            page,
        )
        self.assertIn('id="selected_guide_gc"', page)
        self.assertIn("candidate.gc_percent.toFixed(1)", page)
        gc_guidance = (
            "40&ndash;60% may be favorable but does not "
            "guarantee guide activity."
        )
        guidance_start = page.index(gc_guidance)
        self.assertEqual(1, page.count(gc_guidance))
        self.assertLess(page.index('id="selected_guide_gc"'), guidance_start)
        self.assertLess(guidance_start, page.index('id="selected_pam"'))
        for reference in (
            "https://www.nature.com/articles/srep19675",
            "https://academic.oup.com/nar/article/46/3/1375/4754467",
            "https://academic.oup.com/bioinformatics/article/36/9/2684/"
            "5714741",
        ):
            self.assertIn('href="{}"'.format(reference), page)
        self.assertEqual(3, page.count('rel="noopener noreferrer"'))
        self.assertIn("Guide GC:", page)
        self.assertIn('id="stop_cassette_menu"', page)
        self.assertIn("Reading frame 1:", page)
        self.assertIn("Reading frame 3:", page)
        self.assertNotIn("Reading frame 0:", page)
        self.assertIn("Disruption donor", page)
        self.assertIn("Oriented cassette (forward):", page)
        self.assertIn("Complete donor (forward):", page)
        self.assertNotIn("forward reference", page)
        self.assertIn(
            'class="l_field field_break">Oriented cassette (forward):',
            page,
        )
        self.assertIn(
            'class="l_field field_break">Right homology arm:',
            page,
        )
        self.assertIn('id="donor_total_length"', page)
        self.assertEqual(2, page.count("HR-template construction oligos"))
        self.assertIn('id="donor_hrfw"', page)
        self.assertIn('id="donor_hrrv"', page)
        self.assertIn('id="donor_hr_overlap"', page)
        self.assertIn('id="donor_hr_overlap_reverse"', page)
        self.assertIn('id="donor_hr_overlap_length"', page)
        self.assertIn('id="donor_hr_product_length"', page)
        self.assertIn("Overlap sequence (forward):", page)
        self.assertIn("Overlap sequence (reverse):", page)
        self.assertIn("rev_complement(oligos.overlap)", page)
        self.assertLess(
            page.index("Overlap sequence (forward):"),
            page.index("Overlap sequence (reverse):"),
        )
        self.assertLess(
            page.index("Overlap sequence (reverse):"),
            page.index("Overlap length:"),
        )
        self.assertIn("function hr_oligos(sequence, donor)", page)
        self.assertLess(
            page.index('id="donor_total_length"'),
            page.index("HR-template construction oligos"),
        )
        self.assertLess(
            page.index("HR-template construction oligos"),
            page.index("Insertion-checking primers"),
        )
        self.assertNotIn(donors[0].sequence, page)
        self.assertIn(
            "donor.left_arm + oriented_insert + donor.right_arm",
            page,
        )
        self.assertIn(FIRST, page)
        self.assertNotIn(OLD_FIRST, page)
        self.assertNotIn("computational candidate", page.lower())
        self.assertIn("Insertion-checking primers", page)
        self.assertIn("Edit-spanning PCR", page)
        self.assertIn(
            'id="junction_primer_details" class="annotation_details"',
            page,
        )
        self.assertIn("<summary>Junction-checking primers</summary>", page)
        self.assertIn("Left-junction PCR", page)
        self.assertIn("Right-junction PCR", page)
        details_start = page.index('<details id="junction_primer_details"')
        details_end = page.index("</details>", details_start)
        details = page[details_start:details_end]
        opening_tag = page[details_start:page.index(">", details_start)]
        self.assertNotIn(" open", opening_tag)
        self.assertIn("Left-junction PCR", details)
        self.assertIn("Right-junction PCR", details)
        self.assertLess(page.index("Edit-spanning PCR"), details_start)
        source = (
            "source: packaged PomBase GFF3 and gene viability snapshot"
        )
        source_start = page.index(source)
        self.assertEqual(1, page.count(source))
        self.assertLess(page.index("Cut-site annotation"), source_start)
        self.assertLess(
            source_start,
            page.index('id="cut_annotation_summary"'),
        )
        cut_details_start = page.index('<details id="cut_annotation_details"')
        cut_details_end = page.index("</details>", cut_details_start)
        self.assertNotIn(source, page[cut_details_start:cut_details_end])
        self.assertIn('id="insertion_primer_forward"', page)
        self.assertIn('id="left_junction_reverse"', page)
        self.assertIn('id="right_junction_forward"', page)
        self.assertIn("Expected WT product:", page)
        self.assertIn("Expected disrupted product:", page)
        self.assertIn('fetch("/insertion-primers?"', page)
        self.assertIn("cassette_id: cassette.id", page)
        self.assertIn("coding_strand: annotation.coding_strand", page)
        self.assertIn("var insertion_primer_cache = {};", page)
        self.assertIn('var key = index + ":" + cassette.id;', page)
        self.assertIn("update_insertion_primers(guide_number, cassette);", page)

        css = (PROJECT_ROOT / "css/crispr4p.css").read_text(encoding="utf-8")
        shell = (PROJECT_ROOT / "template/bahler_template.html").read_text(
            encoding="utf-8"
        )
        self.assertIn(".toggle_button::before {", css)
        self.assertIn(".gc_guidance {", css)
        self.assertIn("font-style: italic;", css)
        self.assertIn(
            '.toggle_button[aria-expanded="false"]::before {',
            css,
        )
        self.assertIn(
            'button.setAttribute("aria-expanded", expanded ? "false" : "true")',
            shell,
        )
        self.assertIn(
            'label === "Gene viability" && value === "inviable (essential)"',
            page,
        )
        self.assertIn('row.className += " annotation_warning";', page)
        self.assertIn(".project_workflow_header {", css)
        self.assertIn(".project_workflow .l_field {", css)
        self.assertIn(".project_workflow .annotation_details {", css)
        self.assertEqual(3, css.count("#ECE2CB"))


if __name__ == "__main__":
    unittest.main()
