import json
import unittest
from pathlib import Path
from types import SimpleNamespace

from crispr4p.annotations import GeneAnnotation, GenomeAnnotations
from crispr4p.disruption import (
    StopCassette,
    has_junction_pam,
    load_cassettes,
    recut_sites,
    target_strand,
)
from crispr4p.models import DesignResult
from crispr4p.resources import read_fasta
from crispr4p.service import Crispr4pService
from crispr4p.spedit import reverse_complement
from crispr4p.web_views import annotation_rows, cassette_rows, render_design


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA = PROJECT_ROOT / "data"
CANDIDATES = DATA / "stop_cassettes.json"
FASTA = DATA / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa"
FIRST = "TGATGACCTAGTGACCTAGTAGG"
OLD_FIRST = "TAGTAGCCTAGTGACCTAGTAGG"
ADE6_GUIDE = "ACATTGGCTTACGACGGTCG"
ADE6_REVERSE_GUIDE = "GTGGCGACAGGGACACCTCG"


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
        cls.reference = read_fasta(FASTA)["III"].sequence

    def test_packaged_candidates(self):
        self.assertEqual(tuple(range(1, 11)), tuple(c.id for c in self.cassettes))
        self.assertFalse(any(has_junction_pam(c.sequence) for c in self.cassettes))
        self.assertEqual(FIRST, self.cassettes[0].sequence)
        self.assertEqual(
            "TGATGACTTAGTGACCTAGTAGG",
            self.cassettes[-1].sequence,
        )

        data = json.loads(CANDIDATES.read_text(encoding="utf-8"))
        self.assertEqual(1, data["version"])
        self.assertEqual("ASM294v2.26", data["assembly"])

    def test_candidate_details(self):
        cassette = self.cassettes[0]

        self.assertEqual("TGATGACCTAGTGACCTAGT", cassette.guide)
        self.assertEqual("AGG", cassette.pam)
        self.assertEqual(45.0, cassette.gc_percent)
        self.assertEqual(3, len(cassette.frames))
        self.assertEqual(
            "CCTACTAGGTCACTAGGTCATCA",
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

    def test_target_strand(self):
        plus = gene("plus", "target", "+")
        minus = gene("minus", "other", "-")

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

        self.assertEqual((self.cassettes,), choices)

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

        rows = cassette_rows((self.cassettes,))
        guide_rows = annotation_rows(result.guides, (annotation,), "ade6")
        page = render_design(
            result,
            (annotation,),
            template,
            cassette_choices=(self.cassettes,),
        )

        self.assertEqual(1, len(rows))
        self.assertEqual(10, len(rows[0]))
        self.assertEqual(
            "TGA* TGA* CCT AGT GAC CTA GTA",
            rows[0][0]["frames"][0],
        )
        self.assertEqual("+", guide_rows[0]["coding_strand"])
        self.assertIn("Stop-Cassette Disruption Design", page)
        self.assertIn('id="stop_cassette_menu"', page)
        self.assertIn("Reading frame 1:", page)
        self.assertIn("Reading frame 3:", page)
        self.assertNotIn("Reading frame 0:", page)
        self.assertIn(FIRST, page)
        self.assertNotIn(OLD_FIRST, page)
        self.assertNotIn("Cassette sequence (forward reference)", page)
        self.assertNotIn("computational candidate", page.lower())


if __name__ == "__main__":
    unittest.main()
