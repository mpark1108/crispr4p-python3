import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from crispr4p.annotations import GenomeAnnotations
from crispr4p.web_views import annotation_rows
from tests.test_annotations import GFF, VIABILITY


def guide(cut_left):
    return SimpleNamespace(
        chromosome="I",
        seed="ACATTGGCTTACGACGGTCG",
        gc_percent=55.0,
        pam="TGG",
        pam_coordinates=(cut_left + 4, cut_left + 6),
        cut_coordinates=(cut_left, cut_left + 1),
        strand=1,
    )


class WebAnnotationDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory()
        directory = Path(cls.temporary_directory.name)
        gff_path = directory / "annotations.gff3"
        viability_path = directory / "gene_viability.tsv"
        gff_path.write_text(GFF, encoding="utf-8")
        viability_path.write_text(VIABILITY, encoding="utf-8")
        cls.annotations = GenomeAnnotations.from_files(
            gff_path,
            viability_path,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def candidate(self, cut_left, target_name=None):
        annotation = self.annotations.annotate_cut(
            "I",
            (cut_left, cut_left + 1),
        )
        return annotation_rows(
            (guide(cut_left),),
            (annotation,),
            target_name=target_name,
        )[0]

    def test_primary_target_and_non_coding_overlap_are_grouped_by_gene(self):
        candidate = self.candidate(200, target_name="plus_gene")

        self.assertEqual("I", candidate["chromosome"])
        self.assertEqual(55.0, candidate["gc_percent"])
        self.assertEqual([204, 206], candidate["pam_coordinates"])
        self.assertEqual([200, 201], candidate["cut_coordinates"])
        self.assertEqual("+", candidate["strand"])
        self.assertEqual(2, candidate["gene_count"])
        self.assertEqual("plus", candidate["genes"][0]["gene_id"])

        genes = {gene["gene_id"]: gene for gene in candidate["genes"]}
        primary = genes["plus"]
        self.assertEqual("Primary target", primary["role"])
        self.assertTrue(primary["is_primary"])
        self.assertEqual("viable (non-essential)", primary["viability"])
        self.assertEqual("coding sequence (CDS)", primary["contexts"][0]["region"])
        self.assertEqual(
            {"base": 81, "total": 200, "percent": 40.5},
            primary["contexts"][0]["cds"],
        )

        overlap = genes["antisense"]
        self.assertEqual("Additional overlap", overlap["role"])
        self.assertFalse(overlap["is_protein_coding"])
        self.assertNotIn("warning", overlap)
        self.assertEqual("condition-dependent", overlap["viability"])
        self.assertEqual("non-coding exon", overlap["contexts"][0]["region"])

        # The view model must remain directly safe to serialize for JavaScript.
        json.dumps(candidate)

    def test_additional_protein_coding_overlap_remains_informational(self):
        candidate = self.candidate(200, target_name="antisense_gene")
        genes = {gene["gene_id"]: gene for gene in candidate["genes"]}

        self.assertEqual("Primary target", genes["antisense"]["role"])
        self.assertEqual("Additional overlap", genes["plus"]["role"])
        self.assertTrue(genes["plus"]["is_protein_coding"])
        self.assertNotIn("warning", genes["plus"])

    def test_feature_boundary_is_explicit_in_summary_and_details(self):
        candidate = self.candidate(219, target_name="plus_gene")
        primary = next(gene for gene in candidate["genes"] if gene["is_primary"])
        context = primary["contexts"][0]

        self.assertTrue(context["crosses_boundary"])
        self.assertEqual(
            "coding sequence (CDS) / intron boundary",
            context["region"],
        )
        self.assertIsNone(context["block"])
        self.assertEqual("coding sequence (CDS)", context["left_block"]["label"])
        self.assertEqual("intron", context["right_block"]["label"])

    def test_intergenic_result_includes_nearest_genes(self):
        candidate = self.candidate(550)

        self.assertTrue(candidate["is_intergenic"])
        self.assertEqual([], candidate["genes"])
        self.assertEqual(
            ["plus", "minus"],
            [gene["gene_id"] for gene in candidate["nearest_genes"]],
        )
        self.assertEqual(
            ["viable (non-essential)", "inviable (essential)"],
            [gene["viability"] for gene in candidate["nearest_genes"]],
        )

    def test_candidate_and_annotation_counts_must_remain_aligned(self):
        with self.assertRaisesRegex(ValueError, "counts must match"):
            annotation_rows((guide(200),), ())


if __name__ == "__main__":
    unittest.main()
