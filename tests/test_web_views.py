import hashlib
import unittest
from pathlib import Path

from crispr4p.models import DesignResult, OligoAnalysisResult
from crispr4p.web_views import (
    OligoMatchView,
    render_design_result,
    render_execution_error,
    render_missing_query_error,
    render_oligo_length_error,
    render_oligo_result,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GUIDE = "ACATTGGCTTACGACGGTCG"
SPEDIT_FORWARD = (
    "CTAGAGGTCTCGGACTACATTGGCTTACGACGGTCGGTTTCGAGACCCTTCC"
)
SPEDIT_REVERSE = (
    "GGAAGGGTCTCGAAACCGACCGTCGTAAGCCAATGTAGTCCGAGACCTCTAG"
)


class TestErrorViews(unittest.TestCase):
    def test_missing_query_error_is_unchanged(self) -> None:
        self.assertEqual(
            '<font color="red"><h3>Error: Please fill either Name, '
            'Coordinates, or Oligo Sequence</h3></font>',
            render_missing_query_error(),
        )

    def test_execution_error_is_unchanged(self) -> None:
        self.assertEqual(
            '<font color="red"><h3>ERROR during execution: failure</h3></font>',
            render_execution_error(ValueError("failure")),
        )

    def test_oligo_length_error_is_unchanged(self) -> None:
        self.assertEqual(
            '<font color="red"><h3>Error: Oligo sequence must be 20 bp '
            '(seed only) or 23 bp (seed + PAM). Current length: 19</h3></font>',
            render_oligo_length_error(19),
        )


class TestResultViews(unittest.TestCase):
    def test_oligo_renderer_preserves_exact_html(self) -> None:
        analysis = OligoAnalysisResult(
            oligo_sequence=GUIDE,
            seed=GUIDE,
            pam="NGG",
            n_mismatch=0,
            spedit_forward=SPEDIT_FORWARD,
            spedit_reverse=SPEDIT_REVERSE,
            has_internal_bsai=False,
            match_counts={8: 5, 10: 1, 12: 1, 14: 1, 16: 1, 18: 1, 20: 1},
            full_matches=[
                OligoMatchView(
                    chromosome="III",
                    pam_coordinates=(1316795, 1316797),
                    cut_coordinates=(1316791, 1316792),
                    strand=1,
                    seed=GUIDE,
                    pam="TGG",
                )
            ],
        )
        result = render_oligo_result(analysis)

        self.assertEqual(
            "dd389db45eaca6b3babe396638609b0752b865f5638c14c32aa6c0e7e4c4082e",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )

    def test_design_renderer_preserves_exact_html(self) -> None:
        primer_tuple = (
            GUIDE,
            "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
            "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
            (1316795, 1316797),
            1,
            "TGG",
        )
        legacy_result = (
            [[GUIDE, primer_tuple, 5, 1, 1, 1, 1, 1, 1]],
            ("HR_FORWARD", "HR_REVERSE", "DELETED_DNA"),
            [
                {
                    "PRIMER_LEFT_0_SEQUENCE": "LEFT",
                    "PRIMER_LEFT_0_TM": 59.4,
                    "PRIMER_RIGHT_0_SEQUENCE": "RIGHT",
                    "PRIMER_RIGHT_0_TM": 60.1,
                    "PRIMER_PAIR_0_PRODUCT_SIZE": 255,
                    "negative_result": 1913,
                }
            ],
            "ade6",
            "III",
            "1316337",
            "1317995",
        )
        template_text = (
            PROJECT_ROOT / "template" / "container_table.html"
        ).read_text(encoding="utf-8")

        result = render_design_result(
            DesignResult.from_legacy(legacy_result),
            template_text,
        )

        self.assertEqual(
            "93376b12610f2f10077b8c7b6bd80992a0d9cb65fdb93aed927026330e11ea3a",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()
