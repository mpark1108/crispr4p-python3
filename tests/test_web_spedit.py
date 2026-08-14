import unittest

from crispr4p.web_views import spedit_rows


def make_crispr4p_row(guide: str):
    primer_tuple = (
        guide,
        "LEGACY_FORWARD",
        "LEGACY_REVERSE",
        (100, 103),
        1,
        "TGG",
    )

    return [guide, primer_tuple]


class TestSpeditCandidateData(unittest.TestCase):
    def test_preserves_candidate_order(self) -> None:
        first = "ACATTGGCTTACGACGGTCG"
        second = "TTGATAGCAACAGTGGCGAC"

        candidates = spedit_rows(
            [
                make_crispr4p_row(first),
                make_crispr4p_row(second),
            ]
        )

        self.assertEqual(first, candidates[0]["guide"])
        self.assertEqual(second, candidates[1]["guide"])

    def test_generates_oligos_for_selected_candidate(self) -> None:
        guide = "TTGATAGCAACAGTGGCGAC"

        candidates = spedit_rows(
            [make_crispr4p_row(guide)]
        )

        self.assertEqual(
            (
                "CTAGAGGTCTCGGACT"
                "TTGATAGCAACAGTGGCGAC"
                "GTTTCGAGACCCTTCC"
            ),
            candidates[0]["forward"],
        )

        self.assertEqual(
            (
                "GGAAGGGTCTCGAAAC"
                "GTCGCCACTGTTGCTATCAA"
                "AGTCCGAGACCTCTAG"
            ),
            candidates[0]["reverse"],
        )

    def test_internal_bsai_warning_is_aligned(self) -> None:
        safe = "TTGATAGCAACAGTGGCGAC"
        unsafe = "TTTTGAATGGTCTCAGTTGT"

        candidates = spedit_rows(
            [
                make_crispr4p_row(safe),
                make_crispr4p_row(unsafe),
            ]
        )

        self.assertFalse(candidates[0]["has_internal_bsai"])
        self.assertTrue(candidates[1]["has_internal_bsai"])


if __name__ == "__main__":
    unittest.main()
