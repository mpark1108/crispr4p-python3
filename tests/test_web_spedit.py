import unittest

from webapp import build_spedit_table


def make_crispr4p_row(guide: str):
    """Create the minimum CRISPR4P result-row structure needed by the view."""
    primer_tuple = (
        guide,
        "LEGACY_FORWARD",
        "LEGACY_REVERSE",
        (100, 103),
        1,
        "TGG",
    )

    return [guide, primer_tuple]


class TestSpeditWebOutput(unittest.TestCase):
    def test_table_contains_published_ade6_oligos(self) -> None:
        guide = "TTGATAGCAACAGTGGCGAC"

        output = build_spedit_table(
            [make_crispr4p_row(guide)]
        )

        self.assertIn(guide, output)

        self.assertIn(
            "CTAGAGGTCTCGGACT"
            "TTGATAGCAACAGTGGCGAC"
            "GTTTCGAGACCCTTCC",
            output,
        )

        self.assertIn(
            "GGAAGGGTCTCGAAAC"
            "GTCGCCACTGTTGCTATCAA"
            "AGTCCGAGACCTCTAG",
            output,
        )

    def test_table_labels_sequences_as_52_nt(self) -> None:
        output = build_spedit_table(
            [make_crispr4p_row("TTGATAGCAACAGTGGCGAC")]
        )

        self.assertIn("SpEDIT forward oligo, 52 nt", output)
        self.assertIn("SpEDIT reverse oligo, 52 nt", output)

    def test_table_warns_about_internal_bsai_site(self) -> None:
        guide = "AAAAAGGTCTCAAAAAAAAA"

        output = build_spedit_table(
            [make_crispr4p_row(guide)]
        )

        self.assertIn(
            "Warning: internal BsaI site in guide",
            output,
        )

    def test_empty_results_are_handled(self) -> None:
        output = build_spedit_table([])

        self.assertIn(
            "No guide candidates were available",
            output,
        )


if __name__ == "__main__":
    unittest.main()
