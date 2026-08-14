import unittest

from crispr4p.spedit import (
    has_bsai,
    make_oligos,
    reverse_complement,
)


# Published guide and oligo sequences from the SpEDIT paper's Table 2.
PUBLISHED_SPEDIT_OLIGOS = {
    "ade6": (
        "TTGATAGCAACAGTGGCGAC",
        "CTAGAGGTCTCGGACTTTGATAGCAACAGTGGCGACGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGTCGCCACTGTTGCTATCAAAGTCCGAGACCTCTAG",
    ),
    "ura4": (
        "CCTTGTATAATACCCTCGCC",
        "CTAGAGGTCTCGGACTCCTTGTATAATACCCTCGCCGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGGCGAGGGTATTATACAAGGAGTCCGAGACCTCTAG",
    ),
    "meu27": (
        "TATTAGCCTTTGAAGGATTT",
        "CTAGAGGTCTCGGACTTATTAGCCTTTGAAGGATTTGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACAAATCCTTCAAAGGCTAATAAGTCCGAGACCTCTAG",
    ),
    "clr5": (
        "AGCTTGTGGCTGACCGTTAA",
        "CTAGAGGTCTCGGACTAGCTTGTGGCTGACCGTTAAGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACTTAACGGTCAGCCACAAGCTAGTCCGAGACCTCTAG",
    ),
    "cup1:4xtetO": (
        "ATTTCTTTTGCTTTACGGTC",
        "CTAGAGGTCTCGGACTATTTCTTTTGCTTTACGGTCGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGACCGTAAAGCAAAAGAAATAGTCCGAGACCTCTAG",
    ),
    "cup1-GFP": (
        "GCTCAGGCTAAACGTCGGAA",
        "CTAGAGGTCTCGGACTGCTCAGGCTAAACGTCGGAAGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACTTCCGACGTTTAGCCTGAGCAGTCCGAGACCTCTAG",
    ),
    "epe1": (
        "GGACTTTTAAGATGGATTCC",
        "CTAGAGGTCTCGGACTGGACTTTTAAGATGGATTCCGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGGAATCCATCTTAAAAGTCCAGTCCGAGACCTCTAG",
    ),
}


class TestSpeditOligos(unittest.TestCase):
    def test_published_spedit_examples(self) -> None:
        """Generated oligos must match all examples published by SpEDIT."""
        for name, (
            guide,
            expected_forward,
            expected_reverse,
        ) in PUBLISHED_SPEDIT_OLIGOS.items():
            with self.subTest(name=name):
                forward, reverse = make_oligos(guide)

                self.assertEqual(expected_forward, forward)
                self.assertEqual(expected_reverse, reverse)

    def test_all_published_oligos_are_52_nt(self) -> None:
        """Every generated SpEDIT oligo must contain exactly 52 bases."""
        for name, (guide, _, _) in PUBLISHED_SPEDIT_OLIGOS.items():
            with self.subTest(name=name):
                forward, reverse = make_oligos(guide)

                self.assertEqual(52, len(forward))
                self.assertEqual(52, len(reverse))

    def test_reverse_oligo_is_reverse_complement(self) -> None:
        """The reverse oligo must complement the complete forward oligo."""
        for name, (guide, _, _) in PUBLISHED_SPEDIT_OLIGOS.items():
            with self.subTest(name=name):
                forward, reverse = make_oligos(guide)

                self.assertEqual(reverse_complement(forward), reverse)

    def test_normalizes_lowercase_and_whitespace(self) -> None:
        """Lowercase guides with surrounding whitespace should be accepted."""
        expected = make_oligos(
            "TTGATAGCAACAGTGGCGAC"
        )

        actual = make_oligos(
            "  ttgatagcaacagtggcgac\n"
        )

        self.assertEqual(expected, actual)

    def test_rejects_short_guide(self) -> None:
        """Guides shorter than 20 nt must be rejected."""
        with self.assertRaisesRegex(ValueError, "20-nt"):
            make_oligos("ACGT")

    def test_rejects_long_guide(self) -> None:
        """Guides longer than 20 nt must be rejected."""
        with self.assertRaisesRegex(ValueError, "20-nt"):
            make_oligos("A" * 21)

    def test_rejects_invalid_nucleotide(self) -> None:
        """Only A, C, G, and T are accepted."""
        with self.assertRaisesRegex(
            ValueError,
            "invalid nucleotide",
        ):
            make_oligos("A" * 19 + "N")

    def test_rejects_non_string_guide(self) -> None:
        """Non-string guide values must be rejected."""
        with self.assertRaisesRegex(
            TypeError,
            "must be a string",
        ):
            make_oligos(12345)  # type: ignore[arg-type]

    def test_detects_forward_bsai_site(self) -> None:
        """Detect the forward BsaI recognition sequence GGTCTC."""
        self.assertTrue(
            has_bsai(
                "AAAAAGGTCTCAAAAAAAAA"
            )
        )

    def test_detects_reverse_bsai_site(self) -> None:
        """Detect the reverse-complement BsaI sequence GAGACC."""
        self.assertTrue(
            has_bsai(
                "AAAAAGAGACCAAAAAAAAA"
            )
        )

    def test_bsai_detection_normalizes_input(self) -> None:
        """BsaI detection should normalize case and whitespace."""
        self.assertTrue(
            has_bsai(
                "  aaaaaggtctcaaaaaaaaa\n"
            )
        )

    def test_accepts_guide_without_internal_bsai_site(self) -> None:
        """A standard published guide should not trigger the warning."""
        self.assertFalse(
            has_bsai(
                "TTGATAGCAACAGTGGCGAC"
            )
        )


if __name__ == "__main__":
    unittest.main()
