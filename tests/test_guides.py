import unittest
from unittest.mock import patch

from crispr4p.crispr4p import NGG, PrimerDesign
from crispr4p.guides import is_match, match_guide


SEED = "AACCGGTTAACCGGTTAACC"


def replace_base(sequence, position):
    replacement = "T" if sequence[position] != "T" else "A"
    return sequence[:position] + replacement + sequence[position + 1:]


class TestGuideMatching(unittest.TestCase):
    def setUp(self):
        self.guide = NGG("query", 0, 1, SEED, "NGG")
        self.exact = NGG("I", 1, 1, SEED, "TGG")
        self.distal_mismatch = NGG(
            "I",
            2,
            1,
            replace_base(SEED, 0),
            "TGG",
        )
        self.middle_mismatch = NGG(
            "I",
            3,
            1,
            replace_base(SEED, 8),
            "TGG",
        )
        self.bucket = (
            self.exact,
            self.exact,
            self.distal_mismatch,
            self.middle_mismatch,
        )
        self.index = {SEED[-8:]: self.bucket}

    def test_exact_matching_preserves_legacy_cumulative_counts(self):
        returned_guide, table = match_guide(
            self.guide,
            self.index,
            n_mismatch=0,
        )

        self.assertIs(self.guide, returned_guide)
        self.assertEqual(
            {8: 4, 10: 3, 12: 2, 14: 2, 16: 2, 18: 2, 20: 1},
            {length: len(hits) for length, hits in table.items()},
        )
        self.assertEqual([1, 1, 2, 3], [hit.pos for hit in table[8]])
        self.assertIsNot(self.bucket, table[8])
        self.assertEqual(4, len(self.bucket))

    def test_one_mismatch_is_allowed_only_after_exact_eight_base_lookup(self):
        _, table = match_guide(
            self.guide,
            self.index,
            n_mismatch=1,
        )

        self.assertEqual(
            {8: 4, 10: 3, 12: 3, 14: 3, 16: 3, 18: 3, 20: 3},
            {length: len(hits) for length, hits in table.items()},
        )

    def test_sequence_comparison_preserves_mismatch_thresholds(self):
        self.assertFalse(is_match("AAAA", "AAAT", -1))
        self.assertFalse(is_match("AAAA", "AAAT", 0))
        self.assertTrue(is_match("AAAA", "AAAT", 1))
        self.assertTrue(is_match("AAAA", "AAAT", 2))

    def test_legacy_primer_design_methods_delegate_to_public_functions(self):
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.NGGs = self.index
        sentinel = (self.guide, {8: []})

        with patch(
            "crispr4p.crispr4p.match_guide",
            return_value=sentinel,
        ) as matcher:
            result = designer._single_table_worker(self.guide, 2)

        matcher.assert_called_once_with(self.guide, self.index, 2)
        self.assertIs(sentinel, result)
        self.assertEqual(
            is_match("AAAA", "AAAT", 1),
            PrimerDesign.genomeCompare("AAAA", "AAAT", 1),
        )


if __name__ == "__main__":
    unittest.main()
