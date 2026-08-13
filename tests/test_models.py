import unittest
from dataclasses import FrozenInstanceError

from crispr4p.models import (
    DesignResult,
    GuideCandidate,
    OligoAnalysisResult,
    OligoMatch,
)


LEGACY_RESULT = (
    [["guide"]],
    ("forward", "reverse", "deleted"),
    [{"primer": "checking"}],
    "ade6",
    "III",
    "1316337",
    "1317995",
)

GUIDE = "ACATTGGCTTACGACGGTCG"
GUIDE_PRIMER = (
    GUIDE,
    "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
    "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
    (1316795, 1316797),
    1,
    "TGG",
)
GUIDE_ROW = [GUIDE, GUIDE_PRIMER, 5, 1, 1, 1, 1, 1, 1]
VALID_LEGACY_RESULT = (
    [GUIDE_ROW],
    ("forward", "reverse", "deleted"),
    [{"primer": "checking"}],
    "ade6",
    "III",
    "1316337",
    "1317995",
)


class TestDesignResult(unittest.TestCase):
    def test_names_every_legacy_top_level_field(self) -> None:
        result = DesignResult.from_legacy(LEGACY_RESULT)

        self.assertIs(LEGACY_RESULT[0], result.guide_table)
        self.assertIs(LEGACY_RESULT[1], result.hr_dna)
        self.assertIs(LEGACY_RESULT[2], result.checking_primers)
        self.assertEqual("ade6", result.name)
        self.assertEqual("III", result.chromosome)
        self.assertEqual("1316337", result.start)
        self.assertEqual("1317995", result.end)

    def test_adapter_returns_exact_original_legacy_tuple(self) -> None:
        result = DesignResult.from_legacy(LEGACY_RESULT)

        self.assertIs(LEGACY_RESULT, result.to_legacy())

    def test_result_is_immutable(self) -> None:
        result = DesignResult.from_legacy(LEGACY_RESULT)

        with self.assertRaises(FrozenInstanceError):
            result.name = "ura4"

    def test_rejects_non_tuple_legacy_result(self) -> None:
        with self.assertRaisesRegex(TypeError, "must be a tuple"):
            DesignResult.from_legacy(list(LEGACY_RESULT))

    def test_rejects_wrong_legacy_tuple_length(self) -> None:
        with self.assertRaisesRegex(ValueError, "seven items"):
            DesignResult.from_legacy(LEGACY_RESULT[:-1])

    def test_exposes_immutable_named_guide_candidates(self) -> None:
        result = DesignResult.from_legacy(VALID_LEGACY_RESULT)

        self.assertIsInstance(result.guides, tuple)
        self.assertEqual(1, len(result.guides))
        guide = result.guides[0]
        self.assertIsInstance(guide, GuideCandidate)
        self.assertEqual("III", guide.chromosome)
        self.assertEqual(GUIDE, guide.seed)
        self.assertEqual(GUIDE_PRIMER[1], guide.forward_cloning_oligo)
        self.assertEqual(GUIDE_PRIMER[2], guide.reverse_cloning_oligo)
        self.assertEqual((1316795, 1316797), guide.pam_coordinates)
        self.assertEqual((1316791, 1316792), guide.cut_coordinates)
        self.assertEqual(1, guide.strand)
        self.assertEqual("TGG", guide.pam)
        self.assertEqual(
            {8: 5, 10: 1, 12: 1, 14: 1, 16: 1, 18: 1, 20: 1},
            dict(guide.match_counts),
        )
        self.assertIs(GUIDE_ROW, guide.to_legacy())
        with self.assertRaises(TypeError):
            guide.match_counts[8] = 99

    def test_typed_guides_are_a_snapshot_of_mutable_legacy_rows(self) -> None:
        local_row = list(GUIDE_ROW)
        legacy_result = (local_row,), *VALID_LEGACY_RESULT[1:]
        result = DesignResult.from_legacy(legacy_result)
        guide = result.guides[0]

        local_row[2] = 99

        self.assertEqual(5, guide.match_counts[8])
        self.assertEqual(99, result.guide_table[0][2])

    def test_malformed_nested_rows_remain_accepted_until_typed_access(self):
        result = DesignResult.from_legacy(LEGACY_RESULT)

        self.assertIs(LEGACY_RESULT, result.to_legacy())
        with self.assertRaisesRegex(ValueError, "nine items"):
            _ = result.guides


class TestOligoAnalysisResult(unittest.TestCase):
    def make_result(self):
        return OligoAnalysisResult(
            oligo_sequence="ACATTGGCTTACGACGGTCG",
            seed="ACATTGGCTTACGACGGTCG",
            pam="NGG",
            n_mismatch=0,
            spedit_forward="forward",
            spedit_reverse="reverse",
            has_internal_bsai=False,
            match_counts={8: 5, 20: 1},
            full_matches=[
                OligoMatch(
                    chromosome="III",
                    pam_coordinates=[1316795, 1316797],
                    cut_coordinates=[1316791, 1316792],
                    strand=1,
                    seed="ACATTGGCTTACGACGGTCG",
                    pam="TGG",
                )
            ],
        )

    def test_nested_collections_are_immutable_snapshots(self) -> None:
        result = self.make_result()

        self.assertIsInstance(result.full_matches, tuple)
        self.assertEqual((1316795, 1316797), result.full_matches[0].pam_coordinates)
        self.assertEqual((1316791, 1316792), result.full_matches[0].cut_coordinates)
        with self.assertRaises(TypeError):
            result.match_counts[8] = 99

    def test_top_level_result_is_immutable(self) -> None:
        result = self.make_result()

        with self.assertRaises(FrozenInstanceError):
            result.seed = "A" * 20


if __name__ == "__main__":
    unittest.main()
