import unittest
from dataclasses import FrozenInstanceError

from crispr4p.models import DesignResult


LEGACY_RESULT = (
    [["guide"]],
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


if __name__ == "__main__":
    unittest.main()
