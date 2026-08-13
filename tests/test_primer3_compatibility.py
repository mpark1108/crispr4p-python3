import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

import crispr4p.crispr4p as crp


class TestPrimer3Compatibility(unittest.TestCase):
    def test_prefers_current_primer3_api(self) -> None:
        expected = {"result": "modern"}
        modern = Mock(return_value=expected)
        legacy = Mock()
        bindings = SimpleNamespace(
            design_primers=modern,
            designPrimers=legacy,
        )
        sequence_args = {"SEQUENCE_ID": "test"}
        global_args = {"PRIMER_OPT_SIZE": 20}

        with patch.object(crp, "primer3", bindings):
            actual = crp._design_primers(sequence_args, global_args)

        self.assertIs(expected, actual)
        modern.assert_called_once_with(sequence_args, global_args)
        legacy.assert_not_called()

    def test_falls_back_for_older_primer3_releases(self) -> None:
        expected = {"result": "legacy"}
        legacy = Mock(return_value=expected)
        bindings = SimpleNamespace(designPrimers=legacy)
        sequence_args = {"SEQUENCE_ID": "test"}
        global_args = {"PRIMER_OPT_SIZE": 20}

        with patch.object(crp, "primer3", bindings):
            actual = crp._design_primers(sequence_args, global_args)

        self.assertIs(expected, actual)
        legacy.assert_called_once_with(sequence_args, global_args)


if __name__ == "__main__":
    unittest.main()
