import unittest
from pathlib import Path

from crispr4p.models import DesignResult
from crispr4p.service import Crispr4pService


LEGACY_RESULT = (
    [["guide"]],
    ("forward", "reverse", "deleted"),
    [{"primer": "checking"}],
    "ade6",
    "III",
    "1316337",
    "1317995",
)


class RecordingDesigner:
    instances = []

    def __init__(self, *args, **kwargs):
        self.constructor_args = args
        self.constructor_kwargs = kwargs
        self.run_web_calls = []
        self.__class__.instances.append(self)

    def runWeb(self, **kwargs):
        self.run_web_calls.append(kwargs)
        return LEGACY_RESULT


class TestCrispr4pService(unittest.TestCase):
    def setUp(self) -> None:
        RecordingDesigner.instances = []
        self.service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            precomputed_folder="cache",
            designer_factory=RecordingDesigner,
        )

    def test_gene_design_delegates_without_changing_legacy_result(self) -> None:
        result = self.service.design_gene("ade6", n_mismatch=2)

        self.assertIsInstance(result, DesignResult)
        self.assertIs(LEGACY_RESULT, result.to_legacy())
        self.assertIs(LEGACY_RESULT[0], result.guide_table)
        self.assertEqual("ade6", result.name)
        designer = RecordingDesigner.instances[0]
        self.assertEqual(
            ("genome.fa", "coordinates.txt", "synonyms.txt"),
            designer.constructor_args,
        )
        self.assertEqual(
            {"precomputed_folder": "cache"},
            designer.constructor_kwargs,
        )
        self.assertEqual(
            [
                {
                    "name": "ade6",
                    "cr": None,
                    "start": None,
                    "end": None,
                    "strand": None,
                    "nMismatch": 2,
                }
            ],
            designer.run_web_calls,
        )

    def test_region_design_delegates_all_coordinates_and_options(self) -> None:
        result = self.service.design_region(
            "III",
            "100",
            "200",
            strand="-1",
            n_mismatch=1,
        )

        self.assertIsInstance(result, DesignResult)
        self.assertIs(LEGACY_RESULT, result.to_legacy())
        self.assertEqual("III", result.chromosome)
        self.assertEqual("1316337", result.start)
        self.assertEqual("1317995", result.end)
        self.assertEqual(
            [
                {
                    "name": None,
                    "cr": "III",
                    "start": "100",
                    "end": "200",
                    "strand": "-1",
                    "nMismatch": 1,
                }
            ],
            RecordingDesigner.instances[0].run_web_calls,
        )

    def test_each_request_gets_fresh_query_state(self) -> None:
        self.service.design_gene("ade6")
        self.service.design_region("III", "100", "200")

        self.assertEqual(2, len(RecordingDesigner.instances))
        self.assertIsNot(
            RecordingDesigner.instances[0],
            RecordingDesigner.instances[1],
        )

    def test_project_data_factory_uses_packaged_reference_files(self) -> None:
        service = Crispr4pService.from_project_data(
            precomputed_folder="cache",
            designer_factory=RecordingDesigner,
        )
        service.design_gene("ade6")

        sequence_file, coordinates_file, synonyms_file = (
            Path(value)
            for value in RecordingDesigner.instances[0].constructor_args
        )
        self.assertEqual("data", sequence_file.parent.name)
        self.assertEqual(
            "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            sequence_file.name,
        )
        self.assertEqual("COORDINATES.txt", coordinates_file.name)
        self.assertEqual("SYNONIMS.txt", synonyms_file.name)
        self.assertTrue(sequence_file.is_file())
        self.assertTrue(coordinates_file.is_file())
        self.assertTrue(synonyms_file.is_file())


if __name__ == "__main__":
    unittest.main()
