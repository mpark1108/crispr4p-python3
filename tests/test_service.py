import unittest
from pathlib import Path
from types import SimpleNamespace

from crispr4p.models import DesignResult, OligoAnalysisResult
from crispr4p.service import Crispr4pService, OligoLengthError


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


class RecordingDesigner:
    instances = []

    def __init__(self, *args, **kwargs):
        self.constructor_args = args
        self.constructor_kwargs = kwargs
        self.run_web_calls = []
        self.index_calls = 0
        self.oligo_calls = []
        self.__class__.instances.append(self)

    def runWeb(self, **kwargs):
        self.run_web_calls.append(kwargs)
        return LEGACY_RESULT

    def getNGGsFromGenome(self):
        self.index_calls += 1

    def _single_table_worker(self, query, n_mismatch):
        self.oligo_calls.append((query, n_mismatch))
        hit = SimpleNamespace(
            chromosome="III",
            strand=1,
            seed=GUIDE,
            pam="TGG",
        )
        return query, {
            8: [hit] * 5,
            10: [hit],
            12: [hit],
            14: [hit],
            16: [hit],
            18: [hit],
            20: [hit],
        }

    def getOligoHitCoordinates(self, match):
        return (1316795, 1316797), (1316791, 1316792)


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
        self.service.analyze_oligo(GUIDE)

        self.assertEqual(3, len(RecordingDesigner.instances))
        self.assertIsNot(
            RecordingDesigner.instances[0],
            RecordingDesigner.instances[1],
        )
        self.assertIsNot(
            RecordingDesigner.instances[1],
            RecordingDesigner.instances[2],
        )

    def test_oligo_analysis_returns_structured_existing_results(self) -> None:
        result = self.service.analyze_oligo(
            f"  {GUIDE.lower()}  ",
            n_mismatch=1,
        )

        self.assertIsInstance(result, OligoAnalysisResult)
        self.assertEqual(GUIDE, result.oligo_sequence)
        self.assertEqual(GUIDE, result.seed)
        self.assertEqual("NGG", result.pam)
        self.assertEqual(1, result.n_mismatch)
        self.assertEqual(5, result.match_counts[8])
        self.assertEqual(1, result.match_counts[20])
        self.assertEqual(1, len(result.full_matches))
        self.assertEqual("III", result.full_matches[0].chromosome)
        self.assertEqual(
            (1316795, 1316797),
            result.full_matches[0].pam_coordinates,
        )
        self.assertEqual(
            (1316791, 1316792),
            result.full_matches[0].cut_coordinates,
        )
        designer = RecordingDesigner.instances[0]
        self.assertEqual(1, designer.index_calls)
        query, n_mismatch = designer.oligo_calls[0]
        self.assertEqual(GUIDE, query.seed)
        self.assertEqual("NGG", query.pam)
        self.assertEqual(1, n_mismatch)

    def test_oligo_analysis_preserves_supplied_pam(self) -> None:
        result = self.service.analyze_oligo(GUIDE + "TGG")

        self.assertEqual(GUIDE + "TGG", result.oligo_sequence)
        self.assertEqual(GUIDE, result.seed)
        self.assertEqual("TGG", result.pam)

    def test_oligo_analysis_rejects_invalid_length_before_loading_data(self) -> None:
        with self.assertRaises(OligoLengthError) as raised:
            self.service.analyze_oligo("A" * 19)

        self.assertEqual(19, raised.exception.sequence_length)
        self.assertEqual([], RecordingDesigner.instances)

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
