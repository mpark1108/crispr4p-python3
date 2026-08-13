import contextlib
import hashlib
import io
import os
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from crispr4p.crispr4p import PrimerDesign


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
FASTA = DATA_DIR / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa"
COORDINATES = DATA_DIR / "COORDINATES.txt"
SYNONIMS = DATA_DIR / "SYNONIMS.txt"

ADE6_GUIDE = "ACATTGGCTTACGACGGTCG"
ADE6_COORDINATES = ("III", "1316337", "1317995", "1")


class TestCoreCharacterization(unittest.TestCase):
    """Record outputs that the boundary refactor must continue to produce."""

    @classmethod
    def setUpClass(cls) -> None:
        # Build the real genome PAM index once for this test class. Calling
        # run_ avoids reading or writing the application's pickle cache.
        cls.designer = PrimerDesign(
            str(FASTA),
            str(COORDINATES),
            str(SYNONIMS),
            regression=True,
        )
        chromosome, start, end, _ = ADE6_COORDINATES
        (
            cls.table,
            cls.hr_dna,
            cls.checking_primers,
            cls.guide_matches,
        ) = cls.designer.run_(
            chromosome,
            int(start),
            int(end),
            0,
            "ade6",
        )

    def test_ade6_annotation_coordinates_are_unchanged(self) -> None:
        self.assertEqual(
            ADE6_COORDINATES,
            tuple(self.designer.annotationParser_.getCoordsFromName("ade6")),
        )

    def test_ade6_candidate_order_and_match_counts_are_unchanged(self) -> None:
        self.assertEqual(144, len(self.table))
        self.assertTrue(all(len(row) == 9 for row in self.table))

        expected_first_candidates = [
            ("ACATTGGCTTACGACGGTCG", [5, 1, 1, 1, 1, 1, 1]),
            ("GTGGCGACAGGGACACCTCG", [7, 1, 1, 1, 1, 1, 1]),
            ("AAACAGGTTGTAGGGATCCT", [10, 1, 1, 1, 1, 1, 1]),
            ("TGCCGCTCATTTGCCTGGTA", [12, 1, 1, 1, 1, 1, 1]),
            ("GCATAAGTAACCATGCGATC", [12, 1, 1, 1, 1, 1, 1]),
        ]
        actual_first_candidates = [
            (row[0], row[2:]) for row in self.table[:5]
        ]
        self.assertEqual(expected_first_candidates, actual_first_candidates)
        self.assertEqual(len(self.table), len(self.guide_matches))

    def test_first_ade6_candidate_legacy_primer_tuple_is_unchanged(self) -> None:
        self.assertEqual(
            (
                ADE6_GUIDE,
                "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
                "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
                (1316795, 1316797),
                1,
                "TGG",
            ),
            self.table[0][1],
        )

    def test_ade6_hr_sequences_are_unchanged(self) -> None:
        self.assertEqual([100, 100, 500], [len(value) for value in self.hr_dna])
        self.assertEqual(
            [
                "3efe8b6edeeca79ecd414b61bebcf8f08a8e9100779c35ec6144871f206ea51b",
                "024abcb465fdef8af25533b459abe0419cbc412d8b9ebc3eda4cc05a4174d89d",
                "0c355f0be1819366d73f0604dac12a49db005672b345bc496ca715a286922318",
            ],
            [
                hashlib.sha256(value.encode("ascii")).hexdigest()
                for value in self.hr_dna
            ],
        )

    def test_ade6_checking_primers_are_unchanged(self) -> None:
        self.assertEqual(2, len(self.checking_primers))
        first, second = self.checking_primers

        self.assertEqual("ACTGCGCACTAACTCACTACA", first["PRIMER_LEFT_0_SEQUENCE"])
        self.assertEqual("CGTCGCAGCACATTATTCGG", first["PRIMER_RIGHT_0_SEQUENCE"])
        self.assertEqual(255, first["PRIMER_PAIR_0_PRODUCT_SIZE"])
        self.assertEqual(1913, first["negative_result"])

        self.assertEqual("ACTGCGCACTAACTCACTACA", second["PRIMER_LEFT_1_SEQUENCE"])
        self.assertEqual("CAGCACATTATTCGGGGGGT", second["PRIMER_RIGHT_1_SEQUENCE"])
        self.assertEqual(250, second["PRIMER_PAIR_1_PRODUCT_SIZE"])
        self.assertEqual(1908, second["negative_result"])

    def test_legacy_runweb_result_shape_and_arguments_are_unchanged(self) -> None:
        legacy_result = (
            self.table,
            self.hr_dna,
            self.checking_primers,
            self.guide_matches,
        )

        with patch.object(
            self.designer,
            "run",
            return_value=legacy_result,
        ) as run:
            result = self.designer.runWeb(name="  ADE6  ", nMismatch=0)

        run.assert_called_once_with(
            "III",
            "1316337",
            "1317995",
            0,
            "ade6",
        )
        self.assertEqual(7, len(result))
        self.assertIs(self.table, result[0])
        self.assertIs(self.hr_dna, result[1])
        self.assertIs(self.checking_primers, result[2])
        self.assertEqual(("ade6", "III", "1316337", "1317995"), result[3:])

    def test_legacy_cache_filenames_are_unchanged(self) -> None:
        old_folder = self.designer.precomputed_folder
        try:
            with tempfile.TemporaryDirectory() as temporary_directory:
                self.designer.precomputed_folder = temporary_directory

                gene_path = self.designer._genPrecomputedName(
                    "ade6", 0, "III", "1316337", "1317995"
                )
                coordinate_path = self.designer._genPrecomputedName(
                    None, 2, "III", "100", "200"
                )

                self.assertEqual(
                    "SPCC1322.13_v3_n0.pickle",
                    os.path.basename(gene_path),
                )
                self.assertEqual(
                    "III_100_200_v3_n2.pickle",
                    os.path.basename(coordinate_path),
                )
        finally:
            self.designer.precomputed_folder = old_folder

    def test_legacy_oligo_cli_summary_and_coordinates_are_unchanged(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.designer.runCL(
                ["--oligo", ADE6_GUIDE, "--mismatch", "0"]
            )

        text = output.getvalue()
        self.assertIn(
            f"Querying S. pombe genome for oligo seed: {ADE6_GUIDE}",
            text,
        )
        self.assertRegex(text, re.compile(r"^\s*8\s+\|\s+5$", re.MULTILINE))
        self.assertRegex(text, re.compile(r"^\s*20\s+\|\s+1$", re.MULTILINE))
        self.assertIn("PAM coordinates: 1316795 - 1316797", text)
        self.assertIn("Cut: 1316791 | 1316792", text)
        self.assertIn("Strand: +", text)

    def test_legacy_oligo_cli_rejects_invalid_length(self) -> None:
        output = io.StringIO()
        with self.assertRaisesRegex(SystemExit, "1"):
            with contextlib.redirect_stdout(output):
                self.designer.runCL(["--oligo", "A" * 19])

        self.assertIn(
            "Oligo sequence must be 20 bp (seed only) or 23 bp",
            output.getvalue(),
        )


if __name__ == "__main__":
    unittest.main()
