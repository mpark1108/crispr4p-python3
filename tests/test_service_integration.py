import hashlib
import unittest

from crispr4p.service import Crispr4pService
from crispr4p.web_views import render_oligo_result


GUIDE = "ACATTGGCTTACGACGGTCG"


class TestServiceIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        service = Crispr4pService.from_project_data(
            precomputed_folder="precomputed",
        )
        cls.result = service.analyze_oligo(GUIDE, n_mismatch=0)
        cls.design_result = service.design_gene("ade6", n_mismatch=0)

    def test_real_ade6_oligo_analysis_is_unchanged(self) -> None:
        self.assertEqual(
            {8: 5, 10: 1, 12: 1, 14: 1, 16: 1, 18: 1, 20: 1},
            dict(self.result.match_counts),
        )
        self.assertEqual(1, len(self.result.full_matches))
        match = self.result.full_matches[0]
        self.assertEqual("III", match.chromosome)
        self.assertEqual((1316795, 1316797), match.pam_coordinates)
        self.assertEqual((1316791, 1316792), match.cut_coordinates)
        self.assertEqual(1, match.strand)
        self.assertEqual(GUIDE, match.seed)
        self.assertEqual("TGG", match.pam)
        self.assertEqual(
            "dd389db45eaca6b3babe396638609b0752b865f5638c14c32aa6c0e7e4c4082e",
            hashlib.sha256(
                render_oligo_result(self.result).encode("utf-8")
            ).hexdigest(),
        )

    def test_real_ade6_design_service_result_is_unchanged(self) -> None:
        self.assertEqual("ade6", self.design_result.name)
        self.assertEqual("III", self.design_result.chromosome)
        self.assertEqual("1316337", self.design_result.start)
        self.assertEqual("1317995", self.design_result.end)
        self.assertEqual(144, len(self.design_result.guide_table))
        self.assertEqual(GUIDE, self.design_result.guide_table[0][0])
        self.assertEqual(
            [5, 1, 1, 1, 1, 1, 1],
            self.design_result.guide_table[0][2:],
        )


if __name__ == "__main__":
    unittest.main()
