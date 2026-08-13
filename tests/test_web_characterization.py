import hashlib
import unittest
import urllib.parse
from io import BytesIO
from types import SimpleNamespace
from unittest.mock import Mock, patch

import webapp
from crispr4p.crispr4p import NGG


GUIDE = "ACATTGGCTTACGACGGTCG"
SPEDIT_FORWARD = (
    "CTAGAGGTCTCGGACTACATTGGCTTACGACGGTCGGTTTCGAGACCCTTCC"
)


def make_handler():
    """Construct a handler without opening a network socket."""
    return webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)


class FakeOligoDesigner:
    hit = NGG("III", 0, 1, GUIDE, "TGG")

    def __init__(self, *args, **kwargs):
        self.NGGs = None

    def getNGGsFromGenome(self):
        self.NGGs = {}

    def _single_table_worker(self, query, mismatches):
        return query, {
            8: [self.hit] * 5,
            10: [self.hit],
            12: [self.hit],
            14: [self.hit],
            16: [self.hit],
            18: [self.hit],
            20: [self.hit],
        }

    def getOligoHitCoordinates(self, hit):
        return (1316795, 1316797), (1316791, 1316792)


PRIMER_TUPLE = (
    GUIDE,
    "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
    "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
    (1316795, 1316797),
    1,
    "TGG",
)
DESIGN_TABLE = [[GUIDE, PRIMER_TUPLE, 5, 1, 1, 1, 1, 1, 1]]
HR_DNA = ("HR_FORWARD", "HR_REVERSE", "DELETED_DNA")
CHECKING_PRIMERS = [
    {
        "PRIMER_LEFT_0_SEQUENCE": "LEFT",
        "PRIMER_LEFT_0_TM": 59.4,
        "PRIMER_RIGHT_0_SEQUENCE": "RIGHT",
        "PRIMER_RIGHT_0_TM": 60.1,
        "PRIMER_PAIR_0_PRODUCT_SIZE": 255,
        "negative_result": 1913,
    }
]


class FakeDesignPrimerDesign:
    def __init__(self, *args, **kwargs):
        pass

    def runWeb(self, *args, **kwargs):
        return (
            DESIGN_TABLE,
            HR_DNA,
            CHECKING_PRIMERS,
            "ade6",
            "III",
            "1316337",
            "1317995",
        )


class TestWebRenderingCharacterization(unittest.TestCase):
    def test_oligo_result_html_is_unchanged(self) -> None:
        handler = make_handler()
        fake_crp = SimpleNamespace(
            PrimerDesign=FakeOligoDesigner,
            NGG=NGG,
        )

        with patch.object(webapp, "crp", fake_crp):
            result = handler.run_oligo_model(GUIDE, 0)

        self.assertEqual(
            "dd389db45eaca6b3babe396638609b0752b865f5638c14c32aa6c0e7e4c4082e",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )
        self.assertIn("<tr><td>8 bp</td><td>5</td></tr>", result)
        self.assertIn("<tr><td>20 bp</td><td>1</td></tr>", result)
        self.assertIn("1316795 - 1316797", result)
        self.assertIn("1316791 | 1316792", result)
        self.assertIn(SPEDIT_FORWARD, result)

    def test_design_result_html_is_unchanged(self) -> None:
        handler = make_handler()
        fake_crp = SimpleNamespace(PrimerDesign=FakeDesignPrimerDesign)

        with patch.object(webapp, "crp", fake_crp):
            result = handler.run_design_model("ade6", None, None, None)

        self.assertEqual(
            "93376b12610f2f10077b8c7b6bd80992a0d9cb65fdb93aed927026330e11ea3a",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )
        self.assertIn("<b>Name</b>=ade6", result)
        self.assertIn(GUIDE, result)
        self.assertIn("5'-HR_FORWARD-3'", result)
        self.assertIn("LEFT", result)
        self.assertIn("59 &deg;C", result)
        self.assertIn("255 (bp)", result)
        self.assertIn(SPEDIT_FORWARD, result)


class TestHttpPostCharacterization(unittest.TestCase):
    def make_post_handler(self, parameters):
        payload = urllib.parse.urlencode(parameters).encode("utf-8")
        handler = make_handler()
        handler.headers = {"Content-Length": str(len(payload))}
        handler.rfile = BytesIO(payload)
        handler.run_oligo_model = Mock(return_value="OLIGO RESULT")
        handler.run_design_model = Mock(return_value="DESIGN RESULT")
        handler.serve_form = Mock()
        return handler

    def test_oligo_takes_precedence_and_input_is_normalized(self) -> None:
        handler = self.make_post_handler(
            {
                "name": "ade6",
                "oligo_sequence": f"  {GUIDE.lower()}  ",
                "oligo_mismatch": "not-an-integer",
            }
        )

        handler.process_post()

        handler.run_oligo_model.assert_called_once_with(GUIDE, 0)
        handler.run_design_model.assert_not_called()
        handler.serve_form.assert_called_once_with("OLIGO RESULT")

    def test_gene_query_is_trimmed_and_forwarded(self) -> None:
        handler = self.make_post_handler({"name": "  ade6  "})

        handler.process_post()

        handler.run_design_model.assert_called_once_with(
            "ade6", None, None, None
        )
        handler.run_oligo_model.assert_not_called()
        handler.serve_form.assert_called_once_with("DESIGN RESULT")

    def test_coordinate_query_is_trimmed_and_forwarded(self) -> None:
        handler = self.make_post_handler(
            {
                "chromosome": " III ",
                "coor_lower": " 1316337 ",
                "coor_upper": " 1317995 ",
            }
        )

        handler.process_post()

        handler.run_design_model.assert_called_once_with(
            None, "III", "1316337", "1317995"
        )
        handler.run_oligo_model.assert_not_called()
        handler.serve_form.assert_called_once_with("DESIGN RESULT")

    def test_missing_query_keeps_existing_error_message(self) -> None:
        handler = self.make_post_handler({})

        handler.process_post()

        handler.run_oligo_model.assert_not_called()
        handler.run_design_model.assert_not_called()
        handler.serve_form.assert_called_once_with(
            '<font color="red"><h3>Error: Please fill either Name, '
            'Coordinates, or Oligo Sequence</h3></font>'
        )


if __name__ == "__main__":
    unittest.main()
