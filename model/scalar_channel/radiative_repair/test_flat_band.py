import unittest

from .flat_band import certificates


class FlatBandCertificateTests(unittest.TestCase):
    def test_federated_flat_band_and_controls(self):
        rows = certificates()["federations"]
        expected = [
            (20, 30, 375, 31, 0),
            (80, 120, 6900, 121, 24),
            (320, 480, 114000, 481, 104),
        ]
        actual = [
            (
                row["cells"],
                row["lowest_band_dimension"],
                row["symmetric_dark_lower_bound_seam_constraints"],
                row["disabled_seam_control_lowest_band_dimension"],
                row["full_rank_seam_control_dimension_at_original_gap"],
            )
            for row in rows
        ]
        self.assertEqual(actual, expected)
