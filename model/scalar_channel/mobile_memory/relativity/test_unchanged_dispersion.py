from fractions import Fraction

import numpy as np
from scipy.linalg import expm, sqrtm

from ....contrast_memory_spherical_v1 import carrier
from .unchanged_dispersion import WAIT, _federation_record, run


def test_exact_mass_shell_and_receiver_separation_on_all_carriers():
    result = run()
    assert Fraction(result["charged_kinetic_compression_factor"]) == Fraction(1, 96)
    assert Fraction(result["charged_squared_speed_upper_bound"]) == Fraction(11, 512)
    np.testing.assert_allclose(result["cage_squared_speed_lower_bound"], 4 * np.sqrt(6))
    assert [row["cells"] for row in result["federations"]] == [20, 80, 320]
    for row in result["federations"]:
        assert row["receiver_laplacian_element"] == 0
        assert row["receiver_native_squared_shell_defect"] == row["two_step_paths"] > 0
        upper = Fraction(row["receiver_native_probability_upper_bound"])
        lower = Fraction(row["receiver_spectral_audit_probability_lower_bound"])
        gap = Fraction(row["receiver_spectral_audit_excess_lower_bound"])
        assert lower - upper == gap > 0
        receiver = row["receiver_cell"]
        native = row["native_cell_probabilities"]
        audit = row["slope_matched_spectral_audit_cell_probabilities"]
        parked = row["parked_cell_probabilities"]
        for records in (native, audit, parked):
            assert len(records) == row["cells"]
            np.testing.assert_allclose(sum(records), 1, rtol=0, atol=1e-12)
        assert audit[receiver] - native[receiver] > float(gap)
        assert parked[receiver] == 0


def test_cell_detector_records_against_dense_propagators():
    row = _federation_record(0)
    federation = carrier.federation(0)
    laplacian = np.zeros((federation.cells, federation.cells))
    for a, b in federation.external:
        difference = np.zeros(federation.cells)
        difference[a // 12] = 1
        difference[b // 12] = -1
        laplacian += np.outer(difference, difference)
    audit = sqrtm(24 * np.eye(federation.cells) + 4 * np.sqrt(6) * laplacian)
    wait = float(WAIT)
    for generator, key in (
        (laplacian, "native_cell_probabilities"),
        (audit, "slope_matched_spectral_audit_cell_probabilities"),
    ):
        amplitudes = expm(-1j * wait * generator)[:, row["sender_cell"]]
        np.testing.assert_allclose(abs(amplitudes) ** 2, row[key], atol=1e-12, rtol=0)
