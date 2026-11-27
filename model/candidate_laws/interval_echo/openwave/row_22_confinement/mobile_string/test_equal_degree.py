"""Port-level inversion and held Ramsey records on the re-paired strand."""
import numpy as np
import pytest

from . import commission, equal_degree as control, material as mat, measure


@pytest.mark.parametrize("level", range(3))
def test_material_preparation_and_operator_inversion(level):
    board = mat.federation(level)
    sites = control.ports(level)
    dimers = np.asarray(sites).reshape(16, 2)
    np.testing.assert_array_equal(board.degree[dimers], [[5, 5], [6, 6]] * 8)
    stock = commission.material_stock(level, 16, sites)
    with np.load(control.RECEIVERS, allow_pickle=False) as bank:
        assert bank[f"stock_sha_{level}"].item() == stock["ordered_address_sha256"]
        np.testing.assert_array_equal(bank[f"ports_{level}"], sites)
    for kappa in (8, 0):
        calibrated = mat.calibration_energies(kappa)
        np.testing.assert_allclose(calibrated, [4 * kappa, 1, kappa], rtol=0, atol=1e-12)
        for length in mat.LENGTHS:
            lower, upper = (control.pinned_matrix(board, sites, length, i, calibrated) for i in (0, 1))
            # Bitwise complement on every unfixed physical orientation reverses basis order.
            np.testing.assert_allclose(upper - lower[::-1, ::-1], 4 * kappa * np.eye(len(upper)),
                                       rtol=0, atol=1e-10)
    # Restore the original route only. Its physical-port diagonal must recover
    # the measured reference field and spoil the constant-work conjugacy.
    original = mat.route(level)
    cal = mat.calibration_energies()
    low, high = (control.pinned_matrix(board, original, 4, pin, cal) for pin in (0, 1))
    for pin, physical in enumerate((low, high)):
        reduced = mat.core(4, pin, mat.G)
        if pin == 0:
            offset = physical[0, 0] - reduced[0, 0]
        np.testing.assert_allclose(physical - reduced, offset * np.eye(len(reduced)), atol=1e-10)
    assert np.ptp(np.diag(high - low[::-1, ::-1])) > 1


def test_held_quadratures_and_clock():
    for (level, kappa), (cal, readings) in control.decode().items():
        np.testing.assert_allclose(cal, [4 * kappa, 1, kappa], rtol=0, atol=0.1)
        assert set(readings) == set(mat.LENGTHS)
        for length, pair in readings.items():
            word = pair["sharp_intact"]
            if kappa:
                assert abs(word["line"] - 4 * kappa) < 0.002, (level, kappa, length, word)
                assert word["min_visibility"] > 0.975
                assert word["status"] == "resolved"
            elif word["line"] is None:
                # At zero contact the free spins can lose their return
                # visibility. Such an unwrap cannot select a winding.
                assert word["status"] == "unresolved"
                assert word["low_visibility_holds"]
            else:
                assert abs(word["line"]) < 0.002
    with np.load(control.RECEIVERS, allow_pickle=False) as bank:
        for key in bank.files:
            if key.startswith("static_"):
                counts = bank[key]
                assert counts.shape == (2, 2, len(mat.HOLDS))
                assert np.issubdtype(counts.dtype, np.integer)
                assert np.all(counts <= measure.SHOTS)
                _, _, strength, length = key.split("_")
                if strength == "0":
                    # Independent spins, both pinned interiors identical:
                    # their two-copy return is cos(t/2)^(2m), real for
                    # every hold, including low-visibility winding misses.
                    real = np.cos(mat.HOLDS / 2) ** (2 * (int(length) - 2))
                    p = np.stack(((1 + real) / 2, np.full_like(real, 0.5)))
                    sigma = np.sqrt(measure.SHOTS * p * (1 - p))
                    assert np.all(abs(counts - measure.SHOTS * p) <= 8 * np.maximum(1, sigma))
