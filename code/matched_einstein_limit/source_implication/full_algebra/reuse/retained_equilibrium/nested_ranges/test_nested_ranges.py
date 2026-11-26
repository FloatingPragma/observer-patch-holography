import json

import numpy as np

from . import certificate, kernels, records


def test_masks_come_from_the_same_saved_return_banks():
    raw = records.read()
    assert raw == json.loads((kernels.HERE / "records.json").read_text())
    c = raw["intervals"]
    assert 0 < c["2"]["area"][0] < c["4"]["area"][0]
    assert 0 < c["2"]["volume"][0] < c["4"]["volume"][0]
    assert c["6"]["area"][1] < 0 and c["6"]["volume"][0] < 0 < c["6"]["volume"][1]


def test_exact_late_fiber_and_nested_observables_are_regular():
    saved = json.loads((kernels.HERE / "certificate.json").read_text())
    report = certificate.evaluate()
    assert report["late_range_rank"] == 12 and report["fiber_dimension"] == 6
    assert saved["early_geometry_rank"] == report["early_geometry_rank"] == 4
    assert (
        max(
            saved["right_inverse_residual_infinity_upper"],
            report["right_inverse_residual_infinity_upper"],
        )
        < 0.00374
    )


def test_deficits_have_the_correct_order_and_positivity():
    with np.load(kernels.HERE / "deficits.npz") as saved:
        numerators = saved["numerators"]
    deficit = (numerators[:, :, :, 0] + 1j * numerators[:, :, :, 1]) / kernels.GRID
    np.testing.assert_allclose(deficit, deficit.conj().swapaxes(-1, -2), atol=1e-15)
    # Every successful early return has a smaller score than a timeout.
    assert np.linalg.eigvalsh(deficit).max() < 1e-13
    effects = (
        deficit + np.asarray(kernels.CAPS)[:, None, None, None, None] * np.eye(4) / 2
    )
    assert np.linalg.eigvalsh(np.diff(effects, axis=0)).min() > 0
