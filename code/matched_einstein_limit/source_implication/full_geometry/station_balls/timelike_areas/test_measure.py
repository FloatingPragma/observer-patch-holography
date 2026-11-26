import json
import math

import numpy as np

from . import measure as m


def test_pullback_radial_coordinate_and_independent_surface():
    tilted = m.Slice(5/3, m.S, 10)
    ell = 1/128
    metric = tilted.metric(ell)
    assert np.min(np.linalg.eigvalsh(metric)) > 0.97
    np.testing.assert_allclose((metric @ tilted.n[:, :, None])[:, :, 0], tilted.n,
                               atol=1.5e-15, rtol=0)
    assert math.isclose(tilted.area(ell), (tilted.angular @ np.sqrt(np.linalg.det(metric))) * ell**2,
                        abs_tol=2e-17)
    delta = 1e-7
    assert math.isclose((tilted.volume(ell+delta)-tilted.volume(ell-delta))/(2*delta),
                        tilted.area(ell), rel_tol=2e-9)


def test_quadrature_resolution_and_matched_proper_volume():
    ref = json.loads((m.HERE / "results.json").read_text())
    repeat = m.measure(order=14)
    for name in ("rest", "axial"):
        for first, second in zip(ref["tables"][name], repeat["tables"][name], strict=True):
            target = first["target_volume"]
            for sign in ("plus", "minus"):
                a, b = first[sign], second[sign]
                assert math.isclose(a["volume"], target, rel_tol=2e-13)
                assert math.isclose(b["volume"], target, rel_tol=2e-13)
                assert math.isclose(a["boundary_area"], b["boundary_area"], rel_tol=2e-13)
            if name == "rest":
                assert math.isclose(first["plus"]["boundary_area"], first["minus"]["boundary_area"], rel_tol=2e-13)
            assert abs(first["signed_half_difference_per_s_ell4"]-second["signed_half_difference_per_s_ell4"]) < 1e-7


def test_each_curvature_and_signed_limit():
    saved = json.loads((m.HERE / "results.json").read_text())
    limit = -64*math.pi/45
    for name, gamma in (("rest", 1), ("axial", 5/3)):
        for label, s in (("plus", m.S), ("minus", -m.S)):
            guu = (1+2*gamma**2)*(m.ROUND_RADIUS**-2-3*m.B) + 3*(gamma**2-1)*s
            errors = []
            for row in saved["tables"][name]:
                ell = row["ell"]
                estimate = (4*math.pi*ell**2-row[label]["boundary_area"])/(4*math.pi*ell**4/15)
                errors.append(abs(estimate-guu))
            assert all(a > 3.9*b for a, b in zip(errors, errors[1:]))
            assert errors[-1] < 0.071
    differences = [abs(row["signed_half_difference_per_s_ell4"]-limit) for row in saved["tables"]["axial"]]
    assert all(a > 3.9*b for a, b in zip(differences, differences[1:]))
    assert differences[-1] < 0.0012
