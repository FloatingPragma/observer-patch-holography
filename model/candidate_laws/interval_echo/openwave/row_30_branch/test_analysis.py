"""Sensitivity checks for the diagnostic; the spectra here are mathematical fixtures."""

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_branch.classes import (
    calibrate_half_gap, spatial_operators,
)
from model.candidate_laws.interval_echo.openwave.row_30_branch.experiment import dispersion


def test_two_opposite_kg_branches_recover_one_gap_speed_and_held_curvature():
    k2 = np.array([0., .25, 1., 2.])
    basis = np.eye(len(k2))
    m, c = .2, .5
    shell = np.sqrt(m**2 + c**2*k2)
    panels = {}
    for name, energy in [("upper_222", 6.+shell), ("lower_122", 6.-shell)]:
        panel = dispersion(basis, (k2, basis), energy, basis)
        panel.update(reference_energy=energy[0], has_held_wavenumber=True)
        panels[name] = panel
    calibration = calibrate_half_gap(panels)
    np.testing.assert_allclose([calibration["rest_half_gap"], calibration["speed_squared"]], [m, c**2])
    for panel in panels.values():
        assert panel["passes_with_common_speed_and_half_gap"]
        assert panel["curvature_resolvable_at_enrolled_precision"]
        assert max(row["shift_error_over_forecast"] for row in panel["classes"][1:]) < 1e-12
    panels["upper_222"]["classes"][1]["shift_from_rest"] *= -1.
    assert not calibrate_half_gap(panels)["calibrated"]


def test_harmonic_reduction_measures_paths_through_grounded_out_vertices():
    # The two source sites are joined by a path of two unit resistors in series.
    graph = np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]])
    operators = spatial_operators(graph, [0,2])
    np.testing.assert_allclose(operators["induced"], np.zeros((2,2)))
    np.testing.assert_allclose(operators["harmonic"], [[.5,-.5],[-.5,.5]])
