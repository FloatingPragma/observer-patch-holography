from dataclasses import replace
from math import pi

import numpy as np
from scipy.optimize import brentq

from .measure import Source, expose, interior_volume, match_volume, run


def test_rate_only_changes_redshift_and_preserves_all_geometry_tapes():
    args = (0.2, 64, 2**-44)
    flat = expose(Source(0, 0), *args)
    rate = expose(Source(0.25, 0), *args)
    assert flat.central_ticks_per_emitter_unit > rate.central_ticks_per_emitter_unit
    np.testing.assert_array_equal(flat.radial_echoes, rate.radial_echoes)
    np.testing.assert_array_equal(flat.interior_laps, rate.interior_laps)
    assert flat.rim_lap == rate.rim_lap
    assert interior_volume(replace(rate, rim_lap=0), args[-1]) == interior_volume(rate, args[-1])


def test_clocked_rim_agrees_with_integrated_spatial_metric():
    ell = 0.35
    source = Source(1, 1)
    b = -1 / 8
    target = 4 * pi * ell**3 / 3

    def volume(r):
        return 4*pi*(r**3/3 + 3*b*r**5/5 + 3*b**2*r**7/7 + b**3*r**9/9)

    radius = brentq(lambda r: volume(r) - target, ell, 1.2*ell)
    continuum_rim = 4*pi*(radius*(1+b*radius**2))**2
    errors = []
    for shells in (32, 64, 128):
        measured, area, _ = match_volume(source, ell, shells, tick=2**-44)
        assert abs(measured-target) < 5e-12
        errors.append(abs(area-continuum_rim))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_ordinary_calibration_and_held_response_at_equal_volume():
    rows = run()
    fine = {row["arm"]: row for row in rows if row["ell"] == 0.05}
    kappa = fine["ordinary"]["inferred_coupling"]
    assert abs(kappa-1.5) < 0.001
    assert abs(fine["protected"]["inferred_coupling"] - kappa) < 0.001
    for arm in ("ordinary rate only", "protected rate only"):
        assert fine[arm]["rim_clock"] > 1
        assert abs(fine[arm]["inferred_coupling"]) < 1e-6
        assert abs(fine[arm]["spatial_component_per_density"] - 0.5) < 1e-6
    for arm in ("ordinary", "protected"):
        assert abs(fine[arm]["spatial_component_per_density"]) < 0.001
        errors = [abs(row["inferred_coupling"]-1.5) for row in rows if row["arm"] == arm]
        assert all(small < large/3.8 for large, small in zip(errors, errors[1:]))
