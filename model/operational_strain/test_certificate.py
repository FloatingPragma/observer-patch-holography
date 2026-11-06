"""Independent work identification, complete native effects and hostile records."""
from dataclasses import replace
from fractions import Fraction as F
import json
import subprocess
import sys

import numpy as np
import pytest

from . import certificate as c


@pytest.mark.parametrize("x", [F(1), F(2), F(1, 1000), F(17, 13), F(10000)])
def test_log_enclosures_against_independent_high_precision_series(x):
    import mpmath as mp
    with mp.workdps(100):
        exact = mp.log(mp.mpf(x.numerator)/x.denominator)
        interval = c.log_bounds(x)
        convert = lambda a: mp.mpf(a.numerator)/a.denominator
        assert convert(interval.lo) <= exact <= convert(interval.hi)
        assert interval.hi-interval.lo < F(1, 10**20)


@pytest.mark.parametrize("amplitude,offset,t,alpha", [
    (F(72), F(72), F(1, 4), F(1, 3)),
    (F(432, 7), F(576, 7), F(1, 8), F(1, 2)),
    (F(1, 10000), F(3), F(1, 2), F(1, 6)),
])
def test_three_settings_recover_unknown_scale_and_zero(amplitude, offset, t, alpha):
    rows = c.supplied_power_records(amplitude, offset, t)
    result = c.identify_power_response(rows)
    assert result["contraction"] == c.point(t)
    assert result["amplitude"] == c.point(amplitude)
    assert result["offset"] == c.point(offset)
    assert result["alpha"].contains(alpha)
    assert result["pressure"].contains(amplitude*alpha/3)
    assert result["reference_ratio"].contains(amplitude*alpha/(amplitude+offset))


def test_every_corner_of_bounded_record_error_contains_the_true_pressure():
    from itertools import product
    epsilon = F(1, 10000)
    exact = c.supplied_power_records(72, 72, F(1, 4))
    for signs in product((-1, 1), repeat=3):
        records = tuple(replace(row, energy=row.energy+sign*epsilon, error=epsilon)
                        for row, sign in zip(exact, signs, strict=True))
        result = c.identify_power_response(records)
        assert result["alpha"].contains(F(1, 3))
        assert result["amplitude"].contains(72)
        assert result["offset"].contains(72)
        assert result["pressure"].contains(8)
        assert result["reference_ratio"].contains(F(1, 6))


def test_two_settings_leave_distinct_positive_work_extensions():
    first = c.two_setting_family(144, 90, F(1, 4))
    second = c.two_setting_family(144, 90, F(1, 8))
    r1 = c.supplied_power_records(*first, F(1, 4))
    r2 = c.supplied_power_records(*second, F(1, 8))
    assert r1[:2] == r2[:2]
    assert r1[2].energy == F(153, 2) != r2[2].energy == F(333, 4)
    p1 = c.identify_power_response(r1)["pressure"]
    p2 = c.identify_power_response(r2)["pressure"]
    assert p1.contains(8) and p2.contains(F(72, 7)) and p1.hi < p2.lo


def test_constant_energy_identifies_zero_work_only_under_the_positive_amplitude_family():
    records = c.supplied_power_records(7, 13, 1)
    result = c.identify_power_response(records)
    assert result["pressure"] == c.point(0)
    assert result["alpha"] == c.point(0)
    assert result["amplitude"] is result["offset"] is None
    uncertain = tuple(replace(row, error=F(1, 1000)) for row in records)
    assert c.identify_power_response(uncertain) == {"status": "unresolved_energy_differences"}
    touching = tuple(replace(row, energy=row.energy+sign, error=1)
                     for row, sign in zip(records, (-1, 1, -1), strict=True))
    assert c.identify_power_response(touching) == {"status": "unresolved_energy_differences"}
    assert c.identify_power_response(c.supplied_power_records(7, -7, 1))["reference_ratio"] is None


def test_energy_offset_changes_ratio_without_changing_pressure():
    a = c.identify_power_response(c.supplied_power_records(72, 0, F(1, 4)))
    b = c.identify_power_response(c.supplied_power_records(72, 72, F(1, 4)))
    assert a["pressure"] == b["pressure"]
    assert a["reference_ratio"].contains(F(1, 3))
    assert b["reference_ratio"].contains(F(1, 6))


def test_general_stencil_propagates_meter_error_and_regular_remainder():
    # H(V)=V^3: the central derivative error is exactly h^2, M=6.
    center, h, epsilon = F(3), F(1, 10), F(1, 10000)
    for left_sign in (-1, 1):
        for right_sign in (-1, 1):
            rows = (c.EnergyRecord(center-h, (center-h)**3+left_sign*epsilon, epsilon),
                    c.EnergyRecord(center+h, (center+h)**3+right_sign*epsilon, epsilon))
            result = c.pressure_stencil(rows, derivative_order=3, derivative_bound=6)
            assert result["pressure"].contains(-27)
            assert result["regularity_radius"] == h*h
            assert result["record_error_radius"] == epsilon/h
            shifted = tuple(replace(row, energy=row.energy+F(10000, 7)) for row in rows)
            assert c.pressure_stencil(shifted, derivative_order=3, derivative_bound=6) == result
    # H(V)=V^2 has a certified second derivative bound2 on the full interval.
    rows = (c.EnergyRecord(center-h, (center-h)**2), c.EnergyRecord(center+h, (center+h)**2))
    interval = c.pressure_stencil(rows, derivative_order=2, derivative_bound=2)["pressure"]
    assert interval == c.Interval(-6-h, -6+h)


def test_finite_samples_do_not_bound_an_unrestricted_smooth_positive_derivative():
    nodes = (F(2), F(3), F(4), F(192), F(2))
    coefficients = c.invisible_strain_polynomial(3, nodes)
    evaluate = lambda x: sum((a*x**k for k, a in enumerate(coefficients)), F())
    slope = sum((k*a*F(3)**(k-1) for k, a in enumerate(coefficients) if k), F())
    assert all(evaluate(x) == 0 for x in nodes)
    assert evaluate(F(3)) == 0 and slope == 1
    # At every sample B exp(d u/B)=B, but its reference derivative is d.
    assert all(F(7)*d/F(7)*slope == d for d in (-1000, 0, 1000))


def test_complete_native_census_effects_on_coherent_reference_input():
    report = c.native_census_certificate()
    assert report["raw_input_branches"] == 264
    assert report["effect_matrix_units_checked"] == 13068
    rng = np.random.default_rng(1607)
    amplitudes = rng.normal(size=132)+1j*rng.normal(size=132)
    amplitudes /= np.linalg.norm(amplitudes)
    state = amplitudes.reshape(66, 2)
    received = np.zeros(3)
    index = {x: i for i, x in enumerate(c.census.native.states(binary=True))}
    for row in report["all_tapes"]:
        output = np.zeros((66, 2), complex)
        for source, endpoint in zip(row["input_indices"], row["endpoints"], strict=True):
            output[index[endpoint]] += state[source]/2
        received[row["decoded_initial_charge"]] += np.linalg.norm(output)**2
    expected = np.zeros(3)
    for source, occupation in enumerate(c.census.native.states(binary=True)):
        expected[sum(occupation[:2])] += np.linalg.norm(state[source])**2
    np.testing.assert_allclose(received, expected, atol=8e-16)
    assert abs(sum(received)-1) < 1e-15


def test_weighted_force_interval_has_the_native_expectation():
    p = (F(17, 25), F(23, 75), F(1, 75))
    result = c.sector_work_interval(p, (c.Interval(1, 2), c.Interval(-3, -2), c.point(7)))
    assert result.lo == sum(x*y for x, y in zip(p, (1, -3, 7), strict=True))
    assert result.hi == sum(x*y for x, y in zip(p, (2, -2, 7), strict=True))


@pytest.mark.parametrize("fault", ["state", "order", "volume", "negative_error", "boolean", "float"])
def test_invalid_record_contracts_are_rejected(fault):
    records = c.supplied_power_records(72, 72, F(1, 4))
    with pytest.raises(ValueError):
        if fault == "state":
            c.identify_power_response((records[0], replace(records[1], state_tag="different"), records[2]))
        elif fault == "order":
            c.identify_power_response(tuple(reversed(records)))
        elif fault == "volume":
            c.identify_power_response((records[0], records[1], replace(records[2], volume=12289)))
        elif fault == "negative_error":
            replace(records[0], error=-1)
        elif fault == "boolean":
            replace(records[0], energy=True)
        else:
            replace(records[0], energy=144.0)


def test_conflicting_or_unresolved_responses_are_not_decoded_as_masses_or_pressure():
    with pytest.raises(ValueError):
        c.identify_power_response(tuple(c.EnergyRecord(3*64**j, j) for j in range(3)))
    with pytest.raises(ValueError):
        c.identify_power_response(tuple(c.EnergyRecord(3*64**j, 4-j) for j in range(3)))
    uncertain = c.supplied_power_records(1, 1, F(1, 4), error=1)
    assert "pressure" not in c.identify_power_response(uncertain)


def test_receipt_replays_all_exact_claims_and_preserves_scope():
    result = c.run()
    assert json.loads(c.OUTPUT.read_text()) == result
    assert result["native_scalar_source"]["exact_events_replayed"] == 240
    assert result["scope"]["finite_family_derivative_identified_from_added_records"] is True
    assert result["scope"]["energy_meter_or_work_coupling_constructed_from_native_repair"] is False
    assert result["scope"]["physical_volume_and_strain_intervention_source_selected"] is False
    assert not any(path.startswith("docs/") for path in result["source_sha256"])


def test_runtime_dependencies_are_bound_including_eager_package_initializers():
    script = """
from model.operational_strain import certificate as c
from pathlib import Path
import json, sys
c.run()
pins = c.source_bindings()
loaded = {Path(m.__file__).resolve().relative_to(c.ROOT).as_posix()
          for m in tuple(sys.modules.values()) if getattr(m, '__file__', None)
          and Path(m.__file__).is_file()
          and Path(m.__file__).resolve().is_relative_to(c.ROOT)
          and '.venv/' not in Path(m.__file__).resolve().relative_to(c.ROOT).as_posix()}
print(json.dumps(sorted(loaded-pins.keys())))
"""
    result = subprocess.run([sys.executable, "-c", script], cwd=c.ROOT,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == []
