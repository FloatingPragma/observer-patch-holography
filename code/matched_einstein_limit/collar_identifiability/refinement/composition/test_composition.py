from itertools import product

import numpy as np
import pytest

from . import measure


@pytest.fixture(scope="module")
def results():
    return measure.run()


def test_raw_edge_law_and_product_restrictions(results):
    assert [r["boundary"] for r in results["local_arrays"]] == [3, 6, 12]
    assert all(r["charge_per_u_squared"] == pytest.approx(r["leading"], rel=1e-11)
               for r in results["local_arrays"])
    fresh, formed = .28, .73
    left, right = np.array([.4, .12, .6]), np.array([.3, .7, .6])
    signs = np.array(list(product((-1, 1), repeat=3)))
    y = (formed-fresh)*(left-right)
    s = .5
    # Weakening derivative from the full joint of three separately identified return bits.
    def joint(pairing):
        p = measure.edge_probabilities(left, right, pairing, fresh, formed)
        return np.prod(np.where(signs == 1, p, 1-p), axis=1)

    tangent = (joint(s-1e-6)-joint(s+1e-6))/2e-6
    assert -tangent @ np.log(joint(s)) == pytest.approx(y @ np.arctanh(s*y), abs=3e-10)
    for column in range(3):
        assert joint(s)[signs[:, column] == 1].sum() == pytest.approx((1+s*y[column])/2)


def test_observer_ball_smearing(results):
    for row in results["ball_smearing"]:
        assert row["pooled_square"] == pytest.approx(1)
        assert row["scaled_difference"] == pytest.approx(1/28)


def test_closed_transport_and_repaired_current(results):
    first = results["transport"][0]
    assert first["squared_sum"] == 5/32
    assert first["new_unresolved_variance"] == 27/32
    # Check an independently evaluated antisymmetric current on a non-spherical signed datum.
    chart = np.random.default_rng(17).normal(size=(5, 6, 7))
    outgoing = np.zeros_like(chart)
    for axis in range(3):
        edge = (chart-np.roll(chart, -1, axis=axis))/8
        outgoing += edge-np.roll(edge, 1, axis=axis)
    np.testing.assert_allclose(measure.periodic_step(chart)-chart, -outgoing, atol=2e-15)
    assert measure.predecessor_variance(chart).sum() > 0
    assert all(row["full_second_moment"] == 1 for row in results["transport"])


def test_phase_population_and_visibility(results):
    for row in results["phase_readout"]["two_rail"]:
        assert row["response_per_population"] == pytest.approx(row["population_slope_limit"], abs=2e-13)
        np.testing.assert_allclose(row["half_visibility_quadratures"], np.array(row["quadratures"])/2)
    for row in results["phase_readout"]["three_rail"]:
        assert row["response_per_population"] == pytest.approx(row["population_slope_limit"], abs=2e-13)
        assert row["response_per_population"] > 0
        assert row["zero_visibility_max_bias"] < 1e-14
    # All failure and empty-rail events are retained. A purely imaginary signal has real part zero.
    np.testing.assert_allclose(measure.measured_quadratures(.4j, .5, 1), [0., .4], atol=2e-14)
    tangent = results["phase_readout"]["matched_tangents"]
    assert tangent["same_local_occupation_derivative"] == pytest.approx(.25)
    assert tangent["native_shape_slope"] == pytest.approx(tangent["shape_slope_limit"], abs=2e-13)
    assert tangent["population_slope"] > 0 > tangent["native_shape_slope"]


def test_weighted_number_meter(results):
    data = results["weighted_occupancy_read"]
    p, q = np.asarray(data["formed"]), np.asarray(data["independent"])
    rows = data["gains"]
    assert [r["gain"] for r in (rows[0], rows[-1])] == pytest.approx([.064542391355, .058704049578])
    assert all(r["slope"] < 0 for r in rows)
    for row in rows:
        h = row["hit"]
        assert row["gain"] == pytest.approx(-(q-p) @ np.log((1-h)*p+h*q))
    assert data["flat"][-1]["normalized_charge"] == pytest.approx(1, abs=3e-9)
    # Turning off the material read-stage contrast leaves the same position lottery with no usable gain.
    control = measure.admission.material_panel(1, 1, 0)["records"][-1]
    assert np.max(abs(np.asarray(control["d"]))) < 1e-13


def test_independent_address_and_hopping_reads_on_held_states(results):
    panel = results["weighted_energy_read"]
    ground, staggered, wave = panel["rows"]
    np.testing.assert_allclose([r["weighted_occupation"] for r in panel["rows"]], ground["weighted_occupation"],
                               rtol=1e-14)
    assert ground["weighted_energy"] == pytest.approx(panel["eigenenergy"]*ground["weighted_occupation"], abs=1e-16)
    assert staggered["weighted_energy"] > wave["weighted_energy"] > ground["weighted_energy"]
    assert panel["largest_lottery_offset"] < .5
    for row in panel["rows"]:
        assert abs(row["weighted_energy"]-row["observed_signed_estimator"]) < 1e-16
        assert row["terminal_hit"] == pytest.approx(.5+panel["gamma"]*row["weighted_energy"], abs=1e-15)


def test_duplicating_the_detector_and_copying_its_record(results):
    data = results["weighted_occupancy_read"]
    p, q = np.asarray(data["formed"]), np.asarray(data["independent"])
    delta = q-p
    volume = .35
    for switch in (.01, .001):
        fresh_returns = np.outer(p+switch*delta, p+switch*delta)
        base = np.outer(p, p)
        full = base+volume*(fresh_returns-base)
        first = p + volume*switch*delta
        one = -switch*delta @ np.log(first)
        two = -np.sum((fresh_returns-base)*np.log(full))
        copied = np.diag(first)
        copy_slope = -switch*delta @ np.log(copied.diagonal())
        assert copy_slope == pytest.approx(one)
        assert two/one == pytest.approx(2, abs=2e-7)
