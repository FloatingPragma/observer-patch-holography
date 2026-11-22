import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from . import gateway_acquisition as acquisition


@pytest.fixture(scope="module")
def preparation():
    return acquisition.prepare_layout((1, 3, 10))


def test_every_first_click_and_native_clear_column(preparation):
    # The isolated cell in candidate.generator has one degree-reference term per port.
    # Dropping just these references keeps the full scalar phase of every native contact.
    h = candidate.generator(cell()).toarray() - candidate.KAPPA * np.eye(12)
    pulse = expm(-1j * h)
    tips = list(preparation.gateways)
    vacancy = np.eye(12)
    vacancy[tips, tips] = 0
    clear = vacancy @ pulse
    assert preparation.first_click.shape == (1024, 3, 12)
    for i in range(1024):
        before_read = pulse @ np.linalg.matrix_power(clear, i)
        np.testing.assert_allclose(preparation.first_click[i], abs(before_read[tips]) ** 2,
                                   atol=2e-12, rtol=0)
    np.testing.assert_allclose(preparation.accepted_vectors, np.linalg.matrix_power(clear, 1024),
                               atol=2e-11, rtol=0)
    np.testing.assert_allclose(preparation.success_by_launch + preparation.first_click.sum(axis=(0, 1)),
                               np.ones(12), atol=2e-12, rtol=0)


def test_ordered_attempts_full_histories_and_exhaustion(preparation):
    p = preparation
    rho = np.zeros((12, 12), complex)
    reach = history = 1.0
    clicks = []
    for i, launch in enumerate(acquisition.LAUNCHES):
        conditional = p.first_click[:, :, launch]
        failed, accepted = p.attempt_joint(i + 1)
        np.testing.assert_allclose(failed, reach * conditional, atol=1e-14, rtol=1e-12)
        assert accepted == pytest.approx(reach * p.success_by_launch[launch], abs=1e-14)
        assert p.history_probability(clicks, accepted=True) == pytest.approx(
            history * p.success_by_launch[launch], abs=0, rel=1e-12)
        rho += reach * np.outer(p.accepted_vectors[:, launch], p.accepted_vectors[:, launch].conj())
        r, g = np.unravel_index(np.argmax(conditional), conditional.shape)
        clicks.append((r + 1, p.gateways[g]))
        history *= conditional[r, g]
        reach *= conditional.sum()
    assert len(clicks) == 48 and 0 < history < reach
    assert p.history_probability(clicks, accepted=False) == pytest.approx(history, abs=0)
    assert p.exhaustion == pytest.approx(reach, abs=1e-16)
    assert p.exhaustion + p.attempt_accepted.sum() == pytest.approx(1, abs=2e-12)
    np.testing.assert_allclose(rho, p.accepted_density, atol=2e-13, rtol=0)


def test_physical_wiring_and_recycled_load_paths(preparation):
    board = federation(0)
    home = 13
    wires = acquisition.gateway_wiring(board, home)
    assert tuple(w[0] for w in wires) == preparation.gateways
    for port, external, sign in wires:
        home_port = 12 * home + port
        assert external // 12 != home
        assert ((home_port, external) if sign == 1 else (external, home_port)) in board.external
        for target in range(12):
            path = preparation.routes[port, target]
            assert path[0] == port and path[-1] == target and len(path) <= 4
            assert all(12 * home + b in board.neighbours[12 * home + a] for a, b in zip(path, path[1:]))
    with pytest.raises(ValueError):
        acquisition.scored_appointments(acquisition.prepare_layout((0, 2, 8)), level=0, home=13)


def test_exchange_off_and_one_sharp_load_after_failure():
    p = acquisition.prepare_layout((0, 2, 8), pulse=0)
    expected = np.zeros((1024, 3, 12))
    for i, port in enumerate(p.gateways):
        expected[0, i, port] = 1
    np.testing.assert_array_equal(p.first_click, expected)
    np.testing.assert_array_equal(p.attempt_reach, [1, 1] + [0] * 46)
    np.testing.assert_array_equal(p.attempt_accepted, [0, 1] + [0] * 46)
    assert p.exhaustion == 0
    assert p.history_probability([(1, 0)], accepted=True) == 1
    assert p.history_probability([(2, 0)], accepted=True) == 0
    np.testing.assert_allclose(p.accepted_density, np.diag(np.eye(12)[1]), atol=1e-14, rtol=0)


def test_record_timing_and_free_flight_controls(preparation, monkeypatch):
    # Acquisition must use only signed gateway queries. A port read resolves interior phase.
    with monkeypatch.context() as m:
        def premature_read():
            raise AssertionError("read all home ports before release")
        m.setattr(acquisition, "terminal_port_effects", premature_read)
        p = acquisition.prepare_layout(preparation.gateways)
    diagonal = np.diag(np.diag(p.accepted_density))
    assert np.linalg.norm(diagonal - p.accepted_density) > 0.1
    rows = acquisition.scored_appointments(p)["readings"]
    assert [row["t"] for row in rows] == [0, 1, 7, 37]
    for row in rows:
        for arm in ("coherent", "dephased", "unfiltered_uniform", "stochastic_swap"):
            record = row[arm]
            assert len(record["ports"]) == 12 and min(record["ports"]) >= 0
            assert record["source_exhaustion"] == (0 if arm == "unfiltered_uniform" else p.exhaustion)
            assert sum(record["ports"]) + record["exterior_quiet"] + record["source_exhaustion"] == pytest.approx(1, abs=2e-12)
        assert row["coherent"]["joint_home"] > 0.9999
    np.testing.assert_allclose(rows[0]["coherent"]["ports"], rows[0]["dephased"]["ports"])
    assert rows[-1]["dephased"]["joint_home"] == pytest.approx(0.57014382, abs=1e-7)
    assert rows[-1]["unfiltered_uniform"]["joint_home"] == pytest.approx(0.38850324, abs=1e-7)
    assert rows[-1]["stochastic_swap"]["joint_home"] == pytest.approx(0.15111872, abs=1e-7)


def test_cyclic_source_order_and_shortened_monitoring():
    # An exchange-off diagnostic isolates the ordering and the deterministic retry.
    p = acquisition.prepare_layout((1, 3, 10), pulse=0, start=3, direction=-1, rounds=4)
    assert p.launches == (3, 2, 1, 0, 11, 10, 9, 8, 7, 6, 5, 4) * 4
    assert p.rounds == 4 and p.attempt_joint(1)[0][0, 1] == 1
    assert p.history_probability([(1, 3)], accepted=True) == 1
    port_two = acquisition.scored_appointments(p, exposures=(0,))["readings"][0]["coherent"]["ports"]
    np.testing.assert_allclose(port_two, np.eye(12)[2], atol=1e-14)
    with pytest.raises(ValueError):
        p.history_probability([(5, 3)], accepted=True)

    # Post hoc failure with 32 clear rounds, four ordered sweeps and the same free detector.
    tips = [p for p, _, _ in acquisition.gateway_wiring(federation(0), 4)]
    short = acquisition.prepare_layout(tips, rounds=32)
    record = acquisition.scored_appointments(short, home=4, exposures=(37,))["readings"][0]
    assert record["coherent"]["joint_home"] == pytest.approx(0.922409229351, abs=1e-10)


@pytest.mark.parametrize("invalid", [False, True, np.bool_(True), 1.5, float("nan"), "1"])
def test_reject_invalid_indices_and_flags(preparation, invalid):
    with pytest.raises(ValueError):
        acquisition.prepare_layout((invalid, 3, 10))
    with pytest.raises(ValueError):
        acquisition.prepare_layout((1, 3, 10), pulse=invalid)
    for key in ("start", "direction", "rounds"):
        with pytest.raises(ValueError):
            acquisition.prepare_layout((1, 3, 10), **{key: invalid})
    with pytest.raises(ValueError):
        acquisition.scored_appointments(preparation, level=invalid)
    with pytest.raises(ValueError):
        acquisition.gateway_wiring(federation(0), invalid)
    with pytest.raises(ValueError):
        preparation.attempt_joint(invalid)
    with pytest.raises(ValueError):
        preparation.history_probability([(invalid, 1)], accepted=True)
    with pytest.raises(ValueError):
        preparation.history_probability([(1, invalid)], accepted=True)


@pytest.mark.parametrize("flag", [0, 1, 0.5, float("nan"), "accepted", None])
def test_history_requires_boolean_terminal_outcome(preparation, flag):
    with pytest.raises(ValueError):
        preparation.history_probability([], accepted=flag)


@pytest.mark.parametrize("levels", [(), (0, 1), (0, 1, 1, 2), (0, True, 2), (0, 1, 2.0), (0, 1, 3)])
def test_canonical_survey_rejects_incomplete_or_ambiguous_levels(levels):
    with pytest.raises(ValueError):
        acquisition.gateway_acquisition(levels)


def test_no_partial_or_nonfinite_record_can_pass(preparation):
    with pytest.raises(ValueError):
        preparation.history_probability([], accepted=False)
    with pytest.raises(ValueError):
        preparation.history_probability([(1, 1)] * 48, accepted=True)
    with pytest.raises(ValueError):
        preparation.history_probability([(1025, 1)], accepted=True)
    with pytest.raises(ValueError):
        preparation.attempt_joint(49)
    with pytest.raises(ValueError):
        acquisition.prepare_layout(preparation.gateways, rounds=1025)
    for invalid in ((), (True,), (float("inf"),), (float("nan"),), (-0.5,)):
        with pytest.raises(ValueError):
            acquisition.scored_appointments(preparation, exposures=invalid)
    for bad in (float("nan"), float("inf")):
        with pytest.raises(ArithmeticError):
            acquisition._detector_record(np.full(12, bad), slice(0, 12), np.eye(12), 0)
