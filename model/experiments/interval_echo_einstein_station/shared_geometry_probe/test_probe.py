"""Cross-check the registered reduction, the complete mixed launch and both receipts."""

import hashlib
import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from . import flight, measure


def test_registered_reduction_and_inactive_controls():
    board = measure.federation(0)
    sites = measure.layouts(board)[0]
    x = measure.PAULI[0]
    shared = -(np.kron(x, measure.UNIT) + np.kron(measure.UNIT, x)) / 2
    for s in (-1, 1):
        for q in (-1, 1):
            for arm, coefficient in [
                ("low", 32),
                ("high", 32),
                ("quarter", 32),
                ("collar", 33),
                ("disconnected_low", 0),
                ("changed_high", 0),
            ]:
                actual, _ = measure.generator(board, sites, s, q, arm)
                expected = shared.copy()
                expected[2, 2] -= coefficient
                expected -= np.trace(expected) / 4 * np.eye(4)
                np.testing.assert_allclose(actual, expected, atol=3e-14)
    # Partial trace using forward motion, including the imaginary Y readout.
    h, source = measure.generator(board, sites, 1, -1, "collar")
    u = expm(-4j * h)
    joint = u @ np.kron(measure.UNIT / 2, source) @ u.conj().T
    probe = np.trace(joint.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    target = (1 + np.einsum("amn,nm->a", measure.PAULI, probe).real) / 2
    np.testing.assert_allclose(
        target, measure.probabilities(h, source, 4)[-1], atol=1e-12
    )


def test_flight_against_full_forward_density():
    board = flight.formation.federation(0)
    route = flight.formation.region(board)
    edges = flight.monitor_edges(board, route)
    h = flight.formation.hamiltonian(board, route, 1, 8)
    source = flight.prepared_states(board, route, -1)[0]
    effects = flight.backward_effects(h, edges)[0]
    # Propagate all 960 input basis states forward, sum the final head,
    # and weight every initial position by 1/240 with its full source density.
    propagator = expm_multiply(-1j * h, np.eye(4 * board.n, dtype=complex), traceA=0)
    kernels = propagator.reshape(board.n, 4, board.n, 4)[list(edges[0])]
    for k, vec in enumerate(flight.VECS):
        amp = np.einsum("v,vjxa->jxa", vec.conj(), kernels)
        forward = np.einsum("jxa,ab,jxb->", amp, source, amp.conj()).real / board.n
        backward = np.trace(effects[0, k] @ source).real
        np.testing.assert_allclose(forward, backward, atol=4e-13)


def test_enrollment_and_calibration_receipt():
    rx = json.loads((measure.HERE / "results.json").read_text())
    assert (
        rx["protocol_sha256"]
        == hashlib.sha256((measure.HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert rx["number_of_binary_cohorts"] < 100000
    assert len(rx["trials"]) == 9 * 4 * len(measure.ARMS)
    summary = measure.summarize(rx["trials"], rx["source"], rx["clocks"])
    assert summary == rx["summary"]
    assert summary["clock_max_exposure_error"] < 0.0001
    held = summary["held_energy_calibration"]
    assert (
        sum(held[a]["compared_pointer_components"] for a in ("low", "high", "quarter"))
        == 13992
    )
    assert all(
        held[a]["rejected_pointer_components"] == 0 for a in ("low", "high", "quarter")
    )
    assert held["collar"]["rejected_pointer_components"] == 4680
    mixed = summary["mixed_probe"]
    for arm in ("low", "high"):
        assert (
            mixed[arm]["minimum_mixed_state_distance_by_time"][-1] - measure.EPS
            > 0.4265
        )
    for arm in ("changed_low", "changed_high", "disconnected_low", "disconnected_high"):
        assert max(mixed[arm]["maximum_mixed_state_component_by_time"]) < measure.EPS


def test_flight_receipt_bounds_and_controls():
    rx = json.loads((measure.HERE / "flight.json").read_text())
    assert (
        rx["protocol_sha256"]
        == hashlib.sha256(
            (measure.HERE / "FLIGHT_PROTOCOL.md").read_bytes()
        ).hexdigest()
    )
    endpoint_limits = flight.intervals(np.array([0, flight.TRIALS]))
    assert endpoint_limits[0, 0] == 0 and endpoint_limits[1, 1] == 1
    assert endpoint_limits[0, 1] > 0 and endpoint_limits[1, 0] < 1
    totals = np.zeros((3, 2), int)
    correlation = np.zeros(2, int)
    for row in rx["rows"]:
        hits = np.array(row["plus_counts"])
        limits = flight.intervals(hits)
        np.testing.assert_allclose(limits, row["intervals"], atol=1e-12)
        null = 1 / row["ports"]
        outside = (limits[..., 0] > null) | (limits[..., 1] < null)
        totals += outside.sum(axis=(1, 3))
        correlated = (limits[0, ..., 0] > limits[1, ..., 1]) | (
            limits[0, ..., 1] < limits[1, ..., 0]
        )
        correlation += correlated.sum(axis=(0, 2))
    assert totals[0].tolist() == [144, 0]
    assert totals[2].tolist() == [0, 0]
    assert correlation.tolist() == [115, 0]
