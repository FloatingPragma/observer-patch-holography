import json
from fractions import Fraction
from math import factorial
from unittest.mock import patch

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import clock
from model.executable_model import cell, federation
from model.experiments.interval_echo_clock import publication
from . import confined_sectors
from . import local_release as release
from . import local_transport
from . import pair_return_bound as pair


def test_compact_flight_and_bounded_retries():
    result = release.appointments()
    assert result["cells_checked"] == {20: 20, 80: 80, 320: 320}
    budget = result["bounded_preparation"]
    assert result["prepared"] == {"vacancy": 4, "denominator": 16,
                                  "occupied": dict(zip(("0", "1", "2", "3", "8", "9", "10", "11"),
                                                       (4, 4, 0, 0, 1, 1, 1, 1), strict=True))}
    assert (budget["maximum_half_contacts"], budget["maximum_signed_queries"], budget["maximum_pi_swaps"]) == (
        384, 64, 30)
    assert Fraction(budget["all_failed"]) == Fraction(3, 4)**16
    assert Fraction(budget["prepared"]) + Fraction(budget["all_failed"]) == 1
    for port, route in budget["failure_routes"].items():
        assert route[0] == int(port) and route[-1] == 5
        assert all(tuple(sorted(edge)) in cell().seams for edge in zip(route, route[1:]))
    preps = result["conditional_on_vacancy"]
    # Same initial port intensities, distinct coherence in the inverse instrument.
    np.testing.assert_allclose(preps["coherent"][0]["addresses_and_quiet"][0],
                               np.mean(preps["resolved"][0]["addresses_and_quiet"], axis=0))
    assert np.mean(preps["resolved"][0]["inverse_addresses_and_quiet"], axis=0)[5] == pytest.approx(1/16)
    state, _ = release.prepare()
    state[[4, 5]] = state[[5, 4]]
    with pytest.raises(ArithmeticError):
        release.all_cells(state)
    with pytest.raises(ValueError):
        release.all_cells(np.zeros(12, dtype=int))


def test_reader_and_retry_reject_impossible_records():
    state = np.eye(12)[0]
    for vector, home in ((state, -1), (state, 1), (state[:9], 0), (state.reshape(12, 1), 0)):
        with pytest.raises(ValueError):
            release.record(vector, home)
    _, herald = release.prepare()
    for occupied in (herald["occupied"] | {"0": -3, "1": 11}, herald["occupied"] | {"12": 0},
                     {v: n for v, n in herald["occupied"].items() if v != "2"}):
        with pytest.raises(ValueError):
            release.retry_budget(herald | {"occupied": occupied}, 16)


def test_inverse_reader_from_forward_time_native_pulses():
    strength = np.sqrt(2)
    h0, h1 = release.law.local_term_parts()
    words = [release.law.SRC.index(v) for v in ((0, 0), (1, 0), (0, 1))]
    matchings = json.loads(release.PROTOCOL.read_text(encoding="utf-8"))["matchings"]
    frames = []
    for duration, sequence in ((np.pi/2, matchings), (3*np.pi/2, reversed(matchings))):
        u = expm(-1j*duration*(h0+strength*h1))
        np.testing.assert_allclose(np.delete(u[:, words], words, axis=0), 0, atol=2e-13, rtol=0)
        np.testing.assert_allclose(u[words, words[0]], [1, 0, 0], atol=2e-13, rtol=0)
        native = u[np.ix_(words[1:], words[1:])]
        frame = np.eye(12, dtype=complex)
        for layer in sequence:
            for a, b in layer:
                frame[[a, b]] = native @ frame[[a, b]]
        frames.append(frame)
    forward, inverse = frames
    np.testing.assert_allclose(inverse @ forward, np.exp(-2j*np.pi*strength*len(matchings))*np.eye(12),
                               atol=2e-13, rtol=0)
    # All local port effects, including their off-diagonals, and the exterior quiet effect.
    basis = np.eye(24, dtype=complex)
    for j in range(12):
        for k in range(j+1, 13):
            for phase in (1, 1j):
                state = (basis[j]+phase*basis[k])/np.sqrt(2)
                expected = np.r_[abs(inverse @ state[:12])**2, np.linalg.norm(state[12:])**2]
                np.testing.assert_allclose(release.record(state, 0, reverse=True), expected, atol=2e-13, rtol=0)


def test_pair_certificate_and_compressed_return_on_a_small_graph():
    sector = PairSector(federation(0))
    assert pair.positive_clock_controls(sector)["worst_absolute_error_lower"] == .2817
    for home in (1, 13):
        assert Fraction(pair.home_bound(sector, home)["probability_upper"]) < Fraction(94, 100)
    matrix = 5*np.eye(3) + np.array([[1, -.5, -.5], [-.5, .5, 0], [-.5, 0, .5]])
    idx = [0, 2]
    re, im, den = pair.compressed_polynomial(sparse.csr_matrix(matrix), idx)
    approx = np.asarray(re/den, float) + 1j*np.asarray(im/den, float)
    exact = expm(-1j*(matrix-6*np.eye(3)))[np.ix_(idx, idx)]
    assert np.linalg.norm(approx-exact, 2) < 6**(pair.DEGREE+1)/factorial(pair.DEGREE+1)
    with pytest.raises(ValueError):
        pair.compressed_polynomial(sparse.csr_matrix(matrix/3), idx)


def test_exact_certificates_reject_lossy_conversion():
    for value in (5+1j, np.inf, np.nan, np.int64(2**62), np.uint64(2**63+6)):
        with pytest.raises(ValueError):
            pair.compressed_polynomial(sparse.csr_matrix([[value]]), [0])
    for value in (2+1j, np.inf, 3.5):
        with pytest.raises(ArithmeticError):
            confined_sectors.integer_matrix(np.array([[value]]))
    assert int(confined_sectors.integer_matrix(np.array([[2**63+6]], dtype=np.uint64))[0]) == 2**63+6


def test_transfer_with_a_noninteger_source_phase_on_every_chart():
    result = local_transport.transport()
    assert [sum(c.values()) for c in result["pulse_count_by_directed_interface"].values()] == [60, 240, 960]
    kappa = np.sqrt(2)
    h0, h1 = release.law.local_term_parts()
    u = expm(-1j*np.pi*(h0+kappa*h1))
    idx = [release.law.SRC.index(v) for v in ((1, 0), (0, 1))]
    driven = u[np.ix_(idx, idx)]
    np.testing.assert_allclose(u[:, release.law.SRC.index((0, 0))], np.eye(release.law.DIM)[:, 0],
                               atol=1e-13, rtol=0)
    np.testing.assert_allclose(np.delete(u[:, idx], idx, axis=0), 0, atol=1e-13, rtol=0)
    for chart in result["charts"].values():
        frame = np.eye(24, dtype=complex)[:, list(release.SUPPORT)]
        for a, b in chart["pulses"]:
            frame[[a, b]] = driven @ frame[[a, b]]
        target = np.eye(24)[:, 12+np.array(release.SUPPORT)] * np.exp(-1j*np.pi*kappa*chart["hops_per_rail"])
        np.testing.assert_allclose(frame, target, atol=1e-12, rtol=0)


def test_exception_on_federated_two_load_operator():
    import sympy as sp

    result = confined_sectors.classify()
    assert [(r["one_load_dimension"], r["pair_dimension_at_minus_half"]) for r in result] == [
        (5, 2), (2, 0), (2, 0)]
    assert result[0]["cells"] == {20: 0, 80: 24, 320: 104}
    local = PairSector(cell())
    g = local.hamiltonian(-.5)
    ids = np.flatnonzero(np.isin(local.aa, [0, 1, 2]) | np.isin(local.bb, [0, 1, 2]))
    modes = confined_sectors.confined_at(confined_sectors.integer_matrix(2*g.toarray()),
                                        sp.eye(local.P)[ids.tolist(), :])
    assert modes.shape == (local.P, 2)
    # Act with each real touching seam of the extended 80-cell board, including outward moves.
    b = federation(1)
    for column in range(modes.cols):
        amp = {(60+int(a), 60+int(d)): w for a, d, w in zip(local.aa, local.bb, modes[:, column], strict=True) if w}
        evolved = {}
        for (a, d), w in amp.items():
            for v, other in ((a, d), (d, a)):
                for dest in b.neighbours[v]:
                    if dest != other:
                        target = tuple(sorted((dest, other)))
                        evolved[target] = evolved.get(target, 0) - w
                        evolved[a, d] = evolved.get((a, d), 0) + w
            if d in b.neighbours[a]:
                evolved[a, d] += 2*w
        assert {key: w for key, w in evolved.items() if w} == {key: 12*w for key, w in amp.items()}


@pytest.mark.parametrize("strength", [0., .5, 1., 1.5, 4., 128.])
def test_archived_clock_journeys_across_strengths(strength):
    """Compare the complete publications to the archived runs at strength eight."""
    archived = json.loads((publication.HERE / "publication.json").read_text(encoding="utf-8"))
    with patch.object(clock, "KAPPA", strength):
        clock._primitives.cache_clear()
        clock._dimer_pulse.cache_clear()
        try:
            observed = publication.run()
        finally:
            clock._primitives.cache_clear()
            clock._dimer_pulse.cache_clear()
    for actual, expected in zip(observed["specimens"], archived["specimens"], strict=True):
        for observer in ("home", "traveller"):
            a, b = (row["endpoint"][observer].pop("probabilities") for row in (actual, expected))
            np.testing.assert_allclose(a, b, atol=1e-12, rtol=0)
        assert actual == expected
    assert observed["observed_transition_counts"] == archived["observed_transition_counts"]
