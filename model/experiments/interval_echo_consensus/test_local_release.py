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
from . import pair_return_bound as pair


def test_compact_flight_and_bounded_retries():
    result = release.appointments()
    assert result["cells_checked"] == {20: 20, 80: 80, 320: 320}
    budget = result["bounded_preparation"]
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


def test_pair_certificate_and_compressed_return_on_a_small_graph():
    sector = PairSector(federation(0))
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


def test_archived_clock_journeys_with_exchange_alone():
    """Same discrete transcripts at zero. The reference contains the adopted clock's records."""
    archived = json.loads((publication.HERE / "publication.json").read_text(encoding="utf-8"))
    with patch.object(clock, "KAPPA", 0):
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
