from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as s

from . import certificate as c


def test_committed_receipt_replays_sources_and_all_matrix_unit_controls():
    assert c.run() == json.loads((c.HERE/"results.json").read_text())


def test_basis_tests_alone_miss_destroyed_coherence():
    _, maps, ranks, basis = c.controls()[1]
    bs, cs = c.probe_scores(maps, ranks, basis)
    assert bs == (1, 1, 1)
    assert cs == (s.Rational(1, 2),)
    assert not c.outcome_maps_equal(maps, tuple((p,) for p in c.projectors(ranks)), ranks)


def test_coherent_test_alone_misses_repeatable_corner_twist():
    _, maps, ranks, basis = c.controls()[2]
    bs, cs = c.probe_scores(maps, ranks, basis)
    assert bs == (0, 0, 1) and cs == (1,)
    for branch, p in zip(maps, c.projectors(ranks)):
        assert c.zero(sum((k.H*k for k in branch), s.zeros(3))-p)
        assert all(c.zero(p*k-k) for k in branch)


def test_unread_equality_does_not_hide_wrong_outcome_labels():
    ranks = (2, 1)
    p = c.projectors(ranks)
    maps = ((p[1],), (p[0],))
    assert c.probe_scores(maps, ranks) == ((0, 0, 0), (0,))
    x = s.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert sum((k*x*k.H for branch in maps for k in branch), s.zeros(3)) == sum(
        (k*x*k.H for k in p), s.zeros(3))


def test_complex_phase_context_needs_no_duplicate_coherent_probe():
    _, maps, ranks, basis = c.controls()[4]
    assert c.probe_scores(maps, ranks, basis) == ((1, 1), ())
    _, twisted, _, _ = c.controls()[5]
    assert c.probe_scores(twisted, ranks, basis) == ((0, 0), ())


def test_nonminimal_kraus_representations_are_not_rejected():
    ranks = (2, 1)
    maps = tuple((s.Rational(3, 5)*p, s.Rational(4, 5)*p)
                 for p in c.projectors(ranks))
    assert c.probe_scores(maps, ranks) == ((1, 1, 1), (1,))


@pytest.mark.parametrize("bad", [((s.eye(2)/2,),), ((),), ((s.eye(3),),)])
def test_missing_mass_or_wrong_domain_is_rejected(bad):
    with pytest.raises(ValueError):
        c.probe_scores(bad, (2,))


def test_rounded_exact_inputs_are_rejected():
    with pytest.raises(TypeError):
        c.calibration_bound((2,), [0.0, 0], [0])
    with pytest.raises(TypeError):
        c.probe_scores(((s.eye(2)*s.Float(1.0),),), (2,))
    with pytest.raises(ValueError):
        c.projectors((2,), [[1, 1], [0, 1]])


@pytest.mark.parametrize("ranks,bs,cs", [((0,), [], []), ((True,), [0], []),
    ((2,), [0], [0]), ((2,), [0, 0], []), ((1,), [-1], []), ((1,), [2], [])])
def test_incomplete_or_invalid_probability_contracts_are_rejected(ranks, bs, cs):
    with pytest.raises(ValueError):
        c.calibration_bound(ranks, bs, cs)


def test_outward_roots_and_nontrivial_error_budget():
    for q in (F(0), F(1, 9), F(2), F(7, 13), F(1, 2**100)):
        upper = c.sqrt_upper(q)
        assert upper*upper >= q
        assert upper <= c.sqrt_upper(q, bits=32)
    bound = c.calibration_bound((2, 1), [F(1, 2**40)]*3, [F(1, 2**40)])
    assert 0 < bound < F(1, 10000)
    assert c.calibration_bound((2, 1), [0, 0, 0], [0]) == 0
    assert c.calibration_bound((2, 1), [1, 1, 1], [1]) == 1
    assert c.history_bound(bound, 100, F(1, 100000)) == F(1, 100000)+100*bound
    assert c.history_bound(F(1, 2), 3) == 1


def test_small_exact_rotation_is_bounded_on_entangled_reference_inputs():
    ranks = (2, 1)
    t = s.Rational(1, 1024)
    a, b = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    rotation = s.Matrix([[a, -b, 0], [b, a, 0], [0, 0, 1]])
    ideal = tuple((p,) for p in c.projectors(ranks))
    actual = tuple((rotation*p,) for p in c.projectors(ranks))
    bs, cs = c.probe_scores(actual, ranks)
    delta = float(c.calibration_bound(ranks, [1-x for x in bs], [1-x for x in cs]))
    assert 0 < delta < 1
    rng = np.random.default_rng(719)
    for vector in [np.eye(3).ravel()/np.sqrt(3)] + [
            rng.normal(size=9)+1j*rng.normal(size=9) for _ in range(12)]:
        vector = vector/np.linalg.norm(vector)
        rho = np.outer(vector, vector.conj())
        distance = 0.0
        for aa, ii in zip(actual, ideal):
            k = np.kron(np.array(aa[0]).astype(complex), np.eye(3))
            p = np.kron(np.array(ii[0]).astype(complex), np.eye(3))
            difference = k@rho@k.conj().T-p@rho@p.conj().T
            distance += np.linalg.norm(difference, ord="nuc")/2
        assert distance <= delta+1e-12


def test_nonprojective_effects_satisfy_uniform_calibration_bound():
    _, maps, ranks, basis = c.controls()[6]
    ps = c.projectors(ranks, basis)
    for branch, p in zip(maps, ps):
        effect = sum((k.H*k for k in branch), s.zeros(3))
        assert not c.zero(effect-p)
        assert not c.zero(effect*effect-effect)
    bs, cs = c.probe_scores(maps, ranks, basis)
    delta = float(c.calibration_bound(ranks, [1-x for x in bs], [1-x for x in cs]))
    assert 0 < delta < 0.02
    rng = np.random.default_rng(837)
    for vector in [np.eye(3).ravel()/np.sqrt(3)] + [
            rng.normal(size=9)+1j*rng.normal(size=9) for _ in range(8)]:
        vector = vector/np.linalg.norm(vector)
        rho = np.outer(vector, vector.conj())
        distance = 0.0
        for branch, p in zip(maps, ps):
            ideal = np.kron(np.array(p).astype(complex), np.eye(3))
            actual = np.zeros((9, 9), dtype=complex)
            for k in branch:
                op = np.kron(np.array(k).astype(complex), np.eye(3))
                actual += op@rho@op.conj().T
            distance += np.linalg.norm(actual-ideal@rho@ideal.conj().T, ord="nuc")/2
        assert distance <= delta+1e-12
