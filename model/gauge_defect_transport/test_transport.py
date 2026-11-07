"""Exact completeness, endpoint covariance and fermionic-sign controls."""
import pytest
import sympy as s

from . import transport as tr
from ..defect_field_bridge import full_domain


@pytest.mark.parametrize('completed', [False, True])
def test_all_literal_blocks_have_complete_gauge_covariant_outcomes(completed):
    expected = {-4: 18, -3: 2, 1: 30, 2: 3, 6: 10}
    for charge in tr.NAMES:
        value = tr.literal_transport(charge, completed)
        maps = value['maps']
        rank = int(s.trace(value['source_projector']))
        assert rank == (expected[charge] if completed else value['dimension'])
        assert value['endpoint_checks'] == 72
        identity = s.SparseMatrix(s.eye(maps['forward'].rows))
        for direction in ('forward', 'reverse'):
            k, h = maps[direction], maps[direction+'_hold']
            assert tr.link.zero(k.H*k+h.H*h-identity)
            assert not tr.link.zero(k.H*k-identity)
        w = maps['reversible']
        assert w == w.H
        assert tr.link.zero(w*w-identity)


def test_coherent_rotation_retains_complement_and_angle_is_observable():
    value = tr.literal_transport(6)
    maps = value['maps']
    n = value['normalized'].rows
    source = s.SparseMatrix(2*n, 1, {(n, 0): 1})
    target = maps['forward']*source
    assert target.H*target == s.eye(1)
    for c, a in ((s.Rational(3, 5), s.Rational(4, 5)),
                 (s.Rational(5, 13), s.Rational(12, 13))):
        u = tr.coherent_transport(maps, c, a)
        output = u*source
        assert tr.link.zero(output-(c*source-s.I*a*target))
        assert (output.H*maps['target']*output)[0] == a*a
        complement = s.SparseMatrix.eye(2*n)-maps['source']-maps['target']
        assert tr.link.zero(u*complement-complement)
    # Both complete covariant choices are distinguishable on the same input.
    assert s.Rational(4, 5)**2 != s.Rational(12, 13)**2


def test_link_vacuum_reverse_is_retained_as_a_hold():
    value = tr.literal_transport(6)
    maps = value['maps']
    vacuum_left = s.SparseMatrix(maps['forward'].rows, 1, {(0, 0): 1})
    assert maps['reverse']*vacuum_left == s.zeros(vacuum_left.rows, 1)
    assert maps['reverse_hold']*vacuum_left == vacuum_left


@pytest.mark.parametrize('charge', [-4, -3, 1, 2, 6])
def test_angle_control_uses_complete_gauge_singlet_preparation(charge):
    result = tr.neutralized_control(charge)
    assert result['endpoint_singlet_checks'] == 48
    assert result['same_input_and_detector_rotation_probabilities'] == ['16/25', '144/169']


def test_wrong_normalization_fails_partial_isometry():
    actual = tr.literal_transport(1)['normalized']
    with pytest.raises(ValueError, match='partial isometry'):
        tr.complete_partial_transport(actual/2)


def test_omitting_link_transport_breaks_local_covariance():
    blocks, currents, _, _, _, left, _ = tr._literal_data()
    selected = blocks[1]
    identity = s.SparseMatrix.eye(60*len(selected))
    failures = 0
    for g, lx in zip(currents, left, strict=True):
        matter = s.SparseMatrix(s.kronecker_product(g.extract(selected, selected), s.eye(60)))
        link_action = s.SparseMatrix(s.kronecker_product(s.eye(len(selected)), lx))
        failures += not tr.link.zero((matter+link_action)*identity-identity*link_action)
    assert failures > 0


@pytest.mark.parametrize('flavors', [2, 3])
@pytest.mark.parametrize('endpoints,spectator', [((0, 2), 1), ((0, 1), 2)])
def test_car_sign_with_arbitrary_internal_transport_matches_native_crossings(flavors, endpoints, spectator):
    x, y = endpoints
    wanted = -1 if x < spectator < y else 1
    for origin, target in ((x, y), (y, x)):
        for b in range(flavors):
            for a in range(flavors):
                for label in range(flavors):
                    mask = (1 << (origin*flavors+b)) | (1 << (spectator*flavors+label))
                    output, sign = tr.car_hop_sign(mask, (origin, b), (target, a), flavors)
                    assert sign == wanted
                    assert output == (1 << (target*flavors+a)) | (1 << (spectator*flavors+label))
                    if a == b:
                        coin = 0 if origin == x else 1
                        for law in ('C', 'E'):
                            _, actual, native_sign, _ = full_domain.branch_action(
                                mask, law, coin, sites=3, flavors=flavors, seam=endpoints)
                            assert actual == output
                            assert native_sign == sign


@pytest.mark.parametrize('bad', [s.Matrix([[.5]]), s.Matrix([[s.oo]]),
                                s.Matrix([[s.Symbol('a')]]), s.zeros(1, 2)])
def test_invalid_matrix_rejected(bad):
    with pytest.raises(ValueError):
        tr.complete_partial_transport(bad)


@pytest.mark.parametrize('cosine,sine', [(1, 1), (.6, .8), (s.oo, 0), (s.I, 0)])
def test_invalid_rotation_rejected(cosine, sine):
    with pytest.raises(ValueError):
        tr.coherent_transport(tr.literal_transport(6)['maps'], cosine, sine)
