"""Exact same-operator joins and independent complete-instrument controls."""
import json
from pathlib import Path

import pytest
import sympy as s

from . import bridge as b


def test_same_retained_current_and_Higgs_basis():
    rep = b.representation()
    P = rep['permutation15']
    assert P.T*P == s.eye(15)
    assert len(rep['currents15']) == 12
    assert set(rep['Higgs_channels']) == {'up', 'down', 'charged_lepton'}
    for old, current in zip(rep['source_currents15'], rep['currents15'], strict=True):
        assert b._zero(current*P-P*old)
    W = rep['negative_one_particle_embedding']
    assert W.shape == (180, 90)
    assert W.H*W == 2*s.eye(90)
    assert W[:90, :] == -W[90:, :]


@pytest.mark.parametrize('modes', [1, 2, 3, 4])
def test_dressed_CAR_on_complete_parity_correlated_code(modes):
    W, P, fields = b.parity_code(modes)
    d = P.rows
    Z = s.diag(1, -1)
    assert W.H*W == 2*s.eye(d)
    for a, field in enumerate(fields):
        dressed = s.kronecker_product(Z, field)
        assert dressed*W == W*field
        for other_index, other in enumerate(fields):
            assert field*other.H+other.H*field == (s.eye(d) if a == other_index else s.zeros(d))
            assert field*other+other*field == s.zeros(d)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_native_field_effects_and_reference(law):
    value = b.small_complete_instrument(law, modes=1)
    size = value['encoded_dimension']
    effect = s.zeros(size)
    for row, full, encoded in value['branches']:
        assert full*value['code'] == value['code']*encoded
        effect += encoded.H*encoded
        if row.abort:
            assert not row.inputs
    assert effect == 2*s.eye(size)
    # Exact equality of linear maps extends to any reference, including a
    # coherent two-column reference; no probability-only comparison.
    state = s.zeros(size, 2)
    state[0, 0], state[1, 1], state[-1, 0] = 1, s.I, s.Rational(1, 3)
    for _, full, encoded in value['branches']:
        assert full*(value['code']*state) == value['code']*(encoded*state)


def test_vacuum_exposes_difference_between_fixed_chi_and_Fock_parity():
    W, P, _ = b.parity_code(2)
    Q = s.kronecker_product(s.Matrix([[0, 1], [1, 0]]), s.eye(4))
    assert Q*W[:, 0] == W[:, 0]  # vacuum
    assert Q*W[:, 1] == -W[:, 1]  # one particle
    assert Q*W[:, 3] == W[:, 3]  # two particles
    assert Q*W != -W
    value = b.statistics_certificate()
    assert value['two_occupied_modes_control'] == {
        'native_odd_loop_under_this_code': 1, 'actual_fermionic_mode_exchange': -1}


def test_second_quantizing_Kraus_maps_is_not_TP():
    # On a one-particle line, K0=K1=1/sqrt(2) is a complete instrument.
    # Exterior functor Gamma(K)=1 on the vacuum, so the two lifted vacuum
    # branches have total probability two. Our construction never does it.
    gamma = s.diag(1, 1/s.sqrt(2))
    assert 2*gamma.H*gamma == s.diag(2, 1)
    assert 2*gamma.H*gamma != s.eye(2)


def test_receipt_roundtrip():
    path = Path(b.__file__).with_name('representation_results.json')
    assert b.certificate() == json.loads(path.read_text())
