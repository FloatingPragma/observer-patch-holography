"""Exact common-link occupation, gauge and representation regression controls."""
from fractions import Fraction as F

import pytest

from .certificate import (algebra_certificate, certificate, lift_one,
                          representation_certificate)
from .fock import FockObserver, FockState, coordinates, mode, wedge_insert
from . import spatial_fock
from . import joint_readout


def test_basis_bijection_and_wedge_sign():
    for index in range(7200):
        assert mode(*coordinates(index))==index
    assert wedge_insert((1,4),3)==((1,3,4),-1)
    assert wedge_insert((1,4),4)==(None,0)


def test_zero_pair_state_keeps_exact_arithmetic():
    import sympy as s
    zero=spatial_fock.PairFockState({},7)
    assert zero.norm_squared()==s.Integer(0)
    assert not zero.norm_squared().has(s.Float)
    assert spatial_fock.exchange(zero).norm_squared()==s.Integer(0)


@pytest.mark.parametrize('occupied',[(1,1),(2,1),(-1,),(7200,)])
def test_invalid_occupation_rejected(occupied):
    with pytest.raises(ValueError):
        FockState({(1,occupied):1})


def test_exact_amplitudes_and_source_sector_required():
    observer=FockObserver()
    with pytest.raises(ValueError):
        FockState({(1,()):0.5})
    with pytest.raises(ValueError):
        observer.step(FockState({(0,()):1}))
    with pytest.raises(ValueError):
        FockObserver('multiply-the-classical-control')


def test_same_one_wave_clock_and_reverse():
    observer=FockObserver()
    state=lift_one(observer.base.embedded_ports(1,(0,5)))
    output=observer.step(state)
    assert output==lift_one(observer.base.step(observer.base.embedded_ports(1,(0,5))))
    assert observer.inverse_step(output)==state


def test_control_architecture_difference_and_one_wave_agreement():
    parity=FockObserver(); presence=FockObserver('presence')
    first=mode(5,3,0,0); second=mode(5,4,0,0)
    one=parity.basis_state(1,(first,))
    assert parity.controls(one)==presence.controls(one)
    two=parity.basis_state(1,(first,second))
    assert parity.controls(two)==two
    assert presence.controls(two)!=two
    assert presence.controls(presence.controls(two),inverse=True)==two


def test_exterior_weights_and_statistics_are_separate():
    weights=representation_certificate()
    assert weights['dimension']==15
    assert sum(row['dimension'] for row in weights['multiplets'])==15
    assert algebra_certificate()['even_auxiliary_composites_exchange_sign']==1


def test_complete_same_source_receipt():
    result=certificate()
    execution=result['two_occupation_execution']
    assert sum(map(F,execution['complete_face_probabilities'].values()))==1
    assert execution['output_support']==16200
    assert execution['fixed_positive_Gauss_empty_wave_dimension']==5
    assert result['physical_scope']['SM_species_identified'] is False
    assert result['spatial_Fock_join']['equals_tensor_wave_partial_SWAP'] is False
    assert [r['probability'] for r in result['joint_translation_phase_readout']['complete_joint_records']]==[
        '1/120','59/120','1/4','1/4']


def test_spatial_hopping_empty_and_one_body_blocks():
    i=mode(5,3,0,0)
    state=spatial_fock.positive_basis((1,3),(i,))
    output=spatial_fock.exchange(state)
    assert output.norm_squared()==1
    assert spatial_fock.exchange(output,inverse=True)==state
    assert [b.norm_squared() for _,b in spatial_fock.number_measurement(output,0)]==[F(1,2),F(1,2)]


def test_spatial_data_reject_wrong_register_or_float():
    with pytest.raises(ValueError):
        spatial_fock.PairFockState({((1,0),()):1})
    with pytest.raises(ValueError):
        spatial_fock.PairFockState({((1,3),()):0.5})


def test_joint_translation_carries_the_fermionic_sign():
    mode_index=mode(2,3,0,0)
    state=spatial_fock.PairFockState({((1,3),(mode_index,7200+mode_index)):1})
    moved=joint_readout.translate(state)
    assert moved==spatial_fock.PairFockState({((3,1),(mode_index,7200+mode_index)):-1})
    assert joint_readout.translate(moved)==state
