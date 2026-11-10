"""Exact source-code/current/gauge attachment and adversarial boundary checks."""
import pytest
import sympy as s
from .model import (source_descriptor,currents,encode_basis,native_slot_current,
                    encoded_fock_current,local_encoding_checks,quantum_link_checks,
                    background_covariance_checks,hopping_coefficients,link_generator)
from ..physical_attachment.chiral import representation


def test_complete_local_current_source():
    report=local_encoding_checks()
    assert report['local_F30_current_columns_checked']==12*(1+30+435+3)
    assert report['native_idle_K_on_each_object_cloud']=='0'
    assert source_descriptor()['full_site_Fock_dimension']==str(2**30)
    assert not source_descriptor()['ordinary_qubit_geometric_local_compiler_proved']


def test_bosonic_antisymmetry_countercontrol():
    generator=currents()[0]
    pair=next((a,b) for a in range(30) for b in range(a+1,30) if generator[a,b]!=0)
    proper=encode_basis(pair)
    wrong={k:abs(v) for k,v in proper.items()}
    assert native_slot_current(proper,generator)==encoded_fock_current(pair,generator)
    assert native_slot_current(wrong,generator)!=encoded_fock_current(pair,generator)


def test_quantum_link_ward_on_full_normal_ordered_car():
    result=quantum_link_checks()
    assert result['normal_ordered_full_Fock_Ward_checks']==72
    assert result['wrong_orientation_controls_rejected']==2
    assert not result['quantum_link_identified_with_identity_background']


def test_background_current_and_exact_input_domains():
    assert background_covariance_checks()['exact_nonabelian_background_link_covariance_axes']==3
    with pytest.raises(ValueError):hopping_coefficients(0,1.0)
    with pytest.raises(ValueError):encode_basis((1,1))
    with pytest.raises(ValueError):link_generator(s.I*representation()[2][0],0)
