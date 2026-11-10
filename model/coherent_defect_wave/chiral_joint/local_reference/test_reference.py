"""Complete physical spectra, full-source controls and hostile omissions."""
from copy import deepcopy
from fractions import Fraction
import pytest
import sympy as s

from . import algebra,characters,physical,spectrum,reference


def test_exact_native_particle_hole_and_complete_collective_bound():
    result=algebra.certificate()
    assert result['joint_complex_Lie_dimension']==18
    assert result['quark_collective_color_and_spin_commutators']==11
    assert result['lepton_squared_transfer_operator_upper_bound']==6
    result=characters.certificate()
    assert result['full_quark_Fock_dimension']==2**24
    assert result['distinct_collective_irreps']==28
    assert result['maximal_module']['dimension']==140
    assert result['maximal_module']['multiplicity']==1
    assert result['all_other_color_admissible_quark_Casimir_upper_bound']=='42'
    # Removing the color condition admits C=95/2 and invalidates the bound 135.
    assert Fraction(135-101,2)**2<95*6


def test_omitted_native_yukawa_word_rejected(monkeypatch):
    words=tuple(algebra.yukawa_words())
    monkeypatch.setattr(algebra,'yukawa_words',lambda:words[1:])
    with pytest.raises(AssertionError):algebra.transitions.__wrapped__()


def test_wrong_fundamental_character_and_omitted_maximum_rejected(monkeypatch):
    actual=algebra.quark_clifford()
    with monkeypatch.context() as patch:
        patch.setattr(characters,'quark_clifford',lambda:(*actual[:-1],actual[-1][:-1]))
        with pytest.raises(AssertionError):characters.inventory.__wrapped__()
    rows=deepcopy(characters.inventory())
    rows[12]['modules']=[m for m in rows[12]['modules'] if m['Casimir']!='54']
    monkeypatch.setattr(characters,'inventory',lambda:rows)
    with pytest.raises(AssertionError):characters.certificate()


def test_all_physical_spin_blocks_dark_states_and_global_gap():
    complete=physical.certificate()
    assert complete['full_maximal_module_physical_dimension']==104
    assert len(complete['complete_spin_blocks'])==8
    result=spectrum.certificate()
    assert len(result['complete_maximal_module_energy_records'])==19
    assert result['dark_energy_one_rank']==56
    assert result['global_physical_ground_rank']==result['first_excited_rank']==2
    assert s.Rational(result['global_physical_gap_interval'][0])>0
    assert s.Rational(result['first_excited_energy_interval'][1])<s.Rational(-1113,100)
    assert not result['unconstrained_or_boundary_compensated_local_ground_proved']
    assert not result['interacting_many_cell_vacuum_proved']


def test_incomplete_current_image_and_false_scalar_penalty_rejected():
    with pytest.raises(AssertionError,match='omitted current image'):
        physical.tensor_map(((0,0,0),),(),quark=s.eye(1))
    keys,basis,gram,form,generator=physical.physical_block(3,1)
    bad=s.MutableSparseMatrix(generator);bad[0,0]+=1
    with pytest.raises(ValueError,match='diagonal block'):
        physical.scalar_split(keys,basis,bad)
    bad=s.MutableSparseMatrix(basis)
    blank=next(i for i,key in enumerate(keys) if key[2]==0)
    excited=next(i for i,key in enumerate(keys) if key[2]!=0)
    bad[blank,0]=bad[excited,0]=1
    with pytest.raises(ValueError,match='mixes scalar excitation'):
        physical.scalar_split(keys,bad,generator)


def test_incomplete_rank_and_overlapping_spectral_claims_rejected():
    records=deepcopy(spectrum.certificate()['complete_maximal_module_energy_records'])
    with pytest.raises(ValueError,match='incomplete physical'):
        spectrum.ordered_records(records[:-1],104)
    records[1]['energy_interval']=records[0]['energy_interval']
    with pytest.raises(ValueError,match='strict spectral separation'):
        spectrum.ordered_records(records,104)


def test_actual_stationary_reference_and_half_filled_native_irreps():
    result=reference.certificate()
    assert result['density']=='P_ground/2'
    assert result['exact_total_mean_native_occupation']=='15'
    assert {k:v['exact_ground_mean_occupation'] for k,v in result['native_irrep_ground_occupations'].items()}=={
        'uc':'3','Q':'6','ec':'1','L':'2','dc':'3'}
    assert result['mean_occupations_also_hold_on_each_pure_ground_spin_state']
    assert not result['bare_one_fermion_reference']


def test_missing_declared_and_transitive_custody_fail_closed(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action,'SOURCE_FILES',(*custody.action.SOURCE_FILES,'absent-local-reference.json'))
        with pytest.raises(FileNotFoundError,match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody,'import_closure',lambda seeds:set(seeds)|{'absent-local-reference-import.py'})
    with pytest.raises(FileNotFoundError,match='missing local reference import'):
        custody.input_hashes()
