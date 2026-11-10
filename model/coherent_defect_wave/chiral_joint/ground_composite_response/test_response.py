"""Complete physical images, independent weight diagnostics and hostile omissions."""
import numpy as np
import pytest
import scipy.linalg
import sympy as s

from . import model,spectrum,apparatus
from ..higgs_composites import model as car


def test_complete_actual_ground_response_all_settings_and_dark_states():
    result=spectrum.certificate()
    assert result['complete_joint_record_count']==44
    assert result['positive_probability_records']==40
    rows=result['recorded_settings']
    assert len(rows)==4
    for setting in rows:
        records=setting['complete_target_energy_records']
        assert len(records)==9
        assert sum(row['complete_target_spectral_rank'] for row in records)==26
        dark=[row for row in records if row['target_energy_interval']==['1','1']]
        assert len(dark)==1 and dark[0]['complete_target_spectral_rank']==14
        assert all(s.Rational(row['joint_probability_interval'][0])>0 for row in records)
        assert s.Rational(setting['joint_failure_probability_interval'][0])>0
        assert not setting['conditioned_success_renormalization_applied']
    # Actual full matrices give conjugation as well as spin equality.
    assert all(row['complete_target_energy_records']==rows[0]['complete_target_energy_records'] for row in rows)
    assert all(model.moments(row['kind'],row['probe_spin'])[3]==model.moments('remove',0)[3] for row in rows)


def test_independent_generalized_eigenvector_weight_diagnostic():
    data=model.complete_map('remove',0);records=spectrum.response('remove',0)['complete_target_energy_records']
    hs=np.array(data['source_gram']*data['source_H'],dtype=float)
    gs=np.array(data['source_gram'],dtype=float)
    ht=np.array(data['target_gram']*data['target_H'],dtype=float)
    gt=np.array(data['target_gram'],dtype=float)
    es,vs=scipy.linalg.eigh(hs,gs);et,vt=scipy.linalg.eigh(ht,gt)
    assert abs(es[1]-es[0])<2e-10 and es[2]-es[1]>0.1
    amplitudes=vt.T@(gt@(np.array(data['J'],dtype=float)@vs[:,:2]))
    for row in records:
        interval=list(map(s.Rational,row['target_energy_interval']));center=float(sum(interval)/2)
        ids=np.flatnonzero(np.abs(et-center)<2e-10)
        assert len(ids)==row['complete_target_spectral_rank']
        weight=np.sum(amplitudes[ids,:]**2)/32
        low,high=map(lambda value:float(s.Rational(value)),row['joint_probability_interval'])
        assert low-2e-10<=weight<=high+2e-10
    assert data['J_adjoint']!=data['J'].T


def test_omitted_dark_spectrum_and_wrong_charge_image_rejected():
    data,p,q,_=model.moments('remove',0)
    with pytest.raises(ValueError,match='omits part of the full physical spectrum'):
        model.validate_annihilator(data['target_H'],q.exquo(s.Poly(q.gen-1,q.gen)))
    source=model.source(1);target=model.target(2,0)
    with pytest.raises(ValueError,match='omitted full composite image'):
        model.polynomial_map(source['keys'],(),car.neutral_composite(0))
    with pytest.raises(ValueError,match='omitted full composite image'):
        model.polynomial_map(source['keys'],target['keys'],car.adjoint(car.neutral_composite(0)))
    with pytest.raises(ValueError,match='leaves the retained quark module'):
        model.polynomial_map(source['keys'],target['keys'],car.monomial(((0,False),)))


def test_full_failure_parity_reference_and_preparation_law():
    result=apparatus.certificate()
    assert result['complete_record_count']==44
    assert result['local_parity_preserving_unitary_dilation_dimension']==40
    assert result['ground_selection_is_declared_preparation_law']
    assert not result['arbitrary_input_ground_filter_success_probability_supplied']
    assert all(row['retained_failure_has_no_energy_postselection'] for row in result['settings'])
    assert {(row['reference_success_raw_L'],row['reference_success_twice_Sz'])
            for row in result['settings']}=={(-1,-1),(-1,1),(1,-1),(1,1)}


def test_missing_declared_and_transitive_inputs_fail_closed(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action,'SOURCE_FILES',(*custody.action.SOURCE_FILES,'missing-ground-response.json'))
        with pytest.raises(FileNotFoundError,match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody,'import_closure',lambda seeds:set(seeds)|{'missing-ground-response-import.py'})
    with pytest.raises(FileNotFoundError,match='missing ground response import'):
        custody.input_hashes()
