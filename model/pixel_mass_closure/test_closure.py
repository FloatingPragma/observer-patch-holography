"""Crosschecks for the declared closure, excluded census, and complete memory."""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import census,closure,exclusion,memory,packet,replay


@pytest.fixture(scope='module')
def result():
    return packet.report()


def test_census_indices_and_normalization():
    row=census.coefficients()
    assert row['indices_computed_from_native_15_state_operators']
    assert sum(x['dimension'] for x in row['native_joint_current_blocks'])==15
    assert row['one_copy_indices']=={'SU3':'2','SU2':'2','sum_dimension_Y_squared':'10/3'}
    assert row['SM_three_families_one_complex_Higgs']==['41/10','-19/6','-7']
    assert row['MSSM_three_chiral_families_two_chiral_Higgs']==['33/5','1','-3']


def test_exclusion_is_global_on_stated_positive_ordered_domain():
    row=exclusion.certificate()
    assert row['domain']['P_lower_inclusive']=='1'
    assert float(row['scale_exclusion_margin_at_x_6']['lo'])>1
    assert float(row['heat_numerator_upper_enclosure_at_alpha_2_alpha_3_1_over_6']['hi'])<1/80
    assert row['no_positive_ordered_solution']
    assert not row['excludes_other_transmutation_or_threshold_designs']


@pytest.mark.parametrize('design',closure.DESIGNS)
def test_banach_and_implicit_certificates(result,design):
    r=next(x for x in result['closures'] if x['design']==design)
    assert float(r['lipschitz_upper'])<.00734
    assert r['strict_self_map'] and r['infinite_heat_tails_included']
    checks=r['implicit_and_mass_enclosures']
    assert checks['alpha_u_verification']['endpoint_signs_verified']
    assert checks['alpha_u_verification']['R_u_sign_definite']
    assert checks['mz_verification']['endpoint_signs_verified']
    assert checks['mz_verification']['h_m_sign_definite']
    for field in ('P','alpha','alpha_inv'):
        value=F(result['independent_simultaneous_root_and_quadrature_replay'][design][field])
        box=r['fixed_point_'+field]
        assert F(box['lo'])<value<F(box['hi'])


@pytest.mark.parametrize('design',closure.DESIGNS)
def test_actual_source_law_matches_retained_MCPR(result,design):
    path=closure.ROOT/'code/particles/leptons/derive_charged_mcpr_completion_conditional.py'
    spec=importlib.util.spec_from_file_location('_closure_test_mcpr',path)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    row=next(x for x in result['closures'] if x['design']==design)
    with old.mp.workdps(60):
        p=old.mp.mpf(row['display']['p']);u=old.mp.mpf(row['display']['alpha_u'])
        geom=old.icosahedral_incidence_solve()
        fixed=old.response_fixed_point(p,u,old.register_dimensions(geom))
        target=old.charged_triple(p,u,fixed,geom)['masses_over_v']
        b=closure.C.MpBackend(60)
        actual=closure.mass_coefficients(b,b.num(str(p)),b.num(str(u)))
        assert max(abs(float(x/y)-1) for x,y in zip(actual,target))<2e-14


@pytest.mark.parametrize('design',closure.DESIGNS)
def test_complete_memory_inverse_over_root_box(result,design):
    row=next(x for x in result['closures'] if x['design']==design)
    certificate=row['memory_admissibility']
    assert certificate['cubic_root_count']==3
    assert certificate['quadratic_root_count']==2
    assert len(row['complete_record_controls'])==2
    for control,box in zip(row['complete_record_controls'],certificate['inverse_branches']):
        for field in ('p','q','gain_over_v'):
            assert float(F(box[field][0]))<control[field]<float(F(box[field][1]))
        assert control['retained_records_per_input']==512
        assert control['retained_old_inputs']==12


@pytest.mark.parametrize('design',closure.DESIGNS)
def test_top_is_retained_in_same_mass_response(result,design):
    row=next(x for x in result['closures'] if x['design']==design)
    assert float(row['display']['alpha_inv'])>float(row['same_coordinate_top_omission_control']['inverse_response_without_top'])


def test_no_docs_or_central_circular_evidence():
    paths=packet.source_bindings()
    assert all(not p.startswith('docs/') for p in paths)
    assert 'model/verify.py' not in paths
    assert all(not p.startswith('model/completion/') for p in paths)


def test_no_physical_or_uniqueness_promotion(result):
    assert not result['scope']['physical_particle_map_complete']
    assert not result['scope']['complete_physical_SM_Thomson_response']
    assert not result['scope']['native_current_spectral_response_derived']
    assert not result['scope']['source_law_uniquely_selected']
    assert result['canonical_physical_premises_discharged']==[]


def test_invalid_design_rejected():
    with pytest.raises(ValueError):closure.Point('unknown')


def test_invalid_root_bracket_rejected():
    with pytest.raises(ArithmeticError):memory.bracket(lambda x:x*x+1,-1,1,'1e-5')


@pytest.mark.parametrize('field',['fixed_point_P','fixed_point_alpha_inv'])
def test_independent_replay_rejects_root_mutation(result,field):
    value=deepcopy(result)
    value['closures'][0][field]={'lo':'1','hi':'1.1'}
    with pytest.raises(ValueError):replay.verify_report(value)


def test_validator_rejects_missing_source(result):
    value=deepcopy(result);value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError,match='source bindings'):packet.validate(value)


def test_independent_replay_without_producer_import(result,tmp_path):
    path=tmp_path/'report.json';path.write_text(json.dumps(result))
    script='''
import builtins,json,sys
original=builtins.__import__
def guarded(name,*args,**kwargs):
    if name.endswith(('pixel_mass_closure.closure','pixel_mass_closure.packet','interval_contraction_certificate')):
        raise RuntimeError('producer import forbidden')
    return original(name,*args,**kwargs)
builtins.__import__=guarded
from model.pixel_mass_closure.replay import verify_report
verify_report(json.load(open(sys.argv[1])))
'''
    subprocess.run([sys.executable,'-c',script,str(path)],cwd=closure.ROOT,check=True,capture_output=True,text=True)


def test_cli_check_reads_requested_path(monkeypatch,tmp_path):
    path=tmp_path/'requested.json';path.write_text('{"marker": 7}')
    captured=[]
    monkeypatch.setattr(sys,'argv',['pixel_mass_closure','--check',str(path)])
    monkeypatch.setattr(packet,'validate',captured.append)
    packet.main()
    assert captured==[{'marker':7}]
