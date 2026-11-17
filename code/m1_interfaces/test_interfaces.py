"""Independent finite replay, adversarial receipts and analytic controls."""

from copy import deepcopy
from fractions import Fraction as F
from fnmatch import fnmatchcase
from itertools import product
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest
import yaml

from . import check, model, verify

ROOT=Path(__file__).resolve().parents[2]


@pytest.fixture(scope='module')
def packet():
    return check.strict_load(Path(__file__).with_name('receipt.json'))


def test_complete_receipt(packet):
    result=verify.verify(packet)
    assert result['graphs']==12 and result['cut_cases']==96
    assert result['site_instances']==4*(12**3+16**3+24**3)
    assert result['ordered_read_incidences']>0
    assert result['scale_levels']==12 and result['moment_regimes']==3


@pytest.mark.parametrize('value',[None,True,False,0,-1,13,1.0,'1'])
def test_bad_scale_inputs_fail(value):
    with pytest.raises(ValueError): model.scale(value)
    with pytest.raises(ValueError): check.scale(value)


@pytest.mark.parametrize('mutation',[
    'schema','source','parent_scope','missing_scale','false_threshold','boolean_population',
    'missing_regime','zero_action','hide_alias','wrong_moment','missing_graph','missing_mask',
    'relabel_mask','false_cut','double_count','fake_connectivity','fake_bridge_cost','wrong_boundary',
    'drop_boundary_action','double_boundary_action','wrong_action_normalization',
])
def test_forged_or_truncated_receipts_fail(packet,mutation):
    bad=deepcopy(packet)
    row=bad['graphs']['12:U:clipped']
    if mutation=='schema': bad['accepted']=True
    elif mutation=='source': bad['sources'].popitem()
    elif mutation=='parent_scope': next(iter(bad['parent_claims'].values()))['disposition']='physical_M1_closed'
    elif mutation=='missing_scale': bad['scales'].pop('12')
    elif mutation=='false_threshold': bad['scales']['1']['compactness_error_factor']='0'
    elif mutation=='boolean_population': bad['scales']['1']['q']=True
    elif mutation=='missing_regime': bad['moments'].pop('critical')
    elif mutation=='zero_action': bad['moments']['repaired']['uniform_alpha_L_squared']='0'
    elif mutation=='hide_alias': bad['moments']['critical']['one_alias_eigenvalue_L_squared_bounds']=['0','0']
    elif mutation=='wrong_moment': bad['moments']['raw']['coarse']['axis4']+=1
    elif mutation=='missing_graph': bad['graphs'].pop('24:T:periodic')
    elif mutation=='missing_mask': row['masks'].pop('localized')
    elif mutation=='relabel_mask': row['masks']['residue']=deepcopy(row['masks']['plane'])
    elif mutation=='false_cut': row['masks']['localized']['cut_pairs']+=1
    elif mutation=='double_count': row['masks']['plane']['cut_pairs']*=2
    elif mutation=='fake_connectivity': row['bridge']['nearest_neighbor_connected']=False
    elif mutation=='fake_bridge_cost': row['bridge']['changed_sites']=0
    elif mutation=='wrong_boundary': bad['graphs']['12:U:clipped']=deepcopy(bad['graphs']['12:U:periodic'])
    elif mutation=='drop_boundary_action':
        action=row['uniform_action']
        action['constant_quadratic']=0
        action['affine_boundary_quadratic']=0
        action['affine_quadratic']=action['affine_internal_quadratic']
        action['scaled_affine_quadratic']=str(F(6*action['affine_quadratic'],12*action['second_norm_moment']))
    elif mutation=='double_boundary_action': row['uniform_action']['affine_boundary_quadratic']*=2
    else: row['uniform_action']['scaled_affine_quadratic']='0'
    with pytest.raises(ValueError): verify.verify(bad)


@pytest.mark.parametrize('raw',['{"x":1,"x":2}','{"x":NaN}','{"x":Infinity}'])
def test_noncanonical_json_rejected(tmp_path,raw):
    p=tmp_path/'bad.json';p.write_text(raw,encoding='utf-8')
    with pytest.raises(ValueError): check.strict_load(p)


@pytest.mark.parametrize('case',['moment','boundary'])
def test_optimized_cli_rejects_a_semantic_forgery(packet,tmp_path,case):
    bad=deepcopy(packet)
    if case=='moment':
        bad['moments']['critical']['uniform_alpha_L_squared']='0'
        message='false complete moment or alias certificate'
    else:
        action=bad['graphs']['12:U:clipped']['uniform_action']
        action['constant_quadratic']=0
        action['affine_boundary_quadratic']=0
        action['affine_quadratic']=action['affine_internal_quadratic']
        action['scaled_affine_quadratic']=str(F(6*action['affine_quadratic'],12*action['second_norm_moment']))
        message='false finite cut execution'
    p=tmp_path/'forged.json';p.write_text(json.dumps(bad),encoding='utf-8')
    run=subprocess.run([sys.executable,'-O','-m','m1_interfaces.verify','--receipt',str(p)],cwd=ROOT,
                       env=dict(os.environ,PYTHONPATH=str(ROOT/'code')),capture_output=True,text=True)
    assert run.returncode!=0 and message in run.stderr


def test_verifier_cannot_import_producer():
    script="""
import sys
sys.modules['m1_interfaces.model']=None
sys.modules['m1_interfaces.build']=None
from m1_interfaces import verify,check
print(verify.verify(check.strict_load(verify.HERE/'receipt.json')))
"""
    run=subprocess.run([sys.executable,'-c',script],cwd=ROOT,
                       env=dict(os.environ,PYTHONPATH=str(ROOT/'code')),capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr


def test_duplicate_parent_claim_is_rejected(monkeypatch):
    registry=check.strict_load(ROOT/'claims/claim_registry.yaml')
    parent=next(row for row in registry['claims'] if row['claim_id']=='OPH-BH-CROSSING-READ-AREA-LAW')
    registry['claims'].append(deepcopy(parent))
    monkeypatch.setattr(check,'strict_load',lambda path:registry)
    with pytest.raises(ValueError,match='parent claims must occur exactly once'):
        verify.parent_claims()


def test_ball_columns_against_complete_small_vector_enumeration():
    for r in range(9):
        vs=[v for v in product(range(-r,r+1),repeat=3) if sum(x*x for x in v)<=r*r]
        row=check.ball_columns(r,64)
        assert row==model.ball_columns(r,64)
        assert row['degree']==len(vs)
        assert row['moment']==sum(sum(x*x for x in v) for v in vs)
        assert row['fourth']==sum(sum(x*x for x in v)**2 for v in vs)
        assert row['axis2']==sum(v[0]**2 for v in vs)
        assert row['axis4']==sum(v[0]**4 for v in vs)
        assert row['flux']==sum(max(v[0],0) for v in vs)
        assert row['coordinate_cut']==sum(v[0]*(64-abs(v[1]))*(64-abs(v[2]))
                                          for v in vs if v[0]>0)


def test_complete_axis_slices_put_actual_aliases_inside_exact_bounds(packet):
    # Independent full yz-plane census, rather than the column square-root
    # formulas. Only one-dimensional Fourier probes are evaluated; this is
    # not described as executing the q^3 field or its complete spectrum.
    for kind in ('raw','critical','repaired'):
        row=packet['moments'][kind];r=row['r'];m=row['m']
        yz=np.array(list(product(range(-r,r+1),repeat=2)))
        squares=np.sum(yz*yz,axis=1)
        xs=np.arange(-r,r+1)
        slices=np.array([np.count_nonzero(squares<=r*r-x*x) for x in xs])
        assert sum(slices)==row['fine']['degree']
        assert sum(slices*xs**4)==row['fine']['axis4']
        value=np.dot(slices,1-np.cos(2*np.pi*xs/m))
        value+=sum(2*(1-np.cos(2*np.pi*s/m)) for s in row['remaining_axis_lengths'])
        value*=float(F(row['uniform_alpha_L_squared']))
        lower,upper=map(lambda x:float(F(x)),row['one_alias_eigenvalue_L_squared_bounds'])
        assert lower-1e-10<=value<=upper+1e-10


def test_symbol_bound_including_noncentral_alias_cells():
    # The ball lemma is checked across the Brillouin zone; the universal
    # constant and the all-mode theorem are proved in DERIVATION.md.
    r=8
    vs=np.array([v for v in product(range(-r,r+1),repeat=3) if sum(x*x for x in v)<=r*r])
    for theta in product((-np.pi,-2/r,-1/r,0,1/r,2/r,np.pi),repeat=3):
        p=np.array(theta)
        actual=np.sum(1-np.cos(vs@p))
        lower=r**3/49152*min(r*r*float(p@p),1)
        assert actual+1e-10>=lower


def test_exact_cut_change_and_connected_control(packet):
    for config in check.CONFIGS:
        masks=check.labels(config)
        assert not check.connected(config,masks['localized'])
        assert check.connected(config,masks['connected'])
        for family in ('T','U'):
            for topology in ('periodic','clipped'):
                row=packet['graphs'][f'{config[0]}:{family}:{topology}']
                assert abs(row['bridge']['cut_change'])<=row['bridge']['degree_bound']
                assert row['bridge']['changed_sites']>0


def test_uniform_action_matches_independent_edge_form_and_boundary_control(packet):
    for config in check.CONFIGS:
        for family in ('T','U'):
            for topology in ('periodic','clipped'):
                action=packet['graphs'][f'{config[0]}:{family}:{topology}']['uniform_action']
                # Producer uses the actual shifted-array matrix action. The
                # verifier derives the unordered edge energy and absent-slot
                # diagonal potential independently, without exterior records.
                assert model.graph_cuts(config,family,topology)[3]==action
                assert action['affine_quadratic']==action['affine_internal_quadratic']+action['affine_boundary_quadratic']
                assert (action['constant_quadratic']==0)==(topology=='periodic')
                assert (action['affine_boundary_quadratic']==0)==(topology=='periodic')
                assert F(action['scaled_affine_quadratic'])>0


def test_all_new_lean_theorems_are_audited_and_built():
    source=(ROOT/'Lean/Geometry/M1Interfaces.lean').read_text(encoding='utf-8')
    gate=(ROOT/'Lean/Geometry/M1InterfacesAxiomAudit.lean').read_text(encoding='utf-8')
    names=set(re.findall(r'^theorem (\w+)',source,re.M))
    assert names and {'OPH.M1Interfaces.'+n for n in names}==set(re.findall(r'^audit_m1_interfaces (OPH\.\S+)',gate,re.M))
    assert '#guard_msgs in\naudit_m1_interfaces sorryAx' in gate
    assert '#guard_msgs in\naudit_m1_interfaces Lean.ofReduceBool' in gate
    ci=(ROOT/'.github/workflows/lean.yml').read_text(encoding='utf-8')
    assert re.search(r'^\s*"Geometry.M1InterfacesAxiomAudit"\s*$',ci,re.M)


