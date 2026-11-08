import copy
from fractions import Fraction as F
import json
import math
import subprocess
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import erfc

from . import bounds,finite,packet


@pytest.mark.parametrize('alpha,omega',[(.1,.2),(.1,2.),(.7,1.3),(2.,3.)])
def test_gaussian_filter_multiplier_and_absolute_tail(alpha,omega):
    imaginary=-quad(lambda t:erfc(math.sqrt(alpha)*t)*math.sin(omega*t),0,np.inf,epsabs=1e-11)[0]
    expected=-(-math.expm1(-omega*omega/(4*alpha)))/omega
    assert abs(imaginary-expected)<1e-10
    for t in (.1,1.,2.):
        tail=quad(lambda u:erfc(math.sqrt(alpha)*u),t,np.inf,epsabs=1e-12)[0]
        exact=math.exp(-alpha*t*t)/math.sqrt(math.pi*alpha)-t*erfc(math.sqrt(alpha)*t)
        assert abs(tail-exact)<1e-11
        assert tail<=math.exp(-alpha*t*t)/math.sqrt(math.pi*alpha)+1e-12


def test_source_attachment_and_gap_normalization():
    row=finite.exact_source_attachment()
    assert row['integer_Hadamard_identity_exact']
    assert row['field_parameter_h']=='B/2' and row['full_gap_lower']=='B/8'
    assert not row['scalar_bond_offset_changes_ground_or_gap']


def test_actual_ground_derivative_requires_negative_i_flow_sign():
    full,cut,v=finite.matrices();row=finite.flow((full+cut)/2,v,.0001)
    assert np.linalg.norm(row['vector_derivative']-row['flow_vector'])<1e-12
    assert np.linalg.norm(row['vector_derivative']+row['flow_vector'])>.1
    assert abs(np.vdot(row['ground'],row['vector_derivative']))<1e-13


def test_actual_source_marginal_duality_has_no_tomography_factor():
    full,cut,v=finite.matrices();ef,uf,pf=finite.projector(full);ec,uc,pc=finite.projector(cut)
    delta=finite.marginal_first_two(uf[:,0])-finite.marginal_first_two(uc[:,0])
    eigen,basis=np.linalg.eigh(delta);sign=(basis*np.sign(eigen))@basis.conj().T
    assert abs(np.trace(sign@delta).real/2-np.linalg.norm(delta,ord='nuc')/2)<1e-14
    assert np.linalg.norm(sign,2)<=1+1e-14


def test_rational_certificate_bounds_the_analytic_expression():
    for j in (F(1,3),F(1),F(3)):
        for b in (F(1,4),F(1),F(8)):
            for d in (1,7,103,10000):
                r=bounds.marginal(j,b,d,7,12);g=float(r['gap_lower']);t=float(r['time_cutoff'])
                exact=float(j)*12*(7*t/2*math.exp(-d/3)+(math.sqrt(2*t/(math.pi*g))+1/g)*math.exp(-g*t/2))
                assert exact<=float(r['uncapped_half_trace_upper'])+1e-12
                assert r['gaussian_alpha']==r['gap_lower']/(2*r['time_cutoff'])
                assert not r['exterior_volume_in_bound']


def test_box_planner_has_certified_accuracy_without_state_allocation(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('source allocated in bound planner')
    monkeypatch.setattr(finite.source,'matrix',forbidden)
    for b in (F(1),F(1,4),F(1,16)):
        r=bounds.choose_box(F(1),b,1,27,F(1,10000))
        assert 0<r['half_trace_upper']<=F(1,10000)
        assert r['source_radius']>1 and r['source_sites']==(2*r['source_radius']+1)**3
        assert not r['state_diagonalization_used']


def test_empty_cut_is_exact():
    assert bounds.marginal(F(1),F(1),1,1,0)['half_trace_upper']==0


@pytest.mark.parametrize('args',[(0,1,1,1,1),(1,0,1,1,1),(True,1,1,1,1),(1,1,0,1,1),(1,1,1,0,1),(1,1,1,1,-1),(1,1,1,True,1)])
def test_invalid_or_zero_field_domains_rejected(args):
    with pytest.raises(ValueError):bounds.marginal(*args)


@pytest.mark.parametrize('args',[(1,1,-1,2,1),(1,1,1,1,1),(1,1,1,2,28)])
def test_invalid_box_domains_rejected(args):
    with pytest.raises(ValueError):bounds.box(*args)


@pytest.fixture(scope='module')
def result():return packet.build()


def test_retained_receipt(result):
    assert result==json.loads((packet.HERE/'results.json').read_text())


@pytest.mark.parametrize('key',[k for k,v in packet.claims().items()if not v])
def test_scope_mutation_rejected(result,monkeypatch,key):
    modified=copy.deepcopy(result);modified['claims'][key]=True
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(modified)


@pytest.mark.parametrize('key',['gap','accuracy','binding','type'])
def test_numerical_and_custody_mutation_rejected(result,monkeypatch,key):
    modified=copy.deepcopy(result)
    if key=='gap':modified['source_attachment']['full_gap_lower']='B'
    elif key=='accuracy':modified['box_certificates'][0]['half_trace_upper']='0'
    elif key=='binding':modified['source_bindings'][next(iter(modified['source_bindings']))]='0'*64
    else:modified['all_pass']=1
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(modified)


def test_requested_check_path_is_consumed(tmp_path):
    missing=tmp_path/'missing.json'
    run=subprocess.run([sys.executable,'-m','model.spin_fixed_field_locality','--check',str(missing)],cwd=packet.ROOT,capture_output=True,text=True)
    assert run.returncode and str(missing)in run.stderr


def test_clean_runtime_imports_bound():
    script='''
from pathlib import Path
import sys
from model.spin_fixed_field_locality import packet
packet.build();bindings=packet.bindings();missing=[]
for module in list(sys.modules.values()):
 path=getattr(module,'__file__',None)
 if not path:continue
 path=Path(path).resolve()
 try:key=path.relative_to(packet.ROOT).as_posix()
 except ValueError:continue
 if path.suffix=='.py' and key not in bindings:missing.append(key)
assert not missing,missing
'''
    run=subprocess.run([sys.executable,'-c',script],cwd=packet.ROOT,capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
