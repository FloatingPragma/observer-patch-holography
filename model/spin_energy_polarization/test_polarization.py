import copy
from fractions import Fraction as F
import json
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import algebra,bounds,packet,samples


def test_exact_source_edge_checks_and_current_remainder():
    result=algebra.exact_checks()
    assert len(result['weighted_current_identities'])==3
    assert all(x['exact_weighted_current_commutator']for x in result['weighted_current_identities'])
    assert result['inherited_source_twist_algebra']['exact_single_axis_rotation_energy']
    # Generic norm control genuinely fails to give the missing vanishing limit.
    rows=[bounds.scalar_field_remainder(F(1,n),F(2),F(8))for n in (4,16,256)]
    assert all(r['Q_current_remainder_operator_norm_upper']>128 for r in rows)
    assert all(not r['this_bound_vanishes_as_spacing_tends_to_zero']for r in rows)


def test_weighted_source_rotation_against_independent_finite_matrices():
    pauli=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    left=[np.kron(x/2,np.eye(2))for x in pauli]
    right=[np.kron(np.eye(2),x/2)for x in pauli]
    h=sum(a@b for a,b in zip(left,right))-np.eye(4)/4
    ground=np.array([0.,1.,-1.,0.])/np.sqrt(2)
    for axis in (1,2):
        for weight in (-2.,0.,1.25):
            u=expm(1j*(.31*left[axis]-.17*right[axis]))
            measured=np.vdot(u@ground,weight*h@(u@ground))-np.vdot(ground,weight*h@ground)
            expected=weight*.5*(1-np.cos(.48))
            assert abs(measured-expected)<2e-15
    # This is an independent finite source algebra control; its singlet is
    # never substituted for the ordered infinite-volume reference.


@pytest.mark.parametrize('channel',['Q','P'])
def test_rational_mesh_bounds_enclose_direct_actual_cosine_sum(channel):
    n=4;p=1/n;row=samples.control(channel,n)
    grid=np.arange(-n-2,3*n+3)*p
    x=np.array(np.meshgrid(grid,grid,grid,indexing='ij')).reshape(3,-1).T
    def field(v):
        z=v-np.array([1.,0.,0.])
        return np.prod(np.maximum(0.,1-z*z)**4,axis=1)
    f=field(x);actual=0.;quadratic=0.
    for direction in range(3):
        step=np.eye(3)[direction]*p;other=field(x+step)
        delta=p*p*(f+other)if channel=='Q'else p*(f-other)
        coefficient=(x[:,0]+step[0]/2)/p
        actual+=np.sum(coefficient*(1-np.cos(delta)))
        quadratic+=np.sum(coefficient*delta*delta/2)
    assert abs(quadratic-float(row['quadratic_polarization_per_kappa']))<2e-14
    assert float(row['exact_source_polarization_lower_per_kappa'])<=actual<=float(row['exact_source_polarization_upper_per_kappa'])
    generic=bounds.mesh(channel,F(1,n),F(2),F(1),F(8),F(64),F(512))
    assert abs(actual-float(row['continuum_polarization_per_kappa']))<float(generic['total_error_per_kappa'])


def test_exact_actual_even_translation_law_and_declining_mesh_error():
    for channel in ('Q','P'):
        errors=[]
        for n in (2,4,8,16):
            zero=samples.control(channel,n,(F(0),)*3)
            shifted=samples.control(channel,n)
            assert zero['quadratic_polarization_per_kappa']==0
            assert shifted['quadratic_energy_per_kappa']==zero['quadratic_energy_per_kappa']
            assert shifted['quadratic_polarization_per_kappa']==zero['quadratic_energy_per_kappa']
            errors.append(abs(shifted['quadratic_polarization_per_kappa']-shifted['continuum_polarization_per_kappa']))
        assert all(a>b for a,b in zip(errors,errors[1:]))
    with pytest.raises(ValueError):samples.control('Q',2,(F(1,2),F(0),F(0)))


def test_symmetric_source_kick_bound_decays_and_keeps_its_domain():
    rows=[bounds.mixed_kick(F(1,n),F(2),F(1),F(8),F(64),F(1),F(8),F(64),F(512))for n in (4,16,256)]
    assert all(a['nonlinear_rotation_absolute_error']>b['nonlinear_rotation_absolute_error']for a,b in zip(rows,rows[1:]))
    assert all(r['mixed_quadratic_term_vanishes_by_actual_bond_covariance']for r in rows)
    with pytest.raises(ValueError):bounds.mixed_kick(F(1),F(2),F(1),F(8),F(64),F(1),F(8),F(64),F(512))


def test_complete_configuration_instrument_bound_keeps_zero_record():
    # Exact scalar joint spectral points supply independent control of the
    # square-root complement, including near the maximum creator amplitude.
    for z,w in ((.5+0j,.49+.01j),(.2+.3j,-.1+.4j),(0j,.5j)):
        a=np.array([z,z.conjugate(),np.sqrt(1-2*abs(z)**2)])
        b=np.array([w,w.conjugate(),np.sqrt(1-2*abs(w)**2)])
        assert abs(np.vdot(a,a)-1)<1e-15
        assert np.linalg.norm(a-b)**2<=4*abs(z-w)**2+1e-15
    rows=[bounds.charged_truncation(F(n),6,F(1))for n in (2,4,16)]
    assert all(a['complete_isometry_error_squared_upper']>b['complete_isometry_error_squared_upper']for a,b in zip(rows,rows[1:]))
    assert all(r['complement_retained']and not r['operator_norm_error_claimed']for r in rows)


def test_finite_complete_energy_score_handles_constant_only_source():
    constant=bounds.finite_readout((),F(-3),F(1),F(1))
    assert constant['assign_constant_to_every_record']and constant['exact_increment']==0
    row=bounds.finite_readout((F(1,4),F(-1,8)),F(-3),F(1,1000),F(1,2000))
    assert row['absolute_score_upper']==F(7,2)
    assert row['increment_error_upper']==F(21,2000)
    with pytest.raises(ValueError):bounds.finite_readout((.25,),F(0),F(0),F(0))


@pytest.mark.parametrize('channel',['Q','P'])
def test_Schwartz_tail_decay_and_bound_domains(channel):
    rows=[bounds.schwartz_tail(channel,F(n),6,F(1))for n in (2,4,16)]
    assert all(a['weighted_tail_per_kappa']>b['weighted_tail_per_kappa']for a,b in zip(rows,rows[1:]))
    for bad in (True,1.0,F(-1)):
        with pytest.raises(ValueError):bounds.mesh(channel,bad,F(2),F(1),F(8),F(64),F(512))
    with pytest.raises(ValueError):bounds.schwartz_tail(channel,F(2),2,F(1))
    with pytest.raises(ValueError):bounds.charged_truncation(F(2),3,F(1))


@pytest.fixture(scope='module')
def result():return packet.build()


def test_retained_exact_receipt(result):
    assert result==json.loads((packet.HERE/'results.json').read_text())


@pytest.mark.parametrize('key',[k for k,v in packet.claims().items()if not v])
def test_unproved_physical_attachment_cannot_be_promoted(result,monkeypatch,key):
    value=copy.deepcopy(result);value['claims'][key]=True
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(value)


@pytest.mark.parametrize('part',['bound','record','source','type'])
def test_mutated_exact_control_is_rejected(result,monkeypatch,part):
    value=copy.deepcopy(result)
    if part=='bound':value['general_mesh_bounds_per_kappa'][0]['total_error_per_kappa']='0'
    if part=='record':value['complete_charge_truncation_examples'][0]['complement_retained']=False
    if part=='source':value['source_bindings'][next(iter(value['source_bindings']))]='0'*64
    if part=='type':value['all_pass']=1
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(value)


def test_actual_PATH_is_consumed(tmp_path):
    missing=tmp_path/'missing.json'
    run=subprocess.run([sys.executable,'-m','model.spin_energy_polarization','--check',str(missing)],cwd=packet.ROOT,capture_output=True,text=True)
    assert run.returncode and str(missing)in run.stderr


def test_clean_full_producer_imports_are_bound():
    code='''
from pathlib import Path
import sys
from model.spin_energy_polarization import packet
packet.build(); bound=packet.bindings(); missing=[]
for module in list(sys.modules.values()):
    value=getattr(module,'__file__',None)
    if not value:continue
    p=Path(value).resolve()
    if p.suffix=='.py' and p.is_relative_to(packet.ROOT) and p.relative_to(packet.ROOT).as_posix()not in bound:missing.append(str(p))
assert not missing,missing
'''
    run=subprocess.run([sys.executable,'-W','error','-c',code],cwd=packet.ROOT,capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
