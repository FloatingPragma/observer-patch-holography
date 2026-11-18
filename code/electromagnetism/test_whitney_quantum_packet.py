"""Independent and adversarial controls for the full quantum preparation."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import whitney_quantum_packet as producer
import verify_whitney_quantum_packet as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.load()


def test_current_receipt_fresh_independent_replay(packet):
    result=verifier.verify(packet)
    assert result['accepted'] is True
    assert result['configuration_dimension']==56
    assert result['radiative_dimension']==30
    assert result['initial_projection_norm_squared_lower']=='1/64'
    assert result['quantum_propagation_certified'] is False
    assert result['analytic_domain_proved_by_numeric_replay'] is False
    assert all(min(row['projection_norm_squared_numeric'])>1/64 for row in result['diagnostics'])


def test_producer_rebuild_retains_exact_contract(packet):
    fresh=producer.build()
    assert fresh['source_pins']==packet['source_pins']
    assert fresh['initial_exact']==packet['initial_exact']
    assert fresh['contract']==packet['contract']
    # Do not compare the arbitrary SVD frame across BLAS implementations.
    for a,b in zip(fresh['samples'],packet['samples'],strict=True):
        np.testing.assert_allclose(a['full_q'],b['full_q'],atol=1e-12,rtol=0)
        np.testing.assert_allclose(a['momentum'][30:],b['momentum'][30:],atol=1e-11,rtol=0)


def test_exact_initial_momentum_and_projection_bound(packet):
    exact=packet['initial_exact']
    c=np.array([float(Q(x)) for x in exact['momentum_center_imaginary_in_Qsqrt5']])
    b=np.array([float(Q(x)) for x in exact['momentum_boundary_imaginary_in_Qsqrt5']])
    assert np.max(abs(c+12*b))==0
    p=np.asarray(packet['samples'][0]['momentum'])
    expected=np.r_[np.zeros(43), c@[1,np.sqrt(5)],np.full(12,b@[1,np.sqrt(5)])]
    np.testing.assert_allclose(p,expected,atol=3e-13,rtol=0)
    for row in exact['widths']:
        assert Q(row['A_rational_upper'])<64
        assert Q(row['seed_position_variance'])*Q(row['seed_momentum_variance_hbar1'])==Q(1,4)


def test_time_dependent_gauge_rechart_requires_velocity_term():
    mesh=producer.quantum.geometry(4)
    rng=np.random.default_rng(818)
    q=rng.normal(size=68)/5;v=rng.normal(size=68)/7
    xi=rng.normal(size=13);xi-=xi.mean()
    rate=rng.normal(size=13);rate-=rate.mean()
    psi=q[42:55]+1j*q[55:];pv=v[42:55]+1j*v[55:]
    phase=np.exp(.25j*xi)
    transformed=np.r_[q[:42]+mesh.d@xi,(phase*psi).real,(phase*psi).imag]
    transformed_v=np.r_[v[:42]+mesh.d@rate,(phase*(pv+.25j*rate*psi)).real,(phase*(pv+.25j*rate*psi)).imag]
    a,b=producer.phase_space(q,v,mesh),producer.phase_space(transformed,transformed_v,mesh)
    np.testing.assert_allclose(a['full_q'],b['full_q'],atol=2e-14,rtol=0)
    np.testing.assert_allclose(a['full_velocity'],b['full_velocity'],atol=2e-14,rtol=0)
    np.testing.assert_allclose(a['momentum'],b['momentum'],atol=3e-13,rtol=0)
    # Omitting the gauge-rate term is a physically different tangent.
    wrong=np.r_[v[:42]+mesh.d@rate,(phase*pv).real,(phase*pv).imag]
    bad=producer.phase_space(transformed,wrong,mesh)
    assert np.linalg.norm(bad['full_velocity']-a['full_velocity'])>1e-2


def test_frame_rotation_keeps_physical_preparation(packet):
    mesh=producer.quantum.geometry(4)
    rng=np.random.default_rng(821)
    rotation=np.linalg.qr(rng.normal(size=(30,30)))[0]
    frame=mesh.slice.copy();frame[:42,:30]=frame[:42,:30]@rotation
    changed=replace(mesh,slice=frame)
    source=producer.parent.load()['samples'][1]
    a=producer.phase_space(source['q'],source['velocity'],mesh)
    b=producer.phase_space(source['q'],source['velocity'],changed)
    np.testing.assert_allclose(a['full_q'],b['full_q'],atol=1e-14,rtol=0)
    np.testing.assert_allclose(rotation.T@a['momentum'][:30],b['momentum'][:30],atol=1e-13,rtol=0)
    for sigma in (.25,.5,1):
        aa=producer.overlap_parameters(a['q'],a['momentum'],sigma)
        bb=producer.overlap_parameters(b['q'],b['momentum'],sigma)
        np.testing.assert_allclose(list(aa.values()),list(bb.values()),atol=1e-13,rtol=0)


def test_schur_cotangent_differs_from_naive_principal_block(packet):
    row=packet['samples'][0];mesh=producer.quantum.geometry(4)
    y=np.asarray(row['full_q']);psi=y[42:55]+1j*y[55:]
    metric,_=producer.quantum.coefficients(y[:42],psi,.25,.5,.25,mesh)
    wrong=mesh.slice.T@metric@mesh.slice@np.asarray(row['velocity'])
    assert np.linalg.norm(wrong-np.asarray(row['momentum']))>1e-3


@pytest.mark.parametrize('sigma',[.25,.5,1.])
def test_all_56_seed_modes_and_normalization(sigma):
    q=np.zeros(56);p=np.linspace(-.5,.5,56)
    base=producer.seed_log_half_density(q,q,p,sigma)
    assert base.real==pytest.approx(-14*np.log(2*np.pi*sigma**2))
    for index in range(56):
        point=np.eye(56)[index]*sigma
        change=producer.seed_log_half_density(point,q,p,sigma)-base
        assert change.real==pytest.approx(-.25)
        assert change.imag==pytest.approx(p[index]*sigma)
    assert base.real!=pytest.approx(-28*np.log(2*np.pi*sigma**2))


@pytest.mark.parametrize('sigma',[.25,.5,1.])
def test_neutral_projection_and_radius_against_circle_integral(packet,sigma):
    row=packet['samples'][0];q=np.asarray(row['q']);p=np.asarray(row['momentum'])
    expected=verifier.circle_observables(q,p,sigma)
    actual=producer.overlap_parameters(q,p,sigma)
    assert actual['norm_squared']==pytest.approx(expected['norm_squared'],rel=2e-12)
    assert producer.scalar_radius_moment(q,p,sigma)==pytest.approx(expected['scalar_radius_numeric'],rel=2e-12)
    point=q+np.linspace(-.05,.05,56)
    state=producer.projected_half_density(point,q,p,sigma)
    rotated=producer.projected_half_density(producer.rotate(point,.373),q,p,sigma)
    assert rotated==pytest.approx(state,rel=2e-12,abs=1e-13)
    unprojected=np.exp(producer.seed_log_half_density(point,q,p,sigma))
    assert abs(state-unprojected)>abs(unprojected)*.1
    assert expected['scalar_radius_numeric']!=pytest.approx(26*sigma**2+13,abs=1e-3)


def test_nonzero_charge_projection_discriminant_and_radius():
    q=np.zeros(56);p=q.copy();q[30]=2;p[43]=1
    actual=producer.overlap_parameters(q,p,.5)
    expected=verifier.circle_observables(q,p,.5)
    assert actual['B']==2
    assert actual['norm_squared']==pytest.approx(expected['norm_squared'],rel=2e-12)
    assert producer.scalar_radius_moment(q,p,.5)==pytest.approx(expected['scalar_radius_numeric'],rel=2e-12)
    wrong=np.exp(-actual['A'])*np.i0(np.sqrt(actual['A']**2+actual['B']**2))
    assert abs(wrong-expected['norm_squared'])>1e-2


def test_charge_saturating_overlap_boundary():
    q=np.zeros(56);p=q.copy();q[30]=1;p[43]=2
    value=producer.overlap_parameters(q,p,.5)
    assert value['A']==value['B']==2
    assert value['norm_squared']==pytest.approx(np.exp(-2))
    assert producer.scalar_radius_moment(q,p,.5)==pytest.approx(26*.25)


@pytest.mark.parametrize('path,value',[
    (('contract','configuration_dimension'),5),
    (('contract','configuration_dimension'),56.),
    (('contract','full_56D_preparation'),1),
    (('contract','quantum_propagation_certified'),True),
    (('contract','covariance_evolution_computed'),True),
    (('contract','physical_comparison'),True),
    (('contract','classical_samples_are_quantum_means'),True),
    (('contract','requested_propagation_target','status'),'CERTIFIED'),
    (('contract','requested_propagation_target','residual_upper_bound'),'1/10'),
    (('initial_exact','momentum_center_imaginary_in_Qsqrt5'),['6','2']),
    (('initial_exact','widths',0,'projection_norm_squared_lower'),'1/4'),
    (('initial_exact','widths',0,'seed_position_variance'),'1/4'),
    (('samples',1,'parent_sample_index'),1.),
    (('samples',1,'model_time'),.05),
    (('samples',1,'full_velocity',55),0),
    (('samples',1,'momentum',43),0),
    (('samples',1,'schur_scalar_potential',0),0),
    (('samples',1,'widths',1,'B'),1),
    (('samples',1,'widths',1,'norm_squared'),.5),
    (('samples',1,'widths',1,'scalar_radius_numeric'),5),
    (('orthonormal_coulomb_frame',0,0),1),
    (('samples',),[]),
    (('source_pins',),{}),
])
def test_scientific_and_schema_mutations_fail(packet,path,value):
    changed=deepcopy(packet);item=changed
    for key in path[:-1]:item=item[key]
    item[path[-1]]=value
    with pytest.raises(ValueError):verifier.verify(changed)


def test_source_pins_are_read_fresh(packet,monkeypatch):
    original=Path.read_bytes
    target=verifier.ROOT/'code/electromagnetism/whitney_quantum_packet.py'
    def altered(path):
        data=original(path)
        return data+b'\n' if path==target else data
    monkeypatch.setattr(Path,'read_bytes',altered)
    with pytest.raises(ValueError,match='source pin'):verifier.verify(packet)


@pytest.mark.parametrize('bad',[True,0,-1,float('inf'),float('nan'),1j,'1'])
def test_invalid_widths_fail(bad):
    with pytest.raises((ValueError,TypeError)):
        producer.overlap_parameters(np.zeros(56),np.zeros(56),bad)


@pytest.mark.parametrize('bad',[np.zeros(55),np.zeros(57),np.ones(56,dtype=bool),np.full(56,np.nan),np.zeros(56,dtype=complex)])
def test_invalid_dimensions_and_values_fail(bad):
    with pytest.raises(ValueError):producer.overlap_parameters(bad,np.zeros(56),.5)


@pytest.mark.parametrize('raw',['{"x":1,"x":2}','{"x":NaN}','{"x":1e9999}','{"x":Infinity}'])
def test_strict_json_loader(tmp_path,raw):
    path=tmp_path/'bad.json';path.write_text(raw)
    with pytest.raises(ValueError):verifier.load(path)


def test_standalone_verifier_import_avoids_producer_and_search_path():
    script='''import importlib.util,sys
from pathlib import Path
path=Path(sys.argv[1]);folder=path.parent.resolve()
sys.path[:]=[x for x in sys.path if x and Path(x).resolve()!=folder]
spec=importlib.util.spec_from_file_location('isolated_packet_verifier',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
assert not any(x.endswith('whitney_quantum_packet') and x!='isolated_packet_verifier' for x in sys.modules)
assert not any(x.endswith('whitney_interacting_quantum') for x in sys.modules)
assert module.load()['schema']=='oph.whitney_quantum_packet.v1'
assert module.verify(module.load())['accepted'] is True
assert not any(x.endswith('whitney_interacting_quantum') for x in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',script,str(HERE/'verify_whitney_quantum_packet.py')],capture_output=True,text=True)
    assert result.returncode==0,result.stderr


def test_coulomb_zero_connection_retains_nonzero_electric_field():
    mesh=producer.quantum.geometry(4)
    source=producer.parent.load()['samples'][1]
    mapped=producer.phase_space(source['q'],source['velocity'],mesh)
    electric=-(mapped['full_velocity'][:42]+mesh.d@mapped['transformed_scalar_potential'])
    np.testing.assert_allclose(electric,-np.asarray(source['velocity'])[:42],atol=3e-14,rtol=0)
    assert np.linalg.norm(mapped['full_q'][:42])<1e-13
    assert np.linalg.norm(electric)>.1
    assert np.linalg.norm(electric+mapped['full_velocity'][:42])>.1


def test_full_parent_is_verified_beyond_the_two_prepared_samples(packet,monkeypatch):
    import hashlib
    path=verifier.ROOT/verifier.PARENT_PATH
    raw=verifier.parent.load(path)
    raw['samples'][80]['velocity'][0]+=.01
    forged=producer.canonical(raw)
    old_bytes,old_text=Path.read_bytes,Path.read_text
    def changed_bytes(target):
        return forged if target==path else old_bytes(target)
    def changed_text(target,*args,**kwargs):
        return forged.decode('ascii') if target==path else old_text(target,*args,**kwargs)
    monkeypatch.setattr(Path,'read_bytes',changed_bytes)
    monkeypatch.setattr(Path,'read_text',changed_text)
    changed=deepcopy(packet)
    changed['source_pins'][verifier.PARENT_PATH]=hashlib.sha256(forged).hexdigest()
    with pytest.raises(ValueError):verifier.verify(changed)


@pytest.mark.parametrize('replacement', [0., 1e-12, True, float('nan')])
def test_tiny_reported_charge_is_independently_checked(packet,replacement):
    changed=deepcopy(packet)
    assert changed['samples'][0]['widths'][0]['B'] != 0
    changed['samples'][0]['widths'][0]['B']=replacement
    with pytest.raises(ValueError):verifier.verify(changed)
