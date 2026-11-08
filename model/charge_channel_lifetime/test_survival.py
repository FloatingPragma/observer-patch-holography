"""Independent same-source controls for complete retention and persistence."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from model.protected_defect_catalog import charge
from . import survival, packet


@pytest.fixture(scope='module')
def data():
    return survival.source_data()[0]


@pytest.fixture(scope='module')
def receipt():
    return packet.build()


def test_same_source_width_from_negative_projectors():
    op=charge.operators()
    assert survival.spectral_width_certificate(op)['full_spectral_width_bound']=='14'
    damaged=dict(op);damaged['H']=op['H'].copy();damaged['H'][0,0]+=1
    with pytest.raises(ArithmeticError):survival.spectral_width_certificate(damaged)


def test_local_complete_instrument_factorization_and_charge():
    op=charge.operators()
    survival.instrument_certificate(op)
    for q,a in op['twice_creation'].items():
        assert np.array_equal(op['twice_charge']@a-a@op['twice_charge'],2*q*a)
    bad=dict(op);bad['twice_creation']=dict(op['twice_creation']);bad['twice_creation'][1]=op['local_Z']*2
    with pytest.raises(ArithmeticError):survival.instrument_certificate(bad)


def test_fine_local_X_records_do_not_prepare_pure_branch():
    op=charge.operators();g,_,_=charge.ground_input();g=np.array(g,dtype=float);g/=np.linalg.norm(g)
    for a in op['twice_creation'].values():
        v=a@g;v/=np.linalg.norm(v)
        coefficients=op['walsh'].T@v/16
        dephased_initial_return=float(np.sum(abs(coefficients)**4))
        assert dephased_initial_return<.5
        assert abs(np.dot(v,v)-1)<1e-12


def test_exact_variance_and_two_positive_weights(data):
    assert sum(row['weight'] for row in data.values())==1
    for row in data.values():
        assert row['weight_interval'][0]>F(4,25)
        assert row['variance_interval'][0]>0
        assert len(row['weights'])==256
        assert abs(sum(row['weights'])-1)<=row['gram']


@pytest.mark.parametrize('q',[-1,1])
def test_survival_bounds_enclose_independent_full_eigensystem(q,data):
    op=charge.operators();g,_,_=charge.ground_input();g=np.array(g,dtype=float)
    v=op['twice_creation'][q]@g;v/=np.linalg.norm(v)
    e,u=np.linalg.eigh(op['H']/op['scale']);weights=abs(u.T@v)**2
    for t in (F(),F(1,10),F(1,5),F(2),F(5)):
        lo,hi=survival.survival(data[q],t,32 if t<F(1) else survival.ORDER)
        direct=abs(np.dot(weights,np.exp(-1j*float(t)*e)))**2
        assert float(lo)-2e-14<=direct<=float(hi)+2e-14


@pytest.mark.parametrize('q',[-1,1])
def test_first_threshold_strict_endpoint_bracket(q,data):
    lo,hi=survival.first_threshold(data[q])
    assert 0<lo<hi<F(3,14)
    assert hi-lo<F(1,10**7)
    assert survival.survival(data[q],lo)[0]>F(99,100)
    assert survival.survival(data[q],hi)[1]<F(99,100)
    assert data[q]['variance_interval'][1]*lo**2>F(99,10000)


def test_return_probability_is_original_weight_squared(data):
    op=charge.operators();g,_,_=charge.ground_input();g=np.array(g,dtype=float);g/=np.linalg.norm(g)
    e,u=np.linalg.eigh(op['H']/op['scale'])
    for q in (-1,1):
        a=op['twice_creation'][q]/2;v=a@g;weight=float(v@v)
        for t in (0,.2,2):
            evolved=u@(np.exp(-1j*t*e)*(u.T@v))
            echo=abs(g@op['local_Z']@evolved)**2
            survival_direct=abs(v@evolved/weight)**2
            assert abs(echo-weight**2*survival_direct)<2e-13
        assert abs(weight**2-weight)>.1
        assert float(data[q]['weight_interval'][0])<=weight<=float(data[q]['weight_interval'][1])


def test_fixed_uniform_detector_preserves_charge_contrast():
    op=charge.operators();g,_,_=charge.ground_input();g=np.array(g,dtype=float);g/=np.linalg.norm(g)
    e,u=np.linalg.eigh(op['H']/op['scale']);indices=np.arange(256)
    for q in (-1,1):
        v=op['twice_creation'][q]@g;v/=np.linalg.norm(v)
        for t in (0,.2,2,5):
            evolved=u@(np.exp(-1j*t*e)*(u.T@v))
            means=[np.vdot(evolved,evolved[indices^(1<<x)]).real for x in range(8)]
            ref=[np.dot(g,g[indices^(1<<x)]) for x in range(8)]
            assert abs(sum(s*m for s,m in zip(op['signs'],means))/8-q/4)<2e-13
            assert max(abs((m-r)/2) for m,r in zip(means,ref))>=1/8-1e-13
            assert abs(sum(s*(m-r)/2 for s,m,r in zip(op['signs'],means,ref))-q)<2e-13


def test_complete_echo_and_score_probabilities(receipt):
    rows=receipt['complete_readout']['channels']
    for row in rows:
        for obs in row['observations']:
            rlo,rhi=map(F,obs['original_charge_and_ground_return_interval'])
            nlo,nhi=map(F,obs['same_charge_and_nonreturn_interval'])
            wlo,whi=map(F,row['original_preparation_probability_interval'])
            assert rlo+nlo<=whi and rhi+nhi>=wlo
        score=row['persistent_uniform_site_binary_score']
        assert F(score['original_charge_and_plus_reference_contrast_magnitude_lower'])>F(1,50)
    assert receipt['complete_readout']['persistent_contrast_readout']['complete_charge_setting_and_binary_outcome_count']==32


def test_rebound_and_nonpromotion(receipt):
    result=receipt['complete_readout'];rebound=result['charge_minus_survival_rebound']
    assert F(rebound['5'][0])>F(rebound['2'][1])
    assert all(row['physical_particle_lifetime_or_width'] is None for row in result['channels'])
    assert not result['scope']['full_charge_echo_tape_executed']
    assert not result['scope']['continuum_decay_or_resonance_proved']


@pytest.mark.parametrize('args',[(F(-1),32),(F(1),-1),(1.,32)])
def test_polynomial_rejects_uncontrolled_inputs(args,data):
    with pytest.raises(ValueError):survival.polynomial(data[1]['moments'],*args)


def test_echo_error_preserves_weight_floor():
    bound=survival.echo_error(F(4,25),F(1,10**8),F(1,10**8))
    assert 0<bound<F(1,100000)
    with pytest.raises(ValueError):survival.echo_error(F(),F(1),F(1))


@pytest.mark.parametrize('mutation',['scope','weight','time','binding'])
def test_full_replay_rejects_science_and_binding_mutation(mutation,receipt):
    value=deepcopy(receipt)
    if mutation=='scope':value['complete_readout']['scope']['continuum_decay_or_resonance_proved']=True
    elif mutation=='weight':value['complete_readout']['channels'][0]['complete_spectral_weight_proposal'].pop()
    elif mutation=='time':value['complete_readout']['channels'][0]['first_threshold_time_interval'][0]='0'
    else:value['source_bindings'][survival.WITNESS]='0'*64
    with pytest.raises(ValueError):packet.validate(value)


def test_cli_checks_requested_path(tmp_path,receipt):
    good=tmp_path/'requested.json';good.write_text(packet.canonical(receipt))
    proc=subprocess.run([sys.executable,'-m','model.charge_channel_lifetime','--check',str(good)],cwd=packet.ROOT,capture_output=True,text=True)
    assert proc.returncode==0,proc.stdout+proc.stderr
    good.write_text('{"schema":"forged"}')
    proc=subprocess.run([sys.executable,'-m','model.charge_channel_lifetime','--check',str(good)],cwd=packet.ROOT,capture_output=True,text=True)
    assert proc.returncode!=0


def test_bindings_cover_clean_producer_runtime_imports():
    script='''
import json,sys
from pathlib import Path
from model.charge_channel_lifetime import packet
value=packet.build()
root=packet.ROOT.resolve()
imports=set()
for module in tuple(sys.modules.values()):
    filename=getattr(module,'__file__',None)
    if filename:
        p=Path(filename).resolve()
        if p.suffix=='.py' and p.is_relative_to(root):imports.add(p.relative_to(root).as_posix())
missing=imports-set(value['source_bindings'])
assert not missing,sorted(missing)
assert 'model/protected_defect_catalog/charge.py' in imports
print(json.dumps({'checked':len(imports)}))
'''
    proc=subprocess.run([sys.executable,'-c',script],cwd=packet.ROOT,capture_output=True,text=True)
    assert proc.returncode==0,proc.stdout+proc.stderr
    assert json.loads(proc.stdout)['checked']>10
