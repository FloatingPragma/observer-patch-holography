"""Independent source, completeness, analytic-bound and receipt controls."""
import copy
from fractions import Fraction as F
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import bounds, finite, packet


@pytest.fixture(scope='module')
def result():
    return packet.build()


def test_actual_source_complete_instrument_and_finite_boundary(result):
    row=result['actual_finite_local_source_control']
    assert row['complete_complement_effect']=='21/25 I'
    assert row['all_source_spin_rows_retained']==256
    assert len(row['records'])==2
    assert all(F(x['original_probability_interval'][0])>0 for x in row['records'])
    assert all(x['creator_square_exactly_zero'] for x in row['records'])
    assert not row['finite_source_claimed_to_have_continuum_polynomial_Gram_independence']
    from model.protected_defect_catalog.charge import operators
    op=operators();a=op['twice_creation'][1];b=op['twice_creation'][-1]
    assert np.array_equal(a.T@a+b.T@b+21*np.eye(256,dtype=np.int64),25*np.eye(256,dtype=np.int64))
    assert not np.array_equal(a.T@a+b.T@b,25*np.eye(256,dtype=np.int64))
    assert not np.array_equal(op['twice_charge']@a-a@op['twice_charge'],-2*a)


@pytest.mark.parametrize('seminorm',[F(1,10),F(1),F(2),F(100)])
def test_cutoff_uses_full_analytic_tails(seminorm):
    row=bounds.cutoff(seminorm);tail=row['tail_bounds']
    assert tail['tail_L1_upper']<=F(1,4)
    assert tail['tail_L2_squared_upper']<=F(1,16)
    assert row['configuration_variance_upper']<12
    assert F(37,576)==F(1,16)/36+F(1,16)
    assert row['overlap_lower']==1-F(1,16)
    assert row['cutoff_radius']>=1
    assert not row['sampled_profile_claimed_to_certify_seminorm']


def test_positive_local_band_has_same_probe_and_scaling():
    a=bounds.local_band();b=bounds.local_band(F(3,7))
    assert a['each_charge_original_joint_band_probability_lower']==F(9,200000000)
    assert a['commutator_expectation_lower']<=F(1,10)*F(15,16)/(2*125)
    assert b['energy_window']==[F(3,7)*x for x in a['energy_window']]
    assert b['momentum_annulus']==[F(3,14),F(3,7)]
    assert b['each_charge_original_joint_band_probability_lower']==a['each_charge_original_joint_band_probability_lower']
    assert not a['nonlinear_creator_all_momenta_annular']


@pytest.mark.parametrize('smoothness',[2,3,8,12])
def test_compact_filter_endpoint_and_monotonicity_certificates(smoothness):
    coeff=bounds.smoothstep(smoothness)
    assert coeff[0]==0 and sum(coeff)==1
    for k in range(1,smoothness+1):
        d=bounds.derivative(coeff,k)
        assert d[0]==0 and sum(d)==0
    # Its derivative is c*u^K*(1-u)^K, coefficient by coefficient.
    d=bounds.derivative(coeff,1)
    from math import comb
    c=d[smoothness]
    assert d[:smoothness]==(F(0),)*smoothness
    assert d[smoothness:]==tuple(c*(-1)**j*comb(smoothness,j)for j in range(smoothness+1))
    f=bounds.filter_bounds(F(1,2),F(2),smoothness)
    assert f['spectral_transfer_interval']==[F(1,4),F(4)]
    assert f['time_filter_L1_upper']>0
    assert not f['Fourier_samples_used_as_bounds']


def test_filter_clock_and_tail_scaling():
    f=bounds.filter_bounds(F(1,2),F(2),4)
    assert bounds.tail(f,F(2),F(3))==bounds.tail(f,F(1),F(6))
    assert bounds.tail(f,F(2),F(3),derivative_time=True)==2*bounds.tail(f,F(1),F(6),derivative_time=True)
    a=bounds.commutator(f,F(1),F(1000));b=bounds.commutator(f,F(1),F(2000))
    assert b['commutator_norm_upper']<=a['commutator_norm_upper']
    assert b['one_time_derivative_commutator_upper']==a['one_time_derivative_commutator_upper']/16


@pytest.mark.parametrize('n',[1,2,3,5])
def test_positive_transfer_configuration_budget(n):
    band=bounds.local_band();f=bounds.filter_bounds(*band['dimensionless_energy_window'])
    p=band['each_charge_original_joint_band_probability_lower']/f['time_filter_L1_upper']**2
    epsilon=p**n*4**(n-1)/(4*n*(n-1))if n>1 else F(1,2)
    c=bounds.configuration(n,F(1),F(4),f,p,epsilon)
    assert c['near_tuple_integral_error_upper']<=p**n/4
    assert c['commutator_reordering_error_upper']<=p**n/4
    assert bounds.commutator(f,F(1),F(c['source_support_distance_required']))['commutator_norm_upper']<=epsilon
    assert c['complete_created_energy_support'][0]>0
    assert c['original_preparation_record_count']==3**n
    assert c['rational_good_tuple_original_weight_lower']==p**n/4
    assert not c['source_charge_zero_included']


def test_noncommuting_reordering_constant():
    rng=np.random.default_rng(947)
    for n in (2,3,4):
        matrices=[]
        for _ in range(n):
            x=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
            matrices.append(x/(2*np.linalg.norm(x,2)))
        product=np.eye(3,dtype=complex);effects=np.eye(3,dtype=complex)
        for c in matrices:
            product=c@product
            effects=effects@(c.conj().T@c)
        epsilon=max(np.linalg.norm(a@b-b@a,2)for i,x in enumerate(matrices)for y in matrices[i+1:]
                    for a in (x,x.conj().T)for b in (y,y.conj().T))
        assert np.linalg.norm(product.conj().T@product-effects,2)<=n*(n-1)*4**(1-n)*epsilon+1e-15


def test_cook_derivative_identity_is_exact_without_factorization():
    matrices=[np.array([[1,2],[0,1]]),np.array([[1,0],[3,-1]]),np.array([[2,1],[1,0]])]
    derivatives=[np.array([[0,1],[2,0]]),np.array([[1,1],[0,2]]),np.array([[0,0],[1,1]])]
    def mul(items):
        result=np.eye(2,dtype=np.int64)
        for x in items:result=result@x
        return result
    direct=sum(mul(matrices[:i]+[derivatives[i]]+matrices[i+1:])for i in range(3))
    moved=sum(mul(matrices[:i]+matrices[i+1:]+[derivatives[i]])for i in range(3))
    for i in range(3):
        for j in range(i+1,3):
            comm=derivatives[i]@matrices[j]-matrices[j]@derivatives[i]
            moved+=mul(matrices[:i]+matrices[i+1:j]+[comm]+matrices[j+1:])
    assert np.array_equal(direct,moved)
    row=bounds.cook_bound([F(1,2)]*3,[F(1,100)]*3,{(0,1):F(1,1000),(0,2):F(1,2000),(1,2):F(1,1000)})
    assert row['full_product_derivative_norm_upper']==F(7,800)
    with pytest.raises(ValueError):bounds.cook_bound([1,1],[0,0],{})


@pytest.mark.parametrize('call',[
    lambda:bounds.cutoff(1.0),lambda:bounds.local_band(0),lambda:bounds.local_band(1,F(3,4)),
    lambda:bounds.smoothstep(True),lambda:bounds.filter_bounds(2,1),lambda:bounds.integer_root_ceiling(1,0)])
def test_invalid_unproved_numeric_inputs_fail(call):
    with pytest.raises(ValueError):call()


def test_full_receipt_and_scientific_mutation_rejected(result,monkeypatch):
    assert packet.validate(result)
    changed=copy.deepcopy(result)
    changed['analytic_source_conclusions']['incoming_outgoing_waveoperators_constructed']=True
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(changed)
    changed=copy.deepcopy(result)
    changed['positive_equal_sign_multipackets'][1]['complete_created_energy_support'][0]='0'
    with pytest.raises(ValueError):packet.validate(changed)


def test_requested_path_is_used_and_strict_json_rejects_duplicates(tmp_path,result):
    path=tmp_path/'owner.json';path.write_text(packet.canonical(result))
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    cmd=[sys.executable,'-m','model.spin_multipacket_comparison','--check',str(path)]
    assert subprocess.run(cmd,cwd=packet.ROOT,env=env,capture_output=True,text=True).returncode==0
    path.write_text('{"schema":1,"schema":2}')
    assert subprocess.run(cmd,cwd=packet.ROOT,env=env,capture_output=True,text=True).returncode!=0
    path.write_text(packet.canonical({**result,'all_pass':False}))
    assert subprocess.run(cmd,cwd=packet.ROOT,env=env,capture_output=True,text=True).returncode!=0


def test_source_bindings_cover_fresh_actual_build_imports():
    code='''import json,sys\nfrom pathlib import Path\nfrom model.spin_multipacket_comparison import packet\nv=packet.build()\nr=packet.ROOT.resolve()\npaths=[]\nfor m in tuple(sys.modules.values()):\n p=getattr(m,"__file__",None)\n if p:\n  p=Path(p).resolve()\n  if p.suffix==".py" and p.is_relative_to(r):paths.append(p.relative_to(r).as_posix())\nprint(json.dumps({"imports":paths,"bindings":sorted(v["source_bindings"])}))\n'''
    done=subprocess.run([sys.executable,'-c',code],cwd=packet.ROOT,capture_output=True,text=True,
                        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),check=True)
    data=json.loads(done.stdout)
    assert set(data['imports'])<=set(data['bindings'])
    assert 'model/centered_particle_spectrum/source.py'in data['imports']
    assert not set(data['imports'])<=set(data['bindings'])-{'model/centered_particle_spectrum/source.py'}
