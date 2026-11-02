from copy import deepcopy
from fractions import Fraction as F
import math
import numpy as np
import pytest
from . import source,certificate,angles,experiment,packet


def test_marshall_hamiltonian():
    edges,odd=source.cube();h,scale=source.matrix()
    direct=np.zeros((256,256))
    for mask in range(256):
        for x,y in edges:
            if (mask>>x&1)!=(mask>>y&1):
                direct[mask,mask]-=.5;direct[mask^(1<<x)^(1<<y),mask]+=.5
        for x in range(8):direct[mask^(1<<x),mask]-=(-1)**(x in odd)/8
    signs=np.array([(-1)**sum(mask>>x&1 for x in odd) for mask in range(256)])
    np.testing.assert_array_equal(h/scale,signs[:,None]*direct*signs[None,:])


@pytest.mark.parametrize('n,edges',[(1,()),(2,((0,1),)),(3,((0,1),(1,2)))])
def test_perron_gap_and_power_bound(n,edges):
    h,scale=source.matrix(n,edges,F(1),F(1));e,v=np.linalg.eigh(h/scale)
    bound=source.perron_bound(n,edges,F(1),F(1));assert e[1]-e[0]>=float(F(bound['Hamiltonian_gap_lower']))
    state=np.array(source.rational_power(n,edges,F(1),F(1),3),float);state/=np.linalg.norm(state)
    exact=v[:,0]*np.sign(v[:,0].sum());error=np.linalg.norm(state-exact)
    assert error<=2*math.sqrt(2**n)*float(F(bound['block_power_contraction_upper']))**3
    assert np.min(state)>0


@pytest.mark.parametrize('args',[(3,((0,1),(1,2),(0,2)),F(1),F(1)),(2,((1,0),),F(1),F(1)),(1,(),F(1),F(0)),(True,(),F(1),F(1))])
def test_bad_finite_inputs(args):
    with pytest.raises(ValueError):source.perron_bound(*args)


def test_matrix_rejects_boolean_after_integer_call():
    source.matrix(1,(),F(1),F(1))
    with pytest.raises(ValueError):source.matrix(True,(),F(1),F(1))


def test_complete_spectral_certificate():
    result=certificate.verify(experiment.witness())
    lo,hi=map(float,result['ground_energy_interval']);assert lo<-8.2399123998818<hi
    assert F(result['normalized_ground_vector_error'])<F(1,10**12)


@pytest.mark.parametrize('change',['missing','truncated','float','corrupt','negative'])
def test_reject_untrusted_spectral_witness(change):
    w=experiment.witness()
    if change=='missing':w.pop('energies')
    if change=='truncated':w['vectors']=w['vectors'][:,:-1]
    if change=='float':w['vectors']=w['vectors'].astype(float)
    if change=='corrupt':w['vectors'][0,4]+=1<<40
    if change=='negative':w['vectors'][:,0]*=-1
    with pytest.raises(ValueError):certificate.verify(w)


@pytest.mark.parametrize('a,b,bits',[(1,1,40),(3,7,100),(0,1,40),(1,0,80),(1,10**30,100)])
def test_rational_angle_bound(a,b,bits):
    import mpmath
    row=angles.node(a,b,bits);ctx=mpmath.mp.clone();ctx.prec=bits+100
    target=-ctx.atan2(ctx.sqrt(b),ctx.sqrt(a));convert=lambda x:ctx.mpf(x.numerator)/x.denominator
    assert convert(row['lower'])<=target<=convert(row['upper'])
    assert row['error']<=F(1,1<<bits)
    assert abs(convert(row['angle'])-target)<=convert(row['error'])


def test_exact_walsh_roundtrip():
    from model.centered_native_preparation.tree import walsh
    for level in angles.levels((1,2,3,4,5,6,7,8),60):
        assert walsh(level['coefficients'])==[n['angle'] for n in level['nodes']]


@pytest.mark.parametrize('args',[(0,0,10),(-1,2,10),(1,1,0),(True,1,10)])
def test_reject_bad_angles(args):
    with pytest.raises(ValueError):angles.node(*args)


def test_preparation_program_target():
    word,initial,target,checked,row=experiment.plan()
    assert len(word)==1090 and target.shape==(131072,1)
    assert abs(np.linalg.norm(target)-1)<1e-14
    assert .83547<row['pointer_bit_zero_ideal_probability']<.83548
    assert F(row['complete_ground_isometry_error'])<F(5,10**13)


def test_manifest_includes_initializers_and_dynamic_archive():
    pins=packet.source_bindings()
    assert 'model/full_fock_field/__init__.py' in pins
    assert 'model/full_fock_field/bounds.py' in pins
    assert 'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py' in pins
    assert 'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/__init__.py' in pins
    assert 'evidence/native_geometric_source_20260925/spec.json' in pins


def test_native_record_strict_rejection():
    if not packet.OUTPUT.exists():pytest.skip('run cold producer to test retained archive')
    value,arrays=packet.read(packet.OUTPUT)
    bad=deepcopy(value);bad['theorem']['effective_infinite_ground_marginal_modulus']=True
    with pytest.raises(ValueError):packet.validate(bad,arrays)
    bad=deepcopy(value);bad['all_pass']=1
    with pytest.raises(ValueError):packet.validate(bad,arrays)
    bad=deepcopy(arrays);bad['success']=bad['success'][:-1]
    with pytest.raises(ValueError):packet.validate(value,bad)


def test_near_zero_host_diagnostic_is_compared_absolutely():
    packet.compare({'normalization_residual':1.7319479184152442e-14},
                   {'normalization_residual':1.6542323066914832e-14})
    with pytest.raises(ValueError):
        packet.compare({'normalization_residual':0.},
                       {'normalization_residual':3*packet.REPLAY_ATOL})
    with pytest.raises(ValueError,match='negative'):
        packet.compare({'normalization_residual':0.}, {'normalization_residual':-1e-30})


def test_failure_replay_budget_is_global_in_amplitude():
    zero=np.zeros((25,4),dtype=np.float64)
    packet.compare_arrays(zero,np.full_like(zero,(packet.REPLAY_ATOL/20)**2),'failure_weights')
    distributed=np.full_like(zero,(packet.REPLAY_ATOL/5)**2)
    with pytest.raises(ValueError,match='total failure amplitude'):
        packet.compare_arrays(zero,distributed,'failure_weights')
    with pytest.raises(ValueError):
        packet.compare({'complete_first_failure_probability':0.},
                       {'complete_first_failure_probability':float(distributed.sum())})


def test_replay_agreement_does_not_authorize_a_wrong_target():
    value,arrays=packet.read(packet.OUTPUT)
    altered=dict(arrays);altered['success']=arrays['success'].copy()
    altered['success'][0,0]+=1e-4
    with pytest.raises(ValueError,match='certified source target'):
        packet.validate(value,altered)
