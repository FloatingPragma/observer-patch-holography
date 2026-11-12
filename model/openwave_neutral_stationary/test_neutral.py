"""Independent matrix controls and explicit native source mutations."""
import numpy as np
import pytest
import scipy.linalg as la
import sympy as s
from . import isolation,filter,ground,local_profile,records,translations,replay,dynamics,source


def test_complete_isolation_rejects_insufficient_penalty():
    d=isolation.data()
    assert d['rest'][0]>0 and d['width']<s.Rational(4,10**14)
    assert d['separation']>s.Rational(11,1000)
    with pytest.raises(ValueError,match='does not isolate'):isolation.data(1)
    assert d['old_flux_gap']<d['gap'][0] and 2*d['old_flux_gap']<d['gap'][0]


def test_independent_block_self_energy_and_sylvester_controls():
    rng=np.random.default_rng(92815)
    for n in (5,9,14):
        for rank in (1,2):
            for _ in range(4):
                h0=np.diag([0.]*rank+list(1+rng.random(n-rank)))
                raw=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
                v=raw+raw.conj().T;v[:rank,:rank]=0
                v*=.04/la.norm(v,2)
                values,vec=la.eigh(h0+v);norm=la.norm(v,2)
                d=np.min(np.diag(h0)[rank:]);width=norm**2/(d-norm)
                assert np.min(values)>=-width-1e-13 and values[rank-1]<1e-13
                assert values[rank]>=d-norm-1e-13
                e=vec[:,:rank]@vec[:,:rank].conj().T
                p=np.diag([1.]*rank+[0.]*(n-rank))
                assert la.norm(e-p,2)<=norm/d+1e-13
                # Explicit Sylvester equation with the actual low eigenvalues.
                x=vec[rank:,:rank]
                rhs=-(v@vec[:,:rank])[rank:]
                assert la.norm(h0[rank:,rank:]@x-x@np.diag(values[:rank])-rhs)<1e-12
    altered=np.diag([.01,1.])
    assert la.eigvalsh(altered)[0]>0  # nonzero PVP invalidates the a upper bound


def test_native_projector_against_generalized_eigensolver():
    from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import target
    from model.coherent_defect_wave.chiral_joint.local_reference import physical
    from model.coherent_defect_wave.chiral_joint.local_reference.algebra import LEPTONS
    block=target(2,0);h=block['H'];g=block['gram']
    matrix=lambda m:np.array(m.tolist(),dtype=complex)
    energies,vectors=la.eigh(matrix(g*h),matrix(g));v=vectors[:,0]
    result=local_profile.certificate();lo,hi=map(float,isolation.data()['en'])
    assert lo-1e-12<energies[0]<hi+1e-12
    native=matrix(block['basis'])@v;metric=matrix(physical.metric_on(block['keys']))
    assert abs(np.vdot(native,metric@native)-1)<1e-12
    ec_bits=[i for i,n in enumerate(LEPTONS) if n%15==9]
    total=0
    for row in result['complete_gauge_invariant_native_occupation_records']:
        indicator=np.array([2-sum(mask>>i&1 for i in ec_bits)==row['native_N_ec'] and
                             (scalar!=0)==row['scalar_excited'] for _,mask,scalar in block['keys']])
        probability=np.vdot(native,metric@(indicator*native)).real
        assert abs(probability-row['probability_decimal_interval'][0])<2e-12
        total+=probability
    assert abs(total-1)<1e-12


def test_filter_exact_growth_and_native_spectral_replay():
    f=filter.design()
    assert f['epsilon']==s.Rational(1,10**69)
    assert f['degree']==20951 and f['flattening_power']==8
    assert f['pass_loss']**f['flattening_power']<=f['epsilon']
    assert f['conditional_trace_error']<s.Rational(702,10**15)
    assert f['preparation_success_interval'][1]<s.Rational(3,10**126)
    r=replay.certificate()
    assert r['all_native_block_spectral_ranks_retained']==16
    assert sum(z['native_spectral_rank'] for z in r['energy_records'])==16
    assert filter.smallest_power(s.Rational(3,2),s.Rational(27,8))==3
    # Removing full carrier dimension from conditioning understates leakage.
    assert f['outside_band_probability_upper']>f['epsilon']**2*10**120


def test_all_graded_momenta_and_complete_records():
    t=translations.certificate();r=records.certificate()
    assert t['all_eight_momentum_ranks']==[128]*8
    assert r['complete_record_count']==len(r['records'])==51
    assert r['records'][0]['record']=='complete_preparation_failure'
    assert r['records'][2]['record']=='inside_band_complete_atomic_profile_complement'
    assert {z['cell'] for z in r['records'][3:]}==set(range(8))


def test_native_source_penalty_and_unique_ground():
    from model.openwave_neutral_candidate.source import full_column as original
    key=((1<<9)|(1<<45),tuple(1 if j in (0,7) else 0 for j in range(24)),(0,)*8)
    a=source.full_column(key);b=original(key,periodic=True,coupling=isolation.LAMBDA)
    assert a[key]-b.get(key,0)==2
    assert {k:v for k,v in a.items() if k!=key}=={k:v for k,v in b.items() if k!=key}
    g=ground.certificate()
    assert g['actual_kappa_two_ground_unique'] and g['actual_ground_periodic_momentum']==[0,0,0]
    assert g['actual_ground_gap_decimal_lower']>4.45e-25


def test_localization_and_mass_shell_obstruction():
    r=dynamics.certificate()
    assert not r['unit_model_speed_relativistic_dispersion_compatible']
    assert r['localization_rows'][0]['escape_decimal_upper']<3.61e-12
    assert r['localization_rows'][-1]['escape_decimal_upper']<.001585
    assert not any(z['positive_motion_lower_bound_certified'] for z in r['localization_rows'])
    lo,hi=r['calibration']['neutral_keV_decimal_interval']
    assert 66.1638659224<lo<hi<66.1638659225
    assert not r['calibration']['physical_relativistic_rest_mass_identification']
