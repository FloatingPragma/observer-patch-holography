"""Independent finite-state and complete-pointer centroid controls."""
from fractions import Fraction as F
import numpy as np
import sympy as s
import mpmath as mp
import pytest
from . import bounds as b

def test_actual_native_centroid_and_no_P0_preparation():
    value=b.calculate();lo,hi=map(F,value['finite_clock_conditional_mean_source_excitation_interval'])
    assert F('.6916566792')<lo<hi<F('.6916566827')
    assert hi-lo<F('3.4e-9')
    prep=value['actual_physical_preparation']
    assert not prep['extra_P0_ground_projection_used']
    assert F(prep['original_attempt_bin_success_lower'])>F('.419706')
    assert F(prep['conditional_centroid_preparation_error_upper'])<F('3.012e-12')
    assert value['first_moment_identity']['actual_blank_link_coordinates_replayed']==59
    assert value['clock']['R33_factor_count_upper']==536870400
    assert not value['single_charged_eigenvalue_or_original_Estar_relocalized']
    dlo,dhi=map(F,value['finite_pointer_mean_decoded_source_energy_interval'])
    baseline_error=F(value['clock']['actual_Eg_distance_from_reported_baseline_upper'])
    assert dlo-baseline_error==lo and dhi+baseline_error==hi

def test_complete_flagged_conditioning_bound_for_mixed_complex_states():
    mp.mp.dps=80
    def exact_mpf(x):
        x=s.simplify(x);assert x.is_Rational
        return mp.mpf(int(x.p))/int(x.q)
    def norm1(H):
        assert H==H.H
        a,d=exact_mpf(H[0,0]),exact_mpf(H[1,1]);bb=exact_mpf(s.expand_complex(H[0,1]*H[1,0]))
        return max(abs(a+d),mp.sqrt((a-d)**2+4*bb))
    g=s.Matrix([s.Rational(3,5),s.Rational(4,5)]);tau=g*g.H
    rotation=s.Matrix([[s.Rational(5,13),-s.Rational(12,13)],[s.Rational(12,13),s.Rational(5,13)]])
    checks=0
    for j in range(1,13):
        v=s.Matrix([1,j*s.I]);w=s.Matrix([j+2,1-s.I])
        rho=s.Rational(j,13)*(v*v.H)/(v.H*v)[0]+s.Rational(13-j,13)*(w*w.H)/(w.H*w)[0]
        rho=s.simplify(rho);delta=norm1(rho-tau)/2
        for a,c in ((s.Rational(1,7),s.Rational(5,6)),(s.Rational(4,5),s.Rational(2,3)),(s.Rational(1,100),s.Rational(1,11))):
            K=s.diag(a,c)*rotation;S=s.simplify(K*rho*K.H);T=s.simplify(K*tau*K.H)
            ps,pt=s.trace(S),s.trace(T);p,q=exact_mpf(ps),exact_mpf(pt)
            assert norm1(S-T)+abs(p-q)<=2*delta+mp.mpf('1e-70')
            normalized=s.simplify(S/ps-T/pt)
            assert norm1(normalized)/2<=delta/q+mp.mpf('1e-70')
            # A narrow spectral range can be used after conditioning even
            # when its absolute energy origin is large and negative.
            E=s.diag(-90,-90+s.Rational(1,1000000))
            mean=abs(exact_mpf(s.trace(E*normalized)))
            assert mean<=mp.mpf('1e-6')*delta/q+mp.mpf('1e-70')
            checks+=1
    assert checks==36

def test_zero_first_moment_and_complete_tail_centroid_bound():
    psi=s.Matrix([s.Rational(1,3),s.Rational(2,3),s.Rational(2,3)])
    H=s.diag(2,-s.Rational(1,2),0);residual=H*psi
    assert (psi.T*psi)[0]==1 and (psi.T*H*psi)[0]==0
    delta=residual*psi.T+psi*residual.T;v2=(residual.T*residual)[0]
    assert v2==s.Rational(5,9) and delta**3==v2*delta
    assert (H-delta)*psi==s.zeros(3,1)
    P=s.diag(0,1,1);a=P*psi;other=s.Matrix([0,s.Rational(1,100),-s.Rational(1,200)])
    chi=a+other;norm=(chi.T*chi)[0];moment=abs((chi.T*H*chi)[0])
    selected=abs((a.T*H*a)[0]);assert selected<=v2/2
    rho=b.resolution.root(F(1,100**2)+F(1,200**2),True)
    upper=F(5,18)+2*rho+rho*rho
    assert F(moment)<=upper
    lower=F(89,100);assert F(norm)>lower
    assert F(moment/norm)<=upper/lower

def test_uniform_linear_phase_mean_bound_and_wraparound():
    from model.affine_gauss_charged_clock.model import phase_kernel
    for size in (16,64,256,1024):
        grid=np.arange(size)/size
        for theta in (.25,.25001,1/3,.49999,.5,.74999,.75):
            probabilities=abs(phase_kernel(16*theta,size))**2
            assert abs(probabilities.sum()-1)<1e-13
            actual=np.dot(probabilities,abs(grid-theta))
            assert actual<=float(b.mean_grid_error(size))+1e-14
    size=b.M;grid=np.arange(size)/size
    for theta in (.2500012345,.73123456789):
        probabilities=abs(phase_kernel(16*theta,size))**2
        assert abs(probabilities.sum()-1)<1e-12
        assert np.dot(probabilities,abs(grid-theta))<=float(b.mean_grid_error(size))+1e-13
    for bad in (0,2,17,64.0):
        with pytest.raises(ValueError):b.mean_grid_error(bad)

def test_literal_mixed_density_pointer_mean_and_raw_phase():
    size=32;j=np.arange(size);B=-85.42;r=.0003
    theta=np.array([.27,.51,.73]);energy=B+4*r*theta
    vector=np.array([1.,1j,2.])/np.sqrt(6)
    rho=.7*np.outer(vector,vector.conj())+.3*np.diag([.2,.3,.5])
    t0=np.pi/(2*r);kraus=[]
    for k in range(size):
        # Evaluate the explicit shifted extra interaction, keeping B in
        # the raw-energy definition; no source ground energy is inserted.
        phases=np.exp(-1j*np.outer(j,t0*(energy-B)))
        kraus.append(np.diag(np.exp(2j*np.pi*j*k/size)@phases/size))
    effects=sum(K.conj().T@K for K in kraus)
    np.testing.assert_allclose(effects,np.eye(3),atol=1e-13)
    probabilities=np.array([np.trace(K@rho@K.conj().T).real for K in kraus])
    assert abs(probabilities.sum()-1)<1e-13
    decoded=B+4*r*np.arange(size)/size
    centroid=np.trace(np.diag(energy)@rho).real
    assert abs(np.dot(decoded,probabilities)-centroid)<=4*r*float(b.mean_grid_error(size))+1e-11
    # Actual density need not commute with H; the mean law remains valid.
    assert np.linalg.norm(np.diag(energy)@rho-rho@np.diag(energy))>0

def test_rational_pointer_error_on_actual_phase_domain():
    mp.mp.dps=80;value=b.calculate();gamma=F(value['clock']['R33_controlled_error_upper'])
    limit=mp.mpf(gamma.numerator)/gamma.denominator
    for theta in (mp.mpf(1)/4,mp.mpf(1)/3,mp.mpf(1)/2,mp.mpf(3)/4):
        x=2*mp.pi*theta/b.SUBSTEPS;a=1-x*x/10;d=x/2-x**3/120
        U=((a-1j*d)/(a+1j*d))**b.SUBSTEPS
        error=abs(U**(b.M-1)-mp.exp(-2j*mp.pi*theta*(b.M-1)))
        assert error<=limit
    q,f,o,inside=s.symbols('q f o inside')
    assert s.expand(1-q+q*f+q*o+q*(1-f-o))==1

def test_public_replay_preserves_actual_preparation_and_centroid_scope():
    import copy
    from .certificate import raw_certificate,compare_research,certificate,validate
    raw=raw_certificate();assert compare_research(raw)['all_scientific_fields_equal']
    wrong=copy.deepcopy(raw);wrong['actual_physical_preparation']['extra_P0_ground_projection_used']=True
    with pytest.raises(ValueError):compare_research(wrong)
    value=certificate();validate(value['input_sha256'])
    assert not value['construction_scope']['outer_affine_coupled_source_used']
    assert not value['construction_scope']['actual_prepared_source_state_asserted_exact_or_stationary']
    assert value['clock']['prior_bin_admission_remains_explicit_spectral_service']
    assert 'model/native_charge_centroid/research_reference.json' in value['input_sha256']
