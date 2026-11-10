"""Independent finite controls for the ordered instrument and all branches."""
from fractions import Fraction as F
import numpy as np
import sympy as sp

from . import model as m


def test_exact_sqrt_rounding_and_actual_native_bounds():
    for x in (F(0),F(2),F(1,10**30),F(11,7),F(9,4)):
        lo,hi=m.sqrt_bound(x),m.sqrt_bound(x,True)
        assert lo**2<=x<=hi**2 and hi-lo<=F(1,2**96)
    g=m.ground_bounds();c=m.clock_bounds()
    assert F(g['actual_admitted_bin_ground_gap_lower'])>F('.124237')
    assert F(g['coherent_native_CAR_spatial_pump_conditional_bin_ground_fidelity_lower'])>F('.9997343')
    assert F(c['coherent_spatial_pump_clock_accepted_infidelity_upper'])<F('2.131e-8')
    assert F(c['coherent_spatial_pump_clock_unconditional_true_acceptance_lower'])>F('.419')
    assert not g['outside_J_full_charged_ground_identified']
    assert c['rational_factor_count_upper']==m.SUBSTEPS*(m.M-1)
    assert c['separate_raw_and_scalar_alternative']['factor_count_upper']==1024*c['rational_factor_count_upper']


def test_exact_spatial_CAR_contraction_on_kinematic_boundaries():
    # Two graded native modes, two even port bits and a ternary flux.
    # Vertex0 is the reference vertex; vertex1 lowers its retained flux.
    states=[(w,p0,p1,ell) for w in range(4) for p0 in range(2)
            for p1 in range(2) for ell in (-1,0,1)]
    lookup={v:i for i,v in enumerate(states)};n=len(states)
    operators=[]
    for y in range(2):
        D=sp.zeros(n)
        for j,(w,p0,p1,ell) in enumerate(states):
            ports=[p0,p1]
            if not w&(1<<y) or ports[y] or (y==1 and ell==-1):continue
            ports[y]=1;sign=(-1)**((w&((1<<y)-1)).bit_count())
            D[lookup[(w^(1<<y),*ports,ell-int(y==1))],j]=sign
        operators.append(D)
    d0,d1=operators
    assert d0*d1+d1*d0==sp.zeros(n)
    assert d0*d1.T+d1.T*d0==sp.zeros(n)
    A=sp.Rational(3,5)*d0+sp.Rational(4,5)*d1
    defect=sp.eye(n)-(A*A.T+A.T*A)
    assert defect.is_diagonal() and all(defect[j,j]>=0 for j in range(n))
    # I-A†A=AA†+defect is a literal positive factorization.
    assert sp.eye(n)-A.T*A==A*A.T+defect
    # Retain both valid and invalid flux inputs; the Gauss difference is
    # unchanged because source charge and flux change together.
    gauss=sp.diag(*[ell-(int(p1))*(-1) for w,p0,p1,ell in states])
    assert d0*gauss==gauss*d0 and d1*gauss==gauss*d1


def test_rank_two_projector_and_small_arrival_transfer_without_eigensolver():
    # A complete analytic two-level model with a spectator G1 doublet.
    # Both projector rotation and nonstationary arrival are nonzero.
    phi=.003;c=np.cos(phi);s=np.sin(phi);rotation=np.array([[c,-s],[s,c]])
    energy=np.array([0.,.2]);Hc=np.diag(energy);H=rotation@Hc@rotation.T
    Pc=np.diag([1.,0.]);P=rotation@Pc@rotation.T
    epsilon=.2*abs(s);theta=epsilon/(.2-epsilon)
    assert np.linalg.norm(P-Pc,2)<=theta+1e-15
    assert np.linalg.norm(np.kron(P,np.eye(2))-np.kron(Pc,np.eye(2)),2)<=theta+1e-15
    initial=np.array([1.,1.])/np.sqrt(2)
    target=np.array([1.,-1.])/np.sqrt(2);Q=np.outer(target,target)
    for t in (.25,1/64):
        Uc=np.diag(np.exp(-1j*t*energy));U=rotation@Uc@rotation.T
        chi0=Q@Uc@initial;chi=Q@U@initial
        assert np.linalg.norm(U-Uc,2)<=epsilon*t+1e-14
        lower=np.linalg.norm(Pc@chi0)-theta*np.linalg.norm(chi0)-epsilon*t
        assert lower>0 and np.linalg.norm(P@chi)>=lower-1e-15
        assert np.linalg.norm(chi)<=np.linalg.norm(chi0)+epsilon*t+1e-15
    # Schur degeneracy can be checked literally for these noncommuting
    # finite spin matrices, without identifying physical continuous spin.
    spinx=np.array([[0,1],[1,0]]);spinz=np.diag([1,-1])
    for S in (spinx,spinz):
        assert np.allclose(np.kron(H,np.eye(2))@np.kron(np.eye(2),S),
                           np.kron(np.eye(2),S)@np.kron(H,np.eye(2)))


def test_ordered_backaction_full_fourier_records_and_rare_weights():
    size=32;B=-37.12345;t0=np.pi/8
    energies=np.array([16*3/size,16*11/size]);Hraw=np.diag(B+energies)
    P=np.diag([1.,0.]);Q=np.full((2,2),.5);initial=np.array([1.,0.])
    chi=Q@initial
    assert abs(np.vdot(chi,chi)-.5)<1e-15
    assert abs(np.vdot(P@chi,P@chi)-.25)<1e-15
    assert np.linalg.norm(P@Q-Q@P)>0
    ks=[]
    for k in range(size):
        K=sum((np.exp(2j*np.pi*j*k/size)*np.diag(np.exp(-1j*j*t0*np.diag(Hraw)))
               *np.exp(1j*j*t0*B) for j in range(size)),np.zeros((2,2),complex))/size
        ks.append(K)
    np.testing.assert_allclose(sum((K.conj().T@K for K in ks)),np.eye(2),atol=3e-14)
    probabilities=np.array([np.linalg.norm(K@chi)**2 for K in ks])
    assert abs(probabilities[3]-.25)<1e-14 and abs(probabilities[11]-.25)<1e-14
    assert abs(probabilities.sum()-.5)<1e-14
    # Common autonomous evolution preserves all pointer probabilities.
    free=np.diag(np.exp(-17j*np.diag(Hraw)))
    np.testing.assert_allclose([np.linalg.norm(free@K@chi)**2 for K in ks],probabilities,atol=1e-15)
    for scale in (1.,1e-6,1e-12):
        np.testing.assert_allclose([np.linalg.norm(K@(scale*chi))**2 for K in ks],scale**2*probabilities,atol=scale**2*1e-15)


def test_actual_raw_clock_alias_control_and_all_outcomes():
    c=m.clock_bounds();center=F(c['target_center_above_explicit_pointer_phase'])
    radius=F(c['ground_radius_including_native_raw_phase_uncertainty'])
    low,high=c['accepted_pointer_first'],c['accepted_pointer_last']
    assert abs(16*F(low,m.M)-center)<=m.HALF_WINDOW
    assert abs(16*F(high,m.M)-center)<=m.HALF_WINDOW
    from model.affine_gauss_charged_clock.model import phase_kernel
    for energy in (center-radius,center+radius,F(c['actual_excited_energy_lower_above_pointer_phase']),F(149,10)):
        probabilities=abs(phase_kernel(float(energy)))**2
        assert abs(probabilities.sum()-1)<1e-12
        accepted=probabilities[low:high+1].sum()
        if energy<center+m.HALF_WINDOW:
            assert 1-accepted<=float(F(c['ground_rejection_upper']))
        else:assert accepted<=float(F(c['excited_acceptance_upper']))


def test_all_preparation_failure_records_and_positive_ordered_signal():
    Q,p,w,a=sp.symbols('Q p w a',nonnegative=True)
    law=[1-Q,Q*(1-p),Q*p/2,Q*p*(sp.Rational(1,2)-w),Q*p*w*(1-a),Q*p*w*a]
    assert sp.expand(sum(law))==1
    coherent=law[:4]+[Q*p*w]
    assert sp.expand(sum(coherent))==1
    rows=m.ordered_arrival_bounds()
    assert len(rows)==4
    for row in rows:
        assert F(row['unconditional_arrival_then_true_ground_clock_acceptance_lower'])>0
        assert F(row['clock_accepted_state_infidelity_upper_after_arrival'])<F('.00436')
        assert F(row['actual_arrival_amplitude_upper'])>F(row['comparison_arrival_amplitude_upper'])
        assert not row['same_time_joint_position_energy_projector_claimed']
    apex=next(r for r in rows if r['arrival_vertex']==0 and r['arrival_time']=='1/4')
    assert F(apex['unconditional_arrival_then_true_ground_clock_acceptance_lower'])>F('4.5027e-6')


def test_public_replay_preserves_research_science_and_parent_source():
    import copy,pytest
    from .certificate import raw_certificate,compare_research,certificate,validate_pins
    value=raw_certificate()
    assert compare_research(value)['all_scientific_fields_equal']
    wrong=copy.deepcopy(value)
    wrong['record_law'][2]='Qp/4'
    with pytest.raises(ValueError):compare_research(wrong)
    public=certificate();validate_pins(public['input_sha256'])
    assert public['whole_bin_source_unchanged']
    assert not public['ground_subspace']['outside_J_full_charged_ground_identified']
    assert public['finite_clock']['clock_applied_directly_to_actual_H']
    assert not public['finite_clock']['long_clock_transferred_by_epsilon_times_duration']
    assert public['capacity']['energy_pointer_dimension']==2**20
    assert public['arrival_readout_service']['duration_in_source_clock_units']==1


def test_arrival_readout_generator_requires_no_commutation_with_source():
    Q=sp.ones(2)/2;X=sp.Matrix([[0,1],[1,0]]);I=sp.eye(2)
    L=sp.pi*sp.kronecker_product(Q,I-X)/2
    U=sp.kronecker_product(I-Q,I)+sp.kronecker_product(Q,X)
    # L has only eigenvalues0,pi by this exact polynomial identity.
    assert L*L==sp.pi*L
    assert U==sp.eye(4)-2*L/sp.pi and U.T*U==sp.eye(4)
    H=sp.kronecker_product(sp.diag(0,1),I)
    assert H*L!=L*H
    # Interaction-picture conjugation yields exp(-iH)U at time1;
    # final common H evolution leaves the stored pointer effects fixed.
    pointer=sp.kronecker_product(I,sp.diag(0,1))
    assert H*pointer==pointer*H
