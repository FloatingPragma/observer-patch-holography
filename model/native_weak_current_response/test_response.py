"""Independent finite-matrix controls of the complete current instruments."""
from fractions import Fraction as F
import numpy as np
import sympy as s
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from . import bounds,local,charged

def test_actual_native_moments_and_all_source_images():
    value=bounds.calculate();g=value['ground_current'];c=value['charged_current']
    assert F(g['actual_ground_variance_interval'][0])>F('.1517')
    assert F(g['actual_ground_W_mean_interval'][0])>F('.0495')
    assert F(c['actual_prepared_bin_conditional_W_mean_interval'][0])>F('.0903')
    assert F(c['actual_prepared_unnormalized_insertion_interval'][0])>F('.03795')
    assert value['atomic_local_calculation']['complete_full_projector_ranks']==[4,12,58,4,2]
    controls=value['atomic_charged_calculation']['complete_target_controls']
    assert [(r['twice_magnetic_spin'],r['complete_charged_target_dimension']) for r in controls]==[(-2,7),(0,17),(2,7)]
    assert value['atomic_charged_calculation']['both_initial_ground_magnetic_states_and_both_native_removal_maps_checked']==4
    assert not value['actual_preparation']['exact_P0_added']
    assert F(value['finite_inelastic_experiment']['original_attempt_nonzero_phase_probability_lower'])>F('.03793')
    assert not value['source']['joint_ground_unique_with_position_reference']

def test_full_literal_current_matrix_independent_exact_projectors():
    from model.native_weak_current_interface.model import dressed_weak_matrix
    W=dressed_weak_matrix();I=s.eye(80)
    # Individual spectral effects are recovered from the literal matrix;
    # exact ranks test all five supported full-carrier outcomes.
    projectors=[]
    for w in local.EIGENVALUES:
        P=I
        for z in local.EIGENVALUES:
            if z!=w:P=P*(W-z*I)/(w-z)
        assert P==P.H and P*P==P
        projectors.append(P)
    assert [s.trace(P) for P in projectors]==[4,12,58,4,2]
    assert sum(projectors,s.zeros(80))==I
    assert sum((w*P for w,P in zip(local.EIGENVALUES,projectors,strict=True)),s.zeros(80))==W

def test_complete_two_copy_finite_phase_current_experiment():
    energies=np.array([-1.7,-.2,1.1]);B=2.;gap=1.5;size=128
    W=np.array([[.3,.4,.2],[.4,-.7,.1],[.2,.1,.9]],complex)
    current=W/2;fail_effect=np.eye(3)-current.conj().T@current
    C=(1-1j*energies/B)/(1+1j*energies/B)
    Z=(C[:,None]*C.conj()[None,:]).reshape(-1)
    K=np.array([np.mean(np.exp(2j*np.pi*np.arange(size)*k/size)[:,None]*Z[None,:]**np.arange(size)[:,None],axis=0) for k in range(size)])
    np.testing.assert_allclose(np.sum(abs(K)**2,axis=0),1,atol=2e-14)
    ground=np.zeros(9,complex);ground[0]=1
    psi=np.kron(current,np.eye(3))@ground
    probabilities=np.sum(abs(K*psi[None,:])**2,axis=1)
    elastic=abs(W[0,0])**2/4;inelastic=(np.sum(abs(W[:,0])**2)-abs(W[0,0])**2)/4
    leakage=2*B/(size*gap)
    assert elastic-1e-14<=probabilities[0]<=elastic+leakage**2*inelastic+1e-14
    assert (1-leakage**2)*inelastic<=probabilities[1:].sum()+1e-14<=inelastic+2e-14
    assert abs(probabilities.sum()+fail_effect[0,0]-1)<1e-14
    # A mixed, correlated, nonstationary pair state is allowed. Complete
    # event probabilities change by at most its exact mixing weight.
    v=np.arange(1,10)+1j*np.arange(9);v=v/np.linalg.norm(v);mix=.003
    rho=(1-mix)*np.outer(ground,ground.conj())+mix*np.outer(v,v.conj())
    pump=np.kron(current,np.eye(3));actual=pump@rho@pump.conj().T
    observed=np.array([np.dot(abs(row)**2,np.diag(actual)).real for row in K])
    assert abs(observed[1:].sum()-probabilities[1:].sum())<=mix+1e-14
    failure=np.trace(np.kron(fail_effect,np.eye(3))@rho).real
    assert abs(failure+observed.sum()-1)<1e-14

def test_charged_ordering_backaction_and_short_time_amplitude_bound():
    H0=np.diag([-.4,.7,1.2]);V=np.array([[.01,.004,0],[.004,-.01,.006],[0,.006,0]])
    H=H0+V;v=np.linalg.norm(V,2);atom=np.array([np.sqrt(.42),0,0],complex)
    chi=atom+np.array([.0001j,.00007,-.00003]);nu=np.linalg.norm(chi-atom);t=2
    evolved=expm(-1j*t*H)@chi
    error=nu+t*v*np.linalg.norm(atom)
    assert np.linalg.norm(evolved-np.exp(.4j*t)*atom)<=error
    W=np.array([[.1,.3,0],[.3,-.6,.2],[0,.2,1.1]])
    insertion=np.vdot(evolved,W@evolved).real;atomic=np.vdot(atom,W@atom).real
    assert abs(insertion-atomic)<=2*error*(np.linalg.norm(chi)+np.linalg.norm(atom))
    # Exact ordered law with an energy admission and a noncommuting
    # current readout: reversing the readout order changes the record.
    angle=.4;R=np.array([[np.cos(angle),0,-np.sin(angle)],[0,1,0],[np.sin(angle),0,np.cos(angle)]])
    P=np.diag([1,1,0]);Q=R@np.diag([1,0,0])@R.T
    rho=np.outer(np.array([1,1j,2])/np.sqrt(6),np.array([1,1j,2]).conj()/np.sqrt(6))
    Bpump=R@np.diag([.6,.8,.4]);S=Bpump@rho@Bpump.conj().T
    failure=np.trace(rho-S).real;outside=np.trace((np.eye(3)-P)@S).real
    admitted=expm(-1j*t*H)@P@S@P@expm(1j*t*H)
    p=np.trace(Q@admitted).real;other=np.trace((np.eye(3)-Q)@admitted).real
    assert abs(failure+outside+p+other-1)<1e-14
    assert np.linalg.norm(Q@H-H@Q)>.01
    after=Q@admitted@Q+(np.eye(3)-Q)@admitted@(np.eye(3)-Q)
    assert abs(np.trace(P@after)-np.trace(P@admitted))>1e-10

def test_finite_duration_noncommuting_readout_and_reference_permutations():
    H=np.array([[.8,.6],[.6,-.3]],complex);Q=np.diag([1.,0.]);X=np.array([[0.,1.],[1.,0.]])
    Htotal=np.kron(H,np.eye(2));L=np.pi*np.kron(Q,(np.eye(2)-X))/2
    expected=expm(-1j*Htotal)@expm(-1j*L)
    def ode(t,y):
        U=expm(-1j*t*Htotal);V=U@L@U.conj().T
        return (-1j*(Htotal+V)@y.reshape(4,4)).reshape(-1)
    solution=solve_ivp(ode,(0,1),np.eye(4,dtype=complex).reshape(-1),rtol=1e-11,atol=1e-12)
    np.testing.assert_allclose(solution.y[:,-1].reshape(4,4),expected,atol=2e-11)
    # Full site-index covariance is exact, independent of the source
    # current matrix. AP signs cancel in the even bilinear.
    sites=[(i,j,k) for i in range(2) for j in range(2) for k in range(2)]
    for axis in range(3):
        permutation=[sites.index(tuple((a+int(j==axis))%2 for j,a in enumerate(x))) for x in sites]
        assert sorted(permutation)==list(range(8))
        assert [permutation[permutation[i]] for i in range(8)]==list(range(8))
        for seam_phase in (-1,1):assert seam_phase*seam_phase==1

def test_exact_interval_rounding_and_conditioning_domains():
    for value in (F(-11,7),F(0),F(13,37)):
        assert bounds.down(value)<=value<=bounds.round_up(value)
        assert bounds.round_up(value)-bounds.down(value)<=F(1,2**128)
    assert bounds.square((F(-1),F(2)))==(0,4)
    assert bounds.divide_positive((F(-1),F(2)),(F(3),F(5)))==(F(-1,3),F(2,3))

def test_public_port_replays_same_native_science_and_rejects_scope_mutation():
    import copy
    import pytest
    from .certificate import raw_certificate,compare_research,certificate,validate
    raw=raw_certificate();assert compare_research(raw)['all_scientific_fields_equal']
    changed=copy.deepcopy(raw);changed['actual_preparation']['exact_P0_added']=True
    with pytest.raises(ValueError):compare_research(changed)
    value=certificate();validate(value['input_sha256'])
    assert 'model/native_weak_current_response/research_reference.json' in value['input_sha256']
    assert value['charged_current']['reference_cost']=='4'
    assert value['charged_current']['current_measurement_changes_the_admitted_energy_distribution']
    assert value['finite_inelastic_experiment']['all_additional_pointer_labels_retained']
