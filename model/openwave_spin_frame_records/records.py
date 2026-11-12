"""Native cubic rotations, complete tensor/frame records and excess spin."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as s
from scipy.linalg import block_diag
from model.openwave_lepton_mass_shell.spectrum import native_frames,dirac_hamiltonian
from model.lepton_mass_compiler import analyzer
from model.openwave_kinetic_phase_cycle.compiler import annihilator,rounds_for_error,trace_bound

SIG=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
SIN=(0,1,0,-1)

@lru_cache(maxsize=1)
def cube_certificate():
    points=list(product(range(4),repeat=3));where={p:i for i,p in enumerate(points)}
    P=s.zeros(64);X=s.zeros(64)
    for j,(x,y,z)in enumerate(points):
        P[where[(-y%4,x,z)],j]=1;X[where[(x,-y%4,-z%4)],j]=1
    sig=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    K=s.zeros(128)
    for axis in range(3):
        T=s.zeros(64)
        for j,p in enumerate(points):
            target=list(p);target[axis]=(target[axis]+1)%4;T[j,where[tuple(target)]]=1
        K+=s.kronecker_product(sig[axis],s.I*(T-T.H)/2)
    R=s.diag((1-s.I)/s.sqrt(2),(1+s.I)/s.sqrt(2));RX=-s.I*sig[0]
    U=s.kronecker_product(R,P);UX=s.kronecker_product(RX,X)
    assert (U*K*U.H-K).applyfunc(s.simplify)==s.zeros(128)
    assert UX*K*UX.H==K
    assert (U**4+s.eye(128)).applyfunc(s.simplify)==s.zeros(128)
    eps=s.Matrix([[0,1],[-1,0]])
    assert (R*eps*R.T-eps).applyfunc(s.simplify)==s.zeros(2) and RX*eps*RX.T==eps
    # Full four-mode CAR lift, with the pairing source and all number sectors.
    c=[annihilator(4,j)for j in range(4)];N=[x.H*x for x in c]
    SZ=(N[0]-N[1]+N[2]-N[3])/2
    rot=s.diag(*[s.exp(-s.I*s.pi*SZ[j,j]/2)for j in range(16)])
    V=c[0].H*c[3].H-c[1].H*c[2].H;V+=V.H
    assert (rot*V*rot.H-V).applyfunc(s.simplify)==s.zeros(16)
    parity=s.diag(*[(-1)**mask.bit_count()for mask in range(16)])
    assert (rot**4-parity).applyfunc(s.simplify)==s.zeros(16) and rot**8==s.eye(16)
    return {'periodic_sites':64,'exact_spin_spatial_rows':128,'proper_quarter_turn_covariance':True,
        'proper_x_half_turn_covariance':True,'uniform_pairing_covariance':True,
        'full_Fock_control_dimension':16,'quarter_turns_to_source_parity':4,
        'quarter_turns_to_identity':8,'continuous_spatial_rotation_covariance_claimed':False,
        'all_ancillary_flags_must_be_rotated':True}


def occupied_rotation_phase(mass=None):
    points=list(product(range(4),repeat=3));where={p:i for i,p in enumerate(points)}
    dim=2 if mass is None else 4;R=np.diag(np.exp(np.array([-1,1])*1j*np.pi/4))
    if mass is not None:R=np.kron(np.eye(2),R)
    U=np.zeros((64*dim,64*dim),complex);columns=[]
    for j,address in enumerate(points):
        p=np.array([SIN[n]for n in address]);target=(-address[1]%4,address[0],address[2])
        U[dim*where[target]:dim*(where[target]+1),dim*j:dim*(j+1)]=R
        H=-sum((v*A for v,A in zip(p,SIG)),np.zeros((2,2),complex))if mass is None else dirac_hamiltonian(p,mass)
        values,vectors=np.linalg.eigh(H)
        for col in np.flatnonzero(values<-1e-13):
            v=np.zeros(64*dim,complex);v[dim*j:dim*(j+1)]=vectors[:,col];columns.append(v)
    B=np.array(columns).T;restricted=B.conj().T@U@B
    leakage=float(np.linalg.norm(U@B-B@restricted));phase=np.linalg.det(restricted)
    assert leakage<1e-12 and abs(phase-1)<1e-12
    return {'occupied_one_body_modes':len(columns),'rotation_phase_real':float(phase.real),
        'rotation_phase_imaginary':float(phase.imag),'occupied_subspace_leakage':leakage}


def packet(mass,time,axis='x'):
    beta=np.array([1.,1.])/np.sqrt(2)
    spin=np.array([1.,1.])/np.sqrt(2)if axis=='x'else np.array([1.,0.])
    seed=np.kron(beta,spin);momentum=[];energies=[]
    for n in range(4):
        k=SIN[n];H=dirac_hamiltonian((0,0,k),mass);E=np.sqrt(k*k+mass*mass)
        phi=(np.eye(4)+H/E)@seed/4
        momentum.append(np.exp(-1j*E*time)*phi);energies.append(E)
    momentum=np.array(momentum);weight=float(np.vdot(momentum,momentum).real)
    state=np.array([sum((np.exp(2j*np.pi*n*z/4)*phi for n,phi in enumerate(momentum)),np.zeros(4,complex))/2 for z in range(4)])
    assert abs(np.vdot(state,state).real-weight)<1e-14
    return state/np.sqrt(weight),weight,np.array(energies)


def frame_field_records(mass,time,turns):
    psi,w,_=packet(mass,time);theta=turns*np.pi/2
    R=np.kron(np.eye(2),np.diag(np.exp(np.array([-1,1])*1j*theta/2)))
    psi=psi@R.T;p=float(np.vdot(psi[0],psi[0]).real)
    currents=np.array([np.vdot(psi[0],np.kron(np.eye(2),A)@psi[0]).real for A in SIG])
    n=currents/p;assert np.linalg.norm(n-np.array([np.cos(theta),np.sin(theta),0]))<2e-13
    frame=[{'record':'preparation_failure','probability':1-w}]
    for a in range(3):
        for u in (-1,0,1):
            prob=1-p if u==0 else (p+u*currents[a])/2
            frame.append({'axis':a,'outcome':u,'probability':float(w*prob/3)})
    field=[{'record':label,'probability':prob}for label,prob in
           [('both_failed',(1-w)**2),('first_failed',(1-w)*w),('second_failed',w*(1-w))]]
    for a,b in product(range(3),repeat=2):
        for u in (-1,0,1):
            prob=1-p*p if u==0 else (p*p+u*currents[a]*currents[b])/2
            field.append({'axes':[a,b],'product_outcome':u,'probability':float(w*w*prob/9)})
    assert abs(sum(x['probability']for x in frame)-1)<1e-14
    assert abs(sum(x['probability']for x in field)-1)<1e-14
    amplitude=(2*np.exp(-1j*mass*time)+(1+mass/np.sqrt(1+mass*mass))*np.exp(-1j*np.sqrt(1+mass*mass)*time))/4
    assert abs(p-abs(amplitude)**2/w)<1e-14
    position=np.abs(psi)**2
    return {'source_time':time,'physical_quarter_turns':turns,'plane_success_conditional':p,
        'preparation_success':w,'native_plane_spin_current':(currents/2).tolist(),
        'oriented_frame_axis':n.tolist(),'apolar_tensor_field':(np.outer(currents,currents)/4-np.eye(3)*p*p/12).tolist(),
        'complete_frame_records':frame,'complete_two_copy_tensor_records':field,
        'conditional_all_plane_spinor_position_records':position.tolist()}


def pointer_control(weight):
    plus=np.ones((2,2),complex)/2;H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    rho=np.diag(np.diag(plus));dephased=np.diag(H@rho@H).real
    assert np.max(abs(dephased-.5))<1e-15
    rows=[]
    for turns in (4,8):
        source_phase=-1 if turns==4 else 1
        for rotate_reference in (False,True):
            phase=source_phase*((-1 if turns==4 else 1)if rotate_reference else 1)
            vector=H@np.array([1,phase])/np.sqrt(2);p=abs(vector)**2
            rows.append({'quarter_turns':turns,'compensating_odd_reference_rotated':rotate_reference,
                'complete_probabilities':[1-weight,float(weight*p[0]),float(weight*p[1])]})
    return {'source_and_idle_reference_records':rows,
        'pointer_dephasing_control_probabilities':[1-weight,float(weight*dephased[0]),float(weight*dephased[1])],
        'all_quantum_preparation_failures_retained':True}


def full_spin_records(masses):
    distribution=np.array([1.]);origin=0;variance=0.;occupied=0
    for address in product(range(4),repeat=3):
        p=np.array([SIN[n]for n in address]);n=int(p@p);perp=int(p[0]**2+p[1]**2)
        for mass in masses:
            a=perp/(2*(n+mass*mass));distribution=np.convolve(distribution,[a/2,0,1-a,0,a/2]);origin-=2;variance+=a;occupied+=2
        if n:
            up=(1+p[2]/np.sqrt(n))/2
            for _ in range(39):
                distribution=np.convolve(distribution,[1-up,0,up]);origin-=1;variance+=perp/(4*n);occupied+=1
    values=(np.arange(len(distribution))+origin)/2
    assert abs(distribution.sum()-1)<2e-13 and abs(distribution@values)<1e-11
    assert abs(distribution@(values*values)-variance)<1e-10
    return {'all_twice_spin_bins':list(range(origin,origin+len(distribution))),
        'complete_ground_probability_distribution':distribution.tolist(),
        'occupied_native_reference_modes':occupied,'ground_mean':0,'ground_variance':float(variance),
        'spin_up_excitation_distribution':'same complete probabilities on twice-spin bins shifted by+1',
        'excited_source_mean':'1/2','excited_source_variance':float(variance),
        'unfiltered_source_is_a_spin_eigenstate':False,'full_spatial_Noether_angular_momentum_derived':False}


def tv(records1,records2):return sum(abs(a['probability']-b['probability'])for a,b in zip(records1,records2))/2


def certificate():
    native_frames();cube=cube_certificate()
    boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']
    masses=[float(sum(map(F,box),F())/2)for box in boxes]
    rotations=[occupied_rotation_phase(m)for m in masses]+[occupied_rotation_phase()]
    experiments=[];spin=[]
    for k,m in enumerate(masses):
        for t in (0,1,3):
            rows=[frame_field_records(m,t,j)for j in range(9)]
            f=lambda j:rows[j]['complete_frame_records'];q=lambda j:rows[j]['complete_two_copy_tensor_records']
            assert tv(q(0),q(2))<1e-14 and tv(f(0),f(4))<1e-14
            field_change=tv(q(0),q(1));frame_change=tv(f(0),f(2))
            assert field_change>1e-6 and frame_change>1e-4
            experiments.append({'character':k,'time':t,'records':rows,
                'tensor_quarter_turn_TV':field_change,'frame_half_turn_TV':frame_change,
                'pointer':pointer_control(rows[0]['preparation_success'])})
            psi,w,energies=packet(m,t,'z');prob=abs(psi)**2
            spin_density=[float(np.vdot(v,np.kron(np.eye(2),SIG[2]/2)@v).real)for v in psi]
            assert abs(sum(spin_density)-.5)<2e-15
            assert np.linalg.norm(psi@np.kron(np.eye(2),SIG[2]/2).T-psi/2)<1e-14
            spin.append({'character':k,'time':t,'original_failure':1-w,
                'original_all_position_spinor_probabilities':(w*prob).tolist(),
                'conditional_native_excess_spin_density':spin_density,'conditional_integral':'1/2',
                'original_integrated_excess_spin':w/2,
                'complete_energy_probabilities':(np.sum(abs(np.array([((np.eye(4)+dirac_hamiltonian((0,0,SIN[j]),m)/energies[j])/2)@np.kron(np.array([1.,1.])/np.sqrt(2),np.array([1.,0.]))/2 for j in range(4)]))**2,axis=1)).tolist(),
                'energy_bins':energies.tolist()})
    bounds=[]
    for box in boxes:
        hi=F(box[1]);assert hi<1
        bounds.append({'all_time_original_frame_half_turn_TV_lower':str((1-hi)**2/64),
                       'all_time_original_tensor_quarter_turn_TV_lower':str((1-hi)**4/4096)})
    rounds=rounds_for_error(F(1,10**12),time=F(3),sites=64)
    return {'schema':'oph.native_spin_frame_records.v1','native_cube_symmetry':cube,
        'native_occupied_reference_rotation_controls':rotations,'all_time_period_separation_bounds':bounds,
        'complete_source_field_frame_experiments':experiments,'spin_pure_axial_packets':spin,
        'complete_unfiltered_source_spin_law':full_spin_records(masses),
        'finite_compiler_transfer':{'maximum_source_time':3,'sites':64,'cycles':rounds,
            'unhalved_single_source_trace_error':str(trace_bound(rounds,F(3),sites=64)),
            'two_copy_unhalved_trace_error':'at most twice the single-source bound',
            'all_enlarged_source_modes_rotated':True,'source_parity_response_exact_with_leakage':True},
        'scope':'Actual proper cubic rotations on the same full native kinetic/pairing source. An oriented filtered spin-current frame has period2pi; a two-copy apolar tensor field has periodpi; the odd source amplitude returns at4pi relative to an idle pointer/reference. The tensor is a declared two-copy observable, not the entire microscopic state. Canonical excess spin is1/2; no continuous lattice rotation symmetry, full Noether angular momentum, exchange statistics or physical field identity is inferred.',
        'verdict':'PASS_NATIVE_CUBE_FIELD_FRAME_AND_HALF_SPIN_RECORDS'}
