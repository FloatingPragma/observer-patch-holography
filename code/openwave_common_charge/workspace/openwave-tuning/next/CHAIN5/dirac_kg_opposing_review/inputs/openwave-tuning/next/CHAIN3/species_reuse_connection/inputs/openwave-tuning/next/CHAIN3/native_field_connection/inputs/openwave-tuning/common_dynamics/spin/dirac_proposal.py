"""An explicit spin-dependent replacement of the scalar local walk.

This file declares new dynamics. It does not identify them with W=S R.
All spin modes are co-located at port zero; routed shifts carry Gauss current.
"""
from collections import defaultdict, deque
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent/'vector_field')]
import spin
import autonomous as field
import mode_field

MODEL_ID = 'native-local12-dirac-routed-v3-proposed'
I = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]],complex)
Y = np.array([[0,-1j],[1j,0]],complex)
Z = np.diag([1.,-1.]).astype(complex)
PAULI = dict(I=I,X=X,Y=Y,Z=Z)


def pauli(label):
    """The source compiler uses little-endian tensor labels."""
    answer = np.ones((1,1),complex)
    for letter in label:
        answer = np.kron(PAULI[letter],answer)
    return answer


ALPHA = tuple(pauli(label) for label in ('XX','XY','XZ'))
BETA = pauli('ZI')
SIGMA = tuple(pauli('I'+letter) for letter in 'XYZ')


def diagonalizer_word(axis):
    label = ('XX','XY','XZ')[axis]
    basis=[]
    for j, letter in enumerate(label):
        if letter=='X':
            basis += [spin.compiler.pulse('Y',-math.pi/4,[j])]
        elif letter=='Y':
            basis += [spin.compiler.pulse('X',math.pi/4,[j])]
    return spin.compiler.inverse(basis)+spin.compiler.cnot(0,1)


def ideal_word(word):
    answer=np.eye(4,dtype=complex)
    for event in word:
        targets=event['targets']
        kind=event['kind']
        angle=event['angle']
        # These diagonalizers use H,J,X,Y,Z native Hamiltonian labels.
        if kind=='J':
            projector=np.eye(4,dtype=complex); projector[3,3]=0
            gate=expm(1j*angle*projector)
        else:
            label=['I','I']
            label[targets[0]]=kind if kind in 'XYZI' else 'I'
            generator=(X+Z)/np.sqrt(2) if kind=='H' else PAULI[kind]
            gate=expm(1j*angle*(np.kron(I,generator) if targets[0]==0 else np.kron(generator,I)))
        answer=gate@answer
    return answer


@lru_cache(maxsize=None)
def diagonalizer(axis):
    answer=ideal_word(diagonalizer_word(axis))
    assert np.max(abs(answer@ALPHA[axis]@answer.conj().T-pauli('IZ'))) < 2e-14
    return answer


def mass_word(h,m):
    return spin.compiler.pauli_rotation('ZI',-float(h)*float(m))


def bloch(p,h,m):
    answer=expm(-1j*h*m*BETA)
    for a in range(3):
        answer=expm(-1j*h*float(p[a])*ALPHA[a])@answer
    return answer


def continuum_bound(p,h,m,n):
    """Operator-norm Lie-product bound; terminal TV is capped at one."""
    if h<=0 or type(n) is not int or n<0:
        raise ValueError('positive h and nonnegative integer macro count required')
    values=[abs(float(m))]+list(abs(np.asarray(p,dtype=float)))
    return min(2.,n*h*h*sum(values[a]*values[b] for a in range(4) for b in range(a+1,4)))


@lru_cache(maxsize=None)
def route(cell,axis,direction=1):
    """Actual carrier-seam path; negative travel reverses the positive path."""
    g,B,A,labels=field.data()
    coords={c:i for i,c in enumerate(g.cell_coordinates)}
    c=list(g.cell_coordinates[cell]); c[axis]=(c[axis]+direction)%2
    other=coords[tuple(c)]
    if direction==-1:
        destination,path,current=route(other,axis,1)
        assert destination==cell
        return other,tuple(reversed(path)),-current
    a,b=spin.common.spatial.links()[axis][0]
    def internal(start,end):
        queue=deque([(start,(start,))]); seen={start}
        while queue:
            x,path=queue.popleft()
            if x==end:return path
            for y in g.neighbours[x]:
                if y//12==start//12 and y not in seen:
                    queue.append((y,path+(y,)));seen.add(y)
        raise AssertionError('connected native cell required')
    path=internal(12*cell,12*cell+a)+internal(12*other+b,12*other)
    edges={tuple(edge):j for j,edge in enumerate(g.seams)}
    current=np.zeros(288,dtype=np.int64)
    for x,y in zip(path,path[1:]):
        edge=tuple(sorted((x,y)));j=edges[edge]
        current[j]+=1 if x==edge[0] else -1
    target=np.zeros(96,dtype=np.int64);target[12*cell]=1;target[12*other]-=1
    assert np.array_equal(B@current,target)
    return other,path,current


def initial(cell=0,spinor=None,charge=1):
    if charge not in (-1,1):raise ValueError('declared charge sign required')
    if spinor is None:spinor=np.array([1,0,0,0],complex)
    spinor=np.asarray(spinor,dtype=complex)
    if spinor.shape!=(4,) or abs(np.vdot(spinor,spinor)-1)>1e-12:
        raise ValueError('normalized four-component source preparation required')
    F=charge*field.current_walk.geometry()['P'][:,12*cell]
    return {(cell,bytes((F%16).astype(np.uint8))):spinor.copy()}


def axis_step(state,axis,charge=1):
    result=defaultdict(lambda:np.zeros(4,dtype=complex))
    A=diagonalizer(axis)
    for (cell,F),values in state.items():
        rail=A@values
        for branch in range(4):
            direction=1 if branch<2 else -1
            other,path,current=route(cell,axis,direction)
            result[other,field.add_field(F,charge*current)][branch]+=rail[branch]
    return {key:A.conj().T@value for key,value in result.items()}


def matter_step(state,h,m,charge=1):
    result={key:expm(-1j*h*m*BETA)@value for key,value in state.items()}
    for axis in range(3):result=axis_step(result,axis,charge)
    return result


def plaquette_step(state,beta,face=0):
    """The identical finite Z16 plaquette kick on four-spinor coefficients."""
    a=field.data()[2][:,face];pivot=int(np.flatnonzero(a)[0]);sign=int(a[pivot]);q=16
    fourier=np.exp(2j*np.pi*np.outer(np.arange(q),np.arange(q))/q)/np.sqrt(q)
    kick=(fourier*np.exp(1j*beta*np.cos(2*np.pi*np.arange(q)/q))[None,:])@fourier.conj().T
    orbits=defaultdict(lambda:np.zeros((q,4),complex))
    for (cell,F),values in state.items():
        coordinate=(sign*F[pivot])%q
        canonical=field.add_field(F,-coordinate*a)
        orbits[cell,canonical][coordinate]+=values
    result={}
    for (cell,F),values in orbits.items():
        output=kick@values
        for t in range(q):result[cell,field.add_field(F,t*a)]=output[t]
    return result


def step(state,h,m,alpha=0.,beta=0.,charge=1,face=0,field_substeps=1):
    result=matter_step(state,h,m,charge)
    for _ in range(field_substeps):
        result=field.electric_step(plaquette_step(result,beta,face),alpha)
    return result


def records(state,charge=1):
    B=field.data()[1];reference=np.eye(96,dtype=np.int64)[:,0]*charge
    norm=0.;positions=np.zeros(8);energy=0.;spinmean=np.zeros(3)
    for (cell,F),value in state.items():
        n=np.zeros(96,dtype=np.int64);n[12*cell]=charge
        assert not np.any((n+B@np.frombuffer(F,dtype=np.uint8).astype(np.int64)-reference)%16)
        weight=float(np.vdot(value,value).real);norm+=weight;positions[cell]+=weight
        energy+=weight*field.electric_energy(F)
        spinmean += [float(np.vdot(value,s@value).real)/2 for s in SIGMA]
    assert abs(norm-1)<2e-10
    return dict(norm_squared=norm,root_cell_probabilities=positions.tolist(),electric_energy=energy,
                spin_expectation=spinmean.tolist(),basis_blocks=len(state),modular_Gauss=True,
                discarded_probability=0.,success_normalization=False)


def background_matrix(axis,edge_angles,charge=1):
    """Exact character action of the same path-dependent field translations."""
    A=diagonalizer(axis);answer=np.zeros((32,32),complex)
    for cell in range(8):
        for branch in range(4):
            other,path,current=route(cell,axis,1 if branch<2 else -1)
            projector=np.outer(A.conj().T[:,branch],A[branch])
            answer[4*other:4*other+4,4*cell:4*cell+4]+=np.exp(-1j*charge*np.dot(edge_angles,current))*projector
    return answer


def run():
    paths=[Path(__file__),Path(spin.__file__),Path(field.__file__),Path(mode_field.__file__),
           HERE.parent/'continuum/MANIFEST.json',HERE.parent/'manybody/MANIFEST.json',
           HERE.parent/'clock/MANIFEST.json',HERE.parent/'energy/MANIFEST.json']
    pins={str(p.relative_to(spin.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    eye=np.eye(4)
    generators=list(ALPHA)+[BETA]
    clifford=max(np.max(abs(a@b+b@a-(2*eye if i==j else 0)))
                  for i,a in enumerate(generators) for j,b in enumerate(generators))
    assert clifford==0
    assert np.max(abs(sum((s/2)@(s/2) for s in SIGMA)-.75*eye))==0
    native=[]
    for axis in range(3):
        word=diagonalizer_word(axis)
        executed=spin.recovery.execute_complete(word,2,np.eye(4,dtype=complex))
        error=float(np.linalg.norm(executed.success-diagonalizer(axis),ord=2))
        assert error<2e-11
        native.append(dict(axis=axis,word=word,native_pulses=len(word),operator_error=error,
                           export_failure_norm=float(np.linalg.norm(executed.failure_gram,ord=2))))
    word2=spin.compiler.pauli_rotation('IZ',-math.pi)
    word4=spin.compiler.pauli_rotation('IZ',-2*math.pi)
    two=spin.native_controlled_records(word2,np.array([1,0,0,0],complex),2)
    four=spin.native_controlled_records(word4,np.array([1,0,0,0],complex),2)
    for law in ('C','E'):
        assert abs(two['laws'][law]['bits'][1]-1)<2e-11
        assert abs(four['laws'][law]['bits'][0]-1)<2e-11
    covariance=[]
    for axis in range(3):
        angle=.37;U=expm(-.5j*angle*SIGMA[axis]);p=np.array([.19,-.27,.31])
        J=np.zeros((3,3));a=(axis+1)%3;b=(axis+2)%3;J[a,b]=-1;J[b,a]=1
        O=expm(angle*J);H=sum(p[a]*ALPHA[a] for a in range(3))+.43*BETA
        target=sum((O@p)[a]*ALPHA[a] for a in range(3))+.43*BETA
        error=float(np.linalg.norm(U@H@U.conj().T-target));assert error<1e-13;covariance.append(error)
    bounds=[]
    p=np.array([.3,-.2,.1]);m=.4
    H=sum(p[a]*ALPHA[a] for a in range(3))+m*BETA
    for h,n in ((.2,5),(.1,10),(.05,20)):
        error=float(np.linalg.norm(np.linalg.matrix_power(bloch(p,h,m),n)-expm(-1j*n*h*H),ord=2))
        bound=continuum_bound(p,h,m,n);assert error<=bound+1e-12
        bounds.append(dict(h=h,n=n,error=error,operator_bound=bound,terminal_TV_bound=min(1.,bound)))
    doublers=[]
    for corner in product((0.,math.pi),repeat=3):
        sign=(-1)**sum(x!=0 for x in corner)
        assert np.max(abs(bloch(corner,1.,0.)-sign*eye))<1e-14
        doublers.append(dict(momentum=list(corner),Floquet_sign=sign))
    routes=[]
    for cell in range(8):
        for axis in range(3):
            other,path,current=route(cell,axis)
            routes.append(dict(cell=cell,axis=axis,destination=other,vertices=path,
                               edge_current=[(int(e),int(current[e])) for e in np.flatnonzero(current)]))
    B=field.data()[1];rng=np.random.default_rng(3110)
    # Only characters of finite Z16 are physical here: theta=2*pi*z/16.
    theta=2*np.pi*rng.integers(0,16,size=288)/16
    chi=2*np.pi*rng.integers(0,16,size=96)/16
    G=np.diag(np.repeat(np.exp(1j*chi[::12]),4))
    gauge=[]
    for axis in range(3):
        U=background_matrix(axis,theta)
        transformed=background_matrix(axis,theta+B.T@chi)
        error=float(np.linalg.norm(transformed-G@U@G.conj().T,ord=2))
        assert error<1e-13 and np.linalg.norm(U.conj().T@U-np.eye(32))<1e-13
        gauge.append(error)
    # An omitted intracell current breaks Gauss on actual native port endpoints.
    wrong=route(0,0)[2].copy();e=next(e for e in np.flatnonzero(wrong) if field.data()[0].seams[e][0]//12==field.data()[0].seams[e][1]//12)
    wrong[e]=0
    assert np.any(B@wrong-B@route(0,0)[2])
    experiments=[]
    initial_spin=np.array([1,1j,1,-1j],complex)/2
    for alpha,beta,charge in ((0.,0.,1),(.17,.41,1),(.17,.41,-1)):
        state=initial(spinor=initial_spin,charge=charge)
        states=[]
        for n in range(3):
            states.append(records(state,charge))
            if n<2:state=step(state,.1,.4,alpha,beta,charge)
        experiments.append(dict(alpha=alpha,beta_face0=beta,other_face_couplings=0.,charge=charge,records=states))
    assert abs(experiments[0]['records'][-1]['electric_energy']-experiments[1]['records'][-1]['electric_energy'])>1e-6
    # Pauli identity including the exact oscillator-cutoff defect.
    N=8;qB=.2;a=np.diag(np.sqrt(np.arange(1,N)),1)
    px=np.sqrt(qB/2)*(a+a.T);py=-1j*np.sqrt(qB/2)*(a-a.T)
    comm=px@py-py@px;boundary=np.eye(N);boundary[-1,-1]-=N
    assert np.max(abs(comm-1j*qB*boundary))<1e-14
    kinetic=np.kron(px,ALPHA[0])+np.kron(py,ALPHA[1])+np.kron(np.eye(N),m*BETA)
    rhs=np.kron(px@px+py@py+m*m*np.eye(N),eye)-qB*np.kron(boundary,SIGMA[2])
    pauli_error=float(np.linalg.norm(kinetic@kinetic-rhs));assert pauli_error<1e-13
    # A spectator spin representation supplies a 2pi phase but no spin-B term.
    spectator=np.kron(np.diag([.1,.3]),np.eye(2));spin_up=np.diag(spectator)[::2];spin_down=np.diag(spectator)[1::2]
    assert np.array_equal(spin_up,spin_down)
    source={'model_id':MODEL_ID,'replaces':'W=S R or its scalar mass coin; no unitary equivalence claimed',
            'matter':'four co-located root-port modes per coarse cell; explicit Clifford matrices and optional explicit CAR lift',
            'macro_order':['mass exp(-i h m beta)','axis x','axis y','axis z','declared field substeps','common clock tick'],
            'clock':'t=n tau0 h; c_D=1/tau0; distinct from old 1/(sqrt48 tau0)',
            'central_field_schedule':{'field_substeps':8,'alpha':'1/16','beta':'3/8','alpha_beta':'3/128',
                                      'uniform_curvature_norm_bound':158,'stability_product':'237/64 < 4',
                                      'field_macro_speed':'8 sqrt(2 alpha beta/3) = 1',
                                      'scope':'harmonic full homogeneous face response; finite one-face experiments below use separate declared diagnostic couplings'},
            'architecture_inputs':['Clifford representation','root-port code and actual routes','m','charge sign','field couplings and substeps'],
            'spin_diagonalizers':native,'mass_word_example':mass_word(.1,.4),'routes':routes,
            'edge_increment_generator':'manybody/mode_field.py:field_block(16); paired opposite rails carry nu_a-nu_b',
            'rail_transport_factorization':'onsite plus/minus swap followed by disjoint nearest-cell minus/plus swap; internal paths require declared empty transit modes and current on every traversed seam',
            'full_native_status':'spin words executed; route-current identities and full sparse finite logical dynamics executed; full routed matter-field native pulse tape, transit reset and aggregate error budget are not yet executed',
            'rotation_scope':'continuum generator covariance and native spin-lift readouts; a complete finite carrier rotation apparatus must also move the spatial chart and Gauss dressing. Fixed xyz Trotter order is not exactly rotation covariant.'}
    result={'schema':'oph.local12.dirac_spin_proposal.v1','model_id':MODEL_ID,
            'checks':{'explicit_Clifford_architecture':True,'native_spin_words_all_columns':True,'native_two_pi_four_pi_phase_records':True,
                      'continuum_physical_rotation_covariance':True,'controlled_free_Dirac_limit':True,'full_Brillouin_doublers_retained':True,
                      'actual_carrier_routes_and_Gauss_current':True,'finite_character_gauge_covariance':True,'finite_joint_field_intervention':True,
                      'missing_intracell_current_rejected':True,'Pauli_identity_with_cutoff_defect':True,'spectator_g0_control':True},
            'Clifford_residual':clifford,'rotation_covariance_residuals':covariance,'gauge_residuals':gauge,
            'free_continuum_cases':bounds,'two_pi_native_records':two,'four_pi_native_records':four,
            'doublers':doublers,'finite_joint_experiments':experiments,'Pauli_identity_error':pauli_error,
            'conditional_magnetic_moment':'If the same routed gauge field converges to a canonical covariant derivative, H_D squared = Pi squared +m squared -q Sigma dot B, hence g=2 at leading positive-energy weak-field order. That interacting limit is not established by this finite Z16 experiment.',
            'spectator_control':'An independently supplied SU2 spin lift with spin-independent scalar dynamics has no intrinsic spin-field splitting; its g is zero.',
            'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(spin.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'DIRAC_SOURCE_PROPOSAL.json').write_text(json.dumps(source,indent=2)+'\n')
    (HERE/'DIRAC_SPIN.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'free_bounds':bounds,'native_spin_errors':[r['operator_error'] for r in native],
                      'field_energy_effect':experiments[1]['records'][-1]['electric_energy']-experiments[0]['records'][-1]['electric_energy']}))


if __name__=='__main__':run()
