"""Independent full-Fock controls and kinetic product-formula review."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[3]


def lower(n,mode):
    out=np.zeros((2**n,2**n),complex)
    for mask in range(2**n):
        if mask&(1<<mode):out[mask^(1<<mode),mask]=(-1)**sum(bool(mask&(1<<j)) for j in range(mode))
    return out


def run():
    n=8;c=[lower(n,j) for j in range(n)];dim=2**n
    u=np.array([np.sqrt(7/29),np.exp(.37j)*np.sqrt(22/29)]);W=np.outer(u,u.conj())
    V=np.zeros((dim,dim),complex)
    for start in (0,4):
        for i in range(2):
            for j in range(2):V+=W[i,j]*c[start+i].conj().T@c[start+2+j].conj().T
    V+=V.conj().T
    K=np.diag([.8*(sum((mask>>j)&1 for j in range(4))-sum((mask>>j)&1 for j in range(4,8))) for mask in range(dim)])
    one=np.zeros((n,n));perm=(4,5,6,7,0,1,2,3)
    for j,t in enumerate(perm):one[t,j]=1
    # Exterior signs computed independently as Slater determinants.
    exterior=np.zeros((dim,dim),complex)
    for mask in range(dim):
        occupied=[j for j in range(n) if mask>>j&1];dest=sorted(perm[j] for j in occupied)
        value=round(np.linalg.det(one[np.ix_(dest,occupied)])) if occupied else 1
        exterior[sum(1<<j for j in dest),mask]=value
    assert all(np.linalg.norm(exterior@c[j]@exterior.conj().T-c[perm[j]])<1e-12 for j in range(n))
    leg=np.diag([(-1)**sum((mask>>j)&1 for j in (0,1,4,5)) for mask in range(dim)])
    inverse=exterior@leg;errors=[]
    for index in range(4):
        R=np.diag([1j**(index*sum((mask>>j)&1 for j in (1,3,5,7))) for mask in range(dim)])
        H=K+R.conj().T@V@R
        errors.append(float(np.linalg.norm(inverse@H@inverse.conj().T+H)))
        assert errors[-1]<1e-12
    wrong=float(np.linalg.norm(leg@(K+V)@leg+(K+V)));assert wrong>1
    unsigned=abs(exterior)
    assert np.linalg.norm(unsigned@c[0]@unsigned.T-c[4])>1
    # Independent4-mode physical embedding with complex W and noncommutingK.
    c4=[lower(4,j) for j in range(4)];v=np.zeros((16,16),complex)
    for i in range(2):
        for j in range(2):v+=W[i,j]*c4[i].conj().T@c4[2+j].conj().T
    v+=v.conj().T;k=np.diag([.8*mask.bit_count() for mask in range(16)])
    rotations=[np.diag([1j**(index*(((mask>>1)&1)+((mask>>3)&1))) for mask in range(16)]) for index in range(4)]
    H=[k+R.conj().T@v@R for R in rotations];average=sum(H)/4
    J=np.eye(16)[:,[0,1,4,5]];target=J.conj().T@(k+v)@J
    assert np.linalg.norm(average@J-J@target)<1e-13
    assert np.linalg.norm(k@v-v@k)>1
    p=1/(4-4**(1/3));time=.7;reference=J@expm(-1j*time*target)
    runs=[]
    for rounds in (16,32,64):
        def second(h):
            result=np.eye(16,dtype=complex)
            for i,w in ((0,1/8),(1,1/8),(2,1/8),(3,1/4),(2,1/8),(1,1/8),(0,1/8)):
                result=expm(-1j*h*w*H[i])@result
            return result
        A=second(p*time/rounds);B=second((1-4*p)*time/rounds)
        actual=np.linalg.matrix_power(A@A@B@A@A,rounds)@J
        error=float(np.linalg.norm(actual-reference,2))
        av=actual.ravel()/np.linalg.norm(actual);bv=reference.ravel()/np.linalg.norm(reference)
        overlap=np.vdot(bv,av)
        # The perpendicular residual avoids cancellation in1-|overlap|^2.
        choi=2*np.linalg.norm(av-bv*overlap)
        radius=(7/3)*(4*.8+1);bound=radius**5*time**5/(10*rounds**4)
        assert 2*error<bound
        probabilities=[float(np.linalg.norm(actual[[word|mask for word in (0,1,4,5)],:])**2/4) for mask in (0,2,8,10)]
        assert abs(sum(probabilities)-1)<1e-12
        runs.append({'rounds':rounds,'physical_operator_error':error,'complete_reference_Choi_trace_distance':float(choi),
            'floating_isometry_Choi_norm_defect':float(abs(np.linalg.norm(actual)/2-1)),
            'analytic_complete_trace_bound':bound,'all_ancillary_probabilities':probabilities})
    assert runs[0]['physical_operator_error']/runs[1]['physical_operator_error']>15
    receipt_path=ROOT/'model/openwave_kinetic_phase_cycle/results.json';receipt=json.loads(receipt_path.read_text())
    for p,h in receipt['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    paths=[*receipt['source_sha256'],'model/openwave_kinetic_phase_cycle/results.json',str(Path(__file__).relative_to(ROOT))]
    return {'schema':'oph.independent_kinetic_phase_review.v1','verdict':'PASS',
        'full_Fock_control_dimension':256,'complex_rank_one_W_effect':str(F(7,29)),
        'complete_signed_inverse_residuals':errors,'leg_parity_only_inverse_mutation':wrong,
        'independent_CAR_Slater_permutation':True,'unsigned_permutation_mutation_rejected':True,
        'independent_complete_evolution_runs':runs,
        'mathematical_review':'Spatial inversion reverses the homogeneous Weyl term and leaves onsite pairing; complete-leg parity performs the other sign reversal. Their product supplies actual negative full-H segments. The 24s active-mode restriction is an invariant tensor factorization with free spectators evolving for the same signed time, not a carrier deletion. The fourth-order word identity and full-Fock radius give the stated reference-uniform bound.',
        'scope_review':'The analyzer phase/background and global coherent controls are supplied. Negative-time compilation does not imply laboratory pulse time equals target evolution time, or derive interacting charged minima. Compiler speed alone does not establish a lepton mass spectrum.',
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(paths))}}


if __name__=='__main__':
    value=run();Path(__file__).with_suffix('.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');print(value['verdict'])
