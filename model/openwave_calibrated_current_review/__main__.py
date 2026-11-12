from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product, permutations
import argparse, hashlib, json
import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm
from model.common_source import calibrated_matter as native
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.openwave_calibrated_hadrons.source import sector, parameters, INPUTS
from model.openwave_baryon_weak_formfactors.native import currents, quartet
from model.openwave_baryon_weak_formfactors.certificate import exact_local
from model.openwave_baryon_beta.decay import construction, eigensystem, evolve
from model.openwave_nuclear_refinement.levels import local_census
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TARGETS=('openwave_baryon_weak_formfactors','openwave_baryon_beta')

def pins():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    for name in TARGETS:
        path=ROOT/'model'/name/'results.json'
        receipt=json.loads(path.read_text())
        for p,h in receipt['source_sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h, p
        paths.update(receipt['source_sha256']);paths.add(str(path.relative_to(ROOT)))
    for name in ('openwave_calibrated_hadrons','openwave_nuclear_refinement'):
        path=ROOT/'model'/name/'results.json';receipt=json.loads(path.read_text())
        for p,h in receipt['source_sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
        paths.update(receipt['source_sha256']);paths.add(str(path.relative_to(ROOT)))
    assert not any(p.startswith('docs/') for p in paths)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}

def one_body(A,basis):
    """Direct bit-level CAR construction, independent of target operator helper."""
    where={int(mask):i for i,mask in enumerate(basis)};rows=[];cols=[];values=[]
    for a,b in zip(*np.nonzero(A)):
        a=int(a);b=int(b)
        for col,mask in enumerate(basis):
            mask=int(mask)
            if not (mask>>b)&1:continue
            v=mask^(1<<b)
            if (v>>a)&1:continue
            sign=(-1)**((mask&((1<<b)-1)).bit_count()+(v&((1<<a)-1)).bit_count())
            rows.append(where[v|(1<<a)]);cols.append(col);values.append(sign*A[a,b])
    return sp.coo_matrix((values,(rows,cols)),shape=(len(basis),len(basis))).tocsr()

def maximum(matrix):
    return float(np.max(abs(matrix.data),initial=0)) if sp.issparse(matrix) else float(np.max(abs(np.asarray(matrix)),initial=0))

def independent_currents():
    d=sector('baryon');modes=d['modes'];idx={tuple(mode):i for i,mode in enumerate(modes)}
    subsets,*_=representation();raising=np.asarray(native.weak_raising_current(families=False),dtype=complex)
    third=(raising@raising.conj().T-raising.conj().T@raising)/2
    qi=[[subsets.index((c,w)) for c in range(3)] for w in (3,4)]
    sigmas=(np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
    full=[];reduced=[];errors=[];target=currents()
    for x in range(2):
        fr=[];rr=[]
        for mu,sigma in enumerate(sigmas):
            fs=[];rs=[]
            for k,weak in enumerate((raising,third)):
                A=np.zeros((24,24),complex)
                for s,z,w,v,c in product(range(2),range(2),range(2),range(2),range(3)):
                    A[idx[(0,x,s,w,c)],idx[(0,x,z,v,c)]]=sigma[s,z]*weak[qi[w][c],qi[v][c]]
                J=one_body(A,d['basis']);j=(d['B'].conj().T@J@d['B']).toarray()
                errors.extend((maximum(J-target['full'][x][mu][k]),maximum(J@d['B']-d['B']@j)))
                fs.append(J);rs.append(j)
            fr.append(fs);rr.append(rs)
        full.append(fr);reduced.append(rr)
    assert max(errors)<4e-13
    return d,full,reduced,sigmas,raising,third,max(errors)

def calculate():
    d,full,reduced,sigmas,raising,third,current_error=independent_currents()
    p=parameters();local=exact_local();v=list(map(int,local['integer_neutron_vector']))
    norm=sum(x*x for x in v);assert norm==int(local['integer_vector_norm_squared'])
    jp=sum(full[x][0][0] for x in range(2));jz=sum(full[x][3][0] for x in range(2))
    a=jp@np.array(v);b=jz@np.array(v)
    assert maximum(a.imag)==0 and maximum(b.imag)==0
    a=[int(round(x.real)) for x in a];b=[int(round(x.real)) for x in b]
    assert F(sum(x*y for x,y in zip(a,b)),norm)==F(5,3)
    assert F(sum((3*y-5*x)**2 for x,y in zip(a,b)),9*norm)==F(8,9)
    assert local_census(3)['exact_native_local_nucleon_axial_ratio']=='5/3'
    gamma=F(8,3)*p['U']-F(4,9)*p['A']-F(3,4)*(p['Js']+p['Jw'])-3*p['t']
    axial_error=36*p['t']/gamma;assert F(5,3)-axial_error>F(12753,10000)
    beta=construction();masks=[sum(1<<i for i in row) for row in combinations(range(8),4)]
    def index(w,x,s):return 4*w+2*x+s
    ne=np.array([(mask>>4).bit_count() for mask in masks])
    lz=np.array([sum(1-2*s for w,x,s in product(range(2),repeat=3) if mask>>index(w,x,s)&1) for mask in masks])
    qh=np.rint(np.diag(d['q']).real).astype(int);zh=np.rint(2*np.diag(d['z']).real).astype(int)
    charge=(qh[:,None]-ne).ravel();spin=(zh[:,None]+lz).ravel()
    keep=np.flatnonzero((charge==0)&(spin==1));outside=np.flatnonzero((charge!=0)|(spin!=1))
    assert np.array_equal(keep,beta['keep']) and len(keep)==398 and np.sum(charge==0)==1400
    counts=[int(np.sum((charge==0)&(np.tile(ne,120)==k))) for k in range(3)]
    assert counts==[40,640,720]
    lp=[];V=sp.csr_matrix((8400,8400),dtype=complex)
    for x in range(2):
        row=[]
        for mu,sigma in enumerate(sigmas):
            ops=[]
            for weak in (raising[np.ix_([10,11],[10,11])],third[np.ix_([10,11],[10,11])]):
                A=np.zeros((8,8),complex)
                for w,z,s,r in product(range(2),repeat=4):A[index(w,x,s),index(z,x,r)]=sigma[s,r]*weak[w,z]
                ops.append(one_body(A,masks))
            lplus,l3=ops;qplus,q3=map(sp.csr_matrix,reduced[x][mu]);eta=1 if mu==0 else -1
            V+=eta*(sp.kron(qplus,lplus.conj().T)+sp.kron(qplus.conj().T,lplus)+2*sp.kron(q3,l3))
            row.append(ops)
        lp.append(row)
    V=V.tocsr();outside_error=maximum(V[keep][:,outside]);source_error=maximum(V[keep][:,keep]-beta['V'])
    assert outside_error<1e-12 and source_error<2e-12
    # Exterior fourth power of an integer Hadamard matrix: every minor/4 is exact.
    U=np.zeros((8,8),dtype=np.int64)
    for w,s,k,x in product(range(2),repeat=4):U[index(w,k,s),index(w,x,s)]=(-1)**(k*x)
    occ=[[i for i in range(8) if mask>>i&1] for mask in masks]
    perms=[(row,(-1)**sum(row[i]>row[j] for i in range(4) for j in range(i+1,4))) for row in permutations(range(4))]
    def minor(a,b):
        return sum(sign*int(U[a[0],b[row[0]]])*int(U[a[1],b[row[1]]])*int(U[a[2],b[row[2]]])*int(U[a[3],b[row[3]]]) for row,sign in perms)
    Fnum=np.array([[minor(a,b) for b in occ] for a in occ],dtype=np.int64)
    assert np.array_equal(Fnum@Fnum.T,16*np.eye(70,dtype=np.int64))
    Fourier=Fnum/4;T=sp.kron(sp.eye(120),sp.csr_matrix(Fourier)).tocsr()[keep][:,keep]
    fourier_error=maximum(T-beta['T']);assert fourier_error<2e-14
    # Full occupation diagonal after readout, including colour decoherence.
    q=float((p['delta']-1)/2);A=np.zeros((8,8))
    for w,s in product(range(2),repeat=2):
        ids=[index(w,x,s) for x in range(2)];K=q/2*np.array([[1,-1],[-1,1]])
        A[np.ix_(ids,ids)]=np.eye(2)+K if w else -K
    Hl=one_body(A,masks)+2*q*sp.eye(70)
    rawpos=d['H'].diagonal().real[:,None]+Hl.diagonal().real[None,:]
    rawmom=d['H'].diagonal().real[:,None]+np.diag(Fourier@Hl@Fourier.T).real[None,:]
    for x in range(2):
        for mu in range(4):
            eta=1 if mu==0 else -1
            for k,coef in ((0,1),(1,2)):
                qop=full[x][mu][k];lop=lp[x][mu][k]
                if k==0:
                    a=qop.diagonal()[:,None]*lop.conj().T.diagonal()[None,:]
                    b=qop.diagonal()[:,None]*np.diag(Fourier@lop.conj().T@Fourier.T)[None,:]
                    rawpos+=beta['coupling']*eta*2*a.real;rawmom+=beta['coupling']*eta*2*b.real
                else:
                    rawpos+=beta['coupling']*eta*coef*(qop.diagonal()[:,None]*lop.diagonal()[None,:]).real
                    rawmom+=beta['coupling']*eta*coef*(qop.diagonal()[:,None]*np.diag(Fourier@lop@Fourier.T)[None,:]).real
    weight=abs(d['B']).power(2).T
    readout_error=max(maximum(np.asarray(weight@rawpos).ravel()[keep]-beta['measured_energy']),maximum(np.asarray(weight@rawmom).ravel()[keep]-beta['measured_momentum_energy']))
    assert readout_error<2e-10
    W,_=quartet(d['h']);reference=np.eye(70)[:,masks.index(15)]
    means=[]
    for hstate in (W[:,0],W[:,2]):
        state=np.kron(hstate,reference);means.append(float(np.vdot(state,V@state).real))
    assert abs(means[0]+1)<2e-12 and abs(means[1]-1)<2e-12
    eigen=eigensystem(beta);psi=evolve(beta,512,eigen)
    pair=float(np.sum(abs(psi[beta['ne']==1])**2));two=float(np.sum(abs(psi[beta['ne']==2])**2))
    assert pair>.96 and two>1e-4
    ref=expm(-16j*(beta['H'].toarray()-eigen[2]*np.eye(398)))@beta['initial']
    prop_error=float(np.linalg.norm(ref-evolve(beta,16,eigen)));assert prop_error<2e-8
    failure=evolve(beta,16,eigen,beta['failure']);assert abs(np.vdot(failure,failure).real+beta['success']-1)<2e-11
    return {'schema':'oph.calibrated_native_current_independent_review.v1','verdict':'PASS_NATIVE_CURRENTS_COMPLETE_BETA_AND_EXPLICIT_PHYSICAL_FAILURES',
        'source_sha256':pins(),'direct_bit_CAR_current_and_intertwiner_error':current_error,
        'independent_exact_local_projector_gA':'5/3','exact_axial_leakage_squared':'8/9',
        'certified_axial_interval_rational':[str(F(5,3)-axial_error),str(F(5,3)+axial_error)],
        'scalar_local_repair_obstruction':'Multiplicity-one native (S,I)=(1/2,1/2); an invariant scalar cannot change the fixed native current matrix element.',
        'charge_zero_electron_count_dimensions':counts,'complete_charge_dimension':1400,'spin_projection_dimension':398,
        'independent_full_contact_restriction_error':source_error,'discarded_spin_charge_block_error':outside_error,
        'exact_exterior_Fourier_denominator':4,'Fourier_comparison_error':fourier_error,
        'full_colour_occupation_readout_work_diagonal_error':readout_error,
        'filled_neutral_reference_contact_means':means,
        'contact_expectation_shift_magnitude_GeV_rational':str(F(1,64)*F(INPUTS['electron_unit_GeV'])),
        'conditional_pair_probability_at_512':pair,'conditional_two_pair_probability_at_512':two,
        'preparation_failure_probability':1-beta['success'],'independent_Pade_vs_spectral_error_at_16':prop_error,
        'scope':'Strong-source p/n/pion eigenmasses are explicitly calibrated. The beta contact shifts filled-reference energy expectations by opposite G*m_e; the strong/free endpoint is not an interacting pole endpoint. Finite coherent conversion retains virtual two-pair states and all failed preparations. No parity-violating angular law, physical lifetime or successful axial-ratio match is established.'}

def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=parser.parse_args()
    before=pins();out=calculate();assert before==pins();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),out)
    print(out['verdict'])

if __name__=='__main__':main()
