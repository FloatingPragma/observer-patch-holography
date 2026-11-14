"""Independent Codex-worker review: preparation and continuing apparatus."""
from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.fiber_preparation import FiberPreparation

def require(ok,message):
    if not ok:raise ValueError(message)

def run():
    reports=[]
    for modulus in (4,6):
        prep=FiberPreparation(Construction(3,((0,1),(1,2),(2,0)),.23),modulus)
        full,v=prep.tiny_isometry();inc=np.zeros((6,6),int)
        for j,(a,b) in enumerate(prep.edges):inc[a,j]=1;inc[b,j]=-1
        allwords=np.array(list(itertools.product(range(modulus),repeat=6)),int)
        divergences=(allwords@inc.T)%modulus
        circuit=np.zeros_like(v)
        for position in range(6):
            demand=(np.eye(6,dtype=int)[0]-np.eye(6,dtype=int)[position])%modulus
            fiber=allwords[np.all(divergences==demand,axis=1)]
            require(len(fiber)==modulus,'full Gauss enumeration cardinality differs')
            declared={prep.record_word(position,(a,)) for a in range(modulus)}
            require(declared==set(map(tuple,fiber)),'preparation misses Gauss records')
            for word in fiber:circuit[full.index[(position,tuple(word))],position]=1/np.sqrt(modulus)
        # Full-register arithmetic is reversible even away from zero ancillas:
        # F'=F+solution(c,pos); c'=c-F'_chord, then invert in reverse order.
        rng=np.random.default_rng(modulus);inverse_error=0
        for _ in range(100):
            pos=int(rng.integers(6));a=int(rng.integers(modulus));f=rng.integers(modulus,size=6)
            shifted=(f+np.array(prep.record_word(pos,(a,))))%modulus
            work=(a-shifted[prep.chords[0]])%modulus
            restored_a=(work+shifted[prep.chords[0]])%modulus
            restored_f=(shifted-np.array(prep.record_word(pos,(int(restored_a),))))%modulus
            inverse_error=max(inverse_error,int(np.max(abs(restored_f-f))),int(abs(restored_a-a)))
        # Fourier of |0> has exactly uniform amplitude; zero record registers
        # make the work register uncompute to zero for every coordinate.
        for a in range(modulus):
            word=prep.record_word(2,(a,));require((a-word[prep.chords[0]])%modulus==0,'Fourier work ancilla not erased')
        psi=np.array([1,2j,0,1-1j,0,0],complex);psi/=np.linalg.norm(psi);rho=np.outer(psi,psi.conj())
        _,branches=prep.record_measurement(rho,0);lift=circuit@rho@circuit.conj().T;code=circuit@circuit.conj().T;outside=np.eye(full.n)-code
        branch_error=0.;continuation_error=0.;cross=0.;probability=0.;u=expm(-.12j*full.hamiltonian)
        for row in branches:
            mask=np.array([state[1][0]==row['residue'] for state in full.states]);expected=lift*mask[:,None]*mask[None,:]
            branch_error=max(branch_error,float(np.max(abs(expected-row['poststate_subnormalized']))));probability+=row['probability']
            cross=max(cross,float(np.linalg.norm(code@expected@outside)))
            if row['probability']>1e-14:
                normalized=expected/row['probability'];coherent=u@normalized@u.conj().T
                direct=.8*coherent+.2*sum((k@coherent@k.conj().T for ops in full.operators for _,k in ops))/len(full.operators)
                actual,_=full.step(normalized,.12,.2);continuation_error=max(continuation_error,float(np.max(abs(actual-direct))))
        cut_error=0.
        for bits in range(64):
            region=tuple(i for i in range(6) if bits>>i&1)
            for cut in prep.closed_cut_measurement(rho,region):
                mask=np.array([int(state[0] in region)==cut['regional_charge'] for state in full.states])
                measured=lift*mask[:,None]*mask[None,:]
                compressed=circuit@cut['poststate_subnormalized']@circuit.conj().T
                cut_error=max(cut_error,float(np.max(abs(measured-compressed))))
        require(cut_error<1e-12,'closed cut continuation not in uniform code')
        require(inverse_error==0 and branch_error<1e-12 and continuation_error<1e-12 and abs(probability-1)<1e-12 and cross>.01,'prep/apparatus audit failed')
        reports.append({'modulus':modulus,'all_Gauss_fibers_enumerated':True,'Fourier_tree_isometry_error':float(np.max(abs(circuit-v))),'arbitrary_register_inverse_error':inverse_error,'complete_record_branch_error':branch_error,'all64_closed_cut_branch_intertwiner_error':cut_error,'retained_code_complement_crossblock_norm':cross,'full_native_continuation_error':continuation_error})
    paths=[Path(__file__),*[ROOT/'candidate/oph_openwave'/name for name in ['fiber_preparation.py','construction.py','uniform_native.py','modular_native.py']]]
    return {'schema':'oph.fiber_preparation.independent_review.v1','reviewer':'Independent Codex worker; not an opposite-team Fable audit','outcome':'PASS_WITH_DECLARED_RESOURCES','controls':reports,'scope':'Ideal added modular Fourier and reversible arithmetic gates. Composite modulus works via unit-coefficient tree elimination. Full branch state retains code/complement crossblocks; weights are diagnostics, not another P/Q measurement. No general efficient gate-count or native-only preparation claim.','full_criterion_met':False,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
if __name__=='__main__':
    r=run();Path(__file__).with_name('FIBER_PREPARATION_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
