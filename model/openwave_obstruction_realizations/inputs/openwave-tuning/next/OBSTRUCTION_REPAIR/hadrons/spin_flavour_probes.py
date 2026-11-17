"""Actual CAR state reconstruction and complete spectra for six triquark probes.

Exact limb arithmetic replays compressed source matrices. Probabilities remain
numerical diagnostics; all six highest-weight labels and failures are retained.
"""
from pathlib import Path
from collections import deque
from fractions import Fraction as F
import json
import numpy as np
import scipy.linalg as la
from spin_flavour_source import key,operators,ModRank,primitive,projection,exact_action
from staggered_vacuum import annihilation
HERE=Path(__file__).resolve().parent;CACHE=HERE/'spin_flavour_cache'

from exact_integer_product import exact_product


def state_frame(row,all_rows):
    lam=row['U4_highest_weight'];name=key(lam);saved=CACHE/f'CAR_frame_{name}.npz';block=json.loads((CACHE/f'block_{name}.json').read_text())
    if saved.exists():z=np.load(saved);return z['words'].tolist(),z['R'],block
    words,ops=operators(lam,row['joint_colour_and_flavour_weight_dimension']);n=len(words);rank=ModRank();basis=[];q=deque();rng=np.random.default_rng(314159+sum((i+1)*a for i,a in enumerate(lam)));expected=row['multiplicity']
    def add(v):
        v=primitive(v)
        if rank.add(v):basis.append(v);q.append(v);return True
        return False
    for seed in range(expected+1):
        if len(basis)==expected:break
        if seed==0 and lam==[3,3,3,3]:
            v=np.zeros(n,dtype=np.int64);v[words.index(((1<<12)-1)<<12)]=1
        else:v=projection(rng.integers(-1,2,n,dtype=np.int64),ops,lam,[r['U4_highest_weight'] for r in all_rows])
        add(v)
        while q and len(basis)<expected:
            v=q.popleft()
            for tag in ['M','E','K']:
                add(exact_action(ops[tag],v))
                if len(basis)==expected:break
    assert len(basis)==expected
    R=np.array(basis,dtype=np.int64).T
    assert exact_product(R,R).tolist()==block['G']
    for tag in ['M','K','E']:
        W=np.array([exact_action(ops[tag],v) for v in basis],dtype=np.int64).T
        assert exact_product(R,W).tolist()==block['components'][tag]
    np.savez_compressed(saved,words=np.array(words,dtype=np.int64),R=R)
    print('EXACT FRAME REPLAY',lam,R.shape,flush=True)
    return words,R,block

def eigenframe(R,block,kappa,m):
    G=block['G'];d=len(G);M,K,E=[block['components'][x] for x in ['M','K','E']]
    A=[[m*M[i][j]+F(1,8)*K[i][j]+kappa*F(E[i][j],6) for j in range(d)] for i in range(d)]
    gf=np.array(G,float);sc=1/np.sqrt(np.diag(gf));cond=np.linalg.cond(sc[:,None]*gf*sc[None,:])
    if cond<1e8:
        vals,vec=la.eigh(sc[:,None]*np.array(A,float)*sc[None,:],sc[:,None]*gf*sc[None,:]);V=R@(sc[:,None]*vec)
    else:
        import mpmath as mp
        with mp.workdps(90):
            L=mp.cholesky(mp.matrix(G));Li=L**-1;aa=mp.matrix([[mp.mpf(v.numerator)/v.denominator for v in row] for row in A]);D,U=mp.eigsy(Li*aa*Li.T);vec=Li.T*U
            V=R@np.array(vec.tolist(),float);vals=np.array(list(D),float)
    residual=float(np.linalg.norm(V.T@V-np.eye(d)));assert residual<1e-8
    return vals,V,residual

def main():
    census=json.loads((HERE/'SPIN_FLAVOUR_CENSUS.json').read_text());sectors={r['N']:r for r in census['sectors']};states={}
    for lam in [[3,3,3,3],[3,3,3,0],[3,3,2,1],[3,0,0,0],[2,1,0,0]]:
        sec=sectors[sum(lam)];row=next(x for x in sec['SU4_irreps'] if x['U4_highest_weight']==lam);states[tuple(lam)]=state_frame(row,sec['SU4_irreps'])
    inp=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text());m=F(inp['432']['declared_light_midpoint']);results=[]
    for kappa in [F(1,2),F(1),F(2)]:
        vw,VR,vb=states[3,3,3,3];_,vac,vr=eigenframe(VR,vb,kappa,m);psi=vac[:,0];records=[]
        for lam in [(3,0,0,0),(2,1,0,0)]:
            pw,PR,pb=states[lam];target=tuple(3-x for x in reversed(lam));tw,TR,tb=states[target];ti={w:i for i,w in enumerate(tw)};energies,EV,orth=eigenframe(TR,tb,kappa,m)
            for label in range(PR.shape[1]):
                norm=sum(abs(int(x)) for x in PR[:,label]);emitted=np.zeros(len(tw))
                for op,coef in zip(pw,PR[:,label]):
                    if not coef:continue
                    modes=[j for j in range(24) if op>>j&1]
                    # Weyl reversal of internal labels, in the original CAR
                    # monomial order; anticommutation signs are then literal.
                    modes=[12*(j//12)+3*(3-((j%12)//3))+j%3 for j in modes]
                    for initial,amp in zip(vw,psi):
                        word=initial;sg=1
                        for mode in modes:
                            word,z=annihilation(mode,word)
                            if word is None:break
                            sg*=z
                        if word is not None:emitted[ti[word]]+=float(coef)*sg*amp/norm
                weights=abs(EV.T@emitted)**2;success=float(weights.sum());actual=float(emitted@emitted)
                assert abs(success-actual)<1e-8 and -1e-10<=success<=1+1e-10
                records.append(dict(source_operator_U4_highest_weight=list(lam),polynomial_label=label,
                    target_U4_highest_weight=list(target),exact_l1_normalizer=norm,
                    success_probability_diagnostic=success,failure_probability_diagnostic=1-success,
                    all_target_spectral_probabilities_diagnostic=weights.tolist(),
                    all_target_energies_diagnostic=energies.tolist(),orthogonality_residual=orth,
                    baryon_change=-1,reservoir_baryon_change_on_success=1))
        results.append(dict(kappa=str(kappa),vacuum_U4_highest_weight=[3,3,3,3],vacuum_frame_orthogonality_residual=vr,all_six_highest_weight_probe_records=records))
    out=dict(schema='oph.spin_flavour.probes.v1',actual_CAR_state_maps_replayed_exactly=True,
        all_five_frames_reproduce_every_compressed_source_component_exactly=True,
        all_six_highest_weight_triquark_polynomial_labels_retained=True,
        all_target_spectral_lines_and_failed_probe_probabilities_retained=True,
        full_SU4_multiplet_apparatus_executed=False,native_program_executed=False,
        probabilities_are_diagnostics_not_certified_enclosures=True,results=results)
    (HERE/'SPIN_FLAVOUR_PROBES.json').write_text(json.dumps(out,indent=2)+'\n');print('ALL PROBES COMPLETE',flush=True)
if __name__=='__main__':main()
