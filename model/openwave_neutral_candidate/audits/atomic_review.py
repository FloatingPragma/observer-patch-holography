"""Independent dense Green, parity and min--max checks for finite atoms."""
from itertools import product
from pathlib import Path
import hashlib
import json
import numpy as np
import scipy.linalg as la

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_atomic_refinement'


def dense_atom(n,length):
    a=length/(n+1);points=list(product(range(n),repeat=3));index={p:i for i,p in enumerate(points)}
    lap=6*np.eye(n**3)
    for p,i in index.items():
        for axis in range(3):
            for sign in (-1,1):
                q=list(p);q[axis]+=sign;q=tuple(q)
                if q in index:lap[i,index[q]]=-1
    rho=np.zeros(n**3);rho[index[(n//2,)*3]]=1
    green=la.solve(lap,rho,assume_a='pos')
    potential=-4*np.pi*green/a
    return lap/(2*a*a)+np.diag(potential),potential,lap,points,index


def independent_frame(n,parity):
    points=list(product(range(n),repeat=3));index={p:i for i,p in enumerate(points)}
    center=n//2;columns=[]
    for representative in product(*(range(center+int(odd),n) for odd in parity)):
        vector=np.zeros(n**3)
        orbit={}
        for flip in product((0,1),repeat=3):
            p=tuple(n-1-x if f else x for x,f in zip(representative,flip))
            sign=(-1)**sum(f*odd for f,odd in zip(flip,parity))
            orbit[p]=sign
        for p,sign in orbit.items():vector[index[p]]=sign/np.sqrt(len(orbit))
        columns.append(vector)
    return np.column_stack(columns)


def report():
    from model.openwave_field_observables.orbitals import reduced_hamiltonian
    from model.openwave_atomic_refinement.refinement import spectrum
    paths=[TARGET/p for p in ('refinement.py','__main__.py','PROOF.md','results.json')]
    paths += [ROOT/'model/openwave_field_observables/orbitals.py',Path(__file__)]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    saved=json.loads((TARGET/'results.json').read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in saved['source_sha256'].items())
    h,pot,lap,points,index=dense_atom(5,6.)
    allframes=[];allvalues=[];errors=[];dims=[]
    for parity in product((False,True),repeat=3):
        frame=independent_frame(5,parity)
        assert np.max(abs(frame.T@frame-np.eye(frame.shape[1])))<1e-14
        block=frame.T@h@frame;produced,a=reduced_hamiltonian(5,6.,parity)
        errors.append(float(np.max(abs(block-produced.toarray()))));dims.append(frame.shape[1])
        assert errors[-1]<2e-14
        allframes.append(frame);allvalues.extend(la.eigvalsh(block))
    frame=np.column_stack(allframes)
    assert frame.shape==(125,125) and np.max(abs(frame@frame.T-np.eye(125)))<1e-14
    union_error=float(np.max(abs(np.sort(allvalues)-la.eigvalsh(h))))
    assert union_error<3e-14
    # Fixed-plane mutation loses orthonormality and changes the kinetic form.
    bad=allframes[0].copy();center=index[(2,2,2)];bad[center,0]*=np.sqrt(2)
    assert la.norm(bad.T@bad-np.eye(bad.shape[1]))>0.9
    comparisons=[]
    for ns,nl in ((3,5),(5,7)):
        hs,vs,_,ps,ids=dense_atom(ns,ns+1.)
        hl,vl,_,pl,idl=dense_atom(nl,nl+1.)
        offset=(nl-ns)//2;embedding=np.zeros((nl**3,ns**3))
        for p,i in ids.items():embedding[idl[tuple(x+offset for x in p)],i]=1
        difference=embedding.T@hl@embedding-hs
        assert np.max(abs(difference-np.diag(np.diag(difference))))<1e-14
        assert np.max(np.diag(difference))<0
        es=la.eigvalsh(hs);el=la.eigvalsh(hl)
        assert np.max(el[:len(es)]-es)<0
        comparisons.append({'small':ns,'large':nl,'spacing':1.,
            'largest_restricted_quadratic_form_change':float(np.max(np.diag(difference))),
            'largest_ordered_eigenvalue_change':float(np.max(el[:len(es)]-es))})
    eig,vec=la.eigh(h)
    positive=abs(vec[:,0])+1e-12;ratios=(h@positive)/positive
    assert min(ratios)-1e-13<=eig[0]<=max(ratios)+1e-13
    wrong_residual=float(la.norm(h@vec[:,-1]-eig[-1]*vec[:,-1]))
    assert wrong_residual<1e-13 and eig[-1]-eig[0]>1
    small=spectrum(5,6.)
    assert abs(small['sectors'][0]['energies'][0]-eig[0])<1e-11
    assert all(row['spacing']==1/3 for row in saved['volume_rows'])
    assert all(row['length']==24 for row in saved['spacing_rows'])
    assert [row['grid'] for row in saved['volume_rows']]==[71,95,119]
    assert [row['grid'] for row in saved['spacing_rows']]==[71,95,127,159]
    assert saved['volume_rows'][0]==saved['spacing_rows'][0]
    assert not saved['continuum_spacing_error_certified']
    result={'schema':'oph.independent_atomic_refinement_review.v1','input_sha256':pins,
        'dense_source_dimension':125,'complete_parity_dimensions':dims,
        'largest_compressed_operator_error':max(errors),'complete_spectral_union_error':union_error,
        'fixed_spacing_nested_comparisons':comparisons,
        'incorrect_reflection_plane_normalization_rejected':True,
        'wrong_index_small_residual_control':{'residual':wrong_residual,'index_is_ground':False},
        'independent_positivity_ground_interval':[float(min(ratios)),float(max(ratios))],
        'volume_and_spacing_grids_checked':True,
        'scope':'Dense floating calculations independently verify the finite construction and theorem premises. Neither residuals nor successive-grid differences certify continuum error or an observed hydrogen spectrum.',
        'verdict':'PASS_INDEPENDENT_ATOMIC_REFINEMENT_REVIEW'}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==d for p,d in pins.items())
    return result


if __name__=='__main__':
    result=report();(HERE.parent/'atomic_review.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(result['verdict'])
