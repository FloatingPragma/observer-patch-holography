"""Independent graph-Green and full-parity audit of the finite atom."""
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, eigsh
from model.openwave_atomic_refinement.refinement import spectrum
from model.openwave_field_observables.orbitals import reduced_hamiltonian

ROOT=Path(__file__).resolve().parents[3]


def independent(n,length,sign=1):
    a=length/(n+1)
    rows=[];cols=[];data=[]
    for x in itertools.product(range(n),repeat=3):
        i=np.ravel_multi_index(x,(n,)*3)
        rows.append(i);cols.append(i);data.append(6.)
        for j in range(3):
            for step in (-1,1):
                y=list(x);y[j]+=step
                if 0<=y[j]<n:
                    rows.append(i);cols.append(np.ravel_multi_index(tuple(y),(n,)*3));data.append(-1.)
    lap=sparse.csr_matrix((data,(rows,cols)),shape=(n**3,n**3))
    rho=np.zeros(n**3);rho[np.ravel_multi_index((n//2,)*3,(n,)*3)]=1.
    green=spsolve(lap,rho)
    pot=-sign*4*np.pi*green/a
    return (lap/(2*a*a)+sparse.diags(pot)).tocsr(),green,pot


def frame(n,odd):
    sites=list(itertools.product(*(range(n//2+int(p),n) for p in odd)))
    result=np.zeros((n**3,len(sites)))
    for column,x in enumerate(sites):
        orbit={}
        for flips in itertools.product((0,1),repeat=3):
            y=tuple(n-1-x[j] if flips[j] else x[j] for j in range(3))
            value=(-1)**sum(flips[j]*odd[j] for j in range(3))
            if y in orbit:assert orbit[y]==value
            orbit[y]=value
        for y,value in orbit.items():result[np.ravel_multi_index(y,(n,)*3),column]=value/np.sqrt(len(orbit))
    assert np.max(abs(result.T@result-np.eye(len(sites))))<1e-13
    return result


def run():
    checks=[]
    for n,length in ((5,6.),(7,8.)):
        h,g,p=independent(n,length)
        full=np.linalg.eigvalsh(h.toarray());all_values=[];block_error=0.
        for odd in itertools.product((False,True),repeat=3):
            u=frame(n,odd);block=u.T@(h@u)
            observed,_=reduced_hamiltonian(n,length,odd)
            block_error=max(block_error,float(np.max(abs(block-observed.toarray()))))
            all_values.extend(np.linalg.eigvalsh(block))
        spectral_error=float(np.max(abs(np.sort(all_values)-full)))
        assert block_error<1e-12 and spectral_error<1e-11
        actual=spectrum(n,length)
        expected=[]
        for odd in ((False,False,False),(True,False,False)):
            u=frame(n,odd);expected.append(np.linalg.eigvalsh(u.T@(h@u)))
        for row,want in zip(actual['sectors'],expected):
            assert np.max(abs(np.array(row['energies'])-want[:len(row['energies'])]))<2e-11
        wrong,_,_=independent(n,length,-1)
        assert eigsh(wrong,k=1,which='SA',return_eigenvectors=False)[0]>0
        # Fixed-plane normalization mutation destroys the orthonormal reduction.
        u=frame(n,(False,False,False));bad=u.copy();bad[n**3//2]*=np.sqrt(2)
        mutation=float(np.max(abs(bad.T@bad-np.eye(bad.shape[1]))));assert mutation>.5
        checks.append({'grid':n,'all_parity_dimensions':len(all_values),'full_dimension':n**3,
            'maximum_block_error':block_error,'full_spectrum_error':spectral_error,
            'repulsive_potential_has_positive_lowest_level':True,
            'fixed_plane_normalization_mutation':mutation})
    hs,gs,ps=independent(5,6.);hl,gl,pl=independent(9,10.)
    change=pl.reshape((9,)*3)[2:7,2:7,2:7].ravel()-ps
    assert np.max(change)<0
    for odd in itertools.product((False,True),repeat=3):
        us=frame(5,odd);ul=frame(9,odd)
        es=np.linalg.eigvalsh(us.T@(hs@us))
        el=np.linalg.eigvalsh(ul.T@(hl@ul))
        assert np.all(el[:len(es)]<=es+1e-11)
    receipt=json.loads((ROOT/'model/openwave_atomic_refinement/results.json').read_text())
    for path,digest in receipt['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    paths=[*receipt['source_sha256'],'model/openwave_atomic_refinement/results.json',str(Path(__file__).relative_to(ROOT))]
    return {'schema':'oph.independent_atomic_refinement_review.v1','verdict':'PASS',
        'independent_method':'Direct sparse graph assembly and linear solve; full dense spectra versus eight independently normalized reflection-orbit frames.',
        'finite_grid_checks':checks,'fixed_spacing_all_parity_domain_comparison':True,
        'maximum_nested_potential_change':float(np.max(change)),
        'mutation_controls':['repulsive charge sign','wrong central-plane orbit normalization'],
        'review_conclusion':'Exact Dirichlet domain ordering and parity reduction are correct. Executed negative finite levels and separate spacing/volume diagnostics support hydrogen-like quantized standing waves. Supplied mass/kinetic/coupling, numerical eigenvalue indexing and physical precision/particle attachment are explicit limitations.',
        'continuum_or_outward_rounded_enclosure_claimed':False,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}


if __name__=='__main__':
    output=run();Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
    print(output['verdict'])
