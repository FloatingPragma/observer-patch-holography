"""Gauss-eliminated open colour chains with a fixed two-cell attachment."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
import numpy as np
import scipy.sparse as sp
import sympy as sy
from scipy.linalg import eigh
from model.openwave_nuclear_distance.spatial import carrier
from model.openwave_spin_gradient_nuclei import source as native

@lru_cache(None)
def chain(L):
    assert L>=2
    d=carrier(L);p=native.parameters();U=p['U']-p['kappa_colour']/2;assert U>0
    n=d['n'];N=len(n);E=d['E'];idx=d['index'];angular=d['angular'];counts=d['counts']
    K=sp.csr_matrix((N,N));S=K.copy();I=K.copy();edgeI=K.copy();gradient=K.copy()
    for x in range(L):S+=4*angular((x,),True);I+=4*angular((x,),False)
    for x in range(L-1):
        for ss,f in product(range(2),repeat=2):
            a,b=idx[x,ss,f],idx[x+1,ss,f];K+=E(a,a)+E(b,b)-E(a,b)-E(b,a)
        edgeI+=4*angular((x,x+1),False)
        gradient+=8*(angular((x,),True)+angular((x+1,),True))-4*angular((x,x+1),True)
    prefixes=np.cumsum(n,axis=1)[:,:-1]
    C3=2*np.sum(n*(3-n),axis=1);P3=2*np.sum(prefixes*(3-prefixes),axis=1)
    spans=np.array([max(d['orbitals'][j][0]for j in lab)-min(d['orbitals'][j][0]for j in lab)for lab in d['labels']])
    assert np.array_equal(P3,4*spans)
    up=counts[:,[j for j,o in enumerate(d['orbitals'])if o[2]==0]].sum(axis=1)
    spinup=counts[:,[j for j,o in enumerate(d['orbitals'])if o[1]==0]].sum(axis=1)
    parts={'N':3*sp.eye(N,dtype=np.int64),'D':sp.diags(3-2*up,dtype=np.int64),'B9':sp.diags(np.sum(n*n,axis=1),dtype=np.int64),
      'C3':sp.diags(C3,dtype=np.int64),'P3':sp.diags(P3,dtype=np.int64),'K':K,'S4':S,'I4':I,'EdgeI4':edgeI,'Gradient4':gradient}
    weights={'N':p['mbar'],'D':p['delta']/2,'B9':p['A']/9,'C3':U/3,'P3':p['kappa_colour']/3,
      'K':p['t'],'S4':p['Js']/4,'I4':p['Jw']/4,'EdgeI4':p['eta']/4,'Gradient4':p['spin_gradient']/4}
    H=sum(float(weights[k])*a for k,a in parts.items()).tocsr();ids=np.flatnonzero((up==1)&(spinup==2))
    assert np.max(abs((H-H.T).data),initial=0)<1e-12
    gram=np.array([6//np.prod([factorial(int(v))for v in row])for row in counts],dtype=int)
    reduced={}
    for name,a in parts.items():
        small=a.tocsr()[ids][:,ids].toarray()*np.sqrt(gram[ids,None]*gram[None,ids])
        assert np.max(abs(small-np.rint(small)),initial=0)<1e-10
        reduced[name]=np.rint(small).astype(np.int64)
    exact=sy.zeros(len(ids))
    for name,w in weights.items():exact+=sy.Rational(w.numerator,w.denominator)*sy.Matrix(reduced[name].tolist())
    g=gram[ids];exact=sy.Matrix(len(ids),len(ids),lambda i,j:exact[i,j]/sy.sqrt(int(g[i])*int(g[j])))
    assert np.max(abs(np.array(exact.evalf(18),float)-H[ids][:,ids].toarray()))<2e-11
    if L==2:
        old=native.sector('baryon');assert tuple(old['basis'])==d['full']
        overlap=old['B'].T@d['B'];difference=old['h']@overlap-overlap@H
        assert np.max(abs(difference))<3e-10
    # Uniform position-dependent cost lower bound on this entire N3 source.
    b=p['A']/9-F(2,3)*U;assert b>0
    floor=3*(p['mbar']-p['delta']/2)+6*U+3*b
    return {**d,'parts':parts,'weights':weights,'H':H,'ids':ids,'exact':exact,'spans':spans,'floor':floor,'adjusted_local_colour':U,'parameters':p}

@lru_cache(None)
def native_transfer():
    from model.openwave_six_quark.source import colour48,exact_product
    from model.openwave_nuclear_distance.nuclear import carrier as nuclear_carrier
    b=native.components('baryon');c0=colour48(b['basis'],(0,)).tocsr();c1=colour48(b['basis'],(1,)).tocsr()
    assert (c0+c1-b['parts']['C48']).nnz==0
    assert exact_product(c0-c1,b['Binteger']).nnz==0
    assert exact_product(b['parts']['Cg48'],b['Binteger']).nnz==0
    d=nuclear_carrier();a0=colour48(d['src'],(0,)).tocsr();a1=colour48(d['src'],(1,)).tocsr();ids=d['keep']
    assert (a0+a1)[ids][:,ids].shape==d['parts']['C48'].shape
    assert ((a0+a1)[ids][:,ids]-d['parts']['C48']).nnz==0
    difference=(a0-a1)[ids][:,ids]
    assert exact_product(difference,d['P']).nnz==0
    assert exact_product(d['parts']['Cg48'],d['P']).nnz==0
    # The same identities are false on the retained full coloured complement.
    assert (c0-c1).nnz>0 and b['parts']['Cg48'].nnz>0
    return {'all_native_baryon_singlets':120,'full_native_baryon_occupation_dimension':len(b['basis']),
      'exact_nuclear_Ms0_Mi0_colour_projector_rank':d['colour_weight_dimension'],
      'all_nuclear_magnetic_singlets':2520,'full_native_nuclear_Ms0_Mi0_dimension':len(d['src']),
      'baryon_local_colour_difference_times_exact_basis_nonzeros':0,
      'nuclear_local_colour_difference_times_exact_complete_projector_nonzeros':0,
      'global_colour_times_both_exact_projectors_nonzeros':0,
      'coloured_complement_identity_fails_and_is_retained':True,
      'spin_isospin_multiplet_extension':'Local colour operators are spin/isospin scalars. Exact zero on the complete Ms0/Mi0 colour projector extends to every magnetic partner.'}
