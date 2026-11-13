"""Native one-body charge-conjugation map and complete strong-polynomial symmetry."""
from itertools import product
import numpy as np
from .source import native,ORBITALS,source


def charge_conjugation():
    qg,ag,joins,qcharge,acharge=native()
    modes=tuple((z,*orb,c) for z,orb,c in product(range(2),ORBITALS,range(3)))
    index={v:i for i,v in enumerate(modes)}
    C=np.zeros((96,96),complex)
    for x,f,s,w in ORBITALS:
        for c,d in product(range(3),repeat=2):
            value=joins[w][c,d]
            C[index[1,x,f,s,w,d],index[0,x,f,s,w,c]]=value
            C[index[0,x,f,s,w,c],index[1,x,f,s,w,d]]=value.conjugate()
    assert np.array_equal(C.conj().T,C) and np.array_equal(C@C,np.eye(96))
    Q=np.diag([(qcharge if z==0 else acharge)[w][c] for z,x,f,s,w,c in modes])
    assert np.max(abs(C@Q@C+Q))<1e-15
    signs=[]
    for axis in range(8):
        G=np.zeros((96,96),complex)
        for z,x,f,s,w,c,d in product(range(2),range(2),range(2),range(2),range(2),range(3),range(3)):
            G[index[z,x,f,s,w,c],index[z,x,f,s,w,d]]=(qg if z==0 else ag)[w][axis][c,d]
        sign=1 if np.max(abs(qg[0][axis].T+qg[0][axis]))<1e-15 else -1
        assert np.max(abs(C@G@C-sign*G))<1e-15
        signs.append(sign)
    d=source('meson');lookup={v:i for i,v in enumerate(d['labels'])}
    permutation=[lookup[b,a] for a,b in d['labels']]
    H=d['H0'].toarray()
    assert np.array_equal(H,H[np.ix_(permutation,permutation)])
    return {'one_body_dimension':96,'C_squared_exact':True,
            'electric_charge_reversed':True,'colour_outer_automorphism_signs':signs,
            'complete_meson_dimension':256,'complete_meson_permutation_error':0,
            'light_source_transfer_error':d['light_transfer_error'],
            'polynomial_scope':'Strong action only: species-symmetric masses/kinetic/spin/pair terms; local net-baryon square; colour Casimirs; paired quark/antiquark flavour Casimirs minus c min(Nq,Na). The weak interaction is not asserted C invariant.'}


if __name__=='__main__':
    import json
    print(json.dumps(charge_conjugation()),flush=True)
