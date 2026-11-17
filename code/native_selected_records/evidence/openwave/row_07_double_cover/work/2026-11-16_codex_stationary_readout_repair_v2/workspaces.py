"""Finite native-event unitary extensions and executable V,V,V† microtraces."""
from itertools import product
import hashlib
import numpy as np
from native_source import native
from amplified_native import require

VALUES=tuple(range(-2,4))
LOADS=tuple(product(VALUES,repeat=3))
INDEX={q:i for i,q in enumerate(LOADS)}
MODULUS=11


def pack(q,c,j):return (INDEX[tuple(q)]*2+c)*MODULUS+(j%MODULUS)


def unpack(index):
    qc,j=divmod(int(index),MODULUS);q,c=divmod(qc,2)
    return LOADS[q],c,j if j<=5 else j-MODULUS


def extension(edge,law):
    size=len(LOADS)*2*MODULUS
    permutation=np.full(size,-1,dtype=int);used=set();records=[]
    for q in LOADS:
        for c in (0,1):
            out,event=native(q,edge,c,law)
            i,o=pack(q,c,0),pack(out,c,event[3])
            require(o not in used,'Native isometry collision')
            permutation[i]=o;used.add(o)
            records.append({'input':i,'output':o,'loads':list(q),'event':event})
    missing={}
    for i in range(size):
        if i not in used:
            q,c,_=unpack(i);missing.setdefault((sum(q),c),[]).append(i)
    missing={sector:iter(indices)for sector,indices in missing.items()}
    completion=[]
    for i in range(size):
        if permutation[i]<0:
            q,c,_=unpack(i)
            permutation[i]=next(missing[sum(q),c]);completion.append([i,int(permutation[i])])
    require(np.array_equal(np.sort(permutation),np.arange(size)),'Unitary extension is not bijective')
    require(all((sum(unpack(i)[0]),unpack(i)[1])==(sum(unpack(permutation[i])[0]),unpack(permutation[i])[1])for i in range(size)),
        'Native extension changes total load or source coin')
    inverse=np.argsort(permutation)
    return permutation,inverse,records,completion


def audit():
    cases=[];extensions=[]
    for law in ('native','swap'):
        for edge in ((0,1),(1,2),(0,2)):
            v,vi,records,completion=extension(edge,law)
            domain=np.arange(len(v)*MODULUS)
            qcj1,j2=np.divmod(domain,MODULUS)
            qc,j1=np.divmod(v[qcj1],MODULUS)
            qc=qc^((j1!=0).astype(int))
            qc,j2=np.divmod(v[qc*MODULUS+j2],MODULUS)
            j1=(j1+j2)%MODULUS
            qc,j2=np.divmod(vi[qc*MODULUS+j2],MODULUS)
            full=(qc*MODULUS+j1)*MODULUS+j2
            require(np.array_equal(np.sort(full),domain),'Complete two-workspace circuit is not a permutation')
            extensions.append({'law':law,'edge':list(edge),'dimension':len(v),
                'native_image':records,'orthogonal_complement_permutation':completion,
                'completion_rule':'Sorted unused labels within each fixed total-load and source-coin sector.',
                'full_workspace_dimension':len(full),
                'full_workspace_permutation_sha256':hashlib.sha256(full.astype('<i8').tobytes()).hexdigest()})
            for background in range(-1,3):
                for sign,position in product((-1,1),range(3)):
                    q=tuple(background+sign*int(i==position)for i in range(3))
                    for coin in (0,1):
                        first=v[pack(q,coin,0)];q1,c1,j1=unpack(first)
                        cc=c1^int(j1!=0)
                        second=v[pack(q1,cc,0)];q2,c2,j2=unpack(second)
                        cleared=(j1+j2)%MODULUS
                        final=vi[second];q3,c3,scratch=unpack(final)
                        require(cleared==0 and scratch==0 and q3==q1 and c3==cc,'Workspace microtrace failed')
                        cases.append({'law':law,'edge':list(edge),'background':background,
                            'sign':sign,'input':list(q),'coin':coin,'after_first':[list(q1),c1,j1],
                            'after_second':[list(q2),c2,j2],'first_current_after_add':cleared,
                            'after_inverse':[list(q3),c3,scratch],
                            'forward_accesses':2,'inverse_accesses':1})
    return {'domain':'All signed local loads in [-2,3]^3; signed currents [-5,5] are faithfully encoded modulo11 in the explicit finite extension.',
        'admitted_reversal_sector':'background±one excitation; common source uses background0 after its entry guard.',
        'inverse_access':'Declared coherent V† on the entire displayed finite unitary extension.',
        'extensions':extensions,'microtraces':cases,'count':len(cases)}
