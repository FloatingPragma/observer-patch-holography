"""Integer occupation-coordinate contraction of every native triquark channel.

Each colour-epsilon column is the unnormalized sum of its distinct orbital
permutations. Its Gram weight is 3!/prod(n_i!). E_ij has coordinate n_i+1
for i!=j and n_i on the diagonal, so the Gram numerator is the integer
identity g_new*(n_i+1)=g_old*n_j. This is a combinatorial CAR contraction;
no floating contraction or rounding supplies its entries.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import factorial
import sympy as sy
from model.openwave_hadron_refinement_flow.positivity import inputs,declared_kappa
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients

ORBITALS=tuple(product(range(4),range(2),range(2)))
INDEX={v:i for i,v in enumerate(ORBITALS)}


def gram(label):
    value=6
    for n in Counter(label).values():value//=factorial(n)
    return value


def canonical(label):
    origin=min(ORBITALS[i][0]for i in label)
    return origin,tuple((ORBITALS[i][0]-origin,*ORBITALS[i][1:])for i in label)


def one_body(label,i,j):
    count=Counter(label)
    if not count[j]:return {}
    if i==j:return {label:count[j]}
    coefficient=count[i]+1
    out=list(label);out.remove(j);out.append(i);target=tuple(sorted(out))
    assert gram(target)*coefficient==gram(label)*count[j]
    return {target:coefficient}


def add(out,term,coefficient=1):
    for key,value in term.items():out[key]+=coefficient*value


@lru_cache(None)
def ladder(label,cells,spin,raising):
    out=defaultdict(int)
    for x,other in product(cells,range(2)):
        a=(x,0,other)if spin else(x,other,0)
        b=(x,1,other)if spin else(x,other,1)
        add(out,one_body(label,INDEX[a if raising else b],INDEX[b if raising else a]))
    return dict(out)


@lru_cache(None)
def angular4(label,cells,spin):
    axis=1 if spin else 2
    z=sum(1-2*ORBITALS[i][axis]for i in label if ORBITALS[i][0]in cells)
    out=defaultdict(int);out[label]=z*(z+2)
    for middle,value in ladder(label,cells,spin,True).items():
        add(out,ladder(middle,cells,spin,False),4*value)
    return dict(out)


def hamiltonian(label):
    p,_,_=inputs();k=declared_kappa();dU,dn=degree_coefficients()
    n=[sum(ORBITALS[i][0]==x for i in label)for x in range(4)]
    nu=sum(ORBITALS[i][2]==0 for i in label);prefix=[sum(n[:x])for x in range(1,4)]
    out=defaultdict(F)
    out[label]=3*p['mbar']+p['delta']/2*(3-2*nu)+p['A']/9*sum(v*v for v in n)
    out[label]+=(p['U']-k/2)*F(2,3)*sum(v*(3-v)for v in n)+k*F(2,3)*sum(v*(3-v)for v in prefix)
    for x in range(4):
        add(out,angular4(label,(x,),True),p['Js']/4)
        add(out,angular4(label,(x,),False),p['Jw']/4)
    for x in range(3):
        for s,w in product(range(2),repeat=2):
            i,j=INDEX[x,s,w],INDEX[x+1,s,w]
            for a,b,c in ((i,i,1),(j,j,1),(i,j,-1),(j,i,-1)):add(out,one_body(label,a,b),c*p['t'])
        add(out,angular4(label,(x,x+1),False),p['eta']/4)
        for cell in (x,x+1):add(out,angular4(label,(cell,),True),p['spin_gradient']/2)
        add(out,angular4(label,(x,x+1),True),-p['spin_gradient']/4)
    for x in (1,2):
        add(out,angular4(label,(x,),True),-p['spin_gradient']/4)
        add(out,angular4(label,(x,),False),-p['eta']/4)
        out[label]+=(-p['t']+dn)*n[x]+dU*F(2,3)*n[x]*(3-n[x])
    return dict(out)


@lru_cache(None)
def complete():
    labels=tuple(combinations_with_replacement(range(16),3))
    selected=[lab for lab in labels if max(ORBITALS[i][0]for i in lab)-min(ORBITALS[i][0]for i in lab)<=1]
    shapes=sorted({canonical(lab)[1]for lab in selected});assert len(shapes)==100
    lookup={shape:i for i,shape in enumerate(shapes)}
    rows={lab:lookup[shape]for lab in selected for origin,shape in [canonical(lab)]if origin==1}
    A=sy.zeros(100)
    for label in selected:
        j=lookup[canonical(label)[1]]
        for target,value in hamiltonian(label).items():
            if target in rows:A[rows[target],j]+=sy.Rational(value)*gram(target)
    G=sy.diag(*[gram(next(lab for lab,i in rows.items()if i==j))for j in range(100)])
    assert A==A.T
    angular=[]
    for spin in (True,False):
        S=sy.zeros(100)
        for label,j in rows.items():
            for target,value in angular4(label,(0,1,2,3),spin).items():
                if target in rows:S[rows[target],j]+=sy.Rational(gram(target)*value,4)
        assert S==S.T and A*G.inv()*S==S*G.inv()*A
        angular.append(S)
    return A,G,tuple(angular)
