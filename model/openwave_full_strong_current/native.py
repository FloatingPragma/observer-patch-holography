"""Exact full light-sector CAR/Gauss action on finite-support Fock vectors.

Modes are (species, cell, spin, weak, canonical colour).  Canonical anti-colour
is the native signed join, explicitly authenticated when native words enter.
There is no finite-volume boundary: all hopping/correlated successors survive.
"""
from bisect import bisect_left
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,combinations_with_replacement,permutations,product
from math import factorial
import sympy as sy
from model.openwave_hadron_refinement_flow.positivity import inputs,declared_kappa
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_baryon_dispersion import exact_native as tri
from model.openwave_baryon_dispersion.source import repair_matrix
from model.openwave_native_composites.spectrum import native as native_generators


def clean(value):return {k:v for k,v in value.items()if v}
def add(out,value,scale=1):
    if scale:
        for key,coefficient in value.items():out[key]+=scale*coefficient


def annihilate(state,modes):
    state=list(state);sign=1
    for mode in sorted(modes):
        i=bisect_left(state,mode)
        if i==len(state)or state[i]!=mode:return (),0
        sign*=(-1)**i;state.pop(i)
    return tuple(state),sign


def create(state,modes):
    state=list(state);sign=1
    for mode in reversed(sorted(modes)):
        i=bisect_left(state,mode)
        if i<len(state)and state[i]==mode:return (),0
        sign*=(-1)**i;state.insert(i,mode)
    return tuple(state),sign


def one_body(state,destination,origin):
    middle,sign=annihilate(state,(origin,))
    if not sign:return {}
    final,other=create(middle,(destination,))
    return {final:sign*other}if other else {}


def action_vector(function,vector,*args):
    out=defaultdict(F)
    for state,value in vector.items():add(out,function(state,*args),value)
    return clean(out)


def translate(state,shift):return tuple((sp,x+shift,s,w,c)for sp,x,s,w,c in state)
def dot(left,right):return sum((value*right.get(key,0)for key,value in left.items()),F())


@lru_cache(None)
def angular(state,cells,spin=True,species=None):
    selected=[m for m in state if m[1]in cells and(species is None or m[0]==species)]
    axis=2 if spin else 3;z=sum(1-2*m[axis]for m in selected)
    out=defaultdict(F);out[state]+=F(z*(z+2),4)
    raised=defaultdict(F)
    for mode in selected:
        if mode[axis]==1:
            target=list(mode);target[axis]=0;add(raised,one_body(state,tuple(target),mode))
    for middle,value in raised.items():
        for mode in middle:
            if mode[1]in cells and(species is None or mode[0]==species)and mode[axis]==0:
                target=list(mode);target[axis]=1;add(out,one_body(middle,tuple(target),mode),value)
    return clean(out)


@lru_cache(None)
def colour_E(state,cells,a,b):
    out=defaultdict(F)
    for mode in state:
        species,x,s,w,c=mode
        if x not in cells:continue
        if species==0 and c==b:
            add(out,one_body(state,(species,x,s,w,a),mode))
        elif species==1 and c==a:
            add(out,one_body(state,(species,x,s,w,b),mode),-1)
    return clean(out)


@lru_cache(None)
def colour_casimir(state,cells):
    # Fundamental Fierz identity, with the full antifundamental generator.
    difference=sum((1 if sp==0 else-1)for sp,x,s,w,c in state if x in cells)
    out=defaultdict(F);out[state]-=F(difference*difference,6)
    for a,b in product(range(3),repeat=2):
        for middle,value in colour_E(state,cells,b,a).items():
            add(out,colour_E(middle,cells,a,b),value/2)
    return clean(out)


@lru_cache(None)
def pair_spin(state):
    # 3/2+2 s_q.s_a = 1+spin-swap for every q/anti pair at distance<=1.
    out=defaultdict(F)
    for q in state:
        if q[0]:continue
        for a in state:
            if not a[0]or abs(q[1]-a[1])>1:continue
            out[state]+=2 if q[2]==a[2]else 1
            if q[2]!=a[2]:
                qt=(q[0],q[1],a[2],q[3],q[4]);at=(a[0],a[1],q[2],a[3],a[4])
                for middle,value in one_body(state,at,a).items():add(out,one_body(middle,qt,q),value)
    return clean(out)


@lru_cache(None)
def coefficients():
    p,_,unit=inputs();dU,dn=degree_coefficients()
    return {**p,'U':p['U']+dU,'mbar':p['mbar']+dn,'kappa':declared_kappa(),'unit':unit}


def local_completion(state):
    lo=min(m[1]for m in state);hi=max(m[1]for m in state)
    counts={x:(sum(m[0]==0 and m[1]==x for m in state),sum(m[0]==1 and m[1]==x for m in state))for x in range(lo,hi+1)}
    value=sum(min(row)for row in counts.values())
    for x in range(lo,hi):
        a,b=counts[x],counts[x+1]
        value+=min(a[0]+b[0],a[1]+b[1])-min(a)-min(b)
    return value


@lru_cache(None)
def base_action(state):
    """Every original finite-range term, including electric and mixed pair terms."""
    if not state:return {}
    p=coefficients();out=defaultdict(F);N=len(state)
    lo,hi=min(m[1]for m in state),max(m[1]for m in state)
    sites=tuple(range(lo,hi+1))
    difference=sum(2*m[3]-1 for m in state)
    out[state]+=(p['mbar']+p['t'])*N+p['delta']*difference/2
    out[state]-=F(3,4)*(p['Jw']+p['eta'])*local_completion(state)
    for x in sites:
        d=sum((1 if sp==0 else-1)for sp,y,s,w,c in state if y==x)
        out[state]+=p['A']*d*d/9
        add(out,colour_casimir(state,(x,)),p['U']-p['kappa']/2)
        add(out,angular(state,(x,)),p['Js']-p['spin_gradient'])
        for species in (0,1):add(out,angular(state,(x,),False,species),p['Jw']-p['eta'])
    # Zero boundary flux on the global singlet carrier; all intervening edges
    # retain their complete prefix Casimir. No string energy is discarded.
    for edge in range(lo,hi):
        add(out,colour_casimir(state,tuple(range(lo,edge+1))),p['kappa'])
    for x in range(lo-1,hi+1):
        cells=(x,x+1)
        for species in(0,1):add(out,angular(state,cells,False,species),p['eta'])
        for cell in cells:add(out,angular(state,(cell,)),2*p['spin_gradient'])
        add(out,angular(state,cells),-p['spin_gradient'])
    for mode in state:
        sp,x,s,w,c=mode
        for y in(x-1,x+1):add(out,one_body(state,(sp,y,s,w,c),mode),-p['t'])
    add(out,pair_spin(state),p['pair_triplet'])
    # The inherited three-Nx=3 nuclear term vanishes identically at N=5.
    return clean(out)


def epsilon(shape):
    out=defaultdict(F);divisor=1
    for number in Counter(shape).values():divisor*=factorial(number)
    for colours in permutations(range(3)):
        sign=(-1)**sum(colours[i]>colours[j]for i in range(3)for j in range(i+1,3))
        modes=[(0,x,s,w,c)for(x,s,w),c in zip(shape,colours)]
        if len(set(modes))<3:continue
        order=(-1)**sum(modes[i]>modes[j]for i in range(3)for j in range(i+1,3))
        out[tuple(sorted(modes))]+=F(sign*order,divisor)
    out=clean(out);assert all(value in(-1,1)for value in out.values())
    return out


def pair_column(shape):
    (x,s,w),(y,t,z)=shape
    return {((0,x,s,w,c),(1,y,t,z,c)):F(1)for c in range(3)}


@lru_cache(None)
def triquark_frame():
    labels=tuple(combinations_with_replacement(range(16),3))
    shapes=sorted({tri.canonical(lab)[1]for lab in labels if max(tri.ORBITALS[i][0]for i in lab)-min(tri.ORBITALS[i][0]for i in lab)<=1})
    columns=tuple(epsilon(shape)for shape in shapes);A,G,_=tri.complete()
    assert len(columns)==100 and all(dot(col,col)==G[j,j]for j,col in enumerate(columns))
    return shapes,columns,A,G


@lru_cache(None)
def meson_frame():
    shapes=tuple(((x,s,w),(y,t,z))for x,y in((0,0),(0,1),(1,0))for s,w,t,z in product(range(2),repeat=4))
    columns=tuple(pair_column(shape)for shape in shapes)
    reverse={state:(i,value)for i,col in enumerate(columns)for state,value in col.items()}
    A=sy.zeros(48)
    for j,col in enumerate(columns):
        for shift in(-1,0,1):
            value=action_vector(base_action,{translate(state,shift):c for state,c in col.items()})
            for state,c in value.items():
                if state in reverse:
                    i,b=reverse[state];A[i,j]+=sy.Rational(b*c)
    G=3*sy.eye(48);assert A==A.T
    return shapes,columns,A,G


def inverse_by_blocks(A):
    remaining=set(range(A.rows));out=sy.zeros(A.rows);sizes=[]
    while remaining:
        component={min(remaining)};front=list(component)
        while front:
            i=front.pop()
            for j in tuple(remaining-component):
                if A[i,j]:component.add(j);front.append(j)
        ids=sorted(component);remaining-=component;sizes.append(len(ids))
        block=A.extract(ids,ids).inv()
        for i,a in enumerate(ids):
            for j,b in enumerate(ids):out[a,b]=block[i,j]
    assert A*out==sy.eye(A.rows)
    return out,tuple(sizes)


@lru_cache(None)
def correlated_data(kind):
    shapes,columns,A,G=triquark_frame()if kind=='triquark'else meson_frame()
    inverse,sizes=inverse_by_blocks(A);zeta=repair_matrix()[2]
    J=sy.Rational(zeta)*inverse
    reverse={state:(j,value)for j,col in enumerate(columns)for state,value in col.items()}
    rows=tuple(tuple((i,F(J[i,j]))for i in range(J.rows)if J[i,j])for j in range(J.cols))
    return columns,reverse,rows,A,G,sizes


def correlated_action(state,kind):
    columns,reverse,rows,_,_,_=correlated_data(kind)
    out=defaultdict(F)
    subsets=(combinations(tuple(m for m in state if m[0]==0),3)if kind=='triquark'
             else ((q,a)for q in state if q[0]==0 for a in state if a[0]==1))
    for subset in subsets:
        subset=tuple(sorted(subset));origin=min(m[1]for m in subset)
        if max(m[1]for m in subset)-origin>1:continue
        local=translate(subset,-origin)
        if local not in reverse:continue
        j,amplitude=reverse[local];remaining,sgn=annihilate(state,subset);assert sgn
        for i,coefficient in rows[j]:
            for displacement,laplace in((0,2),(-1,-1),(1,-1)):
                for created,b in columns[i].items():
                    final,other=create(remaining,translate(created,origin+displacement))
                    if other:out[final]+=amplitude*sgn*coefficient*laplace*b*other
    return clean(out)


def full_action(state):
    out=defaultdict(F);add(out,base_action(state))
    add(out,correlated_action(state,'triquark'))
    add(out,correlated_action(state,'meson'))
    return clean(out)


def from_native_mask(mask):
    joins=native_generators()[2]
    assert all((join==[[1,0,0],[0,-1,0],[0,0,1]]).all()for join in joins)
    rows=[];phase=1
    for i in range(96):
        if not(mask>>i&1):continue
        species=int(i>=48);orbital,colour=divmod(i%48,3)
        x,family,s,w=orbital//8,(orbital//4)%2,(orbital//2)%2,orbital%2
        assert family==0
        rows.append((species,x,s,w,colour))
        if species and colour==1:phase=-phase
    assert rows==sorted(rows)
    return tuple(rows),phase


def native_vector(value):
    out=defaultdict(F)
    for mask,c in value.items():
        state,phase=from_native_mask(mask);out[state]+=F(c)*phase
    return clean(out)
