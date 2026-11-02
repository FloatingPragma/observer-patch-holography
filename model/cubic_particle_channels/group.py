"""Exact binary octahedral rotations, including their hidden-frame action."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import numpy as np

from model.field_symmetry import symmetry
from model.cube_eta_protection import cube
from model.centered_particle_spectrum import source
from model.interacting_vacuum import ground


@dataclass(frozen=True)
class Q:
    """An exact element a+b sqrt(2); never silently accept a float."""
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        if type(self.a) not in (int, F) or type(self.b) not in (int, F):
            raise ValueError('rational coefficients required')
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def coerce(x):
        if type(x) is Q: return x
        if type(x) in (int, F): return Q(x)
        raise ValueError('exact coefficient required')

    def __add__(self, x):
        x=self.coerce(x); return Q(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,x): return self+-self.coerce(x)
    def __rsub__(self,x): return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x); return Q(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        if type(x) not in (int,F) or not x: raise ValueError('nonzero rational divisor required')
        return Q(self.a/x,self.b/x)
    def __bool__(self): return bool(self.a or self.b)
    def rational(self):
        if self.b: raise ArithmeticError('uncancelled radical')
        return self.a
    def integer(self):
        x=self.rational()
        if x.denominator!=1: raise ArithmeticError('nonintegral coefficient')
        return x.numerator
    def json(self): return [str(self.a),str(self.b)]


ONE=(Q(1),Q(),Q(),Q())

def multiply(q,r):
    w,x,y,z=q; a,b,c,d=r
    return (w*a-x*b-y*c-z*d,w*b+x*a+y*d-z*c,
            w*c-x*d+y*a+z*b,w*d+x*c-y*b+z*a)

def conjugate(q): return (q[0],-q[1],-q[2],-q[3])
def negative(q): return tuple(-x for x in q)

@lru_cache(maxsize=1)
def elements():
    rows=[]
    for j in range(4):
        for s in (-1,1): rows.append(tuple(Q(s if k==j else 0) for k in range(4)))
    rows.extend(tuple(Q(F(s,2)) for s in signs) for signs in product((-1,1),repeat=4))
    for pair in combinations(range(4),2):
        for signs in product((-1,1),repeat=2):
            row=[Q() for _ in range(4)]
            for j,s in zip(pair,signs,strict=True):row[j]=Q(0,F(s,2))
            rows.append(tuple(row))
    if len(set(rows))!=48:raise ArithmeticError('binary rotation group coverage')
    return tuple(sorted(rows,key=lambda q:tuple((x.a,x.b) for x in q)))

@lru_cache(maxsize=48)
def rotation(q):
    columns=[]
    for j in range(3):
        v=tuple(Q(1 if k==j+1 else 0) for k in range(4))
        r=multiply(multiply(q,v),conjugate(q))
        if r[0]:raise ArithmeticError('nonvector quaternion rotation')
        columns.append([x.integer() for x in r[1:]])
    return np.array(columns,dtype=np.int64).T

@lru_cache(maxsize=48)
def permutation(q):
    sites=cube.GEOMETRY.sites;lookup={x:j for j,x in enumerate(sites)};R=rotation(q)
    return tuple(lookup[tuple(int(v) for v in ((R@(2*np.array(x)-1)+1)//2))] for x in sites)

def frame(site):
    q=ONE
    for j,n in enumerate(site):
        generator=tuple(Q(-1 if k==j+1 else 0) for k in range(4))
        for _ in range(n%4):q=multiply(q,generator)
    return q

def spinor_parts(q):
    """Real and imaginary parts of U(q)=w I-i(x sigma_x+y sigma_y+z sigma_z)."""
    w,x,y,z=q
    return ((w,-y),(y,w)), ((-z,-x),(-x,z))

def determinant3(R):
    return int(R[0,0])*(int(R[1,1])*int(R[2,2])-int(R[1,2])*int(R[2,1]))-int(R[0,1])*(int(R[1,0])*int(R[2,2])-int(R[1,2])*int(R[2,0]))+int(R[0,2])*(int(R[1,0])*int(R[2,1])-int(R[1,1])*int(R[2,0]))

@lru_cache(maxsize=48)
def hidden(q):
    p=permutation(q);sites=cube.GEOMETRY.sites
    blocks=[multiply(multiply(frame(sites[p[x]]),q),conjugate(frame(sites[x]))) for x in range(8)]
    W=blocks[0];sign=[]
    for block in blocks:
        if block==W:sign.append(1)
        elif block==negative(W):sign.append(-1)
        else:raise ArithmeticError('spatial rotation is not a common hidden spin rotation')
    S=np.zeros((8,8),dtype=np.int64);S[p,np.arange(8)]=sign
    return S,W

def parity(permutation):
    return (-1)**sum(x>y for i,x in enumerate(permutation) for y in permutation[i+1:])

def diagonal_sign(q):
    diagonals=tuple((1,y,z) for y,z in product((-1,1),repeat=2));lookup={x:i for i,x in enumerate(diagonals)}
    indices=[]
    for x in diagonals:
        y=rotation(q)@x
        if y[0]<0:y=-y
        indices.append(lookup[tuple(int(v) for v in y)])
    return parity(indices)

IRREPS=('A1','A2','E','T1','T2','G1','G2','H')
DIMENSIONS=(1,1,2,3,3,2,2,4)

@lru_cache(maxsize=48)
def characters(q):
    R=rotation(q);t=int(np.trace(R));a=diagonal_sign(q);w=q[0]
    e=(t*t+int(np.trace(R@R)))//2-1-a*t
    return dict(zip(IRREPS,(Q(1),Q(a),Q(e),Q(t),Q(a*t),2*w,2*a*w,8*w*w*w-4*w),strict=True))

def decompose(values):
    out={}
    for name in IRREPS:
        m=(sum((characters(q)[name]*values(q) for q in elements()),Q())/48).integer()
        if m<0:raise ArithmeticError('negative representation multiplicity')
        if m:out[name]=m
    return out

def spin_character(q,twice_spin):
    if type(twice_spin) is not int or twice_spin<0:raise ValueError('nonnegative twice-spin required')
    a,b=Q(1),2*q[0]
    if twice_spin==0:return a
    for _ in range(1,twice_spin):a,b=b,2*q[0]*b-a
    return b

@lru_cache(maxsize=4)
def walsh_projector_numerator(degree):
    if type(degree) is not int or degree not in range(4):raise ValueError('Walsh degree0..3 required')
    chars=np.array([[(-1)**((k&x).bit_count()) for x in range(8)]
                    for k in range(8) if k.bit_count()==degree],dtype=np.int64)
    return chars.T@chars

def spatial_character(q,degree):
    p=permutation(q);N=walsh_projector_numerator(degree)
    return Q(F(sum(int(N[x,p[x]]) for x in range(8)),8))

@lru_cache(maxsize=1)
def certificate():
    rows=elements();lookup={q:i for i,q in enumerate(rows)}
    for q in rows:
        if multiply(q,conjugate(q))!=ONE:raise ArithmeticError('nonunit rotation')
        R=rotation(q)
        if not np.array_equal(R.T@R,np.eye(3,dtype=np.int64)) or determinant3(R)!=1:
            raise ArithmeticError('improper spatial rotation')
        for r in rows:
            qr=multiply(q,r)
            if qr not in lookup or not np.array_equal(rotation(qr),R@rotation(r)):
                raise ArithmeticError('rotation multiplication fails')
    for a in IRREPS:
        for b in IRREPS:
            product_sum=sum((characters(q)[a]*characters(q)[b] for q in rows),Q())
            if product_sum!=Q(48 if a==b else 0):raise ArithmeticError('irreducible character orthogonality')
    if sum(d*d for d in DIMENSIONS)!=48:raise ArithmeticError('incomplete irreducible list')
    A=np.asarray(ground.band_data()[0],dtype=np.int64)
    for site in cube.GEOMETRY.sites:
        real,imag=spinor_parts(frame(site));actual=symmetry.onsite_frame(site)
        if any(real[i][j].integer()!=actual[i,j].real or imag[i][j].integer()!=actual[i,j].imag
               for i in range(2) for j in range(2)):
            raise ArithmeticError('quaternion frame differs from actual source')
    source_autos={tuple(int(x) for x in p):(s,axes,shift) for axes,shift,p,s in source.automorphisms()}
    for q in rows:
        S,W=hidden(q);p=permutation(q)
        if not np.array_equal(S@A,A@S):raise ArithmeticError('source hopping covariance fails')
        signs,_,_=source_autos[p]
        if not np.array_equal(S[p,np.arange(8)],signs):raise ArithmeticError('wrong signed source frame')
        if 2*W[0]*int(np.trace(A@S))!=Q():raise ArithmeticError('charged band characters differ')
    result={'binary_group_order':48,'proper_spatial_group_order':len({tuple(rotation(q).flat) for q in rows}),
      'full_Fock_central_rotation':'Gamma(-I)=(-1)^N on all65536 states',
      'source_hopping_and_contact_covariant':True,'character_orthogonality_exact':True,
      'ground_character':'A1',
      'hidden_spin_vector':decompose(lambda q:Q(int(np.trace(rotation(hidden(q)[1]))))),
      'charged_one_particle':decompose(lambda q:2*q[0]*sum(x==y for x,y in enumerate(permutation(q)))),
      'charged_each_visible_line':decompose(lambda q:Q(int(np.trace(hidden(q)[0])))*hidden(q)[1][0]),
      'density_by_walsh_degree':{str(k):decompose(lambda q,k=k:spatial_character(q,k)) for k in range(4)},
      'spin_by_walsh_degree':{str(k):decompose(lambda q,k=k:spatial_character(q,k)*int(np.trace(rotation(hidden(q)[1])))) for k in range(4)},
      'SU2_restriction_by_twice_spin':{str(n):decompose(lambda q,n=n:spin_character(q,n)) for n in range(9)}}
    return result
