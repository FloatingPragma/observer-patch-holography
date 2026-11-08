"""Exact-source, literal-control interval replay of the admitted beam.

Every integrated port-control binary64 payload is read as an exact rational,
matching certified_reference.source. All six processor coordinates and every
complement row are enclosed before extracting the successful local block.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,isqrt
from model.certified_reference import source
from model.native_effect_synthesis import intervals as iv
from model.native_fair_local_preparation import coupling
from .literal_angles import ZERO,IDENTITY,add,multiply


def product(a,b):
    out=[[ZERO for _ in range(len(b[0]))]for _ in a]
    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k]==ZERO:continue
            for j in range(len(b[0])):
                if b[k][j]!=ZERO:out[i][j]=add(out[i][j],multiply(a[i][k],b[k][j]))
    return out


def divide(z,n):return iv.divide(z[0],n),iv.divide(z[1],n)
def plus_matrices(a,b):return [[add(x,y)for x,y in zip(r,s,strict=True)]for r,s in zip(a,b,strict=True)]


@lru_cache(maxsize=64)
def primitive(kind,angle):
    controls=source.literal_controls(kind,float(angle));spec=source.specification()
    root=isqrt(5*iv.Q*iv.Q);sqrt5=(root,root+1)
    generator=[]
    for i in range(6):
        row=[]
        for j in range(6):
            components=[]
            for part in (0,1):
                a=b=F()
                for c,g in zip(controls,spec,strict=True):
                    x,y=g[i][j][part];a+=c*x;b+=c*y
                components.append(iv.add(iv.exact(a),iv.mul(iv.exact(b),sqrt5)))
            row.append(tuple(components))
        generator.append(row)
    bound=max(sum(F(max(abs(z[0]),abs(z[1])),iv.Q)for cell in row for z in cell)for row in generator)
    B=(bound.numerator+bound.denominator-1)//bound.denominator
    N=96
    if not 0<B<N+2:raise ArithmeticError('literal source exponential exceeds finite Taylor domain')
    remainder=F(B**(N+1),factorial(N+1))/(1-F(B,N+2))
    radius=(remainder.numerator*iv.Q+remainder.denominator-1)//remainder.denominator
    term=[[IDENTITY if i==j else ZERO for j in range(6)]for i in range(6)]
    total=term
    for n in range(1,N+1):
        term=[[divide(z,n)for z in row]for row in product(generator,term)]
        total=plus_matrices(total,term)
    return tuple(tuple(tuple((a-radius,b+radius)for a,b in z)for z in row)for row in total)


def lifted_success(event):
    full=primitive(event['kind'],F(event['angle']));targets=event['targets'];bits=len(targets)
    out=[[ZERO for _ in range(4)]for _ in range(4)]
    for i in range(4):
        for j in range(4):
            if any(((i>>b)&1)!=((j>>b)&1)for b in range(2)if b not in targets):continue
            ii=sum(((i>>b)&1)<<(bits-1-k)for k,b in enumerate(targets))
            jj=sum(((j>>b)&1)<<(bits-1-k)for k,b in enumerate(targets))
            out[i][j]=full[ii][jj]
    return out


def certificate():
    word=coupling.words()['beam'][1]
    output=[[IDENTITY if i==j else ZERO for j in range(4)]for i in range(4)]
    complement_rows=0
    for event in word:
        output=product(lifted_success(event),output)
        complement_rows+=6-(4 if event['kind']=='J' else 2)
    entries=[]
    for i in (0,3):
        for j in (1,2):
            for k,name in enumerate(('real','imag')):
                lo,hi=output[i][j][k]
                if lo>0 or hi<0:
                    entries.append({'output':i,'input':j,'component':name,
                                    'interval':[str(F(lo,iv.Q)),str(F(hi,iv.Q))]})
    if not entries:raise ArithmeticError('actual literal source-control leakage was not separated from zero')
    return {'pulses':len(word),'unique_literal_processor_responses':primitive.cache_info().currsize,
            'processor_dimension':6,'complete_complement_rows':complement_rows,
            'matrix_exponential_taylor_order':96,'interval_bits':iv.BITS,
            'source_semantics':'exact Q(sqrt5) generators, exact binary64 integrated control rationals',
            'successful_work_number_leakage_entries':entries,
            'complement_outputs_projected_away':False,
            'symbolic_pi_or_ideal_generator_substitution':False}
