"""Exact full-graph compression and complete vacuum energy measure."""
from collections import defaultdict,Counter
from functools import lru_cache
from fractions import Fraction
from math import isqrt
from pathlib import Path
import json
import numpy as np
import sympy as s
from .model import graph
PARTITION=Path(__file__).with_name('equitable_partition.npz')


def generate_partition():
    keys,h,phase,depth,digest=graph();diagonal=h.diagonal()
    colors=[2 if key==(0,0) else int(diagonal[i]) for i,key in enumerate(keys)]
    for _ in range(len(keys)):
        mapping={};nextcolors=[]
        for i in range(h.shape[0]):
            totals=defaultdict(int)
            for offset in range(h.indptr[i],h.indptr[i+1]):
                j=int(h.indices[offset]);totals[colors[j]]+=int(h.data[offset])*int(phase[i])*int(phase[j])
            signature=(colors[i],tuple(sorted(totals.items())))
            if signature not in mapping:mapping[signature]=len(mapping)
            nextcolors.append(mapping[signature])
        if len(set(colors))==len(mapping):break
        colors=nextcolors
    else:raise ValueError('finite partition did not stabilize')
    np.savez_compressed(PARTITION,colors=np.array(nextcolors,dtype=np.int32),source_digest=np.array(digest))


@lru_cache(maxsize=1)
def quotient():
    keys,h,phase,depth,digest=graph()
    with np.load(PARTITION,allow_pickle=False) as archive:
        if str(archive['source_digest'])!=digest:raise ValueError('partition belongs to another source')
        colors=archive['colors'].copy()
    if colors.dtype!=np.int32 or colors.shape!=(len(keys),):raise ValueError('full partition array required')
    classes=sorted(set(map(int,colors)))
    if classes!=list(range(len(classes))):raise ValueError('contiguous nonempty partition required')
    rows={};sizes=Counter(map(int,colors))
    for i in range(h.shape[0]):
        totals=defaultdict(int)
        for offset in range(h.indptr[i],h.indptr[i+1]):
            j=int(h.indices[offset]);totals[int(colors[j])]+=int(h.data[offset])*int(phase[i])*int(phase[j])
        row={k:v for k,v in totals.items() if v};c=int(colors[i])
        if c in rows:
            if rows[c]!=row:raise ValueError('full-source equitable intertwiner failed')
        else:rows[c]=row
    q=s.ImmutableSparseMatrix(len(classes),len(classes),{(i,j):v for i,row in rows.items() for j,v in row.items()})
    g=s.diag(*(sizes[i] for i in classes));assert g*q==q.T*g
    origin=int(colors[keys.index((0,0))]);assert sizes[origin]==1
    blank=[i for i in classes if q[i,i]==0];exc=[i for i in classes if q[i,i]==1]
    assert len(blank)+len(exc)==len(classes)
    assert q.extract(blank,blank)==s.zeros(len(blank))
    assert q.extract(exc,exc)==s.eye(len(exc))
    a=q.extract(blank,exc)*q.extract(exc,blank)
    gb=g.extract(blank,blank);assert gb*a==a.T*gb
    return q,g,blank,exc,a,origin,depth


@lru_cache(maxsize=1)
def cyclic():
    q,g,blank,exc,a,origin,depth=quotient()
    v=s.zeros(len(blank),1);v[blank.index(origin)]=1
    columns=[v]
    # A proposed dimension is checked by exact independence and closure.
    for _ in range(18):columns.append(a*columns[-1])
    extended=s.Matrix.hstack(*columns);reduced,pivots=extended.rref()
    if pivots!=tuple(range(18)):raise ValueError('specified cyclic dimension fails exact rank check')
    c=s.Matrix.hstack(*columns[:18]);coefficients=s.Matrix(reduced[:18,18])
    assert c*coefficients==columns[18]
    companion=s.zeros(18)
    for j in range(17):companion[j+1,j]=1
    companion[:,17]=coefficients
    assert a*c==c*companion
    metric=c.T*g.extract(blank,blank)*c
    assert metric*companion==companion.T*metric
    x=s.Symbol('x');p=s.Poly(x**18-sum(coefficients[j]*x**j for j in range(18)),x)
    moments=[columns[k][blank.index(origin)] for k in range(18)]
    numerator=s.Poly(sum(sum(p.nth(j)*moments[j-k-1] for j in range(k+1,19))*x**k for k in range(18)),x)
    assert s.gcd(p,numerator).degree()==0 and numerator.nth(17)==1
    return p,numerator,c,companion,metric


def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
    values=[x*y for x in a for y in b];return min(values),max(values)
def divide(a,b):
    if b[0]<=0<=b[1]:raise ValueError('zero-crossing divisor interval')
    inverse=sorted((1/b[0],1/b[1]));return mul(a,inverse)
def evaluate(poly,interval):
    result=(s.Rational(0),s.Rational(0))
    for coefficient in poly.all_coeffs():result=add(mul(result,interval),(coefficient,coefficient))
    return result

def sqrt_interval(interval,bits=160):
    bounds=[];scale=1<<bits
    for value in interval:
        value=s.Rational(value);assert value>=0
        n=isqrt(int(value.p)*scale*scale//int(value.q))
        bounds.append((s.Rational(n,scale),s.Rational(n+1,scale)))
    return bounds[0][0],bounds[1][1]

def outward(interval,bits):
    scale=1<<bits
    return [str(s.Rational(s.floor(interval[0]*scale),scale)),str(s.Rational(s.ceiling(interval[1]*scale),scale))]


@lru_cache(maxsize=1)
def spectrum():
    q,g,blank,exc,a,origin,depth=quotient();p,numerator,c,companion,metric=cyclic()
    roots=[]
    for factor,multiplicity in s.factor_list(p)[1]:
        assert multiplicity==1
        for interval,m in s.polys.polytools.intervals(factor,eps=s.Rational(1,10**90)):
            assert m==1 and interval[0]>0;roots.append(tuple(map(s.Rational,interval)))
    roots.sort();assert len(roots)==18
    records=[]
    for root in roots:
        weight=divide(evaluate(numerator,root),evaluate(p.diff(),root))
        assert 0<weight[0]<=weight[1]<=1
        discriminant=sqrt_interval(add((1,1),mul((4,4),root)))
        inverse=divide((1,1),discriminant)
        for sign in (-1,1):
            energy=mul((s.Rational(1,2),)*2,add((1,1),mul((sign,sign),discriminant)))
            split=mul((s.Rational(1,2),)*2,add((1,1),mul((-sign,-sign),inverse)))
            probability=mul(weight,split)
            records.append({'energy_interval':outward(energy,40),'probability_interval':outward(probability,48),
                            'spectral_rank_in_cyclic_source':1})
    records.sort(key=lambda r:s.Rational(r['energy_interval'][0]))
    assert sum(s.Rational(r['probability_interval'][0]) for r in records)<=1
    assert sum(s.Rational(r['probability_interval'][1]) for r in records)>=1
    first=records[0]
    return {'full_connected_dimension':97878,'equitable_dimension':q.rows,
        'scalar_blank_dimension':len(blank),'scalar_excited_dimension':len(exc),
        'squared_cyclic_dimension':18,'full_physical_cyclic_dimension':36,
        'exact_full_graph_intertwiner':True,'cyclic_basis_is_polynomials_of_full_H_on_Gauss_singlet':True,
        'squared_minimal_polynomial':str(p.as_expr()),'squared_polynomial_factorization':str(s.factor(p.as_expr())),
        'vacuum_resolvent_numerator':str(numerator.as_expr()),
        'complete_energy_records':records,
        'outside_cyclic_preparation_sector_record':{'probability_interval':['0','0'],'retained':True},
        'ground_energy_interval':first['energy_interval'],
        'ground_filter_success_interval':first['probability_interval'],
        'ground_filter_failure_interval':[str(1-s.Rational(first['probability_interval'][1])),str(1-s.Rational(first['probability_interval'][0]))],
        'numerical_eigensolver_used':False,'global_other_sector_minimum_proved':False}
