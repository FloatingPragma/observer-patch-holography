#!/usr/bin/env python3
"""Compile the supplied native CAR factors to integer-link Gauss dressing.

The native snapshot retains exact binary coefficients, all 114 factors and
192 modes. The construction uses these values without refitting them.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,hashlib,json,sys,subprocess
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import sympy as sy
try:
    from .audit_core import annihilator,bfs_tree
except ImportError:
    from audit_core import annihilator,bfs_tree
PIN='b37ba4dae0e069dc030ceb0f130e52905fbb5465'

@lru_cache(None)
def normal(word):
    """Complete CAR normal ordering; creators first, each block index-sorted."""
    word=tuple(word)
    for k in range(len(word)-1):
        (i,ci),(j,cj)=word[k:k+2]
        if ci==cj and i==j:return ()
        if ((not ci) and cj) or (ci==cj and i>j):
            result=defaultdict(int)
            swapped=word[:k]+(word[k+1],word[k])+word[k+2:]
            for w,c in normal(swapped):result[w]-=c
            if not ci and cj and i==j:
                for w,c in normal(word[:k]+word[k+2:]):result[w]+=c
            return tuple((w,c) for w,c in result.items() if c)
    return ((word,1),)

def clean(poly):
    out={}
    for w,c in poly.items():
        c=sy.expand(c)
        if c:out[w]=c
    return out

def multiply(left,right):
    out=defaultdict(lambda:sy.Integer(0))
    for a,x in left.items():
        for b,y in right.items():
            for w,c in normal(a+b):out[w]+=x*y*c
    return clean(out)

def dagger(poly):
    out=defaultdict(lambda:sy.Integer(0))
    for w,c in poly.items():
        for v,s in normal(tuple((i,not creation) for i,creation in reversed(w))):out[v]+=sy.conjugate(c)*s
    return clean(out)

def add_scaled(out,poly,weight):
    for w,c in poly.items():out[w]+=weight*c

def exact_number(x):return sy.Rational(F(x))

def bilinear(entries):
    out=defaultdict(lambda:sy.Integer(0))
    for i,j,re,im in entries:out[((int(i),True),(int(j),False))]+=exact_number(re)+sy.I*exact_number(im)
    return clean(out)

def expand_source_terms(terms):
    out=defaultdict(lambda:sy.Integer(0))
    for term in terms:
        kind=term['kind']
        if kind=='positive_linear_fermion_square':
            L=defaultdict(lambda:sy.Integer(0))
            for key,creation in [('annihilation_coefficients',False),('creation_coefficients',True)]:
                for i,re,im in term[key]:L[((int(i),creation),)]+=exact_number(re)+sy.I*exact_number(im)
            p=multiply(dagger(clean(L)),clean(L));weight=exact_number(term['energy_electron_units_diagnostic'])
        else:
            weight=exact_number(term['weight'])
            if kind=='one_body':p=bilinear(term['matrix'])
            elif kind=='one_body_square':
                p=bilinear(term['matrix']);p=multiply(p,p)
            elif kind=='one_body_product':p=multiply(bilinear(term['left_matrix']),bilinear(term['right_matrix']))
            else:raise ValueError('Unsupported term kind '+kind)
        add_scaled(out,p,weight)
    return clean(out)

def charge_changes(word,mode_charges,mode_vertices,n):
    s=[0]*n
    for i,creation in word:s[mode_vertices[i]]+=(1 if creation else -1)*mode_charges[i]
    return s

def tree_flow(s,edges,root=0):
    n=len(s);parent,pe,depth,order=bfs_tree(n,edges,root)
    if sum(s):raise ValueError('Nonneutral monomial: explicit charged compensator required')
    r=list(s);flow=[0]*len(edges)
    for v in reversed(order[1:]):
        a,b=edges[pe[v]]
        flow[pe[v]]=r[v] if a==v else -r[v]
        r[parent[v]]+=r[v];r[v]=0
    assert not any(r)
    return flow

def compile_poly(poly,mode_charges,mode_vertices,n,edges):
    poly=multiply({():sy.Integer(1)},poly)
    if clean(poly)!=dagger(poly):raise ValueError('Input polynomial is not Hermitian')
    if len(mode_charges)!=len(mode_vertices):raise ValueError('Bad mode metadata')
    result=[]
    for word,coefficient in sorted(poly.items()):
        s=charge_changes(word,mode_charges,mode_vertices,n)
        flow=tree_flow(s,edges)
        # B flow = charge change, checked independently from the tree algorithm.
        check=[0]*n
        for value,(a,b) in zip(flow,edges):check[a]+=value;check[b]-=value
        assert check==s
        result.append({'word':[[i,c] for i,c in word], 'coefficient_exact_as_supplied':str(coefficient),
                       'charge_change':s,'link_shift':flow})
    return result

def tests():
    # Verify all length <=4 words on two modes against the full four-state CAR representation.
    import itertools
    ann=[annihilator(2,j) for j in range(2)];letters=[(i,c) for i in range(2) for c in (False,True)]
    count=0
    def matrix(w):
        m=sy.eye(4)
        for i,c in w:m=m*(ann[i].T if c else ann[i])
        return m
    for k in range(5):
        for word in itertools.product(letters,repeat=k):
            rebuilt=sum((s*matrix(w) for w,s in normal(word)),sy.zeros(4))
            assert matrix(word)==rebuilt;count+=1
    # Neutral cross-cell weak-style exchange and pair/hop terms.
    charges=[-3,0,-3,0];vertices=[0,0,1,1]
    O={((1,True),(0,False),(2,True),(3,False)):sy.Rational(1,7)}
    poly=defaultdict(lambda:sy.Integer(0));add_scaled(poly,O,1);add_scaled(poly,dagger(O),1)
    hop={((2,True),(0,False)):sy.Rational(-1,8)}
    add_scaled(poly,hop,1);add_scaled(poly,dagger(hop),1)
    poly=clean(poly)
    result=compile_poly(poly,charges,vertices,2,[(0,1)])
    assert all(abs(x['link_shift'][0])==3 for x in result)
    # A global-charge-violating monomial must be rejected, not silently dressed.
    bad={((0,True),):sy.Integer(1),((0,False),):sy.Integer(1)}
    try:compile_poly(bad,charges,vertices,2,[(0,1)])
    except ValueError:rejected=True
    else:rejected=False
    assert rejected
    return {'normal_order_matrix_words_checked':count,'gauged_compound_and_hop_terms':result,
            'nonneutral_input_rejected':True,'full_native_import_executed_in_this_test':False}

def native_export(repo=None):
    """Compile the complete specified snapshot without a checkout dependency."""
    try:
        from .snapshot import snapshot
    except ImportError:
        from snapshot import snapshot
    src=snapshot()
    ch=src['mode_charge_native_thirds']; vtx=src['mode_vertex']
    poly=expand_source_terms(src['terms'])
    terms=compile_poly(poly,ch,vtx,2,[(0,1)])
    # This serialization is the canonical factor-polynomial commitment.
    return {'upstream_commit':PIN,'source_terms':len(src['terms']), 'mode_count':192,
      'expanded_monomials':len(terms),'terms':terms,'mode_charge_native_thirds':ch,'mode_vertex':vtx,
      'scope':'Exact gauging of supplied dyadic snapshots. This adds a field architecture, changes interacting spectra, and does not assert physical masses or photon identification.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--native',action='store_true');p.add_argument('--output',type=Path);a=p.parse_args()
    out=native_export() if a.native else tests()
    txt=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if a.output:a.output.write_text(txt)
    print(txt,end='')
if __name__=='__main__':main()
