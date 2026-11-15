"""Exact spatial class of the literal carrier's committed current cochain.

Standard library only. This constructs the integer lattice, rather than
inferring integral topology from a numerical rank or a charge label.
"""
from collections import deque
from fractions import Fraction
from itertools import combinations
import cmath
import math


def determinant(matrix):
    a=[list(map(Fraction,row)) for row in matrix]
    det=Fraction(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return Fraction(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];det=-det
        d=a[k][k];det*=d
        a[k]=[x/d for x in a[k]]
        for i in range(k+1,len(a)):
            c=a[i][k]
            a[i]=[x-c*y for x,y in zip(a[i],a[k])]
    return det


def construct(edges):
    edges=tuple(tuple(e) for e in edges);lookup={e:i for i,e in enumerate(edges)}
    raw=[f for f in combinations(range(12),3) if all(e in lookup for e in combinations(f,2))]
    def boundary(face):
        return {lookup[tuple(sorted((a,b)))]:1 if a<b else -1
                for a,b in zip(face,face[1:]+face[:1])}
    oriented={0:raw[0]}
    while len(oriented)<len(raw):
        grew=False
        for j,face in enumerate(raw):
            if j in oriented:continue
            for done in tuple(oriented.values()):
                common=set(boundary(face))&set(boundary(done))
                if len(common)==1:
                    e=next(iter(common))
                    oriented[j]=face if boundary(face)[e]==-boundary(done)[e] else (face[0],face[2],face[1])
                    grew=True;break
        if not grew:raise ValueError('disconnected triangle surface')
    A=[[0]*20 for _ in edges];B=[[0]*len(edges) for _ in range(12)]
    for j in range(20):
        for e,s in boundary(oriented[j]).items():A[e][j]=s
    for e,(a,b) in enumerate(edges):B[a][e]=1;B[b][e]=-1
    assert len(edges)==30 and len(raw)==20
    assert all(sorted(x for x in row if x)==[-1,1] for row in A)
    assert all(sum(B[v][e]*A[e][f] for e in range(30))==0 for v in range(12) for f in range(20))
    tree=[];seen={0};todo=deque([0])
    while todo:
        v=todo.popleft()
        for e,(a,b) in enumerate(edges):
            if v not in (a,b):continue
            w=b if v==a else a
            if w not in seen:seen.add(w);todo.append(w);tree.append(e)
    chords=[e for e in range(30) if e not in tree]
    tree_det=determinant([[B[v][e] for e in tree] for v in range(1,12)])
    cycle_det=determinant([[A[e][f] for f in range(19)] for e in chords])
    assert abs(tree_det)==abs(cycle_det)==1
    dual=[]
    for e in range(30):dual.append((A[e].index(-1),A[e].index(1)))
    loops=[]
    for v in range(12):
        step={}
        for e,(a,b) in enumerate(dual):
            if not B[v][e]:continue
            s=B[v][e];u,w=(a,b) if s==1 else (b,a)
            assert u not in step;step[u]=(w,e,s)
        start=min(step);u=start;loop=[]
        for _ in range(5):u,e,s=step[u];loop.append((e,s))
        assert u==start and len(set(e for e,s in loop))==5
        loops.append(loop)
    return {'seams':edges,'faces':[oriented[j] for j in range(20)],'B':B,'A':A,
            'dual_edges':dual,'dual_port_boundaries':loops,'tree_edges':tree,'chord_edges':chords,
            'tree_incidence_minor_determinant':int(tree_det),
            'nineteen_face_cycle_minor_determinant':int(cycle_det),
            'neutral_charge_rank':11,'integer_cycle_rank':19}


def divergence(edges,current):
    if len(current)!=len(edges):raise ValueError('complete edge cochain required')
    q=[0]*12
    for (a,b),j in zip(edges,current):q[a]+=j;q[b]-=j
    return tuple(q)


def discrete_curve_degree(geometry,current,region):
    """Full lifted C6 paths on the actual dual boundary; no phase unwrapping.

    The invariant is computed from oriented spatial paths, including every
    edge traversal. Fractional edges are not shared-vertex circle maps.
    """
    if any(type(j)is not int for j in current):raise ValueError('dual edge does not close at its shared vertices')
    transitions=0
    for v in region:
        for e,s in geometry['dual_port_boundaries'][v]:
            transitions+=6*s*current[e]
    assert transitions%6==0
    return transitions//6


def raw_pointer_degree(pointer,cycle):
    """Decode only actual measurement currents on the six spatial seams."""
    oriented={tuple(sorted((a,b))):(1 if a<b else -1)
              for a,b in zip(cycle,cycle[1:]+cycle[:1])}
    steps=0
    for record in pointer['records']:
        if record['stage']!='measure':continue
        edge=tuple(record['seam']);j=record['current']
        if edge not in oriented or type(j)is not int or abs(j)>1:raise ValueError('not a native cycle record')
        steps+=oriented[edge]*j
    if pointer['success'] is not True:return None
    if steps%6:raise ValueError('pointer path is not closed')
    return steps//6


def phase_degree(geometry,current,region,deformed=False):
    """Fixed-resolution diagnostic of actual glued spatial circle curves.

    The stated |F_e|<=3 domain implies phase increments below pi/2 with
    64 subdivisions. The exact integral-lattice theorem has no sampling cap.
    """
    if any(type(j)is not int or abs(j)>3 for j in current):raise ValueError('outside fixed phase-resolution domain')
    total=0.
    for v in region:
        contour=[]
        for e,s in geometry['dual_port_boundaries'][v]:
            a,b=geometry['dual_edges'][e]
            ga=.31*math.sin(a+1) if deformed else 0.
            gb=.31*math.sin(b+1) if deformed else 0.
            amp=.23*math.cos(a+b+1) if deformed else 0.
            values=[cmath.exp(1j*(ga+(gb-ga)*k/64+2*math.pi*current[e]*k/64+amp*math.sin(2*math.pi*k/64))) for k in range(65)]
            if s==-1:values.reverse()
            if contour:assert abs(contour[-1]-values[0])<1e-12
            contour.extend(values if not contour else values[1:])
        assert abs(contour[0]-contour[-1])<1e-12
        increments=[cmath.phase(b*a.conjugate()) for a,b in zip(contour,contour[1:])]
        assert max(map(abs,increments))<math.pi/2
        total+=math.fsum(increments)/(2*math.pi)
    return total
