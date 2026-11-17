"""Declared nested angular/radial complexes, with full native vertex blocks."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
import numpy as np
import scipy.sparse as sp
import sympy as sy

from ...spatial import Complex

BASE=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
BASE_FACES=((0,1,2),(0,3,1),(0,2,3),(1,3,2))


def canonical_face(face):
    face=tuple(face)
    forward=[face[i:]+face[:i] for i in range(len(face))]
    rev=face[::-1]
    backward=[rev[i:]+rev[:i] for i in range(len(face))]
    a,b=min(forward),min(backward)
    return (a,1) if a<b else (b,-1)


@dataclass
class Refinement:
    level:int
    n:int
    raw:tuple
    angular_faces:tuple
    labels:tuple
    graph:Complex
    cells:tuple
    B1:sp.csr_matrix
    B2:sp.csr_matrix
    B3:sp.csr_matrix


def shell(level):
    if not isinstance(level,int) or level<0:
        raise ValueError('nonnegative integer refinement level required')
    n=2**level
    raw=[];lookup={};triangles=[]
    def vertex(a,b,c,i,j):
        point=tuple(F((n-i-j)*a[k]+i*b[k]+j*c[k],n) for k in range(3))
        if point not in lookup:
            lookup[point]=len(raw);raw.append(point)
        return lookup[point]
    for a,b,c in BASE_FACES:
        a,b,c=BASE[a],BASE[b],BASE[c]
        for i in range(n):
            for j in range(n-i):
                triangles.append(tuple(vertex(a,b,c,x,y) for x,y in ((i,j),(i+1,j),(i,j+1))))
                if i+j<=n-2:
                    triangles.append(tuple(vertex(a,b,c,x,y) for x,y in ((i+1,j),(i+1,j+1),(i,j+1))))
    count=len(raw)
    labels=tuple((a,j) for j in range(n+1) for a in range(count))
    def at(a,j):return j*count+a
    units=[]
    for p in raw:
        p=tuple(sy.Rational(q.numerator,q.denominator) for q in p)
        norm=sy.sqrt(sum(x*x for x in p))
        units.append(tuple(x/norm for x in p))
    positions=tuple(tuple((1+sy.Rational(j,n))*x for x in units[a]) for a,j in labels)
    faces=[];flookup={};cells=[]
    def face_index(face):
        canonical,sign=canonical_face(face)
        if canonical not in flookup:
            flookup[canonical]=len(faces);faces.append(canonical)
        return flookup[canonical],sign
    for j in range(n):
        for tri in triangles:
            boundary=[(*face_index(tuple(at(a,j+1) for a in tri)),1),
                      (*face_index(tuple(at(a,j) for a in tri)),-1)]
            for u,v in zip(tri,tri[1:]+tri[:1]):
                boundary.append((*face_index((at(u,j),at(v,j),at(v,j+1),at(u,j+1))),1))
            cells.append(tuple((k,sign*factor) for k,sign,factor in boundary))
    edges=tuple(sorted({tuple(sorted((u,v))) for face in faces for u,v in zip(face,face[1:]+face[:1])}))
    graph=Complex(f'dyadic_tetrahedral_shell_level_{level}',positions,edges,tuple(faces))
    elookup={e:i for i,e in enumerate(edges)}
    r=[];c=[];x=[]
    for i,(u,v) in enumerate(edges):r.extend((u,v));c.extend((i,i));x.extend((1,-1))
    b1=sp.csr_matrix((x,(r,c)),shape=(len(positions),len(edges)),dtype=np.int64)
    r=[];c=[];x=[]
    for j,face in enumerate(faces):
        for u,v in zip(face,face[1:]+face[:1]):
            r.append(elookup[tuple(sorted((u,v)))]);c.append(j);x.append(1 if u<v else -1)
    b2=sp.csr_matrix((x,(r,c)),shape=(len(edges),len(faces)),dtype=np.int64)
    r=[];c=[];x=[]
    for j,boundary in enumerate(cells):
        for k,sign in boundary:r.append(k);c.append(j);x.append(sign)
    b3=sp.csr_matrix((x,(r,c)),shape=(len(faces),len(cells)),dtype=np.int64)
    return Refinement(level,n,tuple(raw),tuple(triangles),labels,graph,tuple(cells),b1,b2,b3)


def nested_maps(coarse,fine):
    if fine.level!=coarse.level+1:
        raise ValueError('maps are adjacent dyadic refinements')
    fl={p:i for i,p in enumerate(fine.raw)};fc=len(fine.raw)
    injection=tuple(2*j*fc+fl[coarse.raw[a]] for a,j in coarse.labels)
    fine_edges={e:i for i,e in enumerate(fine.graph.edges)}
    paths=[];r=[];c=[];x=[]
    for k,(u,v) in enumerate(coarse.graph.edges):
        a,j=coarse.labels[u];b,h=coarse.labels[v]
        if j==h:
            raw=tuple((p+q)/2 for p,q in zip(coarse.raw[a],coarse.raw[b]))
            middle=2*j*fc+fl[raw]
        else:
            assert a==b and h==j+1
            middle=(2*j+1)*fc+fl[coarse.raw[a]]
        path=(injection[u],middle,injection[v]);paths.append(path)
        for s,t in zip(path,path[1:]):
            r.append(fine_edges[tuple(sorted((s,t)))]);c.append(k);x.append(1 if s<t else -1)
    edge_map=sp.csr_matrix((x,(r,c)),shape=(len(fine.graph.edges),len(coarse.graph.edges)),dtype=np.int64)
    vertex_map=sp.csr_matrix((np.ones(len(injection),dtype=np.int64),(injection,range(len(injection)))),
                             shape=(fine.graph.nvertices,coarse.graph.nvertices))
    # A total block map for energy sums, distinct from an inverse physical state map.
    coarse_raw=coarse.raw
    nearest={p:min(range(len(coarse_raw)),key=lambda i:(sum((p[k]-coarse_raw[i][k])**2 for k in range(3)),i))
             for p in fine.raw}
    block=tuple((j//2)*len(coarse.raw)+nearest[fine.raw[a]] for a,j in fine.labels)
    return dict(vertex_injection=injection,edge_paths=tuple(paths),edge_chain=edge_map,
                vertex_chain=vertex_map,energy_block_map=block)


def audit(levels=range(4)):
    rows=[];members=[]
    for level in levels:
        m=shell(level);members.append(m)
        assert (m.B1@m.B2).nnz==0 and (m.B2@m.B3).nnz==0
        assert len(m.raw)==2*m.n*m.n+2 and len(m.angular_faces)==4*m.n*m.n
        assert m.graph.nvertices==(2*m.n*m.n+2)*(m.n+1)
        degrees=np.diff(m.B3.tocsr().indptr)
        assert np.all((degrees==1)|(degrees==2))
        totals=np.asarray(m.B3.sum(axis=1)).ravel()
        assert np.all(totals[degrees==2]==0)
        assert int(np.sum(degrees==1))==8*m.n*m.n
        # Reversing one face occurrence destroys closure of an actual cell.
        changed=m.B3.tolil();f,_=m.cells[0][0];changed[f,0]*=-1
        assert (m.B2@changed.tocsr()).nnz>0
        points=np.array(m.graph.positions,dtype=float)
        maxedge=max(float(np.linalg.norm(points[u]-points[v])) for u,v in m.graph.edges)
        bound=float((8*sy.sqrt(6)+1)/m.n)
        assert maxedge<=bound
        rows.append(dict(level=level,subdivisions=m.n,vertices=m.graph.nvertices,
                         native_modes=96*m.graph.nvertices,edges=len(m.graph.edges),
                         faces=len(m.graph.faces),cells=len(m.cells),boundary_faces=int(np.sum(degrees==1)),
                         measured_max_edge=maxedge,proved_cell_diameter_upper=bound,
                         exact_boundary_squared_zero=True,reversed_face_detected=True))
    maps=[]
    for coarse,fine in zip(members[:2],members[1:3]):
        m=nested_maps(coarse,fine)
        defect=fine.B1@m['edge_chain']-m['vertex_chain']@coarse.B1
        assert defect.nnz==0
        assert all(fine.graph.positions[j]==coarse.graph.positions[i] for i,j in enumerate(m['vertex_injection']))
        assert all(m['energy_block_map'][j]==i for i,j in enumerate(m['vertex_injection']))
        mode_map=[96*v+s for v in m['vertex_injection'] for s in range(96)]
        assert len(set(mode_map))==96*coarse.graph.nvertices
        bad=m['edge_chain'].tolil();row=bad[:,0].nonzero()[0][0];bad[row,0]=0
        assert (fine.B1@bad.tocsr()-m['vertex_chain']@coarse.B1).nnz>0
        maps.append(dict(coarse=coarse.level,fine=fine.level,retained_modes=len(mode_map),
                         incidence_intertwiner_exact=True,positions_preserved=True,
                         deleted_path_segment_detected=True))
    return dict(members=rows,maps=maps,scope='declared mesh and CAR/gauge incidence maps; no Hamiltonian or observer-readout convergence')
