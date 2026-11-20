"""Integer certificates for charged rest multiplicity and the full nonconverting input space."""
from fractions import Fraction
from functools import reduce
from math import lcm

import numpy as np
from scipy.sparse.csgraph import shortest_path
import sympy as sy

from ...contrast_memory_spherical_v1 import carrier
from ..carrier import require
from .source import one_body


def mod_rank(matrix, p=101):
    a = np.array(matrix,dtype=np.int64)%p
    k = 0
    for j in range(a.shape[1]):
        rows = np.flatnonzero(a[k:,j])
        if not len(rows):
            continue
        i = k+rows[0]
        a[[k,i]] = a[[i,k]]
        pivot = a[k,j:]*pow(int(a[k,j]),-1,p)%p
        below = np.flatnonzero(a[k+1:,j])+k+1
        a[below,j:] = (a[below,j:]-a[below,j,None]*pivot)%p
        k += 1
        if k==len(a):
            break
    return k


def rest_certificate(Q,level):
    # F0 + sqrt(5)*F1 has range the 5-sqrt(5) eigenspace. Reduce sqrt(5) to 45 mod 101.
    I = np.eye(12,dtype=int)
    other = (Q-4*I) @ (Q-10*I)
    F0,F1 = other@(Q-5*I),-other
    F = (F0+45*F1)%101
    car = carrier.federation(level)
    inside = [set() for _ in range(car.cells)]
    for a,b in car.seams:
        if a//12==b//12:
            inside[a//12].add((a%12,b%12))
    require(all(edges==set(carrier.CELL_SEAMS) for edges in inside)
            and len({i for edge in car.external for i in edge})==2*len(car.external),
            'carrier lost a complete cell or the external matching')
    incident = [[] for _ in range(car.cells)]
    for e,(i,j) in enumerate(car.external):
        incident[i//12].append((e,i%12))
        incident[j//12].append((e,j%12))
    ties = {i:set() for i in range(len(car.external))}
    killed = set()
    constrained = 0
    for triple in incident:
        require(len(triple)==3,'external tripod differs')
        ids,ports = zip(*triple)
        rank = mod_rank(F[list(ports)])
        if rank==3:
            killed.update(ids)
        else:
            pairs = [(u,v) for u in range(3) for v in range(u) if
                     all(np.array_equal(M[ports[u]],-M[ports[v]]) for M in (F0,F1))]
            require(rank==2 and len(pairs)==1,'unclassified local row dependence')
            u,v = pairs[0]
            ties[ids[u]].add(ids[v])
            ties[ids[v]].add(ids[u])
            killed.add(ids[3-u-v])
            constrained += 1
    visited = set()
    free = 0
    for e in ties:
        if e in visited:
            continue
        component = {e}
        stack = [e]
        while stack:
            nxt = ties[stack.pop()]-component
            component.update(nxt)
            stack.extend(nxt)
        visited.update(component)
        free += int(component.isdisjoint(killed))
    return dict(cells=car.cells,external_seams=len(car.external),dependent_tripods=constrained,
                independent_row_relations=free,rest_multiplicity=3*car.cells-len(car.external)+free)


def bright_certificate(Q):
    n = len(Q)
    H = np.kron(Q,np.eye(n,dtype=int))+np.kron(np.eye(n,dtype=int),Q)
    D = np.eye(n*n,dtype=int)[:,np.arange(n)*(n+1)]
    Hc = H-10*np.eye(n*n,dtype=int)
    blocks = [D]
    for _ in range(9):
        blocks.append(Hc@blocks[-1])
    minor_rank = mod_rank(np.column_stack(blocks),1009)
    # Floating point proposes a rational projector. Only the integer identities below certify it.
    values,V = np.linalg.eigh(H)
    P = np.zeros(H.shape)
    for group in np.split(np.arange(len(H)),np.flatnonzero(np.diff(values)>1e-8)+1):
        E = V[:,group]
        U,s,_ = np.linalg.svd(E.T@D,full_matrices=False)
        visible = E @ U[:,s>1e-9]
        P += visible@visible.T
    rationals = {x:Fraction(float(x)).limit_denominator(10**6) for x in np.unique(P.round(12))}
    den = reduce(lcm,(f.denominator for f in rationals.values()),1)
    B = np.array([[int(den*rationals[x]) for x in row] for row in P.round(12)],dtype=object)
    require(np.array_equal(B,B.T) and np.array_equal(B@B,den*B) and
            np.array_equal(H@B,B@H) and np.array_equal(B@D,den*D) and
            B.trace()==minor_rank*den,'uncertified converting subspace')
    distances = shortest_path(carrier.cell().adjacency(),unweighted=True).astype(int)
    bounds = {}
    for d in np.unique(distances):
        entries = set(B.diagonal()[distances.ravel()==d])
        require(len(entries)==1,'distance class not homogeneous')
        bounds[str(d)] = str(Fraction(int(entries.pop()),den))
    return dict(charged_dimension=n*n,converting_span=minor_rank,exact_dark_dimension=n*n-minor_rank,
                projector_denominator=den,neutral_probability_upper_bound_by_port_distance=bounds)


def certificates():
    Q = one_body(carrier.cell())[0].toarray().astype(int)
    x = sy.symbols('x')
    require(sy.Matrix(Q).charpoly(x).as_poly()==sy.Poly((x-10)*(x-4)**5*((x-5)**2-5)**3,x),
            'wrong integer cell polynomial')
    return dict(rest=[rest_certificate(Q,i) for i in range(3)],one_cell=bright_certificate(Q))
