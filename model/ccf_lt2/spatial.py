"""Explicit finite oriented carriers and compact CCF field operators.

The radial shell is a declared S2 x interval carrier. It is not a proof of
the inherited observer spatial metric. Finite field matrices retain all
flux states including walls and Gauss-mismatch outcomes.
"""
from dataclasses import dataclass
from itertools import product
import numpy as np
import scipy.sparse as sp
import sympy as sy


@dataclass(frozen=True)
class Complex:
    name: str
    positions: tuple
    edges: tuple
    faces: tuple

    @property
    def nvertices(self):
        return len(self.positions)

    def incidence(self):
        if len({tuple(sorted(e)) for e in self.edges}) != len(self.edges):
            raise ValueError('duplicate undirected edge')
        if any(u == v or not (0 <= u < self.nvertices and 0 <= v < self.nvertices)
               for u, v in self.edges):
            raise ValueError('invalid edge endpoints')
        B = np.zeros((self.nvertices, len(self.edges)), dtype=np.int64)
        C = np.zeros((len(self.faces), len(self.edges)), dtype=np.int64)
        lookup = {e: j for j, e in enumerate(self.edges)}
        for j, (u, v) in enumerate(self.edges):
            B[u, j], B[v, j] = 1, -1
        for j, face in enumerate(self.faces):
            if len(face) < 3 or len(set(face)) != len(face):
                raise ValueError('face boundary must contain at least three distinct vertices without repetition')
            for u, v in zip(face, face[1:] + face[:1]):
                if (u, v) in lookup:
                    C[j, lookup[(u, v)]] += 1
                elif (v, u) in lookup:
                    C[j, lookup[(v, u)]] -= 1
                else:
                    raise ValueError('face boundary edge is absent')
        if np.any(B @ C.T):
            raise ValueError('face boundary is not closed')
        return B, C

    def edge_geometry(self):
        directions, spacings = [], []
        for u, v in self.edges:
            delta = [sy.sympify(b)-sy.sympify(a) for a, b in zip(self.positions[u], self.positions[v])]
            length = sy.sqrt(sum(x*x for x in delta))
            if length == 0:
                raise ValueError('zero length edge')
            directions.append(tuple(sy.simplify(x/length) for x in delta))
            spacings.append(sy.simplify(length))
        return directions, spacings


def edge():
    return Complex('one_edge_chart', ((0, 0, 0), (1, 0, 0)), ((0, 1),), ())


def triangle():
    return Complex('triangle_chart', ((0, 0, 0), (1, 0, 0), (0, 1, 0)),
                   ((0, 1), (0, 2), (1, 2)), ((0, 1, 2),))


def radial_shell(layers=2):
    """Tetrahedral S2 triangulation times a finite radial interval.

    Each vertex holds a complete native matter block (not a changed internal
    twelve-port cell). No claim of rotational/refinement convergence is made.
    """
    if layers < 2:
        raise ValueError('radial product requires at least two layers')
    base = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    sphere_faces = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
    sphere_edges = tuple((u, v) for u in range(4) for v in range(u+1, 4))
    radii = tuple(1 + sy.Rational(j, layers-1) for j in range(layers))
    positions = tuple(tuple(r*x/sy.sqrt(3) for x in v) for r in radii for v in base)
    edges = [(4*j+u, 4*j+v) for j in range(layers) for u, v in sphere_edges]
    edges += [(4*j+u, 4*(j+1)+u) for j in range(layers-1) for u in range(4)]
    faces = [tuple(4*j+v for v in f) for j in range(layers) for f in sphere_faces]
    faces += [(4*j+u, 4*j+v, 4*(j+1)+v, 4*(j+1)+u)
              for j in range(layers-1) for u, v in sphere_edges]
    graph = Complex('tetrahedral_sphere_radial_product', positions, tuple(sorted(edges)), tuple(faces))
    graph.incidence()
    return graph


def shift_matrix(cap, shift=1):
    size = 2*cap+1
    if abs(shift) >= size:
        return sp.csr_matrix((size, size), dtype=complex)
    return sp.diags(np.ones(size-abs(shift)), -shift, shape=(size, size), dtype=complex, format='csr')


def kron_all(factors):
    out = sp.csr_matrix([[1.]], dtype=complex)
    for factor in factors:
        out = sp.kron(out, factor, format='csr')
    return out


def field_operators(graph, cap, electric=1, magnetic=1):
    """Complete bounded-flux Hamiltonian, no physical-sector projection."""
    B, C = graph.incidence()
    m = len(graph.edges)
    size = 2*cap+1
    dimension = size**m
    if dimension > 100000:
        raise ValueError('explicit field witness exceeds 100000 states')
    eye = sp.eye(size, dtype=complex, format='csr')
    identity = sp.eye(dimension, dtype=complex, format='csr')
    E = [kron_all([sp.diags(np.arange(-cap, cap+1), format='csr') if k == e else eye
                   for k in range(m)]) for e in range(m)]
    loops = [kron_all([shift_matrix(cap, int(value)) for value in row]) for row in C]
    H_e = sum((float(electric)/2*(x@x) for x in E), sp.csr_matrix(identity.shape, dtype=complex))
    H_b = sum((float(magnetic)*(identity-(W+W.conj().T)/2) for W in loops),
              sp.csr_matrix(identity.shape, dtype=complex))
    gauss = [sum((int(B[v, e])*E[e] for e in range(m)), sp.csr_matrix(identity.shape, dtype=complex))
             for v in range(graph.nvertices)]
    return dict(H=H_e+H_b, electric=H_e, magnetic=H_b, E=E, loops=loops,
                gauss=gauss, dimension=dimension, cap=cap)
