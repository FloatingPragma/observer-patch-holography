"""Exact full-sector parity checks, including the source's signed cube."""
import numpy as np
from scipy import sparse

from ..strong_coupling_refinement.algebra import (
    geometry, operators, spin_embedding, validate_graph,
)


def bipartition(sites, edges):
    validate_graph(sites, edges)
    adjacency = [[] for _ in range(sites)]
    for x, y, _ in edges:
        adjacency[x].append(y)
        adjacency[y].append(x)
    colors = [None]*sites
    for start in range(sites):
        if colors[start] is not None:
            continue
        colors[start] = 0
        pending = [start]
        while pending:
            x = pending.pop()
            for y in adjacency[x]:
                if colors[y] is None:
                    colors[y] = 1-colors[x]
                    pending.append(y)
                elif colors[y] == colors[x]:
                    raise ValueError('bipartite graph required')
    return tuple(colors)


def exact_audit(sides):
    sites, edges = geometry(sides)
    colors = bipartition(sites, edges)
    masks, T, doublons = operators(sites, edges)
    P = spin_embedding(sites, masks)
    selected = [x for x, color in enumerate(colors) if color == 0]
    parity = np.array([(-1)**sum(((mask >> (2*x)) & 3).bit_count()
                                for x in selected) for mask in masks], dtype=np.int64)
    S = sparse.diags(parity, dtype=np.int64)
    Q = sparse.eye(len(masks), dtype=np.int64, format='csr')-P@P.T
    B, C = Q@T@P, Q@T@Q
    identities = {
        'hopping_anticommutes_with_sublattice_parity': S@T+T@S,
        'spin_space_has_scalar_parity': S@P-(-1)**len(selected)*P,
        'one_hop_has_one_doublon': sparse.diags(doublons, dtype=np.int64)@B-B,
        'odd_virtual_return_vanishes': B.T@C@B,
    }
    for name, matrix in identities.items():
        matrix.eliminate_zeros()
        if matrix.nnz:
            raise ArithmeticError(name)
    return {
        'sides': list(sides), 'sites': sites, 'edges': len(edges),
        'half_filled_dimension': len(masks), 'spin_dimension': 2**sites,
        'sublattice': selected, 'identities': {name: True for name in identities},
        'virtual_charge_states_retained': True,
    }


def certificate():
    return [exact_audit(sides) for sides in ((2, 1, 1), (3, 1, 1), (2, 2, 1), (2, 2, 2))]
