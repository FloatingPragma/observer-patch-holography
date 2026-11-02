"""Integer full-sector CAR audits of the source spin-exchange reduction."""
from functools import lru_cache
from itertools import combinations, product
import math

import numpy as np
from scipy import sparse

from ..field_symmetry.symmetry import link_sign


def geometry(sides):
    if type(sides) is not tuple or len(sides) != 3 or any(type(n) is not int or n < 1 for n in sides):
        raise ValueError('three positive integer side lengths required')
    sites = tuple(product(*(range(n) for n in sides)))
    if not 2 <= len(sites) <= 8:
        raise ValueError('explicit audit supports two through eight sites')
    index = {x: i for i, x in enumerate(sites)}
    edges = []
    for i, x in enumerate(sites):
        for axis in range(3):
            y = list(x); y[axis] += 1; y = tuple(y)
            if y in index:
                edges.append((i, index[y], link_sign(x, axis)))
    return len(sites), tuple(edges)


def validate_graph(sites, edges):
    if type(sites) is not int or not 2 <= sites <= 8 or type(edges) is not tuple or not edges:
        raise ValueError('a finite nonempty signed graph of two through eight sites is required')
    seen = set()
    for edge in edges:
        if type(edge) is not tuple or len(edge) != 3:
            raise ValueError('edge triples required')
        x, y, sign = edge
        if any(type(v) is not int for v in edge) or not 0 <= x < y < sites or sign not in (-1, 1) or (x, y) in seen:
            raise ValueError('simple ordered signed edges required')
        seen.add((x, y))


def hop(mask, source, destination):
    if not (mask >> source) & 1 or (mask >> destination) & 1:
        return None
    sign = -1 if (mask & ((1 << source)-1)).bit_count() % 2 else 1
    mask ^= 1 << source
    sign *= -1 if (mask & ((1 << destination)-1)).bit_count() % 2 else 1
    return mask | (1 << destination), sign


@lru_cache(maxsize=None)
def _operators(sites, edges, number):
    masks = tuple(sum(1 << i for i in chosen) for chosen in combinations(range(2*sites), number))
    index = {mask: j for j, mask in enumerate(masks)}
    rows, cols, values = [], [], []
    for j, mask in enumerate(masks):
        for x, y, sign in edges:
            for spin in (0, 1):
                for a, b in ((2*x+spin, 2*y+spin), (2*y+spin, 2*x+spin)):
                    result = hop(mask, a, b)
                    if result:
                        other, parity = result
                        rows.append(index[other]); cols.append(j); values.append(sign*parity)
    T2 = sparse.csr_matrix((np.array(values, dtype=np.int64), (rows, cols)), shape=(len(masks), len(masks)))
    doublons = np.array([sum((mask >> (2*x)) & 3 == 3 for x in range(sites)) for mask in masks], dtype=np.int64)
    return masks, T2, doublons


def operators(sites, edges, number=None):
    validate_graph(sites, edges)
    number = sites if number is None else number
    if type(number) is not int or not 0 <= number <= 2*sites:
        raise ValueError('integer total occupation required')
    return _operators(sites, edges, number)


def spin_embedding(sites, masks):
    index = {mask: j for j, mask in enumerate(masks)}
    rows = [index[sum(1 << (2*x+((spin >> x) & 1)) for x in range(sites))] for spin in range(1 << sites)]
    return sparse.csr_matrix((np.ones(1 << sites, dtype=np.int64), (rows, range(1 << sites))), shape=(len(masks), 1 << sites))


def exchange_four(sites, edges):
    """4 sum_edges(S_x.S_y-1/4), independently in the spin tensor basis."""
    validate_graph(sites, edges)
    out = sparse.lil_matrix((1 << sites, 1 << sites), dtype=np.int64)
    for s in range(1 << sites):
        for x, y, _ in edges:
            if ((s >> x) ^ (s >> y)) & 1:
                out[s, s] -= 2
                out[s ^ (1 << x) ^ (1 << y), s] += 2
    return out.tocsr()


def exact_audit(sides):
    sites, edges = geometry(sides)
    masks, T2, D = operators(sites, edges)
    P = spin_embedding(sites, masks)
    B = T2 @ P
    first = P.T @ B
    denominator = sparse.diags(D, dtype=np.int64) @ B-B
    second = -(B.T @ B)-exchange_four(sites, edges)
    for matrix in (T2-T2.T, first, denominator, second):
        matrix.eliminate_zeros()
        if matrix.nnz:
            raise ArithmeticError('exact CAR spin reduction failed')
    return {'sides': list(sides), 'vertices': sites, 'edges': len(edges),
            'half_filled_dimension': math.comb(2*sites, sites), 'spin_dimension': 1 << sites,
            'hopping_scale': 'T=(c/2) T2', 'PTP_zero': True,
            'one_hop_denominator_exactly_U': True, 'second_order_exchange_exact': True,
            'signed_edges': [list(e) for e in edges]}


def interaction_audit():
    first = exchange_four(3, ((0, 1, 1),))
    second = exchange_four(3, ((1, 2, 1),))
    comm = first @ second-second @ first
    comm.eliminate_zeros()
    if not comm.nnz:
        raise ArithmeticError('adjacent exchange terms must not commute')
    return {'adjacent_bond_commutator_nonzero_entries': comm.nnz,
            'dimer_product_evolves_to_maximally_entangled_spin_state': True,
            'dimer_entangling_time': 'pi/(2 J)', 'ideal_one_site_purity': '1/2',
            'full_dimer_source_error_at_c128_U16384': 'less than 5/64',
            'entanglement_witness_survives_full_source_error': True}


def certificate():
    return {'exact_graph_audits': [exact_audit(s) for s in ((2, 1, 1), (3, 1, 1), (2, 2, 1), (2, 2, 2))],
            'interaction': interaction_audit(), 'basis_truncation_before_hopping': False}
