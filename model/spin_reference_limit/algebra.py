"""Complete finite occupation controls for the charge/spin dictionary."""
import numpy as np

from ..strong_coupling_refinement.algebra import hop, validate_graph


def carrier(sites, edges):
    validate_graph(sites, edges)
    if sites > 4:
        raise ValueError('explicit complete-carrier audit supports at most four sites')
    dimension = 4**sites
    T2 = np.zeros((dimension, dimension), dtype=np.int64)
    defect = np.zeros((sites, dimension), dtype=np.int64)
    permutation = np.zeros((dimension, dimension), dtype=np.int64)
    for mask in range(dimension):
        for x in range(sites):
            defect[x, mask] = int(((mask >> (2*x)) & 3).bit_count() != 1)
        for x, y, sign in edges:
            for spin in (0, 1):
                for a, b in ((2*x+spin, 2*y+spin), (2*y+spin, 2*x+spin)):
                    moved = hop(mask, a, b)
                    if moved is not None:
                        other, parity = moved
                        T2[other, mask] += sign*parity
    for charge in range(2**sites):
        for spin in range(2**sites):
            mask = 0
            for x in range(sites):
                q, s = (charge >> x) & 1, (spin >> x) & 1
                local = (1 << s) if not q else (3 if s else 0)
                mask |= local << (2*x)
            permutation[mask, charge*2**sites+spin] = 1
    return T2, defect, permutation


def audit(sites, edges):
    T2, Q, dictionary = carrier(sites, edges)
    N = Q.sum(axis=0)
    P = dictionary[:, :2**sites]
    D0 = np.where(N[:, None] == N[None, :], T2, 0)
    # Any spin matrix is extended as the identity on all charge factors.
    flip = np.eye(2**sites, dtype=np.int64)[:, np.arange(2**sites) ^ 1]
    embedded = dictionary@np.kron(np.eye(2**sites, dtype=np.int64), flip)@dictionary.T
    identities = {
        'dictionary_unitary': dictionary.T@dictionary-np.eye(4**sites, dtype=np.int64),
        'onsite_penalty_annihilates_spin_space': N[:, None]*P,
        'hopping_compression_zero': P.T@T2@P,
        'neutral_hopping_annihilates_spin_space': D0@P,
        'embedded_spin_commutes_with_penalty': N[:, None]*embedded-embedded*N[None, :],
        'embedded_spin_restricts_to_original_spin': P.T@embedded@P-flip,
        'embedding_is_unital_product_preserving': embedded@embedded-np.eye(4**sites, dtype=np.int64),
    }
    if any(np.any(value) for value in identities.values()):
        raise ArithmeticError('complete source charge/spin identity failed')
    if set(N) != set(range(sites+1)):
        raise ArithmeticError('every defect sector must be present')
    return {'sites': sites, 'full_Fock_dimension': 4**sites,
            'spin_dimension': 2**sites, 'defect_sectors': list(range(sites+1)),
            'identities': {key: True for key in identities},
            'source_projected_before_hopping': False}


def certificate():
    return [audit(2, ((0, 1, 1),)), audit(3, ((0, 1, -1), (1, 2, 1)))]
