"""Exact complete-Fock checks of charge conservation and the repair join."""
from functools import lru_cache

import numpy as np

from ..defect_field_bridge import full_domain


def dimensions(sites):
    if type(sites) is not int or not 2 <= sites <= 4:
        raise ValueError('the explicit full-Fock oracle supports two to four sites')
    return 1 << (2*sites)


def defects(sites):
    dim = dimensions(sites)
    return np.array([sum((((mask >> (2*x)) & 3).bit_count()-1)**2
                         for x in range(sites)) for mask in range(dim)], dtype=np.int64)


def protected(sites, support=None):
    support = range(sites) if support is None else support
    return np.array([all(((mask >> (2*x)) & 3).bit_count() == 1 for x in support)
                     for mask in range(dimensions(sites))], dtype=bool)


@lru_cache(maxsize=None)
def branch_matrices(sites, law, seam):
    """Integer matrices equal sqrt(2) times the complete native Kraus maps."""
    dim = dimensions(sites)
    result = []
    rows = full_domain.finite_branches(sites=sites, flavors=2, law=law, seam=seam)
    for coin in (0, 1):
        for raw in range(full_domain.RAW_SLOTS):
            matrix = np.zeros((dim, dim), dtype=np.int64)
            for origin, output, sign in rows[coin, raw]:
                matrix[output, origin] = sign
            matrix.setflags(write=False)
            result.append(matrix)
    return tuple(result)


def twice_unread_dual(observable, sites, law, seam):
    return sum((K.T@observable@K for K in branch_matrices(sites, law, seam)
                if np.any(K)), np.zeros_like(observable))


def spin_embedding(sites, spin):
    """Unital charge tensor spin dictionary; no occupation state is removed."""
    dim = dimensions(sites)
    spin = np.asarray(spin)
    if spin.shape != (1 << sites, 1 << sites):
        raise ValueError('a complete spin tensor matrix is required')
    # Per site: up=(q=0,s=0), down=(0,1), empty=(1,0), double=(1,1).
    dictionary = [(1, 0), (0, 0), (0, 1), (1, 1)]
    labels = []
    for mask in range(dim):
        charge = sum(dictionary[(mask >> (2*x)) & 3][0] << x for x in range(sites))
        label = sum(dictionary[(mask >> (2*x)) & 3][1] << x for x in range(sites))
        labels.append((charge, label))
    return np.array([[spin[a, b] if qa == qb else 0 for qb, b in labels]
                     for qa, a in labels], dtype=spin.dtype)


def annihilator(sites, mode):
    """Independent Jordan-Wigner matrix on every occupation mask."""
    dim = dimensions(sites)
    result = np.zeros((dim, dim), dtype=np.int64)
    for mask in range(dim):
        if (mask >> mode) & 1:
            result[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return result


def hopping_twice(sites, edges):
    """The unit-c source hopping is T2/2, constructed from CAR matrices."""
    dim = dimensions(sites)
    out = np.zeros((dim, dim), dtype=np.int64)
    modes = [annihilator(sites, j) for j in range(2*sites)]
    for x, y, sign in edges:
        for spin in (0, 1):
            a, b = modes[2*x+spin], modes[2*y+spin]
            out += sign*(a.T@b+b.T@a)
    return out


def source_audit(sites, law, seam):
    dim = dimensions(sites)
    N, P = defects(sites), protected(sites, seam)
    identity = np.eye(dim, dtype=np.int64)
    effects = np.zeros(dim, dtype=np.int64)
    outputs = np.zeros(dim, dtype=np.int64)
    for index, K in enumerate(branch_matrices(sites, law, seam)):
        if np.any((N[:, None]-N[None, :])*K):
            raise ArithmeticError('native repair changed the defect number')
        target = identity[:, P] if index % 16 == 0 else np.zeros((dim, sum(P)), dtype=np.int64)
        if not np.array_equal(K[:, P], target):
            raise ArithmeticError('complete protected instrument is not the ideal append')
        # Each branch is a partial injection, so both products are diagonal.
        if np.any(np.count_nonzero(K, axis=0) > 1) or np.any(np.count_nonzero(K, axis=1) > 1):
            raise ArithmeticError('native branch is not a partial injection')
        effects += np.sum(K*K, axis=0)
        outputs += np.sum(K*K, axis=1)
    if not np.all(effects == 2) or not np.all(outputs == 2):
        raise ArithmeticError('complete repair is not trace preserving and unital')
    return {'sites': sites, 'law': law, 'seam': list(seam),
            'full_Fock_dimension': dim, 'coin_raw_slots': 32,
            'every_branch_commutes_with_defect_number': True,
            'complete_protected_append_identity': True,
            'unread_channel_unital_and_trace_preserving': True}


def hopping_audit(law):
    sites, edges = 3, ((0, 1, -1), (1, 2, 1))
    N, P = defects(sites), protected(sites)
    T2 = hopping_twice(sites, edges)
    D02 = np.where(N[:, None] == N[None, :], T2, 0)
    if np.any(T2[np.ix_(P, P)]) or np.any(D02[:, P]) or np.any(D02[P, :]):
        raise ArithmeticError('charge-diagonal hopping must annihilate protected inputs')
    # Exact zero protected sandwiches for every overlapping repair/hopping pair.
    increments = []
    for seam in ((0, 1), (1, 2), (0, 2)):
        for edge in edges:
            term = hopping_twice(sites, (edge,))
            increment = twice_unread_dual(term, sites, law, seam)-2*term
            support = sorted(set(seam+edge[:2]))
            local_P = protected(sites, support)
            if np.any(increment[np.ix_(local_P, local_P)]):
                raise ArithmeticError('repair energy increment has a protected component')
            increments.append(int(np.count_nonzero(increment)))
    if not any(increments):
        raise ArithmeticError('control must retain a nonzero repair energy increment')
    return {'law': law, 'all_64_Fock_states_retained': True,
            'charge_diagonal_hopping_annihilates_protected_space': True,
            'energy_increment_zero_protected_sandwich': True,
            'nonzero_increment_entries': increments}


def record_correction_twice(spin_out, record_out, spin_empty, law='C'):
    """Twice C=I-E-I0+id on a vacuum/one-event archive cylinder.

    The one-event readout is spin_out tensor record_out; spin_empty acts
    on the empty-archive block. The ideal coin is coherent |+>, raw zero.
    """
    A, empty = spin_embedding(2, spin_out), spin_embedding(2, spin_empty)
    record_out = np.asarray(record_out)
    if record_out.shape != (32, 32):
        raise ValueError('complete two-coin by sixteen-raw record required')
    matrices = branch_matrices(2, law, (0, 1))
    result = np.zeros_like(A, dtype=np.result_type(A, empty, record_out))
    nonzero = [j for j, K in enumerate(matrices) if np.any(K)]
    for j in nonzero:
        for k in nonzero:
            result += record_out[j, k]*(matrices[j].T@A@matrices[k])
    result -= twice_unread_dual(empty, 2, law, (0, 1))
    result -= (record_out[0, 0]+record_out[0, 16]+record_out[16, 0]+record_out[16, 16])*A
    return result+2*empty


def certificate():
    record_checks = []
    for law in ('C', 'E'):
        spin = np.arange(16, dtype=np.int64).reshape((4, 4)) % 5-2
        tape = np.arange(1024, dtype=np.int64).reshape((32, 32)) % 7-3
        correction = record_correction_twice(spin, tape, spin.T, law)
        P, N = protected(2), defects(2)
        if (np.any(correction[:, P]) or np.any(correction[P, :])
                or np.any((N[:, None]-N[None, :])*correction)):
            raise ArithmeticError('complete coherent archive correction is not charge supported')
        record_checks.append({'law': law, 'all_32_coin_raw_slots_retained': True,
                              'coherent_record_correction_charge_supported': True})
    return {
        'complete_native_branch_audits': [source_audit(sites, law, (0, sites-1))
                                         for sites in (2, 3, 4) for law in ('C', 'E')],
        'repair_energy_controls': [hopping_audit(law) for law in ('C', 'E')],
        'complete_archive_controls': record_checks,
        'proof_uses_full_source_and_complete_original_instrument': True,
    }
