"""Constants of the binding theorem on the federations and on one cell.

In the sector of two loads (doubled loads included) let `P` be the span of the pairs on a seam and `Q` its
complement. On every federation of complete cells with a matching of external seams

    P H P <= (10 kappa + 6) P,    Q H Q >= 12 kappa Q,    ||Q H P|| <= sqrt(10),

and for `kappa > 3` a state that starts in `P` has `Pr(no pair on a seam) <= 40/(2 kappa - 6)^2` at every value
of the evolution parameter. The proof is in PROOF.md. This module measures the two carrier constants that
enter: `amax`, the largest level of the adjacent block of the exchange (bound 6), and `beta`, the largest
singular value of the coupling block `Q K P` (bound sqrt 10); and the exact single-cell constants.
"""
from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from model.executable_model import carrier as cr
from .sectors import PairSector, doubled_hamiltonian, start_vector


def blocks(car):
    """Adjacent block `A` of the exchange (seams x seams) and coupling `B` (other pairs x seams), from the seams."""
    deg = car.degree
    eidx = {e: i for i, e in enumerate(car.seams)}
    nb = car.neighbours
    Ar, Ac, Av, lines, cols, ridx = [], [], [], [], [], {}
    for i, (a, b) in enumerate(car.seams):
        Ar.append(i); Ac.append(i); Av.append((deg[a] + deg[b] - 2) / 2)
        for mover, other in ((a, b), (b, a)):
            for t in nb[mover]:
                if t == other:
                    continue
                key = (min(t, other), max(t, other))
                if key in eidx:
                    Ar.append(eidx[key]); Ac.append(i); Av.append(-0.5)
                else:
                    lines.append(ridx.setdefault(key, len(ridx))); cols.append(i)
    E = len(car.seams)
    A = sparse.csr_matrix((Av, (Ar, Ac)), shape=(E, E))
    B = sparse.csr_matrix((np.full(len(lines), -0.5), (lines, cols)), shape=(len(ridx), E))
    return A, B


def _largest(M):
    if M.shape[0] <= 2000:
        return float(np.linalg.eigvalsh(M.toarray())[-1])
    return float(eigsh(M, k=1, which='LA', return_eigenvectors=False, v0=start_vector(M.shape[0]))[0])


def constants(car):
    A, B = blocks(car)
    amax = _largest(A)
    beta = math.sqrt(_largest((B.T @ B).tocsr()))
    return {'adjacent_block_largest_level': amax, 'coupling_largest_singular_value': beta,
            'four_beta_squared': 4 * beta ** 2,
            'amax_below_6': amax <= 6 + 1e-12, 'beta_below_sqrt10': beta <= math.sqrt(10) + 1e-12}


def costs(car):
    """The values of `C` on the sector of two loads: pair on a seam, two separate loads, a doubled load."""
    deg = car.degree
    aa, bb = np.triu_indices(car.n, 1)
    adj = PairSector(car).adjacent
    contrast = deg[aa] + deg[bb] - 2 * adj
    C = contrast + (6 - deg[aa]) + (6 - deg[bb])
    doubled = 4 * deg + 4 * (6 - deg)
    return {'pair_on_a_seam': sorted(set(C[adj].tolist())), 'separate_loads': sorted(set(C[~adj].tolist())),
            'doubled_load': sorted(set(doubled.tolist()))}


def operator_inequalities(car, kappa):
    """Measured sides of the three inequalities (unshifted law) on the full two-load sector."""
    ps = PairSector(car)
    K = ps.exchange()
    adj = np.flatnonzero(ps.adjacent)
    rest = np.flatnonzero(~ps.adjacent)
    PHP_top = 10 * kappa + _largest(K[adj][:, adj])
    KQQ = K[rest][:, rest].tocsr()
    # smallest level of the separate-load block (Lanczos; K >= 0)
    w = eigsh(KQQ, k=1, which='SA', return_eigenvectors=False, v0=start_vector(KQQ.shape[0]))
    sep_min = 12 * kappa + float(w[0])
    dbl_min = float(np.linalg.eigvalsh(doubled_hamiltonian(car, kappa).toarray())[0]) + 12 * kappa
    QHQ_bottom = min(sep_min, dbl_min)
    B = K[rest][:, adj]
    QHP_norm = math.sqrt(_largest((B.T @ B).tocsr()))
    return {'kappa': kappa, 'PHP_largest': PHP_top, 'PHP_bound': 10 * kappa + 6,
            'QHQ_smallest': QHQ_bottom, 'QHQ_bound': 12 * kappa, 'QHP_norm': QHP_norm, 'QHP_bound': math.sqrt(10),
            'holds': PHP_top <= 10 * kappa + 6 + 1e-9 and QHQ_bottom >= 12 * kappa - 1e-9
            and QHP_norm <= math.sqrt(10) + 1e-12}


def single_cell():
    """Exact constants of one cell without degree term (all degrees five): `a = 8 kappa + 5`,
    `c = min(10 kappa + 4 - sqrt6, 20 kappa)`, `beta = 2`."""
    car = cr.cell()
    ps = PairSector(car)
    K = ps.exchange().toarray()
    adj = np.flatnonzero(ps.adjacent)
    rest = np.flatnonzero(~ps.adjacent)
    amax = float(np.linalg.eigvalsh(K[np.ix_(adj, adj)])[-1])
    sep_min = float(np.linalg.eigvalsh(K[np.ix_(rest, rest)])[0])
    beta = float(np.linalg.svd(K[np.ix_(rest, adj)], compute_uv=False)[0])
    dbl_min = float(np.linalg.eigvalsh(car.laplacian().toarray() / 2)[0])
    adjm = car.adjacency().toarray()
    two = adjm @ adjm
    # distance orbits of the nonadjacent pairs: quotient [[3, -sqrt5], [-sqrt5, 5]] with bottom 4 - sqrt6
    quotient = np.array([[3, -math.sqrt(5)], [-math.sqrt(5), 5]])
    return {'adjacent_block_largest_level': amax, 'separate_block_smallest_level': sep_min,
            'coupling_norm': beta, 'doubled_block_smallest_level': dbl_min,
            'adjacent_pairs': int(len(adj)), 'distance_two_pairs': int((two[ps.aa[rest], ps.bb[rest]] > 0).sum()),
            'antipodal_pairs': int((two[ps.aa[rest], ps.bb[rest]] == 0).sum()),
            'amax_is_5': abs(amax - 5) < 1e-12, 'beta_is_2': abs(beta - 2) < 1e-12,
            'separate_minimum_is_4_minus_sqrt6': abs(sep_min - (4 - math.sqrt(6))) < 1e-12,
            'quotient_bottom_is_4_minus_sqrt6': abs(np.linalg.eigvalsh(quotient)[0] - (4 - math.sqrt(6))) < 1e-12,
            'formula': {'a': '8 kappa + 5', 'c': 'min(10 kappa + 4 - sqrt6, 20 kappa)', 'beta': '2'},
            'bound_16_over_gamma_squared': {str(k): (min(1.0, 16 / (2 * k - 1 - math.sqrt(6)) ** 2)
                                                     if 2 * k - 1 - math.sqrt(6) > 0 else 1.0) for k in (4, 8, 128)}}


def emitter_bounds():
    """Exact rationals at `kappa = 128`, `g = 64`, for every seam of the spherical family."""
    beta2, gamma, g = 10, Fraction(250), Fraction(64)
    # (2 beta/gamma + sqrt(10)/(2g))^2 = 10 (2/gamma + 1/(2g))^2 with beta = sqrt 10
    no_pair = beta2 * (2 / gamma + 1 / (2 * g)) ** 2
    band = beta2 * (1 / gamma + 1 / (2 * g)) ** 2
    return {'no_pair_doubled_or_fuel': f'{no_pair.numerator}/{no_pair.denominator}',
            'no_pair_doubled_or_fuel_value': float(no_pair),
            'outside_band_or_fuel': f'{band.numerator}/{band.denominator}',
            'outside_band_or_fuel_value': float(band),
            'theorem_at_kappa_8': 40 / (16 - 6) ** 2, 'theorem_at_kappa_128': 40 / 250 ** 2}


def run(levels=(0,)):
    out = {'single_cell': single_cell(), 'emitter': emitter_bounds()}
    for level in levels:
        car = cr.federation(level)
        entry = constants(car)
        entry['costs'] = costs(car) if car.n <= 1000 else None
        entry['bound_with_carrier_constants'] = {
            str(k): min(1.0, entry['four_beta_squared'] / (2 * k - entry['adjacent_block_largest_level']) ** 2)
            for k in (8, 128)}
        out[f'level{level}'] = entry
    return out
