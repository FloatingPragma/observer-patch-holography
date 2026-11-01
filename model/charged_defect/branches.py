"""Exact integer certificate for the same charged-doublet repair lift.

The state space is 66 unordered occupied block pairs, with both internal
doublets retained.  A movement is a graded block swap; a hold is identity.
The raw records are precisely the original native C/E records.  This is a
supplied field-block service, not the bare H66 quantum channel.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import combinations, product

import numpy as np
from scipy import sparse

from ..exchange_sector import statistics
from ..observer import native
from ..particle_transport import algebra


PAIRS = statistics.PAIRS
INDEX = statistics.INDEX
DIMENSION = 4*len(PAIRS)
# Internal order: up/up, up/down, down/up, down/down.  The columns are an
# UNNORMALIZED singlet followed by the three triplets; all checks are exact.
COUPLING = np.array([[0, 1, 0, 0], [1, 0, 1, 0],
                     [-1, 0, 1, 0], [0, 0, 0, 1]], dtype=np.int64)
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                 [0, 1, 0, 0], [0, 0, 0, 1]], dtype=np.int64)


@dataclass(frozen=True)
class Branch:
    law: str
    seam: tuple[int, int]
    coin: int
    raw: int
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    signs: tuple[int, ...]
    abort: bool

    @property
    def label(self):
        return self.law, self.seam, self.coin, self.raw

    def matrix(self):
        """Integer partial injection; physical scalar amplitude is 1/sqrt(2)."""
        return sparse.coo_matrix((self.signs, (self.outputs, self.inputs)),
                                 shape=(DIMENSION, DIMENSION), dtype=np.int64).tocsr()


def branch(law, seam, coin, raw):
    """Transport the two LABELLED internal states, then reorder occupied ports.

    This reconstruction uses native.native for the record and endpoint rule;
    it does not obtain its signs or internal permutation from signed branches.
    """
    # Reuse strict public-domain validation, including zero raw slots.
    statistics.branch(law, seam, coin, raw, -1)
    rows = []
    for x, pair in enumerate(PAIRS):
        occupied = tuple(int(p in pair) for p in range(12))
        updated, (_, current) = native.native(occupied, seam, coin)
        emitted = {0: 0, 1: 1, -1: 2}[current]
        if law == 'E':
            emitted += 4*((occupied[seam[0]]+occupied[seam[1]]) % 2)
        if raw != emitted:
            continue
        destination = tuple(p for p, n in enumerate(updated) if n)
        transported = list(pair)
        if current:
            permutation = {seam[0]: seam[1], seam[1]: seam[0]}
            transported = [permutation.get(p, p) for p in pair]
        crossing = transported[0] > transported[1]
        if tuple(sorted(transported)) != destination:
            raise ArithmeticError('labelled block transport differs from native endpoint')
        for first, second in product(range(2), repeat=2):
            output_spins = (second, first) if crossing else (first, second)
            rows.append((4*x+2*first+second,
                         4*INDEX[destination]+2*output_spins[0]+output_spins[1],
                         -1 if crossing else 1))
    inputs, outputs, signs = tuple(zip(*rows)) if rows else ((), (), ())
    return Branch(law, seam, coin, raw, inputs, outputs, signs,
                  statistics.branch(law, seam, coin, raw, -1).abort)


def all_branches(law):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('C or E law required')
    return tuple(branch(law, seam, coin, raw) for seam in statistics.SEAMS
                 for coin in (0, 1) for raw in range(4 if law == 'C' else 8))


def _zero(matrix):
    matrix.eliminate_zeros()
    return matrix.nnz == 0


def decomposition_certificate(law):
    """Bind every raw integer column to K+ on singlet and K- on triplet."""
    algebra._native_dictionary_chart()
    coupling = sparse.kron(sparse.eye(66, dtype=np.int64), COUPLING, format='csr')
    counts = np.zeros((30, DIMENSION), dtype=np.int64)
    output_counts = np.zeros_like(counts)
    checks, positive, zero_slots, crossing_columns = 0, 0, 0, 0
    for row in all_branches(law):
        matrix = row.matrix()
        positive += bool(row.inputs)
        zero_slots += not row.inputs
        seam_index = statistics.SEAMS.index(row.seam)
        np.add.at(counts[seam_index], list(row.inputs), 1)
        np.add.at(output_counts[seam_index], list(row.outputs), 1)
        signed = statistics.branch(*row.label, -1)
        rows, cols, values = [], [], []
        for source, target, sign in zip(signed.inputs, signed.outputs, signed.signs, strict=True):
            for internal in range(4):
                rows.append(4*target+internal); cols.append(4*source+internal)
                values.append(1 if internal == 0 else sign)
            crossing_columns += sign == -1
        target = sparse.coo_matrix((values, (rows, cols)), shape=(DIMENSION,)*2).tocsr()
        if not _zero(matrix@coupling-coupling@target):
            raise ArithmeticError('charged internal spin decomposition failed')
        if len(set(row.outputs)) != len(row.outputs):
            raise ArithmeticError('branch is not a partial isometry')
        if len(row.inputs) != 4*len(signed.inputs):
            raise ArithmeticError('some internal columns were dropped')
        checks += len(row.inputs)
    if not np.array_equal(counts, np.full_like(counts, 2)):
        raise ArithmeticError('coin-weighted complete instrument is not trace preserving')
    if not np.array_equal(output_counts, np.full_like(counts, 2)):
        raise ArithmeticError('coin-weighted complete instrument is not unital')
    return {'law': law, 'raw_slots_including_zero': positive+zero_slots,
            'positive_slots': positive, 'zero_slots': zero_slots,
            'integer_input_columns_checked': checks,
            'crossing_spatial_columns': crossing_columns,
            'all_30_seams_complete': True, 'all_30_seams_unital': True,
            'source_dimension': DIMENSION,
            'singlet_spatial_dimension': 66, 'triplet_dimension': 198,
            'triplet_noiseless_factor_dimension': 3,
            'decomposition': 'K_plus tensor I_1 direct_sum K_minus tensor I_3',
            'coupling_gram_diagonal': [2, 1, 2, 1]}


def commutant_certificate(law):
    """An exact full commutant classification, not merely exhibited charges.

    Distinct effect signatures force position-diagonal 4x4 blocks.  A connected
    signed transition graph transports those blocks by conjugation with SWAP.
    One odd signed loop forces [X,SWAP]=0, hence every block is the same member
    of M_1 direct_sum M_3.  The displayed tree/odd edge are proof witnesses.
    """
    rows = statistics.branches(law, -1)
    signatures = [tuple(x in row.inputs for row in rows) for x in range(66)]
    if len(set(signatures)) != 66:
        raise ArithmeticError('complete native effects no longer separate positions')
    edges = {x: {} for x in range(66)}
    for number, row in enumerate(rows):
        for x, y, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            if x == y:
                continue
            edges[x].setdefault(y, (number, sign))
            edges[y].setdefault(x, (number, sign))
    root = 0
    queue, parity, tree = deque([root]), {root: 1}, []
    odd_edge = None
    while queue:
        x = queue.popleft()
        for y, (number, sign) in sorted(edges[x].items()):
            if y not in parity:
                parity[y] = parity[x]*sign
                queue.append(y)
                tree.append([x, y, number, sign])
            elif parity[y] != parity[x]*sign and odd_edge is None:
                odd_edge = [x, y, number, sign]
    if len(parity) != 66 or len(tree) != 65 or odd_edge is None:
        raise ArithmeticError('connected odd-loop branch graph required')
    generators = [np.outer(COUPLING[:, 0], COUPLING[:, 0])]
    generators += [np.outer(COUPLING[:, i], COUPLING[:, j])
                   for i, j in product(range(1, 4), repeat=2)]
    # Exact block checks for every generator and all possible branch factors.
    if any(not np.array_equal(g@SWAP, SWAP@g) for g in generators):
        raise ArithmeticError('claimed protected operator fails a branch')
    # Integer linear independence: the coupled-basis image is supported in
    # ten distinct matrix entries, including the singlet diagonal.
    coupled = [COUPLING.T@g@COUPLING for g in generators]
    supports = [tuple(map(tuple, np.argwhere(g))) for g in coupled]
    if any(len(s) != 1 for s in supports) or len(set(supports)) != 10:
        raise ArithmeticError('protected algebra basis is not independent')
    return {'law': law, 'separated_positions': 66,
            'separated_position_pairs': len(tuple(combinations(range(66), 2))),
            'transition_tree': tree, 'odd_chord': odd_edge,
            'tree_parities': [parity[x] for x in range(66)],
            'commutant_complex_dimension': 10,
            'commutant': 'I_66 tensor (M_1 direct_sum M_3)',
            'generated_star_algebra_complex_dimension': 2*66*66,
            'generated_star_algebra': 'M_66 direct_sum (M_66 tensor I_3)',
            'maximal_noiseless_factor_dimension': 3,
            'whole_four_state_internal_factor_protected': False}


def report():
    return {'dimensions': {'spatial': 66, 'internal': 4, 'total': DIMENSION},
            'laws': {law: {'decomposition': decomposition_certificate(law),
                           'commutant': commutant_certificate(law)} for law in ('C', 'E')}}
