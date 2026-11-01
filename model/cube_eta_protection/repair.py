"""Actual four-dimer native block service on the complete65536-state cube.

Only the existing seams whose endpoints lie in ports(0,1,5,7) are admitted.
This realizes the unchanged local rule on four blocks, not twelve blocks
or an invented four-port primitive. Operations are exact low-rank maps plus
the complete outside projector on the full cube Fock space.
"""
from functools import lru_cache
from itertools import combinations, product
import math

import numpy as np
from scipy import sparse
import sympy as sy

from ..particle_transport import algebra as local
from ..observer import native
from ..exchange_sector import statistics
from . import cube


PORTS = (0, 1, 5, 7)
PAIRS = tuple(combinations(range(4), 2))
SEAMS = tuple(pair for pair in PAIRS if tuple(sorted(PORTS[j] for j in pair)) in statistics.SEAMS)


def permute_mask(mask, permutation):
    occupied = [j for j in range(len(permutation)) if mask >> j & 1]
    mapped = [permutation[j] for j in occupied]
    crossings = sum(mapped[i] > mapped[j] for i in range(len(mapped)) for j in range(i+1, len(mapped)))
    return sum(1 << j for j in mapped), (-1)**crossings


def relocate(columns, permutation):
    coo = columns.tocoo()
    mapping = {int(row): permute_mask(int(row), permutation) for row in set(coo.row)}
    rows = [mapping[int(row)][0] for row in coo.row]
    values = [value*mapping[int(row)][1] for row, value in zip(coo.row, coo.data, strict=True)]
    return sparse.coo_matrix((values, (rows, coo.col)), shape=columns.shape).tocsr()


def block_to_cube():
    return tuple(2*(4*x+block)+spin for block in range(4) for x in range(2) for spin in range(2))


@lru_cache(maxsize=1)
def exact_embedding():
    """Unnormalized columns A+sqrt17 B in physical cube CAR order."""
    g0 = np.zeros(16, complex); g1 = np.zeros(16, complex)
    g0[5], g0[10], g0[3], g0[12] = -1j, 1j, 4, 4
    g1[3] = g1[12] = -1
    hole = np.zeros((16, 2), complex)
    hole[1, 0] = hole[2, 1] = 1
    hole[8, 0] = hole[4, 1] = 1j
    parts = [[], []]
    for pair in PAIRS:
        for spins in product(range(2), repeat=2):
            selected = dict(zip(pair, spins, strict=True))
            a = sparse.csr_matrix([[1]], dtype=complex)
            b = sparse.csr_matrix([[0]], dtype=complex)
            for block in reversed(range(4)):
                p = sparse.csr_matrix((hole[:, selected[block]] if block in selected else g0)[:, None])
                q = sparse.csr_matrix((np.zeros(16) if block in selected else g1)[:, None])
                a, b = (sparse.kron(a, p)+17*sparse.kron(b, q),
                        sparse.kron(a, q)+sparse.kron(b, p))
            parts[0].append(a); parts[1].append(b)
    return tuple(relocate(sparse.hstack(part, format='csr'), block_to_cube()) for part in parts)


@lru_cache(maxsize=1)
def embedding():
    # Bind the algebraic columns to the actual inherited normalized ground
    # and hole columns; do not choose eigenvector phases numerically.
    ground, hole = local.ground_column(), local.candidate_columns('hole')
    columns = []
    for pair in PAIRS:
        for spins in product(range(2), repeat=2):
            selected = dict(zip(pair, spins, strict=True)); col = sparse.csr_matrix([[1]], dtype=complex)
            for block in reversed(range(4)):
                factor = hole[:, selected[block]] if block in selected else ground
                col = sparse.kron(col, sparse.csr_matrix(factor[:, None]), format='csr')
            columns.append(col)
    return relocate(sparse.hstack(columns, format='csr'), block_to_cube())


def matrices(law, seam):
    if law not in ('C', 'E') or seam not in SEAMS:
        raise ValueError('existing C/E law and an admitted native four-block seam required')
    slots = []
    edge = tuple(sorted(PORTS[j] for j in seam))
    for coin in (0, 1):
        for raw in range(4 if law == 'C' else 8):
            matrix = np.zeros((24, 24), dtype=np.int64)
            for p, pair in enumerate(PAIRS):
                state = tuple(int(port in (PORTS[pair[0]], PORTS[pair[1]])) for port in range(12))
                output, (_, current) = native.native(state, edge, coin)
                record = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E': record += 4*((state[edge[0]]+state[edge[1]]) % 2)
                if record != raw: continue
                destination = tuple(j for j, port in enumerate(PORTS) if output[port])
                labels = ([seam[1] if j == seam[0] else seam[0] if j == seam[1] else j for j in pair]
                          if current else list(pair))
                for a, b in product(range(2), repeat=2):
                    crossed = labels[0] > labels[1]
                    spin = (b, a) if crossed else (a, b)
                    matrix[4*PAIRS.index(destination)+2*spin[0]+spin[1], 4*p+2*a+b] = -1 if crossed else 1
            slots.append((('native', coin, raw), matrix))
    return tuple(slots)


def apply(state, law, seam):
    """Every native outcome and Q on all65536 input rows; never normalize."""
    if state.shape[0] != cube.DIMENSION:
        raise ValueError('complete cube Fock input rows required')
    v = embedding(); compressed = v.conj().T@state
    out = [(label, v@(matrix@compressed)/math.sqrt(2)) for label, matrix in matrices(law, seam)]
    out.append((('outside',), state-v@compressed))
    return tuple(out)


def certificate():
    a, b = exact_embedding(); eta = cube.operators().eta_minus.real
    if (eta@a).nnz or (eta@b).nnz:
        raise ArithmeticError('actual cube repair support is not an exact eta-lowest-weight space')
    # Close the analytic normalization premise, not merely a rounded overlap.
    radical = sy.sqrt(17)
    ratio = 4-radical
    source_a_squared = (1+4/radical)/2
    source_b_over_a = -sy.Rational(1, 8)/(radical/4*source_a_squared)
    if (sy.simplify(source_b_over_a-ratio) != 0
            or sy.simplify(source_a_squared*(1+ratio**2)-1) != 0):
        raise ArithmeticError('consumed analytic ground formula differs from the exact normalized column')
    if ((a.conj().T@a+17*b.conj().T@b-35904*sparse.eye(24)).nnz
            or (a.conj().T@b+b.conj().T@a+8704*sparse.eye(24)).nnz):
        raise ArithmeticError('exact sqrt17 column Gram normalization failed')
    v = embedding()
    denominator = 4+4*(4-math.sqrt(17))**2  # two holes and two ground norm factors
    binding = float(sparse.linalg.norm(v-(a+math.sqrt(17)*b)/denominator))
    if binding > 3e-13 or sparse.linalg.norm(v.conj().T@v-sparse.eye(24)) > 3e-13:
        raise ArithmeticError('actual normalized cube block dictionary differs from exact columns')
    count = 0
    for law in ('C', 'E'):
        for seam in SEAMS:
            slots = matrices(law, seam)
            if (not np.array_equal(sum(m.T@m for _, m in slots), 2*np.eye(24))
                    or not np.array_equal(sum(m@m.T for _, m in slots), 2*np.eye(24))):
                raise ArithmeticError('native complete four-block law lost completeness or unitality')
            # The abstract matrices must be the physical graded permutation
            # of the actual block columns, including fermionic crossing signs.
            permutation = tuple(2*(4*x+(seam[1] if block == seam[0] else seam[0] if block == seam[1] else block))+spin
                                for x in range(2) for block in range(4) for spin in range(2))
            for label, matrix in slots:
                support = np.flatnonzero(np.any(matrix != 0, axis=0))
                for exact in (a, b):
                    expected = (relocate(exact[:, support], permutation)
                                if label[2] % 4 in (1, 2) else exact[:, support])
                    if (expected-exact@sparse.csr_matrix(matrix[:, support])).nnz:
                        raise ArithmeticError('native repair matrix differs from the physical graded block law')
            count += len(slots)+1
    # Orthogonality follows factorwise for every image column, not only a probe.
    x = cube.eta_image()
    if (a[cube.weights(6)].conj().T@x).nnz or (b[cube.weights(6)].conj().T@x).nnz:
        raise ArithmeticError('a complete cube eta-image column reaches repair support')
    return {'native_block_ports': list(PORTS),
            'admitted_native_seams': [[PORTS[j] for j in pair] for pair in SEAMS],
            'x_directed_dimers': [[j, j+4] for j in range(4)],
            'reference': 'unchanged axis0 dimer ground on each unoccupied block',
            'charged_kind': 'unchanged lower hole doublet, exactly two occupied blocks',
            'repair_support_dimension': 24, 'repair_support_CAR_number': 6,
            'exact_support_eta': 1, 'exact_eta_minus_support_zero': True,
            'complete_outcomes_including_zero_and_outside': count,
            'protected_full_eta_image_columns': 1820,
            'all_native_outcomes_annihilate_protected_image': True,
            'outside_Q_is_identity_on_protected_image': True,
            'complete_reference_preservation': 'the complete channel is identity on all operators supported in the eta image, including arbitrary reference cross terms',
            'exact_column_normalization': 'n=136-32sqrt17; Gram=(35904-8704sqrt17)I=n^2 I',
            'analytic_ground_normalization_proved': True,
            'normalized_source_binding_residual': binding,
            'binding_diagnostic_admission': 3e-13,
            'complete_full_Fock_domain': cube.DIMENSION,
            'twelve_block_source_execution_claimed': False,
            'spin_frame_dressing_or_modified_repair_used': False}
