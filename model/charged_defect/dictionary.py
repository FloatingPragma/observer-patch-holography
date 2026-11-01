"""Attach the exact branch algebra to the unchanged charged dimer spectrum.

The local algebra identities are numerical checks of explicit analytic
columns.  The arbitrary-block statements follow by tensor products and the
exact graded permutation law.  This file does not execute a 48-mode source.
"""
from __future__ import annotations

from itertools import combinations, product
from fractions import Fraction
import math

import numpy as np
from scipy import sparse

from ..particle_transport import algebra
from ..full_fock_field.field import FieldGeometry, FullFock


def embedded_columns(blocks=3, kind='hole', axis=0):
    """Small full-Fock realization, with every internal column retained.

    Formal twelve-block identities are proved factorwise; allocating that
    full space here is deliberately outside the evaluator's bounded domain.
    """
    if type(blocks) is not int or blocks not in (2, 3):
        raise ValueError('two or three complete dimers required')
    if type(kind) is not str or kind not in ('hole', 'particle'):
        raise ValueError('one identical charged doublet kind required')
    ground = algebra.ground_column(axis)
    candidates = algebra.candidate_columns(kind, axis)
    result = []
    for pair in combinations(range(blocks), 2):
        for spins in product(range(2), repeat=2):
            factors = [ground for _ in range(blocks)]
            for occupied, spin in zip(pair, spins, strict=True):
                factors[occupied] = candidates[:, spin]
            column = np.ones(1, dtype=complex)
            for factor in reversed(factors):
                column = np.kron(column, factor)
            result.append(column)
    return np.column_stack(result)


def full_extension(blocks=3, kind='hole', axis=0, seam=(0, 1), law='E'):
    """Complete full-Fock instrument, including an outside-sector outcome.

    Returns each coin/raw outcome (scaled by 1/sqrt(2)) and Q=I-P.  Outside
    states are retained in their own outcome, not discarded or normalized.
    This small exact-definition evaluator is not native source execution.
    """
    if type(blocks) is not int or blocks not in (2, 3):
        raise ValueError('two or three complete dimers required')
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(x) is not int for x in seam)
            or not 0 <= seam[0] < seam[1] < blocks):
        raise ValueError('ordered distinct seam in the declared blocks required')
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('C or E required')
    embedding = sparse.csr_matrix(embedded_columns(blocks, kind, axis))
    pairs = tuple(combinations(range(blocks), 2))
    dimension = 1 << (4*blocks)
    permutation = list(range(blocks))
    permutation[seam[0]], permutation[seam[1]] = seam[1], seam[0]
    targets, phases = zip(*(algebra.block_permutation_address(x, tuple(permutation))
                           for x in range(dimension)), strict=True)
    swap = sparse.coo_matrix((phases, (targets, range(dimension))),
                             shape=(dimension, dimension)).tocsr()
    outputs = []
    for coin in (0, 1):
        for raw in range(4 if law == 'C' else 8):
            support = []
            moving = None
            for x, pair in enumerate(pairs):
                left, right = int(seam[0] in pair), int(seam[1] in pair)
                current = ((left-right)-(2*coin-1)*((left+right) % 2))//2
                record = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E':
                    record += 4*((left+right) % 2)
                if record == raw:
                    support.extend(range(4*x, 4*x+4))
                    if moving is not None and moving != bool(current):
                        raise ArithmeticError('one raw record mixes hold and movement')
                    moving = bool(current)
            selected = embedding[:, support]
            projector = selected@selected.conj().T
            outputs.append((('native', coin, raw),
                            (swap@projector if moving else projector)/math.sqrt(2)))
    projector = embedding@embedding.conj().T
    outputs.append((('outside',), sparse.eye(dimension, format='csr')-projector))
    return tuple(outputs)


def tensor_generator(blocks, axis=0):
    if type(blocks) is not int or blocks not in (2, 3):
        raise ValueError('two or three complete dimers required')
    local = sparse.csr_matrix(algebra.setup(axis).generator)
    return sum((sparse.kron(sparse.eye(16**(blocks-j-1)),
                 sparse.kron(local, sparse.eye(16**j)), format='csr')
                for j in range(blocks)), start=sparse.csr_matrix((16**blocks,)*2))


def report():
    residuals = {}
    for axis in range(3):
        generator = algebra.setup(axis).generator
        ground = algebra.ground_column(axis)
        numbers = np.diag([x.bit_count() for x in range(16)])
        for kind, charge in (('hole', -1), ('particle', 1)):
            columns = algebra.candidate_columns(kind, axis)
            errors = [np.linalg.norm(columns.conj().T@columns-np.eye(2)),
                      np.linalg.norm(ground.conj()@columns),
                      np.linalg.norm(generator@columns-(algebra.GROUND_ENERGY+algebra.LOWER_GAP)*columns),
                      np.linalg.norm((numbers-2*np.eye(16))@columns-charge*columns)]
            residual = float(max(errors))
            if not math.isfinite(residual) or residual > 1e-12:
                raise ArithmeticError('same charged dimer dictionary failed')
            residuals[f'{kind}_axis{axis}'] = residual
    # A genuinely different law: uniform four-site nearest-neighbor field.
    # It contains an inter-dimer hopping term and does not preserve this sector.
    field = FullFock(8)
    geometry = FieldGeometry((4, 1, 1), speed=algebra.SPEED, interaction=algebra.INTERACTION)
    uniform = field.hamiltonian(geometry)
    controls = {}
    for kind in ('hole', 'particle'):
        embedding = embedded_columns(2, kind)
        image = uniform@embedding
        leakage = float(np.linalg.norm(image-embedding@(embedding.conj().T@image)))
        if not math.isfinite(leakage) or leakage < 1e-3:
            raise ArithmeticError('uniform joined-field leakage control disappeared')
        # Here every entry is Gaussian dyadic: W=2V and A=16H are integer.
        # B=(4I-WW*)AW and the leakage is B/128.  The explicit magnitude
        # cap keeps every integer multiply/add exactly within binary64.
        local = np.zeros((16, 2), complex)
        entries = (((1, 0, 1), (2, 1, 1), (4, 1, 1j), (8, 0, 1j))
                   if kind == 'hole' else
                   ((7, 0, -1j), (11, 1, 1j), (13, 1, -1), (14, 0, 1)))
        for row, column, value in entries:
            local[row, column] = value
        w = np.column_stack([np.kron(local[:, j], local[:, i])
                             for i, j in product(range(2), repeat=2)])
        a = 16*uniform.toarray()
        if np.linalg.norm(w/2-embedding) > 1e-14:
            raise ArithmeticError('analytic charged-column chart changed')
        if (np.max(abs(w)) > 1 or np.max(abs(a)) > 100
                or not np.array_equal(w.real, np.rint(w.real))
                or not np.array_equal(w.imag, np.rint(w.imag))
                or not np.array_equal(a.real, np.rint(a.real))
                or not np.array_equal(a.imag, np.rint(a.imag))):
            raise ArithmeticError('exact dyadic leakage premise failed')
        b = 4*a@w-w@(w.conj().T@a@w)
        singlet = np.array([0, 1, -1, 0])
        exact_gram = 32*(np.eye(4, dtype=int)+np.outer(singlet, singlet))
        if not np.array_equal(b.conj().T@b, exact_gram):
            raise ArithmeticError('uniform field singlet/triplet leakage identity failed')
        numerator = sum(int(z.real)**2+int(z.imag)**2 for z in b.ravel())
        fraction = Fraction(numerator, 128**2)
        if fraction != Fraction(3, 256):
            raise ArithmeticError('exact joined-field leakage changed')
        controls[kind] = {'Frobenius_norm': leakage,
                          'exact_squared_Frobenius_norm': str(fraction),
                          'leakage_gram': '(3 P_singlet + P_triplet)/512',
                          'singlet_normalized_squared_leakage': '3/512',
                          'triplet_normalized_squared_leakage': '1/512',
                          'nonzero_nonleaking_input_exists': False}
    extension = {}
    generator = tensor_generator(3)
    for law in ('C', 'E'):
        operations = full_extension(law=law)
        total = sum((v.conj().T@v for _, v in operations),
                    start=sparse.csr_matrix(generator.shape))
        completeness = float(sparse.linalg.norm(total-sparse.eye(4096)))
        energy = max(float(sparse.linalg.norm(generator@v-v@generator)) for _, v in operations)
        if not all(math.isfinite(x) and x < 1e-12 for x in (completeness, energy)):
            raise ArithmeticError('complete full-Fock extension failed')
        extension[law] = {'evaluated_full_fock_dimension': 4096,
                          'outcomes_including_zero_and_outside': len(operations),
                          'completeness_Frobenius_residual': completeness,
                          'energy_commutator_max_Frobenius_residual': energy,
                          'numerical_admission': 1e-12}
    return {
        'candidate': 'same lowest charged dimer doublet on every occupied block',
        'formal_dimer_blocks': 12, 'formal_full_field_modes': 48,
        'formal_full_fock_dimension': 2**48,
        'two_defect_carrier_dimension': 264,
        'same_dimer_gap': algebra.LOWER_GAP,
        'two_defect_reference_subtracted_generator': 2*algebra.LOWER_GAP,
        'gap_formula': 'sqrt(17)/8 - 1/16',
        'charge_per_occupied_block': {'hole': -1, 'particle': 1},
        'unoccupied_block_density': 'exact unchanged half-filled dimer ground',
        'native_current_conversion': 'excess CAR number q times the native defect current',
        'selected_pair_total_charge': {'hole': -2, 'particle': 2},
        'local_dictionary_max_residuals': residuals,
        'local_dictionary_admission': 1e-12,
        'joined_uniform_field_outside_sector_leakage': controls,
        'small_full_extension_checks': extension,
        'full_extension': 'all native coin/raw outcomes on P plus Q=I-P outside outcome',
        'full_extension_energy_conserved': True,
        'two_defect_carrier_energy_spread': 0,
        'distinct_internal_gaps_at_fixed_charge': 1,
        'common_scale_or_chemical_potential_splits_triplet': False,
        'new_absolute_physical_masses_predicted': 0,
        'triplet_is_three_generations': False,
        'binding_or_topological_charge_derived': False,
        'bare_H66_channel_protection': False,
        'full_48_mode_source_executed': False,
    }
