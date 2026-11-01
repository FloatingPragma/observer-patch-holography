"""Same-action SU(2) premises and complete unchanged repair interfaces.

The chain's hidden spin has identical two-site block frames.  Local ground
projectors are singlets, both charged projectors are doublets, and graded
whole-block swaps commute with the sum of those local spins.  These exact
facts prove the full-Fock result without projecting onto the leaky original
charged carrier.  Three-dimensional frame mismatch is retained separately.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations, product

import numpy as np
from scipy import sparse

from ..charged_defect import dictionary
from ..field_symmetry import symmetry
from ..field_stress.energy import GaussianInteger as GI, commutator
from ..full_fock_field.field import FieldGeometry
from ..particle_transport import algebra


GEOMETRY = FieldGeometry((6, 1, 1), speed=1/8, interaction=1.)
BLOCK_PORTS = (0, 1, 5)  # A literal native triangle; local indices are 0,1,2.
NUMERICAL_ADMISSION = 2e-11


@lru_cache(maxsize=1)
def operators():
    return symmetry.operators(GEOMETRY)


def local_certificate():
    """Exact Gaussian-integer premises for arbitrary tensor powers of dimers."""
    local = symmetry.operators(FieldGeometry((2, 1, 1), speed=1/8, interaction=1.))
    singlets = np.zeros((16, 2), complex)
    singlets[5, 0], singlets[10, 0] = -1j, 1j
    singlets[3, 1] = singlets[12, 1] = 1
    if any(np.any(s.matrix()@singlets) for s in local.spin_twice):
        raise ArithmeticError('the two analytic ground-generating columns are not singlets')
    ground = algebra.ground_column()
    if np.linalg.norm(ground-singlets@(singlets.conj().T@ground)/2) > 1e-14:
        raise ArithmeticError('the consumed ground left its exact singlet span')
    projector_checks = 0
    for kind, entries in (
        ('hole', ((1, 0, 1), (2, 1, 1), (4, 1, 1j), (8, 0, 1j))),
        ('particle', ((7, 0, -1j), (11, 1, 1j), (13, 1, -1), (14, 0, 1)))):
        column = np.zeros((16, 2), complex)
        for row, index, value in entries:
            column[row, index] = value
        if np.linalg.norm(column/np.sqrt(2)-algebra.candidate_columns(kind)) > 1e-14:
            raise ArithmeticError('the unchanged charged columns left their exact chart')
        numerator = GI.from_exact(column@column.conj().T)
        for spin in local.spin_twice:
            if not commutator(spin, numerator).is_zero():
                raise ArithmeticError('charged doublet projector is not SU(2) invariant')
            projector_checks += 1
    # Whole-block graded swaps commute with sums of identical even local spins.
    swap = GI.from_exact(algebra.block_swap())
    for spin in local.spin_twice:
        total = GI.from_exact(sparse.kron(sparse.eye(16), spin.matrix())+
                              sparse.kron(spin.matrix(), sparse.eye(16)))
        if not commutator(swap, total).is_zero():
            raise ArithmeticError('unchanged graded block swap breaks the common spin')
    # The global chain frame repeats every two sites, including its signs.
    for site in range(24):
        for axis in range(3):
            if not np.array_equal(symmetry.hidden_spin_matrix((site, 0, 0), axis),
                                  symmetry.hidden_spin_matrix((site % 2, 0, 0), axis)):
                raise ArithmeticError('identical dimer frames fail on the native twelve-block chain')
    return {'analytic_ground_generating_singlets': 2,
            'exact_charged_projector_commutators': projector_checks,
            'exact_complete_256_state_swap_commutators': 3,
            'identical_two_site_spin_frames': True,
            'formal_twelve_block_chain_sites': 24,
            'formal_twelve_block_source_executed': False}


def repair_outcomes(kind='hole', law='E', seam=(0, 1)):
    """Use the predecessor's complete operator definitions without alterations."""
    return dictionary.full_extension(3, kind, 0, seam, law)


def full_certificate():
    ops = operators()
    for spin in ops.spin_twice:
        if not commutator(ops.K, spin).is_zero():
            raise ArithmeticError('the joined full-field generator breaks hidden spin')
    count = 0
    largest = 0.
    for kind, law, seam in product(('hole', 'particle'), ('C', 'E'), combinations(range(3), 2)):
        outcomes = repair_outcomes(kind, law, seam)
        complete = sparse.csr_matrix((4096, 4096), dtype=complex)
        for _, operator in outcomes:
            complete += operator.conj().T@operator
            for spin in ops.spin_twice:
                s = spin.matrix()
                largest = max(largest, float(sparse.linalg.norm(operator@s-s@operator)))
            count += 1
        if sparse.linalg.norm(complete-sparse.eye(4096)) > NUMERICAL_ADMISSION:
            raise ArithmeticError('complete full-Fock repair lost an outcome')
    if not np.isfinite(largest) or largest > NUMERICAL_ADMISSION:
        raise ArithmeticError('actual complete repair operator breaks the common spin')
    # This new common spin does not repair the missing joined-energy account.
    generator = ops.K.matrix(ops.scale)
    pulled = sum((operator.conj().T@generator@operator
                  for _, operator in repair_outcomes()), start=sparse.csr_matrix(generator.shape))
    energy_change = float(sparse.linalg.norm(pulled-generator))
    if not np.isfinite(energy_change) or energy_change < 1e-3:
        raise ArithmeticError('joined-energy disturbance control disappeared')
    return {'sites': 6, 'full_fock_dimension': 4096,
            'native_triangle_ports': list(BLOCK_PORTS),
            'exact_full_field_spin_commutators': 3,
            'complete_repair_outcomes_checked': count,
            'outcome_spin_commutators_checked': 3*count,
            'largest_repair_spin_Frobenius_residual': largest,
            'numerical_admission': NUMERICAL_ADMISSION,
            'joined_energy_pullback_Frobenius_change': energy_change,
            'joined_energy_conserved_by_repair': False,
            'outside_outcome_retained': True}


def frame_obstruction():
    """Exact cube-tiling witness against extending unchanged swaps to 3D.

    Four x-directed dimers have bases (0,y,z), with y,z in {0,1}.  The
    original charged columns are identical in their local frames.  Global
    hidden-spin axes have different signs on these four blocks.  We use the
    actual native seams among ports (0,1,5,7); two movements already destroy
    every linear combination of the three global hidden-spin generators.
    """
    bases = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1))
    ports = (0, 1, 5, 7)
    pairs = tuple(combinations(range(4), 2))
    index = {pair: k for k, pair in enumerate(pairs)}
    spins = []
    for axis in range(3):
        sigma = np.asarray(symmetry.SIGMA[axis], dtype=complex)
        signs = [(-1)**sum(site[j] for j in range(3) if j != axis) for site in bases]
        spins.append(sparse.block_diag([sparse.csr_matrix(signs[a]*np.kron(sigma, np.eye(2))+
                                        signs[b]*np.kron(np.eye(2), sigma))
                                        for a, b in pairs]).toarray())
    gram = np.zeros((3, 3), complex)
    witnesses = []
    for seam in ((0, 1), (0, 2)):
        # coin=0, current=+1 => raw E flag 5, moving left occupied to right.
        matrix = np.zeros((24, 24), complex)
        left, right = seam
        for x, pair in enumerate(pairs):
            if left not in pair or right in pair:
                continue
            labels = [right if p == left else p for p in pair]
            crossing = labels[0] > labels[1]
            y = index[tuple(sorted(labels))]
            for a, b in product(range(2), repeat=2):
                out = (b, a) if crossing else (a, b)
                matrix[4*y+2*out[0]+out[1], 4*x+2*a+b] = -1 if crossing else 1
        commutators = [matrix@s-s@matrix for s in spins]
        for i, j in product(range(3), repeat=2):
            gram[i, j] += np.vdot(commutators[i], commutators[j])
        witnesses.append({'native_seam': [ports[left], ports[right]], 'coin': 0, 'raw_E': 5})
    if np.any(gram.imag) or np.any(gram.real != np.rint(gram.real)):
        raise ArithmeticError('cube obstruction lost exact integer arithmetic')
    integer = gram.real.astype(np.int64)
    if np.any(integer != np.diag(np.diag(integer))) or np.any(np.diag(integer) <= 0):
        raise ArithmeticError('cube movements no longer exclude the full linear spin span')
    return {'geometry': [2, 2, 2], 'x_dimer_base_sites': [list(x) for x in bases],
            'native_block_ports': list(ports), 'witnesses': witnesses,
            'compressed_two_doublet_dimension': 24,
            'exact_commutator_Gram': integer.tolist(),
            'surviving_global_spin_generator_linear_dimension': 0,
            'all_possible_noiseless_algebras_excluded': False,
            'frame_dressed_repair_silently_substituted': False}


def aligned_translation(displacement):
    """All hidden-spin axes match iff the three displacement parities agree."""
    if (type(displacement) is not tuple or len(displacement) != 3
            or any(type(x) is not int for x in displacement)):
        raise ValueError('three integer displacement components required')
    return displacement[0] % 2 == displacement[1] % 2 == displacement[2] % 2


def alignment_certificate():
    # The exact formula S_a(x)=(-1)^sum_{j!=a}x_j sigma_a makes this finite
    # parity calculation a proof for every integer x and displacement.
    rows = []
    for displacement in product(range(2), repeat=3):
        signs = [(-1)**sum(displacement[j] for j in range(3) if j != axis)
                 for axis in range(3)]
        if all(sign == 1 for sign in signs) != aligned_translation(displacement):
            raise ArithmeticError('translation parity characterization failed')
        for site in product(range(4), repeat=3):
            translated = tuple(x+y for x, y in zip(site, displacement, strict=True))
            for axis in range(3):
                if not np.array_equal(symmetry.hidden_spin_matrix(translated, axis),
                                      signs[axis]*symmetry.hidden_spin_matrix(site, axis)):
                    raise ArithmeticError('literal onsite frames violate the alignment criterion')
        rows.append({'displacement_parity': list(displacement),
                     'spin_axis_signs': signs, 'aligned': aligned_translation(displacement)})
    return {'exact_parity_cases': rows,
            'criterion': 'delta_0 = delta_1 = delta_2 modulo 2',
            'all_integer_coordinates_and_displacements': True,
            'even_side_block_tilings_satisfy_alignment': True,
            'odd_odd_odd_displacements_also_align': True,
            'conditional_block_requirements': [
                'identical internal action and coordinate-identified block operations',
                'reference and candidate support projectors commute with every local hidden-spin generator',
                'every complete repair Kraus operator acts through those invariant projectors and aligned block permutations'],
            'new_cube_candidate_or_source_execution_supplied': False}
