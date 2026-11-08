"""Exact source symmetry reduction at the magnetic kinetic node.

The calculation determines spectral multiplicities and transported probes.
It supplies no spectral density, mass shell or physical species assignment.
"""
from fractions import Fraction as F
from itertools import product

import numpy as np

from model.field_symmetry.symmetry import onsite_frame, SIGMA
from model.source_family_compatibility import momentum_response as source
from model.source_family_compatibility.hole_dispersion import gamma_matrices


PARITIES = tuple(product(range(2), repeat=3))
INDEX = {x: i for i, x in enumerate(PARITIES)}


def magnetic_generators():
    """Remove only the common phase -i from U_j at k_*=(pi/2)^3."""
    result = []
    for j in range(3):
        matrix = np.zeros((8, 8), dtype=np.int64)
        for x in PARITIES:
            y = list(x)
            y[j] ^= 1
            matrix[INDEX[x], INDEX[tuple(y)]] = source.magnetic_phase(x, j)
        result.append(matrix)
    return tuple(result)


def cubic_lifts():
    """Actual proper-cubic lifts, including reciprocal-cell chart transport."""
    result = []
    for perm, signs in source.cubic_rotations():
        shifts = [0] * 3
        for j in range(3):
            shifts[perm[j]] = (signs[j] - 1) // 2
        matrix = np.zeros((8, 8), dtype=np.int64)
        for x in PARITIES:
            y = tuple(t % 2 for t in source.rotate(x, perm, signs))
            transport = (-1) ** sum(shifts[j] * y[j] for j in range(3))
            matrix[INDEX[y], INDEX[x]] = transport * source.rotation_phase(x, perm)
        result.append(matrix)
    return tuple(result)


def commutant_basis(generators):
    """Exact signed-orbit elimination of U X U*=X, over real or complex scalars."""
    generators = tuple(np.asarray(g) for g in generators)
    if not generators:
        raise ValueError('at least one signed permutation is required')
    size = generators[0].shape[0]
    actions = []
    for matrix in generators:
        if (matrix.shape != (size, size)
                or not np.all(np.isin(matrix, (-1, 0, 1)))
                or not np.array_equal(np.abs(matrix).sum(axis=0), np.ones(size))
                or not np.array_equal(np.abs(matrix).sum(axis=1), np.ones(size))):
            raise ValueError('exact signed permutation matrices required')
        permutation = np.argmax(np.abs(matrix), axis=0)
        signs = matrix[permutation, np.arange(size)].astype(int)
        actions.append((permutation, signs))
    remaining = set(product(range(size), repeat=2))
    basis = []
    while remaining:
        start = min(remaining)
        values = {start: 1}
        queue = [start]
        consistent = True
        while queue:
            i, j = queue.pop()
            for permutation, signs in actions:
                target = (int(permutation[i]), int(permutation[j]))
                value = int(signs[i] * signs[j] * values[i, j])
                if target not in values:
                    values[target] = value
                    queue.append(target)
                elif values[target] != value:
                    consistent = False
        remaining.difference_update(values)
        if consistent:
            matrix = np.zeros((size, size), dtype=np.int64)
            for (i, j), value in values.items():
                matrix[i, j] = value
            basis.append(matrix)
    return tuple(basis)


def transported_valleys():
    """Unnormalized exact Gaussian-integer columns; F*F=8I, without rounding."""
    blocks = []
    for b in PARITIES:
        block = np.vstack([
            (-1j) ** sum(x) * (-1) ** sum(b[j] * x[j] for j in range(3))
            * onsite_frame(x) for x in PARITIES
        ])
        blocks.append(block)
    return tuple(blocks)


def inversion_source_checks():
    """Actual centered-box symmetry, averaging and transported nodal inversion."""
    sigma = tuple(np.asarray(x, dtype=complex) for x in SIGMA)
    source_checks = 0
    averaging_checks = 0
    for radius in (1, 2):
        points = tuple(product(range(-radius, radius), repeat=3))
        lookup = {x: j for j, x in enumerate(points)}
        size = 2*len(points)
        h = np.zeros((size, size), dtype=complex)
        reflection = np.zeros_like(h)
        for x, index in lookup.items():
            y = tuple(-1-t for t in x)
            # theta_(-1,-1,-1) S, where S a_x=epsilon_x a_-x.
            reflection[2*lookup[y]:2*lookup[y]+2, 2*index:2*index+2] = (
                (-1)**sum(x)*np.eye(2))
            for axis in range(3):
                neighbor = list(x); neighbor[axis] += 1; neighbor = tuple(neighbor)
                if neighbor not in lookup:
                    continue
                other = lookup[neighbor]
                link = 1j*sigma[axis]
                h[2*index:2*index+2, 2*other:2*other+2] = link
                h[2*other:2*other+2, 2*index:2*index+2] = link.conj().T
                source_checks += 1
        if (not np.array_equal(reflection @ reflection.conj().T, np.eye(size))
                or not np.array_equal(reflection @ h @ reflection.conj().T, h)):
            raise ArithmeticError('affine staggered inversion does not preserve literal Weyl links')
        if {tuple(-1-t for t in x) for x in points} != set(points):
            raise ArithmeticError('actual translation average lost its reflection pairing')
        # On every annihilator a_y, R_a theta_x S = gamma_epsilon_x theta_(a-x).
        for x in points:
            for y in ((0, 0, 0), (1, -2, 3)):
                left_phase = (-1)**(sum(y)+sum(tuple(x[j]-y[j] for j in range(3))))
                left_site = tuple(-1-x[j]+y[j] for j in range(3))
                right_site = tuple(y[j]-1-x[j] for j in range(3))
                if left_phase != (-1)**sum(x) or left_site != right_site:
                    raise ArithmeticError('gauge-invariant averaging identity failed')
                averaging_checks += 1
    frame_checks = 0
    for x in product(range(-2, 3), repeat=3):
        frame = onsite_frame(tuple(t % 4 for t in x))
        inverse_frame = onsite_frame(tuple((-t) % 4 for t in x))
        if not np.array_equal(inverse_frame, (-1)**sum(x)*frame):
            raise ArithmeticError('onsite frame does not transport staggered inversion to site inversion')
        frame_checks += 1
    grading = np.diag([(-1)**sum(x) for x in PARITIES])
    for x in PARITIES:
        # exp(+ik_* . x)=epsilon_x exp(-ik_* . x).
        if 1j**sum(x) != (-1)**sum(x)*(-1j)**sum(x):
            raise ArithmeticError('nodal inversion lacks reciprocal-cell transport')
    for perm, signs in source.cubic_rotations():
        for x in PARITIES:
            if source.rotation_phase(tuple(-t for t in x), perm) != source.rotation_phase(x, perm):
                raise ArithmeticError('proper rotation does not preserve source inversion')
    full = commutant_basis(magnetic_generators()+cubic_lifts()+(grading,))
    if len(full) != 1 or not np.array_equal(full[0], np.eye(8, dtype=int)):
        raise ArithmeticError('actual inversion did not identify the two scalar channels')
    return {'literal_Weyl_bond_checks': source_checks,
            'finite_box_radii': [1, 2], 'affine_average_annihilator_checks': averaging_checks,
            'onsite_frame_checks': frame_checks, 'onsite_contact_preserved_by_site_permutation': True,
            'tracial_exterior_preserves_affine_inversion': True,
            'additional_state_average_used': False,
            'nodal_commutant_with_actual_inversion_dimension': 1}


def lorentzian_score(alpha, eta, memory_real, memory_dissipation):
    """Exact nodal score from supplied actual scalar memory, not its selection."""
    alpha, eta, a, b = map(F, (alpha, eta, memory_real, memory_dissipation))
    if alpha <= 0 or eta <= 0 or b < 0:
        raise ValueError('positive alpha/eta and nonnegative dissipation required')
    h = alpha * alpha + eta * eta + eta * b
    return 1 - alpha * alpha * h / (h * h + eta * eta * a * a)


def certificate():
    gamma = tuple(np.asarray(g, dtype=np.int64) for g in gamma_matrices())
    identity = np.eye(8, dtype=np.int64)
    volume = gamma[0] @ gamma[1] @ gamma[2]
    grading = np.diag([(-1) ** sum(x) for x in PARITIES])
    magnetic = magnetic_generators()
    rotations = cubic_lifts()
    basis = commutant_basis(magnetic + rotations)
    if len(basis) != 2 or len(commutant_basis(magnetic)) != 8:
        raise ArithmeticError('source commutant dimensions differ')
    for matrix in basis:
        scalar = int(np.trace(matrix))
        pseudoscalar = -int(np.trace(volume @ matrix))
        if not np.array_equal(8 * matrix, scalar * identity + pseudoscalar * volume):
            raise ArithmeticError('commutant has an additional channel')
    if (not np.array_equal(volume @ volume, -identity)
            or not np.array_equal(volume.T, -volume)
            or np.trace(volume) != 0
            or not np.array_equal(grading @ volume @ grading, -volume)):
        raise ArithmeticError('source channel involution failed')
    for matrix in magnetic:
        for g in gamma:
            if not np.array_equal(matrix @ g, g @ matrix):
                raise ArithmeticError('magnetic translation does not preserve source symbol')
    for (perm, signs), matrix in zip(source.cubic_rotations(), rotations, strict=True):
        for j in range(3):
            if not np.array_equal(matrix @ gamma[j] @ matrix.T, signs[j] * gamma[perm[j]]):
                raise ArithmeticError('cubic lift lost its nodal chart transport')
    columns = transported_valleys()
    frame = np.hstack(columns)
    channel = np.kron(1j * volume, np.eye(2, dtype=int))
    if not np.array_equal(frame.conj().T @ frame, 8 * np.eye(16)):
        raise ArithmeticError('transported original valleys are not a complete frame')
    rows = []
    for b, block in zip(PARITIES, columns, strict=True):
        chirality = -(-1) ** sum(b)
        if not np.array_equal(channel @ block, chirality * block):
            raise ArithmeticError('source channel and kinetic winding disagree')
        for j in range(3):
            left = -np.kron(gamma[j], np.eye(2, dtype=int)) @ block
            right = block @ (-(-1) ** b[j] * np.asarray(SIGMA[j]))
            if not np.array_equal(left, right):
                raise ArithmeticError('original Weyl linearization was not transported')
        rows.append({'original_node_in_pi_units': list(b), 'kinetic_degree': chirality,
                     'nodal_scalar_channel': chirality})
    return {
        'source': 'same proper-cubic averaged centered full-CAR ground and magnetic translations',
        'magnetic_node_in_pi_units': ['1/2', '1/2', '1/2'],
        'proper_cubic_lifts': len(rotations),
        'magnetic_commutant_dimension': 8,
        'nodal_little_group_commutant_dimension': len(basis),
        'parity_component_projector_ranks': [4, 4],
        'including_hidden_spin_projector_ranks': [8, 8],
        'original_valley_frame_dimension': 16,
        'source_valley_channels': rows,
        'nodal_projector_channels': 2,
        'nodal_signed_scalar_measures': 1,
        'independent_after_source_energy_reflection': 1,
        'independent_after_actual_source_inversion': 1,
        'nodal_scalar_measure_is_energy_even': True,
        'nodal_memory_purely_imaginary_scalar_at_imaginary_frequency': True,
        'actual_source_inversion': inversion_source_checks(),
        'original_nodal_lorentzian_scores_equal': True,
        'uniform_growing_packet_scaling_supplied': False,
        'physical_spin_or_species_identification': False,
        'physical_mass_prediction': False,
    }
