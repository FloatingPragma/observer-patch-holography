"""Full-Fock operator crosswalk for the configured native-current CAR field."""
from functools import lru_cache
import hashlib

import sympy as s

from model.common_source import native_matter_attachment as matter
from . import algebra


def mode_index(spin, internal, family):
    if any(type(x) is not int for x in (spin, internal, family)):
        raise TypeError('exact integer field indices required')
    if not (0 <= spin < 2 and 0 <= internal < 15 and 0 <= family < 3):
        raise ValueError('the declared two-spin, fifteen-state, three-family ordering is required')
    return 45*spin+3*internal+family


def _permutation(values):
    values = tuple(values)
    if not values or any(type(x) is not int for x in values) or sorted(values) != list(range(len(values))):
        raise ValueError('a complete finite permutation is required')
    return values


def _mask(mask, modes):
    if type(modes) is not int or modes < 1 or type(mask) is not int or not 0 <= mask < (1 << modes):
        raise ValueError('an occupation mask on the stated complete Fock carrier is required')


def permutation_from_matrix(matrix):
    matrix = s.ImmutableSparseMatrix(matrix)
    if matrix.rows != matrix.cols or matrix.rows == 0:
        raise ValueError('a nonempty square permutation matrix is required')
    values = []
    for column in range(matrix.cols):
        occupied = [row for row in range(matrix.rows) if matrix[row, column] != 0]
        if len(occupied) != 1 or matrix[occupied[0], column] != 1:
            raise ValueError('the exact source map must be a permutation with unit entries')
        values.append(occupied[0])
    return _permutation(values)


def _create(mode, mask):
    if mask >> mode & 1:
        return None
    return mask | (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def _annihilate(mode, mask):
    if not mask >> mode & 1:
        return None
    return mask ^ (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def creation_action(mode, mask, modes):
    _mask(mask, modes)
    if type(mode) is not int or not 0 <= mode < modes:
        raise ValueError('a creation mode in the carrier is required')
    return _create(mode, mask)


def annihilation_action(mode, mask, modes):
    _mask(mask, modes)
    if type(mode) is not int or not 0 <= mode < modes:
        raise ValueError('an annihilation mode in the carrier is required')
    return _annihilate(mode, mask)


def _exterior(mask, permutation):
    image, sign = 0, 1
    while mask:
        bit = mask & -mask
        destination = permutation[bit.bit_length()-1]
        sign *= (-1)**((image >> (destination+1)).bit_count())
        image |= 1 << destination
        mask ^= bit
    return image, sign


def exterior_permutation(mask, permutation):
    """Gamma(P) on a complete occupation column, including its wedge sign."""
    permutation = _permutation(permutation)
    _mask(mask, len(permutation))
    return _exterior(mask, permutation)


@lru_cache(maxsize=1)
def source_currents():
    """Twelve immutable anti-Hermitian M15 matrices in literal native-port order."""
    return tuple(s.ImmutableSparseMatrix(g) for g in matter.representation()['matter_currents'])


def physical_permutation():
    return permutation_from_matrix(matter.field_intertwiner())


def _add(out, mask, coefficient):
    value = s.expand(out.get(mask, 0)+coefficient)
    if value:
        out[mask] = value
    else:
        out.pop(mask, None)


def _one_body(entries, mask):
    out = {}
    for (i, j), value in entries.items():
        removed = _annihilate(j, mask)
        if removed is None:
            continue
        created = _create(i, removed[0])
        if created is not None:
            _add(out, created[0], value*removed[1]*created[1])
    return out


def one_body_action(generator, mask):
    """Literal sum G_ij a_i^dagger a_j on any occupation sector."""
    generator = s.ImmutableSparseMatrix(generator)
    if (generator.rows != generator.cols or generator.rows == 0
            or any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in generator)):
        raise ValueError('an exact finite coefficient matrix is required')
    _mask(mask, generator.rows)
    return _one_body(generator.todok(), mask)


def current_creation_residual(generator, mode, mask):
    """[dGamma(G),a_mode^dagger]-sum_i G_i,mode a_i^dagger, as a full-Fock column."""
    generator = s.ImmutableSparseMatrix(generator)
    # Admission is shared with the independently usable literal action.
    original = one_body_action(generator, mask)
    created = creation_action(mode, mask, generator.rows)
    out = {}
    if created is not None:
        for target, value in _one_body(generator.todok(), created[0]).items():
            _add(out, target, created[1]*value)
    for target, value in original.items():
        final = _create(mode, target)
        if final is not None:
            _add(out, final[0], -value*final[1])
    for i in range(generator.rows):
        final = _create(i, mask)
        if final is not None:
            _add(out, final[0], -generator[i, mode]*final[1])
    return out


def wedge_census(permutation):
    """Exhaust the signed lift and every creation column for at most sixteen modes."""
    permutation = _permutation(permutation)
    n = len(permutation)
    if n > 16:
        raise ValueError('exhaustive finite occupation control is capped at sixteen modes')
    dimension = 1 << n
    images, signs = [0]*dimension, [1]*dimension
    for mask in range(1, dimension):
        bit = mask & -mask
        rest = mask ^ bit
        destination = permutation[bit.bit_length()-1]
        images[mask] = images[rest] | (1 << destination)
        signs[mask] = signs[rest]*(-1)**((images[rest] & ((1 << destination)-1)).bit_count())
    if len(set(images)) != dimension or any(images[m].bit_count() != m.bit_count() for m in range(dimension)):
        raise ArithmeticError('the exterior source map is not grade-preserving and unitary')
    for mask in range(dimension):
        for mode in range(n):
            left = _create(mode, mask)
            right = _create(permutation[mode], images[mask])
            if left is None:
                if right is not None:
                    raise ArithmeticError('source Pauli blocking changed under the lift')
            elif (images[left[0]], left[1]*signs[left[0]]) != (right[0], signs[mask]*right[1]):
                raise ArithmeticError('the signed creation intertwiner failed')
    width = (n+7)//8
    digest = hashlib.sha256(b''.join(image.to_bytes(width, 'little')+bytes([sign == 1])
                                    for image, sign in zip(images, signs))).hexdigest()
    return {'modes': n, 'full_Fock_dimension': dimension,
            'creation_columns_including_zero_columns': n*dimension,
            'all_occupation_sectors_and_creation_columns_checked': True,
            'signed_permutation_columns_sha256': digest}


def casimir_and_weight_records():
    """Casimirs act on the creation-operator multiplet by double commutators."""
    r = matter.representation(); selected = list(r['selected'])
    canonical = matter.canonical_currents()
    hermitian = [matter.exterior_derivative(s.I*matter.spin_lift(g)).extract(selected, selected)
                 for g in canonical]
    color = sum((g*g/4 for g in hermitian[:6]), s.zeros(15))
    h1, h2 = hermitian[6:8]
    color += (h1*h1+h2*h2+h1*h2)/3
    weak = sum((g*g for g in hermitian[8:11]), s.zeros(15))
    for name, actual in (('color_C2', color), ('weak_C2', weak)):
        if not matter.zero(actual-r['operators'][name].extract(selected, selected)):
            raise ArithmeticError('the field-operator Casimir differs from the internal source record')
    matches = algebra.match_components()
    targets = {mask: name for name, record in matches.items() for mask in record['source_masks']}
    permutation = physical_permutation()
    records = []
    for internal, record in enumerate(algebra.source_components()):
        native_modes = [mode_index(spin, internal, family) for spin in range(2) for family in range(3)]
        records.append({**record, 'color_weight': list(record['color_weight']),
                        'native_internal_index': internal,
                        'native_field_modes': native_modes,
                        'action_field_modes': [permutation[index] for index in native_modes],
                        'internal_comparison_component': targets[record['occupation_mask']],
                        'auxiliary_occupation_parity': 1,
                        'field_creation_parity': -1,
                        'physical_species_assignment': None})
    return records


def certificate():
    """Replay source intertwiners and exact finite controls; select no physical species."""
    r = matter.representation()
    p15, _ = matter.action_intertwiner()
    p90 = matter.field_intertwiner()
    permutation = physical_permutation()
    old = matter.existing_action._representation()
    currents = source_currents()
    for port, current in enumerate(currents):
        target = old[2][r['native']['native_to_fixture_port'][port]]
        native90 = s.kronecker_product(s.eye(2), current, s.eye(3))
        target90 = s.kronecker_product(s.eye(2), target, s.eye(3))
        if not matter.zero(target90*p90-p90*native90):
            raise ArithmeticError('the ninety-mode source-current operator frame differs')
    masks = (0, (1 << 15)-1, int('101010101010101', 2), int('000111100110011', 2))
    for current in currents:
        for mask in masks:
            for mode in range(15):
                if current_creation_residual(current, mode, mask):
                    raise ArithmeticError('the actual-current full-Fock commutator failed')
    for mask in (0, (1 << 90)-1, sum(1 << i for i in range(0, 90, 2))):
        image, sign = exterior_permutation(mask, permutation)
        for mode in range(90):
            left = _create(mode, mask); right = _create(permutation[mode], image)
            if left is None:
                if right is not None:
                    raise ArithmeticError('ninety-mode source blocking differs')
            else:
                final, factor = _exterior(left[0], permutation)
                if (final, left[1]*factor) != (right[0], sign*right[1]):
                    raise ArithmeticError('ninety-mode source creation sign differs')
    return {
        'architecture': 'existing configured C2_spin tensor M15 tensor C3_family CAR field',
        'physical_CAR_modes': 90, 'full_Fock_dimension': str(1 << 90),
        'mode_order': 'spin*45+3*native_internal_index+family',
        'all_twelve_native_port_currents_intertwined': True,
        'current_convention': 'G anti-Hermitian; Hermitian charge T=iG; dGamma(T)=sum T_ij a_i^dagger a_j',
        'full_Fock_creation_identity': '[dGamma(T),a_j^dagger]=sum_i T_ij a_i^dagger',
        'CAR_and_creation_identity_scope': 'all occupation sectors by the CAR algebra, without a one-particle truncation',
        'complete_internal_wedge_control': wedge_census(permutation_from_matrix(p15)),
        'actual_current_occupation_column_controls': len(currents)*len(masks)*15,
        'ninety_mode_creation_column_controls': 270,
        'native_to_action_one_particle_permutation': list(permutation),
        'field_operator_charge_records': casimir_and_weight_records(),
        'Casimir_labels_use_creation_operator_double_commutators': True,
        'one_body_Casimir_identified_with_total_many_body_Casimir': False,
        'annihilation_transforms_in_dual_representation': True,
        'annihilation_identified_with_a_created_antiparticle_without_reference': False,
        'auxiliary_even_register_identified_with_odd_field_Fock_space': False,
        'configured_CAR_architecture_selected_by_unread_repair': False,
        'three_family_factor_is_declared': True,
        'mass_or_family_eigenbasis_selected': False,
        'default_two_flavour_interacting_ground_used': False,
        'additional_native_pulse_tape_executed': False,
        'interacting_relativistic_state_or_particle_pole_proved': False,
        'physical_species': [], 'physical_masses': [], 'empirical_mass_input_used': False,
    }
