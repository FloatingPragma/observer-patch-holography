"""Full-CAR static reference, physical preparation and dressed arrivals."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI, source_descriptor
from ..model import sequence, clean
from ..periodic_family.model import (Configuration, apply_state, norm2, kinetic_terms,
    edges, cell_index, shift_cell)
from ..peter_weyl_links import model as source
from ..peter_weyl_links.clebsch import irreps, native_charges


def validate_key(key, config=Configuration()):
    state, reference = key
    source.validate_key(state, config)
    if type(reference) is not int or not 0 <= reference < (1 << 30):
        raise ValueError('full static reference F30 required')
    return state, reference


def empty_key(config=Configuration()):
    return source.empty_key(config), 0


def column(key, config=Configuration()):
    state, reference = validate_key(key, config)
    return {(target, reference): value for target, value in source.column(state, config).items()}


@lru_cache(maxsize=64)
def one_body_columns(matrix):
    out = defaultdict(list)
    for (a, b), value in matrix.todok().items():
        out[b].append((a, value))
    return dict(out)


def one_body_column(mask, matrix, offset=0):
    out = defaultdict(lambda: s.Integer(0))
    if not isinstance(matrix, s.ImmutableSparseMatrix):
        matrix = s.ImmutableSparseMatrix(matrix)
    columns = one_body_columns(matrix)
    occupied = (mask >> offset) & ((1 << matrix.cols) - 1)
    while occupied:
        bit = occupied & -occupied; b = bit.bit_length() - 1
        occupied ^= bit
        for a, value in columns.get(b, ()):
            step = sequence(mask, ((offset + a, True), (offset + b, False)))
            if step is not None:
                out[step[0]] += value * step[1]
    return clean(out)


@lru_cache(maxsize=12)
def reference_current(current):
    return s.ImmutableSparseMatrix(-s.kronecker_product(s.eye(2), representation()[2][current].T))


def gauss_column(key, site, current, config=Configuration()):
    state, reference = validate_key(key, config)
    out = defaultdict(lambda: s.Integer(0), {
        (target, reference): value for target, value in source.gauss_column(state, site, current, config).items()})
    if site == 0:
        for target, value in one_body_column(reference, reference_current(current)).items():
            out[state, target] += value
    return clean(out)


def spin_column(key, axis, config=Configuration()):
    state, reference = validate_key(key, config)
    mask, links, scalars = state
    matrix = s.kronecker_product(PAULI[axis] / 2, s.eye(15))
    out = defaultdict(lambda: s.Integer(0))
    for site in range(config.side ** 3):
        for target, value in one_body_column(mask, matrix, 30 * site).items():
            out[(target, links, scalars), reference] += value
    for target, value in one_body_column(reference, matrix).items():
        out[state, target] += value
    return clean(out)


def preparation(charge, config=Configuration()):
    if charge not in native_charges():
        raise ValueError('one of the five actual native blocks required')
    block = irreps()[charge]['native_indices']; d = len(block)
    _, links, scalars = source.empty_key(config)
    return {((1 << (15 * spin + a), links, scalars), 1 << (15 * (1 - spin) + a)):
            s.Integer((-1) ** spin) / s.sqrt(2 * d)
            for spin in range(2) for a in block}


def preparation_pulse_column(key, charge, config=Configuration()):
    """K=|Phi><vac|+|vac><Phi| on origin CAR and reference, identity elsewhere.

    Both active vectors are even spin/Gauss singlets. The spectator fields
    are unrestricted, so this defines the finite unitary pulse on all states.
    """
    if charge not in native_charges():
        raise ValueError('one of the five actual native blocks required')
    state, reference = validate_key(key, config)
    mask, links, scalars = state
    local = mask & ((1 << 30) - 1)
    spectator = mask ^ local
    # Reference modes follow every source mode in the matrix ordering.
    # The even local pair crosses the intervening source parity exactly once.
    phase = s.Integer((-1) ** spectator.bit_count())
    block = irreps()[charge]['native_indices']; d = len(block)
    if local == 0 and reference == 0:
        return {((spectator | (1 << (15 * spin + a)), links, scalars),
                 1 << (15 * (1 - spin) + a)): phase * (-1) ** spin / s.sqrt(2 * d)
                for spin in range(2) for a in block}
    for spin in range(2):
        for a in block:
            if local == 1 << (15 * spin + a) and reference == 1 << (15 * (1 - spin) + a):
                return {((spectator, links, scalars), 0): phase * (-1) ** spin / s.sqrt(2 * d)}
    return {}


@lru_cache(maxsize=16)
def directed_arrivals(charge, config=Configuration()):
    """Two orthonormal compensated strings from the actual joining summands."""
    if config.side != 2:
        raise ValueError('the declared two-edge nearest-site record uses side two')
    psi = preparation(charge, config)
    result = []
    for term in kinetic_terms(config):
        if set((term.left_cell, term.right_cell)) != {0, 4}:
            continue
        image = apply_state(psi, lambda key: {
            (target, key[1]): value for target, value in
            source.hopping_column(key[0], (term,)).items()})
        assert norm2(image) == abs(term.coefficient) ** 2
        normalized = {key: value / abs(term.coefficient) for key, value in image.items()}
        result.append((term.path, normalized))
    assert len(result) == 2
    assert not (result[0][1].keys() & result[1][1].keys())
    return tuple(result)


def inner(left, right):
    return s.expand(sum(s.conjugate(value) * right.get(key, 0) for key, value in left.items()))


def arrival_projector(state, charge, config=Configuration()):
    out = defaultdict(lambda: s.Integer(0))
    for _, vector in directed_arrivals(charge, config):
        amplitude = inner(vector, state)
        for key, value in vector.items():
            out[key] += value * amplitude
    return clean(out)


def complete_readout(state, charge, config=Configuration()):
    arrived = arrival_projector(state, charge, config)
    other = defaultdict(lambda: s.Integer(0), state)
    for key, value in arrived.items():
        other[key] -= value
    return arrived, clean(other)


def local_arrival_projector(state, charge, config=Configuration()):
    """Gauge-invariant charged matter and one joining flux, all spectators kept."""
    if config.side != 2:
        raise ValueError('the declared local detector uses side two')
    target = 4; block = irreps()[charge]['native_indices']
    incident = [j for j, (xyz, axis) in enumerate(edges(2))
                if target in (cell_index(2, xyz), cell_index(2, shift_cell(2, xyz, axis)))]
    out = {}
    for key, value in state.items():
        (mask, links, scalars), _ = validate_key(key, config)
        local = (mask >> (30 * target)) & ((1 << 30) - 1)
        if local.bit_count() != 1 or (local.bit_length() - 1) % 15 not in block or scalars[target] != 0:
            continue
        excited = [edge for edge in incident if links[edge] != 0]
        if len(excited) != 1 or excited[0] not in (0, 12):
            continue
        edge = excited[0]; sign, a, b = source.labels()[links[edge]]
        if sign == (-1 if edge == 0 else 1) and a in block and b in block:
            out[key] = value
    return out


def descriptor(config=Configuration()):
    return {'schema': 'oph.common_source_native_family_transport.v1',
        'native_site_program': source_descriptor(), 'source': source.descriptor(config),
        'reference': 'static full30-CAR dual-R15 memory at origin, with two ordinary spin components',
        'reference_full_Fock_dimension': str(2 ** 30),
        'reference_gauge_current': 'dGamma(-I_spin tensor T^transpose)',
        'joint_clock': 'exp(-it H_PW) tensor identity_reference',
        'origin_Gauss_constraint': 'G_origin(source) + Q_reference',
        'other_Gauss_constraints': 'G_site(source)',
        'additional_source_Hamiltonian_term': False,
        'same_full_source_and_reference_for_all_five_preparations': True,
        'charged_or_colored_arrival_without_compensating_reference': False,
        'ground_or_free_particle_dispersion_inferred': False,
        'observed_particle_or_mass_dictionary_asserted': False}
