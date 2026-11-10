"""A declared 119-state link architecture, with exact native endpoint currents."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI, source_descriptor
from ...chiral_higgs.model import higgs_memory
from ..model import clean, sequence
from ..periodic_family.model import (
    Configuration, cells, cell_index, shift_cell, edges, kinetic_terms,
    native_blocks, local_yukawa, translate_fermions,
)


@lru_cache(maxsize=1)
def labels():
    """The sign labels flux orientation, and does not duplicate matter species."""
    return ((0, -1, -1),) + tuple(
        (sign, a, b) for sign in (1, -1)
        for block in native_blocks() for a in block for b in block)


@lru_cache(maxsize=1)
def label_index():
    return {label: index for index, label in enumerate(labels())}


@lru_cache(maxsize=1)
def block_by_index():
    return {a: block for block in native_blocks() for a in block}


@lru_cache(maxsize=256)
def link(a, b):
    """U_ab=(|+,a,b><0|+|0><-,a,b|)/sqrt(dim(block))."""
    if a not in block_by_index() or b not in block_by_index()[a]:
        raise ValueError('indices must lie in one native irreducible block')
    index = label_index()
    factor = 1 / s.sqrt(len(block_by_index()[a]))
    return s.ImmutableSparseMatrix(119, 119, {
        (index[1, a, b], 0): factor,
        (0, index[-1, a, b]): factor,
    })


@lru_cache(maxsize=24)
def endpoint_current(current, endpoint):
    """Forward: (-T^transpose,T); reverse: (T,-T^transpose)."""
    if current not in range(12) or endpoint not in (0, 1):
        raise ValueError('one actual native current and one endpoint required')
    t = representation()[2][current]
    index = label_index()
    result = {}
    for col, (sign, a, b) in enumerate(labels()[1:], 1):
        for (i, j), value in t.todok().items():
            if endpoint == 0 and sign == 1 and i == a:
                target, coefficient = (sign, j, b), -value
            elif endpoint == 1 and sign == 1 and j == b:
                target, coefficient = (sign, a, i), value
            elif endpoint == 0 and sign == -1 and j == a:
                target, coefficient = (sign, i, b), value
            elif endpoint == 1 and sign == -1 and i == b:
                target, coefficient = (sign, a, j), -value
            else:
                continue
            result[index[target], col] = coefficient
    return s.ImmutableSparseMatrix(119, 119, result)


@lru_cache(maxsize=1)
def orientation():
    index = label_index()
    entries = {(0, 0): 1}
    for col, (sign, a, b) in enumerate(labels()[1:], 1):
        entries[index[-sign, b, a], col] = 1
    return s.ImmutableSparseMatrix(119, 119, entries)


def reverse_label(label):
    sign, a, b = labels()[label]
    return 0 if sign == 0 else label_index()[-sign, b, a]


@lru_cache(maxsize=5)
def transporter(block):
    """The full internal matrix, including every state of the link carrier."""
    d = len(block)
    entries = {}
    for a, ia in enumerate(block):
        for b, ib in enumerate(block):
            for (row, col), value in link(ia, ib).todok().items():
                entries[119 * a + row, 119 * b + col] = value
    return s.ImmutableSparseMatrix(119 * d, 119 * d, entries)


def validate_key(key, config=Configuration()):
    mask, links, scalars = key
    volume = config.side ** 3
    if type(mask) is not int or not 0 <= mask < 1 << (30 * volume):
        raise ValueError('occupation outside full Fock carrier')
    if (not isinstance(links, tuple) or len(links) != 3 * volume
            or any(type(v) is not int or not 0 <= v < 119 for v in links)):
        raise ValueError('complete 119-state link carrier required')
    if (not isinstance(scalars, tuple) or len(scalars) != volume
            or any(type(v) is not int or not 0 <= v < 5 for v in scalars)):
        raise ValueError('complete five-state scalar carrier required')
    return mask, links, scalars


def empty_key(config=Configuration()):
    return 0, (0,) * (3 * config.side ** 3), (0,) * (config.side ** 3)


@lru_cache(maxsize=32768)
def path_transitions(values, adjoint=False):
    """All W_ab (or W_ab^dagger) columns, on arbitrary link occupations.

    The adjoint output keeps the original a,b labels. The matter word then
    reverses its endpoints. Every axial path in this family has distinct edges.
    """
    if not values or any(type(v) is not int or v not in range(119) for v in values):
        raise ValueError('a nonempty word of valid link states is required')
    transitions = []
    for block in native_blocks():
        d = len(block)
        candidates = []
        for value in values:
            sign, a, b = labels()[value]
            if sign == 0:
                candidates.append(tuple((i, j, label_index()[(-1 if adjoint else 1), i, j])
                                        for i in block for j in block))
            elif sign == (1 if adjoint else -1) and a in block:
                candidates.append(((a, b, 0),))
            else:
                candidates = []
                break
        if not candidates:
            continue
        chains = [(a, b, (target,)) for a, b, target in candidates[0]]
        for entries in candidates[1:]:
            chains = [(a, last, targets + (target,))
                      for a, middle, targets in chains
                      for first, last, target in entries if middle == first]
        factor = d ** (-s.Rational(len(values), 2))
        transitions.extend((a, b, targets, factor) for a, b, targets in chains)
    return tuple(transitions)


def column(key, config=Configuration()):
    """Exact complete-carrier Hamiltonian action for the declared architecture."""
    mask, links, scalars = validate_key(key, config)
    out = defaultdict(lambda: s.Integer(0))
    out[key] += (config.link_penalty * sum(v != 0 for v in links)
                 + config.scalar_penalty * sum(v != 0 for v in scalars))
    for cell, scalar in enumerate(scalars):
        offset = 30 * cell
        local_mask = (mask >> offset) & ((1 << 30) - 1)
        couplings = config.y_up, config.y_down, config.y_lepton
        for (target, h), coefficient in local_yukawa(local_mask, scalar, couplings).items():
            target_scalars = list(scalars)
            target_scalars[cell] = h
            target_mask = mask ^ ((local_mask ^ target) << offset)
            out[target_mask, links, tuple(target_scalars)] += coefficient
    for term in kinetic_terms(config):
        values = tuple(links[e] for e in term.path)
        for adjoint in (False, True):
            paths = path_transitions(values, adjoint)
            for (u, v), spin_weight in PAULI[term.axis].todok().items():
                weight = term.coefficient * spin_weight
                if adjoint:
                    weight = s.conjugate(weight)
                for a, b, targets, factor in paths:
                    p = 30 * term.left_cell + 15 * u + a
                    q = 30 * term.right_cell + 15 * v + b
                    if adjoint:
                        p, q = q, p
                    step = sequence(mask, ((p, True), (q, False)))
                    if step is None:
                        continue
                    target_links = list(links)
                    for edge, target in zip(term.path, targets, strict=True):
                        target_links[edge] = target
                    out[step[0], tuple(target_links), scalars] += weight * factor * step[1]
    return clean(out)


@lru_cache(maxsize=36)
def current_columns(current, endpoint=None):
    if endpoint is None:
        matrix = representation()[2][current]
    elif endpoint == 'scalar':
        matrix = higgs_memory()[1][current]
    else:
        matrix = endpoint_current(current, endpoint)
    result = defaultdict(list)
    for (row, col), value in s.ImmutableSparseMatrix(matrix).todok().items():
        result[col].append((row, value))
    return dict(result)


def gauss_column(key, cell, current, config=Configuration()):
    mask, links, scalars = validate_key(key, config)
    n = config.side
    if cell not in range(n ** 3) or current not in range(12):
        raise ValueError('cell/current outside source')
    out = defaultdict(lambda: s.Integer(0))
    for spin in range(2):
        for b, changes in current_columns(current).items():
            origin = 30 * cell + 15 * spin + b
            if not mask >> origin & 1:
                continue
            for a, coefficient in changes:
                step = sequence(mask, ((30 * cell + 15 * spin + a, True), (origin, False)))
                if step is not None:
                    out[step[0], links, scalars] += coefficient * step[1]
    for target, coefficient in current_columns(current, 'scalar').get(scalars[cell], ()):
        target_scalars = list(scalars)
        target_scalars[cell] = target
        out[mask, links, tuple(target_scalars)] += coefficient
    xyz = cells(n)[cell]
    for axis in range(3):
        incident = ((3 * cell + axis, 0),
                    (3 * cell_index(n, shift_cell(n, xyz, axis, -1)) + axis, 1))
        for edge, endpoint in incident:
            for target, coefficient in current_columns(current, endpoint).get(links[edge], ()):
                target_links = list(links)
                target_links[edge] = target
                out[mask, tuple(target_links), scalars] += coefficient
    return clean(out)


def translation_column(key, axis, config=Configuration()):
    if axis not in range(3):
        raise ValueError('three spatial translations')
    mask, links, scalars = validate_key(key, config)
    mask, phase = translate_fermions(mask, axis, config.side)
    target_links, target_scalars = [0] * len(links), [0] * len(scalars)
    for cell, xyz in enumerate(cells(config.side)):
        target = cell_index(config.side, shift_cell(config.side, xyz, axis))
        target_scalars[target] = scalars[cell]
        for direction in range(3):
            target_links[3 * target + direction] = links[3 * cell + direction]
    return {(mask, tuple(target_links), tuple(target_scalars)): phase}


def descriptor(config=Configuration()):
    volume = config.side ** 3
    return {
        'schema': 'oph.chiral_oriented_quantum_links.v1',
        'native_site_program': source_descriptor()['program'],
        'full_dimension': f'2^{30 * volume} * 119^{3 * volume} * 5^{volume}',
        'native_block_dimensions': [len(block) for block in native_blocks()],
        'matter_charge_conjugation_under_spatial_rotation': False,
        'link_dimension': 119,
        'side': config.side,
        'spacing': str(config.spacing),
        'source_coefficients': {name: str(getattr(config, name)) for name in
            ('y_up', 'y_down', 'y_lepton', 'scalar_penalty', 'link_penalty')},
        'kinetic_coefficient_source': 'periodic_family.kinetic_terms, with U replaced by the declared 119-state operator',
        'link_amplitude_source': 'architecture choice 1/sqrt(native block dimension)',
        'clock': 'exp(-it H_oriented) on this complete carrier',
        'literal_Peter_Weyl_compression': False,
        'native_idle_link_implementation_proved': False,
        'raw_link_spectrum_or_ground_transferred': False,
        'measured_particle_or_GeV_mass_identified': False,
        'interacting_relativistic_refinement_proved': False,
    }
