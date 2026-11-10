"""Actual periodic three-axis source; all Fock, scalar and quantum links retained."""
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import product

import sympy as s

from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI, link_labels, link_generator, source_descriptor
from ...chiral_regulator.model import momentum_matrix
from ...chiral_higgs.model import higgs_memory, tensor_checks
from ..model import sequence, clean


@dataclass(frozen=True)
class Configuration:
    side: int = 2
    spacing: object = s.pi
    y_up: object = s.Integer(1)
    y_down: object = s.Integer(1)
    y_lepton: object = s.Integer(1)
    scalar_penalty: object = s.Integer(1)
    link_penalty: object = s.Integer(1)

    def __post_init__(self):
        if type(self.side) is not int or self.side < 2 or self.side % 2:
            raise ValueError('an even side at least two is required')
        for name in ('spacing', 'y_up', 'y_down', 'y_lepton', 'scalar_penalty', 'link_penalty'):
            value = s.sympify(getattr(self, name))
            if value.has(s.Float) or value.is_real is not True or value.is_finite is not True:
                raise ValueError('finite exact real source coefficients required')
            object.__setattr__(self, name, value)
        if self.spacing.is_positive is not True:
            raise ValueError('positive spacing required')
        if self.scalar_penalty.is_nonnegative is not True or self.link_penalty.is_nonnegative is not True:
            raise ValueError('nonnegative scalar/link penalties required')

    @property
    def length(self):
        return self.side * self.spacing


@lru_cache(maxsize=None)
def cells(side):
    return tuple(product(range(side), repeat=3))


def cell_index(side, cell):
    x, y, z = cell
    if any(type(q) is not int or not 0 <= q < side for q in cell):
        raise ValueError('cell outside periodic cube')
    return (x * side + y) * side + z


def shift_cell(side, cell, axis, steps=1):
    target = list(cell)
    target[axis] = (target[axis] + steps) % side
    return tuple(target)


@lru_cache(maxsize=None)
def edges(side):
    return tuple((cell, axis) for cell in cells(side) for axis in range(3))


def edge_endpoints(side, edge):
    cell, axis = edges(side)[edge]
    return cell, shift_cell(side, cell, axis)


def path(side, anchor, axis, displacement):
    if axis not in range(3) or not 1 <= displacement < side:
        raise ValueError('forward axial displacement must lie in 1..side-1')
    return tuple(3 * cell_index(side, shift_cell(side, anchor, axis, step)) + axis
                 for step in range(displacement))


@dataclass(frozen=True)
class KineticTerm:
    left_cell: int
    right_cell: int
    axis: int
    displacement: int
    path: tuple
    coefficient: object


@lru_cache(maxsize=None)
def kinetic_terms(config=Configuration()):
    """coefficient*sigma_axis multiplies c_left^dagger W_path c_right+h.c."""
    n = config.side
    derivative = momentum_matrix(n, config.length)
    result = []
    for left, anchor in enumerate(cells(n)):
        for axis in range(3):
            for displacement in range(1, n):
                endpoint = shift_cell(n, anchor, axis, displacement)
                coefficient = s.simplify(-derivative[anchor[axis], endpoint[axis]] / 2)
                wrap = (anchor[axis] + displacement) // n
                formula = (-1) ** wrap * (-s.I * s.pi / (2 * config.length)) * (-1) ** displacement / s.sin(s.pi * displacement / n)
                assert s.simplify(coefficient - formula) == 0
                result.append(KineticTerm(left, cell_index(n, endpoint), axis, displacement,
                    path(n, anchor, axis, displacement), coefficient))
    return tuple(result)


@lru_cache(maxsize=1)
def native_blocks():
    groups = {}
    for i, y in enumerate(representation()[2][-1].diagonal()):
        groups.setdefault(y, []).append(i)
    return tuple(tuple(values) for values in groups.values())


@lru_cache(maxsize=1)
def _native_block_by_index():
    return {i: block for block in native_blocks() for i in block}


@lru_cache(maxsize=1)
def _link_index():
    return {label: i for i, label in enumerate(link_labels())}


@lru_cache(maxsize=1)
def yukawa_words():
    fields, _, _, _ = higgs_memory()
    result = []
    for channel, vertex in representation()[3].items():
        for k, tensor in enumerate(vertex['components']):
            field = fields[k].H if vertex['conjugate_Higgs'] else fields[k]
            for (i, j), coefficient in s.ImmutableSparseMatrix(tensor).todok().items():
                for scalar, ops, weight in ((field, ((i, False), (15 + j, False)), coefficient),
                    (field.H, ((15 + j, True), (i, True)), s.conjugate(coefficient))):
                    for (row, col), value in scalar.todok().items():
                        result.append((channel, row, col, ops, s.expand(weight * value)))
    return tuple(result)


@lru_cache(maxsize=65536)
def local_yukawa(mask, scalar, couplings=(s.Integer(1), s.Integer(1), s.Integer(1))):
    weights = dict(zip(('up', 'down', 'charged_lepton'), couplings, strict=True))
    result = defaultdict(lambda: s.Integer(0))
    for channel, row, col, ops, coefficient in yukawa_words():
        if col != scalar or weights[channel] == 0:
            continue
        transition = sequence(mask, ops)
        if transition is not None:
            result[transition[0], row] += weights[channel] * coefficient * transition[1]
    return clean(result)


def empty_key(config=Configuration()):
    return 0, (0,) * (3 * config.side ** 3), (0,) * (config.side ** 3)


def validate_key(key, config=Configuration()):
    mask, links, scalars = key
    volume = config.side ** 3
    if type(mask) is not int or not 0 <= mask < 1 << (30 * volume):
        raise ValueError('occupation outside full Fock carrier')
    if not isinstance(links, tuple) or len(links) != 3 * volume or any(type(v) is not int or not 0 <= v < 60 for v in links):
        raise ValueError('complete sixty-state oriented links required')
    if not isinstance(scalars, tuple) or len(scalars) != volume or any(type(v) is not int or not 0 <= v < 5 for v in scalars):
        raise ValueError('complete five-state quantum scalars required')
    return mask, links, scalars


def column(key, config=Configuration()):
    """Sparse exact action, defined on every full-carrier basis column."""
    mask, links, scalars = validate_key(key, config)
    result = defaultdict(lambda: s.Integer(0))
    result[key] += config.link_penalty * sum(v != 0 for v in links) + config.scalar_penalty * sum(v != 0 for v in scalars)
    for cell, scalar in enumerate(scalars):
        offset = 30 * cell
        local_mask = (mask >> offset) & ((1 << 30) - 1)
        for (target, h), coefficient in local_yukawa(local_mask, scalar,
                (config.y_up, config.y_down, config.y_lepton)).items():
            new_scalars = list(scalars)
            new_scalars[cell] = h
            new_mask = mask ^ ((local_mask ^ target) << offset)
            result[new_mask, links, tuple(new_scalars)] += coefficient
    labels = link_labels()
    link_index = _link_index()
    groups = _native_block_by_index()
    for term in kinetic_terms(config):
        blank = all(links[e] == 0 for e in term.path)
        decoded = [labels[links[e]] for e in term.path]
        reversible = all(a >= 0 for a, _ in decoded) and all(decoded[i][1] == decoded[i + 1][0] for i in range(len(decoded) - 1))
        for (spin_left, spin_right), spin_weight in PAULI[term.axis].todok().items():
            coefficient = term.coefficient * spin_weight
            if blank:
                for b in range(15):
                    origin = 30 * term.right_cell + 15 * spin_right + b
                    if not mask >> origin & 1:
                        continue
                    for a in groups[b]:
                        destination = 30 * term.left_cell + 15 * spin_left + a
                        transition = sequence(mask, ((destination, True), (origin, False)))
                        if transition is None:
                            continue
                        for internal in product(groups[b], repeat=term.displacement - 1):
                            indices = (a, *internal, b)
                            new_links = list(links)
                            for j, edge in enumerate(term.path):
                                new_links[edge] = link_index[indices[j], indices[j + 1]]
                            result[transition[0], tuple(new_links), scalars] += coefficient * transition[1]
            if reversible:
                a, b = decoded[0][0], decoded[-1][1]
                origin = 30 * term.left_cell + 15 * spin_left + a
                destination = 30 * term.right_cell + 15 * spin_right + b
                transition = sequence(mask, ((destination, True), (origin, False)))
                if transition is not None:
                    new_links = list(links)
                    for edge in term.path:
                        new_links[edge] = 0
                    result[transition[0], tuple(new_links), scalars] += s.conjugate(coefficient) * transition[1]
    return clean(result)


def apply_state(state, action):
    result = defaultdict(lambda: s.Integer(0))
    for key, coefficient in state.items():
        for target, value in action(key).items():
            result[target] += coefficient * value
    return clean(result)


def norm2(state):
    return s.simplify(sum(s.conjugate(v) * v for v in state.values()))


@lru_cache(maxsize=None)
def translate_fermions(mask, axis, side):
    occupied = [j for j in range(30 * side ** 3) if mask >> j & 1]
    transformed = []
    wraps = 0
    for mode in occupied:
        cell, internal = divmod(mode, 30)
        xyz = cells(side)[cell]
        wraps += int(xyz[axis] == side - 1)
        transformed.append(30 * cell_index(side, shift_cell(side, xyz, axis)) + internal)
    inversions = sum(a > b for j, a in enumerate(transformed) for b in transformed[j + 1:])
    return sum(1 << i for i in transformed), s.Integer((-1) ** (wraps + inversions))


def translation_column(key, axis, config=Configuration()):
    if axis not in range(3):
        raise ValueError('three spatial translation directions')
    mask, links, scalars = validate_key(key, config)
    n = config.side
    transformed, phase = translate_fermions(mask, axis, n)
    new_links = [0] * len(links)
    new_scalars = [0] * len(scalars)
    for cell, xyz in enumerate(cells(n)):
        target = cell_index(n, shift_cell(n, xyz, axis))
        new_scalars[target] = scalars[cell]
        for edge_axis in range(3):
            new_links[3 * target + edge_axis] = links[3 * cell + edge_axis]
    return {(transformed, tuple(new_links), tuple(new_scalars)): phase}


@lru_cache(maxsize=None)
def _current_columns(current, endpoint=None):
    if endpoint is None:
        matrix = representation()[2][current]
    elif endpoint == 'scalar':
        matrix = higgs_memory()[1][current]
    else:
        matrix = link_generator(representation()[2][current], endpoint)
    columns = defaultdict(list)
    for (row, col), coefficient in s.ImmutableSparseMatrix(matrix).todok().items():
        columns[col].append((row, coefficient))
    return dict(columns)


def gauss_column(key, cell, current, config=Configuration()):
    mask, links, scalars = validate_key(key, config)
    n = config.side
    if cell not in range(n ** 3) or current not in range(12):
        raise ValueError('cell/current outside source')
    result = defaultdict(lambda: s.Integer(0))
    matter = _current_columns(current)
    for spin in range(2):
        for b, transitions in matter.items():
            origin = 30 * cell + 15 * spin + b
            if not mask >> origin & 1:
                continue
            for a, coefficient in transitions:
                step = sequence(mask, ((30 * cell + 15 * spin + a, True), (origin, False)))
                if step is not None:
                    result[step[0], links, scalars] += coefficient * step[1]
    for target, coefficient in _current_columns(current, 'scalar').get(scalars[cell], ()):
        new_scalars = list(scalars)
        new_scalars[cell] = target
        result[mask, links, tuple(new_scalars)] += coefficient
    xyz = cells(n)[cell]
    for axis in range(3):
        for edge, endpoint in ((3 * cell + axis, 0),
                (3 * cell_index(n, shift_cell(n, xyz, axis, -1)) + axis, 1)):
            for target, coefficient in _current_columns(current, endpoint).get(links[edge], ()):
                new_links = list(links)
                new_links[edge] = target
                result[mask, tuple(new_links), scalars] += coefficient
    return clean(result)


def pair_preparation(cell=0, config=Configuration()):
    """Local lepton-pair/scalar pulse; later idle evolution retains every field."""
    if cell not in range(config.side ** 3):
        raise ValueError('cell outside source')
    _, links, scalars = empty_key(config)
    local = local_yukawa(0, 0, (s.Integer(0), s.Integer(0), s.Integer(1)))
    assert norm2(local) == 4
    state = {}
    for (mask, h), coefficient in local.items():
        new_scalars = list(scalars)
        new_scalars[cell] = h
        state[mask << (30 * cell), links, tuple(new_scalars)] = coefficient / 2
    assert norm2(state) == 1
    return state


def phase_projector(state, labels, config=Configuration()):
    if len(labels) != 3 or any(type(k) is not int or not 0 <= k < 2 * config.side for k in labels):
        raise ValueError('one complete translation-phase label per spatial axis')
    result = dict(state)
    n = config.side
    for axis, label in enumerate(labels):
        total = defaultdict(lambda: s.Integer(0))
        power = result
        for m in range(2 * n):
            coefficient = s.simplify(s.exp(-s.I * s.pi * label * m / n)) / (2 * n)
            for key, value in power.items():
                total[key] += coefficient * value
            power = apply_state(power, lambda key: translation_column(key, axis, config))
        result = clean(total)
    return result


def descriptor(config=Configuration()):
    volume = config.side ** 3
    return {'schema': 'oph.chiral_periodic_quantum_family.v1',
            'native_site_program': source_descriptor()['program'],
            'side': config.side, 'spatial_directions': 3, 'cells': volume,
            'oriented_quantum_links': 3 * volume,
            'full_dimension': f'2^{30 * volume} * 60^{3 * volume} * 5^{volume}',
            'full_Fock_modes_per_cell': 30, 'scalar_states_per_cell': 5,
            'spacing': str(config.spacing), 'length': str(config.length),
            'hbar': '1', 'clock': 'exp(-it H_periodic) with the declared full Hamiltonian',
            'coefficients': {'y_up': str(config.y_up), 'y_down': str(config.y_down),
                'y_lepton': str(config.y_lepton), 'scalar_penalty': str(config.scalar_penalty),
                'link_penalty': str(config.link_penalty)},
            'ordered_forward_path_terms_before_adjoint': 3 * volume * (config.side - 1),
            'each_path_has_half_the_AP_spectral_coefficient': True,
            'primitive_quantum_path_normalized_or_replaced_by_identity': False,
            'native_block_dimensions': [len(block) for block in native_blocks()],
            'finite_range_or_finite_exact_lightcone_claim': False,
            'free_ground_or_open_two_cell_ground_transferred': False,
            'interacting_continuum_or_physical_particle_mass': False}
