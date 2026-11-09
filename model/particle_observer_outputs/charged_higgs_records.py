"""Charged addition/removal records of the same finite quantum-Higgs source.

The complete probe carries the compensating gauge representation. Bare
charged field operators are not asserted to be standalone physical pointers.
The calculated finite energy differences are not relativistic particle masses.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations
import json

import sympy as s

from . import higgs_records as neutral


@lru_cache(maxsize=8192)
def _rational_simplify(value):
    """Exact rational-function cancellation without algebraic sign searches."""
    return s.cancel(s.expand(value))


def _basis(number, scalar):
    return tuple((scalar, sum(1 << i for i in occupied))
                 for occupied in combinations(range(18), number))


class _NeutralSector(dict):
    """Literal source block with lazily requested large spectral projectors.

    Charged pole coefficients need only the exact source norm and number
    energies. Constructing all 307 neutral projector columns for that scalar
    calculation is unnecessary. Requests for any such column use the owning
    neutral implementation unchanged.
    """
    def __missing__(self, key):
        self.update(neutral.sector_data(self['Y'], self['pair_cost']))
        if key not in self:
            raise KeyError(key)
        return self[key]


def _neutral_source_block(y, d):
    h, basis = neutral.source.zero_charge_block(y, pair_cost=d)
    vacuum = s.ImmutableSparseMatrix(307, 1, {(0, 0): 1})
    v = h*vacuum
    g2 = _rational_simplify(4*s.trace(y.H*y))
    if g2.is_positive is not True or _rational_simplify((v.H*v)[0]-g2) != 0:
        raise ValueError('a nonzero exact source coupling is required')
    delta = s.sqrt(d*d+4*g2)
    return _NeutralSector(Y=y, pair_cost=d, H=h, basis=basis, vacuum=vacuum,
        bright_vector=v, g_squared=g2, gap=delta, ground_energy=(d-delta)/2,
        excited_energy=(d+delta)/2, dark_energy=d)


@lru_cache(maxsize=6)
def _charged_data(y, d, fermion_energy):
    ground = _neutral_source_block(y, d)
    upper = _basis(1, 0)+_basis(3, 1)+_basis(3, 2)
    lower = _basis(1, 1)+_basis(1, 2)
    index = {state: i for i, state in enumerate(upper)}
    entries = {}
    for col in range(18):
        for state, value in neutral.source.vertex_action(y, 0, 1 << col).items():
            row = index[state]-18
            if row < 0:
                raise ArithmeticError('the actual source vertex left its scalar/pair block')
            entries[row, col] = value
    w = s.ImmutableSparseMatrix(1632, 18, entries)
    gram = neutral._matrix(w.H*w)
    g2 = ground['g_squared']
    left, right = y*y.H, 2*(y.H*y).T
    density = s.MutableSparseMatrix(18, 18, {})
    for f in range(3):
        for g in range(3):
            for weak in range(2):
                for spin in range(2):
                    density[6*f+2*weak+spin, 6*g+2*weak+spin] = left[f, g]
            for spin in range(2):
                density[6*f+4+spin, 6*g+4+spin] = right[f, g]
    density = neutral._matrix(density)
    if not neutral._zero(gram-(g2*s.eye(18)-density)):
        raise ArithmeticError('the literal full-CAR Gram differs from its density complement')
    if neutral._clean(s.trace(density)-2*g2) != 0:
        raise ArithmeticError('the pair density has the wrong total occupation')
    he = {(i, i): fermion_energy+(d if i >= 18 else 0) for i in range(1650)}
    for (row, col), value in w.todok().items():
        he[18+row, col] = value
        he[col, 18+row] = s.conjugate(value)
    plus_h = s.ImmutableSparseMatrix(1650, 1650, he)
    return {'neutral': ground, 'Y': y, 'pair_cost': d, 'fermion_energy': fermion_energy,
            'scalar_energy': d-2*fermion_energy,
            'plus_basis': upper, 'minus_basis': lower,
            'W_one': w, 'gram': gram, 'pair_density': density,
            'plus_H': plus_h, 'minus_energy': d-fermion_energy,
            'plus_dark_energy': fermion_energy+d,
            'plus_dark_rank': 1614, 'plus_bright_rank': 36}


def charged_data(y, pair_cost=1, fermion_energy=0):
    y = s.ImmutableMatrix(y)
    if y.shape != (3, 3):
        raise ValueError('a three-by-three source family coefficient is required')
    y = s.ImmutableMatrix(y.applyfunc(neutral._exact))
    d = neutral._exact(pair_cost)
    ef = neutral._exact(fermion_energy)
    if ef.is_real is not True or d.is_real is not True:
        raise ValueError('the number energies must be real')
    return dict(_charged_data(y, d, ef))


def coherent_charged_data(pass_amplitude='1', phase=s.S.One, pair_cost=1, fermion_energy=0):
    y = neutral.attenuation.family(pass_amplitude, phase)['Y']
    return charged_data(y, pair_cost, fermion_energy)


def creation_removal_maps(data):
    """Every actual CAR component from all 307 inputs to both complete sectors.

    Coherent addition components pair with a dual charged probe; removal
    components pair with a direct probe. The probe components are retained
    quantum coordinates, not individually committed gauge-dependent labels.
    """
    plus_index = {state: row for row, state in enumerate(data['plus_basis'])}
    minus_index = {state: row for row, state in enumerate(data['minus_basis'])}
    additions, removals = [], []
    for mode in range(18):
        add, remove = {}, {}
        for col, (scalar, mask) in enumerate(data['neutral']['basis']):
            created = neutral.source._create(mask, mode)
            if created is not None:
                destination, sign = created
                add[plus_index[scalar, destination], col] = sign
            removed = neutral.source._annihilate(mask, mode)
            if removed is not None:
                destination, sign = removed
                remove[minus_index[scalar, destination], col] = sign
        additions.append(s.ImmutableSparseMatrix(1650, 307, add))
        removals.append(s.ImmutableSparseMatrix(36, 307, remove))
    return tuple(additions), tuple(removals)


def ground_probe_vectors(data):
    """Columns are c_i^dagger|G> and c_i|G>, with their original norms."""
    ground = data['neutral']
    ket = neutral.preparation_unitary(ground)*ground['vacuum']
    additions, removals = creation_removal_maps(data)
    plus = s.ImmutableSparseMatrix.hstack(*(a*ket for a in additions))
    minus = s.ImmutableSparseMatrix.hstack(*(a*ket for a in removals))
    return {'addition': plus, 'removal': minus, 'ground': ket}


def addition_dark_probability(data, state_vector):
    """Original complete-probe weight of the retained 1614-dimensional block."""
    state = s.ImmutableSparseMatrix(state_vector)
    if state.shape != (307, 1) or neutral._clean((state.H*state)[0]-1) != 0:
        raise ValueError('a normalized full-sector pure input is required')
    additions, _ = creation_removal_maps(data)
    full = s.ImmutableSparseMatrix.hstack(*(a*state for a in additions))
    occupied = full[18:, :]
    compressed = data['W_one'].H*occupied
    return neutral._clean((s.trace(occupied.H*occupied)
        -s.trace(compressed.H*data['gram'].inv()*compressed))/18)


def addition_resolvent(data, absolute_energy):
    """The exact 18x18 addition Green matrix, not an isolated pole ansatz.

    The result equals F^dagger(z-H_plus)^(-1)F wherever that inverse exists.
    Absolute energy z and relative gap omega=z-E_ground are distinguished.
    """
    z = neutral._exact(absolute_energy)
    d, ef = data['pair_cost'], data['fermion_energy']
    ground = data['neutral']
    delta, g2 = ground['gap'], ground['g_squared']
    alpha, beta = (1+d/delta)/2, (1-d/delta)/2
    x = z-ef
    denominator = x*(x-d)*s.eye(18)-data['gram']
    if denominator.det() == 0:
        raise ValueError('the chosen energy lies on a charged pole')
    numerator = alpha*(x-d)*s.eye(18)+(beta*x/g2-2/delta)*data['gram']
    return neutral._matrix(numerator*denominator.inv())


def response_lines(data):
    """All addition/removal energies and original complete-probe weights.

    Squared singular values are computed from the declared source Y. A
    source matrix requiring a harder exact root representation can use the
    operator resolvent without evaluating this optional eigenvalue listing.
    """
    ground = data['neutral']
    g2, delta = ground['g_squared'], ground['gap']
    d, ef = data['pair_cost'], data['fermion_energy']
    alpha, beta = (1+d/delta)/2, (1-d/delta)/2
    family_gram = data['Y'].H*data['Y']
    diagonal, offdiagonal = family_gram[0, 0], family_gram[0, 1]
    if all(_rational_simplify(family_gram[i, j]-(diagonal if i == j else offdiagonal)) == 0
           for i in range(3) for j in range(3)):
        # The phase-zero coherent source has this exact singlet/plane split.
        # Its degeneracy is retained without solving an expanded cubic.
        eigenvalues = {}
        for value, multiplicity in ((_rational_simplify(diagonal-offdiagonal), 2),
                                    (_rational_simplify(diagonal+2*offdiagonal), 1)):
            eigenvalues[value] = eigenvalues.get(value, 0)+multiplicity
    else:
        eigenvalues = family_gram.eigenvals()
    if sum(eigenvalues.values()) != 3:
        raise ArithmeticError('the family singular values are incomplete')
    rows = []
    for sigma2, family_rank in eigenvalues.items():
        sigma2 = _rational_simplify(sigma2)
        for channel, factor, internal_spin in (('left_doublet', 1, 4), ('charged_singlet', 2, 2)):
            density = factor*sigma2
            lam = _rational_simplify(g2-density)
            split = s.sqrt(d*d+4*lam)
            total_add = _rational_simplify(alpha+beta*lam/g2)
            contrast = _rational_simplify(-4*lam/delta-d*(alpha-beta*lam/g2))
            multiplicity = int(family_rank)*internal_spin
            weights = {'lower': _rational_simplify((total_add-contrast/split)/2),
                       'upper': _rational_simplify((total_add+contrast/split)/2)}
            lines = []
            for branch, sign in (('lower', -1), ('upper', 1)):
                absolute = ef+(d+sign*split)/2
                lines.append({'branch': branch, 'absolute_energy': _rational_simplify(absolute),
                    'ground_referenced_energy': _rational_simplify(absolute-ground['ground_energy']),
                    'weight_per_normalized_mode': weights[branch],
                    'unconditional_probe_probability': _rational_simplify(multiplicity*weights[branch]/18)})
            removal_weight = _rational_simplify(beta*density/g2)
            if _rational_simplify(sum(weights.values())+removal_weight-1) != 0:
                raise ArithmeticError('a charged mode violated its CAR spectral sum')
            rows.append({'internal_channel': channel, 'sigma_squared': sigma2,
                         'multiplicity': multiplicity, 'one_body_pair_density': density,
                         'W_gram_eigenvalue': lam, 'addition_lines': lines,
                         'removal_absolute_energy': d-ef,
                         'removal_ground_referenced_energy': _rational_simplify(d-ef-ground['ground_energy']),
                         'removal_weight_per_normalized_mode': removal_weight,
                         'removal_unconditional_probe_probability': _rational_simplify(multiplicity*removal_weight/18)})
    total = sum((sum(line['unconditional_probe_probability'] for line in row['addition_lines'])
                 +row['removal_unconditional_probe_probability'] for row in rows), s.S.Zero)
    if _rational_simplify(total-1) != 0:
        raise ArithmeticError('the complete charged probe loses original event weight')
    return rows


def verify_full_vertex(data):
    """Every1650 addition input and every36 removal input uses the literal V."""
    plus_index = {state: row for row, state in enumerate(data['plus_basis'])}
    columns = [{} for _ in range(1650)]
    for (row, col), value in data['plus_H'].todok().items():
        if row != col:
            columns[col][row] = value
    for col, state in enumerate(data['plus_basis']):
        expected = {plus_index[out]: coefficient
                    for out, coefficient in neutral.source.vertex_action(data['Y'], *state).items()}
        actual = columns[col]
        if any(neutral._clean(expected.get(row, 0)-actual.get(row, 0)) != 0
               for row in expected.keys() | actual.keys()):
            raise ArithmeticError('an actual addition-sector CAR input differs')
    for state in data['minus_basis']:
        if neutral.source.vertex_action(data['Y'], *state):
            raise ArithmeticError('the removal sector unexpectedly has a pair transition')
    return {'complete_addition_input_columns': 1650, 'complete_removal_input_columns': 36}


def _serialize(value):
    if isinstance(value, dict):
        return {k: _serialize(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serialize(v) for v in value]
    if isinstance(value, s.Basic):
        return str(value)
    return value


def certificate():
    data = charged_data(s.diag(1, 2, 3), 1, 0)
    source_replay = verify_full_vertex(data)
    additions, removals = creation_removal_maps(data)
    add_effect = sum((a.H*a for a in additions), s.zeros(307))/18
    remove_effect = sum((a.H*a for a in removals), s.zeros(307))/18
    if not neutral._zero(add_effect+remove_effect-s.eye(307)):
        raise ArithmeticError('the gauge-compensated probe is not complete')
    vectors = ground_probe_vectors(data)
    beta = (1-data['pair_cost']/data['neutral']['gap'])/2
    removal_probability = neutral._clean(s.trace(vectors['removal'].H*vectors['removal'])/18)
    if removal_probability != beta/9:
        raise ArithmeticError('the complete removal probability differs')
    coherent = coherent_charged_data()
    dark_probe = neutral.gauge_invariant_dark_probe(data['neutral'])['vector']
    dark_weight = addition_dark_probability(data, dark_probe)
    if dark_weight.is_positive is not True:
        raise ArithmeticError('the retained dark response control vanished')
    return _serialize({'schema': 'oph.finite_charged_higgs_records.v1',
        'literal_source_replay': source_replay,
        'fixture': {'Y': 'diag(1,2,3)', 'fermion_energy': '0', 'scalar_energy': '1', 'measured_inputs': []},
        'complete_dimensions': {'ground_sector': 307, 'addition_sector': 1650, 'removal_sector': 36,
            'addition_bright': 36, 'addition_dark': 1614, 'ambient_full_CAR_scalar': 786432},
        'Gram_identity': {'left_doublet_each_spin': 'g^2 I-Y Y^dagger',
            'charged_singlet_each_spin': 'g^2 I-2(Y^dagger Y)^T',
            'uniform_lower_bound': 'W^dagger W >= g^2 I/2 > 0'},
        'complete_probe': {'addition': 'sum_i c_i^dagger tensor |i dual>/sqrt(18)',
            'removal': 'sum_i c_i tensor |i>/sqrt(18)',
            'charged_reference_coordinates_retained_coherently': True,
            'bare_charged_pointer_claimed_gauge_invariant': False,
            'addition_probability': 1-removal_probability,
            'removal_probability': removal_probability,
            'postselection': False},
        'fixture_response_lines': response_lines(data),
        'coherent_source_response_lines': response_lines(coherent),
        'retained_zero_ground_addition_dark_weight': {'dark_rank': 1614, 'weight': 0},
        'neutral_dark_input_addition_dark_probability': dark_weight,
        'probe_energy_convention': 'zero probe Hamiltonian during the declared field-energy readout',
        'scope': {'actual_interacting_charged_spectral_response': True,
            'number_energy_parameters_are_declared': True,
            'prepared_ground_is_Q_zero_sector_ground': True,
            'global_ambient_ground_or_all_charge_sector_spectrum': False,
            'same_Hamiltonian_clock_as_neutral_readout': True,
            'native_pulse_tape_executed_here': False,
            'relativistic_charged_particle_or_antiparticle_identity': False,
            'physical_mass_or_GeV_assignment': False}})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
