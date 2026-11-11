"""A coupled native R16 sterile extension and its complete lepton control.

The source uses left-handed n^c, the charge conjugate of a physical
right-handed sterile field. It is distinct from an active antineutrino.
"""
from collections import defaultdict
from functools import lru_cache
import sympy as s

from model.defect_mass_join.source import car, one_cell_polynomial, representation, higgs_memory
from model.coherent_defect_wave.chiral_joint.model import sequence
from model.coherent_defect_wave.chiral_joint.neutral_local_spectrum.symmetries import conjugate_polynomial

ACTIVE = (9, 10, 11, 15, 25, 26, 27, 31)
STERILE = (15, 31)


def lift(poly):
    return {(tuple(((i % 15)+16*(i//15), creation) for i, creation in word), row, col): value
            for (word, row, col), value in poly.items()}


def current(index):
    return lift(car.gauge_current(index))


def spin_current(axis):
    from model.coherent_defect_wave.chiral_higgs.model import PAULI
    return car.one_body(s.kronecker_product(PAULI[axis]/2, s.eye(16)))


def number():
    weights = [-1 if i in (9, 15) else 1 if i in (10, 11) else 0 for i in range(16)]
    return car.one_body(s.kronecker_product(s.eye(2), s.diag(*weights)))


@lru_cache(maxsize=1)
def interactions():
    j = [lift(car.neutral_composite(spin)) for spin in (0, 1)]
    n = [car.monomial(((index, False),)) for index in STERILE]
    pair = car.add(car.multiply(n[0], j[1]), car.scale(car.multiply(n[1], j[0]), -1))
    dirac = car.add(pair, car.adjoint(pair))
    majorana_pair = car.multiply(n[0], n[1])
    majorana = car.add(majorana_pair, car.adjoint(majorana_pair))
    return dirac, majorana


def column(poly, mask, scalar):
    if type(mask) is not int or not 0 <= mask < 1 << 32 or scalar not in range(5):
        raise ValueError('complete32-mode occupation and five-state scalar required')
    result = defaultdict(lambda: s.Integer(0))
    for (word, row, col), coefficient in poly.items():
        if col == scalar:
            step = sequence(mask, word)
            if step is not None:
                result[step[0], row] += coefficient*step[1]
    return car.clean(result)


def neutral_slice(poly):
    result = s.zeros(32)
    for (word, row, col), coefficient in poly.items():
        if len(word) == 2 and all(create for _, create in word) and (row, col) in ((0, 0), (0, 2), (2, 0)):
            i, j = (atom[0] for atom in word)
            result[i, j] += coefficient
            result[j, i] -= coefficient
    return result


@lru_cache(maxsize=1)
def operator_certificate():
    dirac, majorana = interactions()
    full = car.add(lift(one_cell_polynomial()), dirac)
    assert dirac and full == car.adjoint(full)
    assert not car.commutator(number(), full)
    for index in range(12):
        assert not car.commutator(current(index), full)
        assert not car.commutator(current(index), majorana)
    for axis in range(3):
        assert not car.commutator(spin_current(axis), full)
        assert not car.commutator(spin_current(axis), majorana)
    assert car.commutator(number(), majorana)
    assert conjugate_polynomial(full) == full
    assert conjugate_polynomial(majorana) == majorana
    for spin in (0, 1):
        for probe in (lift(car.neutral_composite(spin)), car.monomial(((STERILE[spin], False),))):
            assert conjugate_polynomial(probe) == car.scale(car.adjoint(probe), s.I)
    eom_signs = []
    for spin in (0, 1):
        probe = car.monomial(((STERILE[spin], False),))
        derivative = car.commutator(probe, full)
        partner = car.adjoint(lift(car.neutral_composite(1-spin)))
        sign = next(sign for sign in (1, -1) if derivative == car.scale(partner, sign))
        eom_signs.append(sign)
    anomalous = car.anticommutator(car.commutator(car.monomial(((15, False),)), car.add(full, majorana)),
                                  car.monomial(((31, False),)))
    assert anomalous == car.scale(car.monomial(), -1)
    # Direct EOM check on every neutral active/sterile contact coordinate.
    neutral = (10, 15, 26, 31)
    kernel = neutral_slice(full)
    for i in neutral:
        derivative = car.commutator(car.monomial(((i, False),)), full)
        for j in neutral:
            insertion = car.anticommutator(derivative, car.monomial(((j, False),)))
            value = sum(v for (word, row, col), v in insertion.items() if (row, col) in ((0, 0), (0, 2), (2, 0)))
            assert value == kernel[i, j]
    assert kernel.rank() == 32 and kernel.extract(neutral, neutral).rank() == 4
    active_mask = sum(1 << i for i in ACTIVE)
    leakage = sum(value*s.conjugate(value) for (mask, scalar), value in column(full, 0, 0).items()
                  if mask & ~active_mask)
    assert leakage == 24
    return {'native_source_CAR_modes': 32, 'full_source_dimension': str(5*2**32),
        'scalar_states': 5,
        'sterile_left_handed_conjugate_modes': list(STERILE),
        'physical_right_handed_sterile_field': 'nu_R is charge conjugate to the native left-handed n^c coordinate',
        'coupled_operator': 'n^c_up J_down-n^c_down J_up+h.c.',
        'primitive_coupling': '1',
        'full_native_gauge_Ward_checks': 12, 'full_spin_Ward_checks': 3,
        'full_H_preserves_lepton_number': True,
        'full_Dirac_and_optional_Majorana_H_antiunitary_symmetry': True,
        'active_and_sterile_probes_map_to_i_adjoint_under_antiunitary': True,
        'exact_sterile_EOM_partner': 'opposite-spin J^dagger',
        'exact_sterile_EOM_signs': eom_signs,
        'lepton_number_of_native_n_c': -1,
        'quantum_Higgs_not_replaced_by_background_in_H': True,
        'neutral_algebraic_contact_slice': {'orientation': 'phi_1=phi_1^dagger=1, phi_0=0; algebraic diagnostic only',
            'full32_antisymmetric_rank': 32, 'active_sterile_neutral_coordinates': list(neutral),
            'neutral_four_mode_rank': 4,
            'neutral_block': [[str(x) for x in row] for row in kernel.extract(neutral, neutral).tolist()],
            'one_family_Dirac_rank': 1, 'three_identical_copies_Dirac_rank': 3,
            'three_copy_hierarchy_or_mixing_derived': False,
            'physical_vacuum_background_or_rest_mass_inferred': False},
        'empty_quarks_invariant_in_full_all_unit_source': False,
        'all_unit_blank_quark_leakage_squared': str(leakage),
        'optional_unit_sterile_Majorana_operator': 'n^c_up n^c_down+h.c.',
        'sterile_Majorana_gauge_and_spin_invariant': True,
        'sterile_Majorana_preserves_lepton_number': False,
        'unit_sterile_Majorana_anomalous_contact': '{[n_up,H],n_down}=-I on the full32CAR/scalar carrier',
        'full32_all_unit_physical_ground_computed': False,
        'flavor_semantics': {'active_tau_antineutrino': 'conjugate weak-doublet charged-current flavor associated with the tau lepton; a superposition of antineutrino mass modes once a mixing map is derived',
            'right_handed_sterile_neutrino': 'weak and hypercharge singlet; distinct from active tau antineutrino',
            'tau_flavor_assigned_to_a_separate_pole': False,
            'observed_mixing_or_oscillation_inputs': []}}


def control_poly(majorana=False):
    dirac, mass = interactions()
    result = car.add(lift(one_cell_polynomial(('charged_lepton',))), dirac, mass if majorana else {})
    address = {mode: i for i, mode in enumerate(ACTIVE)}
    assert all(all(i in address for i, _ in word) for word, _, _ in result)
    return {(tuple((address[i], create) for i, create in word), row, col): value
            for (word, row, col), value in result.items()}


def restricted(poly):
    address = {mode: i for i, mode in enumerate(ACTIVE)}
    result = {}
    for (word, row, col), value in poly.items():
        if any(mode not in address for mode, _ in word):
            # These words are normal ordered. An inactive annihilator kills
            # the exactly empty spectator factor before any creation acts.
            if any(mode not in address and not create for mode, create in word):
                continue
            raise ValueError('operator leaves the declared8-mode control')
        result[tuple((address[i], create) for i, create in word), row, col] = value
    return result


def matrix(poly, domain, target):
    positions = {key: i for i, key in enumerate(target)}
    entries = {}
    for col, key in enumerate(domain):
        for destination, value in column(poly, *key).items():
            if destination not in positions:
                raise ArithmeticError('a full scalar/CAR image left the declared target')
            entries[positions[destination], col] = value
    return s.ImmutableSparseMatrix(len(target), len(domain), entries)


@lru_cache(maxsize=1)
def physical_basis():
    r = representation()[2]
    scalar = higgs_memory()[1]
    ty = [r[11][i % 16, i % 16] if i % 16 < 15 else 0 for i in ACTIVE]
    tw = [r[10][i % 16, i % 16] if i % 16 < 15 else 0 for i in ACTIVE]
    keys, raised = [], []
    for mask in range(256):
        y = sum(ty[i] for i in range(8) if mask >> i & 1)
        w = sum(tw[i] for i in range(8) if mask >> i & 1)
        for level in range(5):
            if y+scalar[11][level, level] == 0:
                weight = w+scalar[10][level, level]
                if weight == 0:
                    keys.append((mask, level))
                elif weight == 1:
                    raised.append((mask, level))
    raising = restricted(car.add(current(8), car.scale(current(9), s.I)))
    raw = matrix(raising, keys, raised)
    basis = s.Matrix.hstack(*raw.nullspace())
    gram = basis.H*basis
    left = gram.inv()*basis.H
    assert left*basis == s.eye(basis.cols)
    return {'keys': tuple(keys), 'basis': basis, 'gram': gram, 'left': left,
            'full_carrier_dimension': 1280, 'weight_zero_dimension': len(keys),
            'complete_physical_dimension': basis.cols}


@lru_cache(maxsize=2)
def physical_control(majorana=False):
    data = physical_basis()
    poly = control_poly(majorana)
    assert conjugate_polynomial(poly) == poly
    raw = matrix(poly, data['keys'], data['keys'])
    h = s.ImmutableSparseMatrix(data['left']*raw*data['basis'])
    assert data['basis']*h == raw*data['basis']
    assert data['gram']*h == h.H*data['gram']
    x = s.Symbol('x')
    characteristic = s.Poly(h.charpoly(x).as_expr(), x)
    minimal = characteristic.sqf_part()
    check = s.zeros(h.rows)
    for coefficient in minimal.all_coeffs():
        check = h*check+coefficient*s.eye(h.rows)
    assert check == s.zeros(h.rows)
    intervals = characteristic.intervals(eps=s.Rational(1, 10**45))
    return {**data, 'H': h, 'characteristic': characteristic, 'minimal': minimal,
            'energy_intervals_with_ranks': intervals, 'majorana': majorana}
