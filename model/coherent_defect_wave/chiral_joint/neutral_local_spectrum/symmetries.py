"""Actual local Gauss, spin and centered-charge statements on the odd source."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ...chiral_higgs.model import higgs_memory, PAULI
from ..higgs_composites import model as poly
from .model import seed, graph

FULL_MASK = (1 << 30) - 1
ODD_INDEX_MASK = sum(1 << j for j in range(1, 30, 2))
SCALAR_CONJUGATE = (0, 3, 4, 1, 2)


def apply_polynomial(operator, state):
    out = defaultdict(lambda: s.Integer(0))
    for (mask, scalar), coefficient in state.items():
        for key, value in poly.column(operator, mask, scalar).items():
            out[key] += coefficient * value
    return poly.clean(out)


def conjugate_polynomial(operator):
    """Antiunitary A: c_i -> i c_i^dagger and scalar particle/antiparticle swap."""
    out = defaultdict(lambda: s.Integer(0))
    for (word, row, col), coefficient in operator.items():
        phase = s.Integer(1)
        transformed = []
        for mode, create in word:
            phase *= -s.I if create else s.I
            transformed.append((mode, not create))
        for normal, sign in poly.normal_order(tuple(transformed)):
            out[normal, SCALAR_CONJUGATE[row], SCALAR_CONJUGATE[col]] += s.conjugate(coefficient) * phase * sign
    return poly.clean(out)


def conjugate_key(key):
    mask, scalar = key
    exponent = (-mask.bit_count() + 2 * (mask & ODD_INDEX_MASK).bit_count()) % 4
    return (FULL_MASK ^ mask, SCALAR_CONJUGATE[scalar]), exponent


def spin_flip_key(key):
    mask, scalar = key
    up, down = mask & ((1 << 15) - 1), mask >> 15
    parity = (down.bit_count() + up.bit_count() * down.bit_count()) % 2
    return (down | (up << 15), scalar), (-1) ** parity


def raw_labels(key):
    mask, _ = key
    ell = sum(-((mask >> (15 * spin + 9)) & 1)
              + ((mask >> (15 * spin + 10)) & 1)
              + ((mask >> (15 * spin + 11)) & 1) for spin in range(2))
    baryon3 = sum(sum((mask >> (15 * spin + j)) & 1 for j in (2, 3, 5, 6, 7, 8))
                  - sum((mask >> (15 * spin + j)) & 1 for j in (0, 1, 4, 12, 13, 14))
                  for spin in range(2))
    spin2 = (mask & ((1 << 15) - 1)).bit_count() - (mask >> 15).bit_count()
    return ell, baryon3, spin2


@lru_cache(maxsize=1)
def certificate():
    parts, couplings = poly.quantum_yukawa()
    h = poly.add(*parts.values(), poly.monomial(matrix=higgs_memory()[2]))
    h = {key: value.subs({y: 1 for y in couplings.values()}) for key, value in h.items()}
    h = poly.clean(h)
    actual_seed = poly.column(poly.adjoint(poly.neutral_composite(0)), 0, 0)
    assert actual_seed == seed()
    assert sum(value ** 2 for value in actual_seed.values()) == 2
    assert conjugate_polynomial(h) == h
    for current in range(12):
        gauge = poly.gauge_current(current)
        assert not poly.commutator(gauge, h)
        assert not apply_polynomial(gauge, actual_seed)
        assert conjugate_polynomial(gauge) == poly.scale(gauge, -1)
    for axis in range(3):
        spin = poly.spin_current(axis)
        assert not poly.commutator(spin, h)
        assert conjugate_polynomial(spin) == poly.scale(spin, -1)
    raising = s.zeros(30)
    raising[:15, 15:] = s.eye(15)
    assert not apply_polynomial(poly.one_body(raising), actual_seed)
    assert apply_polynomial(poly.spin_current(2), actual_seed) == {key: value / 2 for key, value in actual_seed.items()}
    charges = [s.Integer(0)] * 15
    charges[9], charges[10], charges[11] = -s.Integer(1), s.Integer(1), s.Integer(1)
    ell = poly.one_body(s.kronecker_product(s.eye(2), s.diag(*charges)))
    assert not poly.commutator(ell, h)
    assert apply_polynomial(ell, actual_seed) == actual_seed
    assert conjugate_polynomial(ell) == poly.add(poly.scale(poly.monomial(), 2), poly.scale(ell, -1))
    centered = poly.add(ell, poly.scale(poly.monomial(), -1))
    assert conjugate_polynomial(centered) == poly.scale(centered, -1)
    # The equal-sign coordinate vector is not the neutral Gauss singlet.
    wrong_seed = {key: abs(value) for key, value in actual_seed.items()}
    assert apply_polynomial(poly.gauge_current(8), wrong_seed)
    keys, _, _, _, _ = graph(); domain = set(keys); image = set()
    for key in keys:
        assert raw_labels(key) == (1, 0, 1)
        assert key[0].bit_count() % 2 == 1
        reflected, exponent = conjugate_key(key)
        assert raw_labels(reflected) == (1, 0, -1)
        back, second_exponent = conjugate_key(reflected)
        assert back == key and (second_exponent - exponent) % 4 == 2
        rotated, _ = spin_flip_key(reflected)
        assert rotated in domain
        image.add(rotated)
    assert image == domain
    return {
        'physical_seed': 'D_up^dagger Omega, two native components with relative minus sign',
        'full_H_native_Gauss_operator_identities': 12,
        'full_H_spin_operator_identities': 3,
        'seed_Gauss_constraints': 12,
        'actual_local_spin_Casimir': '3/4',
        'actual_local_twoSz': '1',
        'full_physical_cyclic_space_has_spin': '1/2',
        'raw_occupation_charge_L': '1',
        'centered_charge_L_minus_one': '0',
        'raw_three_baryon_charge': '0',
        'raw_L_charge_conjugation': 'L -> 2-L',
        'centered_charge_conjugation': 'L-1 -> -(L-1)',
        'raw_L_is_identified_with_measured_lepton_charge': False,
        'full_antiunitary_local_H_identity': True,
        'antiunitary_current_action': 'G_a -> -G_a, S_i -> -S_i',
        'antiunitary_square': 'fermion parity, hence -1 on this odd sector',
        'complete_component_label_and_conjugation_checks': len(keys),
        'charge_conjugation_pairs_ground_within_raw_L_one': True,
        'charge_conjugate_ground_is_opposite_spin_partner': True,
        'charge_conjugation_followed_by_pi_spin_rotation_preserves_component': True,
        'distinct_observed_particle_antiparticle_species_inferred': False,
        'wrong_relative_seed_sign_rejected': True,
        'spin_is_identified_with_a_relativistic_little_group': False,
    }
