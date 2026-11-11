"""Exact neutral composite, lepton-number selection and unit Majorana term."""
from functools import lru_cache
import sympy as s

from model.defect_mass_join.source import car, one_cell_polynomial


def lepton_number():
    weights = [-1 if i == 9 else 1 if i in (10, 11) else 0 for i in range(15)]
    return car.one_body(s.kronecker_product(s.eye(2), s.diag(*weights)))


@lru_cache(maxsize=1)
def majorana_operator():
    up, down = (car.neutral_composite(spin) for spin in (0, 1))
    ordered = car.scale(car.add(car.multiply(up, down), car.scale(car.multiply(down, up), -1)), s.Rational(1, 2))
    return ordered, car.add(ordered, car.adjoint(ordered))


@lru_cache(maxsize=1)
def certificate():
    h, number = one_cell_polynomial(), lepton_number()
    assert not car.commutator(number, h)
    for spin in (0, 1):
        j = car.neutral_composite(spin)
        assert car.commutator(number, j) == car.scale(j, -1)
        assert car.commutator(number, car.adjoint(j)) == car.adjoint(j)
        for current in range(12):
            assert not car.commutator(car.gauge_current(current), j)
    pair, perturbation = majorana_operator()
    assert pair and perturbation == car.adjoint(perturbation)
    assert car.commutator(number, pair) == car.scale(pair, -2)
    assert car.commutator(number, perturbation)
    for current in range(12):
        assert not car.commutator(car.gauge_current(current), perturbation)
    for axis in range(3):
        assert not car.commutator(car.spin_current(axis), perturbation)
    # Every CAR word and scalar matrix unit is a contraction. This bound
    # applies to every native occupation, without an assumed cubic identity.
    norm_bound = sum(abs(value) for value in perturbation.values())
    return {'native_neutral_operator': 'J_s=phi_0 c_(s,11)-phi_1 c_(s,10)',
        'native_all_unit_H_preserves_total_lepton_number': True,
        'neutral_composite_lepton_charge': -1,
        'full_native_current_commutators_per_spin': 12,
        'native_invariant_reference_anomalous_JJ_correlator': '0 at every time by exact lepton number',
        'native_Feshbach_elimination_can_generate_DeltaL_two': False,
        'neutral_composite_is_elementary_free_neutrino': False,
        'declared_unit_Majorana_architecture': {
            'pair': 'M=(J_up J_down-J_down J_up)/2',
            'added_Hamiltonian': 'M+M^dagger with primitive coefficient1',
            'literal_global_cap_field_products_used': True,
            'nonzero_normal_ordered_pair_terms': len(pair),
            'nonzero_normal_ordered_Hermitian_terms': len(perturbation),
            'full_gauge_current_commutators': 12,
            'full_spin_current_commutators': 3,
            'pair_lepton_charge': -2,
            'entire_Fock_operator_norm_upper': str(norm_bound),
            'quark_operators_present': False,
            'empirical_Weinberg_coefficient_or_scale_inserted': False,
            'modified_entire_physical_ground_proved': False}}
