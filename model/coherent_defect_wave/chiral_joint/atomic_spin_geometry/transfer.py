"""Full-carrier spectral isolation from a controlled second-order exchange."""
import sympy as s
from .cube import spectrum
from .momentum import certificate as momentum_certificate
from ..atomic_reference_control.lattice import gap_data, hopping_norm
from ..charged_ground_exchange.resolvent import exchange as exchange_certificate
from ..atomic_even_remainder.bounds import choose_quartic_coupling, uniform_quartic_remainder


def certificate():
    source = gap_data(); tau = hopping_norm(); g = source['gap']
    exchange = exchange_certificate(); spin = spectrum(); momentum = momentum_certificate()
    jlo, jhi = map(s.Rational, exchange['exact_J_interval'])
    gaplo, gaphi = map(s.Rational, spin['gap_interval'])
    next_level = spin['complete_spectrum'][2]
    above = s.Rational(next_level['energy_interval'][0]) + 4
    separation = min(gaplo, above)
    assert separation > s.Rational(4, 5)
    limit = choose_quartic_coupling(g, tau, 2*jlo*separation)
    assert 0 < limit < 1
    v = s.cancel(limit*tau)
    remainder = uniform_quartic_remainder(g, tau, limit)
    gap_lower = s.cancel(2*limit**2*jlo*gaplo - 2*remainder)
    gap_upper = s.cancel(2*limit**2*jhi*gaphi + 2*remainder)
    next_lower = s.cancel(2*limit**2*jlo*above - 2*remainder)
    assert gap_lower > 0 and next_lower > 0 and g-2*v > 0
    positive = [row['momenta_in_pi_units'] for row in momentum['complete_first_triplet_translation_records']
                if row['first_triplet_rank']]
    assert sorted(positive) == [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    full_dimension = (5*2**30)**8 * 119**24
    return {'declared_source': 'H_lambda=H0+lambda*T, N=2, kappa=1, all three Yukawa coefficients=1',
            'all_full_Fock_scalar_link_and_physical_Gauss_sectors_retained': True,
            'actual_hopping_architecture_changed_to_admitted_coefficient': True,
            'admitted_positive_hopping_interval': ['0 (excluded)', str(limit)],
            'reported_nonzero_hopping_coefficient': str(limit),
            'hopping_chosen_from_proof_bounds_without_measured_particle_inputs': True,
            'atomic_separation_lower': str(g), 'full_hopping_norm_upper': str(tau),
            'complete_continued_low_band_rank': 256,
            'second_order_scalar_shift_retained': '-24 lambda^2 A_edge',
            'second_order_spin_scale': '2 lambda^2 J_edge',
            'actual_exchange_interval': [str(jlo), str(jhi)],
            'canonical_graph_full_Q_remainder_at_reported_coefficient': str(remainder),
            'remainder_uses_proved_side_two_bipartite_parity': True,
            'canonical_graph_remainder_formula': '25 (lambda tau)^4/g^3',
            'unique_full_physical_interacting_ground_proved': True,
            'first_excited_cluster_rank': 3,
            'full_physical_ground_gap_interval_at_reported_coefficient': [str(gap_lower), str(gap_upper)],
            'first_triplet_to_other_low_band_levels_gap_lower': str(next_lower),
            'low_band_to_all_other_physical_states_gap_lower': str(g-2*v),
            'exact_ground_momentum_in_pi_units': [0, 0, 0],
            'exact_first_excited_source_momenta_in_pi_units': positive,
            'first_three_levels_equal_by_actual_cubic_rotation_and_translation': True,
            'constant_ranks_follow_from_strict_isolation_for_every_positive_admitted_lambda': True,
            'stationary_reference': 'rank-one projector onto the unique physical ground of H_lambda',
            'formal_full_carrier_ground_filter': {
                'input': 'I/full_dimension', 'successful_joint_Gauss_and_energy_projector_rank': 1,
                'full_dimension': str(full_dimension),
                'success_probability': str(s.Rational(1, full_dimension)),
                'failure_probability': str(1-s.Rational(1, full_dimension)),
                'physical_and_energy_projectors_commute': True,
                'efficient_native_filter_implementation_or_full_histogram_executed': False},
            'unit_hopping_ground_statement_transferred': False,
            'thermodynamic_or_relativistic_limit': False,
            'lowest_gap_is_a_zero_momentum_particle_rest_mass': False,
            'three_momenta_are_three_empirical_generations': False,
            'named_species_or_GeV_masses': False}
