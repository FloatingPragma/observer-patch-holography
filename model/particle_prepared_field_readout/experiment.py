"""Compose preparation, locality, representation cutoff and a complete detector.

The large state is specified by finite words and exact integrals; this module
calculates resource and probability enclosures without allocating that state.
All fractions are in one declared Hamiltonian clock, not physical units.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import comb

import sympy as s

from model.particle_field_preparation import preparation
from model.particle_field_causality import bounds
from model.particle_observer_outputs import field_route, readout_contract


SOURCE_INPUTS = (
    'model/particle_field_preparation/preparation.py',
    'model/particle_field_preparation/PROOF.md',
    'model/particle_field_preparation/REFINEMENT.md',
    'model/particle_field_causality/bounds.py',
    'model/particle_field_causality/PROOF.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_observer_outputs/readout_contract.py',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
)


def _integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' requires an integer at least '+str(minimum))
    return value


def _fraction(value, name, positive=False):
    if type(value) not in (int, str, F):
        raise ValueError(name+' requires exact rational data')
    result = F(value)
    if result < 0 or (positive and not result):
        raise ValueError(name+' has the wrong sign')
    return result


def marginal_support(radius, per_link_degree):
    """Support of the actual periodic preparation's entire box marginal.

    Tracing the exterior preserves each local ket word space. The density
    matrix's bra is not an extra ket degree. Boundary charges are retained
    by a conjugate preparation reference, not projected to box singlets.
    """
    radius = _integer(radius, 'box radius', 1)
    degree = _integer(per_link_degree, 'per-link preparation degree', 1)
    cells = (2*radius+1)**3
    links = 3*cells
    return {
        'box_radius': radius, 'box_cells': cells, 'retained_outgoing_links': links,
        'initial_total_word_degree_upper': links*degree,
        'per_link_ket_word_degree_upper': degree,
        'full_Fock_dimension': {'base': 2, 'exponent': 60*cells},
        'initial_fermion_number': 60,
        'initial_fermions': 'all60 modes filled at origin, all other modes empty',
        'gauge_marginal_rank_upper': {'base': 451, 'exponent': links*degree},
        'conjugate_purifier_dimension_upper': {'base': 451, 'exponent': links*degree},
        'input': 'entire open-box marginal of the specified periodic preparation',
        'boundary_outgoing_links_are_electric_spectators': True,
        'box_only_singlet_projection_performed': False,
        'bra_degree_added_to_initial_ket_cutoff': False,
        'simulator_purifier_identified_with_dynamical_periodic_exterior': False,
    }


def small_clock_volume_bound(constants, duration, radius, periodic_side):
    """Rational envelope of the owning locality theorem on an interior torus box.

    The conservative periodic collar prevents a wraparound shortcut. The
    bound compares local final records, retaining an untouched external
    reference. It does not compare the dynamical exterior's complete state.
    """
    t = _fraction(duration, 'duration')
    radius = _integer(radius, 'box radius', 1)
    side = _integer(periodic_side, 'periodic side', 1)
    d = constants['root_support_radius']
    collar = max(d, 2)
    required_side = 2*(radius+2*collar+d)+1
    if side % 2 != 1 or side < required_side:
        raise ValueError('periodic wraparound intersects the certified comparison collar')
    # At t=0 the owning routine supplies its exact collar and crossing norm
    # without constructing an exponential with large rational denominators.
    geometry = bounds.volume_bound(constants, F(0), radius)
    rate = constants['weighted_local_interaction_upper']
    x = 2*rate*t
    if x > F(1, 2):
        raise ValueError('the rational exponential envelope requires 2 L t <= 1/2')
    # exp(x)-1 <= x/(1-x) <= 2x for 0<=x<=1/2.
    distance = geometry['boundary_separation_lower']
    exponential_upper = 2*x
    raw = (geometry['boundary_perturbation_norm_upper']*t
           *2*exponential_upper/F(2**distance))
    return {
        'periodic_side': side, 'required_periodic_side_lower': required_side,
        'box_radius': radius, 'boundary_separation_lower': distance,
        'boundary_perturbation_norm_upper': geometry['boundary_perturbation_norm_upper'],
        'same_plaquette_coefficient': constants['plaquette_norm'],
        'same_duration': t, 'exponential_argument': x,
        'exponential_minus_one_upper': exponential_upper,
        'complete_local_record_TV_upper': min(F(1), raw),
        'same_entire_box_and_external_reference_marginal': True,
        'dynamical_exterior_state_comparison_claimed': False,
        'open_Wilson_gap_assumed': False,
    }


@lru_cache(maxsize=1)
def native_escape_coefficient():
    """Actual m=0 hopping norm out of a completely filled origin cell.

    All six neighbors are distinct. Unitarity of each native M15 link
    cancels its coordinate dependence in the total first-order norm.
    """
    data = field_route.source_data.load()
    gamma, _ = data['gamma']
    exact = [s.Matrix(g).applyfunc(lambda x: s.nsimplify(x, rational=True))
             for g in gamma]
    beta = exact[3]
    rows = []
    total = s.S.Zero
    for direction in range(3):
        for sign in (-1, 1):
            projector = (s.eye(4)+sign*exact[direction])/2
            hopping = -s.Rational(4, 17)*beta*projector
            if projector*projector != projector or projector.H != projector:
                raise ArithmeticError('the literal native spin projector changed')
            if hopping.H*hopping != s.Rational(16, 289)*projector:
                raise ArithmeticError('the literal native hopping coefficient changed')
            coefficient = len(data['charges'])*s.trace(hopping.H*hopping)
            total += coefficient
            rows.append({'direction': direction, 'sign': sign,
                         'spin_rank': int(s.trace(projector)),
                         'internal_dimension': len(data['charges']),
                         'squared_escape_coefficient': str(coefficient)})
    if total != s.Rational(2880, 289):
        raise ArithmeticError('the same-source escape coefficient differs')
    return {'spacing': F(1), 'Legendre_degree': 0,
            'squared_escape_coefficient': F(2880, 289), 'bonds': rows,
            'valid_for_every_normalized_link_state_and_reference': True,
            'origin_occupation_after_one_hopping_insertion': 59}


def occupation_enclosure(duration, finite_H_norm):
    """Complete61-outcome law enclosed by a same-H Taylor remainder.

    The lower bound concerns this declared finite polynomial experiment.
    Its small rare-event probability need not be resolved by the parent
    exact-sign comparison error.
    """
    t = _fraction(duration, 'duration', True)
    norm = _fraction(finite_H_norm, 'finite Hamiltonian norm', True)
    remainder = t*t*norm*norm/2
    if remainder > t/2:
        raise ValueError('the declared tick does not resolve the first hopping term')
    coefficient = native_escape_coefficient()['squared_escape_coefficient']
    if not 9 < coefficient < 16:
        raise ArithmeticError('the rational square-root enclosure differs')
    lo59 = (3*t-remainder)**2
    hi59 = (4*t+remainder)**2
    low_tail = remainder**2
    if hi59+low_tail > 1:
        raise ValueError('the probability enclosure requires a smaller tick')
    rows = [{'origin_occupation': n,
             'local_Fock_projector_rank': comb(60, n),
             'probability_interval': (F(0), low_tail)} for n in range(59)]
    rows.extend([
        {'origin_occupation': 59, 'local_Fock_projector_rank': 60,
         'probability_interval': (lo59, hi59)},
        {'origin_occupation': 60, 'local_Fock_projector_rank': 1,
         'probability_interval': (1-hi59-low_tail, 1-lo59)},
    ])
    return {'duration': t, 'finite_H_norm_on_conserved_N60_sector_upper': norm,
            'state_vector_Taylor_remainder_upper': remainder,
            'complete_outcomes': rows,
            'sum_of_all61_probabilities': F(1),
            'joint_probability_origin_occupation_0_through_58_upper': low_tail,
            'strictly_positive_finite_program_p59_lower': lo59,
            'unused_six_bit_pointer_labels': [61, 62, 63],
            'unused_pointer_probabilities': [F(0), F(0), F(0)],
            'conditioning_or_histogram_renormalization': False,
            'large_field_state_or_propagator_materialized': False}


@lru_cache(maxsize=8)
def experiment_plan(radius=5, target_bad_bits=20):
    """One fixed high-plaquette-cost source, marginal, clock and detector.

    The polynomial is m=0 so its nonzero detector event has a literal
    source hopping coefficient. No accuracy of m=0 relativistic dispersion
    is asserted. The owning sign comparison supplies an absolute error.
    """
    radius = _integer(radius, 'box radius', 5)
    bits = _integer(target_bad_bits, 'all-time gap probability bits', 1)
    degree, spacing, number = 0, F(1), 60
    d, collar = 1, 2
    side = 2*(radius+2*collar+d)+1
    parent = preparation.energy_confinement_plan(
        side, bits, degree=degree, spacing=spacing, occupied_modes=number)
    kappa = parent['declared_plaquette_coefficient']
    constants = bounds.interaction_constants(degree, spacing, kappa)
    support = marginal_support(radius, parent['per_link_word_degree'])
    links = support['retained_outgoing_links']
    k0 = support['initial_total_word_degree_upper']
    cutoff, step = k0+12, 4
    bounded = parent['kinetic_operator_norm_upper']+kappa*links
    # A tensor word of total degree <=K has sum_link degree_link²<=K².
    finite_norm = F(29, 10)*cutoff*cutoff+bounded
    rate = constants['weighted_local_interaction_upper']
    tick = min(F(1), 1/(128*finite_norm*finite_norm), 1/(4*rate))
    volume = small_clock_volume_bound(constants, tick, radius, side)
    integrated = bounded*tick
    if integrated > F(1, 128):
        raise ArithmeticError('the selected tick exceeds its cutoff norm allowance')
    # Reuse the owning complete-record theorem. Its tail helper accepts a
    # short exact sum plus a rigorous geometric remainder, so the actual
    # small clock gives a compact bound without forty large-denominator terms.
    cut = field_route.cutoff_bound(cutoff, k0, step, F(1, 128))
    tail = field_route.exponential_tail_upper(
        integrated, cut['first_possible_omitted_Dyson_order'], extra=1)
    cut_error = min(F(1), 2*tail)
    if cut_error > F(cut['complete_record_TV_error_upper']):
        raise ArithmeticError('the actual small-clock cutoff enclosure exceeds its coarse control')
    cut.update({
        'integrated_interaction_norm_upper': str(integrated),
        'state_vector_error_upper': str(min(F(2), 2*tail)),
        'complete_record_TV_error_upper': str(cut_error),
        'exact_tail_summation_terms_past_first': 1})
    sign = preparation.all_time_final_record_bound(parent, tick)
    sign_error = sign['complete_final_record_TV_poly_vs_exact_sign_upper']
    volume_error = volume['complete_local_record_TV_upper']
    total = min(F(1), sign_error+volume_error+cut_error)
    records = occupation_enclosure(tick, finite_norm)
    parent_p59 = (max(F(0), records['complete_outcomes'][59]['probability_interval'][0]-total),
                  min(F(1), records['complete_outcomes'][59]['probability_interval'][1]+total))
    return {
        'periodic_preparation_and_action': parent,
        'open_marginal_and_reference': support,
        'locality_constants': constants,
        'finite_cutoff': {
            **cut, 'actual_integrated_interaction_norm_upper': integrated,
            'bounded_action_norm_on_N60_sector_upper': bounded,
            'same_plaquette_coefficient': kappa,
            'same_duration': tick,
            'coordinate_dimension_upper': {'base': 1+450*links, 'exponent': cutoff},
            'all_Fock_occupations_retained': True,
            'initial_state_supported_without_postselection': True,
            'total_degree_compression_claimed_exactly_local': False},
        'clock': {'macro_tick': tick, 'number_of_identical_ticks': 1,
                  'meaning': 'exp(-i T P_K H_open P_K), one declared finite unitary service',
                  'preparation_and_final_detection_are_separate_services': True,
                  'physical_seconds_or_GeV_assigned': False},
        'periodic_sign_comparison': sign,
        'periodic_to_open_local_record_comparison': volume,
        'finite_complete_occupation_readout': records,
        'comparison_budget': {
            'periodic_exact_sign_to_periodic_polynomial': sign_error,
            'periodic_polynomial_to_open_polynomial': volume_error,
            'open_polynomial_to_finite_word_code': cut_error,
            'total_complete_local_record_TV_upper': total,
            'external_untouched_references_retained': True,
            'dynamical_periodic_exterior_compared': False,
            'preparation_or_compiler_error_included': False},
        'periodic_exact_sign_p59_interval': parent_p59,
        'periodic_exact_sign_rare_event_positivity_resolved': parent_p59[0] > 0,
        'interacting_relativistic_continuum_claimed': False,
        'physical_particle_mass_or_species_claimed': False,
    }


def implementation_budget(plan, *, preparation_error=F(0), step_error=F(0),
                          detector_error=F(0), initial_state_error=F(0)):
    """Add complete finite-service errors with every failure record retained."""
    extra = readout_contract.composition_error(
        ticks=1, detector_uses=1, preparation_error=preparation_error,
        step_error=step_error, detector_error=detector_error,
        initial_state_error=initial_state_error)
    return min(F(1), plan['comparison_budget']['total_complete_local_record_TV_upper']+extra)


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    return value


def certificate():
    plan = experiment_plan()
    if plan['comparison_budget']['total_complete_local_record_TV_upper'] >= F(1, 10**6):
        raise ArithmeticError('the composed finite experiment exceeds its accuracy target')
    return serial({
        'schema': 'oph.prepared_field_readout.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'experiment': plan,
        'native_hopping_check': native_escape_coefficient(),
        'complete_law': 'p_n=Tr[(P_n tensor I) exp(-i T H_K) rho_box exp(+i T H_K)]',
        'preparation': 'finite canonical purification of the actual periodic box marginal',
        'source_preparation_and_apparatus_may_be_declared_primitives': True,
        'narrower_native_repair_tape_implementation_claimed': False,
        'numerically_materialized_source_or_propagator': False,
        'independent_measured_particle_input_used': False,
    })
