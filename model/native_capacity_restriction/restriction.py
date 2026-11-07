"""Full number-mixture corollary of the native canonical-sector readout.

Population and occupation are declared integers. No physical capacity law,
preparation protocol or dynamical refinement map is selected here.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
import hashlib

from ..observer_time_fields import canonical_sector as owning
from ..native_program import interpreter as native


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('finite nonnegative integer required')
    return value


def _domain(population, occupied, retained):
    _integer(population, 1); _integer(occupied); _integer(retained)
    if max(occupied, retained) > population:
        raise ValueError('occupation and retained size cannot exceed population')


def falling(number, degree):
    """Reuse the owning exact degree-four helper and extend its finite degree."""
    _integer(number); _integer(degree)
    if degree <= 4:
        return owning.falling(number, degree)
    return prod(range(number-degree+1, number+1)) if degree <= number else 0


def pattern_weight(population, occupied, retained, ones):
    """Probability of one specified retained binary pattern with `ones` ones."""
    _domain(population, occupied, retained); _integer(ones)
    if ones > retained:
        raise ValueError('pattern occupation exceeds retained size')
    if retained == 2 and 2 <= population <= 1000000:
        # This is the same two-port source marginal used by the owning
        # canonical incident-clock field calculation.
        bits = (int(ones >= 1), int(ones == 2))
        return owning.initial_weight(population, occupied, bits)
    return F(falling(occupied, ones)*falling(population-occupied, retained-ones),
             falling(population, retained))


def number_law(population, occupied, retained):
    """All retained number sectors, including zeros, in O(retained) ratios."""
    _domain(population, occupied, retained)
    first = max(0, occupied-(population-retained))
    last = min(retained, occupied)
    law = [F()]*(retained+1)
    law[first] = comb(retained, first)*pattern_weight(population, occupied, retained, first)
    for count in range(first, last):
        law[count+1] = law[count]*F((retained-count)*(occupied-count),
                                  (count+1)*(population-occupied-retained+count+1))
    return tuple(law)


def nested_law(population, occupied, middle, retained):
    """Compose both actual restrictions, preserving every intermediate sector."""
    _domain(population, occupied, middle); _integer(retained)
    if retained > middle:
        raise ValueError('nested retained size exceeds the intermediate size')
    if middle == 0:
        return (F(1),)
    weights = number_law(population, occupied, middle)
    output = [F()]*(retained+1)
    for count, weight in enumerate(weights):
        if weight:
            for final, conditional in enumerate(number_law(middle, count, retained)):
                output[final] += weight*conditional
    return tuple(output)


def canonical_sector_distance(population, occupied, retained, coarse_occupied):
    """Exact half trace distance to the normalized canonical coarse sector.

    The same expression is only a lower bound for an arbitrary normalized
    state supported in that sector.
    """
    _domain(population, occupied, retained); _integer(coarse_occupied)
    if coarse_occupied > retained:
        raise ValueError('coarse occupation exceeds the retained size')
    return 1-comb(retained, coarse_occupied)*pattern_weight(
        population, occupied, retained, coarse_occupied)


def number_variance(population, occupied, retained):
    _domain(population, occupied, retained)
    return (F(retained*occupied*(population-occupied)*(population-retained),
              population*population*(population-1)) if population > 1 else F())


def vacuum_distance(population, occupied, retained):
    """Exact local vacuum distance and the uniform native occupancy bound."""
    _domain(population, occupied, retained)
    return (1-pattern_weight(population, occupied, retained, 0),
            min(F(1), F(retained*occupied, population)))


def bernoulli_bound(population, retained):
    """Half trace distance bound to density occupied/population product state."""
    _integer(population, 1); _integer(retained)
    if retained > population:
        raise ValueError('retained size exceeds population')
    return min(F(1), F(retained*(retained-1), 2*population))


def bernoulli_law(retained, density):
    _integer(retained)
    if type(density) is not F or not 0 <= density <= 1:
        raise ValueError('exact rational density in the closed unit interval required')
    return tuple(comb(retained, k)*density**k*(1-density)**(retained-k)
                 for k in range(retained+1))


def federation_population(level):
    _integer(level)
    return 240 << (2*level)


def _canonical_density(population, occupied):
    if population == 0:
        return {(0, 0, 0, 0): F(1)}
    return {(mask, 0, mask, 0): F(1, comb(population, occupied))
            for chosen in combinations(range(population), occupied)
            for mask in (sum(1 << i for i in chosen),)}


def _window_density(population, occupied, retained):
    return {(mask, 0, mask, 0): value for mask in range(1 << retained)
            if (value := pattern_weight(population, occupied, retained, mask.bit_count()))}


def exact_restriction_certificate(max_population=8):
    """Direct sparse-matrix replay against the owning occupation partial trace."""
    _integer(max_population, 1)
    if max_population > 10:
        raise ValueError('enumeration check is limited to ten ports')
    cases = nesting = distances = 0
    for population in range(1, max_population+1):
        for occupied in range(population+1):
            density = _canonical_density(population, occupied)
            for retained in range(population+1):
                keep = tuple(range(population-retained, population))
                actual = native.occupation_partial_trace(density, population, keep)
                expected = _window_density(population, occupied, retained)
                if actual != expected:
                    raise ArithmeticError('canonical formula differs from native partial trace')
                law = number_law(population, occupied, retained)
                if sum(law) != 1:
                    raise ArithmeticError('the retained number law loses probability')
                mean = F(retained*occupied, population)
                if (sum(k*p for k, p in enumerate(law)) != mean or
                        sum((k-mean)**2*p for k, p in enumerate(law)) !=
                        number_variance(population, occupied, retained)):
                    raise ArithmeticError('canonical number moments disagree')
                for coarse_occupied in range(retained+1):
                    comparison = _canonical_density(retained, coarse_occupied)
                    keys = actual.keys() | comparison.keys()
                    distance = sum(abs(actual.get(key, F())-comparison.get(key, F())) for key in keys)/2
                    if distance != canonical_sector_distance(population, occupied, retained, coarse_occupied):
                        raise ArithmeticError('orthogonal-sector distance disagrees')
                    distances += 1
                for smaller in range(retained+1):
                    if nested_law(population, occupied, retained, smaller) != number_law(population, occupied, smaller):
                        raise ArithmeticError('full-sector restriction is not nested')
                    nesting += 1
                vacuum, upper = vacuum_distance(population, occupied, retained)
                product_law = bernoulli_law(retained, F(occupied, population))
                distance = sum(abs(a-b) for a, b in zip(law, product_law, strict=True))/2
                if vacuum > upper or distance > bernoulli_bound(population, retained):
                    raise ArithmeticError('finite-population approximation bound failed')
                cases += 1
    return {'abstract_binary_populations_checked': max_population,
            'actual_partial_trace_cases': cases, 'nested_restriction_cases': nesting,
            'canonical_sharp_sector_distances': distances,
            'all_number_sectors_including_zeros_retained': True,
            'enumeration_populations_are_not_geodesic_federations': True}


def _interval(value, digits=12):
    scale = 10**digits
    lower = value.numerator*scale//value.denominator
    return {'lower': str(F(lower, scale)), 'upper': str(F(lower+1, scale))}


def certificate():
    geometry = native.federation_source.federation(0)
    if geometry.ports != federation_population(0):
        raise ArithmeticError('geodesic population does not match the pinned source')
    trajectory = []
    for level in range(3):
        retained = federation_population(level); population = 4*retained
        distance = canonical_sector_distance(population, population//3, retained, retained//3)
        trajectory.append({'coarse_level': level, 'coarse_ports': retained,
                           'fine_ports': population, 'declared_density': '1/3',
                           'canonical_state_half_trace_distance': _interval(distance)})
    population = federation_population(200); occupied = population//3; retained = 12
    law = number_law(population, occupied, retained)
    comparison = bernoulli_law(retained, F(occupied, population))
    distance = sum(abs(a-b) for a, b in zip(law, comparison, strict=True))/2
    bound = bernoulli_bound(population, retained)
    if sum(law) != 1 or distance > bound:
        raise ArithmeticError('large native population restriction failed')
    exact_hash = hashlib.sha256(('|'.join(map(str, law))).encode()).hexdigest()
    return {'schema': 'oph.native_capacity_restriction.v1',
            'restriction': exact_restriction_certificate(),
            'successive_sharp_sector_controls': trajectory,
            'large_population_local_readout': {
                'declared_level': 200, 'population': str(population),
                'population_integer_bits': population.bit_length(),
                'declared_density': '1/3', 'retained_ports': retained,
                'complete_number_outcomes': len(law), 'total_mass': '1',
                'complete_number_law_sha256': exact_hash,
                'Bernoulli_half_trace_distance_upper': str(bound),
                'exact_distance_satisfies_bound': distance <= bound,
                'global_density_matrix_or_native_history_materialized': False},
            'scope': {
                'same_fair_global_canonical_source': True,
                'all_sector_restriction_and_nesting_exact': True,
                'common_local_complete_channel_preserves_state_error_bound': True,
                'correlated_exterior_factorization_assumed': False,
                'arbitrary_reference_extensions_identified_by_source_marginal': False,
                'finite_exact_preparation_or_dynamic_refinement_intertwiner': False,
                'physical_capacity_equals_port_count': False,
                'source_selected_occupation_vs_capacity_trajectory': False,
                'biased_H66_record_clock_has_proved_fair_federation_lift': False,
                'level200_interpreter_execution_or_resource_admission': False,
                'physical_clock_field_or_particle_mass_selected': False}}
