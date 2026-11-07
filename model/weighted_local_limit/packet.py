"""Reproduce the controlled weighted canonical local-limit certificate."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path

from tools import strict_json
from ..weighted_source_family import packet as parent
from . import limit


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
SOURCE_FILES = tuple(sorted(set((*parent.SOURCE_FILES,
    'model/weighted_source_family/results.json',
    'Lean/Thermodynamics/FiniteConditionalRepair.lean',
    *(f'model/weighted_local_limit/{name}' for name in (
        '__init__.py', 'limit.py', 'packet.py', 'test_limit.py', 'test_packet.py',
        'README.md', 'PROOF.md'))))))


def source_bindings():
    inherited = parent.source_bindings()
    if strict_json.load(parent.OUTPUT)['source_bindings'] != inherited:
        raise ValueError('the weighted family source receipt is stale')
    inherited.update({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                      for name in SOURCE_FILES})
    return dict(sorted(inherited.items()))


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def rounded_interval(interval, places=24):
    """Small receipt presentation with exact outward decimal-grid rounding."""
    lower, upper = interval
    scale = 10**places
    lo = lower.numerator*scale//lower.denominator
    hi = -(-upper.numerator*scale//upper.denominator)
    return [str(F(lo, scale)), str(F(hi, scale))]


def finite_examples():
    rows = []
    cases = ((mean, cells, retained) for mean in (2,4)
             for cells, retained in ((1,1),(2,1),(5,2),(20,1)))
    for mean, cells, retained in cases:
        counts = limit.canonical_counts(cells, retained, mean=mean)
        if sum(counts.values()) != 1:
            raise ArithmeticError('the full canonical local law is not normalized')
        distance = limit.distance_interval(cells, retained, mean=mean)
        bound = limit.local_bound(cells, retained, mean=mean)
        if distance[1] > bound['half_trace_distance_upper']:
            raise ArithmeticError('finite state distance exceeds the proved envelope')
        condition = limit.conditioning_interval(cells, mean=mean)
        if condition[0] <= bound['conditioning_probability_lower']:
            raise ArithmeticError('the finite conditioning control contradicts the theorem')
        pair = limit.adjacent_count_control(cells, mean=mean)
        if pair['event_probability_lower'] <= pair['theorem_conditioning_lower']:
            raise ArithmeticError('the adjacent-count control contradicts the theorem')
        numbers = [sum((p for (a, b), p in counts.items() if a+b == n), F())
                   for n in range(12*retained+1)]
        rows.append({'mean_per_cell':mean, 'cells':cells, 'retained_cells':retained,
                     'all_count_classes':len(counts),
                     'full_local_number_distribution':[str(x) for x in numbers],
                     'full_state_half_trace_distance_interval':rounded_interval(distance),
                     'conditioning_probability_interval':rounded_interval(condition),
                     'conditioning_probability_lower':str(bound['conditioning_probability_lower']),
                     'half_trace_distance_upper':str(bound['half_trace_distance_upper']),
                     'adjacent_counts':[pair['high_count'], pair['low_count']]})
    return rows


def produce():
    bindings = source_bindings()
    word = limit.native_word()
    parameters = {str(mean):limit.parameters(mean=mean) for mean in (2,4)}
    # These exact inequalities are the small rational controls in the proof.
    if not F(5, 1134) > F(1, 240):
        raise ArithmeticError('the conditioning constant failed')
    if sum((F(7, 10)**j/factorial(j)
            for j in range(5)), F()) <= 2:
        raise ArithmeticError('the elementary logarithm envelope failed')
    capacities = []
    for mean, level in ((mean, level) for mean in (2,4) for level in (0,5,20,100,200,1000)):
        bound = limit.level_bound(level, mean=mean)
        capacities.append({key:(str(value) if type(value) is F else value)
                           for key, value in bound.items()})
    huge = limit.level_bound(200)
    if huge['ports'] <= 10**122 or huge['half_trace_distance_upper'] >= F(2, 10**60):
        raise ArithmeticError('the large-capacity certificate failed')
    result = {'schema':'oph.weighted_local_limit.v1',
              'family':{'repeated_native_word':list(word), 'cells_domain':'integer B >= 1',
                        'retained_cells_domain':'integer 0 <= k <= B',
                        'ports':'12*B', 'declared_total_occupation':'mean_per_cell*B',
                        'declared_mean_occupations_per_cell':[2,4],
                        'density_by_mean':{'2':'1/6','4':'1/3'},
                        'fugacity_by_mean':{'2':'(sqrt(1256)-16)/1000','4':'(1+sqrt(6401))/1600'},
                        'fugacity_polynomial_coefficients_by_mean':{'2':[1000,32,-1],'4':[800,-1,-2]},
                        'source_federation_cell_count':'20*4**level'},
              'algebraic_intervals_by_mean':{
                  mean:{key:rounded_interval(value,36) for key,value in values.items()}
                  for mean,values in parameters.items()},
              'all_domain_theorem':{
                  'canonical_is_product_conditioned_on_total':True,
                  'conditioning_probability_lower':'1/(240*B)',
                  'global_relative_entropy_upper':'log(240*B)',
                  'local_relative_entropy_upper':'k*log(240*B)/B',
                  'half_trace_distance_upper':'min(1,sqrt(k*log(240*B)/(2*B)))',
                  'arbitrary_retained_cell_subset':True,
                  'nondivisor_window_sizes_included':True,
                  'all_local_number_sectors_retained':True,
                  'restrictions_to_ports_within_retained_cells':True,
                  'two_capacity_comparison_by_common_product_reference':True,
                  'fixed_local_complete_quantum_instruments_and_independent_reference':True,
                  'unamplified_contractivity_requires_exact_CPTP':True,
                  'bounded_gain_implementation_error_upper':'G*epsilon+eta',
                  'rounded_primitive_maps_assumed_trace_preserving':False,
                  'successful_postselection_normalization_included':False,
                  'arbitrary_correlated_exterior_reference_access_included':False,
                  'all_domain_proof_medium':'prose using finite log-sum and exact elementary bounds'},
              'finite_examples':finite_examples(),
              'large_capacity_certificates':capacities,
              'scope':{'density_and_fugacity_are_declared':True,
                       'weights_from_declared_weighted_source_extension':True,
                       'native_weighted_federation_admitted':False,
                       'source_selected_capacity_or_density':False,
                       'canonical_to_product_limit_is_static':True,
                       'preparation_or_mixing_time_bound':False,
                       'autonomous_local_connected_dynamics_derived':False,
                       'huge_global_carrier_enumerated':False,
                       'physical_volume_or_cosmological_nats_identified':False,
                       'continuum_QFT_or_measured_mass_prediction':False},
              'source_bindings':bindings}
    if source_bindings() != bindings:
        raise ArithmeticError('weighted local-limit sources changed during replay')
    return result


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('weighted local-limit receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()), encoding='utf-8', newline='\n')
        print('Weighted local-limit receipt written.')
    else:
        verify()
        print('Weighted local-limit receipt verified.')
