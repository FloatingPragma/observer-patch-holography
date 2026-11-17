"""Exact finite reconstruction from closed native configuration histories.

Native integer branch counts define every moment before any matrix-state or
Born interpretation. The reconstructed algebra uses the complete-attempt
clock, not the distinct incident-clock similarity dictionary. Complexification
and interpretation of arbitrary reconstructed effects are explicit scope.
"""
from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer import native, quantum_mixing

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
STATES = native.states(binary=True)
INDEX = {state: i for i, state in enumerate(STATES)}
DIMENSION = len(STATES)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def index(value):
    require(type(value) is int and 0 <= value < DIMENSION, 'a native configuration index is required')
    return value


def length(value):
    require(type(value) is int and 1 <= value <= 16, 'a loop length from one through sixteen is required')
    return value


@lru_cache(maxsize=1)
def source_counts():
    """Input rows/output columns from every actual seam and placement coin."""
    result = [[0]*DIMENSION for _ in STATES]
    for i, state in enumerate(STATES):
        for edge in quantum_mixing.seams():
            for coin in (0, 1):
                final, recorded = native.native(state, edge, coin)
                require(recorded == (coin, state[edge[0]]-final[edge[0]]),
                        'source current differs from actual endpoint displacement')
                result[i][INDEX[final]] += 1
    denominator = 2*len(quantum_mixing.seams())
    require(DIMENSION == 66 and denominator == 60, 'the declared source carrier changed')
    require(all(sum(row) == denominator for row in result), 'source branch count is incomplete')
    require(all(result[i][j] == result[j][i] for i in range(DIMENSION)
                for j in range(DIMENSION)), 'complete-attempt source is not symmetric')
    return tuple(map(tuple, result)), denominator


@lru_cache(maxsize=16)
def path_counts(steps):
    """Count all length-L native paths; no trace probability function is called."""
    length(steps)
    counts, _ = source_counts()
    result = np.linalg.matrix_power(np.array(counts, dtype=object), steps)
    return tuple(tuple(int(value) for value in row) for row in result)


def source_spectral_bounds():
    counts, denominator = source_counts()
    lower = min(Q(2*counts[i][i]-denominator, denominator) for i in range(DIMENSION))
    require(lower == Q(2, 3), 'native diagonal-dominance certificate changed')
    # For symmetric nonnegative C, 2ab >= -a²-b² bounds its quadratic form
    # below by sum_i(2 C_ii-row_sum_i) v_i². Upper bound is the row norm.
    return {'lower': str(lower), 'upper': '1',
            'minimum_diagonal': str(min(Q(counts[i][i], denominator) for i in range(DIMENSION))),
            'strict_positive_definiteness_from_integer_counts': True}


@lru_cache(maxsize=1)
def matrix_unit_paths():
    """All matrix units from configuration projectors and connected native steps."""
    counts, denominator = source_counts()
    result = {}
    for origin in range(DIMENSION):
        paths = {origin: (origin,)}
        queue = deque([origin])
        while queue:
            current = queue.popleft()
            for target, count in enumerate(counts[current]):
                if count and target not in paths:
                    paths[target] = (*paths[current], target)
                    queue.append(target)
        require(len(paths) == DIMENSION, 'native configuration support is disconnected')
        for target, path in paths.items():
            coefficient = Q(1)
            for a, b in zip(path, path[1:]):
                coefficient *= Q(counts[a][b], denominator)
            require(coefficient > 0, 'a spanning matrix-unit path has zero weight')
            result[origin, target] = (path, coefficient)
    return result


def generated_word(path):
    """Exact sparse multiplication of P_x0 T P_x1 ... T P_xm.

    This verifies the operator word separately from the BFS coefficient.
    Configuration projectors are products of the twelve occupation projectors
    or their complements, checked in occupation_projector_audit().
    """
    require(type(path) is tuple and path, 'a nonempty configuration path is required')
    for item in path:
        index(item)
    counts, denominator = source_counts()
    row = [Q() for _ in STATES]
    row[path[0]] = Q(1)
    for target in path[1:]:
        value = sum(row[k]*Q(counts[k][target], denominator) for k in range(DIMENSION))
        row = [Q() for _ in STATES]
        row[target] = value
    return {(path[0], j): value for j, value in enumerate(row) if value}


def occupation_projector_audit():
    """All native configuration projectors from disclosed site occupations."""
    for i, state in enumerate(STATES):
        diagonal = tuple(int(all(candidate[p] == state[p] for p in range(12))) for candidate in STATES)
        require(diagonal == tuple(int(i == j) for j in range(DIMENSION)),
                'occupation fields do not separate native configurations')
    return {'site_fields': 12, 'configuration_projectors': DIMENSION,
            'all_projectors_generated_as_products_of_occupations_and_complements': True}


def closed_path_probability(steps, path):
    """A complete native path event, with L unconstrained return steps.

    Start uniformly, follow the displayed configuration path, then take L
    further native attempts and return to the initial configuration. Other
    outcomes are outside this event, not discarded from its probability.
    """
    length(steps)
    require(type(path) is tuple and path, 'a nonempty configuration path is required')
    for item in path:
        index(item)
    counts, denominator = source_counts()
    numerator = path_counts(steps)[path[-1]][path[0]]
    for a, b in zip(path, path[1:]):
        numerator *= counts[a][b]
    return Q(numerator, DIMENSION*denominator**(steps+len(path)-1))


def return_probability(steps):
    power = path_counts(steps)
    _, denominator = source_counts()
    return Q(sum(power[i][i] for i in range(DIMENSION)), DIMENSION*denominator**steps)


def matrix_unit_moment(steps, a, b):
    """Moment defined by scaled, independently normalized closed-history counts."""
    path, coefficient = matrix_unit_paths()[index(a), index(b)]
    return closed_path_probability(steps, path)/(return_probability(steps)*coefficient)


def reconstructed_density(steps):
    """Density representation derived from all matrix-unit count moments."""
    power = path_counts(steps)
    normalization = sum(power[i][i] for i in range(DIMENSION))
    density = tuple(tuple(Q(power[i][j], normalization) for j in range(DIMENSION))
                    for i in range(DIMENSION))
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            require(matrix_unit_moment(steps, a, b) == density[b][a],
                    'closed-history moment does not equal the derived density coordinate')
    return density


def sparse_moment(steps, entries):
    """Canonical complex-linear extension, with exact rational real/imaginary parts.

    Entries are (row,column):(real,imaginary). General signed combinations are
    reconstructed statistics; they are not declared single native events.
    """
    require(type(entries) is dict, 'a sparse exact matrix dictionary is required')
    real, imaginary = Q(), Q()
    for (a, b), pair in entries.items():
        index(a); index(b)
        require(type(pair) is tuple and len(pair) == 2 and
                all(type(value) in (int, Q) for value in pair),
                'rational complex coordinates are required')
        value = matrix_unit_moment(steps, a, b)
        real += Q(pair[0])*value
        imaginary += Q(pair[1])*value
    return real, imaginary


def example_effects(steps=1):
    """A non-diagonal real projector and its phase-conjugate counterparts."""
    a, b = next((a, b) for (a, b), (path, _) in matrix_unit_paths().items()
                if a != b and len(path) == 2)
    plus = {(a,a):(Q(1,2),0), (b,b):(Q(1,2),0),
            (a,b):(Q(1,2),0), (b,a):(Q(1,2),0)}
    minus = {key: ((-real if key[0] != key[1] else real), imag)
             for key,(real,imag) in plus.items()}
    yplus = {(a,a):(Q(1,2),0), (b,b):(Q(1,2),0),
             (a,b):(0,Q(-1,2)), (b,a):(0,Q(1,2))}
    yminus = {key:(real,-imag) for key,(real,imag) in yplus.items()}
    values = {name:sparse_moment(steps, effect) for name,effect in
              [('plus',plus),('minus',minus),('y_plus',yplus),('y_minus',yminus)]}
    require(all(imag == 0 and 0 <= real <= 1 for real,imag in values.values()),
            'reconstructed projection weights are outside probability range')
    require(values['plus'][0] > values['minus'][0] and values['y_plus'] == values['y_minus'],
            'real coherence or conjugation boundary changed')
    return {'indices':[a,b], 'weights':{name:str(value[0]) for name,value in values.items()},
            'real_reconstructed_coherence_resolved': True,
            'complex_conjugation_orientation_selected': False,
            'these_effects_implemented_as_native_single_trial_readouts': False}


def configuration_attachment():
    """Match the native current-only observer and include moment amplification.

    The adapter executes the same complete-attempt source. Its success flag,
    full configuration word and endpoint comparison make every filtered loop
    event in the reconstruction record-readable. The finite allowance is for
    apparatus/preparation error, not for statistical estimation from samples.
    """
    from ..native_thermal_return import certificate as adapter
    counts, _ = source_counts()
    require(tuple(adapter.BASIS) == tuple(STATES), 'native census uses a different configuration chart')
    require(adapter.attempt_rows() == tuple(tuple((j,value) for j,value in enumerate(row) if value)
                                            for row in counts),
            'current-record observer executes a different attempt transfer')
    budget = adapter.bounds(clock='attempt', length=8, seed_scans=8, retries=64,
                            warmup=1024, absorb=True)
    delta = Q(741,10**20)
    require(budget['complete_law_error_upper'] <= delta,
            'census/preparation error exceeds the reconstruction allowance')
    coefficient = min(value for _,value in matrix_unit_paths().values())
    require(coefficient == Q(1,60**4), 'matrix-unit reconstruction amplification changed')
    raw_bound = delta*(1/coefficient+1)/(Q(1,66)-delta)
    coordinate = Q(634,10**11)
    require(raw_bound <= coordinate, 'ratio reconstruction error exceeds its outward bound')
    # For all complex effects ||E||_F <= sqrt(66) < 9, while an entrywise
    # coordinate error <=epsilon gives ||Delta rho||_F <=66 epsilon.
    return {'clock':'complete native attempt', 'census_seed_scans':8, 'census_retries':64,
            'discarded_warmup_attempts':1024, 'initial_collision_sweep':True,
            'maximum_body_attempts':8, 'native_attempts_per_experiment_max':budget['native_attempts_max'],
            'complete_law_error_upper':str(delta), 'minimum_word_coefficient':str(coefficient),
            'return_probability_lower':'1/66', 'coordinate_error_upper':str(coordinate),
            'all_complex_effect_statistic_error_upper':str(66*9*coordinate),
            'finite_sampling_error_included':False, 'raw_reconstructed_matrix_certified_positive':False,
            'source_census_controller_executed_by_this_packet':False,
            'actual_source_transfer_and_controller_bound_compared':True}


def report():
    counts, denominator = source_counts()
    paths = matrix_unit_paths()
    for (a,b),(path,coefficient) in paths.items():
        require(generated_word(path) == {(a,b):coefficient}, 'matrix-unit generation failed')
    levels = []
    for steps in (1,2,4):
        density = reconstructed_density(steps)
        require(sum(density[i][i] for i in range(DIMENSION)) == 1, 'density is not normalized')
        levels.append({'native_attempts':steps, 'return_probability':str(return_probability(steps)),
                       'different_probability':str(1-return_probability(steps)),
                       'matrix_unit_moments_verified':DIMENSION**2,
                       'density_sha256':hashlib.sha256(json.dumps([[str(x) for x in row] for row in density]).encode()).hexdigest(),
                       'reconstructed_effects':example_effects(steps)})
    sources = ['model/native_quantum_reconstruction/certificate.py',
               'model/native_quantum_reconstruction/PROOF.md',
               'model/native_quantum_reconstruction/test_certificate.py',
               'model/observer/native.py','model/observer/quantum_mixing.py',
               'model/common_source/slow_geometry.py',
               'model/native_thermal_return/certificate.py',
               'model/native_clock/configuration.py',
               'model/record_laws/burn_in.py',
               'evidence/native_geometric_source_20260925/spec.json',
               'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
               'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
               'Lean/EventAlgebra/FiniteBuschGleason.lean']
    return {'schema':'oph.native_quantum_reconstruction.v1','all_pass':True,
            'source':{'states':DIMENSION,'seams':denominator//2,'actual_coin_branches':DIMENSION*denominator,
                      'clock':'every native seam attempt, including zero current',
                      'integer_counts_sha256':hashlib.sha256(json.dumps(counts).encode()).hexdigest()},
            'spectral_bounds':source_spectral_bounds(),
            'generation':{**occupation_projector_audit(),'matrix_units':len(paths),
                          'maximum_token_path_length':max(len(path)-1 for path,_ in paths.values()),
                          'full_real_matrix_algebra_generated':True},
            'closed_history_levels':levels,
            'current_only_configuration_attachment':configuration_attachment(),
            'scope':{'native_count_moment_functional_derived':True,
                     'positive_reconstructed_density_derived':True,
                     'born_form_on_reconstructed_complexification_derived':True,
                     'canonical_complexification_forced_as_physical_scalar_field':False,
                     'arbitrary_complex_effect_native_implementation':False,
                     'native_event_clock_identified_with_unitary_real_time':False,
                     'rho_L_prepared_on_native_quantum_working_carrier':False,
                     'complete_configuration_service_executed_by_this_packet':False,
                     'native_source_fields_same_as_incident_clock_dictionary':False},
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    value=report()
    path=HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    if args.check:
        require(strict_json.load(path)==json.loads(json.dumps(value)), 'reconstruction receipt is stale')
    print(json.dumps({'all_pass':True,'matrix_units':value['generation']['matrix_units'],
                      'native_general_effect_implementation':False}))


if __name__=='__main__':
    main()
