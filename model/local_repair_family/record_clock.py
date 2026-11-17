"""Recover, prepare and enact one continuing reference from native records.

The preparation admission uses only the recorded coin intervals.  Its lower
population kernel and upper common-record survival bound cover every true
member in that simultaneous confidence set.  The source simulator alone
receives private port weights and its terminal state.  Training records can
be included in the untouched reference of the preparation theorem; the
subsequent preparation archive is excluded from the product comparison.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import local_bias as local, inference, mixing, observer, modular_clock

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CONTRACT = HERE / 'record_clock_campaign.json'
OUTPUT = HERE / 'record_clock_results.json'
CAMPAIGN = 'OPH-RECORD-DERIVED-CONTINUING-CLOCK-v1'
CONTRACT_SHA256 = 'b67b82ebfc34282a9281147ec45717de410e02a75610d32062153c2ff6be0fa3'
PREPARATION_ATTEMPTS = 1600
BLOCK = 20
SCALE = 10**12
PREPARATION_CAP = F(1, 100000)
COMPILER_CAP = F(1, 1000000)
GAIN_CAP = F(1000001, 1000000)


def _json(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {str(key): _json(item) for key, item in value.items()}
    if type(value) in (list, tuple):
        return [_json(item) for item in value]
    return value


def _digest(value):
    return hashlib.sha256(json.dumps(_json(value), sort_keys=True,
                                    separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def coin_effect_certificate():
    """Exact complete effect identity, including all input coherences.

    Each raw map is a partial injection, so its effect is a diagonal domain
    projection.  Those domains partition all66 inputs for a fixed coin.
    Thus sum_raw K* K = p_coin I, on arbitrary source/reference states.
    """
    columns = 0
    for name in ('C', 'E'):
        for edge in local.SEAMS:
            for coin in (0, 1):
                covered = [0]*66
                for c, _, pairs in local.raw_maps(name, edge):
                    if c != coin:
                        continue
                    if (len({i for i, _ in pairs}) != len(pairs)
                            or len({j for _, j in pairs}) != len(pairs)):
                        raise ArithmeticError('raw record map is not a partial injection')
                    for i, _ in pairs:
                        covered[i] += 1
                if covered != [1]*66:
                    raise ArithmeticError('complete raw effects do not sum to the coin identity')
                columns += len(covered)
    return {'complete_input_columns': columns,
            'identity': 'sum_raw K(seam,law,coin,raw)^* K(seam,law,coin,raw) = p_coin I66',
            'conditional_on_complete_past_and_arbitrary_quantum_reference': True,
            'source_population_or_state_tomography_used': False}


def faithfulness_certificate():
    """The prepared marginal is faithful even outside the confidence event.

    Eleven fixed E/tree/coin0 records separate all66 basis inputs. Each
    associated history is rank one. At most four positive native moves
    reach any target; an empty seam provides arbitrary positive padding.
    Select an input diagonal at least1/66. Its66 different terminal words
    contribute at least I/(66*300**1600) to every normalized input output.
    No tiny eigenvalue is converted to floating point.
    """
    words = []
    states = local.basis()
    for initial in states:
        state, word = initial, []
        for edge in inference.TREE:
            state, record = local.records.recorded_step('E', state, tuple(sorted(edge)), 0)
            word.append(record)
        words.append(tuple(word))
    if len(set(words)) != 66:
        raise ArithmeticError('the complete E/tree words do not separate all native inputs')
    adjacency = [set() for _ in states]
    for edge in local.SEAMS:
        for _, _, pairs in local.raw_maps('E', edge):
            for i, j in pairs:
                adjacency[i].add(j)
    diameter = 0
    for i, state in enumerate(states):
        if not any(state[a] == state[b] == 0 for a, b in local.SEAMS):
            raise ArithmeticError('positive empty-seam padding unavailable')
        visited, frontier, distance = {i}, {i}, 0
        while len(visited) < 66:
            frontier = set().union(*(adjacency[j] for j in frontier))-visited
            if not frontier:
                raise ArithmeticError('native population transport disconnected')
            visited.update(frontier); distance += 1
        diameter = max(diameter, distance)
    if diameter > 4 or len(inference.TREE)+diameter > PREPARATION_ATTEMPTS:
        raise ArithmeticError('positive rank-one histories exceed preparation horizon')
    return {'separating_E_word_length': len(inference.TREE), 'distinct_input_record_words': 66,
            'word_sha256': _digest(words), 'maximum_population_path': diameter,
            'per_native_branch_probability_lower': '1/300',
            'prepared_marginal_eigenvalue_lower': '1/(66*300^1600)',
            'all_training_histories_including_confidence_failures': True,
            'tiny_bound_used_as_a_floating_precision_claim': False}


def preparation_certificate(counts):
    """Certify1600 unread repairs from records, without a weight oracle.

    L<=P entrywise implies L^20<=P^20. Column minima of L^20 give a
    Doeblin minorization for the actual stochastic P, even though L is
    substochastic. For distinct input pairs, upper coin weights bound the
    probability that the complete next record is shared. The existing
    arbitrary-reference mixing argument then applies with those constants.
    """
    inferred = inference.infer_counts(counts)
    if inferred['status'] != 'accept':
        return {'status': 'stop', 'reason': 'unresolved_record_reference'}
    intervals = inferred['probability_intervals']
    lower = [[F() for _ in range(66)] for _ in range(66)]
    common = [[F() for _ in range(66)] for _ in range(66)]
    for name in ('C', 'E'):
        for edge in local.SEAMS:
            lo, hi = intervals[edge]
            for coin, _, pairs in local.raw_maps(name, edge):
                low, high = (lo, hi) if coin else (1-hi, 1-lo)
                low, high = low/60, high/60
                for i, j in pairs:
                    lower[i][j] += low
                    for k, _ in pairs:
                        common[i][k] += high
    if any(not 0 < sum(row) <= 1 for row in lower):
        raise ArithmeticError('record-derived lower kernel is not substochastic')
    denominator = math.lcm(*(x.denominator for row in lower for x in row))
    numerator = tuple(tuple(x.numerator*(denominator//x.denominator) for x in row)
                      for row in lower)
    powered = np.linalg.matrix_power(np.array(numerator, dtype=object), BLOCK)
    minima = tuple(min(int(powered[i, j]) for i in range(66)) for j in range(66))
    beta = F(sum(minima), denominator**BLOCK)
    rate = min(F(1), max(common[i][j] for i in range(66) for j in range(66) if i != j))
    # Compact outward rational bounds, with no floating spectral computation.
    beta = F(beta.numerator*SCALE//beta.denominator, SCALE)
    rate = F(-(-rate.numerator*SCALE//rate.denominator), SCALE)
    if not 0 < beta <= 1 or not 0 <= rate < 1:
        return {'status': 'stop', 'reason': 'record_box_has_no_mixing_certificate'}
    delta = F(65, 2)*rate**PREPARATION_ATTEMPTS + (1-beta)**(PREPARATION_ATTEMPTS//BLOCK)
    delta = F(-(-delta.numerator*SCALE//delta.denominator), SCALE)
    return {'status': 'accept', 'attempts': PREPARATION_ATTEMPTS,
            'law': 'equal tagged C/E mixture; uniform independent native seams',
            'population_block': BLOCK, 'population_minorization_lower': str(beta),
            'common_record_survival_upper': str(rate), 'half_diamond_upper': str(min(F(1), delta)),
            'lower_kernel_sha256': _digest((numerator, denominator)),
            'powered_column_minima_sha256': _digest((minima, denominator**BLOCK)),
            'private_port_weights_used': False,
            'simultaneous_confidence_event_required': True,
            'training_archive_in_preexisting_reference': True,
            'subsequent_preparation_archive_in_product_target': False}


def continued_source_density(terminal_index):
    """Exact population continuation of the basis trajectory, rounded once.

    The input and its unread native marginals are diagonal. Integer matrix
    products preserve the actual source probabilities through all1600 steps;
    no normalization of a rounded evolving state is used to conceal drift.
    The arbitrary quantum/reference guarantee comes from the separate proof.
    """
    if type(terminal_index) is not int or not 0 <= terminal_index < 66:
        raise ValueError('an original terminal native-state index required')
    kernel = mixing.probability_kernel(modular_clock.parameters())
    denominator = math.lcm(*(x.denominator for row in kernel for x in row))
    matrix = np.array([[x.numerator*(denominator//x.denominator) for x in row]
                       for row in kernel], dtype=object)
    population = np.zeros(66, dtype=object)
    population[terminal_index] = 1
    for _ in range(PREPARATION_ATTEMPTS):
        population = population@matrix
    scale = denominator**PREPARATION_ATTEMPTS
    if sum(population) != scale or any(x <= 0 for x in population):
        raise ArithmeticError('exact prepared population lost mass or faithfulness')
    exact = tuple(F(int(x), scale) for x in population)
    rounded = tuple(float(x) for x in exact)
    if sum(abs(F(x)-p) for x, p in zip(rounded, exact, strict=True))/2 > F(1, 10**15):
        raise ArithmeticError('prepared-state rounding exceeds its fixed diagnostic gate')
    density = np.diag(rounded).astype(complex)
    if (not np.isfinite(density).all() or abs(np.trace(density)-1) > 1e-10
            or np.max(np.abs(density-density.conj().T)) > 1e-10):
        raise ArithmeticError('continued source preparation lost finite probability mass')
    return density


def composition_budget(counts, source_budget):
    """Compose record-uniform bounds with an owned source certificate.

    The caller supplies source_budget only from fresh learned_flow execution;
    report/verify own that execution. This arithmetic helper authenticates
    neither a supplied budget nor a training transcript.
    """
    inferred = inference.infer_counts(counts)
    preparation = preparation_certificate(counts)
    if inferred['status'] != 'accept' or preparation['status'] != 'accept':
        return {'status': 'stop', 'reason': 'unresolved_reference_or_preparation'}
    def bound(key):
        value = source_budget.get(key) if type(source_budget) is dict else None
        if type(value) is not str:
            raise ValueError('canonical rational source bounds required')
        number = F(value)
        if str(number) != value or number < 0:
            raise ValueError('canonical nonnegative source bounds required')
        return number
    compiler = bound('complete_instrument_halfdiamond_upper')
    gain = bound('channel_norm_upper')
    if compiler > 1 or not 1 <= gain <= 2:
        raise ValueError('source compiler error or norm gain outside the fixed domain')
    delta = F(preparation['half_diamond_upper'])
    phase = F(3, 32)
    # Filled from the record-only design, independent of accepted transcripts.
    statistical = inference.fixed_budget()['statistical_failure_upper']
    on_event = gain*delta+compiler+phase
    # Admission bounds every possible accepted transcript by the SAME caps.
    # A data-dependent error cannot serve as the expectation over training.
    fixed_good = GAIN_CAP*PREPARATION_CAP+COMPILER_CAP+phase
    complete = fixed_good+statistical*(GAIN_CAP+1)/2
    # Uniform ratio<=4 on12 ports gives pi_min >=1/(66*16).
    # A trace-zero density perturbation of half trace norm delta has operator
    # norm<=delta. The sandwich/log inequality then compares the actual
    # prepared source-marginal modular group to the true stationary one.
    sandwich = 1056*PREPARATION_CAP
    actual_flow = sandwich/(1-sandwich)
    actual_complete = COMPILER_CAP+phase+actual_flow+statistical*(GAIN_CAP+1)/2
    accept = delta <= PREPARATION_CAP and compiler <= COMPILER_CAP and gain <= GAIN_CAP
    if max(complete, actual_complete) > F(1, 8):
        raise ArithmeticError('fixed average protocol cap exceeds its declared design')
    return {'status': 'accept' if accept else 'stop',
            'reason': 'complete_error_cap_satisfied' if accept else 'component_error_cap_exceeded',
            'preparation': preparation, 'preparation_half_trace_upper': str(delta),
            'inferred_modular_halfdiamond_upper': str(phase),
            'compiler_halfdiamond_upper': str(compiler), 'source_channel_norm_upper': str(gain),
            'simultaneous_confidence_failure_upper': str(statistical),
            'complete_on_confidence_event_half_trace_upper': str(on_event),
            'fixed_preparation_half_trace_cap': str(PREPARATION_CAP),
            'fixed_compiler_halfdiamond_cap': str(COMPILER_CAP),
            'fixed_channel_norm_cap': str(GAIN_CAP),
            'fixed_good_event_half_trace_cap': str(fixed_good),
            'complete_joint_protocol_half_trace_upper': str(complete),
            'uniform_stationary_minimum_eigenvalue_lower': '1/1056',
            'actual_prepared_marginal_relative_operator_error_cap': str(sandwich),
            'actual_prepared_marginal_modular_halfdiamond_cap': str(actual_flow),
            'complete_actual_prepared_marginal_clock_half_trace_upper': str(actual_complete),
            'actual_marginal_target': 'actual prepared joint source/R and training tape, with the modular group of its conditional source marginal; no replacement of that joint state',
            'actual_marginal_comparison_is_statewise_not_a_nonlinear_quantum_channel': True,
            'fixed_complete_error_cap': '1/8', 'conditioned_on_acceptance': False,
            'target': 'same training tape and untouched R; true product reference followed by its true modular clock and complete future record instrument',
            'all_preparation_records_including_final_FIFO_traced_in_comparison': True,
            'physical_apparatus_errors_included': False}


def contract():
    if hashlib.sha256(CONTRACT.read_bytes()).hexdigest() != CONTRACT_SHA256:
        raise ValueError('the fixed record-clock protocol contract changed')
    value = strict_json.load(CONTRACT)
    if (inference.TRIALS_PER_SEAM != value['training']['samples_per_seam']
            or PREPARATION_ATTEMPTS != value['continuation']['preparation_attempts']
            or inference.fixed_budget()['probability_radius'] != F(value['training']['confidence_radius'])
            or PREPARATION_CAP != F(value['acceptance']['preparation_half_trace_cap'])
            or COMPILER_CAP != F(value['acceptance']['compiler_halfdiamond_cap'])
            or GAIN_CAP != F(value['acceptance']['source_channel_norm_cap'])):
        raise ArithmeticError('execution and fixed record-clock design disagree')
    observer.contract()
    return value


def admission(counts):
    """Record-only decision, with compiler failures retained as explicit stops."""
    from . import learned_flow
    inferred = inference.infer_counts(counts)
    if inferred['status'] != 'accept':
        return {'status': 'stop', 'reason': inferred['stop_reason'], 'inference': inferred}
    try:
        source = learned_flow.source_budget(inferred['relative_weights'])
    except ArithmeticError:
        return {'status': 'stop', 'reason': 'source_certification_failed', 'inference': inferred}
    budget = composition_budget(counts, source)
    return {**budget, 'inference': inferred, 'source_budget': source}


def continue_from_records(counts, terminal_index, name):
    """Keep the actual terminal carrier on each stop; prepare only on admission.

    terminal_index is private simulator state, never an inference or source
    program argument.  The mathematical channel bound covers arbitrary
    terminal quantum inputs correlated with the training transcript and R.
    """
    from . import learned_flow
    name = local.records.candidate(name)
    if type(terminal_index) is not int or not 0 <= terminal_index < 66:
        raise ValueError('an original terminal native-state index required')
    decision = admission(counts)
    if decision['status'] != 'accept':
        rho = np.zeros((66, 66), complex)
        rho[terminal_index, terminal_index] = 1
        return decision, {('stop', decision['reason']): rho}
    density = continued_source_density(terminal_index)
    return decision, learned_flow.execute(density, decision['inference']['relative_weights'], name)


def _source_bindings(*providers):
    bindings = {}
    for provider in providers:
        for name, digest in provider.items():
            if name in bindings and bindings[name] != digest:
                raise ArithmeticError('record-clock components disagree on source custody')
            bindings[name] = digest
    for name in ('record_clock.py', 'test_record_clock.py', 'record_clock_campaign.json',
                 'README.md', 'observer.py', 'observer_campaign.json', 'continuing.py'):
        path = 'model/local_repair_family/' + name
        bindings[path] = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    return dict(sorted(bindings.items()))


def report():
    from . import acquisition, learned_flow
    frozen = contract()
    source, counts = acquisition.acquire(modular_clock.parameters(), seed=frozen['training']['seed'],
                                        initial_state=local.basis()[frozen['training']['initial_state_index']])
    decision = admission(counts)
    if decision['status'] != 'accept':
        raise ArithmeticError('the retained full native tape does not meet the complete clock cap')
    weights = decision['inference']['relative_weights']
    flow = learned_flow.report(weights)
    # Fresh preparation and clock execution continue the terminal carrier.
    terminal = source['private_replay']['terminal_state_index']
    density = continued_source_density(terminal)
    true_reference = np.diag(list(map(float, local.reference(modular_clock.parameters())))).astype(complex)
    delta = F(decision['preparation_half_trace_upper'])
    if np.linalg.norm(density-true_reference, 'nuc')/2 > float(delta):
        raise ArithmeticError('continued preparation exceeds its record-only certificate')
    rows = {}
    for name in ('C', 'E'):
        actual = learned_flow.execute(density, weights, name)
        ideal_learned = learned_flow.execute(density, weights, name, raw=False)
        true_target = modular_clock.execute(true_reference, name, raw=False)
        if actual.keys() != ideal_learned.keys() or actual.keys() != true_target.keys():
            raise ArithmeticError('record-derived and true-clock complete alphabets disagree')
        residuals = [sum(np.linalg.norm(actual[k]-target[k], 'nuc') for k in actual)/2
                     for target in (ideal_learned, true_target)]
        mass = sum(float(np.trace(block).real) for block in actual.values())
        compiler = F(decision['compiler_halfdiamond_upper'])
        complete = F(decision['complete_on_confidence_event_half_trace_upper'])
        if (not all(math.isfinite(x) for x in (*residuals, mass))
                or residuals[0] > float(compiler) or residuals[1] > float(complete)
                or abs(mass-1) > float(2*compiler)):
            raise ArithmeticError('continued complete output exceeds its owned source/readout gates')
        rows[name] = {'complete_output_blocks': len(actual),
                      'source_vs_inferred_ideal_half_trace_upper': str(compiler),
                      'source_vs_true_reference_clock_half_trace_upper': str(complete),
                      'complete_output_mass_error_upper': str(2*compiler),
                      'complete_quantum_outputs_checked': True,
                      'all_export_aborts_and_zero_raw_slots_retained': True}
    budget = flow['source_budget']
    attempts = source['protocol']['total_events']+PREPARATION_ATTEMPTS+3
    native_costs = acquisition.native_word_costs()
    native_cost_values = [value for row in native_costs.values() for value in row.values()]
    preparation_pulses = PREPARATION_ATTEMPTS*max(native_cost_values)
    preparation_expectation = F(PREPARATION_ATTEMPTS, 60)*sum(native_cost_values)
    future_costs = acquisition.native_word_costs((*modular_clock.REQUESTED_SEAMS, modular_clock.PROBE))
    future_pulses = max(sum(row[name] for row in future_costs.values()) for name in ('C', 'E'))
    compiled_total = (source['resources']['compiled_word_costs']['native_pulses']
                      + preparation_pulses+future_pulses+budget['maximum_source_pulses'])
    return _json({'schema': 'oph.local_repair_family.record_clock.v1', 'campaign_id': CAMPAIGN,
                 'contract_sha256': CONTRACT_SHA256,
                 'source_sha256': _source_bindings(source['source_sha256'], inference.source_bindings(),
                                                  flow['source_sha256'], mixing.source_bindings()),
                 'coin_effects': coin_effect_certificate(), 'unconditional_prepared_faithfulness': faithfulness_certificate(),
                 'identification': inference.exact_certificate(),
                 'acquisition': {k: v for k, v in source.items() if k != 'source_sha256'},
                 'admission': decision,
                 'compiled_clock': {k: v for k, v in flow.items() if k != 'source_sha256'},
                 'continued_complete_readouts': rows,
                 'tagged_half_C_half_E_readout_covered_by_convexity': True,
                 'resources': {'total_native_attempts': attempts,
                               'native_record_encoding_bits': 10,
                               'one_slot_fifo_bits': 10,
                               'exported_native_record_bits_upper': 10*(attempts-1),
                               'all_committed_native_record_bits_upper': 10*attempts,
                               'additional_clock_status_and_failure_metadata': 'separately counted in compiled_clock.source_budget',
                               'training_count_storage_bits': 60*(inference.TRIALS_PER_SEAM.bit_length()),
                               'training_event_counter_bits': source['protocol']['total_events'].bit_length(),
                               'estimated_rational_weights_storage_bits': sum(w.numerator.bit_length()+w.denominator.bit_length()
                                                                             for w in weights),
                               'coherent_clock_pulses': budget['maximum_source_pulses'],
                               'coherent_clock_load_export_calls': budget['coherent_load_export_calls'],
                               'clock_program_bytes': budget['source_program_bytes'],
                               'preparation_compiled_native_pulses_upper': preparation_pulses,
                               'preparation_compiled_native_pulses_expectation': str(preparation_expectation),
                               'future_native_word_pulses_upper': future_pulses,
                               'complete_supplied_realization_pulses_upper': compiled_total,
                               'complete_supplied_realization_load_export_calls_upper': 2*compiled_total,
                               'compiled_native_cost_scope': 'literal supplied word lengths; the31-million-event native replay executes exact transition tables, not trillions of coherent pulses',
                               'decoder_exact_arithmetic_and_workspace_are_declared_classical_services': True,
                               'archive_reconstruction_hash_is_not_a_replacement_for_physical_record_storage': True},
                 'scope': {'observer_supplied_port_weight_table': False,
                           'source_port_weights_selected_or_unique': False,
                           'coherent_actuator_derived_from_passive_repair': False,
                           'physical_seconds_or_energy_normalization': False,
                           'integrated_gate_accepted': False, 'physical_comparison_count_added': 0},
                 'remaining': frozen['remaining']})


def _compare(path, expected):
    actual = strict_json.load(Path(path))
    if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False):
        raise ValueError('record-clock receipt differs from fresh complete source replay')


def verify(path=OUTPUT):
    _compare(path, report())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    value = report()
    if args.check is not None:
        _compare(args.check, value)
    elif args.write:
        OUTPUT.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n',
                          encoding='utf8', newline='\n')
    else:
        print(json.dumps(value, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
