"""Finite selective population realization from at most five native records.

The source, stationary population law and complete raw alphabet are fixed.
The recovered transfer acts on population functions, not native quantum kets.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json

import sympy as s

from model.native_record_response import response as source, transfer

DIMENSION = 66
RAW_SLOTS = 721
SCORE_BOUND = 3600
GRAM_FLOOR = F(55756800, 2129251819009)
SOURCE_FILES = (
    'model/native_record_response/short_records.py',
    'model/native_record_response/SHORT_RECORD_RECONSTRUCTION.md',
    'model/native_record_response/test_short_records.py',
    'Lean/SourceLaw/RecordRealization.lean',
)


def _norm1(matrix):
    return max(sum(abs(v) for v in matrix[:, j]) for j in range(matrix.cols))


def _norminf(matrix):
    return max(sum(abs(v) for v in matrix[i, :]) for i in range(matrix.rows))


@lru_cache(maxsize=1)
def _model():
    d = source.data()
    source.validate(d)
    rows = [s.ones(1, DIMENSION)]
    labels = [('constant',)]
    for e in range(30):
        rows.append(d['D'][e, :]); labels.append(('single', e))
    for e in range(30):
        for f in range(30):
            # Chronological Y_e followed by Y_f, including first-read disturbance.
            rows.append(d['D'][f, :]*d['B'][e]); labels.append(('pair', e, f))
    candidates = s.Matrix.vstack(*rows)
    selected = source.independent_rows_modular(source._modular(candidates))
    if len(selected) != DIMENSION:
        raise ArithmeticError('two-record functions fail to span the actual population carrier')
    X = s.ImmutableMatrix(candidates[list(selected), :].T)
    A = s.ImmutableMatrix(X.to_DM().inv().to_Matrix())
    Pi = s.ImmutableMatrix(s.diag(*d['pi']))
    if A*X != s.eye(DIMENSION) or X*A != s.eye(DIMENSION):
        raise ArithmeticError('exact record frame inverse failed')
    G = s.ImmutableMatrix(X.T*Pi*X)
    Gi = s.ImmutableMatrix(A*s.diag(*(1/p for p in d['pi']))*A.T)
    gamma = F(1/_norminf(Gi))
    if gamma != GRAM_FLOOR:
        raise ArithmeticError('actual reflected Gram floor changed')
    chosen = tuple(labels[i] for i in selected)
    if (sum(x[0] == 'single' for x in chosen), sum(x[0] == 'pair' for x in chosen)) != (11, 54):
        raise ArithmeticError('short-record basis types changed')
    e0 = s.eye(DIMENSION)[:, 0]
    if X*e0 != s.ones(DIMENSION, 1) or (e0.T*G*e0)[0] != 1:
        raise ArithmeticError('the reconstructed ground is not the source constant')
    return d, chosen, X, A, Pi, G, Gi


def feature_labels():
    return _model()[1]


def score(feature, raw_word):
    """Bounded scalar on every complete two-slot record, including zero slots."""
    d = _model()[0]
    indexed = {transfer.label(row): row for row in d['branches']}
    if (type(raw_word) is not tuple or len(raw_word) != 2
            or any(value not in indexed for value in raw_word)
            or feature not in feature_labels()):
        raise ValueError('one admitted feature and two complete raw labels required')
    if feature[0] == 'constant':
        return 1
    first = indexed[raw_word[0]].feature(source.local.SEAMS[feature[1]])
    if feature[0] == 'single':
        return first
    return first*indexed[raw_word[1]].feature(source.local.SEAMS[feature[2]])


def one_record_rank():
    """Rank of the constant and EVERY raw one-record indicator, including zeros."""
    rows = [s.ones(1, DIMENSION)]
    for branch in _model()[0]['branches']:
        row = s.zeros(1, DIMENSION)
        for i, _ in branch.pairs:
            row[i] = branch.weight
        rows.append(row)
    matrix = s.Matrix.vstack(*rows)
    exact = matrix.to_DM().rank()
    modular = len(source.independent_rows_modular(source._modular(matrix)))
    if exact != modular or exact != 42:
        raise ArithmeticError('the complete one-record rank changed')
    return exact


def _branch(label):
    if type(label) is not tuple:
        raise ValueError('one complete raw label required')
    found = [row for row in _model()[0]['branches'] if transfer.label(row) == label]
    if len(found) != 1:
        raise ValueError('unknown complete native raw label')
    return found[0]


def _backward(row):
    return s.ImmutableSparseMatrix(DIMENSION, DIMENSION,
        {(i, j): s.Rational(row.weight) for i, j in row.pairs})


def moment_matrix(label):
    """Five-event reflected moment X.T Pi L_h X, without dropping raw zeros."""
    _, _, X, _, Pi, _, _ = _model()
    L = _backward(_branch(label))
    return s.ImmutableMatrix((X.T.to_DM()*Pi.to_DM()*L.to_DM()*X.to_DM()).to_Matrix())


def gram(gap=0):
    if type(gap) is not int or gap not in (0, 1):
        raise ValueError('zero or one unread bridge attempt required')
    d, _, X, _, Pi, G, _ = _model()
    return G if gap == 0 else s.ImmutableMatrix(X.T*Pi*d['P'].T*X)


def recover_insertion(moment):
    """Recover the coordinate insertion from exact calibrated moment entries."""
    if (not isinstance(moment, s.MatrixBase) or moment.shape != (DIMENSION, DIMENSION)
            or any(value.is_Rational is not True for value in moment)):
        raise ValueError('an exact rational 66 by 66 moment matrix required')
    Gi = _model()[-1]
    return s.ImmutableMatrix((Gi.to_DM()*s.ImmutableMatrix(moment).to_DM()).to_Matrix())


def word_probability(word):
    """Held-out chronological raw word using recovered insertions only."""
    if type(word) is not tuple:
        raise ValueError('a finite immutable raw word required')
    G = gram(); value = s.eye(DIMENSION)[:, 0]
    for label in reversed(word):
        value = recover_insertion(moment_matrix(label))*value
    return F((G[0, :]*value)[0])


def _positive_words():
    """Fixed commands from a specified input define admissible long raw words."""
    d = _model()[0]
    specifications = (
        (('E', 0, 1), ('C', 1, 0), ('E', 7, 1), ('C', 11, 0), ('E', 20, 0), ('C', 29, 1)),
        (('C', 11, 1), ('E', 19, 0), ('C', 3, 1), ('E', 0, 0), ('C', 17, 0), ('E', 22, 1), ('C', 5, 1)),
    )
    words = []
    for directives in specifications:
        at = 0; word = []
        for law, e, coin in directives:
            available = [(row, dict(row.pairs)[at]) for row in d['branches']
                         if row.law == law and row.seam == source.local.SEAMS[e]
                         and row.coin == coin and at in dict(row.pairs)]
            if len(available) != 1:
                raise ArithmeticError('literal command lost its unique raw source output')
            row, at = available[0]
            word.append(transfer.label(row))
        words.append(tuple(word))
    return tuple(words)


def coarse_history(steps=7):
    """All 2**steps original sign-bit words, recovered after five-event calibration."""
    if type(steps) is not int or not 0 <= steps <= 8:
        raise ValueError('an integer complete-history length zero through eight required')
    d, _, X, _, Pi, G, _ = _model()
    L = [s.zeros(DIMENSION), s.zeros(DIMENSION)]
    for row in d['branches']:
        bit = int(row.current is not None and row.current > 0)
        L[bit] += _backward(row)
    recovered = [recover_insertion(s.ImmutableMatrix(X.T*Pi*op*X)) for op in L]
    e0 = s.eye(DIMENSION)[:, 0]
    states = {(): (e0, s.ones(DIMENSION, 1))}
    # Prepend a new event to the future word: chronological matrix order is retained.
    for _ in range(steps):
        states = {(bit,)+word: (recovered[bit]*v, L[bit]*original)
                  for word, (v, original) in states.items() for bit in (0, 1)}
    result = {}
    for word, (v, original) in sorted(states.items()):
        predicted, actual = F((G[0, :]*v)[0]), F((d['pi'].T*original)[0])
        if predicted != actual or predicted < 0:
            raise ArithmeticError('complete held-out sign history differs from the native recursion')
        result[''.join(map(str, word))] = str(predicted)
    if sum(map(F, result.values()), F()) != 1:
        raise ArithmeticError('held-out complete history lost original probability')
    return result


def noise_budget(entry_error, attempts, *, preparation_error=F(0)):
    """Complete held-out TV bound for a feasible population-kernel inference.

    Each candidate and true insertion moment lies in the same +/-epsilon
    interval. The candidate must be nonnegative and complete in the specified
    population coordinates, with the same configured pi and calibrated X.
    This function does not alter an observed record law or fit a quantum map.
    """
    if (type(entry_error) not in (int, F) or entry_error < 0
            or type(preparation_error) not in (int, F) or not 0 <= preparation_error <= 1
            or type(attempts) is not int or not 0 <= attempts <= 2**32):
        raise ValueError('exact nonnegative errors and a bounded integer horizon required')
    _, _, _, A, Pi, _, _ = _model()
    factor = F(RAW_SLOTS*DIMENSION*_norminf(A)*_norm1(A*Pi.inv()))
    return {'moment_entry_error': str(entry_error), 'attempts': attempts,
            'complete_kernel_TV_lipschitz_factor': str(factor),
            'complete_history_TV_upper': str(min(F(1), F(preparation_error)+attempts*factor*entry_error)),
            'requires_nonnegative_complete_candidate_kernel': True,
            'requires_same_calibrated_population_frame_and_pi': True,
            'full_quantum_or_adaptive_policy_bound': False}


def transfer_energy_error(entry_error):
    """Uniform generalized-eigenvalue/log error for four/five-record Gram data."""
    if type(entry_error) not in (int, F) or entry_error < 0:
        raise ValueError('a nonnegative exact entry allowance required')
    delta = DIMENSION*F(entry_error)
    if delta >= GRAM_FLOOR:
        raise ValueError('entry allowance does not certify the positive Gram')
    eigen = 2*delta/(GRAM_FLOOR-delta)
    lower = transfer.LOWER-eigen
    if lower <= 0:
        raise ValueError('entry allowance does not certify a positive transfer')
    return {'Gram_operator_error_upper': str(delta),
            'transfer_eigenvalue_error_upper': str(eigen),
            'transfer_log_energy_error_upper': str(eigen/lower),
            'source_unit': 'native attempt transfer logarithm',
            'particle_Hamiltonian_energy': False}


def sample_budget(entry_error, failure_probability):
    """Conservative independent-trial count, never an executed sampling claim."""
    if (type(entry_error) not in (int, F) or entry_error <= 0
            or type(failure_probability) not in (int, F) or not 0 < failure_probability < 1):
        raise ValueError('positive exact precision and failure probability below one required')
    quantities = DIMENSION**2*(RAW_SLOTS+2)
    target = F(2*quantities)/failure_probability
    log_upper = 0
    while 2**log_upper < target:
        log_upper += 1
    bound = 2*SCORE_BOUND**4*log_upper/F(entry_error)**2
    count = (bound.numerator+bound.denominator-1)//bound.denominator
    return {'independent_stationary_trials_per_menu_sufficient': count,
            'menus': ['four-record reflected Gram', 'five-record complete raw insertions and summed bridge'],
            'max_events_per_trial': 5, 'max_feature_magnitude': SCORE_BOUND,
            'simultaneous_scalar_quantities': quantities,
            'independent_stationary_reset_or_equivalent_sampling_required': True,
            'sampled_trials_executed': False}


def continuing_calibration_budget(entry_error=F(1, 2**40), failure_probability=F(1, 100)):
    """Native burn-in and martingale budget, with admission/evaluation limits."""
    base = sample_budget(entry_error, failure_probability)
    mixing = source.mixing
    p = source.parameters()
    minor = mixing.population_minorization(p, 20)
    rate = mixing.coherence_rate(p)
    if minor['beta'] < F(4, 25) or rate != F(793, 990):
        raise ArithmeticError('actual native reset constants differ')
    burnin = 8192
    exact = F(65, 2)*rate**burnin+F(21, 25)**(burnin//20)
    grid = 2**128
    upper = F(-(-exact.numerator*grid//exact.denominator), grid)
    bias = 2*SCORE_BOUND**2*upper
    if bias > entry_error/2:
        raise ValueError('the fixed native burn-in does not certify this precision')
    stochastic = sample_budget(F(entry_error)/2, failure_probability)
    attempts = burnin+5
    cap = source.interpreter.MAX_ATTEMPTS
    maximum = cap//attempts
    # Admission is actually executed for the maximal whole-trial finite prefix.
    program = source.interpreter.program(commands=[{'op': 'repair',
        'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': maximum*attempts}])
    source.interpreter.admit(program)
    count = stochastic['independent_stationary_trials_per_menu_sufficient']
    return {'burn_in_attempts_per_trial': burnin,
            'twenty_attempt_minorization_lower': '4/25',
            'coherence_contraction_rate': str(rate),
            'reduced_source_reset_half_trace_upper': str(upper),
            'score_conditional_mean_bias_upper': str(bias),
            'total_moment_entry_error': str(entry_error),
            'martingale_trials_per_menu_sufficient': count,
            'stationary_independent_resets_assumed': False,
            'burn_in_archive_retained_but_not_available_to_later_controls': True,
            'full_source_and_new_burn_in_archive_product_reset_claimed': False,
            'raw_attempts_per_five_record_trial': attempts,
            'raw_attempts_for_two_menus_upper': 2*attempts*count,
            'maximum_native_program_attempts': cap,
            'maximum_complete_trials_in_one_admitted_program': maximum,
            'maximal_whole_trial_prefix_admitted': True,
            'required_sampler_fits_one_admitted_program': 2*attempts*count <= cap,
            'explicit_evaluator_step_cap': source.interpreter.MAX_EVALUATED_STEPS,
            'explicit_evaluator_history_block_cap': source.interpreter.MAX_EVALUATED_BRANCHES,
            'sampling_experiment_executed': False,
            'longer_service_or_multiple_program_execution_established': False}


def resource_stop_control():
    """Actual complete interpreter stop, before any long-burn-in probe completes."""
    d = _model()[0]
    program = source.interpreter.program(commands=[{'op': 'repair',
        'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 8197}])
    density = {(i, 0, i, 0): F(p) for i, p in enumerate(d['pi'])}
    result = source.interpreter.execute(program, density)
    if (result.status != 'RESOURCE_STOP' or result.attempted_steps != 1
            or len(result.branches) != RAW_SLOTS
            or sum(source.interpreter.trace(branch.density) for branch in result.branches) != 1
            or any(not branch.halted for branch in result.branches)):
        raise ArithmeticError('the native complete resource-stop boundary changed')
    return {'status': result.status, 'attempts_executed': result.attempted_steps,
            'retained_terminal_blocks': len(result.branches), 'total_original_probability': '1',
            'completed_calibration_trials': 0, 'input_included_actual_stationary_source': True}


@lru_cache(maxsize=1)
def _certificate_json():
    d, features, X, A, Pi, G, Gi = _model()
    Xd, Ad, Pid = X.to_DM(), A.to_DM(), Pi.to_DM()
    summed = s.zeros(DIMENSION)
    checked, zero = 0, 0
    for row in d['branches']:
        L = _backward(row).to_DM().convert_to(Xd.domain)
        B = Xd.transpose()*Pid*L*Xd
        if Ad.transpose()*B*Ad != Pid*L:
            raise ArithmeticError('a complete native selective map failed exact recovery')
        summed += L.to_Matrix()
        checked += 1; zero += not row.pairs
    if checked != RAW_SLOTS or zero != 421 or summed != d['P'].T:
        raise ArithmeticError('complete insertion sum or zero domain changed')
    if G*Gi != s.eye(DIMENSION):
        raise ArithmeticError('exact reflected Gram inverse failed')
    empty = next(transfer.label(row) for row in d['branches'] if row.label is None and row.kind == 'repair')
    positive = _positive_words()
    words = (*positive, (empty,)+positive[0], positive[0]+positive[1][:2])
    held = []
    for word in words:
        predicted = word_probability(word)
        actual = F(sum(transfer.forward_word(d, word)))
        if predicted != actual:
            raise ArithmeticError('a held-out chronological word differs from native recursion')
        held.append({'word': [list(v) for v in word], 'probability': str(predicted)})
    matrix_digest = hashlib.sha256(json.dumps([[str(v) for v in row] for row in X.tolist()],
        separators=(',', ':')).encode()).hexdigest()
    return json.dumps({'schema': 'oph.short_native_population_realization.v1',
        'dimension': DIMENSION, 'complete_one_record_rank': one_record_rank(),
        'minimum_population_observability_depth': 2, 'complete_raw_slots': checked, 'zero_raw_slots_retained': zero,
        'feature_basis': [list(v) for v in features], 'feature_frame_sha256': matrix_digest,
        'constant_single_pair_counts': [1, 11, 54], 'maximum_feature_events': 2,
        'Gram_calibration_events': 4, 'selective_map_calibration_events': 5,
        'exact_all_selective_map_inverse_checks': checked,
        'Gram_lower_bound': str(GRAM_FLOOR), 'inverse_frame_1_norm': str(_norm1(A)),
        'inverse_frame_infinity_norm': str(_norminf(A)),
        'held_out_original_word_probabilities': held,
        'complete_seven_event_sign_law': coarse_history(),
        'transfer_eigenvalue_interval': [str(transfer.LOWER), '1'],
        'logarithm_degree_16_tail': str(transfer.logarithm_tail(16)),
        'noise_example': noise_budget(F(1, 2**60), 16),
        'transfer_energy_error_example': transfer_energy_error(F(1, 2**40)),
        'sampling_example': sample_budget(F(1, 2**40), F(1, 100)),
        'continuing_source_calibration': continuing_calibration_budget(),
        'actual_complete_resource_stop': resource_stop_control(),
        'scope': {'configured_finite_stationary_population_source': True,
                  'all_fixed_schedule_finite_raw_word_laws_recovered': True,
                  'arbitrary_hidden_dimension_identified_from_five_events': False,
                  'full_quantum_instrument_tomography': False,
                  'physical_clock_or_particle_energy_identified': False,
                  'coherent_CAR_gauge_Higgs_action_selected': False,
                  'empirical_sampling_executed': False,
                  'exact_H66_numbers_evaluated_in_Lean_kernel': False}}, sort_keys=True)


def certificate():
    """Detached public values; the cached representation contains JSON text only."""
    return json.loads(_certificate_json())
