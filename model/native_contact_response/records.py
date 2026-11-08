"""Literal native repair records attached to the finite contact response."""
from collections import defaultdict
from fractions import Fraction as F

from model.native_correlation_preparation import preparation
from model.native_program import interpreter
from . import contact as c


def program(lag=0, warmup=0):
    c.integer(lag); c.integer(warmup)
    if lag+warmup+4 > interpreter.MAX_ATTEMPTS:
        raise ValueError('the complete native probe exceeds the attempt ceiling')
    pair = [{'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(edge)},
             'law': 'E', 'repeat': 1} for edge in c.PROBE]
    prefix = ([{'op': 'repair', 'selector': {'kind': 'uniform'},
                'law': 'mixture', 'repeat': warmup}] if warmup else [])
    middle = ([{'op': 'repair', 'selector': {'kind': 'uniform'},
                'law': 'mixture', 'repeat': lag}] if lag else [])
    return interpreter.program(biased=False, commands=prefix+pair+middle+pair,
                               capacity=1, overflow='export_oldest')


def scored_step(value, command, density, side=None):
    """A linear scored summary of the complete original outcome list.

    This is not a replacement for the retained instrument or its host stops.
    """
    if side not in (None, 'incoming', 'outgoing'):
        raise ValueError('the declared record cut is required')
    result = defaultdict(F)
    outcomes = interpreter._outcomes(value, command, density)
    for outcome in outcomes:
        record = outcome.record
        if side is None:
            score = 1
        elif record.kind == 'repair' and record.law == 'E' and record.label is not None:
            score = record.label if side == 'incoming' else record.label-2*record.current
        else:
            score = 0
        for key, weight in outcome.density.items():
            result[key] += score*weight
    return {key: weight for key, weight in result.items() if weight}, len(outcomes)


def exact_scores(max_lag=4):
    c.integer(max_lag, 0, 12)
    value = program(max_lag)
    pair = value['commands'][:2]
    unread = {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 1}
    density = {(i, 0, i, 0): F(1, 66) for i in range(66)}
    signed = density
    for command in pair:
        signed, slots = scored_step(value, command, signed, 'outgoing')
        if slots != 25:
            raise ArithmeticError('the full native fixed-seam alphabet changed')
    h = c.data()['probe']; g = c.data()['two']; transfer = c.sp.eye(66)-g/60
    vector = h; rows = []
    for lag in range(max_lag+1):
        terminal = signed
        for command in pair:
            terminal, slots = scored_step(value, command, terminal, 'incoming')
        score = sum(weight for (i, a, j, b), weight in terminal.items() if i == j and a == b)
        expected = F(2, 33)*c.fraction((h.T*vector)[0]/4)
        if score != expected:
            raise ArithmeticError('complete native scored records differ from the contact response')
        rows.append({'lag': lag, 'raw_score': score, 'normalized_response': score/F(2, 33)})
        if lag < max_lag:
            signed, slots = scored_step(value, unread, signed)
            if slots != 721:
                raise ArithmeticError('the complete native uniform alphabet changed')
            vector = transfer*vector
    return rows


def preparation_certificate(precision_bits=10):
    # The actual local one-body gap is (5-sqrt5)/60 > 1/30.
    budget = preparation.preparation_budget(12, 2, F(1, 30), precision_bits)
    value = program(4, budget['warmup_attempts'])
    return {**budget, 'admitted_program_digest': interpreter.digest(value),
            'attempts_upper': value['declared_resources']['attempts_upper'],
            'all_warmup_records_retained_but_ignored_by_score': True,
            'claim_is_signed_reference_score_not_full_state_preparation': True}


def sampling_budget(ports, precision, failure):
    """Hoeffding budget for independent complete trials, with no success filter."""
    c.integer(ports, 4)
    precision = c.rational(precision); failure = c.rational(failure)
    if precision <= 0 or not 0 < failure < 1:
        raise ValueError('positive precision and a failure probability in (0,1) are required')
    variance = F(8, ports*(ports-1))
    ratio = F(2, 1)/failure
    ceiling = (ratio.numerator+ratio.denominator-1)//ratio.denominator
    log_upper = ceiling.bit_length()
    required = F(2*log_upper, 1)/(variance*variance*precision*precision)
    count = (required.numerator+required.denominator-1)//required.denominator
    return {'ports': ports, 'q': 2, 'raw_variance': variance,
            'normalized_sampling_error': precision, 'failure_probability_upper': failure,
            'independent_complete_trials': count,
            'independence_and_repeated_preparation_are_required': True}


def certificate():
    value = program(4)
    stopped = interpreter.execute(value, {(0, 0, 0, 0): F(1)}, max_branches=1)
    if stopped.status != 'RESOURCE_STOP' or stopped.attempted_steps != 0:
        raise ArithmeticError('the literal expanded resource-stop control changed')
    return {'program_digest': interpreter.digest(value), 'attempts_upper': 8,
            'FIFO_capacity': 1, 'maximum_exported_records_on_complete_history': 7,
            'complete_probe_pair_slots': 625, 'uniform_attempt_slots': 721,
            'exact_scored_records': exact_scores(),
            'expanded_control_status': stopped.status,
            'source_caps_changed': False, 'scored_factorization_is_expanded_execution': False,
            'preparation': preparation_certificate(),
            'independent_sampling': [sampling_budget(v, F(1, 100), F(1, 1000))
                                     for v in (12, 240, 960)]}
