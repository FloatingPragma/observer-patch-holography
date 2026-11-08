"""Complete native four-E observations of the connected resolvent."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb

from model.native_contact_response import contact, records as local_records
from model.native_pair_response import records as pair_records
from model.native_correlation_sectors import sectors
from . import dressing

source = pair_records.source


def program(level=0, lag=0, warmup=0, cell=0):
    p = dressing.problem(level, cell)
    contact.integer(lag); contact.integer(warmup)
    if lag+warmup+4 > source.MAX_ATTEMPTS:
        raise ValueError('four complete native probes and warmup exceed the attempt ceiling')
    offset = p.core[0]; edges = tuple(tuple(offset+i for i in e) for e in contact.PROBE)
    value = pair_records.program(level, 2, lag, edges, edges)
    if warmup:
        prefix = [{'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': warmup}]
        value = source.program(federation=value['carrier'], commands=prefix+value['commands'],
                               capacity=1, overflow='export_oldest')
    return value


def exact_scores(max_lag=1):
    """Literal level-zero kernels; local canonical multiplicities are exact."""
    contact.integer(max_lag, 0, 1)
    value = program(lag=max_lag); p = dressing.problem(0); v = len(p.neighbours); total = comb(v, 2)
    density = {}
    # The four fixed repairs only read the cell. Each exterior class can be
    # represented by one untouched exterior mask with its exact multiplicity.
    for inside_count in range(3):
        for inside in combinations(p.core, inside_count):
            outside = tuple(range(12, 14-inside_count))
            mask = sum(1 << i for i in (*inside, *outside))
            density[mask, 0, mask, 0] = F(comb(v-12, 2-inside_count), total)
    signed = density; fixed = value['commands'][:2]
    for command in fixed:
        signed, slots = local_records.scored_step(value, command, signed, 'outgoing')
        if slots != 25: raise ArithmeticError('complete E alphabet changed')
    h = contact.source.contrast(*contact.PROBE, v)
    expected = {(sum(1 << i for i in pair), 0, sum(1 << i for i in pair), 0): F(x, total)
                for pair, x in h.items()}
    if signed != expected:
        raise ArithmeticError('reflected past records failed to prepare the signed primitive statistic')
    family = sectors.source.federation_family(0, (1,)*12); vector = h; scale = 2*len(family.seams)
    rows = []
    for lag in range(max_lag+1):
        terminal = signed
        for command in fixed:
            terminal, slots = local_records.scored_step(value, command, terminal, 'incoming')
        score = source.trace(terminal)
        normalized = F(sectors.dot(h, vector), 4*scale**lag)
        variance = F(4, total)
        if score != variance*normalized:
            raise ArithmeticError('complete connected native records disagree with sparse source response')
        rows.append({'lag': lag, 'raw_score': score, 'normalized_response': normalized})
        if lag < max_lag:
            unread = {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 1}
            signed, slots = local_records.scored_step(value, unread, signed)
            if slots != 24*len(family.seams)+1: raise ArithmeticError('uniform alphabet changed')
            vector = sectors.counts_step(family, vector, 2)
    return rows


def abel(level, shift, cutoff):
    contact.integer(level, 0, 20); contact.integer(cutoff, 0, 4096)
    s = contact.rational(shift)
    if s <= 0: raise ValueError('strictly positive resolvent shift required')
    m = 630*4**level; discount = F(2*m)/(2*m+s)
    # The evaluated finite-cutoff API has a declared host bit budget. Its
    # analytic bound remains a formula beyond this evaluator's domain.
    if (cutoff > source.MAX_ATTEMPTS-4
            or (cutoff+1)*max(discount.numerator.bit_length(), discount.denominator.bit_length()) > 12000):
        raise ValueError('finite evaluated cutoff or rational bit budget exceeded')
    return {'level': level, 'shift': s, 'cutoff': cutoff, 'discount': discount,
            'discarded_tail_upper': discount**(cutoff+1),
            'raw_variance': F(8, (240*4**level)*(240*4**level-1)),
            'working_attempts_per_lag_program': cutoff+4,
            'native_program_admitted_at_this_level': level <= 2,
            'all_lags_require_separate_complete_trials_or_an_admitted_selector': True}


def certificate():
    value = program(lag=1)
    stopped = source.execute(value, {(3, 0, 3, 0): F(1)}, max_branches=1)
    if stopped.status != 'RESOURCE_STOP': raise ArithmeticError('expanded host stop control changed')
    controls = []
    for level in range(3):
        val = program(level, lag=32)
        controls.append({'level': level, 'program_digest': source.digest(val),
                         'attempts_upper': val['declared_resources']['attempts_upper'],
                         'uniform_attempt_slots': 24*(630*4**level)+1})
    finite = abel(0, F(24), 512)
    if finite['discarded_tail_upper'] >= F(1, 10000):
        raise ArithmeticError('the evaluated Abel tail no longer fits the declared budget')
    return {'fixed_E_slots': 25, 'fixed_pair_slots': 625, 'fixed_pair_zero_slots': 589,
            'literal_connected_scores': exact_scores(), 'admitted_program_controls': controls,
            'FIFO_capacity': 1, 'complete_records_exported_except_final_slot': True,
            'expanded_execution_control': stopped.status,
            'signed_kernel_factorization_is_completed_expanded_execution': False,
            'finite_Abel_example': finite,
            'independent_trials_per_lag_example': local_records.sampling_budget(240, F(1, 1000), F(1, 513000)),
            'canonical_q2_reference_or_scored_preparation_error_required': True,
            'arbitrary_reference_joint_error_assumed_zero': False,
            'successful_records_or_surviving_histories_postselected': False}
