"""Join independently replayed preparation and exchange-response contracts."""
from fractions import Fraction as F
import json

from ..native_correlation_preparation import preparation
from ..native_exchange_clock import exchange


def _same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(
        right, sort_keys=True, allow_nan=False)


def compose(prepared, clock):
    """Compose already authenticated owner packets; never infer a physical clock."""
    rows = []
    prep_rows = prepared['federation_preparation_budgets']
    clock_rows = clock['evaluated_federations']
    if not _same([row['level'] for row in prep_rows],
                 [row['level'] for row in clock_rows]):
        raise ValueError('preparation and response refinement families differ')
    for left, right in zip(prep_rows, clock_rows):
        graph = left['graph_certificate']
        if (not _same(graph['vertices'], right['ports'])
                or not _same(graph['seams'], right['seams'])
                or not _same(left['probe_seams'], right['probe_seams'])):
            raise ValueError('preparation and response source or probe differs')
        if F(right['native_time']) != exchange.native_time(right['seams'], right['attempts']):
            raise ValueError('response clock is not the counted native exchange time')
        budgets, occupations = left['budgets'], right['occupations']
        if not _same([item['occupation'] for item in budgets],
                     [item['global_occupation'] for item in occupations]):
            raise ValueError('preparation and response occupation differs')
        for budget, observation in zip(budgets, occupations):
            if (not _same(budget['ports'], right['ports'])
                    or not _same(budget['stationary_variance'], observation['variance'])
                    or not _same(budget['normalization'], observation['variance'])):
                raise ValueError('preparation and response normalization differs')
            variance = F(observation['variance'])
            prep_error = F(budget['unnormalized_reference_score_trace_norm_upper'])
            if prep_error/variance != F(budget['normalized_score_error_upper']):
                raise ValueError('preparation normalization does not replay')
            warmup, lag = budget['warmup_attempts'], right['attempts']
            preparation.s.integer(warmup); preparation.s.integer(lag)
            total = warmup + lag + 4
            admitted = total <= preparation.s.source.interpreter.MAX_ATTEMPTS
            digest = None
            if admitted:
                program = preparation.prepared_program(left['level'],
                    budget['occupation'], warmup, lag)
                if program['declared_resources']['attempts_upper'] != total:
                    raise ArithmeticError('composed native program cost differs')
                digest = preparation.s.source.interpreter.digest(program)
            error = exchange.preparation_composition(variance, prep_error,
                F(right['response_budget']['best_response_error_upper']))
            rows.append({'level': left['level'], 'ports': right['ports'],
                'occupation': budget['occupation'], 'variance': str(variance),
                'warmup_attempts': warmup, 'lag_attempts': lag,
                'counted_exchange_time': right['native_time'],
                'normalized_preparation_error_upper': str(prep_error/variance),
                'normalized_clock_error_upper': right['response_budget']['best_response_error_upper'],
                'combined_normalized_reference_score_error_upper': str(error),
                'local_exchange_response_interval': right['local_exchange_response_interval'],
                'complete_attempts': total, 'single_program_admitted': admitted,
                'complete_program_digest': digest, 'expanded_history_executed': False})
    return {'schema': 'oph.native_field_response.join.v1', 'responses': rows,
            'observation': 'signed four-E score with unconditioned warmup and body records summed',
            'physical_particle_identification': False,
            'full_state_or_complete_instrument_preparation': False}
