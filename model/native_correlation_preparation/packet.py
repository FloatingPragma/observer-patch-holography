"""Replay source-bound preparation of complete native correlation scores."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_correlation_sectors import packet as parent
from . import preparation as p

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple(sorted(set((*parent.SOURCE_FILES, *p.mixing.SOURCE_PATHS,
    *(f'model/native_correlation_preparation/{name}' for name in (
        '__init__.py', 'preparation.py', 'packet.py', 'test_preparation.py',
        'test_packet.py', 'PROOF.md', 'README.md'))))))


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = set(SOURCE_FILES) | set(parent.source_bindings())
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths) if not name.startswith('docs/')}


def sparse_effect_controls():
    rows = []
    for level in (0, 1, 2):
        family, first, second, h = p.s.probe_data(level)
        alternate = next(e for e in family.seams if e != second and not set(e) & set(first))
        future = (first, alternate)
        hj = p.s.contrast(*future, len(family.weights))
        moved = hj; denominator = 1
        for lag in range(4):
            effect = p.cross_effect(family, (first, second), future, lag)
            means = {}
            for q in (2, len(family.weights)//3, len(family.weights)//2):
                actual = p.uniform_mean(effect, len(family.weights), q)
                expected = p.s.norm_factor(len(family.weights), 2, q)*p.s.dot(h, moved)/denominator
                if actual != expected:
                    raise ArithmeticError('cross-effect mean differs from the exact primitive response')
                means[str(q)] = str(actual)
            payload = [[list(a), str(v)] for a, v in sorted(effect.items())]
            rows.append({'level': level, 'lag_attempts': lag,
                         'past_seams': [list(first), list(second)],
                         'future_seams': [list(e) for e in future],
                         'nonzero_coefficients': len(effect),
                         'degree_upper': max(map(len, effect), default=0),
                         'effect_sha256': hashlib.sha256(canonical(payload).encode()).hexdigest(),
                         'stationary_cross_scores': means})
            moved = p.s.counts_step(family, moved, 2)
            denominator *= 2*len(family.seams)
    return rows


def produce():
    before = source_bindings()
    result = {
        'schema': 'oph.native_correlation_preparation.v1',
        'source_bindings': before,
        'complete_native_H66_effect_controls': [p.source_effect_control(0, 2, law)
                                                for law in ('C', 'E', 'mixture')],
        'all_input_H66_prepared_response': p.prepared_h66_control(),
        'sparse_federation_cross_effect_controls': sparse_effect_controls(),
        'federation_preparation_budgets': [p.source_budgets(level) for level in (0, 1, 2)],
        'all_domain_theorem': {
            'scored_effect': 'H_past1 H_past2[phi_i * P^t(phi_j)]',
            'degree_upper': 4,
            'centered_polynomial_rank': 'binomial(V,min(4,q,V-q))-1',
            'variance': '4*(q)_2*(V-q)_2/(V)_4',
            'reference_score_trace_norm_error': 'sqrt(rank*variance)*(1-gamma)^warmup',
            'normalized_score_error': 'sqrt(rank/variance)*(1-gamma)^warmup',
            'generic_polynomial_normalized_error': 'sqrt(rank*centered_L2_squared_upper)/normalization*(1-gamma)^warmup',
            'generic_effect_requires_independent_instrument_attachment': True,
            'normalized_stop_error_addition': 'stop_probability/variance',
            'clock': 'counted independent uniform seam attempts with fair placement coin',
            'reference': 'uniform density in the declared binary q sector',
            'normalized_domain': '2 <= q <= V-2',
            'generic_polynomial_domain': 'real squarefree degree <= 4 on the declared q slice'},
        'scope': {
            'arbitrary_initial_quantum_state_and_untouched_reference': True,
            'complete_native_zero_addresses_retained': True,
            'preparation_records_ignored_without_conditioning_in_scored_observation': True,
            'cross_probes_with_distinct_past_and_future_seams': True,
            'global_density_trace_distance_bound_claimed': False,
            'joint_preparation_archive_product_law_claimed': False,
            'normalized_bound_uniform_in_capacity_or_fixed_small_occupation': False,
            'sufficient_budget_is_minimum_preparation_time': False,
            'all_budgets_fit_native_single_program_limit': False,
            'expanded_tree_success_or_postselected_timeout_claimed': False,
            'biased_H66_or_adaptive_selector_attachment': False,
            'physical_vacuum_or_particle_identification': False,
            'new_canonical_premise_discharged': False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during preparation replay')
    return result


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('native correlation-preparation receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native correlation-preparation receipt written.')
    else:
        verify()
        print('Native correlation-preparation receipt verified.')
