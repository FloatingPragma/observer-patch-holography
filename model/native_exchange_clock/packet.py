"""Source-bound exact replay of the native exchange-clock comparison."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_correlation_sectors import packet as owning, records
from . import exchange as e

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple(f'model/native_exchange_clock/{name}' for name in (
    '__init__.py', 'exchange.py', 'packet.py', 'test_exchange.py', 'PROOF.md', 'README.md'))


def source_bindings():
    names = set(owning.source_bindings()) | set(SOURCE_FILES)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sorted(names)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n'


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: encode(item) for key, item in value.items()}
    if type(value) in (tuple, list):
        return [encode(item) for item in value]
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError('only exact JSON scalars and fractions enter this receipt')


def check_source_generator():
    """Compare all fair H66 columns with independently enumerated native coins."""
    source = e.sector.source
    family = source.local_family(source.local.Parameters())
    masks = source.sector_masks(12, 2)
    seams = len(family.seams)
    for mask in masks:
        subset = tuple(i for i in range(12) if mask >> i & 1)
        native = source.population_column(family, mask)
        g = e.generator(family, {subset: 1}, 2)
        for target in set(native) | {sum(1 << i for i in y) for y in g}:
            y = tuple(i for i in range(12) if target >> i & 1)
            expected = F(target == mask) - F(g.get(y, 0), 2 * seams)
            if native.get(target, F(0)) != expected:
                raise ArithmeticError('local exchange generator differs from native endpoint counts')
        counts = e.sector.counts_step(family, {subset: 1}, 2)
        predicted = {y: 2 * seams * (y == subset) - g.get(y, 0)
                     for y in set(g) | {subset}}
        if {y: value for y, value in predicted.items() if value} != counts:
            raise ArithmeticError('local exchange generator differs from owning transfer counts')
    return {'fair_H66_columns': len(masks), 'native_seam_coin_columns': len(masks)*2*seams,
            'source_law': 'fair native endpoint coins; independent uniform selector',
            'complete_raw_attachment': owning.attachment()}


def produce():
    before = source_bindings()
    attachment = check_source_generator()
    rows = []
    for level in (0, 1, 2):
        family, first, second, h = e.sector.probe_data(level)
        ports, seams, maximum_degree = e.graph_data(family)
        if seams != 630 * 4**level or maximum_degree != 6:
            raise ArithmeticError('actual fair federation graph normalization changed')
        if owning.scaling.transfer_floor(level) != F(10, 21):
            raise ArithmeticError('owning positive transfer floor changed')
        attempts = seams // 16
        row = e.response_intervals(family, h, 2, attempts, 8)
        occupations = []
        for occupied in (2, ports // 3, ports // 2):
            variance = 4 * e.sector.norm_factor(ports, 2, occupied)
            program = records.probe_program(level, occupied, attempts)
            if program['declared_resources']['attempts_upper'] != attempts + 4:
                raise ArithmeticError('complete native probe program lost attempts')
            occupations.append({'global_occupation': occupied, 'variance': variance,
                'program_digest': e.sector.source.interpreter.digest(program),
                'complete_attempts': attempts + 4,
                'score_preparation_error_for_half_trace_epsilon': f'2*epsilon/({variance})'})
        row.update({'level': level, 'ports': ports, 'seams': seams,
                    'maximum_degree': maximum_degree,
                    'probe_seams': [list(first), list(second)],
                    'occupations': occupations,
                    'target_time_comparison': e.rounded_time_budget(seams, F(1, 32), 6, 2),
                    'half_filled_full_generator_obstruction':
                        e.full_generator_obstruction(ports, seams, ports // 2)})
        rows.append(row)
    value = {'schema': 'oph.native_exchange_clock.v1', 'source_bindings': before,
        'native_attachment': attachment, 'evaluated_federations': rows,
        'all_domain_bounds': {
            'positive_native_transfer_floor': '10/21',
            'counted_time': 't_n=n/(2M)',
            'full_sector_endpoint_error': 'min(1,3/(5*n)) for n>=1; zero at n=0',
            'floor_time_error': 'min(1,(3/5+kappa/2)/n), n=floor(2*M*t)>=1',
            'degree_r_spectral_ceiling': '2*r*Delta',
            'degree_r_generator_error': 'r*r*Delta*Delta/(M-r*Delta), M>r*Delta',
            'full_half_filled_generator_norm_lower': 'M/16',
            'normalized_prepared_score_composition': 'eta/v_q + endpoint_error',
            'matrix_response_norm_allowance': 'endpoint_error * norm(probe_Gram)',
            'deterministic_vs_poisson_mean_one_record_length_TV': '1-exp(-1)>1/2'},
        'scope': {
            'exact_ideal_native_population_transfer': True,
            'default_weighted_H66_source_covered': False,
            'complete_four_E_score_backaction_included': True,
            'all_occupation_endpoint_response_bound': True,
            'all_occupation_generator_convergence': False,
            'source_attempt_clock_is_read_from_specification': True,
            'Poisson_scheduler_required': False,
            'sparse_moment_evaluation_is_expanded_history_execution': False,
            'complete_record_law_or_quantum_instrument_equivalence': False,
            'raw_quantum_repair_is_unitary': False,
            'physical_time_or_chiral_action_selected': False,
            'continuum_or_particle_mass_identification': False,
            'new_canonical_premise_discharged': False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during exact replay')
    return encode(value)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('native exchange-clock receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native exchange-clock receipt written.')
    else:
        verify()
        print('Native exchange-clock receipt verified.')
