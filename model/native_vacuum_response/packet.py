"""Deterministic source-bound certificate for finite transfer-ground response."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_correlation_preparation import packet as preparing
from ..native_exchange_clock import packet as clock
from . import response as r

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple(f'model/native_vacuum_response/{name}' for name in (
    '__init__.py', 'response.py', 'packet.py', 'test_response.py', 'PROOF.md', 'README.md'))


def source_bindings():
    names = set(preparing.source_bindings()) | set(clock.source_bindings()) | set(SOURCE_FILES)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def local_case():
    source = r.exchange.sector.source
    family = source.local_family(source.local.Parameters(port_weights=(1,)*12))
    first, second = next((a, b) for a in family.seams for b in family.seams if not set(a)&set(b))
    h = r.exchange.sector.contrast(first, second, 12)
    return family, first, second, h


def phase_blindness():
    """Every actual H66 raw slot retains the phase-blind reference marginal."""
    source = r.exchange.sector.source.interpreter
    command = {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 1}
    program = source.program(biased=False, commands=[command], capacity=1, overflow='export_oldest')
    # Two configuration basis vectors entangled with a two-dimensional R.
    original = {(0, 0, 0, 0): F(1, 2), (1, 1, 1, 1): F(1, 2),
                (0, 0, 1, 1): F(1, 2), (1, 1, 0, 0): F(1, 2)}
    phased_real = {(0, 0, 0, 0): F(1, 2), (1, 1, 1, 1): F(1, 2)}
    phased_imag = {(0, 0, 1, 1): F(-1, 2), (1, 1, 0, 0): F(1, 2)}
    values = [source._outcomes(program, command, density)
              for density in (original, phased_real, phased_imag)]

    def marginal(density):
        out = {}
        for (x, a, y, b), value in density.items():
            if x == y:
                out[(a, b)] = out.get((a, b), F())+value
        return {key: value for key, value in out.items() if value}

    zero = 0; work_changed = 0
    for ordinary, real, imag in zip(*values, strict=True):
        if not (ordinary.record == real.record == imag.record):
            raise ArithmeticError('the complete source raw address changed')
        if marginal(ordinary.density) != marginal(real.density) or marginal(imag.density):
            raise ArithmeticError('a bare native record acquired diagonal-phase sensitivity')
        zero += int(not ordinary.density)
        work_changed += int(ordinary.density != real.density or bool(imag.density))
    if len(values[0]) != 721 or work_changed == 0:
        raise ArithmeticError('the complete raw alphabet or terminal-state distinction changed')
    # Check every basis pair in each literal partial map, not just the probe.
    family, _, _, _ = local_case()
    count = 0
    for edge in family.seams:
        for law in ('C', 'E'):
            for coin, label, pairs in source.local.raw_maps(law, edge, True):
                domain = [x for x, _ in pairs]; image = [y for _, y in pairs]
                if len(domain) != len(set(domain)) or len(image) != len(set(image)):
                    raise ArithmeticError('a fine source map ceased to be a partial injection')
                count += 1
    return {'program_digest': source.digest(program), 'complete_raw_slots': len(values[0]),
            'zero_slots_on_entangled_probe': zero, 'literal_partial_maps_checked': count,
            'terminal_work_changed_slots': work_changed,
            'reference_dimension': 2, 'complete_record_reference_marginals_equal': True,
            'terminal_work_equality_asserted': False}


def produce():
    before = source_bindings()
    family, first, second, h = local_case()
    local = r.continuation(30, 5, 2, F(1, 64), 12)
    values = r.exact_native_lags(family, h, 2, local['lag_attempts'])
    moments = r.exchange.moments(family, h, 2, 8)
    taylor = r.unitary_taylor(moments, F(1, 64))
    phase = phase_blindness()
    value = {'schema': 'oph.native_vacuum_response.v1', 'source_bindings': before,
        'local_source': {'ports': 12, 'seams': 30, 'occupation': 2,
            'probe_seams': (first, second), 'continuation': local,
            'exact_normalized_native_lag_scores': values,
            'reconstructed_complex_ground_response': r.reconstruct(local['power_coefficients'], values),
            'exchange_taylor_order': 8, 'exchange_taylor_value': taylor,
            'exchange_taylor_error_upper': (F(20, 64)**9)/362880,
            'finite_intervention_witness': r.response_witness(family, h),
            'complete_record_phase_control': phase},
        'federation_budgets': [r.source_budget(level, 80*4**level) for level in (0, 1, 2)],
        'all_domain_theorem': {
            'ground': 'normalized constant configuration vector, distinct from I/binomial(V,q)',
            'generator': '-2*M*log(P_q), P_q=I-G_q/(2*M)',
            'positive_energy_ground_matrix_response': True,
            'counted_time_convention': 'native attempt lag divided by2*M',
            'blocked_transfer': 'b=floor(M/(2*r*Delta)), Q=P^b>=I/2, M>=2*r*Delta',
            'error': 'epsilon_m+A_m*eta+s with sampling failure probability at most p',
            'stop_allowance': 'A_m*p_stop/v_q only with a proved complete stop bound',
            'native_occupation_domain': '2<=q<=V-2; bounded primitive degree r=2',
            'matrix_operator_error_requires_probe_Gram_factor': True},
        'scope': {'bare_native_records_and_backaction_attached': True,
            'complex_statistic_is_one_native_effect': False,
            'transfer_ground_prepared_on_working_carrier': False,
            'unitary_real_time_executed_by_native_repair': False,
            'displayed_repetition_counts_executed': False,
            'expanded_host_resource_limits_removed': False,
            'all_source_levels_admit_the_sufficient_warmup': False,
            'arbitrary_reference_score_preparation_preserved': True,
            'same_intervention_dictionary_has_explicit_obstruction': True,
            'all_encodings_or_coherent_controls_excluded': False,
            'Higgs_reference_identification': False,
            'physical_relativistic_vacuum_or_mass_identified': False,
            'canonical_physical_premise_closed': False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during exact replay')
    return clock.encode(value)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('transfer-ground response receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native transfer-ground response receipt written.')
    else:
        verify()
        print('Native transfer-ground response receipt verified.')
