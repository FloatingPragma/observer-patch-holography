"""Reproduce the complete noisy observer and its native memory attachment.

Receipts keep compact diagnostics; replay reconstructs every quantum block,
source phase branch, visible record and raw native outcome. No trajectory
dump or selected successful branch replaces those computations.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import continuation_bounds, encoded_observer, infinite_history, protocol, smooth_history
from . import open_bounds, open_environment as environment, open_native, open_protocol as mixed


STEPS = 512
SOURCE_PATHS = (
    'model/common_source/open_observer.py', 'model/common_source/open_protocol.py',
    'model/common_source/open_environment.py', 'model/common_source/open_bounds.py',
    'model/common_source/open_native.py', 'model/common_source/protocol.py',
    'model/common_source/smooth_history.py', 'model/common_source/infinite_history.py',
    'model/common_source/continuation_bounds.py', 'model/common_source/refinement.py',
    'model/common_source/scalar_coupling.py', 'model/common_source/smooth_protocol.py',
    'model/common_source/experiment.py', 'model/common_source/gauge_matter.py',
    'model/common_source/quantum_native.py', 'model/common_source/native_routing.py',
    'model/common_source/encoded_memory.py', 'model/common_source/encoded_observer.py',
    'model/anharmonic_transfer/environment.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py',
    'model/observer/source_scalar/compiler.py', 'model/observer/source_scalar/field.py',
    'model/observer/source_scalar/bounds.py', 'model/observer/source_scalar/cutoff.py',
    'model/observer/source_discrimination/circuit.py',
    'model/observer/source_discrimination/experiment.py', 'model/observer/native.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/requirements.txt',
)


def canonical_environment():
    """Four supplied phase variances, after the four Hamiltonian segments."""
    return environment.Environment(tuple(environment.PhaseNoise(cursor, 0, F(1, 64), 8)
                                         for cursor in (1, 4, 6, 9)))


def _distance(left, right):
    if set(left) != set(right):
        raise ArithmeticError('complete record menus differ')
    return math.fsum(float(np.abs(np.linalg.eigvalsh(
        (left[key]-right[key]+(left[key]-right[key]).conj().T)/2)).sum())/2 for key in left)


def _conditional_purity(history):
    return math.fsum(float(np.trace(b.density@b.density).real)/b.probability
                     for b in history.branches if b.probability > 0)


def _diagnostics(*values):
    if any(not math.isfinite(value) or not 0 <= value <= 1+2e-9 for value in values):
        raise ArithmeticError('complete history diagnostics must be finite and in their probability/distance range')


def correlated_draw_control(program, env):
    """A deliberately different source channel: reuse one draw across all gaps.

    Every source phase unitary is still executed, but the hidden alternatives
    are kept correlated. This is the negative control for treating a cached
    channel as a persistent stochastic sample instead of an independent draw.
    """
    env = environment.validated(program, env)
    if not env.events or any((e.site, e.variance, e.order) !=
                            (env.events[0].site, env.events[0].variance, env.events[0].order)
                            for e in env.events):
        raise ValueError('correlated control needs identical repeated phase menus')
    event = env.events[0]
    environment.event_metadata(program, event)
    channel = environment.source.channel(program.dimension, 1., float(event.variance), event.order)
    occupation = (np.arange(program.work_dimension)//(4*program.dimension**event.site)) % program.dimension
    prepared = protocol.execute(program, stop_after=0).branches[0].state.ravel()
    boundaries = {e.after for e in env.events}
    totals = {}
    for weight, unitary in zip(channel.weights, channel.unitaries, strict=True):
        diagonal = np.diag(unitary)[occupation]
        rows = {(): prepared*diagonal if 0 in boundaries else prepared}
        for cursor, operation in enumerate(program.operations, 1):
            outputs = {}
            for record, vector in rows.items():
                selected = record[operation.record] if operation.kind == 'feedback' else None
                maps, _ = protocol._operation_maps(program, operation, selected, 'source')
                for bit, matrix in enumerate(maps):
                    key = (*record, bit) if operation.kind in protocol.READS else record
                    value = matrix@vector
                    outputs[key] = value*diagonal if cursor in boundaries else value
            rows = outputs
        for record, vector in rows.items():
            totals.setdefault(record, np.zeros((program.work_dimension,)*2, complex))
            totals[record] += float(weight)*np.outer(vector, vector.conj())
    return totals


def _reference_control(env):
    program = protocol.continuation_program(2, STEPS)
    amplitudes = np.zeros((16, 2), complex)
    amplitudes[0, 0], amplitudes[15, 1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    initial = np.outer(amplitudes.ravel(), amplitudes.ravel().conj())
    history = mixed.execute(program, env, initial=initial)
    final = sum(b.density for b in history.branches)
    reference = np.trace(final.reshape(16, 2, 16, 2), axis1=0, axis2=2)
    residual = float(np.linalg.norm(reference-np.eye(2)/2, 2))
    _diagnostics(residual)
    if residual > 2e-10 or dict(history.branches[0].resources)['work_preparations'] != 0:
        raise ArithmeticError('mixed observer lost its reference or re-prepared its supplied input')
    return {'reference_dimension': 2, 'visible_outcomes': len(history.branches),
            'reference_marginal_residual_diagnostic': residual,
            'supplied_entangled_input_replaces_vacuum': True,
            'vacuum_occupation_certificate_applied': False}


def _prepared_noiseless_bound(program, infinite):
    finite = smooth_history.enclose(program, dimension=program.dimension)
    _, upper = smooth_history.distance_interval(infinite, finite)
    temporal = F(continuation_bounds.program_bound(program)['ideal_temporal_instrument_error_upper'])
    rounded = infinite_history._interval(F(0), min(F(1), upper+temporal))[1]
    return F(rounded)


def run():
    env = canonical_environment()
    rows, source_channels = [], {}
    for current in (-1, 0, 1):
        infinite = smooth_history.enclose(protocol.continuation_program(2, STEPS, native_current=current))
        for dimension in (2, 4):
            program = protocol.continuation_program(dimension, STEPS, native_current=current)
            source = mixed.execute(program, env)
            gaussian = mixed.execute(program, env, backend='continuous', noise_backend='gaussian')
            finite_gaussian = mixed.execute(program, env, backend='product', noise_backend='gaussian')
            noiseless = mixed.execute(program)
            prefix = mixed.execute(program, env, stop_after=5)
            restored = mixed.checkpoint_from_json(program, env, mixed.checkpoint_to_json(prefix))
            resumed = mixed.execute(program, env, checkpoint=restored)
            continuation_error = mixed.complete_distance(source, resumed)
            _diagnostics(continuation_error)
            if continuation_error > 2e-12 or source.expanded_tree_resources != resumed.expanded_tree_resources:
                raise ArithmeticError('mixed checkpoint continuation changed the history or costs')
            prefix_error = max(abs(a-b) for a, b in zip(
                mixed.prefix_probabilities(prefix, 2).values(),
                mixed.prefix_probabilities(source, 2).values(), strict=True))
            _diagnostics(prefix_error)
            if prefix_error > 2e-10:
                raise ArithmeticError('continued noise or feedback rewrote a retained record')
            densities = {b.record: b.density for b in source.branches}
            reused = correlated_draw_control(program, env)
            source_gaussian = mixed.complete_distance(source, gaussian)
            finite_difference = mixed.complete_distance(source, finite_gaussian)
            noise_effect = mixed.complete_distance(source, noiseless)
            correlation_effect = _distance(densities, reused)
            purity = _conditional_purity(source)
            _diagnostics(source_gaussian, finite_difference, noise_effect, correlation_effect,
                         purity, source.probability_mass)
            prepared = _prepared_noiseless_bound(program, infinite)
            budget = open_bounds.budget(program, env, prepared_noiseless_error=prepared)
            # Each source menu is reproduced at every use but retained once.
            # Full intermediate moment recurrences are likewise recomputed.
            budget['matched_phase_noise_occupation_certificate'].pop('stages')
            for event in budget['environment_budget']['events']:
                source_metadata = event.pop('source')
                key = f"d{dimension}:site{event['site']}:v{event['variance_exact']}:M{event['quadrature_order']}"
                source_channels[key] = source_metadata
                event['source_channel_id'] = key
            ideal = F(budget['ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper'])
            # Numerical comparisons remain distinct from the ideal theorem.
            if (source_gaussian > float(ideal)+2e-10 or finite_difference > 2e-9
                    or noise_effect <= 1e-5 or correlation_effect <= 1e-5
                    or purity >= .9999):
                raise ArithmeticError('complete noisy-history comparison or intervention control failed')
            counts = {key: max(dict(b.resources)[key] for b in source.branches) for key in mixed.RESOURCE_KEYS}
            for noisy, quiet in zip(source.branches, noiseless.branches, strict=True):
                for key, count in budget['environment_budget']['resources_per_realization'].items():
                    if dict(noisy.resources)[key]-dict(quiet.resources)[key] != count:
                        raise ArithmeticError('source environment was not charged on its realized path')
            rows.append({
                'dimension': dimension, 'native_current': current,
                'complete_output_probabilities': [{'record': list(b.record), 'probability': b.probability}
                                                   for b in source.branches],
                'source_mass_diagnostic': source.probability_mass,
                'mean_conditional_work_purity_diagnostic': purity,
                'source_to_continuous_gaussian_complete_distance_diagnostic': source_gaussian,
                'source_to_product_gaussian_complete_distance_diagnostic': finite_difference,
                'complete_noise_effect_diagnostic': noise_effect,
                'reused_draw_complete_distance_diagnostic': correlation_effect,
                'resume_complete_distance_diagnostic': continuation_error,
                'earlier_record_preservation_error_diagnostic': prefix_error,
                'maximum_realized_path_resources': counts,
                'expanded_visible_tree_resources': dict(source.expanded_tree_resources),
                'prepared_noiseless_certificate_recomputed': True,
                'budget': budget,
            })
    native_rows = []
    for law in ('C', 'E'):
        for dimension in (2, 4):
            active = open_native.execute(env, law, dimension, STEPS, native_mode='protected_qubit')
            passive = open_native.execute(env, law, dimension, STEPS, native_mode='passive_qubit')
            program = protocol.continuation_program(dimension, STEPS)
            required = tuple(key for key in encoded_observer.required_service_keys(program)
                             if key != 'native_terminal_decoder_rotation_and_readout')
            required += tuple(f'environment_{index}:site_{event.site}@{event.after}'
                              for index, event in enumerate(env.events))
            finite_bound = next(row['budget']['ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper']
                                for row in rows if row['dimension'] == dimension)
            native_rows.append({'law': law, 'dimension': dimension,
                                'active_complete_instrument_audit': open_native.audit(active),
                                'passive_complete_instrument_audit': open_native.audit(passive),
                                'memory_comparison': open_native.memory_witness(active, passive),
                                'protected_complete_physical_service_coverage': continuation_bounds.service_error_budget(
                                    F(finite_bound), required_stage_keys=required),
                                'active_resources': open_native.resources(active),
                                'passive_resources': open_native.resources(passive)})
    refinement = []
    template = protocol.continuation_program()
    for dimension in (64, 256, 1024):
        exact_phase = continuation_bounds.complete_error_budget(
            dimension, template.operations, dimension**2)
        order = open_bounds.gaussian_refinement(dimension, tuple(e.variance for e in env.events))
        combined = min(F(1), F(exact_phase['complete_ideal_product_to_infinite_instrument_error_upper'])
                       +F(order['ideal_exact_parameter_quadrature_error_upper']))
        # Keep enormous analytic denominators out of the retained receipt;
        # both fields are rounded outward on the established 12-place grid.
        refinement.append({'dimension': dimension, 'quadrature_orders': order['orders'],
                           'ideal_gaussian_rule_error_upper': infinite_history._interval(
                               F(0), F(order['ideal_exact_parameter_quadrature_error_upper']))[1],
                           'complete_noisy_ideal_refinement_error_upper': infinite_history._interval(F(0), combined)[1],
                           'larger_source_execution_claimed': False})
    root = Path(__file__).resolve().parents[2]
    return {'schema': 'oph.common_source.open_observer.v1',
            'source_sha256': {path: hashlib.sha256((root/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'source_environment_channels': source_channels,
            'environment': env.specification(), 'work_histories': rows,
            'native_attachments': native_rows, 'retained_reference_control': _reference_control(env),
            'analytic_noisy_refinement': refinement,
            'scope': {'canonical_program_grammar_reused': True,
                      'all_visible_records_and_mixed_work_outputs_retained': True,
                      'fresh_environment_draws_are_not_accessible_to_feedback': True,
                      'initial_native_work_independence_required_for_native_wrapper': True,
                      'model_time_and_phase_variance_are_supplied': True,
                      'noise_law_or_rate_selected_from_OPH': False,
                      'noisy_native_memory_or_archive_fault_tolerance_claimed': False,
                      'floating_or_physical_service_error_certified': False,
                      'simultaneous_Hamiltonian_Lindblad_flow_claimed': False},
            'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument('--write', type=Path)
    files.add_argument('--check', type=Path)
    args = parser.parse_args()
    expected = strict_json.load(args.check) if args.check else None
    result = run()
    if args.check:
        compare_receipt(expected, result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    print(canonical({'all_pass': result['all_pass'],
                     'complete_work_histories': len(result['work_histories']),
                     'complete_native_memory_comparisons': len(result['native_attachments'])}))


if __name__ == '__main__':
    main()
