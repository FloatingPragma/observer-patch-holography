"""Own the bounded regional-clock receipt; never rerun training or preparation."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from . import actuator as a, independent as independent
from ..local_repair_family import inference, modular_clock
from ..native_program import members

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')


def retained_weights():
    receipt = members.checked_record_receipt()
    counts = {tuple(x['seam']): (x['trials'], x['coin1_count']) for x in receipt['acquisition']['blocks']}
    decoded = inference.infer_counts(counts)
    if decoded['status'] != 'accept':
        raise ArithmeticError('retained record reference is unresolved')
    return decoded['relative_weights']


def outward(value):
    scale = 10**24
    return str(F(-(-value.numerator*scale//value.denominator), scale))


def report():
    weights = retained_weights()
    # Authenticate the exact Q(sqrt(5)) controls; numerical exponentials below
    # are diagnostics of these controls, not replacement matrix primitives.
    compiler = modular_clock.common.full.compiler
    compiler.source_summary()
    results = []
    for ports in a.REGIONS:
        tape = a.word(weights, ports)
        admitted = a.admitted(weights, ports)
        analytic = independent.angle_certificate(tape, weights, ports)
        analytic = {k: outward(F(v)) if k.endswith('_upper') else v for k, v in analytic.items()}
        array = a.executed(weights, ports)
        # Fresh independent full-address target, with rational log intervals.
        diagonal = np.ones(512, complex)
        for code, ratio in zip(independent.CODES, independent.target(weights, ports), strict=True):
            lo, hi = independent.log_interval(ratio)
            diagonal[code] = np.exp(-1j*float((lo+hi)/2))
        residual = float(np.linalg.norm(array-np.diag(diagonal), 'fro'))
        if residual > float(a.OPERATOR_CAP):
            raise ArithmeticError('independent numerical full-word reconstruction failed')
        complement_rows = admitted['declared_resources']['export_complement_outcomes_per_run']
        results.append({'ports': list(ports), 'algebra': independent.algebra_certificate(weights, ports),
                        'analytic_certificate': analytic,
                        'numeric_full512_Frobenius_diagnostic': residual,
                        'numeric_error_is_physical_precision': False,
                        'resources': admitted['declared_resources'],
                        'word_sha256': hashlib.sha256(json.dumps(admitted, sort_keys=True).encode()).hexdigest(),
                        'retained_first_processor_failure_labels': complement_rows,
                        'retained_final_address_complement_labels': 446,
                        'processor_failure_matrices_materialized': False,
                        'exact_source_address_abort_probability_upper': outward(F(analytic['analytic_operator_error_upper'])**2)})
    r, s = map(set, a.REGIONS)
    union, overlap = tuple(sorted(r | s)), tuple(sorted(r & s))
    qr, qs = (independent.target(weights, p) for p in a.REGIONS)
    qu, qi = (independent.target(weights, p) for p in (union, overlap))
    if any(x*y != u*i for x, y, u, i in zip(qr, qs, qu, qi, strict=True)):
        raise ArithmeticError('overlap product-potential identity failed')
    # Nonnested expectation composition has a different answer on fixed N.
    f = tuple(F(x[7]) for x in independent.STATES)
    true_weights = modular_clock.parameters().port_weights
    er = lambda x: independent.expectation_diagonal(true_weights, a.REGIONS[0], x)
    es = lambda x: independent.expectation_diagonal(true_weights, a.REGIONS[1], x)
    ei = lambda x: independent.expectation_diagonal(true_weights, overlap, x)
    defect = max(abs(x-y) for x, y in zip(er(es(f)), ei(f), strict=True))
    if defect == 0 or ei(er(f)) != ei(f) or er(ei(f)) != ei(f):
        raise ArithmeticError('nested/nonnested expectation control failed')
    paths = [*Path(__file__).parent.glob('*.py'), Path(__file__).with_name('README.md'),
             Path(__file__).with_name('contract.json')]
    paths += [ROOT/p for p in ('model/completion/program_class_v2.json',
        'model/local_repair_family/record_clock_results.json', 'model/native_program/members.py',
        'model/observer/source_scalar/compiler.py', 'code/source_selection_model/response.json',
        'code/m1_source_realization/receipt.json')]
    # The fixed contract is a decision boundary, not a mutable report label.
    contract = json.loads(Path(__file__).with_name('contract.json').read_text())
    if (contract['regions'] != [list(p) for p in a.REGIONS]
            or contract['baseline_upstream'] != '1f14865fe560f5788a990a1b25884d6f91d12858'
            or len(contract['validation_upstream']) != 40
            or any(c not in '0123456789abcdef' for c in contract['validation_upstream'])
            or F(contract['fixed_operator_error_cap']) != a.OPERATOR_CAP
            or contract['common_lift_gauge_port'] != a.BASE_PORT
            or contract['CP02_accepted'] is not False):
        raise ArithmeticError('regional implementation differs from its frozen contract')
    prior = members.checked_record_receipt()['resources']
    total_pulses = sum(x['resources']['native_pulses'] for x in results)
    return {'schema': 'oph.regional_clock.v1', 'task': 'OPH–ARITHMON–REGIONAL-CLOCK-IMPLEMENTATION-0',
            'baseline_upstream': '1f14865fe560f5788a990a1b25884d6f91d12858',
            'validation_upstream': contract['validation_upstream'],
            'status': 'SOURCE_BOUND_REGIONAL_CLOCK_IMPLEMENTED',
            'scope': 'approximate complete-channel implementation using the accepted supplied services',
            'reference': 'learned conditional product on H66; stationary true product is a separate comparator',
            'dimensionless_parameter': '1 per exhibited word', 'regions': results,
            'overlap': list(overlap), 'union': list(union),
            'exact_join': 'V_R(s)V_S(t)=V_S(t)V_R(s); V_R(s)V_S(s)=V_union(s)V_overlap(s)',
            'join_scope': 'product-potential lifts with one common port-7 gauge; overlap action is trivial',
            'nonnested_expectation_tower_counterexample_on_n7': str(defect),
            'nonnested_expectation_reference': 'true stationary product; not the input to phase compilation',
            'two_orders_halfdiamond_upper': outward(2*sum(F(x['analytic_certificate']['analytic_operator_error_upper']) for x in results)),
            'fixed_training_design_failure_upper': str(inference.fixed_budget()['statistical_failure_upper']),
            'learned_vs_true_pair_lift_error_per_unit_total_parameter': '3/32',
            'finite_prepared_reference_claim': 'reduced modular comparison only; no exact global regional invariance',
            'extra_native_attempts': 0, 'extra_pointer_services': 0,
            'appended_requests_budget': {
                'source_pulses': total_pulses, 'load_export_calls': 2*total_pulses,
                'processor_blank_initializations': total_pulses,
                'fresh_invocation_address_blank_bits_upper': 9,
                'continued_invocation_source_resets': 0,
                'source_address_encoding_calls_per_request': 1,
                'source_address_final_exports_per_request': 1,
                'combined_literal_program_bytes': sum(x['resources']['program_description_bytes'] for x in results),
                'request_order_bit': 1, 'request_status_bits': 4,
                'additional_failure_metadata_bits_upper': 1+sum(
                    x['resources']['source_position_counter_bits']+3+9+2 for x in results),
                'processor_complement_environment_dimension_upper': 6,
                'address_complement_environment_dimension': 446,
                'prior_protocol_native_attempts_retained': prior['total_native_attempts'],
                'prior_protocol_record_bits_retained_upper': prior['all_committed_native_record_bits_upper'],
                'prior_plus_two_requests_pulses_upper': prior['complete_supplied_realization_pulses_upper']+total_pulses,
                'prior_plus_two_requests_load_export_calls_upper': prior['complete_supplied_realization_load_export_calls_upper']+2*total_pulses},
            'training_or_preparation_reexecuted': False,
            'CP02_accepted': False, 'physical_time_energy_or_field_gravity_join': False,
            'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    if args.check:
        retained = json.loads(OUTPUT.read_text())
        # Host BLAS diagnostics may differ; all rational claims and bindings
        # must match exactly, and each replay has already enforced its gate.
        for value in (retained, result):
            for row in value['regions']:
                residual = row.pop('numeric_full512_Frobenius_diagnostic')
                if not np.isfinite(residual) or not 0 <= residual <= float(a.OPERATOR_CAP):
                    raise ArithmeticError('invalid retained numerical diagnostic')
        if retained != result:
            raise ArithmeticError('regional receipt differs from independent replay')
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'regions': [x['ports'] for x in result['regions']],
                      'pulses': [x['resources']['native_pulses'] for x in result['regions']]}, sort_keys=True))


if __name__ == '__main__':
    main()
