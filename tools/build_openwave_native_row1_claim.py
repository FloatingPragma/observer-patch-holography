"""Build the source-bound finite native spatial-charge claim register."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_row1_native'
OUTPUT = ROOT / 'claims/openwave_native_row1.json'
PRIMARY_CHECKS = {
    'complete_native_source_histories', 'actual_source_current_counter_readout',
    'actual_native_pointer_all_unknown_inputs',
    'external_joined_executor_matches_independent_composition',
    'all_unconditional_apparatus_laws_complete', 'all_primary_regions',
    'all_required_mutations_rejected', 'actual_joined_timeout_stages_retained',
    'counter_complements_retained',
}
SPATIAL_CHECKS = {
    'native_kernel_balanced_form', 'complete_native_histories',
    'native_record_spatial_degree_identity',
    'nonlinear_homotopy_and_shared_vertex_gauge_invariance',
    'no_operation_events_retained',
    'half_zero_open_and_forgotten_current_mutations_rejected',
}
MUTATIONS = {'forgotten_native_current_rejected', 'larger_charge_alias_rejected',
             'unconditional_failure_omission_rejected'}
REVIEW_CHECKS = {
    'JOINED_READOUT_INDEPENDENT_REVIEW.json': {
        'all_frozen_source_histories_reconstructed_independently',
        'all_holds_and_moves_retained', 'full_conditional_probability_laws_exact',
        'reset_route_and_measurement_timeouts_independently_split',
        'controller_runs_without_terminal_load_or_background_audit_fields',
        'fine_history_coin_factors_exact', 'timeout_has_no_decoded_charge',
        'source_prefix_clock_explicit', 'source_hashes_unchanged',
    },
    'JOINED_READOUT_REVIEW.json': {
        'all456_source_histories_independently_replayed',
        'full_joined_execution_native_records_and_combined_probability',
        'all19008_unconditional_laws_reconstructed_with_separate_timeouts',
        'no_terminal_load_oracle_dependency_mutation_passes',
        'source_idle_and_global_attempt_clock_explicit',
        'zero_and_all_failure_outputs_remain_distinct',
        'maximum_attempt_bound627', 'reviewed_sources_unchanged',
    },
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inside(package, name):
    path = (package / name).resolve()
    require(path.is_relative_to(package.resolve()), 'path escapes public package')
    return path


def named_checks(checks, expected, label):
    require(isinstance(checks, dict) and set(checks) == expected,
            label + ': incomplete or unknown check names')
    require(all(value is True for value in checks.values()), label + ': failed/nonboolean check')


def passed(receipt, label, expected):
    checks = receipt.get('checks')
    named_checks(checks, expected, label)
    require(receipt.get('official_rows_changed') == [], label + ': unexpected official award')


def build(package=PACKAGE):
    package = Path(package)
    manifest = json.loads((package/'manifest.json').read_text())
    outputs = json.loads((package/'OUTPUTS.json').read_text())
    whitelist = manifest['public_file_whitelist']
    require(len(whitelist) == len(set(whitelist)), 'duplicate whitelist file')
    require(set(outputs['public_whitelist']) == set(whitelist), 'output whitelist differs')
    for filename in whitelist:
        require(inside(package, filename).is_file(), 'missing public file: ' + filename)
    for ledger in (manifest, outputs):
        for filename, expected in ledger['files_sha256'].items():
            require(filename in whitelist, 'hashed file outside whitelist')
            require(sha(inside(package, filename)) == expected, 'stale public file: ' + filename)
    primary = json.loads((package/'PRIMARY_RESULT.json').read_text())
    passed(primary, 'primary', PRIMARY_CHECKS)
    named_checks(primary['spatial_checks'], SPATIAL_CHECKS, 'spatial')
    require(primary['public_manifest_sha256'] == sha(package/'manifest.json'), 'primary manifest changed')
    require(primary['full_criterion_met'] is False, 'platform award cannot be inferred locally')
    require(primary['source_histories'] == 456 and primary['joined_executor_comparisons'] == 5472,
            'primary source domain changed')
    require(primary['source_events'] == 888 and primary['source_noops'] == 564,
            'source outcomes omitted')
    named_checks(primary['mutations'], MUTATIONS, 'mutations')
    failures = primary['actual_joined_failure_traces']
    require({r['execution']['reason'] for r in failures} ==
            {'discovery_timeout', 'routing_timeout', 'measurement_timeout'}, 'missing failure stage')
    for failure in failures:
        execution = failure['execution']
        require(execution['success'] is False and execution['pointer']['success'] is False,
                'timeout reported success')
        require(all(execution[key] is None for key in
                    ('relative_charge', 'spatial_degree', 'pointer_winding'))
                and execution['pointer']['winding'] is None, 'timeout decoded as charge')
    templates = primary['apparatus_templates']
    require(len(templates) == 36 and
            {(t['initial_pointer_port'], t['counter_class']) for t in templates} ==
            {(port, label) for port in range(12) for label in (0, 1, 5)},
            'apparatus input columns incomplete')
    decode = {0: 0, 1: 1, 5: -1}
    for template in templates:
        law = template['unconditional_law']
        require(set(law) == {'success', 'timeout'}, 'apparatus outcome names changed')
        success, timeout = Fraction(law['success']), Fraction(law['timeout'])
        require(0 <= success <= 1 and 0 <= timeout <= 1 and success + timeout == 1,
                'apparatus success/timeout law is incomplete')
        w = decode[template['counter_class']]
        require(template['decoded_spatial_winding'] == w and
                template['decoded_relative_charge'] == -w, 'apparatus charge sign changed')
    histories = primary['source_to_apparatus_factorized_joint_law']
    require(len(histories) == 456, 'factorized source histories incomplete')
    setting_laws, seen_histories = {}, set()
    for history in histories:
        setting, coins = history['setting'], tuple(history['source_coin_record'])
        require((setting, coins) not in seen_histories and all(c in (0, 1) for c in coins),
                'duplicate or invalid source coin history')
        seen_histories.add((setting, coins))
        probability = Fraction(history['source_probability'])
        require(probability == Fraction(1, 2 ** len(coins)), 'source coin law changed')
        setting_laws[setting] = setting_laws.get(setting, Fraction(0)) + probability
        readings = history['contour_readings']
        require(len(readings) == 12 and {tuple(r['region_ports']) for r in readings} ==
                {(port,) for port in range(12)}, 'source contour readings omitted')
        for reading in readings:
            port = reading['region_ports'][0]
            w = history['initial_G'][port] - history['present_n'][port]
            require(w in (-1, 0, 1) and reading['relative_charge'] == -w and
                    reading['regional_counter_class'] == w % 6 and
                    history['counter'][port] == w % 6, 'current/charge map changed')
            indices = reading['apparatus_template_indices']
            require(len(indices) == 12 and len(set(indices)) == 12 and
                    all(type(i) is int and 0 <= i < len(templates) for i in indices),
                    'missing or invalid pointer input column')
            require({templates[i]['initial_pointer_port'] for i in indices} == set(range(12))
                    and all(templates[i]['counter_class'] == w % 6 for i in indices),
                    'source to apparatus law mismatched')
    require(set(setting_laws) == set(range(132)) and
            all(p == 1 for p in setting_laws.values()), 'source settings not normalized')
    origins = json.loads((package/'provenance/ORIGINS.json').read_text())
    custody = origins['pinned_source_custody']
    for original, record in custody.items():
        require(record['bundled_path'] in whitelist, 'unbundled review source: ' + original)
        require(sha(inside(package, record['bundled_path'])) == record['sha256'],
                'changed review source: ' + original)
    reviews = []
    for filename in ('JOINED_READOUT_INDEPENDENT_REVIEW.json', 'JOINED_READOUT_REVIEW.json'):
        path = package/'provenance'/filename
        review = json.loads(path.read_text())
        passed(review, filename, REVIEW_CHECKS[filename])
        for original, expected in review['source_sha256'].items():
            require(original in custody and custody[original]['sha256'] == expected,
                    'review target mismatch: ' + original)
        reviews.append({'path': 'code/openwave_row1_native/provenance/' + filename,
                        'sha256': sha(path), 'checks_passed': len(review['checks']),
                        'reviewer_sha256': review['reviewer_sha256']})
    require(reviews[0]['reviewer_sha256'] != reviews[1]['reviewer_sha256'], 'duplicate reviewer')
    packet = json.loads((package/'source_packet.json').read_text())
    original_spatial = json.loads((package/'provenance/NATIVE_DUAL_WINDING.json').read_text())
    require(packet['spatial_result'] == original_spatial, 'source packet changed audited histories')
    sources = {'code/openwave_row1_native/' + name: sha(package/name) for name in whitelist}
    sources['tools/build_openwave_native_row1_claim.py'] = sha(Path(__file__))
    return {
        'schema': 'oph.native_spatial_charge_claim.v1',
        'status': 'PASS_SOURCE_BOUND_FINITE_NATIVE_SPATIAL_CHARGE_AND_COMPLETE_READOUT',
        'finite_bridge_complete': True,
        'official_openwave_status': 'untested',
        'official_rows_changed': [],
        'row': 1,
        'statement': 'Literal native currents define an integral spatial dual-cell class with boundary degree G-n. A current-driven finite observer and auxiliary native pointer read its regional degree with every hold and timeout retained. Integer closed-map degrees exclude a split of one unit into two closed half-unit classes.',
        'unit_ancestry': 'The integer event alphabet and nonnegative total-one endpoint law supply the indivisible source unit. Continuity also holds for fractional currents; integer degree alone does not test preservation of the whole source unit.',
        'domain': 'Twelve-port icosahedral carrier; one nonnegative load at known initial G; zero initial current tape; prescribed source schedules and fair coins; retained reference; six-state counters injective on regional degree {-1,0,1}; source idles during addressed terminal apparatus.',
        'topology': 'The integral quotient Z^30/im(primal face boundaries) maps isomorphically to the neutral integer twelve-port charge lattice. The spatial circle-map representative is a diagnostic realization of this class; it is not an independently evolving electromagnetic field.',
        'readout': 'Unknown one-load auxiliary pointer is discovered, reset and moved through native events. No coherent gates or ideal terminal projector. Up to 627 joined native attempts; each configured preparation/region/pointer state has a complete unconditional law.',
        'measurement_scope': 'The pointer reads current-driven counters. Sampled phase contours test an observer encoding of the same records, not an independently evolving field. Missing moving records are detected against source loads; the meter has no internal record-completeness oracle.',
        'source_histories': primary['source_histories'],
        'joined_executor_comparisons': primary['joined_executor_comparisons'],
        'retained_source_holds': primary['source_noops'],
        'worst_primary_meter_timeout': primary['worst_primary_meter_timeout'],
        'independent_full_join_reviews': reviews,
        'independent_full_join_review_count': len(reviews),
        'interpretation_boundary': 'Finite static mechanism and native observer execution are established on this source. Physical electromagnetic normalization, transfer to the interacting coherent source, and platform interpretation/admission require separate evidence. No mass, spin or force claim is made.',
        'opposing_interpretation': 'Independent opposing reviews confirm the finite identities and support at most a partial physical interpretation because the integer unit is supplied by the source law and the phase is an observer encoding. Numerical agreement is not an agreement to award validation.',
        'law_control_addendum': 'code/openwave_row1_protocol/law_controls.py enumerates fresh development histories and executes fractional, integer-fragmenting and omitted-record laws with fixed sampling; its separately pinned receipt does not alter frozen primary provenance.',
        'reproduce': 'python3 code/openwave_row1_native/run.py --output /tmp/openwave-native-row1.json',
        'source_sha256': sources,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
    else:
        require(json.loads(OUTPUT.read_text()) == result, 'claim register differs from validated sources')
    print(result['status'])


if __name__ == '__main__':
    main()
