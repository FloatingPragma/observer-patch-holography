"""Review the proposed executable class without silently adopting it.

The class contract binds existing exposure, the source interpreter, finite
observer membership and an independent review. Admission is distinct from
physical success. Exact comparison helpers implement the declared interval
and ancestry rules; they provide no experimental likelihood or theory error.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/program_class_v2.json')
EXPOSURE_SNAPSHOTS = {
    'tracking/observation_ledger.json':
        'model/completion/historical/exposure/observation_ledger.json',
    'claims/physical_identification_registry.json':
        'model/completion/historical/exposure/physical_identification_registry.json',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def reviewed_contract_digest(packet):
    """Adoption signatures do not sign themselves; every class clause is signed.

    The interpreter binding is included. Replacing a member receipt, target,
    calibration, resource ceiling or explanatory scope invalidates the review.
    """
    return digest({key: value for key, value in packet.items() if key != 'adoption'})


def bound_file(root, name, expected):
    path = (root/name).resolve()
    require(path.is_relative_to(root) and path.is_file(), f'missing or external bound file: {name}')
    require(hashlib.sha256(path.read_bytes()).hexdigest() == expected, f'stale reviewed source: {name}')


def bound_exposure(root, name, expected):
    """Check original exposure bytes, independently of mutable live commentary.

    Only the two maintained registers use snapshots. All other exposure pins,
    interpreter inputs and review sources keep their exact live-file checks.
    The snapshot hash is the original class pin, never a regenerated value.
    """
    bound_file(root, EXPOSURE_SNAPSHOTS.get(name, name), expected)


def exact_interval(value):
    require(type(value) in (tuple, list) and len(value) == 2, 'two exact interval endpoints required')
    require(all(type(x) in (str, int, F) for x in value), 'floating or Boolean interval endpoint forbidden')
    lower, upper = map(F, value)
    require(lower <= upper, 'reversed interval')
    return lower, upper


def compare_intervals(prediction, observation, *, complete_error, same_quantity, independent_output):
    """A compatibility diagnostic only; no joint statistical interpretation.

    The three Boolean premises must have their own reviewed evidence. This
    arithmetic helper neither establishes them nor promotes calibrated output.
    """
    require(all(type(x) is bool for x in (complete_error, same_quantity, independent_output)),
            'explicit Boolean evidence status required')
    p, q = exact_interval(prediction), exact_interval(observation)
    if not (complete_error and same_quantity and independent_output):
        return 'INCONCLUSIVE'
    if p[1] < q[0] or q[1] < p[0]:
        return 'CANDIDATE_DISCREPANCY'
    if q[0] <= p[0] <= p[1] <= q[1]:
        return 'RETROSPECTIVE_INTERVAL_COMPATIBILITY'
    return 'INCONCLUSIVE'


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.executable_program_class.v2', 'wrong class schema')
    require(packet.get('proposed_program_id') == 'OPH-EXECUTABLE-COMPLETION-2026-10-29-v2',
            'a changed class needs a separately reviewed version')
    adoption = packet.get('adoption', {})
    require(adoption.get('status') in ('PROPOSED_FOR_INDEPENDENT_REVIEW', 'ADOPTED_AFTER_REVIEW'),
            'invalid adoption status')
    hypothesis = strict_json.load(root/'claims/executable_model_hypothesis.json')
    require(packet['hypothesis']['mandatory_obligations'] == [r['id'] for r in hypothesis['obligations']],
            'mandatory physical target removed or changed')
    require(packet['hypothesis']['source_selected_physical_outputs_required'] is True and
            packet['hypothesis']['imported_field_action_counts_as_physical_success'] is False,
            'supplied benchmarks cannot become source-derived physical success')
    require({r['id'] for r in packet['carriers']} ==
            {'local_integer', 'binary_native_federation', 'independent_native_product'},
            'the source class silently drops an existing native carrier family')
    federation = next(r for r in packet['carriers'] if r['id'] == 'binary_native_federation')
    require(federation['local_total_two_tensor_invariance_claimed'] is False,
            'joining seams do not preserve local total-two factors')
    require(packet['refinements']['dynamical_intertwining_assumed'] is False,
            'static refinement cannot be promoted to dynamical intertwining')
    parameters = packet['laws']['local_alternative_family']
    parent = strict_json.load(root/parameters['contract'])['parameter_boundary']
    require(parameters['port_weights'] == {'count': 12, 'min': 1, 'max': 32,
                                         'max_ratio': '4', 'primitive_gcd': 1},
            'local alternative family changed')
    require(parameters['seam_weights'] == {'count': 30, 'min': 1, 'max': 16} and
            parameters['probability_max_denominator'] == parent['probability_max_denominator'] and
            F(parameters['repair_probability_min']) == F(parent['repair_probability_min']),
            'existing alternative-law bounds changed')
    require(packet['services']['source_word_kinds'] == ['X', 'Y', 'Z', 'I', 'H', 'J'],
            'an unreviewed source primitive was introduced')
    require(packet['services']['coherent_services_derived_from_passive_repair'] is False,
            'coherence reachability boundary was erased')
    resources = packet['resource_envelope']
    require(resources['ast_nodes_max'] == 4096 and resources['source_word_literals_max'] == 65536 and
            resources['repeat_count_max'] == 2**32 and resources['fifo_capacity_min'] == 1 and
            resources['fifo_capacity_max'] == 8 and
            int(resources['source_pulses_per_finite_request_max']) == 2**48 and
            int(resources['native_attempts_per_finite_request_max']) == 2**32 and
            int(resources['archive_bits_per_finite_request_max']) == 2**40,
            'resource class changed without a successor contract')
    require(resources['host_limits_are_physical_clock'] is False and
            resources['contract_does_not_authorize_large_runs'] is True,
            'resource envelope changed scientific or execution scope')
    q3 = packet['observer_acceptance']['q3']
    require(q3['historical_status'] == 'UNEXECUTED_NOT_OBSTRUCTED' and
            q3['proposed_active_v2_disposition'] == 'NOT_ADOPTED' and
            q3['evidence_preserved'] is True, 'historical q3 evidence was rewritten')
    require(packet['observer_acceptance']['full_spatial_field_intertwiner_required_for_CP01'] is False,
            'CP-06 physical field work must not be invented as a CP-01 prerequisite')
    prospective = packet['targets']['prospective_arm']
    require(prospective['status'] == 'DISABLED_IN_V2' and prospective['eligible_datasets'] == [] and
            prospective['prospective_claim_allowed'] is False,
            'known target exposure cannot be erased by this program version')
    require(packet['targets']['physical_comparison_rule']['exposure'] == 'ALL_EXISTING_FIXTURES_RETROSPECTIVE',
            'physical comparison exposure changed')
    bindings = packet['exposure']['source_hashes']
    require(isinstance(bindings, dict) and bindings, 'missing exposure/source custody')
    for name, expected in bindings.items():
        bound_exposure(root, name, expected)
    for name in packet['exposure']['registers']:
        require(name in bindings, f'unbound prior exposure register: {name}')
    observations = strict_json.load(root/'tracking/observation_ledger.json')['rows']
    require(len(observations) == 71, 'observation target inventory changed; review successor contract')
    frozen = strict_json.load(root/EXPOSURE_SNAPSHOTS['tracking/observation_ledger.json'])['rows']
    require({r['id']: r['target'] for r in observations} ==
            {r['id']: r['target'] for r in frozen},
            'live observation identities or physical targets differ from frozen exposure')
    require(packet['targets']['structural_menu']['selection'].startswith('all71'),
            'physical target menu silently narrowed')
    require(packet['predecessor']['old_finite_campaigns_preserved'] is True and
            packet['predecessor']['prior_failures_preserved'] is True, 'prior exclusions were reset')

    blockers = []
    implementations = packet.get('review_implementation_bindings', {})
    if not implementations:
        blockers.append('Bind the contract checker, independent tests and explanatory scope.')
    else:
        required = {'model/completion/program_class_v2.py',
                    'model/completion/test_program_class_v2.py',
                    'docs/EXECUTABLE_PROGRAM_CLASS_V2_PROPOSAL.md'}
        require(required <= set(implementations), 'incomplete contract-review implementation custody')
        for name, expected in implementations.items():
            bound_file(root, name, expected)
    interpreter = packet.get('interpreter_binding')
    if not interpreter:
        blockers.append('Bind the completed typed interpreter, member receipt and independent checks.')
    else:
        require(interpreter.get('source_hashes') and interpreter.get('member_receipt') and
                interpreter.get('independent_check'), 'incomplete interpreter binding')
        for name, expected in interpreter['source_hashes'].items():
            bound_file(root, name, expected)
        require(interpreter['member_receipt'] in interpreter['source_hashes'] and
                interpreter['independent_check'] in interpreter['source_hashes'],
                'member or independent check is outside source custody')
    review = adoption.get('independent_review')
    if not review:
        blockers.append('Independent review must verify the class, every CP-01 criterion and the CP-13 model map.')
    else:
        require(review.get('verdict') == 'PASS' and review.get('path') and review.get('sha256'),
                'independent review did not accept the proposed contract')
        bound_file(root, review['path'], review['sha256'])
        reviewed = strict_json.load(root/review['path'])
        expected_digest = reviewed_contract_digest(packet)
        require(review.get('reviewed_contract_digest') == expected_digest and
                reviewed.get('reviewed_contract_digest') == expected_digest and
                reviewed.get('verdict') == 'PASS',
                'independent review does not bind this exact class contract')
        review_sources = reviewed.get('source_bindings', {})
        require({'model/completion/finite_attachment_map.json',
                 'model/native_program/results.json'} <= set(review_sources),
                'independent review must bind the actual theorem map and observer receipt')
        for name, expected in review_sources.items():
            bound_file(root, name, expected)
    if adoption['status'] == 'ADOPTED_AFTER_REVIEW':
        require(not blockers, 'cannot adopt an incomplete or unreviewed contract')
    return {'schema': packet['schema'], 'adoption': adoption['status'],
            'target_observation_rows': len(observations), 'bound_sources': len(bindings),
            'reviewed_contract_digest': reviewed_contract_digest(packet),
            'adoption_blockers': blockers, 'CP00_acceptance_established_by_this_checker': False,
            'CP01_acceptance_established_by_this_checker': False,
            'physical_obligations_discharged': [], 'prospective_claim_allowed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='check proposal, custody and adoption safeguards')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
