"""Renew explicit live custody while preserving the accepted-main semantics.

The November 7 and November 11 authorities, records and checkers are retained
bytes. This successor reverses only the typed r2063 navigation delta, then
checks the accepted-main scientific semantics against bd5e6ba80.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.completion import current_custody_renewal as predecessor
from tools.build_completion_program import integration_source_paths, validate_integration_sources

ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = 'model/completion/custody_authority_r2063.json'
AUTHORITY_SHA256 = '0f47735b000df58b670f9ae4b217f5dcfe8e404868adb3b502f9aa8b37d41f81'
RECORD = 'model/completion/custody_renewal_r2063.json'
REVIEW = 'model/completion/program_class_v2_review.json'
CLASS = 'model/completion/program_class_v2.json'
MAP = 'model/completion/finite_attachment_map.json'
PROGRAM = 'claims/completion_program.json'
INTERFACES = 'model/completion/component_interfaces.json'
TRACKING = 'tracking/premise_register.json'
OWNER = 'model/matched_einstein_limit/results.json'
CHILD = 'model/completion/spectrum_einstein_campaign.json'
SELF = 'model/completion/r2063_custody_renewal.py'
TEST = 'model/completion/test_custody_renewal.py'
CURRENT_TOOL = 'tools/build_completion_program_current.py'
LIVE_INTERFACE_FILES = {
    'model/completion/current_theorem_navigation.py',
    'model/completion/theorem_navigation_current.json',
    'model/completion/native_measurement_parent_adapter.py',
    'model/completion/test_native_measurement_parent_adapter.py',
    'tools/renew_canonical_live_receipts.py',
    'tools/canonical_live_receipt_owners_r2063.json',
    'tools/test_renew_canonical_live_receipts.py',
}
LIVE_CANONICAL_TARGETS = {
    'model/completion/source_flavor_comparison.json',
    'model/defect_current_bridge/results.json',
    'model/mcpr_native_memory_bridge/results.json',
    'model/common_source/chiral_candidate/results.json',
    'model/common_source/native_repair_mass_results.json',
    'model/completion/campaign.json',
    'model/majorana_defect_bridge/results.json',
    'model/native_flavor_orbit/results.json',
    'model/native_program/results.json',
    'model/particle_comparison/results.json',
}
INTEGRATION_CONFIGURATIONS = {
    'model/pair_cyclic_control/campaign.json',
    'model/recovery_elision/campaign.json',
    'model/biased_code_universality/campaign.json',
}
DATE = '2026-11-11'
BASE = 'bd5e6ba80bb3e24c09c6e2289e9b72a627a039fa'
ADDITIONS = {INTERFACES: {SELF, AUTHORITY, CURRENT_TOOL, predecessor.SELF, predecessor.AUTHORITY},
             MAP: set(), REVIEW: {SELF, AUTHORITY, RECORD, CURRENT_TOOL,
                                  predecessor.SELF, predecessor.AUTHORITY, predecessor.RECORD}}
ADDITIONS[INTERFACES] |= LIVE_INTERFACE_FILES
ADDITIONS[REVIEW] |= LIVE_INTERFACE_FILES
INTERFACE_MAINTENANCE_DEPENDENCIES = {
    'model/common_source/chiral_candidate/hierarchy_controls.py',
    'model/common_source/chiral_candidate/history.py',
    'model/common_source/chiral_candidate/run.py',
    'model/common_source/chiral_candidate/spacetime_attachment.py',
    'model/common_source/hierarchy_diagnosis/__init__.py',
    'model/common_source/hierarchy_diagnosis/audit.py',
    'model/common_source/hierarchy_diagnosis/charged_sector.py',
    'model/common_source/hierarchy_diagnosis/collective_response.py',
    'model/common_source/hierarchy_diagnosis/native_covariance.py',
    'model/common_source/hierarchy_diagnosis/native_exterior_attachment.py',
    'model/common_source/hierarchy_diagnosis/native_family_dynamics.py',
    'model/common_source/hierarchy_diagnosis/nonlinear_shapes.py',
    'model/common_source/hierarchy_diagnosis/quadratic_response.py',
    'model/common_source/mass_response/dispersion.py',
    'model/common_source/mass_response/experiment.py',
    'model/common_source/mass_response/gap_certificate.py',
    'model/common_source/mass_response/hubbard.py',
    'model/common_source/mass_response/native_attachment.py',
    'model/common_source/mass_response/run.py',
    'model/common_source/particle_identification/clifford.py',
    'model/common_source/particle_identification/kinetic_response.py',
    'model/common_source/particle_identification/response.py',
    'model/common_source/particle_identification/run.py',
    'model/common_source/replay.py',
    'model/exchange_sector/many_body.py',
    'model/federated_observer/scaling.py',
    'model/gauge_defect_transport/__init__.py',
    'model/gauge_defect_transport/transport.py',
    'model/native_causal_dataflow/__init__.py',
    'model/native_causal_dataflow/dataflow.py',
    'model/native_coexistent_effects/__init__.py',
    'model/native_coexistent_effects/menus.py',
    'model/native_coexistent_effects/packet.py',
    'model/native_correlation_preparation/__init__.py',
    'model/native_correlation_preparation/packet.py',
    'model/native_correlation_preparation/preparation.py',
    'model/native_correlation_sectors/__init__.py',
    'model/native_correlation_sectors/packet.py',
    'model/native_correlation_sectors/records.py',
    'model/native_correlation_sectors/sectors.py',
    'model/native_effect_coverage/__init__.py',
    'model/native_effect_coverage/coverage.py',
    'model/native_effect_coverage/detector.py',
    'model/native_effect_coverage/packet.py',
    'model/native_effect_synthesis/__init__.py',
    'model/native_effect_synthesis/intervals.py',
    'model/native_effect_synthesis/packet.py',
    'model/native_effect_synthesis/synthesis.py',
    'model/native_energy_identification/__init__.py',
    'model/native_energy_identification/energy.py',
    'model/native_energy_identification/packet.py',
    'model/native_energy_measurement/__init__.py',
    'model/native_energy_measurement/measurement.py',
    'model/native_energy_measurement/packet.py',
    'model/native_energy_measurement/reference.py',
    'model/native_energy_measurement/spectral.py',
    'model/native_exchange_clock/__init__.py',
    'model/native_exchange_clock/exchange.py',
    'model/native_exchange_clock/packet.py',
    'model/native_fair_local_preparation/__init__.py',
    'model/native_fair_local_preparation/coupling.py',
    'model/native_fair_local_preparation/packet.py',
    'model/native_fair_local_preparation/preparation.py',
    'model/native_fair_local_preparation/witness.py',
    'model/native_fixed_algebra/occupation_service.py',
    'model/native_interface/tomography.py',
    'model/native_local_response/__init__.py',
    'model/native_local_response/bounds.py',
    'model/native_local_response/cone.py',
    'model/native_local_response/packet.py',
    'model/native_noncentral_reference/__init__.py',
    'model/native_noncentral_reference/packet.py',
    'model/native_noncentral_reference/reference.py',
    'model/native_phase_experiment/__init__.py',
    'model/native_phase_experiment/experiment.py',
    'model/native_phase_experiment/independent.py',
    'model/native_phase_experiment/lineage.py',
    'model/native_phase_experiment/provenance.py',
    'model/native_phase_experiment/verify.py',
    'model/native_quantum_gauss/__init__.py',
    'model/native_quantum_gauss/algebra.py',
    'model/native_record_causality/__init__.py',
    'model/native_record_causality/causality.py',
    'model/native_record_response/__init__.py',
    'model/native_record_response/response.py',
    'model/native_record_response/short_records.py',
    'model/native_record_response/transfer.py',
    'model/native_sector_readout/__init__.py',
    'model/native_sector_readout/packet.py',
    'model/native_sector_readout/readout.py',
    'model/native_transfer_control/__init__.py',
    'model/native_transfer_control/control.py',
    'model/native_transfer_control/packet.py',
    'model/native_transfer_control/spectral.py',
    'model/native_vacuum_response/__init__.py',
    'model/native_vacuum_response/packet.py',
    'model/native_vacuum_response/response.py',
    'model/particle_gauge_band/__init__.py',
    'model/particle_gauge_band/quantum_link_completion.py',
    'model/protected_reference/__init__.py',
    'model/protected_reference/gauge_join.py',
    'model/protected_reference/reference.py',
    'model/protected_reference/reversal.py',
    'model/weighted_source_family/__init__.py',
    'model/weighted_source_family/family.py',
    'tools/replay_model_receipts.py',
}
ADDITIONS[INTERFACES] |= INTERFACE_MAINTENANCE_DEPENDENCIES
ADDITIONS[REVIEW] |= INTERFACE_MAINTENANCE_DEPENDENCIES

CURRENT_BINDINGS = {INTERFACES, MAP, OWNER, CHILD, SELF, TEST, AUTHORITY, CURRENT_TOOL,
                    'code/a5_closure/manifests/source_current_order_sensitive_inventory.json'}
CURRENT_BINDINGS |= LIVE_INTERFACE_FILES
RENEWABLE_TARGETS = {INTERFACES, MAP, OWNER, CHILD, SELF, TEST, AUTHORITY, CURRENT_TOOL, RECORD,
                     'Lean/lakefile.lean', 'tools/build_completion_program.py', 'tools/test_completion_program.py',
                     'model/completion/README.md'}
RENEWABLE_TARGETS |= LIVE_CANONICAL_TARGETS | LIVE_INTERFACE_FILES | INTERFACE_MAINTENANCE_DEPENDENCIES | {'model/completion/test_theorem_map.py'}
RENEWABLE_TARGETS |= INTEGRATION_CONFIGURATIONS


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_authority(root=ROOT):
    require(file_digest(Path(root) / AUTHORITY) == AUTHORITY_SHA256, 'immutable custody authority changed')
    value = strict_json.load(Path(root) / AUTHORITY)
    require(value['immutable_base_commit'] == BASE, 'immutable accepted-main baseline changed')
    require(value['renewable_targets'] == sorted(RENEWABLE_TARGETS), 'reviewed renewal target boundary changed')
    return value


def at_path(value, path):
    for key in path:
        value = value[key]
    return value


def program_pointers(program, authority=None):
    authority = load_authority() if authority is None else authority
    rows = []
    def walk(value, path=()):
        if isinstance(value, dict):
            if 'path' in value and 'sha256' in value:
                rows.append((path, value['path']))
            for key in sorted(value):
                walk(value[key], path + (key,))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, path + (index,))
    walk(program)
    schema = {'locations': len(rows), 'targets': len({p for _, p in rows}), 'sha256': digest(rows)}
    require(schema == authority['program_pointer_schema'], 'immutable live pointer schema changed')
    return rows


def semantic_bodies(documents, authority):
    program = copy.deepcopy(documents[PROGRAM])
    for path, _ in program_pointers(program, authority):
        del at_path(program, path)['sha256']
    finite = {k: copy.deepcopy(v) for k, v in documents[MAP].items() if k != 'source_bindings'}
    for key in ('claim_count', 'lean_module_count', 'contract_coverage_sha256'):
        finite['existing_navigation_coverage'][key] = '<current navigation>'
    klass = copy.deepcopy(documents[CLASS])
    reviewed = copy.deepcopy(documents[REVIEW])
    transition = authority['permitted_class_interpreter_custody_successor']
    target = transition['interpreter_receipt']
    require(target == 'model/native_program/results.json', 'class interpreter renewal target changed')
    require(klass['interpreter_binding']['source_hashes'][target] == transition['current_interpreter_receipt_sha256'],
            'current class interpreter receipt changed')
    from model.completion.program_class_v2 import reviewed_contract_digest
    current_digest = reviewed_contract_digest(klass)
    require(current_digest == transition['current_reviewed_contract_digest']
            and klass['adoption']['independent_review']['reviewed_contract_digest'] == current_digest
            and reviewed['reviewed_contract_digest'] == current_digest,
            'current class fingerprint changed outside its interpreter custody')
    require(reviewed.pop('current_interpreter_custody_successor', None) == transition
            and transition['old_review_signed_current_bytes'] is False
            and transition['scientific_clauses_and_gate_statuses_exactly_preserved'] is True,
            'current interpreter custody cannot rewrite original review attribution')
    klass['interpreter_binding']['source_hashes'][target] = transition['retained_interpreter_receipt_sha256']
    klass['adoption']['independent_review']['reviewed_contract_digest'] = transition['retained_reviewed_contract_digest']
    reviewed['reviewed_contract_digest'] = transition['retained_reviewed_contract_digest']
    klass['adoption']['independent_review']['sha256'] = '<current review>'
    return {
        'program': program,
        'finite': finite,
        'class': klass,
        'interfaces': {k: v for k, v in documents[INTERFACES].items() if k != 'source_bindings'},
        'review': {k: v for k, v in reviewed.items()
                   if k not in ('source_bindings', 'current_source_custody_renewal')},
    }


def tracking_transition(register, authority):
    """Reverse the exact reviewed r2063 delta, then apply the retained authority.

    Only descriptions, listed evidence links and the three explicit assumption
    unions change. Every other byte of semantic JSON reaches the original
    bd5e6ba80 conservation check, including all canonical discharge statuses.
    """
    result = copy.deepcopy(register)
    transition = authority['permitted_r2063_navigation_transition']
    for change in transition['groups']:
        group = result['model_supply']['groups'][change['group_index']]
        require(group['id'] == change['group_id'], 'r2063 navigation group identity changed')
        suffix = change['constructed_append']
        require(bool(suffix) and group['constructed'].endswith(suffix),
                'unreviewed r2063 navigation description change')
        group['constructed'] = group['constructed'][:-len(suffix)]
        for addition in change['source_additions']:
            require(group['sources'].count(addition) == 1, 'r2063 navigation source append changed')
            group['sources'].remove(addition)
        replacement = change.get('remaining_replacement')
        if replacement is not None:
            require(group['remaining'] == replacement['after'], 'unreviewed r2063 reference-domain change')
            group['remaining'] = replacement['before']
    additions = transition['assumptions_added']
    for index in transition['binding_indices']:
        assumptions = result['model_supply']['contract_coverage']['bindings'][index]['assumptions']
        require(all(assumptions.count(a) == 1 for a in additions), 'r2063 navigation assumption union changed')
        result['model_supply']['contract_coverage']['bindings'][index]['assumptions'] = [
            a for a in assumptions if a not in additions]
    predecessor.tracking_transition(result, authority)
    coverage = register['model_supply']['contract_coverage']
    return {'claim_count': len(coverage['claim_ids']), 'lean_module_count': len(coverage['lean_files']),
            'contract_coverage_sha256': digest(coverage)}


def authenticated_documents(root, authority):
    documents = {p: strict_json.load(root / p) for p in (PROGRAM, MAP, CLASS, REVIEW, INTERFACES)}
    for name, value in semantic_bodies(documents, authority).items():
        require(digest(value) == authority['baseline_semantics'][name],
                'accepted-main scientific semantics changed: ' + name)
    for name, expected in authority['historical_file_sha256'].items():
        require(file_digest(root / name) == expected, 'historical evidence changed: ' + name)
    validate_canonical_live_inputs(root, authority)
    validate_integration_configurations(root, authority)
    tracking = strict_json.load(root / TRACKING)
    navigation = tracking_transition(tracking, authority)
    owner = strict_json.load(root / OWNER)
    require(digest({k: v for k, v in owner.items() if k != 'source_bindings'})
            == authority['matched_einstein_scientific_sha256'], 'owning Einstein scientific payload changed')
    child = strict_json.load(root / CHILD)
    row = next(r for r in child['component_results'] if r['id'] == 'matched_einstein_limit')
    row['evidence'][OWNER] = '<live owning receipt>'
    require(digest(child) == authority['live_child_evidence']['semantic_sha256'],
            'spectrum-Einstein campaign changed outside its one live owning receipt')
    return documents, navigation


def validate_integration_configurations(root, authority):
    rows = authority['permitted_integration_configuration_renewals']
    require(set(rows) == INTEGRATION_CONFIGURATIONS, 'integration configuration inventory changed')
    for path, row in rows.items():
        require(file_digest(root / row['retained_path']) == row['retained_sha256'],
                'retained integration configuration changed: ' + path)
        require(file_digest(root / path) == row['replay_verified_sha256'],
                'current integration configuration changed: ' + path)
        original = strict_json.load(root / row['retained_path'])
        current = strict_json.load(root / path)
        require(set(row['allowed_source_pin_targets']) == {
            'model/completion/campaign.json', 'model/native_program/results.json'},
            'integration source-pin transition expanded')
        for target in row['allowed_source_pin_targets']:
            original['source_bindings'][target] = file_digest(root / target)
        require(digest(original) == digest(current), 'integration scientific configuration changed: ' + path)
        validate_integration_sources(root, current)


def validate_canonical_live_inputs(root, authority):
    """Check the retained science and explicitly executed current input bytes.

    This is a custody check of pinned replay results, not a fresh execution or
    an admission of arbitrary changes to a source implementation.
    """
    from tools import renew_canonical_live_receipts as live
    rows = authority['permitted_canonical_live_renewals']
    require(set(rows) == LIVE_CANONICAL_TARGETS, 'canonical live renewal inventory changed')
    for path, row in rows.items():
        require(file_digest(root / row['retained_path']) == row['retained_sha256'],
                'retained canonical input changed: ' + path)
        require(file_digest(root / path) == row['replay_verified_sha256'],
                'current canonical input differs from executed replay: ' + path)
        old = strict_json.load(root / row['retained_path'])
        current = strict_json.load(root / path)
        if path == 'model/completion/campaign.json':
            # This historical campaign configuration uses a list of source
            # declarations, not a replay hash map. Only its one control pointer
            # was explicitly made live by the retained clock contract.
            pointer, = row['reviewed_receipt_pointers']
            require(file_digest(root / pointer['retained_path']) == pointer['retained_sha256']
                    and at_path(old, pointer['json_path']) == pointer['retained_sha256'],
                    'retained campaign control authority changed')
            at_path(old, pointer['json_path'][:-1])[pointer['json_path'][-1]] = file_digest(root / pointer['target'])
            require(digest(old) == digest(current), 'campaign changed outside admitted control custody')
            continue
        projected = live.carry_metadata(path, old, current, row.get('reviewed_source_scope_append'),
                                        row.get('reviewed_receipt_pointers', ()), root=root)
        require(live.canonical(projected) == live.canonical(current),
                'canonical input scientific payload changed: ' + path)
        for here, relative, expected, target in live.bindings(root, path, current):
            require(file_digest(target) == expected, 'current canonical source changed: ' + path + ':' + relative)
    require(set(authority['maintenance_interface_dependency_append']) == INTERFACE_MAINTENANCE_DEPENDENCIES,
            'maintenance interface dependency inventory changed')
    for path, expected in authority['maintenance_interface_dependency_append'].items():
        require(file_digest(root / path) == expected, 'maintenance interface dependency changed: ' + path)
    for path, expected in authority['reviewed_live_interface_sha256'].items():
        require(file_digest(root / path) == expected, 'reviewed custody interface changed: ' + path)
    if root == ROOT:
        from model.completion import current_theorem_navigation
        current_theorem_navigation.validate()


def check_binding_scope(document, path, authority):
    bindings = document['source_bindings']
    added = ADDITIONS[path]
    require(digest(sorted(set(bindings) - added)) == authority['baseline_binding_path_sha256'][path],
            'live binding path coverage changed: ' + path)
    static = {p: v for p, v in bindings.items() if p not in RENEWABLE_TARGETS}
    require(digest(static) == authority['baseline_static_binding_sha256'][path],
            'immutable nonrenewable source custody changed: ' + path)


def check_bindings(document, path, root, authority):
    check_binding_scope(document, path, authority)
    bindings = document['source_bindings']
    require(ADDITIONS[path] <= set(bindings), 'live binding path coverage changed: ' + path)
    for name, expected in bindings.items():
        source = (root / name).resolve()
        require(source.is_relative_to(root) and source.is_file(), 'missing or escaping live source: ' + name)
        require(file_digest(source) == expected, 'live source binding changed: ' + name)


def validate(packet=None, root=ROOT):
    root = Path(root).resolve()
    authority = load_authority(root)
    packet = strict_json.load(root / RECORD) if packet is None else packet
    require(packet.get('schema') == 'oph.current_custody_renewal.v3' and packet.get('date') == DATE,
            'wrong successor custody record')
    require(packet.get('authority') == {'path': AUTHORITY, 'sha256': AUTHORITY_SHA256, 'base_commit': BASE},
            'immutable successor authority changed')
    require(packet.get('new_mathematical_acceptance') is False and packet.get('gate_statuses_changed') is False
            and packet.get('canonical_physical_premises_discharged') == [], 'custody cannot promote scientific acceptance')
    require(packet.get('historical_evidence_rewritten') is False,
            'custody cannot rewrite historical evidence')
    require(packet.get('program_pointer_schema') == authority['program_pointer_schema'],
            'immutable program pointer authority changed')
    require(packet.get('classification_counts') == authority['classification_counts'], 'classification authority changed')
    documents, navigation = authenticated_documents(root, authority)
    require(packet.get('navigation') == navigation, 'successor navigation certificate changed')
    require(all(documents[MAP]['existing_navigation_coverage'][k] == v for k, v in navigation.items()),
            'finite-map live navigation changed')
    require(packet.get('conserved_semantics') == authority['baseline_semantics'], 'immutable conservation authority changed')
    require(packet.get('live_child_evidence') == authority['live_child_evidence'], 'live child evidence authority changed')
    child = strict_json.load(root / CHILD)
    require(next(r for r in child['component_results'] if r['id'] == 'matched_einstein_limit')['evidence'][OWNER]
            == file_digest(root / OWNER), 'live campaign owning receipt changed')
    for name in (INTERFACES, MAP, REVIEW):
        check_bindings(documents[name], name, root, authority)
    validate_integration_sources(root, documents[INTERFACES])
    for path, target in program_pointers(documents[PROGRAM], authority):
        require(at_path(documents[PROGRAM], path)['sha256'] == file_digest(root / target),
                'live program custody target changed: ' + target)
    require(set(packet.get('current_bindings', {})) == CURRENT_BINDINGS, 'current binding coverage changed')
    for path, expected in packet['current_bindings'].items():
        require(file_digest(root / path) == expected, 'current custody binding changed: ' + path)
    expected_pointer = {'date': DATE, 'path': RECORD, 'role': 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE',
                        'sha256': file_digest(root / RECORD)}
    require(documents[REVIEW]['current_source_custody_renewal'] == expected_pointer,
            'review does not bind successor custody')
    require(documents[CLASS]['adoption']['independent_review']['sha256'] == file_digest(root / REVIEW),
            'adoption custody does not bind current review')
    return {'date': DATE, 'immutable_base_commit': BASE, 'classification_counts': authority['classification_counts'],
            'program_pointer_schema': authority['program_pointer_schema'], 'navigation': navigation,
            'new_mathematical_acceptance': False, 'scientific_status_changed': False,
            'historical_evidence_rewritten': False}


def write_json(root, name, value):
    (root / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def renew(root=ROOT):
    """Explicit maintenance operation; validation alone never repairs custody."""
    root = Path(root).resolve()
    authority = load_authority(root)
    documents, navigation = authenticated_documents(root, authority)
    # A live renewal is not permission to repin arbitrary altered science.
    # The unchanged dependency set is checked before any output is written.
    for name in (INTERFACES, MAP, REVIEW):
        check_binding_scope(documents[name], name, authority)
        changed = {p for p, expected in documents[name]['source_bindings'].items()
                   if file_digest(root / p) != expected}
        require(changed <= RENEWABLE_TARGETS,
                'unreviewed source change cannot be renewed: ' + ', '.join(sorted(changed - RENEWABLE_TARGETS)))
    # Execute the owning exact producer, not just its hash calculator.
    from model.matched_einstein_limit import packet as producer
    require(root == ROOT, 'producer-backed renewal requires the actual repository root')
    fresh = producer.report()
    require(digest({k: v for k, v in fresh.items() if k != 'source_bindings'})
            == authority['matched_einstein_scientific_sha256'], 'fresh owning science changed')
    (root / OWNER).write_text(producer.canonical(fresh))
    child = strict_json.load(root / CHILD)
    next(r for r in child['component_results'] if r['id'] == 'matched_einstein_limit')['evidence'][OWNER] = file_digest(root / OWNER)
    write_json(root, CHILD, child)
    for key, value in navigation.items():
        documents[MAP]['existing_navigation_coverage'][key] = value
    for name in (INTERFACES, MAP):
        keys = list(documents[name]['source_bindings'])
        keys += sorted(ADDITIONS[name] - set(keys))
        if name == INTERFACES:
            require(integration_source_paths(root, documents[name]) <= set(keys),
                    'new transitive dependency requires explicit custody review')
        documents[name]['source_bindings'] = {p: file_digest(root / p) for p in keys}
        write_json(root, name, documents[name])
    record = {
        'schema': 'oph.current_custody_renewal.v3', 'date': DATE,
        'authority': {'path': AUTHORITY, 'sha256': AUTHORITY_SHA256, 'base_commit': BASE},
        'scope': 'Live source custody and the exact owner-approved r2063 navigation delta; no historical proof replay or new gate acceptance.',
        'new_mathematical_acceptance': False, 'gate_statuses_changed': False,
        'canonical_physical_premises_discharged': [], 'historical_evidence_rewritten': False,
        'conserved_semantics': authority['baseline_semantics'],
        'live_child_evidence': authority['live_child_evidence'],
        'classification_counts': authority['classification_counts'],
        'program_pointer_schema': authority['program_pointer_schema'], 'navigation': navigation,
        'current_bindings': {p: file_digest(root / p) for p in sorted(CURRENT_BINDINGS)},
        'owning_producer': {'module': 'model.matched_einstein_limit.packet',
                            'scientific_payload_sha256': authority['matched_einstein_scientific_sha256'],
                            'scientific_payload_exactly_preserved': True,
                            'fresh_Lean_kernel_execution_claimed_by_this_record': False},
    }
    write_json(root, RECORD, record)
    review = documents[REVIEW]
    keys = list(review['source_bindings'])
    keys += sorted(ADDITIONS[REVIEW] - set(keys))
    review['source_bindings'] = {p: file_digest(root / p) for p in keys}
    review['current_source_custody_renewal'] = {
        'date': DATE, 'path': RECORD, 'role': 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE',
        'sha256': file_digest(root / RECORD)}
    write_json(root, REVIEW, review)
    documents[CLASS]['adoption']['independent_review']['sha256'] = file_digest(root / REVIEW)
    write_json(root, CLASS, documents[CLASS])
    for path, target in program_pointers(documents[PROGRAM], authority):
        at_path(documents[PROGRAM], path)['sha256'] = file_digest(root / target)
    write_json(root, PROGRAM, documents[PROGRAM])
    return validate(root=root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--renew', action='store_true', help='explicit producer-backed live custody maintenance')
    args = parser.parse_args()
    print(json.dumps(renew() if args.renew else validate(), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
