#!/usr/bin/env python3
"""Current completion interface with an explicit joint dependency owner.

The original builder remains byte-identical because scientific producer
closures import its source-custody helper. This successor preserves every
original load_program check, adding only the typed dependency contract owner.
Its retained rendering helpers do not execute any scientific producer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

if __package__:
    from tools import build_completion_program as retained
else:  # direct execution from tools/
    import build_completion_program as retained

ROOT, PROGRAM, HYPOTHESIS, OUTPUT = retained.ROOT, retained.PROGRAM, retained.HYPOTHESIS, retained.OUTPUT
require, local_file, bound_evidence = retained.require, retained.local_file, retained.bound_evidence
validate_integration_sources, validate_next_wave = retained.validate_integration_sources, retained.validate_next_wave
render, render_issue, issue_link, readiness = retained.render, retained.render_issue, retained.issue_link, retained.readiness
JOINT_DEPENDENCY_ID = 'OPH-JOINT-FIELD-SPECTRUM-CONTRACT-v1'
JOINT_DEPENDENCY_PATH = 'model/particle_joint_reconstruction/dependencies.json'
JOINT_DEPENDENCY_CHECK = 'model/particle_joint_reconstruction/test_dependencies.py'
JOINT_DEPENDENCY_GATES = {'CP-13', 'CP-04', 'CP-05', 'CP-06', 'CP-07', 'CP-09'}


def joint_dependency_owner(root: Path, stages: list) -> dict:
    """Typed owner for the existing dependency graph, not a finite campaign.

    Custody/scope validation remains separate from its mathematical producer.
    No unknown campaign identifier, new gate mapping or acceptance is admitted.
    """
    anchors = [(stage['id'], row) for stage in stages for row in stage.get('progress_evidence', [])
               if row.get('id') == 'joint_field_spectrum_dependency_contract']
    require(len(anchors) == 1 and anchors[0][0] == 'CP-13', 'joint dependency owner anchor missing or moved')
    anchor = anchors[0][1]
    require(anchor.get('campaign_id') == JOINT_DEPENDENCY_ID
            and anchor.get('path') == JOINT_DEPENDENCY_PATH
            and anchor.get('independent_check') == JOINT_DEPENDENCY_CHECK,
            'joint dependency owner identity changed')
    bound_evidence(root, [anchor])
    graph = json.loads(local_file(root, JOINT_DEPENDENCY_PATH).read_text())
    require(graph.get('schema') == 'oph.joint_reconstruction_dependencies.v1'
            and graph.get('primary_program') == 'oph.open_wilson_gauss_observer.v1',
            'joint dependency schema or primary source changed')
    ids = ('fixed_program', 'physical_preparation', 'local_records', 'joint_refinement',
           'relativistic_field', 'chiral_species_attachment', 'particle_observables',
           'independent_comparison', 'pure_gauge_gap', 'geometric_gravity')
    nodes = graph.get('nodes', [])
    require([r.get('id') for r in nodes] == list(ids)
            and [r.get('number') for r in nodes] == list(range(1, 11)),
            'joint dependency node coverage changed')
    finite_domains = (
        'The selected open Wilson family, with sixty CAR modes per cell, product link-word cutoff and the original energy zero.',
        'The exact Gauss-polynomial state of each selected finite member and its retained canonical conjugate reference.',
        'The actual selected Wilson action and prepared state, the complete origin-cell boundary-charge instrument, and a remote physical plaquette phase on the same input and clock.',
    )
    for index, node in enumerate(nodes):
        status = 'finite_constructed' if index < 3 else 'conditional_criterion' if index == 8 else 'physical_proof_required'
        dependencies = [] if index == 0 else ['relativistic_field'] if index >= 8 else [ids[index - 1]]
        require(node.get('status') == status and node.get('depends_on') == dependencies,
                'joint dependency status or required edge changed')
        require(node.get('domain') and node.get('acceptance') and node.get('evidence'),
                'joint dependency scope missing')
        if index < 3:
            require(node['domain'] == finite_domains[index], 'joint finite source domain changed')
        for evidence in node['evidence']:
            local_file(root, evidence.get('path'))
            require(bool(evidence.get('scope')), 'joint dependency evidence scope missing')
    extension = graph.get('conditional_candidate', {})
    require(extension.get('id') == 'oph.conditional_chiral_gauge_higgs_extension.v1'
            and extension.get('status') == 'conditional_extension'
            and all(extension.get(key) is False for key in (
                'discharges_primary_physical_dependencies', 'calibrated_tree_coefficients_are_independent_predictions',
                'full_theory_positive_global_gap_claimed')), 'joint conditional extension promoted')
    require(nodes[8].get('required_shared_nodes') == list(ids[:5])
            and nodes[8].get('required_for_full_particle_spectrum') is False
            and nodes[8].get('implies_full_theory_positive_global_gap') is False
            and nodes[9].get('required_for_full_particle_spectrum') is False
            and nodes[9].get('proper_geometry_or_Einstein_limit_proved') is False,
            'joint dependency physical boundary changed')
    assumptions = {a for node in nodes for a in node.get('conditional_assumptions', [])}
    require(assumptions == {
        'chiral_measure_and_mirror_control', 'common_vacuum_reference_clock', 'extension_carrier_and_action',
        'neutrino_EFT_completion_and_running', 'nonzero_connected_interaction', 'physical_current_spin_detector_dictionary',
        'pure_YM_dense_decay_and_noncollapse', 'pure_YM_source_and_sector_attachment',
        'radiative_infrared_resonance_matching', 'renormalized_four_dimensional_correlations'},
        'joint conditional assumption coverage changed')
    maps = graph.get('requirement_maps', {})
    expected_maps = {key: {str(i) for i in range(1, 11)} for key in
                     ('field_common_theory_ten', 'native_full_spectrum_ten', 'defect_identification_ten')}
    expected_maps['yang_mills_Y1_Y6'] = {f'Y{i}' for i in range(1, 7)}
    require(set(maps) == set(expected_maps), 'joint requirement-map coverage changed')
    for key, targets in maps.items():
        require(set(targets) == expected_maps[key]
                and all(values and all(v in ids or v.startswith('outside_target:') for v in values)
                        for values in targets.values()), 'joint requirement destination changed')
    for stage in stages:
        require(not any(row.get('campaign_id') == JOINT_DEPENDENCY_ID for row in stage.get('accepted_evidence', [])),
                'joint dependency contract cannot accept a completion gate')
    return {'id': JOINT_DEPENDENCY_ID, 'supports': sorted(JOINT_DEPENDENCY_GATES),
            'kind': 'DEPENDENCY_CONTRACT_NOT_COMPONENT_CAMPAIGN'}



def load_program(root: Path = ROOT, *, allow_unlinked: bool = False) -> tuple[dict, dict]:
    root = root.resolve()
    program = json.loads(local_file(root, PROGRAM).read_text())
    hypothesis = json.loads(local_file(root, HYPOTHESIS).read_text())
    require(program.get("schema") == "oph.completion_program.v1", "wrong program schema")
    require(hypothesis.get("schema") == "oph.executable_model_hypothesis.v1", "wrong hypothesis schema")
    require(program.get("hypothesis") == HYPOTHESIS, "wrong hypothesis parent")
    require(hypothesis.get("claim_class") == "working_hypothesis"
            and hypothesis.get("status") == "PROPOSED_NOT_ESTABLISHED",
            "hypothesis cannot be silently promoted to established physics")
    for flag in ("global_positive_mass_gap_claimed", "modular_flow_automatically_physical_time",
                 "all_particles_stable_defects_established", "einstein_premises_realized"):
        require(hypothesis.get(flag) is False, f"unsupported hypothesis promotion: {flag}")
    require(hypothesis.get("common_process_required") is True
            and hypothesis.get("comparison_count_added") == 0, "wrong joint comparison scope")
    local_file(root, hypothesis.get("contract"))
    obligations = hypothesis.get("obligations", [])
    obligation_ids = {row["id"] for row in obligations}
    require(len(obligations) == 5 and obligation_ids == {f"EMH-{i}" for i in range(1, 6)},
            "missing or duplicate hypothesis obligation")
    require(all(row.get("status") == "NOT_ESTABLISHED" for row in obligations),
            "physical obligation promotion needs a separately reviewed contract revision")
    stages = program.get("stages", [])
    ids = [stage["id"] for stage in stages]
    require(len(ids) == len(set(ids)) and set(ids) == {f"CP-{i:02d}" for i in range(14)},
            "missing or duplicate completion gate")
    by_id = {stage["id"]: stage for stage in stages}
    campaigns = program.get("component_campaigns", [])
    require(isinstance(campaigns, list), "invalid component campaigns")
    campaign_ids = set()
    for campaign in campaigns:
        identifier = campaign.get("id")
        require(isinstance(identifier, str) and identifier and identifier not in campaign_ids,
                "missing or duplicate component campaign")
        campaign_ids.add(identifier)
        require(campaign.get("scope") in {"FINITE_COMPONENT_FAMILY", "FINITE_COMPONENT_DIAGNOSTICS"}
                and campaign.get("whole_program_frozen") is False,
                "component campaign cannot freeze the whole program")
        require(bool(campaign.get("supports")) and set(campaign["supports"]) <= set(ids),
                "unmapped component campaign")
        bound_evidence(root, [campaign.get("contract")])
        contract = json.loads(local_file(root, campaign["contract"]["path"]).read_text())
        require(contract.get("campaign_id") == identifier
                and contract.get("scope") == campaign["scope"]
                and contract.get("whole_program_frozen") is False
                and contract.get("physical_comparison_count_added") == 0,
                "component campaign contract scope changed")
        if contract.get('schema') == 'oph.integration_simplification.v1':
            validate_integration_sources(root, contract)
        if campaign["scope"] == "FINITE_COMPONENT_DIAGNOSTICS":
            require(contract.get("boundary_status") == "FROZEN_DIAGNOSTIC_CONTRACT"
                    and contract.get("integrated_gate_accepted") is False
                    and contract.get("native_family_changed") is False,
                    "component diagnostics cannot accept gates or change the native family")
    campaign_by_id = {campaign["id"]: campaign for campaign in campaigns}
    dependency_owner = joint_dependency_owner(root, stages)
    require(dependency_owner['id'] not in campaign_by_id, 'dependency owner cannot be recast as a component campaign')
    campaign_by_id[dependency_owner['id']] = dependency_owner
    campaign_ids.add(dependency_owner['id'])
    claims = json.loads(local_file(root, "claims/claim_registry.yaml").read_text())["claims"]
    claim_ids = {row["claim_id"] for row in claims}
    foundations = program.get("foundations", [])
    require(bool(foundations) and len({row["id"] for row in foundations}) == len(foundations),
            "missing or duplicate retained foundations")
    for foundation in foundations:
        require(foundation.get("claim_id") in claim_ids, "foundation lacks a canonical claim")
        require(set(foundation.get("supports", [])) <= set(ids)
                and bool(foundation.get("achieved")) and bool(foundation.get("boundary")),
                "unscoped retained foundation")
        for path in foundation["evidence"]:
            local_file(root, path)
    issues = []
    for stage in stages:
        deps = stage.get("requires", [])
        require(isinstance(deps, list) and len(deps) == len(set(deps))
                and set(deps) <= set(ids) and stage["id"] not in deps,
                f"invalid dependencies: {stage['id']}")
        require(stage.get("status") in {"OPEN", "IN_PROGRESS", "ACCEPTED", "OBSTRUCTED"},
                f"invalid gate status: {stage['id']}")
        require(bool(stage.get("acceptance")) and bool(stage.get("obstruction"))
                and bool(stage.get("deliverable")) and bool(stage.get("evidence")),
                f"incomplete acceptance/obstruction contract: {stage['id']}")
        require(set(stage.get("obligations", [])) <= obligation_ids
                and bool(stage.get("obligations")), "unmapped hypothesis obligation")
        for path in stage["evidence"]:
            local_file(root, path)
        if stage.get("progress_evidence"):
            bound_evidence(root, stage["progress_evidence"])
            require(all(row.get("campaign_id") in campaign_ids
                        and stage["id"] in campaign_by_id[row["campaign_id"]]["supports"]
                        for row in stage["progress_evidence"]),
                    "component progress lacks a registered campaign for this gate")
        if stage["status"] == "ACCEPTED":
            bound_evidence(root, stage.get("accepted_evidence"))
            require(all(by_id[dep]["status"] == "ACCEPTED" for dep in deps),
                    f"gate accepted before predecessors: {stage['id']}")
        issue = stage.get("issue")
        if issue is not None:
            require(type(issue) is int and issue > 0, "invalid GitHub issue")
            issues.append(issue)
        else:
            require(allow_unlinked, f"unlinked completion issue: {stage['id']}")
    require(len(issues) == len(set(issues)), "duplicate completion issue ownership")
    if "next_wave" in program:
        validate_next_wave(root, program["next_wave"], by_id, campaigns)
    visiting, visited = set(), set()

    def visit(id):
        require(id not in visiting, f"completion dependency cycle at {id}")
        if id in visited:
            return
        visiting.add(id)
        for dep in by_id[id]["requires"]:
            visit(dep)
        visiting.remove(id)
        visited.add(id)

    for id in ids:
        visit(id)
    umbrella = program.get("umbrella_issue")
    require((allow_unlinked and umbrella is None)
            or (type(umbrella) is int and umbrella > 0 and umbrella not in issues),
            "invalid umbrella issue")
    frozen = program.get("boundary_status") == "FROZEN"
    require(program.get("boundary_status") in {"NOT_FROZEN", "FROZEN"}, "invalid boundary state")
    require(frozen == (by_id["CP-00"]["status"] == "ACCEPTED"),
            "frozen class requires an accepted boundary contract")
    require(program.get("decision") in {"NOT_DECIDED", "FALSIFIED_PROGRAM", "STOPPED_RESOURCE", "COMPLETED_DOMAIN"},
            "invalid program decision")
    falsified = program.get("decision") == "FALSIFIED_PROGRAM"
    require(program.get("whole_program_falsified") is falsified, "inconsistent whole-program verdict")
    if program["decision"] != "NOT_DECIDED":
        bound_evidence(root, program.get("decision_evidence"))
    if falsified:
        require(frozen, "whole-program falsification requires a frozen quantified class")
        require(program.get("coverage") == "ALL_FROZEN_REALIZATIONS"
                and program.get("decision_kind") in {"INDISPENSABLE_CONTRADICTION", "UNIVERSAL_REALIZATION_NO_GO", "UNAVOIDABLE_EMPIRICAL_EXCLUSION"},
                "missing whole-program coverage or valid falsifier")
    if program["decision"] == "COMPLETED_DOMAIN":
        require(frozen and all(stage["status"] == "ACCEPTED" for stage in stages),
                "completion requires every accepted gate")
    return program, hypothesis



def ledger_view(root: Path = ROOT) -> dict:
    program, hypothesis = load_program(root)
    parents = {}
    for name in (PROGRAM, HYPOTHESIS, hypothesis["contract"], "claims/claim_registry.yaml"):
        parents[name] = hashlib.sha256(local_file(root, name).read_bytes()).hexdigest()
    claims = {row["claim_id"]: row for row in json.loads(local_file(root, "claims/claim_registry.yaml").read_text())["claims"]}
    foundations = []
    for foundation in program["foundations"]:
        row = dict(foundation)
        row["canonical_claim_status"] = claims[row["claim_id"]]["status"]
        row["source_sha256"] = {name: hashlib.sha256(local_file(root, name).read_bytes()).hexdigest()
                                for name in row["evidence"]}
        foundations.append(row)
    return dict(verification="Program structure and source custody only; not scientific gate acceptance or a new physical comparison",
                parents=parents, hypothesis=hypothesis, program=program, retained_foundations=foundations)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--allow-unlinked", action="store_true", help="local drafting only")
    args = parser.parse_args()
    program, hypothesis = load_program(allow_unlinked=args.allow_unlinked)
    expected = render(program, hypothesis)
    if args.check:
        require((ROOT / OUTPUT).read_text() == expected, "completion program view is stale")
        print("Completion program DAG and issue bindings verified; no physical gate promoted")
    else:
        (ROOT / OUTPUT).write_text(expected)
        print(OUTPUT)


if __name__ == "__main__":
    main()
