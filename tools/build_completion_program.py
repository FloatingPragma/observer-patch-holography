#!/usr/bin/env python3
"""Validate and render the executable completion DAG; no scientific replay."""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO_URL = "https://github.com/FloatingPragma/observer-patch-holography"
PROGRAM = "claims/completion_program.json"
HYPOTHESIS = "claims/executable_model_hypothesis.json"
OUTPUT = "docs/COMPLETION_PROGRAM.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def local_file(root: Path, name: str) -> Path:
    require(isinstance(name, str) and bool(name), "missing local artifact")
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()) and path.is_file(),
            f"missing or escaping artifact: {name}")
    return path


def bound_evidence(root: Path, rows: list) -> None:
    require(isinstance(rows, list) and bool(rows), "missing bound decision evidence")
    for row in rows:
        require(isinstance(row, dict), "untyped decision evidence")
        path = local_file(root, row.get("path"))
        require(hashlib.sha256(path.read_bytes()).hexdigest() == row.get("sha256"),
                f"decision evidence custody changed: {row.get('path')}")
        local_file(root, row.get("independent_check"))
        require(bool(row.get("scope")), "missing evidence scope")


def integration_source_paths(root: Path, contract: dict) -> set[str]:
    """Repository Python closure, including function imports and eager parents."""
    document = local_file(root, contract.get('domain_and_proof'))
    paths = {document.relative_to(root).as_posix(), 'model/requirements.txt'}
    paths.update(p.relative_to(root).as_posix() for p in document.parent.glob('*.py'))
    paths.update(contract.get('independent_checks', []))
    receipts = contract.get('source_receipts', [])
    require(isinstance(receipts, list) and bool(receipts),
            'integration source custody requires owning source receipts')
    for name in receipts:
        receipt = json.loads(local_file(root, name).read_text())
        inherited = receipt.get('source_bindings', receipt.get('source_sha256'))
        require(isinstance(inherited, dict) and bool(inherited),
                f'integration source receipt has no source manifest: {name}')
        paths.add(name)
        paths.update(inherited)
    pending = list(paths)
    visited = set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        candidates = [stem+'.py', stem+'/__init__.py']
        parts = stem.split('/')
        candidates.extend('/'.join(parts[:end])+'/__init__.py'
                          for end in range(1, len(parts)))
        for candidate in candidates:
            if (root/candidate).is_file() and candidate not in paths:
                paths.add(candidate)
                pending.append(candidate)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removesuffix('.py').replace('/', '.')
        package = (module.removesuffix('.__init__') if module.endswith('.__init__')
                   else module.rpartition('.')[0])
        add_module(module)
        for node in ast.walk(ast.parse(local_file(root, name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    return paths


def validate_integration_sources(root: Path, contract: dict) -> None:
    paths = integration_source_paths(root, contract)
    bindings = contract.get('source_bindings')
    require(isinstance(bindings, dict) and paths <= set(bindings),
            'integration source custody omits a transitive dependency')
    for name, expected in bindings.items():
        require(hashlib.sha256(local_file(root, name).read_bytes()).hexdigest() == expected,
                f'integration source custody changed: {name}')
    for name in contract['source_receipts']:
        receipt = json.loads(local_file(root, name).read_text())
        inherited = receipt.get('source_bindings', receipt.get('source_sha256'))
        require(all(bindings.get(path) == expected for path, expected in inherited.items()),
                f'integration source custody differs from owning receipt: {name}')
    require(contract.get('integrated_gate_accepted') is False
            and contract.get('source_law_changed') is False,
            'integration simplification cannot accept a gate or change its source law')


def validate_next_wave(root: Path, wave: dict, by_id: dict, campaigns: list) -> None:
    """A work assignment may reuse evidence, but cannot widen its contract."""
    require(isinstance(wave, dict) and bool(wave.get("id")), "missing next-wave assignment")
    require(wave.get("status") == "PLANNED_PROVISIONAL_PARALLEL_WORK"
            and wave.get("scope") == "SOURCE_BOUND_COMPONENT_WORK_WITH_OPEN_ACCEPTANCE_JOINS",
            "next wave cannot silently become accepted evidence")
    require(wave.get("changes_frozen_campaign") is False
            and wave.get("successor_class_frozen") is False
            and wave.get("physical_comparison_count_added") == 0,
            "next wave cannot change a frozen class or add a physical comparison")
    campaign_by_id = {row["id"]: row for row in campaigns}
    interface = wave.get("shared_interface", {})
    require(isinstance(interface, dict), "untyped next-wave shared interface")
    baseline = campaign_by_id.get(interface.get("baseline_campaign_id"))
    require(baseline is not None
            and interface.get("baseline_contract") == baseline["contract"]["path"],
            "next wave lacks the registered common baseline")
    require(interface.get("same_process_required") is True
            and interface.get("unilateral_primitive_changes_allowed") is False
            and isinstance(interface.get("requirements"), list)
            and bool(interface["requirements"])
            and all(isinstance(row, str) and row.strip() for row in interface["requirements"]),
            "next wave lacks a common interface and revision boundary")
    evidence = wave.get("evidence")
    bound_evidence(root, evidence)
    evidence_ids = [row.get("id") for row in evidence]
    require(all(isinstance(identifier, str) and identifier for identifier in evidence_ids)
            and len(evidence_ids) == len(set(evidence_ids)),
            "missing or duplicate next-wave evidence identifier")
    require(all(row.get("campaign_id") in campaign_by_id for row in evidence),
            "next-wave evidence lacks a registered campaign")
    owners = wave.get("owners", [])
    require(isinstance(owners, list) and 3 <= len(owners) <= 5
            and all(isinstance(row, dict) for row in owners)
            and len({row.get("gate") for row in owners}) == len(owners),
            "next wave requires three to five distinct issue owners")
    reviews = wave.get("embedded_reviews", [])
    require(isinstance(reviews, list) and bool(reviews)
            and all(isinstance(row, dict) for row in reviews)
            and len({row.get("gate") for row in reviews}) == len(reviews),
            "missing or duplicate embedded review owner")
    require(not ({row.get("gate") for row in owners} & {row.get("gate") for row in reviews}),
            "embedded review cannot be another primary work assignment")
    for row in owners + reviews:
        gate = row.get("gate")
        require(gate in by_id and row.get("issue") == by_id[gate]["issue"],
                "next-wave issue differs from its canonical gate owner")
        refs = row.get("evidence", [])
        require(isinstance(refs, list) and bool(refs) and set(refs) <= set(evidence_ids),
                "next-wave work lacks bound starting evidence")
    for row in owners:
        require(all(isinstance(row.get(field), str) and row[field].strip()
                    for field in ("title", "why_now", "deliverable", "may_start_now", "join_condition")),
                "next-wave owner lacks a concrete deliverable or parallel boundary")
    require(all(isinstance(row.get("role"), str) and row["role"].strip() for row in reviews),
            "embedded review lacks its scope")
    require(all(isinstance(wave.get(field), str) and wave[field].strip()
                for field in ("rationale", "acceptance_policy", "deferred_reason")),
            "next-wave rationale or acceptance policy missing")
    deferred = wave.get("deferred_gates", [])
    require(isinstance(deferred, list) and len(deferred) == len(set(deferred))
            and set(deferred) <= set(by_id)
            and not (set(deferred) & {row["gate"] for row in owners + reviews}),
            "invalid deferred next-wave gate")


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


def issue_link(number: int | None, label: str | None = None) -> str:
    return f"[{label or '#' + str(number)}]({REPO_URL}/issues/{number})" if number else "Issue pending"


def readiness(stage: dict, program: dict) -> str:
    by_id = {row["id"]: row for row in program["stages"]}
    if stage["status"] in {"ACCEPTED", "OBSTRUCTED"}:
        return stage["status"]
    missing = [dep for dep in stage["requires"] if by_id[dep]["status"] != "ACCEPTED"]
    if stage["status"] == "IN_PROGRESS":
        return ("In progress; acceptance awaiting " + ", ".join(missing)
                if missing else "In progress")
    return "Awaiting " + ", ".join(missing) if missing else "Ready to start"


def render_next_wave(program: dict) -> list[str]:
    wave = program.get("next_wave")
    if not wave:
        return []
    evidence = {row["id"]: row for row in wave["evidence"]}
    owner_count = {3: "three", 4: "four", 5: "five"}[len(wave["owners"])]
    lines = [f"## Next {owner_count} parallel workstreams", "",
             f"Assignment: `{wave['id']}`. This plans provisional work; it does not "
             "select or freeze a successor law, change the retained campaign or accept a gate.", "",
             wave["rationale"], "",
             "| Owner | Next concrete result | Why this comes next |",
             "| --- | --- | --- |"]
    for owner in wave["owners"]:
        lines.append(f"| {issue_link(owner['issue'])} / {owner['gate']}: {owner['title']} "
                     f"| {owner['deliverable']} | {owner['why_now']} |")
    lines += ["", "Parallel starting work and the required joins:", ""]
    for owner in wave["owners"]:
        refs = ", ".join(f"[{identifier.replace('_', ' ')}](../{evidence[identifier]['path']})"
                         for identifier in owner["evidence"])
        lines.append(f"- **{owner['gate']} / {issue_link(owner['issue'])}:** "
                     f"{owner['may_start_now']} {owner['join_condition']} Starting evidence: {refs}.")
    lines += ["", "Shared source interface:", ""]
    lines += [f"- {requirement}" for requirement in wave["shared_interface"]["requirements"]]
    lines += ["", wave["acceptance_policy"], "", "Review stays inside these workstreams:", ""]
    lines += [f"- {review['gate']} / {issue_link(review['issue'])}: {review['role']}"
              for review in wave["embedded_reviews"]]
    lines += ["", wave["deferred_reason"], ""]
    return lines


def mermaid(program: dict) -> str:
    lines = ["```mermaid", "flowchart TD"]
    for stage in program["stages"]:
        label = stage["id"] + ": " + stage["title"]
        lines.append(f"    {stage['id'].replace('-', '')}[{json.dumps(label)}]")
    for stage in program["stages"]:
        for dep in stage["requires"]:
            lines.append(f"    {dep.replace('-', '')} --> {stage['id'].replace('-', '')}")
    for stage in program["stages"]:
        if stage.get("issue"):
            lines.append(f"    click {stage['id'].replace('-', '')} \"{REPO_URL}/issues/{stage['issue']}\"")
    return "\n".join([*lines, "```"])


def render(program: dict, hypothesis: dict) -> str:
    lines = ["# Executable OPH completion program", "",
             f"Generated by `tools/build_completion_program.py` from [`{PROGRAM}`](../{PROGRAM}).",
             "", f"Coordination: {issue_link(program['umbrella_issue'])}. "
             "Physical target: [the executable hypothesis](OPH_EXECUTABLE_HYPOTHESIS.md).", "",
             "**One native repair/memory process must supply the observer clocks, "
             "defect/maintenance-energy mechanism, interacting field/particle description "
             "and realized Einstein premises together.** Component benchmarks retain their "
             "own contracts; they cannot be pooled into that realization.", "",
             f"Program: `{program['program_id']}`. Boundary: `{program['boundary_status']}`. "
             f"Decision: `{program['decision']}`. No integrated gate is accepted by opening an issue.", "",
             "## Existing foundations to extend", "",
             "Physical construction gates remain open, but the program starts from retained "
             "positive constructions. The rows below link to existing canonical claims; "
             "they are reuse/navigation bindings, not new proof-discharge assertions. "
             "Preserve these capabilities while constructing the missing joins.", "",
             "| Foundation | Achieved component | Remaining join | Builds into |",
             "| --- | --- | --- | --- |"]
    for foundation in program["foundations"]:
        lines.append(f"| {foundation['title']} | {foundation['achieved']} | {foundation['boundary']} | {', '.join(foundation['supports'])} |")
    lines += ["", "In particular, retain the modular-clock candidate and the "
              "conditional 3+1 field construction based on the native response "
              "structure. The slow band's rank three and the native screen's "
              "spectral dimension two are distinct facts; the supplied radial "
              "extension remains part of that field construction's contract.", ""]
    if program.get("component_campaigns"):
        lines += ["## Current bounded investigations", "",
                  "These source-bound component campaigns record actual progress without "
                  "accepting dependent gates or freezing the larger program. Their "
                  "exact exclusions do not cover untested histories or refinements.", ""]
        for campaign in program["component_campaigns"]:
            contract = campaign["contract"]
            lines += [f"- [{campaign['id']}](../{contract['path']}): {contract['scope']} "
                      f"Owning gates: {', '.join(campaign['supports'])}."]
        lines += [""]
    lines += render_next_wave(program)
    class_accepted = next(stage for stage in program['stages']
                          if stage['id'] == 'CP-00')['status'] == 'ACCEPTED'
    boundary_summary = (
        "CP-00 has accepted the reviewed executable law/resource/interpretation class. "
        "The CP-01 observer and CP-13 theorem map retain their declared finite "
        "scope and separate acceptance criteria. Dependent physical gates require their own constructions and predecessor "
        "mappings. CP-12 audits every stage within that class; a whole-program exclusion "
        "must cover every admitted realization of its mandatory target conjunction."
        if class_accepted else
        "CP-00 must identify the complete law/resource/interpretation class "
        "whose exclusion would stop this program. The frozen finite campaign "
        "does not supply that larger boundary. Component investigations "
        "can continue provisionally; their results cannot discharge a dependent gate "
        "without the predecessor mappings. CP-12 audits every stage in parallel "
        "within each declared scope; a whole-program decision requires its frozen class.")
    lines += ["## Next work and dependencies", "",
             boundary_summary, "",
             mermaid(program), "",
             "| Gate | Research task | Work status / acceptance prerequisites | Issue |", "| --- | --- | --- | --- |"]
    for stage in program["stages"]:
        lines.append(f"| {stage['id']} | {stage['title']} | {readiness(stage, program)} | {issue_link(stage['issue'])} |")
    lines += ["", "The arrows are research acceptance dependencies, not a new theorem "
              "entailment graph. The existing [claim graph](../claims/dependency_graph.json) "
              "continues to record theorem/claim dependencies.", "", "## Decisions and stopping", "",
              "| Verdict | Evidence needed | Consequence |", "| --- | --- | --- |",
              "| Candidate rejection | A replayable failure of the specified law/state/access/readout contract | Stop that candidate and retain its result |",
              "| Class obstruction | A proof covering every realization in an explicitly quantified subclass | Stop that class; list remaining admitted classes |",
              "| Whole-program falsification | Verified contradiction in indispensable premises, exclusion of every frozen realization of a required join, or an unavoidable independently excluded prediction | Stop the frozen program; any premise/class change starts a new version |",
              "| Inconclusive | Open maps, incomplete coverage or an unsuccessful finite search | No confirmation and no impossibility claim |",
              "| Resource/value stop | A declared budget or return threshold is reached | A research stop, not scientific falsification |", "",
              ("The executable class is **frozen** by the reviewed CP-00 contract. "
               "No whole-program exclusion is claimed; it would need a proof covering "
               "every admitted realization of the mandatory target conjunction. "
               "A changed class requires a successor version with prior evidence preserved."
               if program['boundary_status']=='FROZEN' else
               "The boundary is **not frozen**, so no whole-program exclusion is claimed. "
               "Freezing CP-00 must specify the admitted class and its decision contract. "
               "Changing it never erases an earlier verdict."), "",
              "Every obstruction records target revision, mandatory premises, quantified "
              "class, evidence and independent check, coverage, surviving escape routes "
              "and the exact consequence. Issue closure alone is not scientific acceptance. "
              "The registry validator checks structure and artifact custody, not the truth "
              "of a mathematical proof or a physical interpretation.", "",
              "## Gate contracts", ""]
    for stage in program["stages"]:
        lines += [f"### {stage['id']}: {stage['title']}", "",
                  f"Owner: {issue_link(stage['issue'])}. Status: `{stage['status']}`. "
                  f"Requires: {', '.join(stage['requires']) or 'none'}. "
                  f"Hypothesis obligations: {', '.join(stage['obligations'])}.", "",
                  stage["deliverable"], "", "Acceptance:", ""]
        lines += [f"- {row}" for row in stage["acceptance"]]
        lines += ["", "Obstruction/stopping test: " + stage["obstruction"], "", "Retained starting evidence:", ""]
        lines += [f"- [{Path(path).name}](../{path})" for path in stage["evidence"]]
        if stage.get("accepted_evidence"):
            lines += ["", "Reviewed acceptance evidence:", ""]
            lines += [f"- [{Path(row['path']).name}](../{row['path']}): {row['scope']}"
                      for row in stage['accepted_evidence']]
        if stage.get("progress_evidence"):
            lines += ["", "Current component progress (not gate acceptance):", ""]
            lines += [f"- [{Path(row['path']).name}](../{row['path']}): {row['scope']}"
                      for row in stage["progress_evidence"]]
        lines += ["", "Earlier related work: " + ", ".join(issue_link(number) for number in stage["legacy_issues"]) + ".", ""]
    lines += ["## Working organization", "",
              program["legacy_issue_policy"] + " Live observer code, "
              "tests and retained results stay in `model/`; the existing particle code "
              "and theorem library keep their scientific ownership. Do not move or "
              "rewrite frozen packets to make the navigation look simpler.", "",
              "New work identifies its CP gate, fixed contract, nearest existing result, "
              "actual model/theorem map, falsifier and preservation checks. Prefer the "
              "smallest exact counterexample or native witness before larger runs. "
              "A supplied benchmark is evidence only for its stated conditional task. "
              "Each accepted gate needs a retained, source-bound result and the explicit "
              "join to its accepted predecessors.", "",
              "The generated [postdiction ledger](POSTDICTION_LEDGER.md) carries the "
              "hypothesis and program custody without adding a physical comparison. "
              "The [falsification program](OPH_FALSIFICATION_PROGRAM.md) retains its "
              "claim/branch-specific rejection rules and the stronger scope required "
              "for whole-program decisions.", ""]
    return "\n".join(lines)


def render_issue(stage: dict, program: dict) -> str:
    by_id = {row["id"]: row for row in program["stages"]}
    lines = ["## Goal", "", stage["deliverable"], "",
             f"Completion owner: {issue_link(program['umbrella_issue'])}. "
             f"Gate `{stage['id']}` in program `{program['program_id']}`.", "",
             f"[Common hypothesis]({REPO_URL}/blob/main/docs/OPH_EXECUTABLE_HYPOTHESIS.md) · "
             f"[Dependency DAG and gate contracts]({REPO_URL}/blob/main/{OUTPUT})", "",
             "## Dependencies", ""]
    lines += ([f"- [ ] {issue_link(by_id[dep]['issue'])}: {by_id[dep]['title']}" for dep in stage["requires"]]
              or ["No predecessor. This issue freezes the program boundary before integrated acceptance."])
    lines += ["", "Existing component work can proceed provisionally; integrated gate "
              "acceptance requires the actual maps to accepted predecessors. The same "
              "law, state/reference, observer access, clocks, energy and refinement "
              "domain must compose. A similarly named benchmark is not that map.", ""]
    for owner in program.get("next_wave", {}).get("owners", []):
        if owner["gate"] == stage["id"]:
            lines += ["## Next planned work", "", owner["deliverable"], "",
                      owner["may_start_now"], "", owner["join_condition"], "",
                      "This assignment does not change the frozen component campaign or "
                      "accept this gate. Wider primitives require the shared versioned contract.", ""]
    lines += ["## Acceptance", ""]
    lines += [f"- [ ] {row}" for row in stage["acceptance"]]
    lines += ["", "## Obstruction and stopping", "", stage["obstruction"], "",
              "Route every exclusion to CP-12 with target version, mandatory premises, "
              "quantified class, evidence/independent checks, coverage, surviving "
              "escape routes and consequence. Candidate failure, class no-go, theorem "
              "error, empirical rejection, inconclusive search and resource stop remain "
              "distinct. A whole-program verdict requires coverage of every frozen "
              "admissible realization; a changed class is a new program version.", "",
              "## Retained starting evidence", ""]
    lines += [f"- [{path}]({REPO_URL}/blob/main/{path})" for path in stage["evidence"]]
    lines += ["", "Related earlier work: " + ", ".join(issue_link(number) for number in stage["legacy_issues"]) + ".", "",
              "Opening or closing this issue does not establish a physical result. "
              "Acceptance needs the source-bound evidence, independent checks and "
              "predecessor joins recorded in the program registry. No retained failure "
              "or empirical exposure is reset by this reorganization.", ""]
    return "\n".join(lines)


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
