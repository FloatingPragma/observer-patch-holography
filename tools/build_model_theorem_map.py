"""Validate and render the model/theorem navigation crosswalk (stdlib only).

The crosswalk is the ``model_supply`` section of
``tracking/premise_register.json``; this tool renders it to
``model/PREMISES.md``. It checks references, inventory coverage and
generated-file freshness. Its contract section additionally enumerates every
registered claim and tracked Lean source, retains canonical inputs, and checks
typed, noncircular connections between narrowly scoped model clauses. It does
not infer claim-consumer edges, logical entailment, or a total of unproved
premises. The premise rows and the claim registry remain authoritative.

A source path may end in ``#fragment`` to link a section of its file; the
fragment must name an explicit ``<a id="fragment"></a>`` anchor in that file.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import date
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from urllib.parse import quote

from check_lean_theorem_count import DECL_RE, strip_comments

ROOT = Path(__file__).resolve().parents[1]
REGISTER = Path("tracking/premise_register.json")
SECTION = "model_supply"
CLAIMS = Path("claims/claim_registry.yaml")  # This repository uses JSON syntax.
OUTPUT = Path("model/PREMISES.md")
SCHEMA = "oph.model_theorem_map.v1"
SCOPES = {"active", "branch", "downstream", "metrology"}
STATUSES = {"conditional", "open", "imported"}
MODEL_IDS = {f"M{i:02d}" for i in range(1, 13)}
ROLES = ("consumed", "open", "boundary")
TOP_KEYS = {"schema", "as_of", "coverage_note", "groups", "contract_coverage"}
CONTRACT_SCHEMA = "oph.model_contract_coverage.v1"
DOMAIN_FACETS = ("carrier", "law", "preparation", "records", "controls", "clock", "regime")
GROUP_KEYS = {
    "id", "title", "links", "scope", "status", "premises", "constructed",
    "remaining", "closure", "sources", "assumptions", "claims",
}
GROUP_REQUIRED = GROUP_KEYS - {"assumptions", "claims"}
FRAGMENT = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*")


def lean_inventory(root: Path) -> dict[str, int]:
    """Enumerate public source declarations, not proof validity or finite scope.

    Use the same comment stripping and declaration grammar as the existing
    theorem-count gate. Git coverage includes new workspace sources before
    staging, but excludes ignored dependency/build trees; fixture roots
    without a Git repository use their own Lean directory.
    """
    root = Path(root).resolve()
    try:
        result = subprocess.run(["git", "-C", str(root), "ls-files", "-z",
                                 "--cached", "--others", "--exclude-standard", "--", "Lean"],
                                check=True, capture_output=True)
        paths = {root / item.decode() for item in result.stdout.split(b"\0")
                 if item.endswith(b".lean")}
    except (OSError, subprocess.CalledProcessError):
        paths = list((root / "Lean").rglob("*.lean"))
    return {path.relative_to(root).as_posix(): len(DECL_RE.findall(strip_comments(
                path.read_text(encoding="utf-8")))) for path in
            sorted(paths, key=lambda path: path.relative_to(root).parts)
            if path.is_file() and not any(part.startswith(".") or part == "lake-packages"
                                        for part in path.relative_to(root).parts)}


def _contract_sources(value: object, where: str, root: Path, errors: list[str]) -> None:
    if not isinstance(value, list) or not value:
        errors.append(f"{where}: expected nonempty exact source anchors")
        return
    for i, source in enumerate(value):
        spot = f"{where}[{i}]"
        if not isinstance(source, dict) or set(source) != {"path", "anchor"}:
            errors.append(f"{spot}: expected path and anchor")
            continue
        name, anchor = source["path"], source["anchor"]
        if not _text(name) or not _text(anchor):
            errors.append(f"{spot}: nonempty path and anchor required")
            continue
        path = PurePosixPath(name)
        target = (root / name).resolve()
        if (path.is_absolute() or ".." in path.parts or "\\" in name or
                not target.is_relative_to(root) or not target.is_file()):
            errors.append(f"{spot}: source must be an existing repository file")
        elif anchor not in target.read_text(encoding="utf-8"):
            errors.append(f"{spot}: exact source anchor absent")


def validate_contracts(data: object, claims: list[dict], root: Path) -> list[str]:
    """Check exhaustive custody and typed joins; never certify entailment.

    A binding records a reviewed, restricted clause, not discharge of its
    whole claim. All canonical assumptions/imports survive the projection.
    A join must name its consumed fact and every changed domain coordinate.
    Its source argument remains subject to ordinary mathematical review.
    """
    errors: list[str] = []
    keys = {"schema", "claim_ids", "lean_files", "domains", "bindings", "connections"}
    if not isinstance(data, dict):
        return ["contract_coverage: expected an object"]
    _keys(data, keys, keys, "contract_coverage", errors)
    if data.get("schema") != CONTRACT_SCHEMA:
        errors.append(f"contract_coverage.schema: expected {CONTRACT_SCHEMA}")
    ids = _string_list(data.get("claim_ids"), "contract_coverage.claim_ids", errors)
    registry = {row["claim_id"]: row for row in claims}
    if ids != list(registry):
        errors.append("contract_coverage: exact claim inventory/order changed; review every new or removed contract")
    files = _string_list(data.get("lean_files"), "contract_coverage.lean_files", errors)
    if files != list(lean_inventory(root)):
        errors.append("contract_coverage: exact Lean source inventory changed; review its model mapping")
    domains = data.get("domains")
    if not isinstance(domains, dict) or not domains:
        errors.append("contract_coverage.domains: expected nonempty named domains")
        domains = {}
    for name, domain in domains.items():
        if not _text(name) or not isinstance(domain, dict):
            errors.append("contract_coverage.domains: nonempty names and objects required")
            continue
        if set(domain) != set(DOMAIN_FACETS) or any(not _text(v) for v in domain.values()):
            errors.append(f"domain {name}: all seven nonempty domain coordinates required")
    bindings = data.get("bindings")
    if not isinstance(bindings, list):
        errors.append("contract_coverage.bindings: expected a list")
        bindings = []
    by_id: dict[str, dict] = {}
    binding_keys = {"id", "domain", "claims", "establishes", "requires", "assumptions",
                    "additional_assumptions", "imports", "scope", "limitations", "sources"}
    for i, binding in enumerate(bindings):
        where = f"binding[{i}]"
        if not isinstance(binding, dict):
            errors.append(f"{where}: expected an object")
            continue
        _keys(binding, binding_keys, binding_keys, where, errors)
        name = binding.get("id")
        if not _text(name) or name in by_id:
            errors.append(f"{where}: unique nonempty binding id required")
        else:
            by_id[name] = binding
        if not _text(binding.get("domain")) or binding["domain"] not in domains:
            errors.append(f"{where}: unknown model domain")
        named = _string_list(binding.get("claims"), f"{where}.claims", errors)
        if not named or set(named) - registry.keys():
            errors.append(f"{where}: nonempty known canonical claims required")
        for field in ("establishes", "requires", "assumptions", "imports"):
            _string_list(binding.get(field), f"{where}.{field}", errors)
        if not binding.get("establishes"):
            errors.append(f"{where}: a restricted established clause is required")
        for field, original in (("assumptions", "assumptions"), ("imports", "imported_results")):
            expected = sorted({item for cid in named if cid in registry
                               for item in registry[cid][original]})
            if binding.get(field) != expected:
                errors.append(f"{where}.{field}: every canonical input must remain, without discharge")
        extra_fields = sorted({field for cid in named if cid in registry for field in registry[cid]
                               if field.endswith("_assumptions")})
        extra = {field: sorted({value for cid in named if cid in registry
                                for value in registry[cid].get(field, [])}) for field in extra_fields}
        if binding.get("additional_assumptions") != extra:
            errors.append(f"{where}.additional_assumptions: all separate canonical assumption fields must remain")
        for field in ("scope", "limitations"):
            if not _text(binding.get(field)):
                errors.append(f"{where}.{field}: nonempty restricted scope required")
        _contract_sources(binding.get("sources"), f"{where}.sources", root, errors)
    connections = data.get("connections")
    if not isinstance(connections, list):
        errors.append("contract_coverage.connections: expected a list")
        connections = []
    delivered: set[tuple[str, str]] = set()
    adjacency = {name: set() for name in by_id}
    for i, connection in enumerate(connections):
        where = f"connection[{i}]"
        if not isinstance(connection, dict):
            errors.append(f"{where}: expected an object")
            continue
        required = {"from", "to", "fact", "coordinate_maps", "sources"}
        _keys(connection, required, required, where, errors)
        first, last = connection.get("from"), connection.get("to")
        if not _text(first) or not _text(last) or first not in by_id or last not in by_id:
            errors.append(f"{where}: unknown binding endpoint")
            continue
        fact = connection.get("fact")
        produced = by_id[first].get("establishes")
        required_facts = by_id[last].get("requires")
        if (not isinstance(produced, list) or not isinstance(required_facts, list) or
                fact not in produced or fact not in required_facts):
            errors.append(f"{where}: produced and required fact must agree exactly")
        if not _text(fact):
            continue
        if (last, fact) in delivered:
            errors.append(f"{where}: duplicate input delivery")
        delivered.add((last, fact)); adjacency[first].add(last)
        first_domain, last_domain = by_id[first].get("domain"), by_id[last].get("domain")
        a = domains.get(first_domain) if _text(first_domain) else None
        b = domains.get(last_domain) if _text(last_domain) else None
        if isinstance(a, dict) and isinstance(b, dict):
            changed = {field for field in DOMAIN_FACETS if a.get(field) != b.get(field)}
            maps = connection.get("coordinate_maps")
            if (not isinstance(maps, dict) or set(maps) != changed or
                    any(not _text(value) for value in maps.values())):
                errors.append(f"{where}: explicit maps required for exactly changed coordinates {sorted(changed)}")
        _contract_sources(connection.get("sources"), f"{where}.sources", root, errors)
    for name, binding in by_id.items():
        facts = binding.get("requires")
        for fact in facts if isinstance(facts, list) and all(_text(x) for x in facts) else []:
            if (name, fact) not in delivered:
                errors.append(f"binding {name}: missing connection for {fact}")
    pending = {name: set(targets) for name, targets in adjacency.items()}
    while pending:
        sinks = {name for name, targets in pending.items() if not targets}
        if not sinks:
            errors.append("contract_coverage: circular connection cannot establish its own input")
            break
        pending = {name: targets - sinks for name, targets in pending.items() if name not in sinks}
    return errors


def _text(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _string_list(value: object, where: str, errors: list[str]) -> list[str]:
    if not isinstance(value, list) or any(not _text(item) for item in value):
        errors.append(f"{where}: expected a list of nonempty strings")
        return []
    if len(set(value)) != len(value):
        errors.append(f"{where}: duplicate entries")
    return value


def _keys(value: dict, allowed: set[str], required: set[str],
          where: str, errors: list[str]) -> None:
    missing, extra = required - value.keys(), value.keys() - allowed
    if missing:
        errors.append(f"{where}: missing keys {sorted(missing)}")
    if extra:
        errors.append(f"{where}: unknown keys {sorted(extra)}")


def _registry_sets(premise_register: object, claim_registry: object,
                   errors: list[str]) -> tuple[set[str], set[str], set[str]]:
    """Check the inventory shapes needed here, without revalidating all science."""
    rows = premise_register.get("rows") if isinstance(premise_register, dict) else None
    claims = claim_registry.get("claims") if isinstance(claim_registry, dict) else None
    if not isinstance(rows, list) or not rows:
        errors.append("premise register: expected nonempty rows")
        rows = []
    if not isinstance(claims, list) or not claims:
        errors.append("claim registry: expected nonempty claims")
        claims = []
    premise_ids: set[str] = set()
    for i, row in enumerate(rows):
        if not isinstance(row, dict) or not _text(row.get("id")):
            errors.append(f"premise register row {i}: missing id")
            continue
        if row["id"] in premise_ids:
            errors.append(f"premise register: duplicate id {row['id']}")
        premise_ids.add(row["id"])
    claim_ids: set[str] = set()
    assumption_names: set[str] = set()
    for i, claim in enumerate(claims):
        if not isinstance(claim, dict) or not _text(claim.get("claim_id")):
            errors.append(f"claim registry row {i}: missing claim_id")
            continue
        cid = claim["claim_id"]
        if cid in claim_ids:
            errors.append(f"claim registry: duplicate id {cid}")
        claim_ids.add(cid)
        for field, value in claim.items():
            if field == "assumptions" or field.endswith("_assumptions"):
                assumption_names.update(_string_list(value, f"{cid}.{field}", errors))
        _string_list(claim.get("assumptions"), f"{cid}.assumptions", errors)
        _string_list(claim.get("imported_results"), f"{cid}.imported_results", errors)
        deps = claim.get("premise_dependencies")
        if not isinstance(deps, dict):
            errors.append(f"{cid}.premise_dependencies: expected an object")
            continue
        for role in ROLES:
            ids = _string_list(deps.get(role), f"{cid}.{role}", errors)
            unknown = set(ids) - premise_ids
            if unknown:
                errors.append(f"{cid}.{role}: unknown premises {sorted(unknown)}")
    return premise_ids, claim_ids, assumption_names


def validate(data: object, premise_register: object, claim_registry: object,
             root: Path) -> list[str]:
    """Return validation errors; accepts in-memory inputs and a fixture root."""
    errors: list[str] = []
    premise_ids, claim_ids, assumption_names = _registry_sets(
        premise_register, claim_registry, errors)
    if not isinstance(data, dict):
        return errors + ["map: expected an object"]
    _keys(data, TOP_KEYS, TOP_KEYS, "map", errors)
    if data.get("schema") != SCHEMA:
        errors.append(f"map.schema: expected {SCHEMA}")
    try:
        stamp = data.get("as_of")
        if not isinstance(stamp, str) or date.fromisoformat(stamp).isoformat() != stamp:
            raise ValueError
    except ValueError:
        errors.append("map.as_of: expected YYYY-MM-DD")
    if not _text(data.get("coverage_note")):
        errors.append("map.coverage_note: expected nonempty string")
    groups = data.get("groups")
    if not isinstance(groups, list) or not groups:
        return errors + ["map.groups: expected a nonempty list"]
    seen: set[str] = set()
    used_pr: set[str] = set()
    used_m: set[str] = set()
    root = Path(root).resolve()
    for i, group in enumerate(groups):
        where = f"group[{i}]"
        if not isinstance(group, dict):
            errors.append(f"{where}: expected an object")
            continue
        _keys(group, GROUP_KEYS, GROUP_REQUIRED, where, errors)
        gid = group.get("id")
        if not isinstance(gid, str) or not re.fullmatch(r"RE-[0-9]{2,}", gid):
            errors.append(f"{where}.id: expected RE- followed by at least two digits")
        elif gid in seen:
            errors.append(f"{where}: duplicate group id {gid}")
        else:
            seen.add(gid)
        for key in ("title", "constructed", "remaining", "closure"):
            if not _text(group.get(key)):
                errors.append(f"{where}.{key}: expected nonempty string")
        if not isinstance(group.get("scope"), str) or group["scope"] not in SCOPES:
            errors.append(f"{where}.scope: expected one of {sorted(SCOPES)}")
        if not isinstance(group.get("status"), str) or group["status"] not in STATUSES:
            errors.append(f"{where}.status: expected one of {sorted(STATUSES)}; closed is not supported")
        for field, known, used in (("links", MODEL_IDS, used_m),
                                  ("premises", premise_ids, used_pr)):
            values = _string_list(group.get(field), f"{where}.{field}", errors)
            if field == "links" and not values:
                errors.append(f"{where}.links: must identify at least one model criterion")
            unknown = set(values) - known
            if unknown:
                errors.append(f"{where}.{field}: unknown references {sorted(unknown)}")
            used.update(values)
        for field, known in (("claims", claim_ids), ("assumptions", assumption_names)):
            if field in group:
                values = _string_list(group[field], f"{where}.{field}", errors)
                if set(values) - known:
                    errors.append(f"{where}.{field}: unknown references {sorted(set(values) - known)}")
        sources = group.get("sources")
        if not isinstance(sources, list) or not sources:
            errors.append(f"{where}.sources: expected a nonempty list")
            continue
        source_keys: set[tuple[str, str]] = set()
        for j, source in enumerate(sources):
            spot = f"{where}.sources[{j}]"
            if not isinstance(source, dict):
                errors.append(f"{spot}: expected an object")
                continue
            _keys(source, {"path", "label", "anchor"}, {"path", "label"}, spot, errors)
            if not _text(source.get("label")):
                errors.append(f"{spot}.label: expected nonempty string")
            anchor = source.get("anchor")
            if "anchor" in source and not _text(anchor):
                errors.append(f"{spot}.anchor: expected nonempty exact substring")
            path = source.get("path")
            if not _text(path):
                errors.append(f"{spot}.path: expected a repository-relative path")
                continue
            file_path, has_fragment, fragment = path.partition("#")
            relative = PurePosixPath(file_path)
            if relative.is_absolute() or ".." in relative.parts or "\\" in path:
                errors.append(f"{spot}.path: must stay inside the repository")
                continue
            target = (root / file_path).resolve()
            if not target.is_relative_to(root) or not target.is_file():
                errors.append(f"{spot}.path: missing file or outside repository: {path}")
                continue
            if has_fragment:
                if not FRAGMENT.fullmatch(fragment):
                    errors.append(f"{spot}.path: fragment must be lowercase words joined by hyphens")
                    continue
                try:
                    if f'<a id="{fragment}"></a>' not in target.read_text(encoding="utf-8"):
                        errors.append(f"{spot}.path: anchor #{fragment} absent from {file_path}")
                except (OSError, UnicodeError) as exc:
                    errors.append(f"{spot}: cannot read source: {exc}")
            key = (path, anchor if isinstance(anchor, str) else "")
            if key in source_keys:
                errors.append(f"{spot}: duplicate source path/anchor")
            source_keys.add(key)
            if _text(anchor):
                try:
                    if anchor not in target.read_text(encoding="utf-8"):
                        errors.append(f"{spot}.anchor: exact substring absent from {path}")
                except (OSError, UnicodeError) as exc:
                    errors.append(f"{spot}: cannot read source: {exc}")
    if premise_ids - used_pr:
        errors.append(f"map: omitted canonical premises {sorted(premise_ids - used_pr)}")
    if MODEL_IDS - used_m:
        errors.append(f"map: missing model criteria {sorted(MODEL_IDS - used_m)}")
    if isinstance(claim_registry, dict) and isinstance(claim_registry.get("claims"), list):
        errors.extend(validate_contracts(data.get("contract_coverage"),
                                         claim_registry["claims"], root))
    return errors


def inventory(data: dict, premise_register: dict, claim_registry: dict,
              root: Path | None = None) -> dict:
    """Count named registry metadata, not unproved mathematical obligations."""
    claims, groups = claim_registry["claims"], data["groups"]
    ordinary = {s for c in claims for s in c["assumptions"]}
    fields = sorted({k for c in claims for k in c if k.endswith("_assumptions")})
    extras = {}
    for field in fields:
        names = {s for c in claims for s in c.get(field, [])}
        extras[field] = {
            "claims_with_field": sum(field in c for c in claims),
            "occurrences": sum(len(c.get(field, [])) for c in claims),
            "unique_names": len(names), "also_ordinary": len(names & ordinary),
        }
    return {
        "claims": len(claims), "canonical_premises": len(premise_register["rows"]),
        "covered_premises": len({p for g in groups for p in g["premises"]}),
        "covered_model_criteria": len({m for g in groups for m in g["links"]}),
        "ordinary_assumptions": len(ordinary),
        "imported_names": len({s for c in claims for s in c["imported_results"]}),
        "extra_assumption_fields": extras,
        "role_counts": {role: len({p for c in claims for p in
                                   c["premise_dependencies"][role]}) for role in ROLES},
        "classifications": dict(sorted(Counter(
            c["premise_dependencies"].get("classification", "unspecified")
            for c in claims).items())),
        "groups": len(groups),
        "mapped_assumptions": len({a for g in groups for a in g.get("assumptions", [])}),
        "mapped_claims": len({c for g in groups for c in g.get("claims", [])}),
        "scope_counts": dict(sorted(Counter(g["scope"] for g in groups).items())),
        "contract_rows": claims,
        "premise_rows": premise_register["rows"],
        "lean_inventory": lean_inventory(ROOT if root is None else root),
    }


DISCHARGE_ORDER = (
    ("derived_from_native_records", "derived from native records"),
    ("constructed_on_declared_inputs", "constructed on declared inputs"),
    ("obstructed", "obstructed"),
    ("declared", "declared input"),
    ("not_addressed", "not addressed"),
)


def render_discharge(counts: dict) -> list[str]:
    """Render what the current executable implementation does with each premise."""
    # The premise-register builder requires a discharge record on every live row;
    # a fixture row without one renders as not addressed.
    default = {"status": "not_addressed", "evidence": [], "note": ""}
    rows = [dict(row, discharge=row.get("discharge", default)) for row in counts["premise_rows"]]
    labels = dict(DISCHARGE_ORDER)
    lines = [
        "## What the current implementation does with each premise", "",
        "The executable construction in this directory is the current implementation of the",
        "specification's dynamical details; the twelve-port carrier is fixed, the details are not.",
        "For each canonical premise the register records whether this implementation derives it",
        "from native records, constructs an instance on declared inputs, proves an obstruction,",
        "takes it as a declared input, or does not address it. The status certifies no physical",
        "premise, and a theorem-library result holds on its premises whatever the status.", "",
        "| Status | Premises |", "| --- | ---: |",
    ]
    for status, label in DISCHARGE_ORDER:
        lines.append(f"| {label} | {sum(1 for r in rows if r['discharge']['status'] == status)} |")
    lines += ["", "| Premise | Current implementation | Basis | Evidence |", "| --- | --- | --- | --- |"]
    for row in rows:
        if row["discharge"]["status"] == "not_addressed":
            continue
        evidence = ", ".join(
            f"[{Path(path).name}](../{quote(path, safe='/#')})" for path in row["discharge"]["evidence"]
        )
        lines.append(
            f"| {row['id']} {_cell(row['name'])} | {labels[row['discharge']['status']]} |"
            f" {_cell(row['discharge']['note'])} | {evidence} |"
        )
    not_addressed = [row["id"] for row in rows if row["discharge"]["status"] == "not_addressed"]
    lines += ["", "Not addressed by the current implementation: " + ", ".join(not_addressed) + ".", ""]
    return lines


def render_contracts(data: dict, counts: dict) -> list[str]:
    """Render a compact view; retain exhaustive inventories in the checked JSON."""
    coverage = data["contract_coverage"]
    bindings = coverage["bindings"]
    rows = counts["contract_rows"]
    bound_claims = {cid for binding in bindings for cid in binding["claims"]}
    lean_evidence = {path for row in rows for path in row.get("evidence", [])
                     if path in counts["lean_inventory"]}
    lines = ["", '<a id="finite-contract-coverage"></a>', "",
             "## Finite-contract coverage and compatible model connections", "",
             "The complete machine inventory is the `model_supply.contract_coverage` section",
             "of the [premise register](../tracking/premise_register.json). It retains every",
             "canonical claim ID and tracked Lean module, all seven domain coordinates, every",
             "binding's canonical assumptions/imports and every sourced connection map. The",
             "[claim registry](../claims/claim_registry.yaml) owns statements, evidence and",
             "premise roles; the [Lean sources](../Lean/) retain their quantified hypotheses.", "",
             "Coverage includes non-theorem claims and is a superset of finite contracts.",
             "Lexical declaration counts decide neither finiteness nor proof validity; Lean",
             "builds remain responsible for kernel checks. Registry or source-module additions",
             "require explicit inventory updates. All unbound claims have the same status:",
             "**Open: no same-model instantiation recorded**. This does not invalidate a",
             "theorem on its own fixture or add a physical premise to an explicit non-consumer.", "",
             "Each binding is a restricted clause, not discharge of its whole registry row.",
             "The JSON retains that row's complete assumption/import union, including inputs",
             "outside the clause. The checker requires exact produced/consumed facts, acyclic",
             "connections and a sourced map for every changed carrier, law, preparation,",
             "record access, control, clock or regime coordinate. These are coverage and",
             "contract checks: mathematical entailment requires review; strings and source",
             "anchors are not proofs. No entry silently discharges a canonical premise.", "",
             "| Coverage | Count |", "| --- | ---: |",
             f"| Canonical claims retained | {len(rows)} |",
             f"| Claims with restricted bindings | {len(bound_claims)} |",
             f"| Claims with open common-model mapping | {len(rows)-len(bound_claims)} |",
             f"| Tracked Lean source modules retained | {len(counts['lean_inventory'])} |",
             f"| Lean modules referenced by canonical evidence | {len(lean_evidence)} |",
             f"| Public theorem/lemma declarations | {sum(counts['lean_inventory'].values())} |",
             f"| Restricted bindings / sourced connections | {len(bindings)} / {len(coverage['connections'])} |", "",
             "### Restricted model bindings", "",
             "The domain key resolves to the full seven-coordinate contract in the machine",
             "inventory. The linked sources give each argument; the boundary column retains",
             "its principal restrictions. Full assumption names and source anchors are not",
             "duplicated here.", "",
             "| Binding and domain | Established clause | Boundary | Sources |",
             "| --- | --- | --- | --- |"]
    for binding in bindings:
        sources = "; ".join(f"[{Path(s['path']).name}](../{quote(s['path'], safe='/')})"
                            for s in binding["sources"])
        lines.append("| " + " | ".join(map(_cell, [
            f"{binding['id']} (`{binding['domain']}`)", binding["scope"],
            binding["limitations"], sources])) + " |")
    if coverage["connections"]:
        lines += ["", "### Sourced domain connections", "",
                  "Every listed changed coordinate has an explicit argument and source anchor",
                  "in the machine inventory. Omitted coordinates must agree exactly.", "",
                  "| Producer → consumer | Required fact | Changed coordinates |",
                  "| --- | --- | --- |"]
    for connection in coverage["connections"]:
        lines.append("| " + " | ".join(map(_cell, [
            f"{connection['from']} → {connection['to']}", connection["fact"],
            ", ".join(connection["coordinate_maps"]) or "None"])) + " |")
    return lines


def _cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def _code(value: str) -> str:
    # Anchors may themselves contain Markdown backticks.
    runs = re.findall(r"`+", value)
    delimiter = "`" * (max(map(len, runs), default=0) + 1)
    return f"{delimiter} {value} {delimiter}"


def render(data: dict, counts: dict) -> str:
    """Deterministic Markdown; call after validate succeeds."""
    lines = [
        "# OPH reverse-engineering table: model and theorem premises", "",
        f"As of {data['as_of']}. Generated from [the navigation map](../tracking/premise_register.json)",
        "by `python tools/build_model_theorem_map.py`. Edit the input, then regenerate.", "",
        data["coverage_note"], "",
        "These groups are workgroup navigation crosswalks, not claim-consumer edges or changes",
        "to proof status. A group links related construction, branch, downstream or metrology",
        "work; that association does not make every listed premise an antecedent of every",
        "listed theorem. Existing claim-registry edges and `explicit_non_consumer` declarations",
        "remain authoritative. Physical attachments excluded by a theorem are not hidden premises",
        "of that theorem. Group statuses describe the remaining connection, not whether a theorem",
        "is proved; the constructed results are recorded separately below. This schema contains",
        "no proof-discharge edges; the implementation status of each premise is the separate",
        "[discharge table](#what-the-current-implementation-does-with-each-premise).", "",
        "## Inventory and limits", "",
        "Counts below come from the [canonical premise register](../tracking/premise_register.json)",
        "and [claim registry](../claims/claim_registry.yaml). Ordinary assumption names are",
        "explained in the [assumption dictionary](../claims/assumption_dictionary.md).", "",
        "| Inventory | Count |", "| --- | ---: |",
        f"| Canonical premise IDs represented in this crosswalk | {counts['covered_premises']}/{counts['canonical_premises']} |",
        f"| Model criteria represented | {counts['covered_model_criteria']}/12 |",
        f"| Registered claims | {counts['claims']} |",
        f"| Unique ordinary `assumptions` names | {counts['ordinary_assumptions']} |",
        f"| Unique `imported_results` names | {counts['imported_names']} |",
        f"| Workgroups | {counts['groups']} |",
        f"| Assumption names explicitly linked by optional group fields | {counts['mapped_assumptions']} |",
        f"| Claims explicitly linked by optional group fields | {counts['mapped_claims']} |", "",
        "Canonical-ID coverage is a reference check. Ordinary assumptions and imported names",
        "are inventoried only; they have not all been individually, semantically bound to the",
        "model by this page. Names can overlap or encode several obligations. These inventories",
        "cannot be added to obtain a certified total of unproved premises, and coverage is not",
        "a global completeness or proof claim.", "",
        "| Existing claim-registry role | Unique canonical IDs in that role |",
        "| --- | ---: |",
    ]
    lines.extend(f"| `{role}` | {counts['role_counts'][role]} |" for role in ROLES)
    lines += ["", "These role sets overlap. They record metadata on existing claims; `consumed`",
              "does not mean discharged, and `boundary` does not mean an unproved antecedent.", "",
              "Additional assumption fields are counted separately, not silently folded into",
              "the ordinary-assumption inventory:", "",
              "| Field | Claims containing field | Occurrences | Unique names | Also ordinary |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for field, values in counts["extra_assumption_fields"].items():
        lines.append(f"| `{field}` | {values['claims_with_field']} | {values['occurrences']} | {values['unique_names']} | {values['also_ordinary']} |")
    if not counts["extra_assumption_fields"]:
        lines.append("| None present | 0 | 0 | 0 | 0 |")
    lines += ["", "Existing claim classifications: " + "; ".join(
        f"`{key}`: {value}" for key, value in counts["classifications"].items()) + ".", "",
        *render_discharge(counts),
        "## Remaining-input and dependency rule", "",
        "For a selected scope G, let I(g) denote each group's required inputs on its",
        "declared family, state, observable algebra, clock, regime and quantitative error domain.",
        "The remaining-input set is U(G) = union over g in G of I(g), with genuinely identical",
        "inputs counted once. This is a bookkeeping rule, not a machine-computed logical set:",
        "the prose, ordinary assumptions and imports require semantic review. Active",
        "construction, optional branches, downstream physical attachments and metrology are",
        "separate scopes; unioning them does not turn all of them into theorem antecedents.", "",
        "A result H implies C discharges an input P only when C establishes P on that same",
        "domain, with the required errors, and H has independent supplied or established",
        "support. Substitution retains any outstanding inputs of H. Circular references,",
        "a change of regulator/state/code, or an imported calibration relabeled as a derivation",
        "do not close an obligation. Establishing a conditional theorem and realizing its",
        "hypotheses in the candidate are distinct tasks. This tool checks neither entailment",
        "nor acyclicity of scientific arguments and changes no existing dependency edges.", "",
        "## Workgroup crosswalk", "",
    ]
    if {f"RE-{i:02d}" for i in range(5, 13)} <= {g["id"] for g in data["groups"]}:
        lines += [
            "For the active path, RE-05–RE-12 form the connected action, inertia, quantum",
            "representation, state, history, refinement and clock cluster. Downstream rows",
            "are not additional prerequisites for the supplied finite scalar theorem.", "",
        ]
    lines += ["| Group | Model criteria | Scope | Work status | Canonical navigation IDs |",
              "| --- | --- | --- | --- | --- |"]
    for group in data["groups"]:
        lines.append(f"| [{_cell(group['id'] + ': ' + group['title'])}](#{group['id'].lower()}) | {', '.join(group['links'])} | {group['scope']} | {group['status']} | {', '.join(group['premises']) or 'None'} |")
    for group in data["groups"]:
        lines += ["", f"<a id=\"{group['id'].lower()}\"></a>", "",
                  f"### {group['id']}: {group['title']}", "",
                  f"**Constructed:** {group['constructed']}", "",
                  f"**Remaining input:** {group['remaining']}", "",
                  f"**Closure condition:** {group['closure']}", ""]
        if group.get("assumptions"):
            lines += ["Explicit assumption-name navigation: " + ", ".join(
                _code(a) for a in group["assumptions"]) + ".", ""]
        if group.get("claims"):
            lines += ["Existing claim-row navigation (no new edges): " + ", ".join(
                _code(c) for c in group["claims"]) + ".", ""]
        lines += ["Sources:", ""]
        for source in group["sources"]:
            label = source["label"].replace("[", "\\[").replace("]", "\\]")
            item = f"- [{label}](../{quote(source['path'], safe='/#')})"
            if "anchor" in source:
                item += "; exact source anchor: " + _code(source["anchor"])
            lines.append(item + ".")
    lines.extend(render_contracts(data, counts))
    return "\n".join(lines) + "\n"


def _unique_object(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="read-only validation and freshness check")
    args = parser.parse_args(argv)
    try:
        premises, claims = (read_json(ROOT / path) for path in (REGISTER, CLAIMS))
        if not isinstance(premises, dict) or SECTION not in premises:
            print(f"{REGISTER}: missing {SECTION} section", file=sys.stderr)
            return 1
        data = premises[SECTION]
        errors = validate(data, premises, claims, ROOT)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            return 1
        counts = inventory(data, premises, claims, ROOT)
        content = render(data, counts)
        target = ROOT / OUTPUT
        if args.check:
            if not target.is_file() or target.read_text(encoding="utf-8") != content:
                print(f"{OUTPUT} is stale; run python tools/build_model_theorem_map.py", file=sys.stderr)
                return 1
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content, encoding="utf-8")
        print(f"Model/theorem map {'checked' if args.check else 'generated'}: "
              f"{counts['groups']} groups; {counts['covered_premises']}/"
              f"{counts['canonical_premises']} canonical IDs; "
              f"{counts['covered_model_criteria']}/12 model criteria; "
              f"{len(counts['contract_rows'])} claim contracts; "
              f"{len(counts['lean_inventory'])} Lean modules. Typed contract checks, not proof discharge.")
        return 0
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"Cannot build model/theorem map: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
