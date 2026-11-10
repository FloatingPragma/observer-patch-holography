#!/usr/bin/env python3
"""Validate the OPH claim-registry seed files.

The registry file is JSON-compatible YAML, so this validator avoids an external
YAML dependency while keeping the requested `.yaml` public path.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from closed_lanes import SUPERSEDED_LANE_SUCCESSORS, is_known, live_owners  # noqa: E402

PAPER_EXTERNAL_REGISTRY_PATTERNS = [
    "claims/claim_registry",
    "\\ophid{claims/",
]

REQUIRED_CLAIM_FIELDS = {
    "claim_id",
    "statement",
    "owner_paper",
    "tier",
    "assumptions",
    "imported_results",
    "oph_specific_delta",
    "novelty_type",
    "evidence",
    "falsifier",
    "scope_if_false",
    "status",
    "claim_class",
    "gates",
    "premise_dependencies",
}

CLAIM_REGISTRY_SCHEMA = 3
PREMISE_REGISTER_RELATIVE = Path("tracking") / "premise_register.json"
PREMISE_DEPENDENCY_FIELDS = ("consumed", "open", "boundary")
PREMISE_DEPENDENCY_CLASSIFICATIONS = {
    "explicit_edges",
    "explicit_non_consumer",
}
PREMISE_ID = re.compile(r"^PR-[0-9]{2}$")
PREMISE_DEPENDENCY_PROJECTION_SHA256 = (
    "9fcda40ea8eb2c6ee7ea5b2ecee2bf46e20a903f93b82d581a34abc4358be419"
)

# Controlled claim classification (issue #512). `status` stays descriptive
# free text; `claim_class` is the machine-checked epistemic class.
CLAIM_CLASS_VOCABULARY = {
    # Definitions, contracts, programs, interpretive frames; no theorem or
    # executed artifact is asserted.
    "declared_structure",
    # A proved implication whose premises include declared (not source-derived)
    # inputs; no claim that the physical world instantiates the branch.
    "conditional_implication",
    # The assertion that a declared physical branch is entered or nonempty.
    "branch_entry",
    # An executed pipeline artifact (audit, archive, selector run, protocol
    # record); the claim is the artifact's existence and content.
    "emitted_artifact",
    # A quantitative pipeline landing compared against measurement, with a
    # declared endpoint boundary or error budget.
    "empirical_implementation",
    # Physically established end to end; requires every live gate closed.
    "physical_establishment",
}

# Classes whose public wording may assert physical establishment. Every other
# class fails the wording gate if its statement claims establishment.
PROMOTED_CLASSES = {"physical_establishment"}

ESTABLISHMENT_WORDING = re.compile(
    r"physically established|experimentally confirmed|physically confirmed"
    r"|empirically established|is established as physical",
    re.IGNORECASE,
)

PINNED_GITHUB_EVIDENCE = re.compile(
    r"https://github\.com/[^/]+/[^/]+/blob/[0-9a-f]{40}/.+"
)

# Artifact medium of one evidence path, keyed by suffix. An unlisted suffix
# fails closed: a new artifact kind must say whether it carries a
# machine-checked proof, an executed run, an emitted record, or prose.
EVIDENCE_MEDIA_BY_SUFFIX = {
    ".lean": "machine_checked_proof",
    ".py": "executed_program",
    ".cpp": "executed_program",
    ".json": "emitted_record",
    ".npz": "emitted_record",
    ".npy": "numeric_data",
    ".tex": "paper_prose",
    ".md": "markdown_prose",
}

# Media that carry a human-readable argument instead of a machine-checked or
# executed artifact.
PROSE_EVIDENCE_MEDIA = {"paper_prose", "markdown_prose"}

# Classes whose row asserts a theorem. `declared_structure` fixes definitions
# and `emitted_artifact` asserts an artifact's existence and content, so
# neither one owes a proof medium.
THEOREM_ASSERTING_CLASSES = {
    "conditional_implication",
    "branch_entry",
    "empirical_implementation",
    "physical_establishment",
}

PROOF_MEDIUM_FIELD = "proof_medium"
PROOF_MEDIUM_VOCABULARY = {"paper_prose"}

# Directory roots whose owner file is a paper source: the TeX paper folders and
# the research notes in docs/research. Every other owner is a protocol or run
# record and must declare that medium explicitly.
PAPER_OWNER_ROOTS = ("paper", "extra", "cosmology", "flagship", "docs/research")
OWNER_MEDIUM_FIELD = "owner_medium"
OWNER_MEDIUM_VOCABULARY = {"protocol_record"}

REQUIRED_TOPICAL_GATE_OWNERS_FIELD = "required_topical_gate_owners"

# Keys that name gates without being a gate list. Each one states why the key
# cannot hide a gate from the invariants that read `gates`.
GATE_SIDE_CHANNEL_ALLOWLIST = {
    REQUIRED_TOPICAL_GATE_OWNERS_FIELD: (
        "validated as a sorted subset of `gates` by "
        "check_topical_gate_owner_projection"
    ),
}

# Historical topical lanes for claim rows whose own statement leaves the
# corresponding physical attachment open. Every such lane is closed; the open
# issue that owns the remaining work is `closed_lanes.live_owner(lane)`.  This is deliberately an explicit
# claim-by-claim policy rather than a keyword rule: exact finite helpers and
# scoped no-go results remain valid without inheriting every downstream lane.
# Every claim named here also declares `required_topical_gate_owners`, and
# check_topical_gate_owner_projection compares the two surfaces in both
# directions, so a one-sided edit to either the registry row or this policy
# fails instead of silently dropping topical ownership.
REQUIRED_V3_TOPIC_GATES_BY_CLAIM: dict[str, frozenset[int]] = {
    "OPH-UNIFIED-TYPED-SPINE": frozenset({740}),
    "OPH-GR-D6-CAPACITY": frozenset({742}),
    "OPH-GR-D6-HORIZON-RECORD": frozenset({742}),
    "OPH-GR-DS-SHOCK-SIGN-ATTACHMENT": frozenset({742}),
    "OPH-GR-DS-DISCRETE-SHOCK-SPECTRUM": frozenset({742}),
    "OPH-COSMO-SCREEN-SPECTRUM": frozenset({742}),
    "OPH-WZ-STRICT-1L-POLE-MAP": frozenset({743, 745}),
    "OPH-SM-Q1-LOCAL-G6": frozenset({743, 745}),
    "OPH-SM-Q2E-CHIRAL-MEASURE-CRITERION": frozenset({743}),
    "OPH-SM-Q2H-POSITIVE-HAMILTONIAN-SOUNDNESS": frozenset({743}),
    "OPH-SM-Q3-BV-RESTORATION": frozenset({743}),
    "OPH-SM-Q4-OS-OBSERVABLE-SECTOR": frozenset({743}),
    "OPH-SM-Q4-RESONANCE-CONTINUATION": frozenset({743}),
    "OPH-A5-PRIMITIVE-PORT-SPIN6": frozenset({742}),
    "OPH-A5-SEAM-CURRENT-EDGE30": frozenset({742}),
    "OPH-MODAL-MAXWELL-FACTORIZATION-BOUNDARY": frozenset({754}),
    "OPH-FINITE-CONSERVATION-WARD-PRECURSOR": frozenset({743}),
    "OPH-QUARK-REGISTER-CLEBSCH": frozenset({745}),
    "OPH-KOIDE-CIRCULANT-IDENTITY": frozenset({745}),
    "OPH-W5-STABILISER-POTENTIAL-BOUNDARY": frozenset({745}),
    "OPH-SM-ROUTE-IDENTIFICATION": frozenset({740}),
    "OPH-Q-N-RESERVE-CANDIDATES": frozenset({742}),
    "OPH-HIER-EW": frozenset({740, 742, 745}),
    "OPH-SCREEN-24-CLOCK-DETERMINANT": frozenset({745}),
    "OPH-ALPHA-PIXEL": frozenset({744}),
    "OPH-DM-CONT": frozenset({742}),
    "OPH-YM-GAP": frozenset({743, 744}),
    "OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX": frozenset({743}),
    "OPH-FINITE-UNITARY-SCATTERING-LIMIT-NO-GO": frozenset({743}),
}


class DuplicateKeyError(ValueError):
    """Raised when a JSON-compatible YAML object repeats a key."""


def _unique_object(pairs: list[tuple[str, object]]) -> dict:
    result: dict = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateKeyError(f"duplicate object key {key!r}")
        result[key] = value
    return result


def load_json(path: Path) -> dict:
    try:
        return json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object
        )
    except (json.JSONDecodeError, DuplicateKeyError) as exc:
        raise SystemExit(f"{path}: invalid JSON-compatible YAML: {exc}") from exc


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, strict=True)
        try:
            headers = reader.fieldnames
            require(
                bool(headers) and all(header.strip() for header in headers),
                f"{path}: missing or empty CSV header",
            )
            duplicates = sorted(
                header for header, count in Counter(headers).items() if count > 1
            )
            require(not duplicates, f"{path}: duplicate CSV headers: {duplicates}")
            rows = []
            for row in reader:
                location = f"{path}: CSV line {reader.line_num}"
                require(None not in row, f"{location}: overflow CSV cells")
                missing = [header for header, value in row.items() if value is None]
                require(not missing, f"{location}: missing CSV cells: {missing}")
                empty = [header for header, value in row.items() if not value.strip()]
                require(not empty, f"{location}: empty CSV cells: {empty}")
                rows.append(row)
            return rows
        except csv.Error as exc:
            raise SystemExit(
                f"{path}: invalid CSV at line {reader.line_num}: {exc}"
            ) from exc


def check_exact_claim_id_projection(
    *,
    surface: str,
    projected_ids: list[str],
    canonical_ids: set[str],
    one_row_per_claim: bool,
) -> None:
    """Require exact ID coverage by a projection of the canonical registry.

    Checking only ``projected_ids <= canonical_ids`` is not bidirectional: a
    newly registered claim can silently disappear from a matrix or DAG while
    every remaining row still names a valid ID. Surfaces defined to have one
    row per claim also reject duplicates. The falsification matrix may retain
    several independently scoped falsification rows for one canonical claim.
    """
    duplicates = sorted(
        claim_id for claim_id, count in Counter(projected_ids).items() if count > 1
    )
    require(
        not one_row_per_claim or not duplicates,
        f"{surface}: duplicate claim ids: {duplicates}",
    )

    projected = set(projected_ids)
    missing = sorted(canonical_ids - projected)
    extra = sorted(projected - canonical_ids)
    require(
        not missing and not extra,
        f"{surface}: claim ids do not exactly match the canonical registry; "
        f"missing={missing}, extra={extra}",
    )


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def check_dependency_graph_edges(raw_edges: object, nodes: set[str]) -> None:
    """Validate that the scientific claim graph is a simple directed DAG.

    A claim dependency points from an antecedent to a downstream consumer.
    Self-loops, duplicate endpoint pairs, and directed cycles therefore make
    the ancestry ill-typed even when every endpoint is a registered claim.
    """

    require(
        isinstance(raw_edges, list),
        "claims/dependency_graph.json: edges must be a list",
    )
    adjacency: dict[str, list[str]] = {node: [] for node in nodes}
    seen_pairs: set[tuple[str, str]] = set()
    for edge in raw_edges:
        require(
            isinstance(edge, dict),
            f"dependency graph edge must be an object: {edge!r}",
        )
        source = edge.get("from")
        target = edge.get("to")
        require(
            source in nodes,
            f"dependency graph edge source is not a declared node: {edge}",
        )
        require(
            target in nodes,
            f"dependency graph edge target is not a declared node: {edge}",
        )
        require(edge.get("role"), f"dependency graph edge lacks role: {edge}")
        require(
            source != target,
            f"dependency graph self-loop is forbidden: {source}",
        )
        pair = (source, target)
        require(
            pair not in seen_pairs,
            f"dependency graph repeats directed edge: {source} -> {target}",
        )
        seen_pairs.add(pair)
        adjacency[source].append(target)

    state: dict[str, int] = {}
    stack: list[str] = []
    stack_index: dict[str, int] = {}

    def visit(node: str) -> None:
        state[node] = 1
        stack_index[node] = len(stack)
        stack.append(node)
        for target in adjacency[node]:
            if state.get(target, 0) == 0:
                visit(target)
            elif state[target] == 1:
                start = stack_index[target]
                cycle = stack[start:] + [target]
                raise SystemExit(
                    "dependency graph contains directed cycle: "
                    + " -> ".join(cycle)
                )
        stack.pop()
        stack_index.pop(node)
        state[node] = 2

    for node in sorted(nodes):
        if state.get(node, 0) == 0:
            visit(node)


def evidence_is_available(root: Path, evidence: str) -> bool:
    """Accept repository files or immutable GitHub source links.

    Cross-repository producers cannot exist in a clean checkout of this
    repository. Their evidence links must therefore pin a full commit hash;
    branch and tag URLs are mutable and remain invalid.
    """

    return (root / evidence).exists() or bool(PINNED_GITHUB_EVIDENCE.fullmatch(evidence))


def release_id_from_tex(root: Path) -> str:
    text = (root / "paper" / "release_info.tex").read_text(encoding="utf-8")
    match = re.search(r"\\newcommand\{\\OPHPaperReleaseID\}\{([^}]+)\}", text)
    require(match is not None, "paper/release_info.tex does not define OPHPaperReleaseID")
    return match.group(1)


def check_claim_views(root: Path) -> None:
    """The claim context covers every claim, and the views render from the registry."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import build_claim_views

    registry = load_json(root / "claims" / "claim_registry.yaml")
    ids = [claim["claim_id"] for claim in registry["claims"]]
    context = registry.get("claim_context")
    require(isinstance(context, dict), "claim registry: missing claim_context section")
    require(list(context) == ids, "claim_context must list every claim once, in registry order")
    for claim_id, entry in context.items():
        require(
            {"prior_work", "owner_section", "falsifiers", "depends_on"} <= set(entry)
            <= {"prior_work", "owner_section", "falsifiers", "depends_on", "matrix_summary"},
            f"{claim_id}: claim_context fields are prior_work, owner_section, falsifiers, depends_on and an optional matrix_summary",
        )
        prior = entry["prior_work"]
        require(
            isinstance(prior, dict) and set(prior) == {"summary", "overlap", "imported_step", "citations"},
            f"{claim_id}: prior_work must carry summary, overlap, imported_step and citations",
        )
        require(isinstance(prior["citations"], list), f"{claim_id}: prior_work.citations must be a list")
        falsifiers = entry["falsifiers"]
        require(
            isinstance(falsifiers, dict)
            and {"physical_identification", "phenomenological"} <= set(falsifiers)
            <= {"physical_identification", "phenomenological", "mathematical"},
            f"{claim_id}: falsifiers must carry physical_identification and phenomenological",
        )
        require(
            isinstance(entry["depends_on"], list)
            and all(isinstance(edge, dict) and set(edge) == {"claim_id", "role"} and edge["claim_id"] in context
                    for edge in entry["depends_on"]),
            f"{claim_id}: depends_on must list claim_id/role records of registered claims",
        )
    stale = build_claim_views.stale_views(root)
    require(
        not stale,
        "claim views differ from the registry: " + ", ".join(stale)
        + "; run python3 tools/build_claim_views.py",
    )

    # Every established result a claim uses is cited from the verified
    # bibliography or marked as a result of the project itself, and the
    # papers' credited-prior-work sections are current.
    import build_prior_work

    try:
        bibliography = build_prior_work.load_bibliography(root / "claims" / "prior_work.json")
        build_prior_work.check_citations(registry["claims"], context, bibliography)
        stale_prior = build_prior_work.stale(root)
    except build_prior_work.PriorWorkError as error:
        require(False, f"prior work: {error}")
    require(
        not stale_prior,
        "prior-work renderings differ: " + ", ".join(stale_prior)
        + "; run python3 tools/build_prior_work.py",
    )


def check_readme_pair(root: Path) -> None:
    """Both READMEs carry the same top-level heading count."""
    english = (root / "README.md").read_text(encoding="utf-8")
    french = (root / "README_FR.md").read_text(encoding="utf-8")
    english_headings = [line for line in english.splitlines() if line.startswith("## ")]
    french_headings = [line for line in french.splitlines() if line.startswith("## ")]
    require(
        len(english_headings) == len(french_headings),
        "README.md and README_FR.md top-level heading counts differ: "
        f"{len(english_headings)} and {len(french_headings)}",
    )


def check_standalone_papers(root: Path) -> None:
    for folder in ["paper", "extra"]:
        for path in (root / folder).glob("*.tex"):
            text = path.read_text(encoding="utf-8")
            for pattern in PAPER_EXTERNAL_REGISTRY_PATTERNS:
                require(
                    pattern not in text,
                    f"{path.relative_to(root)} references the external claim registry; papers must remain standalone",
                )


def check_extra_holds_papers(root: Path) -> None:
    """Keep extra/ for TeX papers; Markdown research notes live in docs/research/."""
    folder = root / "extra"
    if not folder.is_dir():
        return
    stray = sorted(
        path.relative_to(root).as_posix()
        for path in folder.rglob("*.md")
        if path.relative_to(folder).as_posix() != "README.md"
    )
    require(
        not stray,
        f"{', '.join(stray)}: extra/ holds TeX papers; move Markdown research notes to docs/research/",
    )


def dictionary_tokens(root: Path) -> set[str]:
    """Assumption tokens with a canonical dictionary row (backtick-quoted)."""
    text = (root / "claims" / "assumption_dictionary.md").read_text(encoding="utf-8")
    return set(re.findall(r"^\|\s*`([^`]+)`", text, re.MULTILINE))


def premise_ids(root: Path) -> set[str]:
    """Load the canonical premise IDs without inferring any semantic edges."""
    register = load_json(root / PREMISE_REGISTER_RELATIVE)
    rows = register.get("rows", [])
    require(
        isinstance(rows, list) and rows,
        f"{PREMISE_REGISTER_RELATIVE}: premise register has no rows",
    )
    identifiers = [row.get("id") for row in rows if isinstance(row, dict)]
    require(
        len(identifiers) == len(rows)
        and all(isinstance(identifier, str) and PREMISE_ID.fullmatch(identifier) for identifier in identifiers),
        f"{PREMISE_REGISTER_RELATIVE}: malformed premise ids",
    )
    require(
        len(identifiers) == len(set(identifiers)),
        f"{PREMISE_REGISTER_RELATIVE}: duplicate premise ids",
    )
    return set(identifiers)


def check_premise_dependencies(claim: dict, known_premises: set[str]) -> None:
    """Require an explicit direct-edge classification for every claim.

    `consumed` is an antecedent used by the scoped claim, `open` is a missing
    premise needed for the stronger physical claim named on that same row, and
    `boundary` is a premise whose derivation space the claim's no-go or audit
    directly constrains without consuming it. Claim-to-claim transitive ancestry
    stays in the dependency graph. Shared evidence and thematic similarity never
    create premise edges.
    """
    claim_id = claim["claim_id"]
    dependencies = claim.get("premise_dependencies")
    require(
        isinstance(dependencies, dict),
        f"{claim_id}: premise_dependencies must be an object",
    )
    classification = dependencies.get("classification")
    require(
        classification in PREMISE_DEPENDENCY_CLASSIFICATIONS,
        f"{claim_id}: premise dependency classification must be one of "
        f"{sorted(PREMISE_DEPENDENCY_CLASSIFICATIONS)}",
    )
    expected_keys = {"classification", *PREMISE_DEPENDENCY_FIELDS}
    if classification == "explicit_non_consumer":
        expected_keys.add("rationale")
    require(
        set(dependencies) == expected_keys,
        f"{claim_id}: premise_dependencies keys must equal {sorted(expected_keys)}",
    )

    fields: dict[str, list[str]] = {}
    for field in PREMISE_DEPENDENCY_FIELDS:
        values = dependencies[field]
        require(
            isinstance(values, list)
            and all(isinstance(value, str) and PREMISE_ID.fullmatch(value) for value in values),
            f"{claim_id}: premise_dependencies.{field} must be a list of PR-xx ids",
        )
        require(
            len(values) == len(set(values)),
            f"{claim_id}: premise_dependencies.{field} repeats an id",
        )
        require(
            values == sorted(values),
            f"{claim_id}: premise_dependencies.{field} must be sorted",
        )
        unknown = sorted(set(values) - known_premises)
        require(
            not unknown,
            f"{claim_id}: premise_dependencies.{field} cites unknown premise ids {unknown}",
        )
        fields[field] = values

    for index, left in enumerate(PREMISE_DEPENDENCY_FIELDS):
        for right in PREMISE_DEPENDENCY_FIELDS[index + 1 :]:
            overlap = sorted(set(fields[left]) & set(fields[right]))
            require(
                not overlap,
                f"{claim_id}: premise dependency fields {left} and {right} overlap at {overlap}",
            )

    edge_count = sum(len(values) for values in fields.values())
    if classification == "explicit_edges":
        require(
            edge_count > 0,
            f"{claim_id}: explicit_edges classification requires at least one direct premise edge",
        )
    else:
        require(
            edge_count == 0,
            f"{claim_id}: explicit_non_consumer must carry no premise edges",
        )
        rationale = dependencies["rationale"]
        require(
            isinstance(rationale, str) and rationale.strip(),
            f"{claim_id}: explicit_non_consumer requires a nonempty rationale",
        )


def premise_dependency_projection_sha256(claims: list[dict]) -> str:
    """Bind the independently reviewed direct-edge classification.

    The validator cannot infer scientific premise semantics from prose.  This
    projection pin makes every reviewed classification, edge type, and
    claim-specific non-consumer rationale an explicit review event instead of
    allowing a syntactically valid mutation to pass unnoticed.
    """
    projection = {
        claim["claim_id"]: claim["premise_dependencies"] for claim in claims
    }
    raw = (
        json.dumps(
            projection,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def check_non_consumer_rationale_uniqueness(claims: list[dict]) -> None:
    """Reject boilerplate shared across independently classified non-consumers."""
    rationale_owners: dict[str, list[str]] = {}
    for claim in claims:
        dependencies = claim["premise_dependencies"]
        if dependencies["classification"] != "explicit_non_consumer":
            continue
        normalized = " ".join(dependencies["rationale"].split()).casefold()
        rationale_owners.setdefault(normalized, []).append(claim["claim_id"])
    duplicates = {
        rationale: owners
        for rationale, owners in rationale_owners.items()
        if len(owners) > 1
    }
    require(
        not duplicates,
        "explicit_non_consumer rationales must be claim-specific; repeated "
        f"rationales={duplicates}",
    )


def evidence_medium(evidence: str) -> str | None:
    """Classify one evidence path's artifact medium from the path itself.

    The medium is a property of the artifact kind, so it is read from the
    suffix instead of from prose. An unlisted suffix returns None and its
    caller fails closed: a new artifact kind must be classified before it can
    back a claim.
    """
    suffix = Path(evidence.split("#", 1)[0]).suffix.casefold()
    return EVIDENCE_MEDIA_BY_SUFFIX.get(suffix)


def check_evidence_medium(claim: dict) -> bool:
    """Require a theorem proved in prose to say so on its own row.

    A theorem-asserting row whose entire evidence list is prose has no
    machine-checked or executed artifact behind it. That row must declare
    `proof_medium: paper_prose`, and every other row must not: the field would
    otherwise become an unread annotation that weakens a Lean-backed or
    run-backed row by attrition. The declaration is returned so the summary
    line can count the prose-medium theorem rows.
    """
    claim_id = claim["claim_id"]
    evidence = claim["evidence"]
    require(
        isinstance(evidence, list) and evidence,
        f"{claim_id}: evidence must be a nonempty list of artifact paths",
    )
    media: list[str] = []
    for item in evidence:
        medium = evidence_medium(item)
        require(
            medium is not None,
            f"{claim_id}: evidence artifact medium is unclassified; add the "
            f"suffix to EVIDENCE_MEDIA_BY_SUFFIX: {item}",
        )
        media.append(medium)

    prose_only = all(medium in PROSE_EVIDENCE_MEDIA for medium in media)
    theorem_asserting = claim["claim_class"] in THEOREM_ASSERTING_CLASSES
    declared = claim.get(PROOF_MEDIUM_FIELD)
    if theorem_asserting and prose_only:
        require(
            declared is not None,
            f"{claim_id}: claim_class {claim['claim_class']!r} asserts a theorem "
            f"whose evidence is prose only, so the row must declare "
            f"{PROOF_MEDIUM_FIELD}: paper_prose",
        )
        require(
            declared in PROOF_MEDIUM_VOCABULARY,
            f"{claim_id}: {PROOF_MEDIUM_FIELD} {declared!r} is not in the "
            f"controlled vocabulary {sorted(PROOF_MEDIUM_VOCABULARY)}",
        )
        return True
    require(
        declared is None,
        f"{claim_id}: {PROOF_MEDIUM_FIELD} is admissible only on a "
        f"theorem-asserting row whose evidence list is prose only; this row is "
        f"claim_class {claim['claim_class']!r} with media {sorted(set(media))}",
    )
    return False


def check_owner_medium(claim: dict) -> None:
    """Keep a non-paper owner from passing as a paper source.

    Owner files under the paper roots are standalone sources. An owner outside
    them is a protocol or run record, which is admissible only where the claim
    is the artifact itself, so that row must declare the medium and a
    paper-owned row must not.
    """
    claim_id = claim["claim_id"]
    owner = claim["owner_paper"]
    paper_owned = any(owner.startswith(root + "/") for root in PAPER_OWNER_ROOTS)
    declared = claim.get(OWNER_MEDIUM_FIELD)
    if paper_owned:
        require(
            declared is None,
            f"{claim_id}: owner {owner} is a paper source, so the row must not "
            f"declare {OWNER_MEDIUM_FIELD}",
        )
        return
    require(
        declared in OWNER_MEDIUM_VOCABULARY,
        f"{claim_id}: owner {owner} is outside the paper roots "
        f"{list(PAPER_OWNER_ROOTS)}, so the row must declare "
        f"{OWNER_MEDIUM_FIELD}: protocol_record",
    )
    require(
        claim["claim_class"] == "emitted_artifact",
        f"{claim_id}: {OWNER_MEDIUM_FIELD} {declared!r} is admissible only for "
        f"claim_class 'emitted_artifact', not {claim['claim_class']!r}",
    )


def check_topical_gate_owner_projection(
    claims: list[dict], *, check_policy_keys: bool
) -> None:
    """Bind the historical topical-lane policy to the rows it governs.

    The policy map and the registry rows are compared in both directions:
    a row declaring an owner the policy does not name fails, a row the policy
    names without the matching declaration fails, and every declared owner must
    appear in that row's own `gates`. On the live tree a policy key that names
    no registered claim fails as well, so deleting a governed row cannot leave
    the policy pointing at nothing.
    """
    field = REQUIRED_TOPICAL_GATE_OWNERS_FIELD
    for claim in claims:
        claim_id = claim["claim_id"]
        required = REQUIRED_V3_TOPIC_GATES_BY_CLAIM.get(claim_id)
        declared = claim.get(field)
        if required is None:
            require(
                declared is None,
                f"{claim_id}: {field} is declared while the historical topical-lane "
                f"policy names no gate for this claim",
            )
            continue
        require(
            isinstance(declared, list)
            and all(
                isinstance(gate, int) and not isinstance(gate, bool) and gate > 0
                for gate in declared
            ),
            f"{claim_id}: {field} must be a list of positive GitHub issue numbers",
        )
        require(
            declared == sorted(set(declared)),
            f"{claim_id}: {field} must be sorted and must not repeat",
        )
        require(
            set(declared) == set(required),
            f"{claim_id}: {field} {declared} does not equal the historical topical-lane "
            f"policy {sorted(required)}",
        )
        outside = sorted(set(declared) - set(claim["gates"]))
        require(
            not outside,
            f"{claim_id}: {field} names gates that are absent from `gates`: "
            f"{outside}",
        )
    if check_policy_keys:
        registered = {claim["claim_id"] for claim in claims}
        unregistered = sorted(set(REQUIRED_V3_TOPIC_GATES_BY_CLAIM) - registered)
        require(
            not unregistered,
            "the historical topical-lane policy names claims that the registry does "
            f"not contain: {unregistered}",
        )


def check_required_v3_topic_gates(claim: dict) -> None:
    """Keep each claim's historical topical lanes recorded."""
    claim_id = claim["claim_id"]
    required_gates = REQUIRED_V3_TOPIC_GATES_BY_CLAIM.get(claim_id, frozenset())
    missing = sorted(required_gates - set(claim["gates"]))
    require(
        not missing,
        f"{claim_id}: missing required historical topical lanes {missing}",
    )


def check_gates(claim: dict) -> None:
    """Validate the historical lanes of a claim and its promotion contract.

    `gates` records the lane issues under which a claim's open work was
    scoped. The V3 and V4 lanes are closed; each gate must be an open issue or
    a historical lane that `tools/closed_lanes.py` maps to the open issue that
    owns its remaining work, and a gate never names a superseded lane.
    """
    claim_id = claim["claim_id"]
    gates = claim["gates"]
    side_channels = sorted(
        key
        for key in claim
        if key != "gates"
        and "gate" in key.casefold()
        and key not in GATE_SIDE_CHANNEL_ALLOWLIST
    )
    require(
        not side_channels,
        f"{claim_id}: gate lists belong in `gates`, which every gate invariant "
        f"reads; found side-channel keys {side_channels}",
    )
    require(
        isinstance(gates, list)
        and all(isinstance(g, int) and not isinstance(g, bool) and g > 0 for g in gates),
        f"{claim_id}: gates must be a list of positive GitHub issue numbers",
    )
    require(
        len(gates) == len(set(gates)),
        f"{claim_id}: gates must not repeat GitHub issue numbers",
    )
    superseded = {
        gate: SUPERSEDED_LANE_SUCCESSORS[gate]
        for gate in gates
        if gate in SUPERSEDED_LANE_SUCCESSORS
    }
    require(
        not superseded,
        f"{claim_id}: gates name superseded lanes; record their successor lanes "
        f"{superseded}",
    )
    unknown = sorted(gate for gate in gates if not is_known(gate))
    require(
        not unknown,
        f"{claim_id}: gates {unknown} are neither open issues nor historical lanes "
        "with a live owner in tools/closed_lanes.py",
    )
    check_required_v3_topic_gates(claim)
    if claim["claim_class"] in PROMOTED_CLASSES:
        require(
            not gates,
            f"{claim_id}: claim_class {claim['claim_class']!r} asserts physical "
            f"establishment while gates {gates} still record open work",
        )


def check_wording(claim: dict) -> None:
    """Public wording must not be stronger than the claim's class."""
    if claim["claim_class"] in PROMOTED_CLASSES:
        return
    match = ESTABLISHMENT_WORDING.search(claim["statement"])
    if match is not None:
        raise SystemExit(
            f"{claim['claim_id']}: statement wording {match.group(0)!r} asserts "
            f"establishment but claim_class is {claim['claim_class']!r}"
        )


def main(root: Path = ROOT) -> None:
    registry = load_json(root / "claims" / "claim_registry.yaml")
    require(
        registry.get("schema_version") == CLAIM_REGISTRY_SCHEMA,
        f"claim registry schema_version must equal {CLAIM_REGISTRY_SCHEMA}",
    )
    require(
        registry.get("release_id") == release_id_from_tex(root),
        f"registry release_id {registry.get('release_id')!r} does not match paper/release_info.tex",
    )
    claims = registry.get("claims", [])
    require(isinstance(claims, list) and claims, "claim registry has no claims")
    check_standalone_papers(root)
    check_extra_holds_papers(root)
    defined_tokens = dictionary_tokens(root)
    known_premises = premise_ids(root)

    seen: set[str] = set()
    owner_paths: set[str] = set()
    prose_medium_claims: list[str] = []
    for claim in claims:
        missing = REQUIRED_CLAIM_FIELDS - set(claim)
        require(not missing, f"{claim.get('claim_id', '<missing>')}: missing fields {sorted(missing)}")
        claim_id = claim["claim_id"]
        require(claim_id not in seen, f"duplicate claim_id {claim_id}")
        seen.add(claim_id)
        require(claim["statement"].strip(), f"{claim_id}: empty statement")
        require(claim["assumptions"], f"{claim_id}: empty assumptions")
        require(claim["imported_results"], f"{claim_id}: empty imported_results")
        require(claim["oph_specific_delta"].strip(), f"{claim_id}: empty OPH delta")
        require(claim["falsifier"].strip(), f"{claim_id}: empty falsifier")
        require(
            claim["claim_class"] in CLAIM_CLASS_VOCABULARY,
            f"{claim_id}: claim_class {claim.get('claim_class')!r} is not in the "
            f"controlled vocabulary {sorted(CLAIM_CLASS_VOCABULARY)}",
        )
        check_gates(claim)
        check_premise_dependencies(claim, known_premises)
        check_wording(claim)
        if check_evidence_medium(claim):
            prose_medium_claims.append(claim_id)
        check_owner_medium(claim)
        owner = root / claim["owner_paper"]
        require(owner.exists(), f"{claim_id}: owner paper does not exist: {claim['owner_paper']}")
        owner_paths.add(claim["owner_paper"])
        undefined = sorted(set(claim["assumptions"]) - defined_tokens)
        require(
            not undefined,
            f"{claim_id}: assumption tokens without a canonical dictionary row: {undefined}",
        )
        for evidence in claim["evidence"]:
            require(
                evidence_is_available(root, evidence),
                f"{claim_id}: evidence path does not exist and pinned source URL is invalid: {evidence}",
            )

    check_non_consumer_rationale_uniqueness(claims)
    check_topical_gate_owner_projection(
        claims, check_policy_keys=root.resolve() == ROOT.resolve()
    )

    if root.resolve() == ROOT.resolve():
        actual_projection_sha256 = premise_dependency_projection_sha256(claims)
        require(
            actual_projection_sha256 == PREMISE_DEPENDENCY_PROJECTION_SHA256,
            "claim-registry premise-dependency projection differs from the "
            "independently reviewed classification; review the semantic change "
            "and update the pinned projection deliberately",
        )

    for matrix_name, required_columns, one_row_per_claim in [
        ("novelty_matrix.csv", {"claim_id", "closest_prior_work", "oph_specific_delta", "novelty_type", "falsifier"}, True),
        ("falsification_matrix.csv", {"claim_id", "mathematical_falsifier", "physical_identification_falsifier", "phenomenological_falsifier", "scope_if_false"}, False),
    ]:
        matrix_path = root / "claims" / matrix_name
        rows = load_csv(matrix_path)
        require(rows, f"{matrix_path}: no rows")
        require(required_columns.issubset(rows[0].keys()), f"{matrix_path}: missing required columns")
        check_exact_claim_id_projection(
            surface=str(matrix_path.relative_to(root)),
            projected_ids=[row["claim_id"] for row in rows],
            canonical_ids=seen,
            one_row_per_claim=one_row_per_claim,
        )
        if matrix_name == "falsification_matrix.csv":
            registry_scope = {claim["claim_id"]: claim["scope_if_false"] for claim in claims}
            drifted_scope = sorted(
                {
                    row["claim_id"]
                    for row in rows
                    if row["scope_if_false"] != registry_scope[row["claim_id"]]
                }
            )
            require(
                not drifted_scope,
                f"{matrix_path}: scope_if_false differs from the claim registry for {drifted_scope}",
            )
        if matrix_name == "novelty_matrix.csv":
            registry_novelty = {claim["claim_id"]: claim["novelty_type"] for claim in claims}
            drifted = sorted(
                row["claim_id"]
                for row in rows
                if row["novelty_type"] != registry_novelty[row["claim_id"]]
            )
            require(
                not drifted,
                f"{matrix_path}: novelty_type differs from the claim registry for {drifted}",
            )

    graph = load_json(root / "claims" / "dependency_graph.json")
    raw_nodes = graph.get("nodes", [])
    require(
        isinstance(raw_nodes, list) and all(isinstance(node, str) for node in raw_nodes),
        "claims/dependency_graph.json: nodes must be a list of claim ids",
    )
    check_exact_claim_id_projection(
        surface="claims/dependency_graph.json nodes",
        projected_ids=raw_nodes,
        canonical_ids=seen,
        one_row_per_claim=True,
    )
    nodes = set(raw_nodes)
    check_dependency_graph_edges(graph.get("edges", []), nodes)

    # The ledger fields and the rendered views are checked on the repository
    # itself; synthetic fixture registries exercise the structural gates above.
    if root.resolve() == ROOT.resolve():
        check_claim_views(root)
        check_readme_pair(root)

    gated = [claim for claim in claims if claim["gates"]]
    gate_count = len({gate for claim in gated for gate in claim["gates"]})
    owners = live_owners({gate for claim in gated for gate in claim["gates"]})
    classification_counts = Counter(
        claim["premise_dependencies"]["classification"] for claim in claims
    )
    premise_edge_count = sum(
        len(claim["premise_dependencies"][field])
        for claim in claims
        for field in PREMISE_DEPENDENCY_FIELDS
    )
    print(
        f"claim registry OK: {len(seen)} claims, {len(owner_paths)} owner papers, "
        f"{gate_count} historical lanes across {len(gated)} gated claims "
        f"(live owners {', '.join(f'#{n}' for n in owners)}), "
        f"{len(prose_medium_claims)} prose-medium theorem rows, "
        f"{premise_edge_count} direct premise edges "
        f"({classification_counts['explicit_edges']} classified consumers/boundaries, "
        f"{classification_counts['explicit_non_consumer']} explicit non-consumers)"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT)
