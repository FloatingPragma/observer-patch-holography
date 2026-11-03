"""Mutation controls for the Phase-0 fail-closed gates.

Each control changes an isolated temporary fixture and invokes the same
production checker used by the mandatory suite. The suite runs the canonical
public, provenance, and null-model checks before this file; self-contained
claim, public-surface, and release fixtures additionally prove their clean
case here. A control passes only when the named false-green mutation reaches
and is rejected by its intended gate. No control reads the network or changes
a tracked repository file.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from tools import check_claim_registry, check_reader_style


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_BLOCK_START = "<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->"
PUBLIC_BLOCK_END = "<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->"
CLAIM_CHECKER = ROOT / "tools/check_claim_registry.py"
EXTERNAL_CHECKER = ROOT / "tools/check_external_data_provenance.py"
NULL_CHECKER = ROOT / "tools/check_null_models.py"
PUBLIC_BUILDER = ROOT / "tools/build_public_quantitative_section.py"
PUBLIC_CHECKER = ROOT / "tools/check_public_surface_claims.py"
RELEASE_CHECKER = ROOT / "tools/check_github_release_channel.py"


def test_paper_style_gate_rejects_control_plane_identifiers_and_progress() -> None:
    samples = (
        "PR-65 labels the premise.",
        "PR-CC is listed in the paper.",
        "OL-C6 labels the observation.",
        "lane #743 tracks this construction.",
        "issue B19 owns this result.",
        "GitHub issue 730 owns this result.",
        "See https://github.com/example/project/issues/730.",
        "The physical attachment is work in progress.",
        "The continuum route stays open.",
        "The physical attachment is open.",
        "This is an open realization map.",
        "The status matrix lists seven exits.",
        "See the claim-status table.",
    )
    for sample in samples:
        assert any(
            pattern.search(sample)
            for pattern, _label in check_reader_style.PAPER_TRACKING_PATTERNS
        ), sample


def test_paper_style_gate_allows_scientific_open_and_identifier_lookalikes() -> None:
    samples = (
        "Pr-141 is a praseodymium isotope.",
        "The interval remains open.",
        "The channel stays open.",
        "The order remains partial.",
        "Each SIMD lane carries one word.",
        "Phase II is the deconfined phase.",
        "Issue 5 of the journal contains the erratum.",
        r"The cardinality is \# 3.",
        "The linguistic corpus contains 200 documents.",
    )
    for sample in samples:
        assert not any(
            pattern.search(sample)
            for pattern, _label in (
                check_reader_style.PROGRESS_PATTERNS
                + check_reader_style.PAPER_TRACKING_PATTERNS
            )
        ), sample


def test_prose_gate_scans_register_docs_and_holds_out_allowlist() -> None:
    scanned = {
        path.relative_to(ROOT).as_posix()
        for path in check_reader_style.iter_paths(
            check_reader_style.READER_GLOBS
            + check_reader_style.STATUS_GLOBS
            + check_reader_style.PAPER_GLOBS
            + check_reader_style.REGISTER_GLOBS
        )
    }
    for required in (
        "docs/registers/OBSERVATION_LEDGER_V3.md",
        "docs/registers/PREMISE_REGISTER_V3.md",
        "docs/EMPIRICAL_STATUS.md",
        "code/instrument_specs/OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md",
    ):
        assert (ROOT / required).is_file(), required
        assert required in scanned, required
    for allowlisted in (
        "docs/policies/STYLE_GUIDE.md",
        "docs/FROZEN_PREDICTION_LADDER.md",
    ):
        assert (ROOT / allowlisted).is_file(), allowlisted
        assert allowlisted not in scanned, allowlisted


def test_prose_gate_keeps_register_identifiers_legal_inside_register_rows() -> None:
    # The register docs are scanned for banned vocabulary alone.  A row id is
    # the subject of a register row, so the reader-identifier patterns are
    # confined to the reader globs.
    register_paths = {
        path.relative_to(ROOT).as_posix()
        for path in check_reader_style.iter_paths(check_reader_style.REGISTER_GLOBS)
    }
    reader_paths = {
        path.relative_to(ROOT).as_posix()
        for path in check_reader_style.iter_paths(check_reader_style.READER_GLOBS)
    }
    assert "docs/registers/OBSERVATION_LEDGER_V3.md" in register_paths
    assert "docs/registers/OBSERVATION_LEDGER_V3.md" not in reader_paths
    # The confinement is load-bearing: the register docs carry internal
    # identifiers the reader patterns match, and the gate passes on them.
    ledger = (ROOT / "docs/registers/OBSERVATION_LEDGER_V3.md").read_text(encoding="utf-8")
    assert any(
        pattern.search(ledger)
        for pattern, _label in check_reader_style.READER_IDENTIFIER_PATTERNS
    )


def test_main_paper_relevance_diagnostic_preserves_live_rg_routes() -> None:
    text = " ".join(
        (ROOT / "paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex")
        .read_text(encoding="utf-8")
        .split()
    )
    for required in (
        r"relevance alone does not generate \(\mathcal O\)",
        "vanishing operator mixing",
        "not an exhaustive protection theorem",
        "does not select chirality by itself",
    ):
        assert required in text, required
    for forbidden in (
        "only if symmetry forbids that direction",
        "is generated under refinement unless",
        "requires chiral matter content",
    ):
        assert forbidden not in text, forbidden


def test_main_paper_keeps_generalized_entropy_stationarity_as_interface() -> None:
    text = " ".join(
        (ROOT / "paper/tex_fragments/PAPER.tex")
        .read_text(encoding="utf-8")
        .split()
    )
    for required in (
        "proved finite MaxEnt envelope identities",
        "declared fixed-cap generalized-entropy stationarity interface",
        "does not, by itself, identify a variation of the state constraints with a geometric area variation",
        "declared joint state--geometry balance",
    ):
        assert required in text, required
    for forbidden in (
        "MaxEnt selection implies that for variations",
        "derived fixed-cap generalized-entropy stationarity theorem",
        "MaxEnt-selected fixed-cap generalized-entropy stationarity theorem",
    ):
        assert forbidden not in text, forbidden


def test_einstein_surfaces_keep_local_stress_and_ward_supplies_conditional() -> None:
    main_text = " ".join(
        (ROOT / "paper/tex_fragments/PAPER.tex")
        .read_text(encoding="utf-8")
        .split()
    )
    wrapper_text = " ".join(
        (
            ROOT
            / "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex"
        )
        .read_text(encoding="utf-8")
        .split()
    )
    supplement_text = " ".join(
        (ROOT / "paper/tex_fragments/DERIVATION_TECHNICAL_SUPPLEMENT_PORT.tex")
        .read_text(encoding="utf-8")
        .split()
    )
    publication_text = f"{main_text} {wrapper_text} {supplement_text}"
    for required in (
        "Charge compatibility and tensor tomography alone imply neither the face cancellation nor its continuum limit",
        "declared fixed-cap generalized-entropy stationarity interface",
        "conditional same-source local-stress and weak-Ward bridge",
        "conditional local null-stress reading under effective-locality and same-source matching",
        "assume the effective-locality and same-source bridge that identifies the proved half-line endpoint derivative",
        "finite null tomography does not supply them",
    ):
        assert required in publication_text, required
    for forbidden in (
        "Positive null translations and modular charges reconstruct the local conserved stress tensor",
        "Positive null translations and modular charges reconstruct a local conserved stress tensor",
        "closure packets construct the local conserved stress tensor",
        "exact identification of \\(P_\\Omega\\) with the local null-stress charge",
        "half-line generator/charge identification proved inside the null modular bridge itself",
        "half-line generator/charge identification",
        "null modular data reconstruct the stress tensor",
    ):
        assert forbidden not in publication_text, forbidden
    assert "derived fixed-cap generalized-entropy stationarity theorem" not in publication_text
    assert "generator/charge identification internally" not in publication_text

    novelty = (ROOT / "claims/novelty_matrix.csv").read_text(encoding="utf-8")
    falsification = (ROOT / "claims/falsification_matrix.csv").read_text(
        encoding="utf-8"
    )
    assert "T_ab is constructed from OPH modular charges" not in novelty
    dependency_graph = (ROOT / "claims/dependency_graph.json").read_text(
        encoding="utf-8"
    )
    assert "the MI linearity clause consumed by null tomography" not in dependency_graph
    assert "constructed modular stress channel" not in dependency_graph
    assert (
        "Every stated component implication and common-domain premise holds, but the typed Einstein composition fails"
        in falsification
    )
    assert "The D3--D4 Lorentz branch falls" not in next(
        line
        for line in falsification.splitlines()
        if line.startswith("OPH-GR-E2E-BRANCH-ENTRY,")
    )

    registry = json.loads(
        (ROOT / "claims/claim_registry.yaml").read_text(encoding="utf-8")
    )
    by_id = {row["claim_id"]: row for row in registry["claims"]}
    expected = [f"PR-{number}" for number in range(23, 29)]
    for claim_id in (
        "OPH-GR-D5A-ABSOLUTE-EINSTEIN",
        "OPH-GR-E2E-BRANCH-ENTRY",
    ):
        dependencies = by_id[claim_id]["premise_dependencies"]
        assert dependencies["classification"] == "explicit_edges"
        assert dependencies["consumed"] == expected
        assert "null_stress_bridge" in by_id[claim_id]["assumptions"]
        assert "weak_ward_conservation" in by_id[claim_id]["assumptions"]

    stress_assumptions = by_id["OPH-GR-D4C-LOCAL-STRESS"]["assumptions"]
    assert "null_stress_bridge" in stress_assumptions
    assert "weak_ward_conservation" in stress_assumptions


def _run(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *arguments],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def _combined(result: subprocess.CompletedProcess[str]) -> str:
    return result.stdout + result.stderr


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_claim_fixture(root: Path) -> Path:
    for relative in ("paper", "extra", "claims", "code", "tracking"):
        (root / relative).mkdir(parents=True, exist_ok=True)

    (root / "paper/release_info.tex").write_text(
        "\\newcommand{\\OPHPaperReleaseID}{r-test}\n",
        encoding="utf-8",
    )
    (root / "paper/owner.tex").write_text(
        "Fixture owner paper.\n",
        encoding="utf-8",
    )
    (root / "code/witness.py").write_text(
        "# fixture witness\n",
        encoding="utf-8",
    )
    (root / "claims/assumption_dictionary.md").write_text(
        "# Assumption dictionary\n\n"
        "| Token | Meaning |\n"
        "|---|---|\n"
        "| `DECLARED_FIXTURE` | Isolated mutation-control premise. |\n",
        encoding="utf-8",
    )

    registry_path = root / "claims/claim_registry.yaml"
    _write_json(
        registry_path,
        {
            "schema_version": 3,
            "release_id": "r-test",
            "claims": [
                {
                    "claim_id": "OPH-FIXTURE-GATE",
                    "statement": "A conditional fixture implication.",
                    "owner_paper": "paper/owner.tex",
                    "tier": "conditional",
                    "assumptions": ["DECLARED_FIXTURE"],
                    "imported_results": ["none"],
                    "oph_specific_delta": "Fixture delta.",
                    "novelty_type": "mutation control",
                    "evidence": ["code/witness.py"],
                    "falsifier": "The declared premise fails.",
                    "scope_if_false": "This fixture only.",
                    "status": "conditional_fixture",
                    "claim_class": "conditional_implication",
                    "gates": [1026],
                    "premise_dependencies": {
                        "classification": "explicit_edges",
                        "consumed": ["PR-01"],
                        "open": [],
                        "boundary": [],
                    },
                }
            ],
        },
    )
    (root / "claims/novelty_matrix.csv").write_text(
        "claim_id,closest_prior_work,oph_specific_delta,novelty_type,falsifier\n"
        "OPH-FIXTURE-GATE,none,fixture delta,mutation control,premise fails\n",
        encoding="utf-8",
    )
    (root / "claims/falsification_matrix.csv").write_text(
        "claim_id,mathematical_falsifier,physical_identification_falsifier,"
        "phenomenological_falsifier,scope_if_false\n"
        "OPH-FIXTURE-GATE,premise fails,not applicable,not applicable,"
        "This fixture only.\n",
        encoding="utf-8",
    )
    _write_json(
        root / "claims/dependency_graph.json",
        {"nodes": ["OPH-FIXTURE-GATE"], "edges": []},
    )
    _write_json(root / "tracking/premise_register.json", {"rows": [{"id": "PR-01"}]})
    return registry_path


def _copy(relative: str, target_root: Path) -> None:
    source = ROOT / relative
    target = target_root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def _write_public_fixture(root: Path) -> None:
    manifest_relative = Path("claims/public_surface_quantitative_claims.json")
    registry_relative = Path("claims/claim_registry.yaml")
    _copy(manifest_relative.as_posix(), root)
    _copy(registry_relative.as_posix(), root)

    manifest = json.loads(
        (root / manifest_relative).read_text(encoding="utf-8")
    )
    inputs: set[str] = set()
    for row in manifest["rows"]:
        inputs.add(row["producer"]["script"])
        inputs.add(row["producer"]["artifact"])
        reference = row.get("reference")
        if reference is not None:
            inputs.add(reference["artifact"])
        for support in row.get("supporting_artifacts", []):
            inputs.add(support["artifact"])
    for relative in sorted(inputs):
        _copy(relative, root)
    for surface in manifest["comparison_table_surfaces"]:
        _copy(surface["path"], root)
    ledger = root / "code/particles/runs/status/postdiction_ledger.json"
    ledger.parent.mkdir(parents=True, exist_ok=True)
    _write_json(
        ledger,
        {
            "artifact": "oph_postdiction_ledger",
            "sections": {
                "fixture": [
                    {"id": row["lane"]["row_id"]}
                    for surface in manifest["comparison_table_surfaces"]
                    for row in surface["rows"]
                    if row.get("lane", {}).get("kind") == "ledger"
                ]
            },
        },
    )

    for surface in manifest["surfaces"]:
        path = root / surface["path"]
        path.write_text(
            "# Mutation-control surface\n\n"
            f"{PUBLIC_BLOCK_START}\n"
            f"{PUBLIC_BLOCK_END}\n",
            encoding="utf-8",
        )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_release_fixture(root: Path) -> tuple[Path, Path, Path, Path, str]:
    paper = root / "paper/paper.pdf"
    book = root / "book/reverse-engineering-reality-book.pdf"
    paper.parent.mkdir(parents=True)
    book.parent.mkdir(parents=True)
    paper.write_bytes(b"%PDF-1.4\nfixture paper\n")
    book.write_bytes(b"%PDF-1.4\nfixture book\n")

    manifest = root / "paper/paper_release_manifest.json"
    _write_json(
        manifest,
        {
            "release_id": "r-test",
            "book": {
                "built_for_release_id": "r-test",
                "pdf_path": "book/reverse-engineering-reality-book.pdf",
                "sha256": _sha256(book),
                "size_bytes": book.stat().st_size,
            },
            "papers": {
                "fixture": {
                    "pdf_path": "paper/paper.pdf",
                    "sha256": _sha256(paper),
                    "size_bytes": paper.stat().st_size,
                }
            },
            "supplemental_papers": {},
            "extra_papers": {},
        },
    )
    # The book PDF is written into the fixture and recorded in the manifest, but it is not
    # a release asset: the book publishes on its own cadence to the reader-facing site.
    assets = []
    for path in (paper, manifest):
        assets.append(
            {
                "name": path.name,
                "digest": f"sha256:{_sha256(path)}",
                "size": path.stat().st_size,
            }
        )

    release = root / "release.json"
    latest = root / "latest.json"
    tag = root / "tag.json"
    commit = "a" * 40
    _write_json(
        release,
        {
            "tag_name": "r-test",
            "draft": False,
            "prerelease": False,
            "assets": assets,
        },
    )
    _write_json(latest, {"tag_name": "r-test"})
    _write_json(tag, {"tag_name": "r-test", "commit_sha": commit})
    return manifest, release, latest, tag, commit


def _release_command(
    root: Path,
    manifest: Path,
    release: Path,
    latest: Path,
    tag: Path,
    commit: str,
) -> tuple[str, ...]:
    return (
        str(RELEASE_CHECKER),
        "--repo-root",
        str(root),
        "--manifest",
        str(manifest),
        "--release-json",
        str(release),
        "--latest-json",
        str(latest),
        "--tag-json",
        str(tag),
        "--expected-commit",
        commit,
    )


def test_claim_gate_rejects_physical_promotion_with_open_work(
    tmp_path: Path,
) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    clean = _run(str(CLAIM_CHECKER), str(root))
    assert clean.returncode == 0, _combined(clean)

    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["claim_class"] = "physical_establishment"
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert "asserts physical establishment while gates [1026] still record open work" in _combined(
        mutant
    )


def test_claim_gate_rejects_gate_side_channel_keys(tmp_path: Path) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["gates"] = []
    registry["claims"][0]["github_gates"] = [42]
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert "side-channel keys ['github_gates']" in _combined(mutant)


def test_claim_gate_rejects_scope_drift_from_registry(tmp_path: Path) -> None:
    root = tmp_path / "claims"
    _write_claim_fixture(root)
    path = root / "claims/falsification_matrix.csv"
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "This fixture only.", "a paraphrase of the scope"
        ),
        encoding="utf-8",
    )
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert "scope_if_false differs from the claim registry for ['OPH-FIXTURE-GATE']" in (
        _combined(mutant)
    )


def test_claim_gate_rejects_novelty_type_drift_from_registry(tmp_path: Path) -> None:
    root = tmp_path / "claims"
    _write_claim_fixture(root)
    path = root / "claims/novelty_matrix.csv"
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "mutation control", "stronger unregistered label"
        ),
        encoding="utf-8",
    )
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert "novelty_type differs from the claim registry for ['OPH-FIXTURE-GATE']" in (
        _combined(mutant)
    )


@pytest.mark.parametrize("matrix", ["novelty_matrix.csv", "falsification_matrix.csv"])
@pytest.mark.parametrize("cell_break", ["\n", "\r\n"])
@pytest.mark.parametrize(
    ("mutation", "diagnostic"),
    [
        ("matching", None),
        ("overflow", "overflow CSV cells"),
        ("missing", "missing CSV cells"),
        ("empty", "empty CSV cells"),
        ("duplicate_header", "duplicate CSV headers"),
        ("empty_header", "missing or empty CSV header"),
        ("unterminated_quote", "invalid CSV"),
    ],
)
def test_claim_matrix_csv_structure_is_checked_before_field_use(
    tmp_path: Path, matrix: str, cell_break: str, mutation: str, diagnostic: str | None
) -> None:
    root = tmp_path / "claims"
    _write_claim_fixture(root)
    path = root / "claims" / matrix
    with path.open(encoding="utf-8", newline="") as handle:
        headers, body = list(csv.reader(handle))
    # Commas, escaped quotes, and multiline cells are legitimate scientific text.
    body[-1] += ', including "quoted" detail' + cell_break + 'and a second line'
    if matrix == "falsification_matrix.csv":
        # scope_if_false must equal the registry text, so the registry carries
        # the same scientific text and only the CSV structure is under test.
        registry_path = root / "claims" / "claim_registry.yaml"
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
        registry["claims"][0]["scope_if_false"] = body[-1]
        _write_json(registry_path, registry)
    if mutation == "overflow":
        body.append("silently discarded scope")
    elif mutation == "missing":
        body.pop()
    elif mutation == "empty":
        body[-1] = ""
    elif mutation == "duplicate_header":
        headers.append(headers[-1])
        body.append("overwritten scope")
    elif mutation == "empty_header":
        headers[-1] = ""
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(headers)
    if mutation == "unterminated_quote":
        output.write('"unterminated scientific text\n')
    else:
        writer.writerow(body)
    # newline="" keeps the multiline cell's embedded "\n" out of the platform's
    # line-ending translation, so the CSV carries the same bytes the registry
    # JSON escapes and the two agree on every platform.
    path.write_text(output.getvalue(), encoding="utf-8", newline="")

    result = _run(str(CLAIM_CHECKER), str(root))
    if diagnostic is None:
        assert result.returncode == 0, _combined(result)
    else:
        assert result.returncode != 0, _combined(result)
        assert matrix in _combined(result)
        assert diagnostic in _combined(result)


def _live_claims() -> list[dict[str, Any]]:
    registry = json.loads(
        (ROOT / "claims/claim_registry.yaml").read_text(encoding="utf-8")
    )
    return registry["claims"]


def test_claim_gate_requires_a_proof_medium_on_a_prose_only_theorem_row(
    tmp_path: Path,
) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["evidence"] = ["paper/owner.tex"]
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert "must declare proof_medium: paper_prose" in _combined(mutant)


def test_claim_gate_rejects_a_proof_medium_on_an_artifact_backed_row(
    tmp_path: Path,
) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["proof_medium"] = "paper_prose"
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert (
        "proof_medium is admissible only on a theorem-asserting row whose "
        "evidence list is prose only"
    ) in _combined(mutant)


@pytest.mark.parametrize(
    ("evidence", "diagnostic"),
    [
        (["code/witness.txt"], "evidence artifact medium is unclassified"),
        ([], "evidence must be a nonempty list of artifact paths"),
    ],
)
def test_claim_gate_rejects_unclassified_and_empty_evidence(
    tmp_path: Path, evidence: list[str], diagnostic: str
) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    (root / "code/witness.txt").write_text("fixture witness\n", encoding="utf-8")
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["evidence"] = evidence
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert diagnostic in _combined(mutant)


def test_live_prose_medium_rows_fail_closed_without_the_declaration() -> None:
    declared = [
        claim
        for claim in _live_claims()
        if check_claim_registry.PROOF_MEDIUM_FIELD in claim
    ]
    assert declared, "the registry declares no prose-medium theorem row"
    for claim in declared:
        stripped = {
            key: value
            for key, value in claim.items()
            if key != check_claim_registry.PROOF_MEDIUM_FIELD
        }
        with pytest.raises(SystemExit) as failure:
            check_claim_registry.check_evidence_medium(stripped)
        assert "must declare proof_medium: paper_prose" in str(failure.value)


def test_live_artifact_backed_rows_reject_the_prose_declaration() -> None:
    undeclared = [
        claim
        for claim in _live_claims()
        if check_claim_registry.PROOF_MEDIUM_FIELD not in claim
    ]
    assert undeclared
    for claim in undeclared:
        mutant = dict(claim)
        mutant[check_claim_registry.PROOF_MEDIUM_FIELD] = "paper_prose"
        with pytest.raises(SystemExit) as failure:
            check_claim_registry.check_evidence_medium(mutant)
        assert "admissible only on a theorem-asserting row" in str(failure.value)


def test_claim_gate_requires_owner_medium_for_a_nonpaper_owner(
    tmp_path: Path,
) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    (root / "code/protocol.md").write_text("Fixture protocol.\n", encoding="utf-8")
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["owner_paper"] = "code/protocol.md"
    registry["claims"][0]["claim_class"] = "emitted_artifact"
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert "must declare owner_medium: protocol_record" in _combined(mutant)


def test_claim_gate_rejects_owner_medium_on_a_paper_owned_row(
    tmp_path: Path,
) -> None:
    root = tmp_path / "claims"
    registry_path = _write_claim_fixture(root)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry["claims"][0]["owner_medium"] = "protocol_record"
    _write_json(registry_path, registry)
    mutant = _run(str(CLAIM_CHECKER), str(root))
    assert mutant.returncode != 0
    assert (
        "is a paper source, so the row must not declare owner_medium"
        in _combined(mutant)
    )


def _governed_claim(claims: list[dict[str, Any]]) -> dict[str, Any]:
    field = check_claim_registry.REQUIRED_TOPICAL_GATE_OWNERS_FIELD
    return next(claim for claim in claims if field in claim)


def _ungoverned_claim(claims: list[dict[str, Any]]) -> dict[str, Any]:
    field = check_claim_registry.REQUIRED_TOPICAL_GATE_OWNERS_FIELD
    return next(claim for claim in claims if field not in claim)


@pytest.mark.parametrize(
    ("edit", "diagnostic"),
    [
        ("drop_declaration", "must be a list of positive GitHub issue numbers"),
        ("declare_unnamed_row", "policy names no gate for this claim"),
        ("shrink_declaration", "does not equal the historical topical-lane policy"),
        ("drop_owner_from_gates", "names gates that are absent from `gates`"),
        (
            "delete_governed_row",
            "policy names claims that the registry does not contain",
        ),
    ],
)
def test_topical_owner_projection_rejects_one_sided_edits(
    edit: str, diagnostic: str
) -> None:
    field = check_claim_registry.REQUIRED_TOPICAL_GATE_OWNERS_FIELD
    claims = [dict(claim) for claim in _live_claims()]
    governed = _governed_claim(claims)
    if edit == "drop_declaration":
        del governed[field]
    elif edit == "declare_unnamed_row":
        _ungoverned_claim(claims)[field] = list(governed[field])
    elif edit == "shrink_declaration":
        governed[field] = governed[field][:-1]
    elif edit == "drop_owner_from_gates":
        governed["gates"] = [
            gate for gate in governed["gates"] if gate not in governed[field]
        ]
    else:
        claims = [claim for claim in claims if claim is not governed]
    with pytest.raises(SystemExit) as failure:
        check_claim_registry.check_topical_gate_owner_projection(
            claims, check_policy_keys=True
        )
    assert diagnostic in str(failure.value)


def test_live_topical_owner_projection_matches_the_registry() -> None:
    field = check_claim_registry.REQUIRED_TOPICAL_GATE_OWNERS_FIELD
    claims = _live_claims()
    declared = {
        claim["claim_id"]: claim[field] for claim in claims if field in claim
    }
    assert declared == {
        claim_id: sorted(gates)
        for claim_id, gates in (
            check_claim_registry.REQUIRED_V3_TOPIC_GATES_BY_CLAIM.items()
        )
    }
    check_claim_registry.check_topical_gate_owner_projection(
        claims, check_policy_keys=True
    )


def test_external_provenance_gate_rejects_a_forged_artifact_pin(
    tmp_path: Path,
) -> None:
    registry_path = tmp_path / "external_registry.json"
    registry = json.loads(
        (ROOT / "code/audit/external_data_provenance_registry.json").read_text(
            encoding="utf-8"
        )
    )
    registry["entries"][0]["artifact_sha256"] = "0" * 64
    _write_json(registry_path, registry)

    mutant = _run(str(EXTERNAL_CHECKER), "--registry", str(registry_path))
    assert mutant.returncode != 0
    assert "$.entries[0] SHA-256 mismatch" in _combined(mutant)


def test_public_surface_gate_rejects_a_manual_generated_result(
    tmp_path: Path,
) -> None:
    root = tmp_path / "public"
    root.mkdir()
    _write_public_fixture(root)
    generated = _run(str(PUBLIC_BUILDER), "--root", str(root))
    assert generated.returncode == 0, _combined(generated)

    readme = root / "README.md"
    marker = PUBLIC_BLOCK_END
    text = readme.read_text(encoding="utf-8")
    assert marker in text
    readme.write_text(
        text.replace(marker, "false-green result\n" + marker, 1),
        encoding="utf-8",
    )

    mutant = _run(str(PUBLIC_CHECKER), "--root", str(root))
    assert mutant.returncode != 0
    assert "generated quantitative claim block is stale" in _combined(mutant)


def test_null_model_gate_rejects_a_tampered_scorecard(
    tmp_path: Path,
) -> None:
    scorecard = tmp_path / "null_model_scorecard.md"
    canonical = ROOT / "tracking/null_model_scorecard.md"
    shutil.copy2(canonical, scorecard)
    scorecard.write_text(
        scorecard.read_text(encoding="utf-8") + "\nfalse-green result\n",
        encoding="utf-8",
    )

    mutant = _run(
        str(NULL_CHECKER),
        "--root",
        str(ROOT),
        "--output",
        str(scorecard),
        "--check",
    )
    assert mutant.returncode != 0
    assert "has drifted from the null-model inputs" in _combined(mutant)


def test_offline_release_gate_rejects_byte_and_tag_false_greens(
    tmp_path: Path,
) -> None:
    root = tmp_path / "release"
    manifest, release, latest, tag, commit = _write_release_fixture(root)
    command = _release_command(root, manifest, release, latest, tag, commit)
    clean = _run(*command)
    assert clean.returncode == 0, _combined(clean)

    payload = json.loads(release.read_text(encoding="utf-8"))
    paper_asset = next(
        asset for asset in payload["assets"] if asset["name"] == "paper.pdf"
    )
    paper_asset["digest"] = "sha256:" + ("0" * 64)
    _write_json(release, payload)
    byte_mutant = _run(*command)
    assert byte_mutant.returncode != 0
    assert "public digest" in _combined(byte_mutant)

    paper_asset["digest"] = f"sha256:{_sha256(root / 'paper/paper.pdf')}"
    _write_json(release, payload)
    _write_json(tag, {"tag_name": "r-test", "commit_sha": "b" * 40})
    tag_mutant = _run(*command)
    assert tag_mutant.returncode != 0
    assert "public release tag commit differs from the requested commit" in _combined(
        tag_mutant
    )


def test_readme_pair_requires_the_charged_lepton_limitation(tmp_path: Path) -> None:
    english = tmp_path / "README.md"
    french = tmp_path / "README_FR.md"
    english.write_text(
        "# README\n\n## Scope\n\n"
        "No tested mass law reproduces the measured charged-lepton hierarchy.\n",
        encoding="utf-8",
    )
    french.write_text(
        "# README\n\n## Portée\n\n"
        "Aucune loi de masse testée ne reproduit la hiérarchie mesurée des leptons chargés.\n",
        encoding="utf-8",
    )
    check_claim_registry.check_readme_pair(tmp_path)

    french.write_text(
        "# README\n\n## Portée\n\n## Détails\n\n"
        "Aucune loi de masse testée ne reproduit la hiérarchie mesurée des leptons chargés.\n",
        encoding="utf-8",
    )
    with pytest.raises(SystemExit) as failure:
        check_claim_registry.check_readme_pair(tmp_path)
    assert "top-level heading counts" in str(failure.value)

    french.write_text("# README\n\n## Portée\n\nAucune limite.\n", encoding="utf-8")
    with pytest.raises(SystemExit) as failure:
        check_claim_registry.check_readme_pair(tmp_path)
    assert "charged-lepton mass-hierarchy limitation" in str(failure.value)
