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

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from tools import check_reader_style


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


def test_main_paper_relevance_diagnostic_preserves_live_rg_routes() -> None:
    text = " ".join(
        (ROOT / "paper/tex_fragments/PAPER.tex")
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
    (root / "code/witness.txt").write_text(
        "fixture witness\n",
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
                    "evidence": ["code/witness.txt"],
                    "falsifier": "The declared premise fails.",
                    "scope_if_false": "This fixture only.",
                    "status": "conditional_fixture",
                    "claim_class": "conditional_implication",
                    "gates": [42],
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
        "fixture only\n",
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
    assert "asserts physical establishment while gates [42] are open" in _combined(
        mutant
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
