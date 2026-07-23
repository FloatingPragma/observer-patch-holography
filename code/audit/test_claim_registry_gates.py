"""Regression fixtures for the claim-registry validator (issue #512).

Each fixture builds a minimal repository tree, perturbs one defining
antecedent, and requires the validator to fail closed. The unperturbed
tree must pass, so a fixture cannot pass vacuously.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKER = REPO_ROOT / "tools" / "check_claim_registry.py"

spec = importlib.util.spec_from_file_location("check_claim_registry", CHECKER)
checker = importlib.util.module_from_spec(spec)
sys.modules["check_claim_registry"] = checker
spec.loader.exec_module(checker)


def write_fixture_repo(root: Path) -> None:
    (root / "paper").mkdir()
    (root / "claims").mkdir()
    (root / "extra").mkdir()
    (root / "code").mkdir()
    (root / "paper" / "release_info.tex").write_text(
        "\\newcommand{\\OPHPaperReleaseID}{r-test}\n", encoding="utf-8"
    )
    (root / "paper" / "owner.tex").write_text("owner paper\n", encoding="utf-8")
    (root / "code" / "witness.py").write_text("print('witness')\n", encoding="utf-8")
    (root / "claims" / "assumption_dictionary.md").write_text(
        "# Dictionary\n\n| Assumption | Meaning | Primary owner |\n|---|---|---|\n"
        "| `declared_token` | A declared fixture assumption. | Owner paper |\n",
        encoding="utf-8",
    )
    registry = {
        "schema_version": 1,
        "release_id": "r-test",
        "claims": [
            {
                "claim_id": "FIX-1",
                "statement": "Fixture claim.",
                "owner_paper": "paper/owner.tex",
                "tier": "definition",
                "assumptions": ["declared_token"],
                "imported_results": ["none"],
                "oph_specific_delta": "Fixture delta.",
                "novelty_type": "fixture",
                "evidence": ["code/witness.py"],
                "falsifier": "Fixture falsifier.",
                "scope_if_false": "Fixture scope.",
                "status": "declared_basis",
            }
        ],
    }
    (root / "claims" / "claim_registry.yaml").write_text(
        json.dumps(registry), encoding="utf-8"
    )
    (root / "claims" / "novelty_matrix.csv").write_text(
        "claim_id,closest_prior_work,oph_specific_delta,novelty_type,falsifier\n"
        "FIX-1,none,delta,fixture,falsifier\n",
        encoding="utf-8",
    )
    (root / "claims" / "falsification_matrix.csv").write_text(
        "claim_id,mathematical_falsifier,physical_identification_falsifier,"
        "phenomenological_falsifier,scope_if_false\n"
        "FIX-1,m,p,ph,scope\n",
        encoding="utf-8",
    )
    (root / "claims" / "dependency_graph.json").write_text(
        json.dumps({"nodes": ["FIX-1"], "edges": []}), encoding="utf-8"
    )


def edit_registry(root: Path, mutate) -> None:
    path = root / "claims" / "claim_registry.yaml"
    registry = json.loads(path.read_text(encoding="utf-8"))
    mutate(registry)
    path.write_text(json.dumps(registry), encoding="utf-8")


def test_clean_fixture_passes(tmp_path):
    write_fixture_repo(tmp_path)
    checker.main(tmp_path)


def test_broken_evidence_path_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"][0]["evidence"].append("code/missing_witness.py"),
    )
    with pytest.raises(SystemExit, match="evidence path does not exist"):
        checker.main(tmp_path)


def test_undefined_assumption_token_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"][0]["assumptions"].append("undeclared_token"),
    )
    with pytest.raises(SystemExit, match="without a canonical dictionary row"):
        checker.main(tmp_path)


def test_unknown_dependency_node_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    graph = tmp_path / "claims" / "dependency_graph.json"
    graph.write_text(
        json.dumps({"nodes": ["FIX-1", "FIX-GHOST"], "edges": []}),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="unknown nodes"):
        checker.main(tmp_path)


def test_release_id_drift_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(tmp_path, lambda r: r.update(release_id="r-stale"))
    with pytest.raises(SystemExit, match="does not match"):
        checker.main(tmp_path)


def test_live_repository_registry_passes():
    checker.main(REPO_ROOT)
