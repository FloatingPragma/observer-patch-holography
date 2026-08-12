"""Regression fixtures for the claim-registry validator (issue #512).

Each fixture builds a minimal repository tree, perturbs one defining
antecedent, and requires the validator to fail closed. The unperturbed
tree must pass, so a fixture cannot pass vacuously.
"""

import csv
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
    (root / "tracking" / "open_issues").mkdir(parents=True)
    (root / "tracking" / "open_issues" / "open_problem_ledger.json").write_text(
        json.dumps(
            {
                "open_issue_count": 1,
                "rows": [{"number": 42, "title": "open fixture gate"}],
                "closed_out_of_scope_records": [
                    {"number": 7, "title": "closed fixture issue"}
                ],
            }
        ),
        encoding="utf-8",
    )
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
                "claim_class": "conditional_implication",
                "gates": [42],
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


def test_duplicate_registry_object_key_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    path = tmp_path / "claims" / "claim_registry.yaml"
    text = path.read_text(encoding="utf-8")
    path.write_text(
        text.replace(
            '"statement": "Fixture claim.",',
            '"statement": "Fixture claim.", "statement": "Shadowed claim.",',
        ),
        encoding="utf-8",
    )

    with pytest.raises(SystemExit, match="duplicate object key 'statement'"):
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
    with pytest.raises(SystemExit, match="do not exactly match"):
        checker.main(tmp_path)


@pytest.mark.parametrize("matrix_name", ["novelty_matrix.csv", "falsification_matrix.csv"])
def test_registry_claim_missing_from_matrix_fails_closed(tmp_path, matrix_name):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"].append(
            {
                **r["claims"][0],
                "claim_id": "FIX-2",
                "statement": "Second canonical fixture claim.",
            }
        ),
    )
    # Deliberately leave the selected matrix with only FIX-1.  Make the other
    # projections complete so the selected omission is the first failure.
    other = (
        "falsification_matrix.csv"
        if matrix_name == "novelty_matrix.csv"
        else "novelty_matrix.csv"
    )
    path = tmp_path / "claims" / other
    rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
    rows.append({**rows[0], "claim_id": "FIX-2"})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    graph = tmp_path / "claims" / "dependency_graph.json"
    graph.write_text(
        json.dumps({"nodes": ["FIX-1", "FIX-2"], "edges": []}),
        encoding="utf-8",
    )

    with pytest.raises(
        SystemExit,
        match=rf"{matrix_name}: claim ids do not exactly match.*FIX-2",
    ):
        checker.main(tmp_path)


def test_registry_claim_missing_from_dependency_graph_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"].append(
            {
                **r["claims"][0],
                "claim_id": "FIX-2",
                "statement": "Second canonical fixture claim.",
            }
        ),
    )
    for matrix_name in ["novelty_matrix.csv", "falsification_matrix.csv"]:
        path = tmp_path / "claims" / matrix_name
        rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
        rows.append({**rows[0], "claim_id": "FIX-2"})
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    with pytest.raises(
        SystemExit,
        match=r"dependency_graph.json nodes: claim ids do not exactly match.*FIX-2",
    ):
        checker.main(tmp_path)


def test_duplicate_novelty_claim_id_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    matrix_name = "novelty_matrix.csv"
    path = tmp_path / "claims" / matrix_name
    rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
    rows.append(dict(rows[0]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(SystemExit, match=rf"{matrix_name}: duplicate claim ids.*FIX-1"):
        checker.main(tmp_path)


def test_multiple_falsification_rows_for_one_claim_are_permitted(tmp_path):
    write_fixture_repo(tmp_path)
    path = tmp_path / "claims" / "falsification_matrix.csv"
    rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
    rows.append(
        {
            **rows[0],
            "mathematical_falsifier": "independent second mathematical falsifier",
        }
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    checker.main(tmp_path)


def test_duplicate_dependency_node_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    graph = tmp_path / "claims" / "dependency_graph.json"
    graph.write_text(
        json.dumps({"nodes": ["FIX-1", "FIX-1"], "edges": []}),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="duplicate claim ids.*FIX-1"):
        checker.main(tmp_path)


def test_dependency_edge_endpoint_must_be_a_declared_node(tmp_path):
    write_fixture_repo(tmp_path)
    # An edge cannot smuggle in an endpoint omitted from the declared node list.
    graph = tmp_path / "claims" / "dependency_graph.json"
    graph.write_text(
        json.dumps(
            {
                "nodes": ["FIX-1"],
                "edges": [{"from": "FIX-1", "to": "FIX-GHOST", "role": "fixture"}],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="edge target is not a declared node"):
        checker.main(tmp_path)


def test_release_id_drift_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(tmp_path, lambda r: r.update(release_id="r-stale"))
    with pytest.raises(SystemExit, match="does not match"):
        checker.main(tmp_path)


def edit_snapshot(root: Path, mutate) -> None:
    path = root / "tracking" / "open_issues" / "open_problem_ledger.json"
    snapshot = json.loads(path.read_text(encoding="utf-8"))
    mutate(snapshot)
    path.write_text(json.dumps(snapshot), encoding="utf-8")


def test_closed_issue_referenced_as_open_gate_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(tmp_path, lambda r: r["claims"][0]["gates"].append(7))
    with pytest.raises(SystemExit, match="closed on GitHub but referenced"):
        checker.main(tmp_path)


def test_gate_missing_from_github_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(tmp_path, lambda r: r["claims"][0]["gates"].append(999))
    with pytest.raises(SystemExit, match="missing from the GitHub issue snapshot"):
        checker.main(tmp_path)


@pytest.mark.parametrize(
    ("claim_id", "remaining_gates", "missing_owner"),
    [
        ("OPH-DM-CONT", [], 742),
        ("OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX", [730], 743),
        ("OPH-YM-GAP", [743], 744),
        ("OPH-UNIFIED-TYPED-SPINE", [728, 729, 730, 740], 741),
        ("OPH-HIER-EW", [736, 740, 742], 745),
    ],
)
def test_named_v3_topic_gate_owner_cannot_silently_disappear(
    claim_id, remaining_gates, missing_owner
):
    with pytest.raises(
        SystemExit,
        match=rf"{claim_id}: missing required V3 topical gate owners.*{missing_owner}",
    ):
        checker.check_required_v3_topic_gates(
            {"claim_id": claim_id, "gates": remaining_gates}
        )


def test_claim_promoted_while_gate_open_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"][0].update(claim_class="physical_establishment"),
    )
    with pytest.raises(SystemExit, match="while gates .* are open"):
        checker.main(tmp_path)


def test_uncontrolled_claim_class_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"][0].update(claim_class="basically_proved"),
    )
    with pytest.raises(SystemExit, match="not in the controlled vocabulary"):
        checker.main(tmp_path)


def test_wording_stronger_than_class_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_registry(
        tmp_path,
        lambda r: r["claims"][0].update(
            statement="This branch result is physically established."
        ),
    )
    with pytest.raises(SystemExit, match="asserts establishment but claim_class"):
        checker.main(tmp_path)


def test_stale_snapshot_count_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_snapshot(tmp_path, lambda s: s.update(open_issue_count=58))
    with pytest.raises(SystemExit, match="hard-coded or stale count"):
        checker.main(tmp_path)


def test_snapshot_citing_unknown_claim_fails_closed(tmp_path):
    write_fixture_repo(tmp_path)
    edit_snapshot(
        tmp_path,
        lambda s: s["rows"][0].update(blocker="Blocked on OPH-GHOST-CLAIM."),
    )
    with pytest.raises(SystemExit, match="cites unknown claim id"):
        checker.main(tmp_path)


def test_live_repository_registry_passes():
    checker.main(REPO_ROOT)
