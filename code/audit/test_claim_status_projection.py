"""Fail-closed tests for the generated #512 cross-surface status projection."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_scoreboard.py"

spec = importlib.util.spec_from_file_location("build_scoreboard", BUILDER)
scoreboard = importlib.util.module_from_spec(spec)
sys.modules["build_scoreboard"] = scoreboard
spec.loader.exec_module(scoreboard)


def fixture_sources():
    registry = {
        "release_id": "r-test",
        "claims": [
            {
                "claim_id": "FIX-1",
                "claim_class": "conditional_implication",
                "status": "proved_on_declared_branch",
                "gates": [42],
            }
        ],
    }
    snapshot = {
        "repo": "example/project",
        "rows": [
            {"number": 42, "title": "open gate"},
            {"number": 554, "title": "declare selector menus"},
        ],
    }
    accounting = {
        "schema_version": 1,
        "physical_identifications": [
            {
                "identification_id": "fixture_physical_identification",
                "label": "fixture mathematical object / physical object",
                "status": "undischarged",
                "claim_ids": ["FIX-1"],
                "source_anchors": ["Fixture assumption"],
                "blocking_issues": [
                    {
                        "number": 42,
                        "disposition": "open_work_item",
                        "role": "fixture physical source binding",
                    }
                ],
            }
        ],
        "selector_menus": [
            {
                "selector_id": "fixture_selector_inventory",
                "label": "fixture selector menu inventory",
                "source": "docs/SELECTION_LEDGER.md",
                "status": "undeclared",
                "menu_size": None,
                "blocking_issues": [554],
            }
        ],
        "compression_bound": {
            "status": "not_computable",
            "reason": "menu_sizes_undeclared",
            "blocking_issues": [554],
            "p_acc_upper_bound": None,
        },
    }
    return registry, snapshot, accounting


@pytest.mark.parametrize(
    ("surface", "mutation"),
    [
        (
            "claim class",
            lambda r, s, a: r["claims"][0].update(claim_class="branch_entry"),
        ),
        ("descriptive status", lambda r, s, a: r["claims"][0].update(status="open")),
        ("claim gate", lambda r, s, a: r["claims"][0].update(gates=[])),
        (
            "live issue state",
            lambda r, s, a: s["rows"].append(
                {"number": 99, "title": "additional open gate"}
            ),
        ),
        ("release", lambda r, s, a: r.update(release_id="r-next")),
        (
            "physical identification",
            lambda r, s, a: a["physical_identifications"][0].update(
                label="changed identification"
            ),
        ),
        (
            "selector menu",
            lambda r, s, a: a["selector_menus"][0].update(
                label="changed selector menu"
            ),
        ),
    ],
)
def test_projection_digest_binds_every_status_input(surface, mutation):
    registry, snapshot, accounting = fixture_sources()
    before = scoreboard.projection_digest(
        scoreboard.canonical_projection_payload(registry, snapshot, accounting)
    )
    mutation(registry, snapshot, accounting)
    after = scoreboard.projection_digest(
        scoreboard.canonical_projection_payload(registry, snapshot, accounting)
    )
    assert after != before, surface


def test_projection_replacement_is_exact_and_idempotent():
    registry, snapshot, accounting = fixture_sources()
    projection = scoreboard.render_projection(registry, snapshot, accounting)
    seed = (
        "# Surface\n\n"
        f"{scoreboard.PROJECTION_START}\nold\n{scoreboard.PROJECTION_END}\n"
    )
    updated = scoreboard.replace_projection(seed, projection, Path("surface.md"))
    assert projection in updated
    assert scoreboard.replace_projection(
        updated, projection, Path("surface.md")
    ) == updated


@pytest.mark.parametrize(
    "bad",
    [
        "# no markers\n",
        (
            f"{scoreboard.PROJECTION_START}\n{scoreboard.PROJECTION_END}\n"
            f"{scoreboard.PROJECTION_START}\n{scoreboard.PROJECTION_END}\n"
        ),
    ],
)
def test_missing_or_duplicate_projection_markers_fail_closed(bad):
    registry, snapshot, accounting = fixture_sources()
    projection = scoreboard.render_projection(registry, snapshot, accounting)
    with pytest.raises(SystemExit, match="exactly one generated claim-status block"):
        scoreboard.replace_projection(bad, projection, Path("surface.md"))


def test_proof_spine_and_compression_scorecard_carry_identical_live_projection():
    registry, snapshot, accounting = scoreboard.source_documents()
    projection = scoreboard.render_projection(registry, snapshot, accounting)
    for relative in ["docs/PROOF_SPINE.md", "docs/COMPRESSION_SCORECARD.md"]:
        path = REPO_ROOT / relative
        current = path.read_text(encoding="utf-8")
        assert scoreboard.replace_projection(
            current, projection, Path(relative)
        ) == current


def test_scoreboard_leads_with_enumerated_physical_identifications():
    registry, snapshot, accounting = fixture_sources()
    page = scoreboard.render(registry, snapshot, accounting)
    physical = page.index("## Undischarged physical identifications")
    compression = page.index("## Compression-bound state")
    class_table = page.index("| Class | Claims |")
    assert physical < compression < class_table
    assert (
        "1. **`fixture_physical_identification`: "
        "fixture mathematical object / physical object.**"
    ) in page
    assert "[#42](https://github.com/example/project/issues/42)" in page


def test_undeclared_menu_renders_exact_noncomputable_state():
    registry, snapshot, accounting = fixture_sources()
    page = scoreboard.render(registry, snapshot, accounting)
    state = "compression_bound: not_computable - menu sizes undeclared (#554)"
    assert f"`{state}`" in page
    assert state in scoreboard.render_projection(registry, snapshot, accounting)


def test_live_physical_registry_contains_required_audit_rows_and_dispositions():
    _, _, accounting = scoreboard.source_documents()
    by_id = {
        row["identification_id"]: row
        for row in accounting["physical_identifications"]
    }
    assert {
        "alpha_in_thomson",
        "beta_ew_common_load_carrier",
        "capacity_horizon",
        "de_sitter_capacity_transfer_shock_sign",
        "port_bracket_physical_current",
        "repair_generator_yang_mills_hamiltonian",
        "g_si_clock",
    } <= set(by_id)

    alpha_blockers = {
        row["number"]: row["disposition"]
        for row in by_id["alpha_in_thomson"]["blocking_issues"]
    }
    assert alpha_blockers == {
        318: "open_work_item",
        545: "open_work_item",
        425: "parked_computational_blocker",
    }
    assert (
        by_id["repair_generator_yang_mills_hamiltonian"]["source_anchors"]
        == ["Assumption 20"]
    )
    de_sitter_blockers = {
        row["number"]: row["disposition"]
        for row in by_id["de_sitter_capacity_transfer_shock_sign"][
            "blocking_issues"
        ]
    }
    assert de_sitter_blockers == {
        334: "open_work_item",
        505: "open_work_item",
        589: "open_work_item",
        595: "open_work_item",
        608: "open_work_item",
    }
    assert set(
        by_id["de_sitter_capacity_transfer_shock_sign"]["source_anchors"]
    ) == {
        "de_Sitter_fixed_total_capacity_horizon_dictionary",
        "capacity_ledger_to_observer_mass_dictionary",
        "de_Sitter_shock_coefficient_scale_dictionary",
        "Einstein_branch_D5",
        "independent_physical_scale_receipt",
    }


def test_de_sitter_claim_split_preserves_status_boundaries():
    registry, _, _ = scoreboard.source_documents()
    by_id = {claim["claim_id"]: claim for claim in registry["claims"]}
    expected_classes = {
        "OPH-GR-DS-MU2-IDENTITY": "conditional_implication",
        "OPH-GR-DS-CAPACITY-TRANSFER": "conditional_implication",
        "OPH-GR-DS-SHOCK-SIGN-ATTACHMENT": "conditional_implication",
        "OPH-GR-DS-DISCRETE-SHOCK-SPECTRUM": "conditional_implication",
        "OPH-REPAIR-SHOCK-DOMAIN-BOUNDARY": "declared_structure",
    }
    assert {
        claim_id: by_id[claim_id]["claim_class"]
        for claim_id in expected_classes
    } == expected_classes
    assert not any(
        claim["claim_class"] == "physical_establishment"
        for claim in registry["claims"]
    )
    assert {
        "DS-GAUGE",
        "DS-LAPLACIAN",
    } == set(
        by_id["OPH-GR-DS-DISCRETE-SHOCK-SPECTRUM"]["assumptions"]
    )
    assert (
        "no unconstrained interior maximum"
        in by_id["OPH-GR-DS-CAPACITY-TRANSFER"]["statement"]
    )
    discrete = by_id["OPH-GR-DS-DISCRETE-SHOCK-SPECTRUM"]
    assert discrete["gates"] == [608]
    assert "draft_only_no_live_gate" not in discrete["status"]


def test_numeric_p_acc_is_forbidden_while_any_menu_is_undeclared():
    registry, snapshot, accounting = fixture_sources()
    accounting["compression_bound"]["p_acc_upper_bound"] = 0.5
    with pytest.raises(SystemExit, match="numeric P_acc is forbidden"):
        scoreboard.render(registry, snapshot, accounting)


def test_undeclared_menu_cannot_hide_a_numeric_size():
    registry, snapshot, accounting = fixture_sources()
    accounting["selector_menus"][0]["menu_size"] = 3
    with pytest.raises(SystemExit, match="cannot carry a numeric menu_size"):
        scoreboard.render(registry, snapshot, accounting)


def test_declared_menu_requires_a_positive_integer_size():
    registry, snapshot, accounting = fixture_sources()
    menu = accounting["selector_menus"][0]
    menu.update(status="declared", blocking_issues=[])
    accounting["compression_bound"].update(
        status="not_computed",
        reason="calculation_not_run",
        blocking_issues=[],
    )
    with pytest.raises(SystemExit, match="requires a positive integer menu_size"):
        scoreboard.render(registry, snapshot, accounting)


def test_physical_identification_requires_live_canonical_claim_and_issue():
    registry, snapshot, accounting = fixture_sources()
    accounting["physical_identifications"][0]["claim_ids"] = ["FIX-GHOST"]
    with pytest.raises(SystemExit, match="unknown canonical claim ids"):
        scoreboard.render(registry, snapshot, accounting)

    registry, snapshot, accounting = fixture_sources()
    snapshot["rows"] = [
        row for row in snapshot["rows"] if row["number"] != 42
    ]
    with pytest.raises(SystemExit, match="absent from the canonical open-issue ledger"):
        scoreboard.render(registry, snapshot, accounting)


def test_injectable_root_checker_rejects_mutated_generated_page(tmp_path):
    registry, snapshot, accounting = fixture_sources()
    for relative in ["claims", "tracking/open_issues", "docs"]:
        (tmp_path / relative).mkdir(parents=True, exist_ok=True)
    (tmp_path / "claims" / "claim_registry.yaml").write_text(
        json.dumps(registry),
        encoding="utf-8",
    )
    (
        tmp_path / "tracking" / "open_issues" / "open_problem_ledger.json"
    ).write_text(json.dumps(snapshot), encoding="utf-8")
    (tmp_path / "claims" / "physical_identification_registry.json").write_text(
        json.dumps(accounting),
        encoding="utf-8",
    )
    (tmp_path / "docs" / "SELECTION_LEDGER.md").write_text(
        "# Fixture selector ledger\n",
        encoding="utf-8",
    )
    (tmp_path / "docs" / "PROOF_SPINE.md").write_text(
        "# Fixture proof spine\n\n"
        f"{scoreboard.PROJECTION_START}\n"
        f"{scoreboard.PROJECTION_END}\n",
        encoding="utf-8",
    )

    scoreboard.write_generated_surfaces(tmp_path)
    scoreboard.check_generated_surfaces(tmp_path)

    generated = tmp_path / "tracking" / "claims_scoreboard.md"
    generated.write_text(
        generated.read_text(encoding="utf-8").replace(
            "Status: `undischarged`",
            "Status: `false-green`",
            1,
        ),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="has drifted from its sources"):
        scoreboard.check_generated_surfaces(tmp_path)
