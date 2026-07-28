"""Regression and adversarial tests for the Workstream C rule tables."""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))
sys.path.insert(0, str(ROOT / "checkers"))

import rule_engine_a as engine_a  # noqa: E402
import rule_engine_b as engine_b  # noqa: E402
import check_rule_engines as checker  # noqa: E402
import enumerate_diagram_universe as enumerator  # noqa: E402
import vertex_format  # noqa: E402


@pytest.fixture(scope="module")
def table_a() -> dict:
    return engine_a.build_table()


@pytest.fixture(scope="module")
def table_b() -> dict:
    return engine_b.build_table()


def coefficient_of(table: dict, fields: list[str], structure: str) -> dict[tuple, str]:
    target = sorted(fields)
    for entry in table["entries"]:
        if entry["fields"] == target and entry["structure"] == structure:
            return {
                tuple(tuple(p) for p in m["powers"]): m["prefactor"]
                for m in entry["coefficient"]["monomials"]
            }
    raise AssertionError(f"record {target} {structure} is absent")


def test_engines_agree_exactly(table_a: dict, table_b: dict) -> None:
    assert table_a["entries"] == table_b["entries"]
    assert table_a["table_digest"] == table_b["table_digest"]
    assert table_a["entry_count"] == table_b["entry_count"] == 123


def test_emitted_files_match_builders(table_a: dict, table_b: dict) -> None:
    stored_a = json.loads(engine_a.OUT_PATH.read_text(encoding="utf-8"))
    stored_b = json.loads(engine_b.OUT_PATH.read_text(encoding="utf-8"))
    assert stored_a["table_digest"] == table_a["table_digest"]
    assert stored_b["table_digest"] == table_b["table_digest"]
    assert stored_a["entries"] == table_a["entries"]
    assert stored_b["entries"] == table_b["entries"]


def test_vector_masses(table_a: dict) -> None:
    assert coefficient_of(table_a, ["Wp", "Wm"], "vector_bilinear_mass") == {
        (("g2", 2), ("v", 2)): "1/4"
    }
    assert coefficient_of(table_a, ["W3", "B"], "vector_bilinear_mass") == {
        (("g1", 1), ("g2", 1), ("v", 2)): "-1/4"
    }


def test_tadpole_retained_minimum_not_imposed(table_a: dict) -> None:
    tadpole = coefficient_of(table_a, ["h"], "scalar_tadpole")
    assert tadpole == {(("mu2", 1), ("v", 1)): "1", (("lam", 1), ("v", 3)): "-1"}


def test_yang_mills_cubic_and_quartic(table_a: dict) -> None:
    cubic = coefficient_of(table_a, ["W3", "Wm", "Wp"], "yang_mills_three_point")
    assert cubic == {(("I", 1), ("g2", 1)): "-1"}
    assert coefficient_of(table_a, ["Wm", "Wm", "Wp", "Wp"], "yang_mills_four_point_12_34") == {
        (("g2", 2),): "1/2"
    }
    assert coefficient_of(table_a, ["Wm", "Wm", "Wp", "Wp"], "yang_mills_four_point_13_24") == {
        (("g2", 2),): "-1/2"
    }
    assert coefficient_of(table_a, ["W3", "W3", "Wm", "Wp"], "yang_mills_four_point_12_34") == {
        (("g2", 2),): "-1"
    }
    assert coefficient_of(table_a, ["W3", "W3", "Wm", "Wp"], "yang_mills_four_point_13_24") == {
        (("g2", 2),): "1"
    }


def test_goldstone_derivative_couplings_present(table_a: dict) -> None:
    assert coefficient_of(table_a, ["Gm", "Gp", "W3"], "scalar_scalar_gauge") == {
        (("I", 1), ("g2", 1)): "-1/2"
    }
    assert coefficient_of(table_a, ["G0", "h", "W3"], "scalar_scalar_gauge") == {
        (("g2", 1),): "1/2"
    }


def test_mixing_cancelled_by_solved_gauge_fixing(table_a: dict) -> None:
    assert not any(e["structure"] == "vector_scalar_mixing" for e in table_a["entries"])
    assert table_a["gauge_fixing"]["mixing_cancelled"] == [
        {"fields": ["B", "G0"], "structure": "vector_scalar_mixing"},
        {"fields": ["G0", "W3"], "structure": "vector_scalar_mixing"},
        {"fields": ["Gm", "Wp"], "structure": "vector_scalar_mixing"},
        {"fields": ["Gp", "Wm"], "structure": "vector_scalar_mixing"},
    ]


def test_ghost_and_goldstone_share_xi_poles(table_a: dict) -> None:
    g0 = coefficient_of(table_a, ["G0", "G0"], "scalar_bilinear_mass")
    assert g0[(("g2", 2), ("v", 2), ("xi2", 1))] == "-1/8"
    assert g0[(("g1", 2), ("v", 2), ("xi1", 1))] == "-1/8"
    z_ghost = coefficient_of(table_a, ["c3", "c3_bar"], "ghost_scalar_mass")
    assert z_ghost == {(("g2", 2), ("v", 2), ("xi2", 1)): "-1/4"}
    mixed = coefficient_of(table_a, ["c3_bar", "cB"], "ghost_scalar_mass")
    assert mixed == {(("g1", 1), ("g2", 1), ("v", 2), ("xi2", 1)): "1/4"}


def test_right_handed_hypercharges_flip_census_conjugates(table_a: dict) -> None:
    assert coefficient_of(table_a, ["uR_bar", "uR", "B"], "fermion_vector_current") == {
        (("g1", 1),): "2/3"
    }
    assert coefficient_of(table_a, ["eR_bar", "eR", "B"], "fermion_vector_current") == {
        (("g1", 1),): "-1"
    }


def test_charged_yukawa_and_dagger_records(table_a: dict) -> None:
    charged = coefficient_of(table_a, ["dL_bar", "uR", "Gm"], "fermion_scalar_yukawa")
    assert charged[(("Yu[1][2]", 1),)] == "1"
    dagger = coefficient_of(table_a, ["uR_bar", "dL", "Gp"], "fermion_scalar_yukawa")
    assert dagger[(("Yud[2][1]", 1),)] == "1"
    mass = coefficient_of(table_a, ["uL_bar", "uR"], "fermion_bilinear_mass")
    assert mass[(("Yu[3][1]", 1), ("sqrt2", 1), ("v", 1))] == "-1/2"


def test_checker_passes_and_controls_fire() -> None:
    table_a = json.loads(checker.TABLE_A.read_text(encoding="utf-8"))
    table_b = json.loads(checker.TABLE_B.read_text(encoding="utf-8"))
    universe = json.loads(checker.UNIVERSE.read_text(encoding="utf-8"))
    assert checker.verify_tables(table_a, table_b) == []
    assert checker.verify_universe(table_a, universe) == []
    controls = checker.run_controls(table_a, universe)
    for name, result in controls.items():
        assert result["expected_failure"] and result["failed"], name


def test_universe_covers_all_blocks() -> None:
    universe = enumerator.build_universe()
    for name, block in universe["external_blocks"].items():
        assert block["counts"]["bubbles"] > 0, name
        assert block["counts"]["seagulls"] > 0, name
        assert block["counts"]["tadpole_loops"] > 0, name
    stored = json.loads(enumerator.OUT_PATH.read_text(encoding="utf-8"))
    assert stored["universe_digest"] == universe["universe_digest"]


def test_conservation_holds_for_every_record(table_a: dict) -> None:
    for entry in table_a["entries"]:
        assert vertex_format.conservation_violations(entry["fields"]) == []


def test_format_rejects_charge_violation() -> None:
    assert vertex_format.conservation_violations(["Wp", "Wp", "W3"]) == ["electric charge 2"]


def test_monomial_canonical_folding() -> None:
    assert vertex_format.monomial(1, ("sqrt2", 2)) == ("2", ())
    assert vertex_format.monomial(1, ("sqrt2", -1)) == ("1/2", (("sqrt2", 1),))
    assert vertex_format.monomial(1, ("I", 2)) == ("-1", ())
    assert vertex_format.monomial(1, ("I", -1)) == ("-1", (("I", 1),))
    assert vertex_format.monomial(Fraction(1, 3), ("g2", 1), ("g2", 1)) == (
        "1/3", (("g2", 2),)
    )


def test_polynomial_cancels_exactly() -> None:
    combined = vertex_format.polynomial(
        vertex_format.monomial(1, ("g2", 1)),
        vertex_format.monomial(-1, ("g2", 1)),
    )
    assert combined == {"monomials": []}


def test_record_rejects_empty_polynomial() -> None:
    with pytest.raises(ValueError):
        vertex_format.record(["h"], {"monomials": []}, "scalar_tadpole")
