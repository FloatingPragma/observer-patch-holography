"""Regression and adversarial tests for the Workstream D counterterm packet."""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))
sys.path.insert(0, str(ROOT / "checkers"))

import counterterm_producer as producer  # noqa: E402
import check_counterterms as checker  # noqa: E402


def load_packet() -> dict:
    return json.loads(producer.OUT_PATH.read_text(encoding="utf-8"))


def test_packet_is_deterministic() -> None:
    packet = producer.build_packet()
    stored = load_packet()
    assert packet["ct_terms_digest"] == stored["ct_terms_digest"]
    assert packet["ct_matrix"] == stored["ct_matrix"]
    assert packet["packet_sha256"] == stored["packet_sha256"]


def test_rank_and_global_symmetry_nullspace() -> None:
    packet = load_packet()
    assert packet["ct_matrix"]["columns"] == 156
    assert packet["ct_matrix"]["rank"] == 154
    null = packet["ct_matrix"]["nullspace"]
    assert len(null) == 2
    supports = [set(v) for v in null]
    quark = {f"{m}{d}[{i}][{i}]" for m in ("dZQL", "dZuR", "dZdR") for d in ("", "d") for i in (1, 2, 3)}
    lepton = {f"{m}{d}[{i}][{i}]" for m in ("dZLL", "dZeR") for d in ("", "d") for i in (1, 2, 3)}
    assert sorted(map(len, supports)) == [12, 18]
    assert any(s == quark for s in supports)
    assert any(s == lepton for s in supports)


def test_gauge_poles_derive_from_census_betas() -> None:
    packet = load_packet()
    poles = packet["uv_poles"]["gauge_sector"]["poles"]
    assert Fraction(poles["dg1"]["pole_coefficient"]) == Fraction(41, 12)
    assert Fraction(poles["dg2"]["pole_coefficient"]) == Fraction(-19, 12)
    assert Fraction(poles["dg3"]["pole_coefficient"]) == Fraction(-7, 2)
    assert "dg1" not in packet["uv_poles"]["open_pole_values"]["generators"]
    assert "dlam" in packet["uv_poles"]["open_pole_values"]["generators"]


def test_no_gauge_parameter_generator() -> None:
    packet = load_packet()
    assert not any(g.startswith("dxi") for g in packet["generators"])


def test_w_mass_ct_terms() -> None:
    table = json.loads(producer.TABLE_PATH.read_text(encoding="utf-8"))
    entry = next(e for e in table["entries"]
                 if e["fields"] == ["Wm", "Wp"] and e["structure"] == "vector_bilinear_mass")
    terms = producer.ct_terms_for_record(entry)
    by_gen: dict[str, Fraction] = {}
    for term in terms:
        by_gen[term["generator"]] = by_gen.get(term["generator"], Fraction(0)) + Fraction(term["prefactor"])
    assert by_gen["dg2"] == Fraction(1, 2)
    assert by_gen["dv"] == Fraction(1, 2)
    assert by_gen["dZW"] == Fraction(1, 4)
    assert set(by_gen) == {"dg2", "dv", "dZW"}


def test_yukawa_ct_terms_carry_contraction_indices() -> None:
    table = json.loads(producer.TABLE_PATH.read_text(encoding="utf-8"))
    entry = next(e for e in table["entries"]
                 if e["fields"] == ["Gm", "dL_bar", "uR"] and e["structure"] == "fermion_scalar_yukawa")
    terms = producer.ct_terms_for_record(entry)
    generators = {t["generator"] for t in terms}
    assert "dYu[1][2]" in generators
    assert "dZH" in generators
    zq_terms = [t for t in terms if t["generator"] == "dZQLd[3][1]"]
    assert sorted((t["bar_gen"], t["plain_gen"]) for t in zq_terms) == [(3, 1), (3, 2), (3, 3)]
    zq = next(t for t in zq_terms if t["plain_gen"] == 2)
    assert zq["residual_powers"] == [["Yu[1][2]", 1]]
    assert Fraction(zq["prefactor"]) == Fraction(1, 2)
    zu = next(t for t in terms if t["generator"] == "dZuR[2][3]" and t["bar_gen"] == 1)
    assert zu["plain_gen"] == 3 and zu["residual_powers"] == [["Yu[1][2]", 1]]


def test_unreachable_directions_are_symmetry_forbidden() -> None:
    packet = load_packet()
    unreachable = {(tuple(d["fields"]), d["structure"]) for d in packet["general_local_basis"]["unreachable"]}
    assert (("G0", "h"), "scalar_bilinear_mass") in unreachable
    assert (("G0",), "scalar_tadpole") in unreachable
    assert (("B", "Wm", "Wp"), "yang_mills_three_point") in unreachable
    assert (("Gm", "W3", "Wp"), "scalar_scalar_gauge") not in unreachable
    assert (("B", "Gm", "Wp"), "scalar_gauge_gauge") not in unreachable
    assert (("Gm", "W3", "Wp"), "scalar_gauge_gauge") in unreachable


def test_controls_fire() -> None:
    packet = load_packet()
    for name, control in packet["controls"].items():
        assert control["expected_failure"] and control["failed"], name


def test_checker_passes() -> None:
    verdict = checker.check()
    assert verdict["status"] == "PASS", verdict["problems"]


def test_checker_rejects_mutated_packet(tmp_path, monkeypatch) -> None:
    packet = load_packet()
    packet["ct_matrix"]["rank"] = 153
    mutated = tmp_path / "renormalization_ct_1.json"
    mutated.write_text(json.dumps(packet), encoding="utf-8")
    monkeypatch.setattr(checker, "PACKET_PATH", mutated)
    verdict = checker.check()
    assert verdict["status"] == "FAIL"
    assert any("rank" in p for p in verdict["problems"])
