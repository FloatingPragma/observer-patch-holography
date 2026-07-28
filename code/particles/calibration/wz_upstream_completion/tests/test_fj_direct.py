"""Regression tests for the FJ direct-engine stage: reduction core,
mass-basis specialization, and the emitted vector-block payload."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import loop_reduction as lr  # noqa: E402
import fj_spectrum  # noqa: E402

PAYLOAD = json.loads((ROOT / "outputs" / "fj_direct_vector_blocks.json").read_text(encoding="utf-8"))


def test_reduction_ward_identity_fermion() -> None:
    msq = sp.Symbol("msq")
    kkp = lr.k2 - lr.kp
    p_proj = sp.expand(4 * (2 * lr.kp * (lr.kp - lr.p2) - lr.p2 * (kkp - msq)))
    reduced = lr.reduce_two_point(p_proj, 1, 1).subs({lr.m1sq: msq, lr.m2sq: msq})
    reduced = reduced.subs({sp.Symbol("A0m2"): sp.Symbol("A0m1"), sp.Symbol("A0pm2"): sp.Symbol("A0pm1")})
    assert sp.simplify(reduced) == 0


def test_reduction_scalar_qed_pole() -> None:
    msq = sp.Symbol("msq")
    g_sb = sp.expand(4 * lr.k2 - 4 * lr.kp + lr.p2)
    p_sb = sp.expand((2 * lr.kp - lr.p2) ** 2)
    ident = {sp.Symbol("A0m2"): sp.Symbol("A0m1"), sp.Symbol("A0pm2"): sp.Symbol("A0pm1")}
    gred = sp.expand(lr.reduce_two_point(g_sb, 1, 1).subs({lr.m1sq: msq, lr.m2sq: msq}).subs(ident))
    pred = sp.expand(lr.reduce_two_point(p_sb, 1, 1).subs({lr.m1sq: msq, lr.m2sq: msq}).subs(ident))
    gtot = sp.expand(gred - 2 * lr.d_sym * sp.Symbol("A0m1"))
    ptot = sp.expand(pred - 2 * lr.p2 * sp.Symbol("A0m1"))
    assert sp.simplify(ptot) == 0
    pi_t, _ = lr.transverse_longitudinal(gtot, ptot)
    pole = lr.uv_pole(sp.expand(pi_t), {lr.m1sq: msq})
    assert sp.simplify(pole + lr.p2 / 3) == 0


def test_spectrum_checks_pass() -> None:
    records = fj_spectrum.specialized_records()
    spectrum, checks = fj_spectrum.spectrum_and_checks(records)
    assert checks == []
    assert sp.simplify(spectrum["Z"] - (fj_spectrum.g1 ** 2 + fj_spectrum.g2 ** 2) * fj_spectrum.v ** 2 / 4) == 0
    assert spectrum["A"] == 0 and spectrum["cA"] == 0


def test_payload_controls_all_passed() -> None:
    for name, control in PAYLOAD["controls"].items():
        assert control.get("passed", True), name


def test_photon_block_transversality_and_sectors() -> None:
    aa = PAYLOAD["blocks"]["AA"]
    assert sp.simplify(sp.sympify(aa["longitudinal_pole"])) == 0
    control = PAYLOAD["controls"]["photon_fermion_sector"]
    assert sp.simplify(sp.sympify(control["engine"]) - sp.sympify(control["expected"])) == 0


def test_blocks_present_with_diagram_counts() -> None:
    for name, minimum in (("AA", 20), ("AZ", 20), ("ZZ", 30), ("WpWm", 40)):
        block = PAYLOAD["blocks"][name]
        assert block["diagram_count"] >= minimum, (name, block["diagram_count"])


def test_tadpole_insertions_present() -> None:
    for name in ("ZZ", "WpWm"):
        kinds = {d["kind"] for d in PAYLOAD["blocks"][name]["diagrams"]}
        assert "tadpole_insertion" in kinds, name
    tree = PAYLOAD["one_point_h"]["contributions"]["tree"]
    assert "mu2" in tree and "lam" in tree


def test_ckm_structure_in_ww_block() -> None:
    ww_pole = PAYLOAD["blocks"]["WpWm"]["transverse_pole"]
    assert "V11" in ww_pole and "Vc32" in ww_pole or "Vc23" in ww_pole
