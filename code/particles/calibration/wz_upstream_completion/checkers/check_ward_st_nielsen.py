#!/usr/bin/env python3
"""Workstream G: Ward, Slavnov-Taylor and Nielsen replay checker.

The checker reads emitted coefficient records and block payloads; it
never imports a loop engine.  Every quantity it certifies is
recomputed here from the payload strings with its own parser and
sympy algebra:

* per-diagram re-summation: the stored block poles equal the sums of
  the stored per-diagram poles;
* photon Ward replay: the COMPLETE longitudinal AA expression (all
  instantiated loop functions, not only the pole) sums to zero on the
  minimum-restored surface mu2 = lam v^2; the offset is itself of
  tadpole order, so the identity closes there at strict one loop
  while the FJ chart keeps the minimum unimposed in the records;
* charge universality and Nielsen replay: the combination
  dZ_e = -(1/2) dPi_T^AA/dp2 - (g1/(2 g2)) (2/m_Z^2) Pi_T^AZ(0) is
  xi-independent (the Nielsen statement for the charge) and equals
  (b1 + b2)/2 e^2 from the census betas;
* custodial replay: the G0G0 and G+G- poles coincide at g1 -> 0,
  unit mixing, degenerate quark masses;
* FJ-equivalence replay: the subtraction insertion minus v-shift is
  recomputed here from the receipt fields for every block;
* counterterm reachability: every block pole is polynomial in p2 of
  degree at most one (renormalizable structure), and no pole
  component lies on a structurally unreachable direction of the
  counterterm packet (the engine emits no G0-h or CP-odd block).

The remaining full gauge-condition correlator protection identity
(the G_Z G_Z combination with the ghost two-point functions) is
recorded as the declared residual item of this target.
"""

from __future__ import annotations

import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCALAR_PATH = ROOT / "outputs" / "fj_direct_scalar_blocks.json"
CONVERTED_PATH = ROOT / "outputs" / "fj_converted_equivalence.json"
CT_PATH = ROOT / "outputs" / "renormalization_ct_1.json"
MATCHING_PATH = ROOT / "outputs" / "eft_matching_1.json"
OUT_PATH = ROOT / "outputs" / "ward_st_nielsen_check.json"

g1, g2, v = sp.symbols("g1 g2 v", positive=True)
lam, mu2v, xi, p2, d = sp.symbols("lam mu2 xi p2 d", positive=True)


def namespace() -> dict[str, sp.Expr]:
    ns = {"g1": g1, "g2": g2, "v": v, "lam": lam, "mu2": mu2v, "xi": xi, "p2": p2, "d": d}
    for i in (1, 2, 3):
        for name in (f"mfu{i}", f"mfd{i}", f"mfe{i}"):
            ns[name] = sp.Symbol(name)
        for j in (1, 2, 3):
            ns[f"V{i}{j}"] = sp.Symbol(f"V{i}{j}")
            ns[f"Vc{i}{j}"] = sp.Symbol(f"Vc{i}{j}")
    return ns


NS = namespace()


def parse(text: str, at_minimum: bool = False) -> sp.Expr:
    """Parse an emitted value; at_minimum evaluates on the
    minimum-restored surface mu2 -> lam v^2, including inside the
    loop-function arguments.  The FJ chart keeps the minimum
    unimposed, and the (lam v^2 - mu2) offset is itself of tadpole
    order, so symmetry identities close on this surface at strict one
    loop while the offset bookkeeping lives in the tadpole sector."""

    opaque: dict[str, str] = {}

    def hide(match: re.Match) -> str:
        token = f"OPQ{len(opaque)}"
        opaque[token] = match.group(0)
        return token

    hidden = re.sub(r"[ABC][0-9]*[0p]*\[[^\]]*\]", hide, text)
    local = dict(NS)
    for token, original in opaque.items():
        local[token] = sp.Symbol(canonical_loop_symbol(original, at_minimum))
    expr = sp.sympify(hidden, locals=local)
    if at_minimum:
        expr = expr.subs(mu2v, lam * v ** 2)
    return expr


def canonical_mass_text(text: str, at_minimum: bool) -> str:
    value = sp.sympify(text, locals=NS)
    if at_minimum:
        value = value.subs(mu2v, lam * v ** 2)
    return str(sp.simplify(sp.expand(value)))


def canonical_loop_symbol(name: str, at_minimum: bool = False) -> str:
    """Exact argument symmetries of the loop basis: B0 and C22 are
    symmetric under mass exchange, and I(2,1) with masses (a, b) is
    I(1,2) with (b, a); the canonical form normalizes the argument
    text, applies the evaluation surface, sorts symmetric pairs, and
    renames C12 to C21 on a swap."""

    match = re.match(r"(B0|C21|C12|C22)\[([^|]*)\|([^\]]*)\]$", name)
    if not match:
        single = re.match(r"(A0p?)\[([^\]]*)\]$", name)
        if single:
            head, arg = single.groups()
            return f"{head}[{canonical_mass_text(arg, at_minimum)}]"
        return name
    head, a, b = match.groups()
    a = canonical_mass_text(a, at_minimum)
    b = canonical_mass_text(b, at_minimum)
    if head in ("B0", "C22"):
        first, second = sorted((a, b))
        return f"{head}[{first}|{second}]"
    if a <= b:
        return f"{head}[{a}|{b}]"
    swapped = "C12" if head == "C21" else "C21"
    return f"{swapped}[{b}|{a}]"


def check() -> dict[str, Any]:
    vector = json.loads(VECTOR_PATH.read_text(encoding="utf-8"))
    scalar = json.loads(SCALAR_PATH.read_text(encoding="utf-8"))
    converted = json.loads(CONVERTED_PATH.read_text(encoding="utf-8"))
    ct = json.loads(CT_PATH.read_text(encoding="utf-8"))
    matching = json.loads(MATCHING_PATH.read_text(encoding="utf-8"))
    problems: list[str] = []
    replays: dict[str, Any] = {}

    # 1. per-diagram pole re-summation for every block.
    for source, projections in ((vector, ("g_pole", "p_pole")), (scalar, ("pole",))):
        for name, block in source["blocks"].items():
            if "g_pole" in block:
                g_sum = sum(parse(diagram["g_pole"]) for diagram in block["diagrams"])
                p_sum = sum(parse(diagram["p_pole"]) for diagram in block["diagrams"])
                if sp.simplify(g_sum - parse(block["g_pole"])) != 0:
                    problems.append(f"{name}: stored g-pole differs from the diagram sum")
                if sp.simplify(p_sum - parse(block["p_pole"])) != 0:
                    problems.append(f"{name}: stored p-pole differs from the diagram sum")
            else:
                total = sum(parse(diagram["pole"]) for diagram in block["diagrams"])
                divide = p2 if block["normalization"].startswith("Sigma") else sp.Integer(1)
                if sp.simplify(total / divide - parse(block["pole"])) != 0:
                    problems.append(f"{name}: stored pole differs from the diagram sum")
    replays["per_diagram_resummation"] = "all blocks"

    # 2. complete photon longitudinal cancellation, loop functions and
    # finite parts included.
    aa_longitudinal = sp.Integer(0)
    for diagram in vector["blocks"]["AA"]["diagrams"]:
        aa_longitudinal += parse(diagram["p"], at_minimum=True)
    aa_longitudinal = sp.simplify(sp.expand(aa_longitudinal))
    if aa_longitudinal != 0:
        problems.append(f"photon longitudinal expression survives: {aa_longitudinal}")
    replays["photon_ward_full_expression"] = str(aa_longitudinal)

    # 3. charge universality and the Nielsen statement for the charge.
    e_sq = g1 ** 2 * g2 ** 2 / (g1 ** 2 + g2 ** 2)
    mz_sq = (g1 ** 2 + g2 ** 2) * v ** 2 / 4
    aa_slope = sp.expand(parse(vector["blocks"]["AA"]["transverse_pole"])).coeff(p2)
    az_zero = sp.simplify(parse(vector["blocks"]["AZ"]["transverse_pole"]).subs(p2, 0))
    dz_e = sp.simplify(-sp.Rational(1, 2) * aa_slope - (g1 / (2 * g2)) * 2 * az_zero / mz_sq)
    b1 = sp.Rational(Fraction(matching["gauge_betas"]["coefficients"]["b1"]))
    b2 = sp.Rational(Fraction(matching["gauge_betas"]["coefficients"]["b2"]))
    if dz_e.has(xi):
        problems.append("charge combination is xi-dependent: Nielsen replay fails")
    if sp.simplify(dz_e - (b1 + b2) / 2 * e_sq) != 0:
        problems.append(f"charge combination differs from census (b1+b2)/2: {dz_e}")
    replays["charge_universality_nielsen"] = str(dz_e)

    # 4. custodial replay.
    custodial: dict[sp.Symbol, sp.Expr] = {g1: sp.Integer(0)}
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            custodial[sp.Symbol(f"V{i}{j}")] = sp.Integer(1 if i == j else 0)
            custodial[sp.Symbol(f"Vc{i}{j}")] = sp.Integer(1 if i == j else 0)
        custodial[sp.Symbol(f"mfd{i}")] = sp.Symbol(f"mfu{i}")
    g0_pole = parse(scalar["blocks"]["G0G0"]["pole"]).subs(custodial)
    gp_pole = parse(scalar["blocks"]["GpGm"]["pole"]).subs(custodial)
    if sp.simplify(sp.expand(g0_pole - gp_pole)) != 0:
        problems.append("custodial replay fails")
    replays["custodial"] = "G0G0 equals GpGm in the custodial limit"

    # 5. FJ-equivalence replay: the subtraction is recomputed from the
    # receipt's own insertion and v-shift fields, never trusted from
    # its difference field.
    for name, receipt in converted["equivalence_receipts"].items():
        recomputed = sp.expand(parse(receipt["insertion"]) - parse(receipt["v_shift_of_tree_record"]))
        if sp.simplify(recomputed) != 0:
            problems.append(f"FJ equivalence residual in {name}")
    replays["fj_equivalence"] = sorted(converted["equivalence_receipts"])

    # 6. counterterm reachability of the poles.
    unreachable = {(tuple(direction["fields"]), direction["structure"])
                   for direction in ct["general_local_basis"]["unreachable"]}
    if (("G0", "h"), "scalar_bilinear_mass") not in unreachable:
        problems.append("counterterm packet lost the G0-h unreachable direction")
    emitted_blocks = set(vector["blocks"]) | set(scalar["blocks"])
    if "G0h" in emitted_blocks:
        problems.append("an engine emitted a CP-odd G0-h block")
    for source in (vector, scalar):
        for name, block in source["blocks"].items():
            pole_text = block.get("transverse_pole", block.get("pole"))
            pole = sp.expand(parse(pole_text))
            if sp.degree(sp.Poly(pole, p2)) > 1:
                problems.append(f"{name}: pole is not linear in p2")
    replays["ct_reachability"] = "all block poles linear in p2; no unreachable-direction block emitted"

    verdict = {
        "schema": "ward_st_nielsen_check.v1",
        "target": "WARD_ST_NIELSEN_1",
        "status": "PASS" if not problems else "FAIL",
        "replays": replays,
        "problems": problems,
        "residual": {
            "gauge_condition_correlator_protection": (
                "the G_Z G_Z combination over {ZZ longitudinal, Z-G0 "
                "mixing, G0G0, ghost two-point functions} is the "
                "declared residual projection of this target; its "
                "inputs (A-G0 and Z-G0 mixings) are emitted, the ghost "
                "two-point blocks are the missing input"
            ),
        },
    }
    return verdict


def main() -> int:
    verdict = check()
    OUT_PATH.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": verdict["status"], "problems": verdict["problems"][:4]}))
    return 0 if verdict["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
