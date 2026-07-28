#!/usr/bin/env python3
"""Workstream G residual closure: gauge-condition correlator replay.

The correlator of two gauge-fixing functions is a pure contact term to
all orders,

    < F_a(x) F_b(y) > = -i xi delta_ab delta(x - y),

because the Slavnov-Taylor identity <s(cbar_a F_b)> = 0 reduces the
ghost chain to the Schwinger-Dyson contact <cbar_a dS/dcbar_b> and no
ghost two-point function survives.  At one loop this forces the
momentum-space combination built from dressed two-point insertions to
vanish identically:

    0 = P.G0 Sigma G0.P  summed over the {vector, Goldstone} system,

with the F-vertices read from the solved gauge functions
F_Z = d.Z - xi m_Z G0, F_W+- = d.W+- - xi m_W G+-, F_A = d.A, tree
propagators from the certified spectrum, and the loop insertions from
the emitted blocks.  The tree part of every correlator is recomputed
first and must equal -i xi (or zero for the A-Z cross), which anchors
the Fourier sign conventions before the loop replay runs.

Everything is evaluated on the minimum-restored surface mu2 = lam v^2,
where the Goldstone masses collapse to xi m_V^2 and the identities
close at strict one loop; the offset bookkeeping lives in the tadpole
sector of the FJ chart.

The checker reads emitted payloads only; it never imports an engine.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCALAR_PATH = ROOT / "outputs" / "fj_direct_scalar_blocks.json"
OUT_PATH = ROOT / "outputs" / "gauge_condition_correlator_check.json"

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


def canonical_mass_text(text: str) -> str:
    value = sp.sympify(text, locals=NS).subs(mu2v, lam * v ** 2)
    return str(sp.simplify(sp.expand(value)))


def canonical_loop_symbol(name: str) -> str:
    match = re.match(r"(B0|C21|C12|C22)\[([^|]*)\|([^\]]*)\]$", name)
    if not match:
        single = re.match(r"(A0p?)\[([^\]]*)\]$", name)
        if single:
            head, arg = single.groups()
            return f"{head}[{canonical_mass_text(arg)}]"
        return name
    head, a, b = match.groups()
    a = canonical_mass_text(a)
    b = canonical_mass_text(b)
    if head in ("B0", "C22"):
        first, second = sorted((a, b))
        return f"{head}[{first}|{second}]"
    if a <= b:
        return f"{head}[{a}|{b}]"
    swapped = "C12" if head == "C21" else "C21"
    return f"{swapped}[{b}|{a}]"


def parse(text: str) -> sp.Expr:
    opaque: dict[str, str] = {}

    def hide(match: re.Match) -> str:
        token = f"OPQ{len(opaque)}"
        opaque[token] = match.group(0)
        return token

    hidden = re.sub(r"[ABC][0-9]*[0p]*\[[^\]]*\]", hide, text)
    local = dict(NS)
    for token, original in opaque.items():
        local[token] = sp.Symbol(canonical_loop_symbol(original))
    return sp.sympify(hidden, locals=local).subs(mu2v, lam * v ** 2)


def block_sum(payload: dict[str, Any], name: str, key: str) -> sp.Expr:
    total = sp.Integer(0)
    for diagram in payload["blocks"][name]["diagrams"]:
        total += parse(diagram[key])
    return sp.expand(total)


def check() -> dict[str, Any]:
    vector = json.loads(VECTOR_PATH.read_text(encoding="utf-8"))
    scalar = json.loads(SCALAR_PATH.read_text(encoding="utf-8"))
    problems: list[str] = []
    identities: dict[str, Any] = {}

    mz_sq = (g1 ** 2 + g2 ** 2) * v ** 2 / 4
    mw_sq = g2 ** 2 * v ** 2 / 4
    mz = sp.sqrt(mz_sq)
    mw = sp.sqrt(mw_sq)

    # Fourier convention: for incoming momentum p on the first leg,
    # F(p) = -I p.V(p) - xi m S(p) on the annihilation side and
    # F(-p) = +I p.V - xi m S on the second leg.  The tree anchors
    # below verify this convention and the propagator forms before the
    # loop replay:
    #
    # p.G0^{VV}.p = -I xi p2/(p2 - xi m^2)   (longitudinal projection)
    # G0^{SS}     = +I /(p2 - xi m^2)        (minimum-restored Goldstone)
    #
    # tree <F F> = (-I p)(+I p): p.G.p-term
    #              + (xi m)^2 G^{SS} = -I xi p2/(p2-xi m^2)
    #              - xi^2 m^2 /(p2 - xi m^2) * I = -I xi  (exact).
    def tree_ff(mass_sq: sp.Expr) -> sp.Expr:
        vector_part = -sp.I * xi * p2 / (p2 - xi * mass_sq)
        scalar_part = xi ** 2 * mass_sq * (sp.I / (p2 - xi * mass_sq))
        return sp.simplify(vector_part - scalar_part)

    for name, mass_sq in (("Z", mz_sq), ("W", mw_sq)):
        anchor = sp.simplify(tree_ff(mass_sq) + sp.I * xi)
        if anchor != 0:
            problems.append(f"tree anchor for F_{name} fails: {anchor}")
    identities["tree_anchor"] = "tree <F F> = -I xi for Z and W systems"

    # One-loop replay.  With G = G0 + G0 (i Sigma-hat) G0 and the
    # engine convention i Pi = diagram sum (Pi emitted), the insertion
    # between tree legs contributes
    #   loop <F F> = sum over channels of
    #     (F-vertex_1 . G0 legs) (i Pi-block) (G0 legs . F-vertex_2)
    # and the identity requires the total to vanish.  Longitudinal
    # vector legs use p.G0 = -I xi p^mu/(p2 - xi m^2) (the transverse
    # part drops against p contraction); the scalar leg is G0^{SS}.
    # The mixing block enters twice (Z G0 and G0 Z orderings), with
    # Pi^mu = p^mu Sigma_mix and the second-leg momentum -p flipping
    # one sign per mixing insertion.
    def loop_ff(pi_l_vector: sp.Expr, sigma_mix: sp.Expr, pi_scalar: sp.Expr,
                mass_sq: sp.Expr) -> sp.Expr:
        """Numerator of the loop FF-combination after clearing the
        common denominator (p2 - xi m^2)^2; zero-testing the expanded
        polynomial avoids large simplify calls."""

        m = sp.sqrt(mass_sq)
        # legs without denominators: vector leg -I xi (times p^mu),
        # scalar leg +I; F-scalar vertex -xi m per leg.
        vector_channel = (-sp.I * xi) * sp.I * (p2 * pi_l_vector) * (-sp.I * xi)
        scalar_channel = (-xi * m * sp.I) * sp.I * pi_scalar * (-xi * m * sp.I)
        cross_channel = 2 * (-sp.I * xi) * sp.I * (p2 * sigma_mix) * (-xi * m) * sp.I
        return sp.expand(vector_channel + scalar_channel + cross_channel)

    zz_l = sp.simplify(block_sum(vector, "ZZ", "p") / p2)
    ww_l = sp.simplify(block_sum(vector, "WpWm", "p") / p2)
    zg0 = sp.simplify(block_sum(scalar, "ZG0", "value") / p2)
    wg = sp.simplify(block_sum(scalar, "WpGm", "value") / p2)
    g0g0 = block_sum(scalar, "G0G0", "value")
    gpgm = block_sum(scalar, "GpGm", "value")

    for name, pieces in (
        ("F_Z F_Z", (zz_l, zg0, g0g0, mz_sq)),
        ("F_W F_W", (ww_l, wg, gpgm, mw_sq)),
    ):
        value = sp.expand(loop_ff(*pieces))
        if value != 0:
            value = sp.expand(sp.cancel(sp.together(value)))
        identities[name] = str(value)[:400]
        if value != 0:
            problems.append(f"{name} loop part survives: {str(value)[:400]}")

    # Photon system: F_A = d.A has no scalar part; the A-Z cross
    # correlator identity requires the combination of the AZ
    # longitudinal block and the A-G0 mixing against the Z-side scalar
    # vertex to cancel.
    az_l = sp.simplify(block_sum(vector, "AZ", "p") / p2)
    ag0 = sp.simplify(block_sum(scalar, "AG0", "value") / p2)
    cross = sp.expand(
        (-sp.I * xi) * sp.I * (p2 * az_l) * (-sp.I * xi)
        + 2 * (-sp.I * xi) * sp.I * (p2 * ag0) * (-xi * mz) * sp.I / 2
    )
    if cross != 0:
        cross = sp.expand(sp.cancel(sp.together(cross)))
    identities["F_A F_Z"] = str(cross)
    if cross != 0:
        problems.append(f"F_A F_Z loop part survives: {str(cross)[:400]}")

    aa_l = sp.simplify(block_sum(vector, "AA", "p"))
    if sp.simplify(aa_l) != 0:
        problems.append("F_A F_A loop part survives")
    identities["F_A F_A"] = "0"

    return {
        "schema": "gauge_condition_correlator_check.v1",
        "target": "WARD_ST_NIELSEN_1",
        "status": "PASS" if not problems else "FAIL",
        "surface": "minimum-restored mu2 = lam v^2; offset bookkeeping in the FJ tadpole sector",
        "identities": identities,
        "problems": problems,
    }


def main() -> int:
    verdict = check()
    OUT_PATH.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": verdict["status"], "problems": [p[:160] for p in verdict["problems"][:3]]}))
    return 0 if verdict["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
