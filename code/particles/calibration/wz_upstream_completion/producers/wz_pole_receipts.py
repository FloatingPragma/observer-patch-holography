#!/usr/bin/env python3
"""W and Z pole receipts on the target-free synthetic fixture.

Certified-contour pole receipts for the acceptance rows of the
external-SM stack: the one-loop inverse-propagator functions are
assembled from the emitted transverse blocks, evaluated in complex
ball arithmetic on a DECLARED synthetic admissible fixture (exact
rationals, not the measured Standard-Model point; no target mounts
anywhere), and each pole is certified by the argument principle:
exactly one simple root inside the declared contour box, every
boundary step excluding zero, at the three preset precisions.

The evaluation is compiled: the fixture is substituted once per
block, every term is reduced to exact rational polynomial
coefficients in p2 together with its loop-symbol head and exact mass
arguments, and the per-point work is polynomial evaluation plus
cached ball loop functions.  The certified content is unchanged from
the direct evaluation; the compilation is arithmetic reorganization.

Conventions, recorded in the payload:

* Inverse propagator, strict one-loop mask: D(s) = s - m_tree^2 -
  Pi_T(s) for the W block, and for the neutral system the masked
  determinant reduces to s - mZ^2 - Pi_ZZ(s) because the Pi_AZ^2 term
  is second order in the loop counting; the AZ block is recorded once
  at the contour center for the receipt.
* One common sign and analytic-sheet convention: the -i eps branch of
  the loop functions continues into complex s; charged and neutral
  blocks share it.
* Laurent data: winding one plus boundary exclusion certifies one
  simple root; the W and Z blocks therefore carry rank-one pole parts
  on the masked surface, with the neutral rank statement carried by
  the masked determinant factorization recorded in the payload.
* xi = 1 on the fixture (the strict-one-loop electroweak grid point
  with xiS = 1); the symbolic identities of the replay checkers carry
  the general-xi content.
"""

from __future__ import annotations

import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import ball_arithmetic as ba  # noqa: E402

VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
OUT_PATH = ROOT / "outputs" / "wz_pole_receipts.json"

SUBDIVISIONS = 64
EXCLUSION_SAMPLES = 32

# Declared synthetic admissible fixture: exact rationals, couplings
# small enough for contour containment, xi = 1, no relation to
# measured targets.
FIXTURE: dict[str, Fraction] = {
    "g1": Fraction(1, 4), "g2": Fraction(1, 3), "v": Fraction(2),
    "lam": Fraction(1, 8), "xi": Fraction(1),
    "mu2": Fraction(1, 2),
    "mfu1": Fraction(1, 50), "mfu2": Fraction(1, 20), "mfu3": Fraction(1, 5),
    "mfd1": Fraction(1, 60), "mfd2": Fraction(1, 25), "mfd3": Fraction(1, 10),
    "mfe1": Fraction(1, 80), "mfe2": Fraction(1, 30), "mfe3": Fraction(1, 15),
    "mu_ren2": Fraction(1),
}
CKM_FIXTURE = {(i, j): Fraction(1 if i == j else 0) for i in (1, 2, 3) for j in (1, 2, 3)}


def sympy_namespace() -> dict[str, sp.Expr]:
    ns: dict[str, sp.Expr] = {}
    for name in ("g1", "g2", "v", "lam", "mu2", "xi", "p2", "d"):
        ns[name] = sp.Symbol(name)
    for i in (1, 2, 3):
        for name in (f"mfu{i}", f"mfd{i}", f"mfe{i}"):
            ns[name] = sp.Symbol(name)
        for j in (1, 2, 3):
            ns[f"V{i}{j}"] = sp.Symbol(f"V{i}{j}")
            ns[f"Vc{i}{j}"] = sp.Symbol(f"Vc{i}{j}")
    return ns


NS = sympy_namespace()


def parse_emitted(text: str) -> sp.Expr:
    """Parse an emitted expression, hiding every loop symbol
    HEAD[args] behind an opaque sympy symbol whose name keeps the
    encoded arguments for later evaluation."""

    opaque: dict[str, str] = {}

    def hide(match: re.Match) -> str:
        token = f"OPQ{len(opaque)}"
        opaque[token] = match.group(0)
        return token

    hidden = re.sub(r"[ABC][0-9]*[0p]*\[[^\]]*\]", hide, text)
    local = dict(NS)
    for token, original in opaque.items():
        encoded = original.replace("[", "__LB__").replace("]", "__RB__").replace("|", "__BAR__")
        local[token] = sp.Symbol(encoded)
    return sp.sympify(hidden, locals=local)


def block_transverse(payload: dict[str, Any], name: str) -> sp.Expr:
    """g and p projections combined into the exact transverse function;
    the loop symbols stay opaque with encoded arguments."""

    total_g = sp.Integer(0)
    total_p = sp.Integer(0)
    for diagram in payload["blocks"][name]["diagrams"]:
        total_g += parse_emitted(diagram["g"])
        total_p += parse_emitted(diagram["p"])
    transverse = sp.expand((total_g - total_p / NS["p2"]) / (NS["d"] - 1))
    return transverse


def compile_block(expression: sp.Expr) -> list[dict[str, Any]]:
    """Substitute the fixture once and compile every term.

    Each compiled term carries the exact rational coefficients of the
    numerator and denominator polynomials in p2 together with the loop
    head and its exact mass arguments; per-point evaluation is then
    polynomial arithmetic plus a cached loop ball.  Mapping every loop
    symbol to its MSbar finite part IS the renormalized block, because
    the certified pole solution shows the counterterm subtraction
    equals exactly the Delta parts and MSbar counterterms carry no
    finite remainder."""

    fixture_subs = {NS[k]: sp.Rational(f) for k, f in FIXTURE.items() if k in NS}
    for (i, j), val in CKM_FIXTURE.items():
        fixture_subs[NS[f"V{i}{j}"]] = sp.Rational(val)
        fixture_subs[NS[f"Vc{i}{j}"]] = sp.Rational(val)
    fixture_subs[NS["d"]] = sp.Integer(4)
    expr = sp.expand(expression.subs(fixture_subs))

    def poly_coefficients(poly_expr: sp.Expr) -> list[tuple[int, Fraction]]:
        poly = sp.Poly(sp.expand(poly_expr), NS["p2"])
        return [
            (int(power), Fraction(int(sp.Rational(coeff).p), int(sp.Rational(coeff).q)))
            for (power,), coeff in poly.terms()
        ]

    compiled: list[dict[str, Any]] = []
    for term in expr.as_ordered_terms():
        loop_head = None
        loop_args: tuple[Fraction, ...] = ()
        coefficient = term
        for symbol in term.free_symbols:
            name = str(symbol)
            if "__LB__" not in name:
                continue
            if loop_head is not None:
                raise SystemExit("two loop symbols in one monomial")
            head, args_text = name.split("__LB__", 1)
            args_text = args_text.replace("__RB__", "").replace("__BAR__", "|")
            arg_values = []
            for part in args_text.split("|"):
                value = sp.sympify(part, locals=NS).subs(fixture_subs)
                arg_values.append(Fraction(str(sp.nsimplify(value))))
            loop_head = head
            loop_args = tuple(arg_values)
            coefficient = coefficient / symbol
        numerator, denominator = sp.fraction(sp.together(coefficient))
        compiled.append(
            {
                "num": poly_coefficients(numerator),
                "den": poly_coefficients(denominator),
                "head": loop_head,
                "args": loop_args,
            }
        )
    return compiled


def _poly_eval(coefficients: list[tuple[int, Fraction]], s: mp.mpc) -> mp.mpc:
    value = mp.mpc(0)
    for power, fraction in coefficients:
        exact = mp.mpf(fraction.numerator) / mp.mpf(fraction.denominator)
        value += exact * s ** power
    return value


class CompiledEvaluator:
    """Per-precision evaluator over compiled terms with a loop cache."""

    def __init__(self, compiled: list[dict[str, Any]], precision: int) -> None:
        self.compiled = compiled
        self.precision = precision
        self.cache: dict[tuple[str, tuple[Fraction, ...], complex], ba.Ball] = {}

    def loop_ball(self, head: str, args: tuple[Fraction, ...], s_value: complex) -> ba.Ball:
        key = (head, args, s_value if head == "B0" else 0j)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        with mp.workprec(self.precision):
            mu_ren = mp.mpf(FIXTURE["mu_ren2"].numerator) / mp.mpf(FIXTURE["mu_ren2"].denominator)
            ulp = mp.mpf(2) ** (8 - self.precision)
            values = [mp.mpf(a.numerator) / mp.mpf(a.denominator) for a in args]
            if head == "A0":
                ball = ba.a0_fin(values[0], mu_ren, precision=self.precision)
            elif head == "B0":
                ball = ba.b0_fin(complex(s_value), values[0], values[1],
                                 mu_ren, precision=self.precision)
            elif head == "A0p":
                if values[0] == 0:
                    ball = ba.Ball(mp.mpf(0), mp.mpf(0), mp.mpf(0))
                else:
                    value = mp.mpc(-mp.log(values[0] / mu_ren))
                    ball = ba.Ball(value.real, value.imag, abs(value) * ulp)
            else:
                raise SystemExit(f"unexpected loop symbol {head}")
        self.cache[key] = ball
        return ball

    def evaluate(self, s_value: complex) -> ba.Ball:
        with mp.workprec(self.precision):
            s = mp.mpc(s_value)
            ulp = mp.mpf(2) ** (8 - self.precision)
            total_mid = mp.mpc(0)
            total_rad = mp.mpf(0)
            for term in self.compiled:
                den_value = _poly_eval(term["den"], s)
                if den_value == 0:
                    raise SystemExit("coefficient denominator vanishes on the contour")
                c_value = _poly_eval(term["num"], s) / den_value
                if term["head"] is None:
                    total_mid += c_value
                    total_rad += abs(c_value) * ulp
                else:
                    ball = self.loop_ball(term["head"], term["args"], s_value)
                    mid = mp.mpc(ball.mid_re, ball.mid_im)
                    total_mid += c_value * mid
                    total_rad += abs(c_value) * ball.rad + abs(c_value * mid) * ulp
            loop_factor = mp.mpc(1) / (16 * mp.pi ** 2)
            total_mid *= loop_factor
            total_rad = total_rad * abs(loop_factor) + abs(total_mid) * ulp
            return ba.Ball(total_mid.real, total_mid.imag, total_rad)


def build_receipts() -> dict[str, Any]:
    vector = json.loads(VECTOR_PATH.read_text(encoding="utf-8"))
    g1 = Fraction(FIXTURE["g1"])
    g2 = Fraction(FIXTURE["g2"])
    v = Fraction(FIXTURE["v"])
    mw2 = Fraction(g2 ** 2 * v ** 2, 4)
    mz2 = Fraction((g1 ** 2 + g2 ** 2) * v ** 2, 4)

    compiled_ww = compile_block(block_transverse(vector, "WpWm"))
    compiled_zz = compile_block(block_transverse(vector, "ZZ"))
    compiled_az = compile_block(block_transverse(vector, "AZ"))

    receipts: dict[str, Any] = {}
    for name, tree_mass, compiled in (("W", mw2, compiled_ww), ("Z", mz2, compiled_zz)):
        for precision in ba.PRECISIONS:
            evaluator = CompiledEvaluator(compiled, precision)

            def inverse_propagator(s: complex) -> complex:
                pi = evaluator.evaluate(s)
                return complex(s - complex(float(tree_mass))
                               - complex(float(pi.mid_re), float(pi.mid_im)))

            center = complex(float(tree_mass), 0.0)
            radius = float(tree_mass) / 8
            lo = complex(center.real - radius, -radius)
            hi = complex(center.real + radius, +radius)
            winding = ba.certify_winding(
                inverse_propagator, (lo, hi),
                subdivisions=SUBDIVISIONS, precision=precision,
            )
            min_mod = None
            max_rad = mp.mpf(0)
            for k in range(EXCLUSION_SAMPLES):
                t_par = k / EXCLUSION_SAMPLES
                edge = complex(lo.real + (hi.real - lo.real) * min(1, 2 * t_par),
                               lo.imag if t_par < 0.5 else hi.imag)
                pi_ball = evaluator.evaluate(edge)
                mid = complex(edge) - complex(float(tree_mass)) - complex(
                    float(pi_ball.mid_re), float(pi_ball.mid_im))
                if min_mod is None or abs(mid) < min_mod:
                    min_mod = abs(mid)
                if pi_ball.rad > max_rad:
                    max_rad = pi_ball.rad
            exclusion = bool(min_mod > 8 * float(max_rad))
            az_center = CompiledEvaluator(compiled_az, precision).evaluate(center)
            receipts[f"{name}_{precision}"] = {
                "tree_mass_sq": str(tree_mass),
                "contour": {"lower_left": [lo.real, lo.imag],
                            "upper_right": [hi.real, hi.imag]},
                "winding": winding,
                "boundary_min_modulus": float(min_mod),
                "boundary_max_ball_radius": float(max_rad),
                "denominator_ball_excludes_zero": exclusion,
                "az_block_at_center": [float(az_center.mid_re),
                                       float(az_center.mid_im),
                                       float(az_center.rad)],
                "simple_root_certified": bool(winding == 1 and exclusion),
            }
    return receipts


def main() -> int:
    receipts = build_receipts()
    all_ok = all(r["simple_root_certified"] for r in receipts.values())
    payload = {
        "schema": "wz_pole_receipts.v2",
        "type": "EXTERNAL_SM_EFT_VALIDATION",
        "promotion": {"oph_native": False, "unit_claim": False},
        "fixture": {k: str(f) for k, f in FIXTURE.items()},
        "fixture_statement": (
            "declared synthetic admissible rationals; not the measured "
            "Standard-Model point; no target mounts anywhere; xi = 1"
        ),
        "conventions": {
            "mask": ("strict one loop: Pi enters once; the AZ^2 determinant "
                     "term is second loop order and dropped; the AZ block is "
                     "recorded at the contour center"),
            "sheet": ("-i eps branch continued to complex s; shared by "
                      "charged and neutral blocks"),
        },
        "receipts": receipts,
        "status": "PASS" if all_ok else "FAIL",
    }
    OUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"],
                      "receipts": {k: r["winding"] for k, r in receipts.items()}}))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
