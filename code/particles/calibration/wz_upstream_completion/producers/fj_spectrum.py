#!/usr/bin/env python3
"""Mass-basis specialization of the certified rule table for the
direct FJ engine.

Everything here is a mechanical transformation of the dual-derived
rule table; nothing is entered by hand:

* The evaluation subfamily xi1 = xi2 = xi is declared and recorded
  (the independent-xi stress grid stays with the converted-engine
  stage).
* The neutral vectors and neutral ghosts rotate to the mass basis
  with A = (g1 W3 + g2 B)/gz and Z = (g2 W3 - g1 B)/gz, gz =
  sqrt(g1^2 + g2^2); the rotation is applied by record re-expansion
  and the vanishing of the A mass, the A-Z mass mixing, and the
  massless photon-ghost row is CHECKED, not assumed.
* The propagator spectrum is read from the rotated bilinear records
  with the vacuum minimum NOT imposed: scalar masses keep the exact
  (lam v^2 - mu2) offset of the FJ chart.
* The lepton Yukawa chart is the exact diagonal weak-basis choice
  Ye[i][j] -> delta_ij sqrt2 me_i / v; the quark chart introduces the
  mass-basis symbols mu_i, md_i and the mixing matrix V[i][j] through
  the declared unitary field redefinitions, with the neutral-current
  diagonality following from unitarity."""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

TABLE_PATH = ROOT / "outputs" / "rule_table_engine_a.json"

g1, g2, g3, v, lam, mu2, xi = sp.symbols("g1 g2 g3 v lam mu2 xi", positive=True)
GZ = sp.sqrt(g1 ** 2 + g2 ** 2)
CW = g2 / GZ
SW = g1 / GZ

SYMBOLS: dict[str, sp.Expr] = {
    "g1": g1, "g2": g2, "g3": g3, "v": v, "lam": lam, "mu2": mu2,
    "xi1": xi, "xi2": xi, "sqrt2": sp.sqrt(2), "I": sp.I,
}

# Mass-basis fermion chart symbols.
MU = {i: sp.Symbol(f"mfu{i}", positive=True) for i in (1, 2, 3)}
MD = {i: sp.Symbol(f"mfd{i}", positive=True) for i in (1, 2, 3)}
ME = {i: sp.Symbol(f"mfe{i}", positive=True) for i in (1, 2, 3)}
VCKM = {(i, j): sp.Symbol(f"V{i}{j}") for i in (1, 2, 3) for j in (1, 2, 3)}
VCKMC = {(i, j): sp.Symbol(f"Vc{i}{j}") for i in (1, 2, 3) for j in (1, 2, 3)}

NEUTRAL_ROTATION = {
    "W3": (("Z", CW), ("A", SW)),
    "B": (("Z", -SW), ("A", CW)),
    "c3": (("cZ", CW), ("cA", SW)),
    "cB": (("cZ", -SW), ("cA", CW)),
    "c3_bar": (("cZ_bar", CW), ("cA_bar", SW)),
    "cB_bar": (("cZ_bar", -SW), ("cA_bar", CW)),
}


def coefficient_to_sympy(monomials: list[dict[str, Any]]) -> sp.Expr:
    total = sp.Integer(0)
    for monomial in monomials:
        term = sp.Rational(Fraction(monomial["prefactor"]))
        for name, power in monomial["powers"]:
            if name in SYMBOLS:
                term *= SYMBOLS[name] ** power
            else:
                term *= sp.Symbol(name) ** power
        total += term
    return total


def rotate_record(fields: list[str], coefficient: sp.Expr) -> list[tuple[tuple[str, ...], sp.Expr]]:
    """Expand every neutral gauge-basis field into its mass-basis
    combination, returning the resulting (fields, coefficient) list."""

    expansions: list[tuple[tuple[str, ...], sp.Expr]] = [((), coefficient)]
    for field in fields:
        new_expansions = []
        options = NEUTRAL_ROTATION.get(field, ((field, sp.Integer(1)),))
        for prefix, coeff in expansions:
            for label, weight in options:
                new_expansions.append((prefix + (label,), coeff * weight))
        expansions = new_expansions
    collected: dict[tuple[str, ...], sp.Expr] = {}
    for labels, coeff in expansions:
        key = tuple(sorted(labels))
        collected[key] = sp.simplify(collected.get(key, sp.Integer(0)) + coeff)
    return [(labels, coeff) for labels, coeff in collected.items() if coeff != 0]


def specialized_records() -> list[dict[str, Any]]:
    """The certified table with xi1 = xi2 = xi and the neutral sector
    rotated to the mass basis; coefficients are sympy expressions."""

    table = json.loads(TABLE_PATH.read_text(encoding="utf-8"))
    out: list[dict[str, Any]] = []
    merged: dict[tuple, sp.Expr] = {}
    structures: dict[tuple, str] = {}
    for entry in table["entries"]:
        coefficient = coefficient_to_sympy(entry["coefficient"]["monomials"])
        for labels, coeff in rotate_record(entry["fields"], coefficient):
            key = (labels, entry["structure"])
            merged[key] = sp.simplify(merged.get(key, sp.Integer(0)) + coeff)
            structures[key] = entry["structure"]
    for (labels, structure), coefficient in sorted(merged.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        coefficient = sp.simplify(coefficient)
        if coefficient == 0:
            continue
        out.append({"fields": list(labels), "structure": structure, "coefficient": coefficient})
    return out


def spectrum_and_checks(records: list[dict[str, Any]]) -> tuple[dict[str, sp.Expr], list[str]]:
    """Read the propagator mass-squares from the rotated bilinears and
    verify the diagonalization facts."""

    def bilinear(fields: list[str], structure: str) -> sp.Expr:
        for record in records:
            if record["fields"] == sorted(fields) and record["structure"] == structure:
                return record["coefficient"]
        return sp.Integer(0)

    checks: list[str] = []
    mz_sq = sp.simplify((g1 ** 2 + g2 ** 2) * v ** 2 / 4)
    mw_sq = sp.simplify(g2 ** 2 * v ** 2 / 4)

    if sp.simplify(bilinear(["A", "A"], "vector_bilinear_mass")) != 0:
        checks.append("photon mass does not vanish after rotation")
    if sp.simplify(bilinear(["A", "Z"], "vector_bilinear_mass")) != 0:
        checks.append("A-Z mass mixing survives the rotation")
    zz = bilinear(["Z", "Z"], "vector_bilinear_mass")
    if sp.simplify(zz - mz_sq / 2) != 0:
        checks.append(f"ZZ mass record is {zz}, expected mZ^2/2")
    ww = bilinear(["Wm", "Wp"], "vector_bilinear_mass")
    if sp.simplify(ww - mw_sq) != 0:
        checks.append(f"W mass record is {ww}, expected mW^2")

    if sp.simplify(bilinear(["cA", "cA_bar"], "ghost_scalar_mass")) != 0:
        checks.append("photon ghost acquires a mass")
    for pair in (["cA_bar", "cZ"], ["cA", "cZ_bar"]):
        if sp.simplify(bilinear(pair, "ghost_scalar_mass")) != 0:
            checks.append("neutral ghost mass mixing survives the rotation")
    cz = bilinear(["cZ", "cZ_bar"], "ghost_scalar_mass")
    if sp.simplify(cz + xi * mz_sq) != 0:
        checks.append(f"Z-ghost mass record is {cz}, expected -xi mZ^2")

    # Scalar masses from L = c2 S S (real) or c GpGm (complex pair):
    # real scalar: m^2 = -2 c2; complex pair: m^2 = -c.
    m_h_sq = sp.simplify(-2 * bilinear(["h", "h"], "scalar_bilinear_mass"))
    m_g0_sq = sp.simplify(-2 * bilinear(["G0", "G0"], "scalar_bilinear_mass"))
    m_gp_sq = sp.simplify(-bilinear(["Gm", "Gp"], "scalar_bilinear_mass"))
    if sp.simplify(m_h_sq - (3 * lam * v ** 2 - mu2)) != 0:
        checks.append("h mass differs from 3 lam v^2 - mu2")
    if sp.simplify(m_g0_sq - (lam * v ** 2 - mu2 + xi * mz_sq)) != 0:
        checks.append("G0 mass differs from (lam v^2 - mu2) + xi mZ^2")
    if sp.simplify(m_gp_sq - (lam * v ** 2 - mu2 + xi * mw_sq)) != 0:
        checks.append("charged Goldstone mass differs from (lam v^2 - mu2) + xi mW^2")

    spectrum = {
        "Wp": mw_sq, "Wm": mw_sq, "Z": mz_sq, "A": sp.Integer(0), "Gl": sp.Integer(0),
        "h": m_h_sq, "G0": m_g0_sq, "Gp": m_gp_sq, "Gm": m_gp_sq,
        "cp": xi * mw_sq, "cm": xi * mw_sq, "cZ": xi * mz_sq, "cA": sp.Integer(0),
    }
    for i in (1, 2, 3):
        spectrum[f"u{i}"] = MU[i] ** 2
        spectrum[f"d{i}"] = MD[i] ** 2
        spectrum[f"e{i}"] = ME[i] ** 2
        spectrum[f"nu{i}"] = sp.Integer(0)
    return spectrum, checks


def chart_declaration() -> dict[str, str]:
    return {
        "xi_subfamily": "xi1 = xi2 = xi; the independent-xi stress grid belongs to the converted engine stage",
        "neutral_rotation": "A = (g1 W3 + g2 B)/gz, Z = (g2 W3 - g1 B)/gz, gz = sqrt(g1^2+g2^2); ghosts rotate identically",
        "lepton_chart": "Ye[i][j] = delta_ij sqrt2 mfe_i / v (exact weak-basis freedom; neutrinos massless)",
        "quark_chart": (
            "uL = ULu uL', dL = ULd dL', uR = URu uR', dR = URd dR' with "
            "Yu = (sqrt2/v) ULu diag(mu) URu^dag, Yd = (sqrt2/v) ULd diag(md) URd^dag; "
            "V = ULu^dag ULd; neutral currents stay diagonal by unitarity; "
            "the charged current carries V and the charged Goldstone vertices "
            "carry (mfu_i V[i][j], V[i][j] mfd_j) mass-weighted chiral parts"
        ),
        "minimum": "not imposed: scalar masses keep the exact (lam v^2 - mu2) offset per the FJ contract",
    }


def main() -> int:
    records = specialized_records()
    spectrum, checks = spectrum_and_checks(records)
    if checks:
        for line in checks:
            print("CHECK FAILED:", line)
        return 1
    print(json.dumps({
        "status": "OK",
        "records": len(records),
        "spectrum": {k: str(sp.simplify(val)) for k, val in sorted(spectrum.items())},
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
