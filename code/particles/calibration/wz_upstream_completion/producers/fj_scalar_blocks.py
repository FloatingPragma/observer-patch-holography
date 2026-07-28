#!/usr/bin/env python3
"""Workstream E stage 2: Goldstone, Higgs and vector-scalar mixing
blocks of the direct FJ engine.

Extends the validated vector-block machinery to scalar externals
(hh, G0G0, G+G-) and to the single-index mixing blocks (W+G-, Z G0,
A G0), with the same conventions: Pi = -I x weight x reduced, weight
carrying the loop measure i, propagator prefactors, loop signs, color
and symmetry factors; scalar-type vertices carry i x coefficient x
multiplicity inside the terms, derivative vertices carry the real
momentum form.  A scalar block is a single invariant Pi(p2); a mixing
block is emitted as Sigma_mix(p2) with Pi^mu = p^mu Sigma_mix.

Fermion couplings for the scalar sector follow the declared mass-basis
chart: h couples -mf/v (P_L + P_R) per template; G0 couples
+- I mf/v (P_R - P_L) with the up/down sign split of the table
records; G+ couples (sqrt2/v)(mfu_i V[i][j] P_L - V[i][j] mfd_j P_R)
on the u-bar d side, with the conjugate on the d-bar u side.  The
custodial control checks G0G0 = G+G- exactly in the limit g1 -> 0,
V = 1, mfu = mfd, closing the chart against the boson sectors.

The one-loop A-G0 mixing is computed because the Slavnov-Taylor
replay of Workstream G consumes it together with the AZ longitudinal
part.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import fj_spectrum  # noqa: E402
import loop_reduction as lr  # noqa: E402
from fj_direct_engine import (  # noqa: E402
    Engine, vertex_terms, contract, expand_q, vector_propagator_pieces,
    multiplicity_factor, conj, VECTOR_LINE, SCALAR_LINE, one_point_h, instantiate,
)

OUT_PATH = ROOT / "outputs" / "fj_direct_scalar_blocks.json"

g1, g2, v = fj_spectrum.g1, fj_spectrum.g2, fj_spectrum.v
lam, mu2, xi = fj_spectrum.lam, fj_spectrum.mu2, fj_spectrum.xi
d = lr.d_sym
p2 = lr.p2
k2 = lr.k2
kp = lr.kp

MU, MD, ME = fj_spectrum.MU, fj_spectrum.MD, fj_spectrum.ME
VCKM, VCKMC = fj_spectrum.VCKM, fj_spectrum.VCKMC


def reduce_scalar(terms: list[tuple[sp.Expr, tuple]], internal: set[str],
                  a_power: int, b_power: int) -> sp.Expr:
    closed = contract(expand_q(terms), internal)
    total = sp.Integer(0)
    for coeff, car_mu, car_nu in closed:
        if car_mu is not None or car_nu is not None:
            raise ValueError("open index in a scalar block")
        total += coeff
    return lr.reduce_two_point(sp.expand(total), a_power, b_power)


def reduce_mixed(terms: list[tuple[sp.Expr, tuple]], internal: set[str],
                 a_power: int, b_power: int) -> sp.Expr:
    """p_mu contraction of a single-index block; returns p2 Sigma."""

    closed = contract(expand_q(terms), internal)
    total = sp.Integer(0)
    for coeff, car_mu, car_nu in closed:
        if car_nu is not None:
            raise ValueError("second open index in a mixing block")
        if car_mu is None:
            raise ValueError("missing open index in a mixing block")
        total += coeff * lr.dot(car_mu, "p")
    return lr.reduce_two_point(sp.expand(total), a_power, b_power)


def scalar_vertex_terms(record: dict[str, Any], slots: list[tuple[str, str | None, tuple[str, int]]]):
    """Vertex factor for the scalar-sector structures; falls back to
    the shared builder for the gauge structures."""

    structure = record["structure"]
    c = record["coefficient"]
    mult = multiplicity_factor(record["fields"])
    if structure in ("scalar_potential", "ghost_scalar_mass"):
        return [(sp.I * c * mult, ())]
    return vertex_terms(record, slots)


class ScalarBlock:
    """One scalar or mixing block; ext_vector is None for scalar
    externals, else the vector label carrying the open index mu."""

    def __init__(self, engine: Engine, ext1: str, ext2: str) -> None:
        self.e = engine
        self.ext1 = ext1
        self.ext2 = ext2
        self.mixed = ext1 in VECTOR_LINE
        self.diagrams: list[dict[str, Any]] = []

    def emit(self, kind: str, internal: list[str], weight: sp.Expr, reduced: sp.Expr,
             mass_1: sp.Expr, mass_2: sp.Expr | None) -> None:
        value = sp.expand(-sp.I * weight * reduced)
        if value == 0:
            return
        mass_map = {lr.m1sq: sp.simplify(mass_1)}
        if mass_2 is not None:
            mass_map[lr.m2sq] = sp.simplify(mass_2)
        self.diagrams.append({
            "kind": kind, "internal": internal,
            "canonical": value, "mass_map": mass_map,
        })

    def trilinear(self, fields: list[str]):
        order = ("yang_mills_three_point", "scalar_gauge_gauge", "scalar_scalar_gauge",
                 "scalar_potential")
        for structure in order:
            for record in self.e.by_multiset.get(tuple(sorted(fields)), []):
                if record["structure"] == structure:
                    return record
        return None

    # -- boson bubbles ----------------------------------------------------

    def boson_bubbles(self) -> None:
        import itertools
        for x, y in itertools.combinations_with_replacement(sorted(VECTOR_LINE | SCALAR_LINE), 2):
            v1 = self.trilinear([self.ext1, x, y])
            v2 = self.trilinear([self.ext2, conj(x), conj(y)])
            if v1 is None or v2 is None:
                continue
            symmetry = sp.Rational(1, 2) if (x == y and x in ("Z", "A", "h", "G0")) else sp.Integer(1)
            x_vec = x in VECTOR_LINE
            y_vec = y in VECTOR_LINE
            pieces_x = vector_propagator_pieces(self.e.spectrum[x]) if x_vec \
                else [{"mass": self.e.spectrum[x], "den_power": 1, "g": None, "kk": None}]
            pieces_y = vector_propagator_pieces(self.e.spectrum[y]) if y_vec \
                else [{"mass": self.e.spectrum[y], "den_power": 1, "g": None, "kk": None}]
            ext1_slot = (self.ext1, "mu" if self.mixed else None, ("p", 1))
            v1_terms = scalar_vertex_terms(v1, [
                ext1_slot,
                (x, "r1i" if x_vec else None, ("k", -1)),
                (y, "r2i" if y_vec else None, ("q", 1)),
            ])
            v2_terms = scalar_vertex_terms(v2, [
                (self.ext2, None, ("p", -1)),
                (conj(x), "r1j" if x_vec else None, ("k", 1)),
                (conj(y), "r2j" if y_vec else None, ("q", -1)),
            ])
            if v1_terms is None or v2_terms is None:
                continue
            for piece_x in pieces_x:
                for piece_y in pieces_y:
                    prop_x = self.prop_terms(piece_x, "r1i", "r1j", "k") if x_vec else [(sp.Integer(1), ())]
                    prop_y = self.prop_terms(piece_y, "r2i", "r2j", "q") if y_vec else [(sp.Integer(1), ())]
                    combined = [
                        (c1 * c2 * cx * cy, f1 + f2 + fx + fy)
                        for c1, f1 in v1_terms for c2, f2 in v2_terms
                        for cx, fx in prop_x for cy, fy in prop_y
                    ]
                    internal = set()
                    if x_vec:
                        internal |= {"r1i", "r1j"}
                    if y_vec:
                        internal |= {"r2i", "r2j"}
                    reducer = reduce_mixed if self.mixed else reduce_scalar
                    reduced = reducer(combined, internal, piece_x["den_power"], piece_y["den_power"])
                    weight = sp.I * (-sp.I if x_vec else sp.I) * (-sp.I if y_vec else sp.I) * symmetry
                    self.emit("bubble", [x, y], weight, reduced, piece_x["mass"], piece_y["mass"])

    def prop_terms(self, piece: dict[str, Any], idx_a: str, idx_b: str, vec: str):
        out = []
        if piece["g"]:
            out.append((piece["g"], (("g", idx_a, idx_b),)))
        if piece["kk"]:
            out.append((piece["kk"], (("mom", vec, idx_a), ("mom", vec, idx_b))))
        return out

    # -- ghost bubbles ----------------------------------------------------

    def ghost_bubbles(self) -> None:
        # Ordered ghost pairs (X, Y) as in the vector blocks: V1 =
        # {ext1, Xbar, Y}, V2 = {ext2, Ybar, X}; line 1 carries the Y
        # line (mass on D1), line 2 the X line.
        ghosts = ("cp", "cm", "cZ", "cA")
        for x in ghosts:
            for y in ghosts:
                structure_1 = "ghost_gauge_derivative" if self.mixed else "ghost_scalar_mass"
                c1 = self.e.lookup([self.ext1, x + "_bar", y], structure_1)
                c2 = self.e.lookup([self.ext2, y + "_bar", x], "ghost_scalar_mass")
                if c1 == 0 or c2 == 0:
                    continue
                if self.mixed:
                    terms = [(c1 * (sp.I * c2), (("mom", "q", "mu"),))]
                    reduced = reduce_mixed(terms, set(), 1, 1)
                else:
                    terms = [((sp.I * c1) * (sp.I * c2), ())]
                    reduced = reduce_scalar(terms, set(), 1, 1)
                weight = sp.I * (sp.I * sp.I) * (-1)
                self.emit(f"ghost_bubble_{x}{y}", [y, x], weight, reduced,
                          self.e.spectrum[y], self.e.spectrum[x])

    # -- fermion bubbles --------------------------------------------------

    def scalar_scalar_trace(self, a1, b1, a2, b2, m1, m2):
        return sp.expand(2 * (a1 * b2 + b1 * a2) * (k2 - kp) + 2 * (a1 * a2 + b1 * b2) * m1 * m2)

    def mixed_trace_terms(self, a1, b1, a2, b2, m1, m2):
        # gamma^mu (a1 P_L + b1 P_R)(ksl + m1)(a2 P_L + b2 P_R)(qsl + m2)
        return [
            (2 * m2 * (a1 * b2 + b1 * a2), (("mom", "k", "mu"),)),
            (2 * m1 * (a1 * a2 + b1 * b2), (("mom", "q", "mu"),)),
        ]

    def yukawa_pair(self, scalar: str, template: str, gen: int) -> tuple[sp.Expr, sp.Expr]:
        """(P_L, P_R) couplings of a neutral scalar to a mass-basis
        template generation, from the chart."""

        mass = {"u": MU, "d": MD, "e": ME}[template][gen]
        if scalar == "h":
            return (-mass / v, -mass / v)
        if scalar == "G0":
            sign = 1 if template == "u" else -1
            return (-sign * sp.I * mass / v, sign * sp.I * mass / v)
        raise ValueError(scalar)

    def fermion_bubbles(self) -> None:
        weight_unit = sp.I * (-1) * (sp.I ** 2) * (sp.I ** 2)
        neutral_scalars = {self.ext1, self.ext2} <= {"h", "G0"}
        if neutral_scalars:
            for template, colors in (("u", 3), ("d", 3), ("e", 1)):
                for gen in (1, 2, 3):
                    a1, b1 = self.yukawa_pair(self.ext1, template, gen)
                    a2, b2 = self.yukawa_pair(self.ext2, template, gen)
                    mass = {"u": MU, "d": MD, "e": ME}[template][gen]
                    reduced = lr.reduce_two_point(
                        self.scalar_scalar_trace(a1, b1, a2, b2, mass, mass), 1, 1)
                    self.emit(f"fermion_{template}{gen}", [f"{template}{gen}"],
                              weight_unit * colors, reduced, mass ** 2, mass ** 2)
        if {self.ext1, self.ext2} == {"Gp", "Gm"}:
            s2 = sp.sqrt(2)
            for i in (1, 2, 3):
                for j in (1, 2, 3):
                    # G+ (ubar_i d_j): (sqrt2/v)(mfu_i V P_L - V mfd_j P_R)
                    a1 = s2 / v * MU[i] * VCKM[(i, j)]
                    b1 = -s2 / v * VCKM[(i, j)] * MD[j]
                    # G- (dbar_j u_i): conjugate couplings
                    a2 = -s2 / v * VCKMC[(i, j)] * MD[j]
                    b2 = s2 / v * MU[i] * VCKMC[(i, j)]
                    reduced = lr.reduce_two_point(
                        self.scalar_scalar_trace(a1, b1, a2, b2, MU[i], MD[j]), 1, 1)
                    self.emit(f"fermion_u{i}d{j}", [f"u{i}", f"d{j}"],
                              weight_unit * 3, reduced, MU[i] ** 2, MD[j] ** 2)
            for i in (1, 2, 3):
                a1, b1 = sp.Integer(0), -s2 / v * ME[i]
                a2, b2 = -s2 / v * ME[i], sp.Integer(0)
                reduced = lr.reduce_two_point(
                    self.scalar_scalar_trace(a1, b1, a2, b2, sp.Integer(0), ME[i]), 1, 1)
                self.emit(f"fermion_nu{i}e{i}", [f"nu{i}", f"e{i}"],
                          weight_unit, reduced, sp.Integer(0), ME[i] ** 2)
        if self.mixed and {self.ext1} <= {"Z", "A"} and self.ext2 == "G0":
            for template, colors in (("u", 3), ("d", 3), ("e", 1)):
                l_coup, r_coup = self.e.neutral_fermion_coupling(template, self.ext1)
                if l_coup == 0 and r_coup == 0:
                    continue
                for gen in (1, 2, 3):
                    a2, b2 = self.yukawa_pair("G0", template, gen)
                    mass = {"u": MU, "d": MD, "e": ME}[template][gen]
                    terms = self.mixed_trace_terms(l_coup, r_coup, a2, b2, mass, mass)
                    reduced = reduce_mixed(terms, set(), 1, 1)
                    self.emit(f"fermion_{template}{gen}", [f"{template}{gen}"],
                              weight_unit * colors, reduced, mass ** 2, mass ** 2)
        if self.mixed and self.ext1 == "Wp" and self.ext2 == "Gm":
            s2 = sp.sqrt(2)
            coupling = g2 / s2
            for i in (1, 2, 3):
                for j in (1, 2, 3):
                    a1, b1 = coupling * VCKM[(i, j)], sp.Integer(0)
                    a2 = -s2 / v * VCKMC[(i, j)] * MD[j]
                    b2 = s2 / v * MU[i] * VCKMC[(i, j)]
                    terms = self.mixed_trace_terms(a1, b1, a2, b2, MU[i], MD[j])
                    reduced = reduce_mixed(terms, set(), 1, 1)
                    self.emit(f"fermion_u{i}d{j}", [f"u{i}", f"d{j}"],
                              weight_unit * 3, reduced, MU[i] ** 2, MD[j] ** 2)
            for i in (1, 2, 3):
                a1, b1 = coupling, sp.Integer(0)
                a2, b2 = -s2 / v * ME[i], sp.Integer(0)
                terms = self.mixed_trace_terms(a1, b1, a2, b2, sp.Integer(0), ME[i])
                reduced = reduce_mixed(terms, set(), 1, 1)
                self.emit(f"fermion_nu{i}e{i}", [f"nu{i}", f"e{i}"],
                          weight_unit, reduced, sp.Integer(0), ME[i] ** 2)

    # -- seagulls ---------------------------------------------------------

    def seagulls(self) -> None:
        if self.mixed:
            return
        for loop_field in sorted(VECTOR_LINE | SCALAR_LINE):
            partner = conj(loop_field)
            if partner < loop_field:
                continue
            fields = tuple(sorted([self.ext1, self.ext2, loop_field, partner]))
            symmetry = sp.Rational(1, 2) if loop_field == partner else sp.Integer(1)
            for record in self.e.by_multiset.get(fields, []):
                mult = multiplicity_factor(record["fields"])
                if record["structure"] == "scalar_potential" and loop_field in SCALAR_LINE:
                    terms = [(sp.I * record["coefficient"] * mult, ())]
                    reduced = lr.reduce_tadpole(sum(c for c, _f in terms))
                    weight = sp.I * sp.I * symmetry
                    self.emit("seagull", [loop_field, partner], weight, reduced,
                              self.e.spectrum[loop_field], None)
                if record["structure"] == "scalar_scalar_gauge_gauge" and loop_field in VECTOR_LINE:
                    for piece in vector_propagator_pieces(self.e.spectrum[loop_field]):
                        terms = []
                        base = sp.I * record["coefficient"] * mult
                        if piece["g"]:
                            terms.append((base * piece["g"] * d, ()))
                        if piece["kk"]:
                            terms.append((base * piece["kk"] * k2, ()))
                        total = sp.expand(sum(c for c, _f in terms))
                        if piece["den_power"] == 1:
                            reduced = lr.reduce_tadpole(total)
                        else:
                            reduced = lr.reduce_two_point(total, piece["den_power"], 0)
                        weight = sp.I * (-sp.I) * symmetry
                        self.emit("vector_seagull", [loop_field, partner], weight, reduced,
                                  piece["mass"], None)

    # -- FJ tadpole insertion ---------------------------------------------

    def tadpole_insertion(self, t_total: sp.Expr) -> None:
        if self.mixed:
            head_record = self.trilinear([self.ext1, self.ext2, "h"])
            if head_record is None or head_record["structure"] != "scalar_scalar_gauge":
                return
            head_terms = scalar_vertex_terms(head_record, [
                (self.ext1, "mu", ("p", 1)),
                (self.ext2, None, ("p", -1)),
                ("h", None, ("p", 0)),
            ])
            total = sp.Integer(0)
            for coeff, factors in head_terms:
                if factors and factors[0][0] == "mom":
                    total += coeff * (lr.p2 if factors[0][1] == "p" else sp.Integer(0))
            scalar_value = sp.expand(-sp.I * total * (sp.I / (0 - self.e.mh_sq)) * (sp.I * t_total) / lr.p2 * lr.p2)
            self.diagrams.append({
                "kind": "tadpole_insertion", "internal": ["h"],
                "canonical": scalar_value, "mass_map": {},
            })
            return
        head = self.e.lookup(sorted([self.ext1, self.ext2, "h"]), "scalar_potential")
        if head == 0:
            return
        mult = multiplicity_factor(sorted([self.ext1, self.ext2, "h"]))
        scalar_value = sp.expand(-sp.I * (sp.I * head * mult) * (sp.I / (0 - self.e.mh_sq)) * (sp.I * t_total))
        self.diagrams.append({
            "kind": "tadpole_insertion", "internal": ["h"],
            "canonical": scalar_value, "mass_map": {},
        })

    def compute(self, t_total: sp.Expr) -> dict[str, Any]:
        self.boson_bubbles()
        self.ghost_bubbles()
        self.fermion_bubbles()
        self.seagulls()
        self.tadpole_insertion(t_total)
        pole = sp.Integer(0)
        per_diagram = []
        for diagram in self.diagrams:
            dp = lr.uv_pole(diagram["canonical"], diagram["mass_map"])
            pole += dp
            per_diagram.append({
                "kind": diagram["kind"], "internal": diagram["internal"],
                "value": str(instantiate(diagram["canonical"], diagram["mass_map"])),
                "pole": str(sp.simplify(dp)),
            })
        divide = lr.p2 if self.mixed else sp.Integer(1)
        return {
            "diagram_count": len(per_diagram),
            "diagrams": per_diagram,
            "pole": str(sp.simplify(pole / divide)),
            "_pole": sp.simplify(pole / divide),
            "normalization": "Sigma_mix with Pi^mu = p^mu Sigma_mix" if self.mixed else "Pi(p2)",
        }


def main() -> int:
    engine = Engine()
    t_total, _t_pole, _report = one_point_h(engine)
    payload: dict[str, Any] = {
        "schema": "fj_direct_scalar_blocks.v1",
        "target": "FJ_DIRECT_1",
        "units": "loop measure i/(16 pi^2) stripped; Delta is the single 1/eps pole unit",
        "blocks": {},
    }
    results = {}
    for name, (ext1, ext2) in (
        ("hh", ("h", "h")), ("G0G0", ("G0", "G0")), ("GpGm", ("Gp", "Gm")),
        ("WpGm", ("Wp", "Gm")), ("ZG0", ("Z", "G0")), ("AG0", ("A", "G0")),
    ):
        block = ScalarBlock(engine, ext1, ext2)
        summary = block.compute(t_total)
        results[name] = summary
        payload["blocks"][name] = {k: v for k, v in summary.items() if not k.startswith("_")}
        print(json.dumps({"block": name, "diagrams": summary["diagram_count"]}))

    # Custodial control: at g1 -> 0, V = 1, mfu = mfd the neutral and
    # charged Goldstone blocks coincide exactly.
    custodial_map = {g1: 0}
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            custodial_map[VCKM[(i, j)]] = sp.Integer(1) if i == j else sp.Integer(0)
            custodial_map[VCKMC[(i, j)]] = sp.Integer(1) if i == j else sp.Integer(0)
        custodial_map[MD[i]] = MU[i]
    g0_pole = results["G0G0"]["_pole"].subs(custodial_map)
    gp_pole = results["GpGm"]["_pole"].subs(custodial_map)
    difference = sp.simplify(sp.expand(g0_pole - gp_pole))
    payload["controls"] = {
        "custodial_goldstone_poles": {
            "difference": str(difference),
            "passed": difference == 0,
        },
        "ag0_mixing_present_for_st_replay": {
            "pole": results["AG0"]["pole"],
            "statement": (
                "consumed by the WARD_ST_NIELSEN_1 replay together with "
                "the AZ longitudinal part and the ghost mixing"
            ),
        },
    }
    OUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"custodial_passed": bool(difference == 0), "difference": str(difference)[:200]}))
    return 0 if difference == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
