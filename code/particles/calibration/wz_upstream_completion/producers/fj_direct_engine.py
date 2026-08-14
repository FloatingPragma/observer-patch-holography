#!/usr/bin/env python3
"""Workstream E stage 1: direct FJ one-loop engine, vector blocks.

Computes the one-loop two-point functions of the four mass-basis
vector blocks {W+W-, ZZ, AZ, AA} and the h one-point block directly
from the mass-basis specialization of the certified rule table, with
explicit tadpole topologies kept per the FJ contract.  Every vertex
factor, propagator, combinatorial weight and loop numerator is
assembled mechanically from the records; the reduction to the scalar
basis and the UV poles use the validated loop_reduction core.

Conventions (recorded in the payload):

* Routing: p enters vertex 1 and leaves vertex 2; line 1 carries k
  from V1 to V2 (D1 = k^2 - m1sq); line 2 carries q = k - p from V2
  to V1 (D2 = q^2 - m2sq).  Vertex momenta are incoming: at V1
  (p, -k, +q), at V2 (-p, +k, -q).
* Every loop integral carries the measure factor i/(16 pi^2) with the
  16 pi^2 stripped from the emitted coefficients; a diagram value is
  Pi = -I x weight x (reduced numerator), where the weight collects
  the measure i, the propagator prefactors (-i per vector piece, +i
  per scalar or ghost line, i^2 per fermion line pair inside the
  standard trace), the closed-loop sign (-1 for fermion and ghost
  loops), the color factor (N_c = 3 for quark templates), and the
  symmetry factor (1/2 for a self-conjugate identical pair).
* Vertex factors are i x (momentum-space functional derivative of the
  record term): metric-type records carry i x coefficient x
  multiplicity; single-derivative records (scalar_scalar_gauge on the
  sorted scalar pair, ghost_gauge_derivative on the incoming
  antighost momentum) carry coefficient x multiplicity x momentum,
  the i cancelling against the derivative's -i; the Yang-Mills cubic
  tensor is the role sum sgn(sigma) (p_A)^{idx B} g^{idx A, idx C}
  with real prefactor.
* R_xi vector propagators split exactly by partial fractions; on a
  massless line the squared denominator stays in the reduction basis.
* Fermion loops use the standard four-dimensional trace identities;
  epsilon-tensor terms cannot survive a two-point projection, which
  is the frozen BMHV statement for this stage.
* The quark sector uses the mass-basis chart of fj_spectrum with
  V[i][j] symbols; leptons are diagonal with massless neutrinos.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import fj_spectrum  # noqa: E402
import loop_reduction as lr  # noqa: E402

OUT_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
MATCHING_PATH = ROOT / "outputs" / "eft_matching_1.json"

g1, g2, v = fj_spectrum.g1, fj_spectrum.g2, fj_spectrum.v
lam, mu2, xi = fj_spectrum.lam, fj_spectrum.mu2, fj_spectrum.xi
d = lr.d_sym
p2 = lr.p2
k2 = lr.k2
kp = lr.kp

SELF_CONJUGATE = {"Z", "A", "h", "G0", "Gl"}
CONJ = {"Wp": "Wm", "Wm": "Wp", "Gp": "Gm", "Gm": "Gp"}
VECTOR_LINE = {"Wp", "Wm", "Z", "A"}
SCALAR_LINE = {"h", "G0", "Gp", "Gm"}


def conj(label: str) -> str:
    if label in CONJ:
        return CONJ[label]
    if label.endswith("_bar"):
        return label[:-4]
    if label in SELF_CONJUGATE:
        return label
    return label + "_bar"


def multiplicity_factor(fields: list[str]) -> int:
    out = 1
    for label in set(fields):
        for n in range(2, fields.count(label) + 1):
            out *= n
    return out


# ---------------------------------------------------------------------------
# Index contraction engine
# ---------------------------------------------------------------------------


def contract(terms: list[tuple[sp.Expr, tuple]], internal: set[str]) -> list[tuple[sp.Expr, str | None, str | None]]:
    """Contract all internal indices of ("mom", vec, idx) / ("g", i, j)
    factor products; return project()-ready external terms."""

    results: list[tuple[sp.Expr, str | None, str | None]] = []
    for coefficient, factors in terms:
        work = list(factors)
        coeff = coefficient
        changed = True
        while changed:
            changed = False
            for i, factor in enumerate(work):
                if factor[0] != "g":
                    continue
                _, a, b = factor
                if a in internal and b in internal and a == b:
                    coeff *= d
                    work.pop(i)
                    changed = True
                    break
                idx = a if a in internal else (b if b in internal else None)
                if idx is None:
                    continue
                other = b if idx == a else a
                partner = None
                for j, cand in enumerate(work):
                    if j != i and idx in cand[1:]:
                        partner = j
                        break
                if partner is None:
                    raise ValueError(f"dangling internal index {idx}")
                cand = work[partner]
                if cand[0] == "mom":
                    repl = ("mom", cand[1], other)
                else:
                    far = cand[2] if cand[1] == idx else cand[1]
                    if far == idx:
                        coeff *= d
                        repl = None
                    else:
                        repl = ("g", other, far)
                work = [f for n, f in enumerate(work) if n not in (i, partner)]
                if repl is not None:
                    work.append(repl)
                changed = True
                break
            if changed:
                continue
            for i, factor in enumerate(work):
                if factor[0] != "mom" or factor[2] not in internal:
                    continue
                idx = factor[2]
                for j in range(i + 1, len(work)):
                    cand = work[j]
                    if cand[0] == "mom" and cand[2] == idx:
                        coeff *= lr.dot(factor[1], cand[1])
                        work = [f for n, f in enumerate(work) if n not in (i, j)]
                        changed = True
                        break
                if changed:
                    break
        car_mu = car_nu = None
        metric = False
        for factor in work:
            if factor[0] == "mom":
                if factor[2] == "mu":
                    car_mu = factor[1]
                elif factor[2] == "nu":
                    car_nu = factor[1]
                else:
                    raise ValueError(f"unresolved index in {factor}")
            else:
                if {factor[1], factor[2]} == {"mu", "nu"}:
                    metric = True
                else:
                    raise ValueError(f"unresolved metric {factor}")
        if metric:
            results.append((sp.expand(coeff), None, None))
        else:
            results.append((sp.expand(coeff), car_mu, car_nu))
    return results


def expand_q(terms: list[tuple[sp.Expr, tuple]]) -> list[tuple[sp.Expr, tuple]]:
    out: list[tuple[sp.Expr, tuple]] = []

    def rec(coeff: sp.Expr, done: tuple, todo: tuple) -> None:
        if not todo:
            out.append((coeff, done))
            return
        head, tail = todo[0], todo[1:]
        if head[0] == "mom" and head[1] == "q":
            rec(coeff, done + (("mom", "k", head[2]),), tail)
            rec(-coeff, done + (("mom", "p", head[2]),), tail)
        else:
            rec(coeff, done + (head,), tail)

    for coefficient, factors in terms:
        rec(coefficient, (), tuple(factors))
    return out


def reduce_projections(terms: list[tuple[sp.Expr, tuple]], internal: set[str],
                       a_power: int, b_power: int) -> tuple[sp.Expr, sp.Expr]:
    closed = contract(expand_q(terms), internal)
    g_proj = lr.project(closed, "g")
    p_proj = lr.project(closed, "p")
    return (lr.reduce_two_point(g_proj, a_power, b_power),
            lr.reduce_two_point(p_proj, a_power, b_power))


def vector_propagator_pieces(mass_sq: sp.Expr) -> list[dict[str, Any]]:
    if mass_sq == 0:
        return [
            {"mass": sp.Integer(0), "den_power": 1, "g": sp.Integer(1), "kk": sp.Integer(0)},
            {"mass": sp.Integer(0), "den_power": 2, "g": sp.Integer(0), "kk": -(1 - xi)},
        ]
    return [
        {"mass": mass_sq, "den_power": 1, "g": sp.Integer(1), "kk": -1 / mass_sq},
        {"mass": xi * mass_sq, "den_power": 1, "g": sp.Integer(0), "kk": 1 / mass_sq},
    ]


# ---------------------------------------------------------------------------
# Engine
# ---------------------------------------------------------------------------


class Engine:
    def __init__(self) -> None:
        self.records = fj_spectrum.specialized_records()
        self.spectrum, checks = fj_spectrum.spectrum_and_checks(self.records)
        if checks:
            raise SystemExit("fj engine: specialization checks failed: " + "; ".join(checks))
        self.by_multiset: dict[tuple, list[dict[str, Any]]] = {}
        for record in self.records:
            self.by_multiset.setdefault(tuple(sorted(record["fields"])), []).append(record)
        self.mh_sq = self.spectrum["h"]

    def lookup(self, fields: list[str], structure: str) -> sp.Expr:
        for record in self.by_multiset.get(tuple(sorted(fields)), []):
            if record["structure"] == structure:
                return record["coefficient"]
        return sp.Integer(0)

    def neutral_fermion_coupling(self, template: str, vector: str) -> tuple[sp.Expr, sp.Expr]:
        left = self.lookup([f"{template}L_bar", f"{template}L", vector], "fermion_vector_current")
        right = self.lookup([f"{template}R_bar", f"{template}R", vector], "fermion_vector_current")
        return sp.simplify(left), sp.simplify(right)


# vertex slot: (label, lorentz_index_or_None, (vector, sign)) with the
# momentum given as incoming.


def vertex_terms(record: dict[str, Any], slots: list[tuple[str, str | None, tuple[str, int]]]):
    """Feynman factor of one trilinear record on the given slots."""

    structure = record["structure"]
    c = record["coefficient"]
    mult = multiplicity_factor(record["fields"])
    if structure == "yang_mills_three_point":
        assignment = {label: (idx, mom) for label, idx, mom in slots}
        terms = []
        fields = record["fields"]
        for perm in itertools.permutations(range(3)):
            sign = sp.Integer(1)
            for x in range(3):
                for y in range(x + 1, 3):
                    if perm[x] > perm[y]:
                        sign = -sign
            ordered = [fields[perm.index(r)] for r in range(3)]
            f_a, f_b, f_c = ordered
            idx_a, mom_a = assignment[f_a]
            idx_b, _ = assignment[f_b]
            idx_c, _ = assignment[f_c]
            vec, msign = mom_a
            terms.append((c * mult_free(record) * sign * msign, (("mom", vec, idx_b), ("g", idx_a, idx_c))))
        return terms
    if structure == "scalar_gauge_gauge":
        index_slots = [idx for label, idx, _mom in slots if label in VECTOR_LINE]
        if len(index_slots) != 2:
            return None
        return [(sp.I * c * mult, (("g", index_slots[0], index_slots[1]),))]
    if structure == "scalar_scalar_gauge":
        scalars = sorted(label for label, _i, _m in slots if label in SCALAR_LINE)
        if len(scalars) != 2 or scalars[0] == scalars[1]:
            return None
        vector_slot = next(idx for label, idx, _m in slots if label in VECTOR_LINE)
        mom_first = next(m for label, _i, m in slots if label == scalars[0])
        mom_second = next(m for label, _i, m in slots if label == scalars[1])
        return [
            (c * mult * mom_first[1], (("mom", mom_first[0], vector_slot),)),
            (-c * mult * mom_second[1], (("mom", mom_second[0], vector_slot),)),
        ]
    return None


def mult_free(record: dict[str, Any]) -> int:
    # the cubic record has three distinct fields, so the multiplicity
    # factor is one; kept explicit in the provenance record.
    return multiplicity_factor(record["fields"])


class BlockComputer:
    def __init__(self, engine: Engine, ext1: str, ext2: str) -> None:
        self.e = engine
        self.ext1 = ext1
        self.ext2 = ext2
        self.diagrams: list[dict[str, Any]] = []

    def emit(self, kind: str, internal: list[str], weight: sp.Expr,
             g_red: sp.Expr, p_red: sp.Expr, mass_1: sp.Expr, mass_2: sp.Expr | None) -> None:
        g_final = sp.expand(-sp.I * weight * g_red)
        p_final = sp.expand(-sp.I * weight * p_red)
        if g_final == 0 and p_final == 0:
            return
        mass_map = {lr.m1sq: sp.simplify(mass_1)}
        if mass_2 is not None:
            mass_map[lr.m2sq] = sp.simplify(mass_2)
        self.diagrams.append({
            "kind": kind, "internal": internal,
            "g_canonical": g_final, "p_canonical": p_final,
            "mass_map": mass_map,
        })

    # -- boson bubbles ----------------------------------------------------

    def trilinear(self, fields: list[str]):
        for structure in ("yang_mills_three_point", "scalar_gauge_gauge", "scalar_scalar_gauge"):
            for record in self.e.by_multiset.get(tuple(sorted(fields)), []):
                if record["structure"] == structure:
                    return record
        return None

    def boson_bubbles(self) -> None:
        for x, y in itertools.combinations_with_replacement(sorted(VECTOR_LINE | SCALAR_LINE), 2):
            v1 = self.trilinear([self.ext1, x, y])
            v2 = self.trilinear([self.ext2, conj(x), conj(y)])
            if v1 is None or v2 is None:
                continue
            symmetry = sp.Rational(1, 2) if (x == y and x in SELF_CONJUGATE) else sp.Integer(1)
            x_vec = x in VECTOR_LINE
            y_vec = y in VECTOR_LINE
            pieces_x = vector_propagator_pieces(self.e.spectrum[x]) if x_vec \
                else [{"mass": self.e.spectrum[x], "den_power": 1, "g": None, "kk": None}]
            pieces_y = vector_propagator_pieces(self.e.spectrum[y]) if y_vec \
                else [{"mass": self.e.spectrum[y], "den_power": 1, "g": None, "kk": None}]
            v1_terms = vertex_terms(v1, [
                (self.ext1, "mu", ("p", 1)),
                (x, "r1i" if x_vec else None, ("k", -1)),
                (y, "r2i" if y_vec else None, ("q", 1)),
            ])
            v2_terms = vertex_terms(v2, [
                (self.ext2, "nu", ("p", -1)),
                (conj(x), "r1j" if x_vec else None, ("k", 1)),
                (conj(y), "r2j" if y_vec else None, ("q", -1)),
            ])
            if v1_terms is None or v2_terms is None:
                continue
            for piece_x in pieces_x:
                for piece_y in pieces_y:
                    prop_terms_x = self.prop_terms(piece_x, "r1i", "r1j", "k") if x_vec else [(sp.Integer(1), ())]
                    prop_terms_y = self.prop_terms(piece_y, "r2i", "r2j", "q") if y_vec else [(sp.Integer(1), ())]
                    combined = [
                        (c1 * c2 * cx * cy, f1 + f2 + fx + fy)
                        for c1, f1 in v1_terms for c2, f2 in v2_terms
                        for cx, fx in prop_terms_x for cy, fy in prop_terms_y
                    ]
                    internal = set()
                    if x_vec:
                        internal |= {"r1i", "r1j"}
                    if y_vec:
                        internal |= {"r2i", "r2j"}
                    g_red, p_red = reduce_projections(combined, internal, piece_x["den_power"], piece_y["den_power"])
                    weight = sp.I * (-sp.I if x_vec else sp.I) * (-sp.I if y_vec else sp.I) * symmetry
                    self.emit("bubble", [x, y], weight, g_red, p_red, piece_x["mass"], piece_y["mass"])

    def prop_terms(self, piece: dict[str, Any], idx_a: str, idx_b: str, vec: str):
        out = []
        if piece["g"]:
            out.append((piece["g"], (("g", idx_a, idx_b),)))
        if piece["kk"]:
            out.append((piece["kk"], (("mom", vec, idx_a), ("mom", vec, idx_b))))
        return out

    # -- ghost bubbles ----------------------------------------------------

    def ghost_bubbles(self) -> None:
        # Ordered ghost pairs (X, Y): V1 = {ext1, Xbar, Y}, V2 =
        # {ext2, Ybar, X}; line 1 (momentum k) is the Y line, line 2
        # (momentum q) the X line; the antighost momenta are +q at V1
        # and +k at V2.  Diagonal pairs reproduce the neutral-block
        # loops; the mixed pairs carry the charged-block ghost sector.
        ghosts = ("cp", "cm", "cZ", "cA")
        for x in ghosts:
            for y in ghosts:
                c1 = self.e.lookup([self.ext1, x + "_bar", y], "ghost_gauge_derivative")
                c2 = self.e.lookup([self.ext2, y + "_bar", x], "ghost_gauge_derivative")
                if c1 == 0 or c2 == 0:
                    continue
                terms = [(c1 * c2, (("mom", "q", "mu"), ("mom", "k", "nu")))]
                g_red, p_red = reduce_projections(terms, set(), 1, 1)
                weight = sp.I * (sp.I * sp.I) * (-1)
                self.emit(f"ghost_bubble_{x}{y}", [y, x], weight, g_red, p_red,
                          self.e.spectrum[y], self.e.spectrum[x])

    # -- fermion bubbles --------------------------------------------------

    def fermion_terms(self, a, b, c_, d_, mass_1, mass_2):
        chir = 2 * (a * c_ + b * d_)
        flip = 2 * (a * d_ + b * c_) * mass_1 * mass_2
        return [
            (chir, (("mom", "k", "mu"), ("mom", "q", "nu"))),
            (chir, (("mom", "k", "nu"), ("mom", "q", "mu"))),
            (-chir * (k2 - kp), (("g", "mu", "nu"),)),
            (flip, (("g", "mu", "nu"),)),
        ]

    def fermion_bubbles(self) -> None:
        weight = sp.I * (-1) * (sp.I ** 2) * (sp.I ** 2)  # measure, loop, vertices, propagators
        if self.ext1 in ("Z", "A") and self.ext2 in ("Z", "A"):
            for template, colors, masses in (
                ("u", 3, fj_spectrum.MU), ("d", 3, fj_spectrum.MD),
                ("e", 1, fj_spectrum.ME), ("nu", 1, None),
            ):
                l1, r1 = self.e.neutral_fermion_coupling(template, self.ext1)
                l2, r2 = self.e.neutral_fermion_coupling(template, self.ext2)
                if (l1 == 0 and r1 == 0) or (l2 == 0 and r2 == 0):
                    continue
                for gen in (1, 2, 3):
                    mass = masses[gen] if masses else sp.Integer(0)
                    terms = self.fermion_terms(l1, r1, l2, r2, mass, mass)
                    g_red, p_red = reduce_projections(terms, set(), 1, 1)
                    self.emit(f"fermion_{template}{gen}", [f"{template}{gen}"],
                              weight * colors, g_red, p_red, mass ** 2, mass ** 2)
        if {self.ext1, self.ext2} == {"Wp", "Wm"}:
            coupling = g2 / sp.sqrt(2)
            for i in (1, 2, 3):
                for j in (1, 2, 3):
                    vij = fj_spectrum.VCKM[(i, j)]
                    vcij = fj_spectrum.VCKMC[(i, j)]
                    terms = self.fermion_terms(coupling * vij, 0, coupling * vcij, 0,
                                               fj_spectrum.MU[i], fj_spectrum.MD[j])
                    g_red, p_red = reduce_projections(terms, set(), 1, 1)
                    self.emit(f"fermion_u{i}d{j}", [f"u{i}", f"d{j}"], weight * 3,
                              g_red, p_red, fj_spectrum.MU[i] ** 2, fj_spectrum.MD[j] ** 2)
            for i in (1, 2, 3):
                terms = self.fermion_terms(coupling, 0, coupling, 0, sp.Integer(0), fj_spectrum.ME[i])
                g_red, p_red = reduce_projections(terms, set(), 1, 1)
                self.emit(f"fermion_nu{i}e{i}", [f"nu{i}", f"e{i}"], weight,
                          g_red, p_red, sp.Integer(0), fj_spectrum.ME[i] ** 2)

    # -- seagulls ---------------------------------------------------------

    def seagulls(self) -> None:
        for loop_field in sorted(VECTOR_LINE | SCALAR_LINE):
            partner = conj(loop_field)
            if partner < loop_field:
                continue
            fields = tuple(sorted([self.ext1, self.ext2, loop_field, partner]))
            symmetry = sp.Rational(1, 2) if loop_field == partner else sp.Integer(1)
            for record in self.e.by_multiset.get(fields, []):
                mult = multiplicity_factor(record["fields"])
                if record["structure"] == "scalar_scalar_gauge_gauge" and loop_field in SCALAR_LINE:
                    terms = [(sp.I * record["coefficient"] * mult, (("g", "mu", "nu"),))]
                    closed = contract(expand_q(terms), set())
                    g_red = lr.reduce_tadpole(lr.project(closed, "g"))
                    p_red = lr.reduce_tadpole(lr.project(closed, "p"))
                    weight = sp.I * sp.I * symmetry
                    self.emit("seagull", [loop_field, partner], weight, g_red, p_red,
                              self.e.spectrum[loop_field], None)
                if record["structure"].startswith("yang_mills_four_point") and loop_field in VECTOR_LINE:
                    self.vector_seagull(record, loop_field, partner, symmetry)

    def vector_seagull(self, record: dict[str, Any], loop_field: str, partner: str, symmetry) -> None:
        fields = record["fields"]
        tag = record["structure"].rsplit("_", 2)
        pairing = {"12": (0, 1), "34": (2, 3), "13": (0, 2), "24": (1, 3), "14": (0, 3), "23": (1, 2)}
        first = pairing[tag[1]]
        second = pairing[tag[2]]
        # Functional differentiation: sum over every label-preserving
        # assignment of the record positions to the four leg slots.
        # The pairing tensor acts on record POSITIONS, so assignments
        # with identical field labels generate DISTINCT slot tensors
        # (the degenerate pairing labels of the table merge value
        # partitions; the assignment sum restores the full Feynman
        # rule).  The assignment sum replaces the multiplicity factor.
        import itertools as it
        slot_names = ["mu", "nu", "li", "lj"]
        slot_labels = [self.ext1, self.ext2, loop_field, partner]
        bases = []
        for perm in it.permutations(range(4)):
            if all(fields[k] == slot_labels[perm[k]] for k in range(4)):
                s = [slot_names[perm[k]] for k in range(4)]
                bases.append((("g", s[first[0]], s[first[1]]),
                              ("g", s[second[0]], s[second[1]])))
        if not bases:
            return
        for piece in vector_propagator_pieces(self.e.spectrum[loop_field]):
            terms = []
            for base in bases:
                if piece["g"]:
                    terms.append((sp.I * record["coefficient"] * piece["g"], base + (("g", "li", "lj"),)))
                if piece["kk"]:
                    terms.append((sp.I * record["coefficient"] * piece["kk"],
                                  base + (("mom", "k", "li"), ("mom", "k", "lj"))))
            closed = contract(expand_q(terms), {"li", "lj"})
            g_proj = lr.project(closed, "g")
            p_proj = lr.project(closed, "p")
            if piece["den_power"] == 1:
                g_red = lr.reduce_tadpole(g_proj)
                p_red = lr.reduce_tadpole(p_proj)
            else:
                g_red = lr.reduce_two_point(g_proj, piece["den_power"], 0)
                p_red = lr.reduce_two_point(p_proj, piece["den_power"], 0)
            weight = sp.I * (-sp.I) * symmetry
            self.emit("vector_seagull", [loop_field, partner], weight, g_red, p_red, piece["mass"], None)

    # -- FJ tadpole insertions --------------------------------------------

    def tadpole_insertion(self, t_total: sp.Expr) -> None:
        head = self.e.lookup([self.ext1, self.ext2, "h"], "scalar_gauge_gauge")
        if head == 0:
            return
        mult = multiplicity_factor(sorted([self.ext1, self.ext2, "h"]))
        # Pi_insert = (head vertex) x (i/(0 - mh^2)) x (i T_h), reduced
        # to the final convention Pi = -I x (i c mult) x (i/-mh^2) x (i T).
        scalar_value = sp.expand(-sp.I * (sp.I * head * mult) * (sp.I / (0 - self.e.mh_sq)) * (sp.I * t_total))
        self.diagrams.append({
            "kind": "tadpole_insertion", "internal": ["h"],
            "g_canonical": sp.expand(scalar_value * d), "p_canonical": sp.expand(scalar_value * p2),
            "mass_map": {},
        })

    def compute(self, t_total: sp.Expr) -> None:
        self.boson_bubbles()
        self.ghost_bubbles()
        self.fermion_bubbles()
        self.seagulls()
        self.tadpole_insertion(t_total)


# ---------------------------------------------------------------------------
# One-point block
# ---------------------------------------------------------------------------


def one_point_h(engine: Engine) -> tuple[sp.Expr, sp.Expr, dict[str, str]]:
    """T_h with T = -I x sum of diagram amplitudes; returns the total
    (over instantiated A0[mass] symbols), its UV pole, and the report."""

    def a0(mass_sq: sp.Expr) -> sp.Expr:
        return sp.Symbol(f"A0[{sp.simplify(mass_sq)}]")

    contributions: dict[str, sp.Expr] = {}
    pole = sp.Integer(0)
    tree = engine.lookup(["h"], "scalar_tadpole")
    contributions["tree"] = tree

    def add_loop(name: str, weight: sp.Expr, coefficient: sp.Expr, mass_sq: sp.Expr,
                 d_factor: sp.Expr = sp.Integer(1)) -> None:
        nonlocal pole
        value = sp.expand(-sp.I * weight * coefficient * d_factor.subs(d, 4 - 2 * lr.eps).subs(lr.eps, 0) * a0(mass_sq))
        exact = sp.expand(-sp.I * weight * coefficient * d_factor * a0(mass_sq))
        contributions[name] = exact
        pole_part = sp.expand(-sp.I * weight * coefficient * d_factor * mass_sq * lr.Delta)
        pole_series = sp.series(sp.together(pole_part.subs(d, 4 - 2 * lr.eps).coeff(lr.Delta)), lr.eps, 0, 1).removeO()
        pole += sp.simplify(pole_series.subs(lr.eps, 0))

    for scalar, sym in (("h", sp.Rational(1, 2)), ("G0", sp.Rational(1, 2))):
        c3 = engine.lookup(sorted(["h", scalar, scalar]), "scalar_potential")
        if c3 != 0:
            mult = multiplicity_factor(sorted(["h", scalar, scalar]))
            add_loop(f"scalar_loop_{scalar}", sp.I * sp.I * sym, sp.I * c3 * mult, engine.spectrum[scalar])
    c3 = engine.lookup(["Gm", "Gp", "h"], "scalar_potential")
    if c3 != 0:
        add_loop("scalar_loop_Gpm", sp.I * sp.I, sp.I * c3, engine.spectrum["Gp"])

    for name, pair, sym in (("Wpm", ["Wm", "Wp"], sp.Integer(1)), ("Z", ["Z", "Z"], sp.Rational(1, 2))):
        cv = engine.lookup(sorted(["h"] + pair), "scalar_gauge_gauge")
        if cv == 0:
            continue
        mult = multiplicity_factor(sorted(["h"] + pair))
        for piece in vector_propagator_pieces(engine.spectrum[pair[0]]):
            g_part = (piece["g"] or 0) * d
            kk_part = (piece["kk"] or 0) * k2
            reduced = lr.reduce_tadpole(sp.expand(g_part + kk_part))
            # reduced is coefficient x A0m1 in canonical form
            coefficient = reduced.coeff(lr.BASIS[(1, 0)]).subs(lr.m1sq, piece["mass"])
            if coefficient == 0:
                continue
            add_loop(f"vector_loop_{name}_{piece['mass']}", sp.I * (-sp.I) * sym,
                     sp.I * cv * mult * coefficient, piece["mass"])

    for ghost in ("cp", "cm", "cZ"):
        cg = engine.lookup(sorted(["h", ghost, ghost + "_bar"]), "ghost_scalar_mass")
        if cg != 0:
            add_loop(f"ghost_loop_{ghost}", sp.I * sp.I * (-1), sp.I * cg, engine.spectrum[ghost])

    for template, colors, masses in (("u", 3, fj_spectrum.MU), ("d", 3, fj_spectrum.MD), ("e", 1, fj_spectrum.ME)):
        for gen in (1, 2, 3):
            mf = masses[gen]
            # h coupling -mf/v (P_L + P_R); Tr[(a P_L + b P_R)(ksl+m)] = 2 m (a + b)
            add_loop(f"fermion_loop_{template}{gen}", sp.I * sp.I * (-1) * colors,
                     sp.I * (-mf / v) * 2 * mf, mf ** 2)

    total = sp.expand(sum(contributions.values()))
    report = {name: str(value) for name, value in contributions.items()}
    return total, sp.simplify(pole), report


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------


def instantiate(expression: sp.Expr, mass_map: dict) -> sp.Expr:
    m1v = mass_map.get(lr.m1sq)
    m2v = mass_map.get(lr.m2sq)
    subs: dict[sp.Expr, sp.Expr] = {}
    if m1v is not None:
        subs[sp.Symbol("A0m1")] = sp.Symbol(f"A0[{m1v}]")
        subs[sp.Symbol("A0pm1")] = sp.Symbol(f"A0p[{m1v}]")
    if m2v is not None:
        subs[sp.Symbol("A0m2")] = sp.Symbol(f"A0[{m2v}]")
        subs[sp.Symbol("A0pm2")] = sp.Symbol(f"A0p[{m2v}]")
    if m1v is not None and m2v is not None:
        subs[sp.Symbol("B0")] = sp.Symbol(f"B0[{m1v}|{m2v}]")
        subs[sp.Symbol("C21")] = sp.Symbol(f"C21[{m1v}|{m2v}]")
        subs[sp.Symbol("C12")] = sp.Symbol(f"C12[{m1v}|{m2v}]")
        subs[sp.Symbol("C22")] = sp.Symbol(f"C22[{m1v}|{m2v}]")
    out = expression.subs(subs)
    return sp.expand(out.subs(mass_map))


def block_summary(diagrams: list[dict[str, Any]]) -> dict[str, Any]:
    g_pole = sp.Integer(0)
    p_pole = sp.Integer(0)
    per_diagram = []
    typed_poles: list[tuple] = []
    for diagram in diagrams:
        gp = lr.uv_pole(diagram["g_canonical"], diagram["mass_map"]) if diagram["mass_map"] else \
            lr.uv_pole(diagram["g_canonical"], {})
        pp = lr.uv_pole(diagram["p_canonical"], diagram["mass_map"]) if diagram["mass_map"] else \
            lr.uv_pole(diagram["p_canonical"], {})
        g_pole += gp
        p_pole += pp
        typed_poles.append((diagram["kind"] if not (diagram["kind"] == "bubble" and all(f in ("h", "G0", "Gp", "Gm") for f in diagram["internal"])) else "bubble_scalar", gp, pp))
        per_diagram.append({
            "kind": diagram["kind"],
            "internal": diagram["internal"],
            "g": str(instantiate(diagram["g_canonical"], diagram["mass_map"])),
            "p": str(instantiate(diagram["p_canonical"], diagram["mass_map"])),
            "g_pole": str(sp.simplify(gp)),
            "p_pole": str(sp.simplify(pp)),
        })
    pi_t_pole = sp.simplify(sp.expand((g_pole - p_pole / p2) / 3))
    pi_l_pole = sp.simplify(p_pole / p2)
    return {
        "diagram_count": len(per_diagram),
        "diagrams": per_diagram,
        "g_pole": str(sp.simplify(g_pole)),
        "p_pole": str(sp.simplify(p_pole)),
        "transverse_pole": str(pi_t_pole),
        "longitudinal_pole": str(pi_l_pole),
        "_pi_t_pole": pi_t_pole,
        "_typed_poles": typed_poles,
    }


def main() -> int:
    engine = Engine()
    t_total, t_pole, t_report = one_point_h(engine)
    payload: dict[str, Any] = {
        "schema": "fj_direct_vector_blocks.v1",
        "target": "FJ_DIRECT_1",
        "units": "loop measure i/(16 pi^2) stripped; Delta is the single 1/eps pole unit",
        "one_point_h": {"contributions": t_report, "uv_pole": str(t_pole)},
        "blocks": {},
    }
    blocks = {}
    for name, (ext1, ext2) in (("AA", ("A", "A")), ("AZ", ("A", "Z")),
                               ("ZZ", ("Z", "Z")), ("WpWm", ("Wp", "Wm"))):
        computer = BlockComputer(engine, ext1, ext2)
        computer.compute(t_total)
        summary = block_summary(computer.diagrams)
        blocks[name] = summary
        payload["blocks"][name] = {k: v for k, v in summary.items() if not k.startswith("_")}
        summary["_name"] = name
        print(json.dumps({"block": name, "diagrams": summary["diagram_count"],
                          "T_pole": summary["transverse_pole"][:120]}))

    # Validation controls.  The photon block decomposes into three
    # exactly known sectors; the transversality of the photon is exact.
    from fractions import Fraction
    from vertex_format import FIELD_TABLE
    aa = blocks["AA"]
    e_sq = g1 ** 2 * g2 ** 2 / (g1 ** 2 + g2 ** 2)
    fermion_pole = sp.Integer(0)
    scalar_pole = sp.Integer(0)
    gauge_pole = sp.Integer(0)
    for diagram, gp, pp in aa["_typed_poles"]:
        t_part = sp.simplify((gp - pp / p2) / 3)
        if diagram.startswith("fermion"):
            fermion_pole += t_part
        elif diagram in ("bubble_scalar", "seagull"):
            scalar_pole += t_part
        else:
            gauge_pole += t_part
    # Sum of N_c Q^2 over the Dirac templates, from the component table.
    charge_sum = sp.Integer(0)
    for template, colors in (("uL", 3), ("dL", 3), ("eL", 1)):
        q = sp.Rational(Fraction(FIELD_TABLE[template]["charge"]))
        charge_sum += colors * q ** 2 * 3  # three generations
    controls = {}
    controls["photon_longitudinal_pole_vanishes"] = {
        "value": aa["longitudinal_pole"],
        "passed": sp.simplify(sp.sympify(aa["longitudinal_pole"])) == 0,
    }
    expected_fermion = sp.simplify(-sp.Rational(4, 3) * charge_sum * e_sq * p2)
    controls["photon_fermion_sector"] = {
        "expected": str(expected_fermion),
        "engine": str(sp.simplify(fermion_pole)),
        "passed": sp.simplify(fermion_pole - expected_fermion) == 0,
    }
    expected_scalar = sp.simplify(-e_sq * p2 / 3)
    scalar_p2 = sp.simplify(sp.expand(scalar_pole).coeff(p2) * p2)
    controls["photon_scalar_sector_p2"] = {
        "expected": str(expected_scalar),
        "engine": str(scalar_p2),
        "passed": sp.simplify(scalar_p2 - expected_scalar) == 0,
    }
    expected_gauge = sp.simplify(2 * (sp.Rational(13, 6) - xi / 2) * e_sq * p2)
    gauge_p2 = sp.simplify(sp.expand(gauge_pole).coeff(p2) * p2)
    controls["photon_gauge_sector_p2"] = {
        "expected": str(expected_gauge),
        "engine": str(gauge_p2),
        "reference": "R_xi vector wave-function formula [13/6 - xi/2] C_A with the machine C_A = 2 of one charged pair",
        "passed": sp.simplify(gauge_p2 - expected_gauge) == 0,
    }
    # Charge-universality census binding: with dZ_AA = -dPi_T^AA/dp2
    # at 0 and dZ_ZA = -(2/mZ^2) Pi_T^AZ(0), the combination
    # dZ_e = -(1/2) dZ_AA - (sw/(2 cw)) dZ_ZA must be xi-independent
    # and equal (b1 + b2)/2 e^2 from the census betas of the matching
    # packet.  The AZ block is NOT transverse in R_xi (the charged
    # gauge-fixing functions break em covariance; the block obeys a
    # Slavnov-Taylor identity replayed in Workstream G), and exactly
    # its zero-momentum value completes the gauge-independent charge.
    matching = json.loads(MATCHING_PATH.read_text(encoding="utf-8"))
    b1 = sp.Rational(Fraction(matching["gauge_betas"]["coefficients"]["b1"]))
    b2 = sp.Rational(Fraction(matching["gauge_betas"]["coefficients"]["b2"]))
    mz_sq = (g1 ** 2 + g2 ** 2) * v ** 2 / 4
    sw_over_cw = g1 / g2
    aa_slope = sp.expand(blocks["AA"]["_pi_t_pole"]).coeff(p2)
    az_zero = sp.simplify(blocks["AZ"]["_pi_t_pole"].subs(p2, 0))
    # Sign dictionary anchored on QED: the fermion sector alone must
    # give dZ_e = +16/3 e^2, which fixes dZ_AA = +dPi/dp2 and
    # dZ_ZA = +2 Pi^AZ(0)/mZ^2 in the engine Pi convention.
    dz_aa = aa_slope
    dz_za = 2 * az_zero / mz_sq
    dz_e = sp.simplify(-sp.Rational(1, 2) * dz_aa - (sw_over_cw / 2) * dz_za)
    expected_dz_e = sp.simplify((b1 + b2) / 2 * e_sq)
    controls["charge_universality_census_binding"] = {
        "dZ_AA_pole": str(sp.simplify(dz_aa)),
        "dZ_ZA_pole": str(sp.simplify(dz_za)),
        "dZ_e_pole": str(dz_e),
        "expected": str(expected_dz_e),
        "xi_independent": not dz_e.has(xi),
        "passed": bool(sp.simplify(dz_e - expected_dz_e) == 0 and not dz_e.has(xi)),
    }
    controls["az_longitudinal_is_st_constrained"] = {
        "value": payload["blocks"]["AZ"]["longitudinal_pole"],
        "statement": (
            "nonzero by R_xi structure: the charged gauge-fixing "
            "functions break em covariance, so the AZ block obeys the "
            "Slavnov-Taylor identity with ghost-mixing and A-G0 "
            "functions instead of naive transversality; the identity "
            "is replayed by the WARD_ST_NIELSEN_1 checker"
        ),
    }
    payload["controls"] = controls
    all_passed = all(c.get("passed", True) for c in controls.values())
    OUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"controls_passed": all_passed}))
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
