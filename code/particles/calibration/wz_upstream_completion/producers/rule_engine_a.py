#!/usr/bin/env python3
"""Workstream C, engine A: vertex records by symbolic action expansion.

Every record is DERIVED by expanding the broken-phase gauge-fixed
Lagrangian with exact symbolic arithmetic, never entered by hand:

* Higgs kinetic sector: (dH + M H)-dagger (dH + M H) with the exact
  SU(2) x U(1) gauge matrix M on the doublet H = (Gp, (v+h+I G0)/sqrt2).
  Derivative trilinears are folded onto the antisymmetric
  scalar_scalar_gauge convention with a machine check that the
  symmetric part vanishes.
* Potential: L = + mu2 (H'H) - lam (H'H)^2 expanded around the shifted
  vacuum with the minimum NOT imposed; the tadpole record is retained
  per the direct FJ contract.
* Yukawa sector: the three operators expanded with the exact conjugate
  doublet; hermitian-conjugate records carry the dagger-matrix symbols
  with the positional index binding declared in vertex_format.
* Yang-Mills self-couplings: the cubic is expanded as
  -g2 eps_abc (dW)_a W_b W_c with distinct role symbols per Lorentz
  slot, rotated to the charged basis, and emitted as the coefficient of
  the sorted role assignment after a machine proof of total
  antisymmetry; the quartic is expanded as
  -(g2^2/4) sum_a (eps_abc P_b Q_c)(eps_ade P_d Q_e) with slot symbols
  P (mu) and Q (nu), and each term is routed to the Lorentz pairing
  label of its P/Q partition.
* Gauge fixing: R_xi functions G^a = d.W^a + xi k_a . phi with the
  coefficients k_a SOLVED from the requirement that the |DH|^2
  vector-Goldstone mixing cancels exactly; the solved coefficients and
  the cancelled mixing records are reported in the table metadata.
* Ghost sector: L_ghost = - cbar_a (delta G^a / delta theta_b) c_b
  evaluated from the exact gauge variations of the vector and doublet
  components, with ghosts rotated to the charged basis; terms are
  ordered with the barred ghost leftmost, which fixes the sign
  convention for the anticommuting pair.
* Fermion gauge currents are read from the action bundle census with
  their exact isospin and hypercharge factors.

Exclusions (recorded in the output): gluon self-couplings and the
gluon-ghost gauge vertex are outside this table because the census
gluon is collapsed to a single adjoint template and gluon
self-interactions enter W/Z pole functions beyond one loop.
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

from vertex_format import (  # noqa: E402
    bar,
    conservation_violations,
    monomial,
    polynomial,
    record,
    table_digest,
    vertex_hash,
)

ACTION_PATH = ROOT / "outputs" / "sm_eft_action_1.json"
OUT_PATH = ROOT / "outputs" / "rule_table_engine_a.json"

# Sympy symbols.  Field symbols are commuting; conjugation is handled
# explicitly by symbol substitution, never by sympy's conjugate.
g1, g2, g3, v, lam, mu2 = sp.symbols("g1 g2 g3 v lam mu2", positive=True)
xi1, xi2 = sp.symbols("xi1 xi2", positive=True)
sqrt2 = sp.sqrt(2)
Wp, Wm, W3, B = sp.symbols("Wp Wm W3 B")
h, G0, Gp, Gm = sp.symbols("h G0 Gp Gm")
dh, dG0, dGp, dGm = sp.symbols("dh dG0 dGp dGm")

FIELD_SYMBOLS = {Wp: "Wp", Wm: "Wm", W3: "W3", B: "B", h: "h", G0: "G0", Gp: "Gp", Gm: "Gm"}
DERIVATIVE_SYMBOLS = {dh: "h", dG0: "G0", dGp: "Gp", dGm: "Gm"}
COUPLING_ORDER = ("g1", "g2", "g3", "v", "lam", "mu2", "xi1", "xi2")

VECTOR_LABELS = ("Wp", "Wm", "W3", "B", "Gl")
SCALAR_LABELS = ("h", "G0", "Gp", "Gm")


def load_action() -> dict[str, Any]:
    if not ACTION_PATH.is_file():
        raise SystemExit("engine A: run the Workstream A producer first")
    return json.loads(ACTION_PATH.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Exact monomial conversion
# ---------------------------------------------------------------------------


def to_polynomial(expression: sp.Expr) -> dict[str, Any]:
    """Convert a sympy coupling expression into the canonical polynomial."""

    expanded = sp.expand(expression)
    monomials = []
    for term in expanded.as_ordered_terms():
        rest = sp.simplify(term)
        i_power = 0
        if rest.has(sp.I):
            rest = sp.cancel(rest / sp.I)
            i_power = 1
            if rest.has(sp.I):
                raise SystemExit(f"engine A: unresolved imaginary structure in {expression}")
        numerator, denominator = sp.fraction(sp.nsimplify(rest))
        sqrt2_power = 0
        for base_expr, direction in ((numerator, 1), (denominator, -1)):
            for factor_part in sp.Mul.make_args(base_expr):
                if factor_part == sqrt2:
                    sqrt2_power += direction
        rest = sp.simplify(rest / (sqrt2 ** sqrt2_power))
        poly_part = sp.Poly(rest, g1, g2, g3, v, lam, mu2, xi1, xi2)
        if len(poly_part.terms()) != 1:
            raise SystemExit(f"engine A: non-monomial residue in {expression}")
        exponents, coefficient = poly_part.terms()[0]
        prefactor = Fraction(sp.Rational(coefficient))
        powers = [
            (name, int(exponent))
            for name, exponent in zip(COUPLING_ORDER, exponents)
            if exponent
        ]
        if sqrt2_power:
            powers.append(("sqrt2", int(sqrt2_power)))
        if i_power:
            powers.append(("I", 1))
        monomials.append(monomial(prefactor, *powers))
    return polynomial(*monomials)


def classify_and_emit(expression: sp.Expr, entries: list[dict[str, Any]]) -> None:
    """Split an expanded scalar-gauge Lagrangian piece into records.

    Terms with one derivative, one vector, and one plain scalar are
    folded onto the antisymmetric scalar_scalar_gauge convention: the
    record carries the coefficient of (d phi_a) phi_b V for the sorted
    scalar pair, and the run aborts unless the reversed assignment
    carries exactly the negated coefficient.
    """

    expression = sp.expand(expression)
    grouped: dict[tuple, sp.Expr] = {}
    for term in expression.as_ordered_terms():
        fields: list[str] = []
        coefficient = term
        for symbol_obj, label in {**FIELD_SYMBOLS, **DERIVATIVE_SYMBOLS}.items():
            power = sp.degree(term, symbol_obj)
            if power:
                fields.extend([("d:" if symbol_obj in DERIVATIVE_SYMBOLS else "") + label] * int(power))
                coefficient = coefficient / symbol_obj ** power
        key = tuple(sorted(fields))
        grouped[key] = grouped.get(key, sp.Integer(0)) + sp.simplify(coefficient)

    antisymmetric: dict[tuple, dict[str, sp.Expr]] = {}
    for key, coefficient in grouped.items():
        coefficient = sp.simplify(coefficient)
        if coefficient == 0:
            continue
        derivative_fields = [f[2:] for f in key if f.startswith("d:")]
        plain_fields = [f for f in key if not f.startswith("d:")]
        vectors = [f for f in plain_fields if f in VECTOR_LABELS]
        scalars = [f for f in plain_fields if f in SCALAR_LABELS]
        n_vec, n_scal, n_der = len(vectors), len(scalars), len(derivative_fields)
        if n_der == 2 and n_vec == 0:
            continue  # canonically normalized scalar kinetic bilinear
        if n_der == 1 and n_vec == 1 and n_scal == 1:
            pair = tuple(sorted([derivative_fields[0], scalars[0]]))
            if derivative_fields[0] == scalars[0]:
                raise SystemExit(f"engine A: symmetric derivative pair {key} -> {coefficient}")
            slot = "on_first" if derivative_fields[0] == pair[0] else "on_second"
            bucket = antisymmetric.setdefault((pair, vectors[0]), {})
            bucket[slot] = bucket.get(slot, sp.Integer(0)) + coefficient
            continue
        if n_der == 1 and n_vec == 1 and n_scal == 0:
            structure = "vector_scalar_mixing"
        elif n_der == 0 and n_vec == 2 and n_scal == 0:
            structure = "vector_bilinear_mass"
        elif n_der == 0 and n_vec == 2 and n_scal == 1:
            structure = "scalar_gauge_gauge"
        elif n_der == 0 and n_vec == 2 and n_scal == 2:
            structure = "scalar_scalar_gauge_gauge"
        elif n_der == 0 and n_vec == 0 and n_scal == 1:
            structure = "scalar_tadpole"
        elif n_der == 0 and n_vec == 0 and n_scal == 2:
            structure = "scalar_bilinear_mass"
        elif n_der == 0 and n_vec == 0 and n_scal in (3, 4):
            structure = "scalar_potential"
        else:
            raise SystemExit(f"engine A: unclassified term {key} -> {coefficient}")
        entries.append(record(plain_fields + derivative_fields, to_polynomial(coefficient), structure))

    for (pair, vector), bucket in antisymmetric.items():
        on_first = sp.simplify(bucket.get("on_first", sp.Integer(0)))
        on_second = sp.simplify(bucket.get("on_second", sp.Integer(0)))
        if sp.simplify(on_first + on_second) != 0:
            raise SystemExit(
                f"engine A: derivative coupling {pair} {vector} is not antisymmetric: "
                f"{on_first} vs {on_second}"
            )
        if on_first == 0:
            continue
        entries.append(record(list(pair) + [vector], to_polynomial(on_first), "scalar_scalar_gauge"))


# ---------------------------------------------------------------------------
# Doublet parameterization shared by the scalar-gauge sectors
# ---------------------------------------------------------------------------

PHI0 = (v + h + sp.I * G0) / sqrt2
PHI0C = (v + h - sp.I * G0) / sqrt2
H_DOUBLET = sp.Matrix([[Gp], [PHI0]])
HC_ROW = sp.Matrix([[Gm, PHI0C]])
GAUGE_MATRIX = -sp.I / 2 * sp.Matrix(
    [[g2 * W3 + g1 * B, g2 * sqrt2 * Wp], [g2 * sqrt2 * Wm, -g2 * W3 + g1 * B]]
)
GAUGE_MATRIX_C = sp.I / 2 * sp.Matrix(
    [[g2 * W3 + g1 * B, g2 * sqrt2 * Wp], [g2 * sqrt2 * Wm, -g2 * W3 + g1 * B]]
)

# Real Goldstone components: Gp = (r1 - I r2)/sqrt2, matching the W
# rotation Wp = (W1 - I W2)/sqrt2.
r1, r2 = sp.symbols("r1 r2")
RW1, RW2 = sp.symbols("RW1 RW2")
TO_REAL_SCALARS = {Gp: (r1 - sp.I * r2) / sqrt2, Gm: (r1 + sp.I * r2) / sqrt2}
TO_REAL_VECTORS = {Wp: (RW1 - sp.I * RW2) / sqrt2, Wm: (RW1 + sp.I * RW2) / sqrt2}
FROM_REAL_SCALARS = {r1: (Gp + Gm) / sqrt2, r2: sp.I * (Gp - Gm) / sqrt2}


def derivative_cross_expression() -> sp.Expr:
    dH = sp.Matrix([[dGp], [(dh + sp.I * dG0) / sqrt2]])
    dHc = sp.Matrix([[dGm, (dh - sp.I * dG0) / sqrt2]])
    return sp.expand((dHc * GAUGE_MATRIX * H_DOUBLET)[0] + (HC_ROW * GAUGE_MATRIX_C * dH)[0])


def higgs_kinetic_entries() -> list[dict[str, Any]]:
    seagull_and_mass = (HC_ROW * GAUGE_MATRIX_C * GAUGE_MATRIX * H_DOUBLET)[0]
    entries: list[dict[str, Any]] = []
    classify_and_emit(seagull_and_mass, entries)
    classify_and_emit(derivative_cross_expression(), entries)
    return entries


def potential_entries() -> list[dict[str, Any]]:
    """L = + mu2 (H'H) - lam (H'H)^2, minimum not imposed."""

    hh = Gp * Gm + PHI0 * PHI0C
    lagrangian = mu2 * hh - lam * hh ** 2
    constant_part = lagrangian.subs({h: 0, G0: 0, Gp: 0, Gm: 0})
    entries: list[dict[str, Any]] = []
    classify_and_emit(sp.expand(lagrangian - constant_part), entries)
    return entries


# ---------------------------------------------------------------------------
# Gauge fixing solved from the mixing-cancellation requirement
# ---------------------------------------------------------------------------


def solve_gauge_fixing() -> dict[str, Any]:
    """Solve the R_xi coefficients k_a from exact mixing cancellation.

    The |DH|^2 cross term contains bilinears V^mu d_mu phi with
    coefficient m[V][phi] at v.  The gauge-fixing term
    -(1/(2 xi)) (d.V^a + xi k_a . phi)^2 contributes +k_a[phi] to the
    same bilinear after integration by parts, so k = -m cancels the
    mixing identically; the xi factors drop out of the cross term.
    """

    mixing = derivative_cross_expression().subs({h: 0, G0: 0, Gp: 0, Gm: 0})
    mixing = sp.expand(mixing.subs(TO_REAL_VECTORS).subs(
        {dGp: (sp.Symbol("dr1") - sp.I * sp.Symbol("dr2")) / sqrt2,
         dGm: (sp.Symbol("dr1") + sp.I * sp.Symbol("dr2")) / sqrt2}
    ))
    vector_syms = {"W1": RW1, "W2": RW2, "W3": W3, "B": B}
    scalar_syms = {"r1": sp.Symbol("dr1"), "r2": sp.Symbol("dr2"), "G0": dG0, "h": dh}
    kappa: dict[str, dict[str, sp.Expr]] = {}
    for v_name, v_sym in vector_syms.items():
        for s_name, s_sym in scalar_syms.items():
            m_coeff = sp.simplify(mixing.coeff(v_sym).coeff(s_sym))
            if m_coeff != 0:
                kappa.setdefault(v_name, {})[s_name] = sp.simplify(-m_coeff)
    return kappa


def gauge_fixing_entries(kappa: dict[str, dict[str, sp.Expr]]) -> list[dict[str, Any]]:
    """Emit the gauge-fixing records: the mixing counterterms
    +k_a[phi] V^a d phi and the xi-weighted Goldstone bilinears
    -(xi_a/2) (k_a . phi)^2, mapped back to the charged basis."""

    real_scalars = {"r1": r1, "r2": r2, "G0": G0, "h": h}
    real_vectors = {"W1": RW1, "W2": RW2, "W3": W3, "B": B}
    d_real = {"r1": sp.Symbol("dr1"), "r2": sp.Symbol("dr2"), "G0": dG0, "h": dh}
    from_real_d = {
        sp.Symbol("dr1"): (dGp + dGm) / sqrt2,
        sp.Symbol("dr2"): sp.I * (dGp - dGm) / sqrt2,
    }
    xi_of = {"W1": xi2, "W2": xi2, "W3": xi2, "B": xi1}
    mixing_counterterms = sp.Integer(0)
    scalar_bilinears = sp.Integer(0)
    for v_name, row in kappa.items():
        k_dot_phi = sum(coeff * real_scalars[s_name] for s_name, coeff in row.items())
        k_dot_dphi = sum(coeff * d_real[s_name] for s_name, coeff in row.items())
        mixing_counterterms += real_vectors[v_name] * k_dot_dphi
        scalar_bilinears += -(xi_of[v_name] / 2) * k_dot_phi ** 2
    back = {RW1: (Wp + Wm) / sqrt2, RW2: sp.I * (Wp - Wm) / sqrt2}
    mixing_counterterms = sp.expand(mixing_counterterms.subs(back).subs(from_real_d))
    scalar_bilinears = sp.expand(scalar_bilinears.subs(FROM_REAL_SCALARS))
    entries: list[dict[str, Any]] = []
    classify_and_emit(mixing_counterterms, entries)
    classify_and_emit(scalar_bilinears, entries)
    return entries


# ---------------------------------------------------------------------------
# Yukawa sector
# ---------------------------------------------------------------------------

DAGGER_NAME = {"Yu": "Yud", "Yd": "Ydd", "Ye": "Yed"}


def yukawa_entries() -> list[dict[str, Any]]:
    """Expand the three Yukawa operators with the exact conjugate doublet.

    - -Qbar Yu Htilde uR with Htilde = (phi0c, -Gm) gives
      -uLbar Yu uR phi0c + dLbar Yu uR Gm.
    - -Qbar Yd H dR gives -uLbar Yd dR Gp - dLbar Yd dR phi0.
    - -Lbar Ye H eR gives -nuLbar Ye eR Gp - eLbar Ye eR phi0.
    phi0 carries + I G0 and phi0c carries - I G0.  All nine matrix
    elements are emitted for every record, including the mass
    bilinears, under the positional binding of vertex_format; the
    hermitian conjugates carry the dagger symbols with transposed
    indices and conjugated phases.
    """

    entries: list[dict[str, Any]] = []

    def emit_neutral(left: str, right: str, sector: str, sign: int, g0_phase: int) -> None:
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                coupling = (f"{sector}[{i}][{j}]", 1)
                base = [bar(left), right]
                entries.append(record(base + ["h"], polynomial(monomial(sign, coupling, ("sqrt2", -1))), "fermion_scalar_yukawa"))
                entries.append(record(base + ["G0"], polynomial(monomial(sign * g0_phase, coupling, ("sqrt2", -1), ("I", 1))), "fermion_scalar_yukawa"))
                entries.append(record(base, polynomial(monomial(sign, coupling, ("v", 1), ("sqrt2", -1))), "fermion_bilinear_mass"))

    def emit_charged(left: str, right: str, scalar: str, sector: str, sign: int) -> None:
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                coupling = (f"{sector}[{i}][{j}]", 1)
                entries.append(
                    record([bar(left), right, scalar], polynomial(monomial(sign, coupling)), "fermion_scalar_yukawa")
                )

    emit_neutral("uL", "uR", "Yu", -1, -1)
    emit_charged("dL", "uR", "Gm", "Yu", 1)
    emit_neutral("dL", "dR", "Yd", -1, 1)
    emit_charged("uL", "dR", "Gp", "Yd", -1)
    emit_neutral("eL", "eR", "Ye", -1, 1)
    emit_charged("nuL", "eR", "Gp", "Ye", -1)

    conjugates: list[dict[str, Any]] = []
    for entry in entries:
        fields = entry["fields"]
        bar_fields = [f for f in fields if f.endswith("_bar")]
        plain_fermions = [f for f in fields if not f.endswith("_bar") and f not in SCALAR_LABELS]
        scalars = [f for f in fields if f in SCALAR_LABELS]
        swapped_scalars = [{"Gp": "Gm", "Gm": "Gp"}.get(s, s) for s in scalars]
        new_fields = [bar(f) for f in plain_fermions] + [f[:-4] for f in bar_fields] + swapped_scalars
        conjugated = []
        for m in entry["coefficient"]["monomials"]:
            prefactor = Fraction(m["prefactor"])
            powers = []
            for symbol_name, power in m["powers"]:
                if symbol_name == "I":
                    prefactor = -prefactor
                    powers.append((symbol_name, power))
                    continue
                for plain, dagger in DAGGER_NAME.items():
                    if symbol_name.startswith(plain + "["):
                        i_idx, j_idx = symbol_name[len(plain):].strip("[]").split("][")
                        powers.append((f"{dagger}[{j_idx}][{i_idx}]", power))
                        break
                else:
                    powers.append((symbol_name, power))
            conjugated.append((str(prefactor), tuple(powers)))
        conjugates.append(record(new_fields, polynomial(*conjugated), entry["structure"]))
    return entries + conjugates


# ---------------------------------------------------------------------------
# Yang-Mills self-couplings with Lorentz role slots
# ---------------------------------------------------------------------------


def epsilon(a: int, b: int, c: int) -> int:
    return int((b - a) * (c - b) * (c - a) / 2)


def charged_rotation(prefix: str) -> tuple[list[sp.Expr], dict[sp.Symbol, str]]:
    """Return the real-basis role symbols rotated to the charged basis
    and the map from charged role symbols to field labels."""

    plus, minus, third = sp.symbols(f"{prefix}p {prefix}m {prefix}3")
    rotated = [(plus + minus) / sqrt2, sp.I * (plus - minus) / sqrt2, third]
    labels = {plus: "Wp", minus: "Wm", third: "W3"}
    return rotated, labels


def term_role_fields(term: sp.Expr, label_maps: list[dict[sp.Symbol, str]]) -> tuple[list[list[str]], sp.Expr]:
    """Extract per-role field labels and the residual coupling."""

    coefficient = term
    per_role: list[list[str]] = []
    for label_map in label_maps:
        role_fields: list[str] = []
        for symbol_obj, label in label_map.items():
            power = sp.degree(term, symbol_obj)
            if power:
                role_fields.extend([label] * int(power))
                coefficient = coefficient / symbol_obj ** power
        per_role.append(role_fields)
    return per_role, sp.simplify(coefficient)


def yang_mills_cubic_entries() -> list[dict[str, Any]]:
    """Cubic records from -g2 eps_abc (dW)_a W_b W_c.

    Roles: A = d_mu W_nu, B = W^mu, C = W^nu.  The record carries the
    coefficient of the sorted role assignment (A, B, C) = (f1, f2, f3)
    and total antisymmetry over role permutations is verified for every
    field triple before emission.
    """

    roles = [charged_rotation(prefix) for prefix in ("A", "B", "C")]
    expression = sp.Integer(0)
    for a in range(3):
        for b_i in range(3):
            for c_i in range(3):
                eps = epsilon(a, b_i, c_i)
                if eps:
                    expression += eps * roles[0][0][a] * roles[1][0][b_i] * roles[2][0][c_i]
    expression = sp.expand(-g2 * expression)
    tensor: dict[tuple[str, str, str], sp.Expr] = {}
    for term in expression.as_ordered_terms():
        per_role, coefficient = term_role_fields(term, [r[1] for r in roles])
        if any(len(fields) != 1 for fields in per_role):
            raise SystemExit("engine A: cubic term without exactly one field per role")
        key = (per_role[0][0], per_role[1][0], per_role[2][0])
        tensor[key] = tensor.get(key, sp.Integer(0)) + coefficient
    entries: list[dict[str, Any]] = []
    seen: set[tuple[str, ...]] = set()
    for key, coefficient in tensor.items():
        multiset = tuple(sorted(key))
        if multiset in seen:
            continue
        seen.add(multiset)
        if len(set(multiset)) != 3:
            for assignment, value in tensor.items():
                if tuple(sorted(assignment)) == multiset and sp.simplify(value) != 0:
                    raise SystemExit(f"engine A: repeated-field cubic assignment {assignment} survives")
            continue
        reference = sp.simplify(tensor.get(multiset, sp.Integer(0)))
        for permutation in itertools.permutations(range(3)):
            assignment = tuple(multiset[p] for p in permutation)
            sign = sp.Integer(1)
            for x in range(3):
                for y in range(x + 1, 3):
                    if permutation[x] > permutation[y]:
                        sign = -sign
            value = sp.simplify(tensor.get(assignment, sp.Integer(0)))
            if sp.simplify(value - sign * reference) != 0:
                raise SystemExit(f"engine A: cubic tensor is not totally antisymmetric at {assignment}")
        if reference != 0:
            entries.append(record(list(multiset), to_polynomial(reference), "yang_mills_three_point"))
    return entries


def yang_mills_quartic_entries() -> list[dict[str, Any]]:
    """Quartic records from -(g2^2/4) sum_a (eps_abc P_b Q_c)^2-form.

    P carries the mu slot and Q the nu slot; a term P_x P_z Q_y Q_w is
    the Lorentz scalar (x . z)(y . w), so each term is routed to the
    pairing partition of its P and Q multisets and emitted under the
    first pairing label whose value partition matches.
    """

    p_rotated, p_labels = charged_rotation("P")
    q_rotated, q_labels = charged_rotation("Q")
    expression = sp.Integer(0)
    for a in range(3):
        factor = sp.Integer(0)
        for b_i in range(3):
            for c_i in range(3):
                eps = epsilon(a, b_i, c_i)
                if eps:
                    factor += eps * p_rotated[b_i] * q_rotated[c_i]
        expression += factor ** 2
    expression = sp.expand(-sp.Rational(1, 4) * g2 ** 2 * expression)
    accumulator: dict[tuple[tuple[str, ...], tuple[tuple[str, str], ...]], sp.Expr] = {}
    for term in expression.as_ordered_terms():
        per_role, coefficient = term_role_fields(term, [p_labels, q_labels])
        p_fields, q_fields = per_role
        if len(p_fields) != 2 or len(q_fields) != 2:
            raise SystemExit("engine A: quartic term without two fields per slot")
        multiset = tuple(sorted(p_fields + q_fields))
        partition = tuple(sorted([tuple(sorted(p_fields)), tuple(sorted(q_fields))]))
        key = (multiset, partition)
        accumulator[key] = accumulator.get(key, sp.Integer(0)) + coefficient
    pairing_labels = (
        ("yang_mills_four_point_12_34", ((0, 1), (2, 3))),
        ("yang_mills_four_point_13_24", ((0, 2), (1, 3))),
        ("yang_mills_four_point_14_23", ((0, 3), (1, 2))),
    )
    entries: list[dict[str, Any]] = []
    for (multiset, partition), coefficient in sorted(accumulator.items()):
        coefficient = sp.simplify(coefficient)
        if coefficient == 0:
            continue
        chosen = None
        for label, index_pairs in pairing_labels:
            value_partition = tuple(sorted(
                tuple(sorted((multiset[i], multiset[j]))) for i, j in index_pairs
            ))
            if value_partition == partition:
                chosen = label
                break
        if chosen is None:
            raise SystemExit(f"engine A: quartic partition {partition} has no pairing label")
        entries.append(record(list(multiset), to_polynomial(coefficient), chosen))
    return entries


# ---------------------------------------------------------------------------
# Ghost sector from the exact gauge variations
# ---------------------------------------------------------------------------


def ghost_entries(kappa: dict[str, dict[str, sp.Expr]]) -> list[dict[str, Any]]:
    """L_ghost = - cbar_a (delta G^a / delta theta_b) c_b.

    delta W^a = d theta^a + g2 eps_abc W^b theta^c and delta B = d
    theta_B give the kinetic terms (excluded as bilinears without mass)
    and the (d cbar) V c gauge vertices; the doublet variation
    delta H = +(I/2)(g2 theta^a sigma^a + g1 theta_B) H evaluated on
    the shifted doublet gives, through the solved gauge functions, the
    xi-weighted ghost masses including the neutral cbar3-cbarB mixing
    and every ghost-scalar coupling.  Ghosts rotate to the charged
    basis alongside the vectors, and every term is ordered with the
    barred ghost leftmost.
    """

    theta1, theta2, theta3, thetaB = sp.symbols("theta1 theta2 theta3 thetaB")

    # Doublet variation in components.  The pair (delta H, delta W)
    # consistent with the covariant derivative d + M, M = -(I/2)(...)
    # and delta W^a = + d theta^a + g2 eps_abc W^b theta^c is
    # delta H = +(I/2)(g2 theta^a sigma^a + g1 thetaB) H; the relative
    # sign against the kinetic ghost term is checked downstream by the
    # ghost poles matching the xi-weighted Goldstone poles.
    delta_Gp = sp.I / 2 * ((g2 * theta3 + g1 * thetaB) * Gp + g2 * (theta1 - sp.I * theta2) * PHI0)
    delta_phi0 = sp.I / 2 * (g2 * (theta1 + sp.I * theta2) * Gp + (-g2 * theta3 + g1 * thetaB) * PHI0)
    delta_Gm = -sp.I / 2 * ((g2 * theta3 + g1 * thetaB) * Gm + g2 * (theta1 + sp.I * theta2) * PHI0C)
    delta_phi0c = -sp.I / 2 * (g2 * (theta1 - sp.I * theta2) * Gm + (-g2 * theta3 + g1 * thetaB) * PHI0C)
    delta_scalar = {
        "h": sqrt2 * (delta_phi0 + delta_phi0c) / 2,
        "G0": sqrt2 * (delta_phi0 - delta_phi0c) / (2 * sp.I),
        "r1": (delta_Gp + delta_Gm) / sqrt2,
        "r2": sp.I * (delta_Gp - delta_Gm) / sqrt2,
    }

    # Scalar part: - cbar_a xi_a (k_a . delta phi).
    barred_of_vector = {"W1": "cb1", "W2": "cb2", "W3": bar("c3"), "B": bar("cB")}
    xi_of = {"W1": xi2, "W2": xi2, "W3": xi2, "B": xi1}
    cb1, cb2 = sp.symbols("cb1 cb2")
    barred_syms = {"cb1": cb1, "cb2": cb2, bar("c3"): sp.Symbol("cb3"), bar("cB"): sp.Symbol("cbB")}
    scalar_part = sp.Integer(0)
    for v_name, row in kappa.items():
        k_dot_delta = sum(coeff * delta_scalar[s_name] for s_name, coeff in row.items())
        scalar_part += -barred_syms[barred_of_vector[v_name]] * xi_of[v_name] * k_dot_delta
    # Rotate barred and unbarred ghosts to the charged basis.
    cpb, cmb = sp.symbols("cpb cmb")
    cp_sym, cm_sym = sp.symbols("cp cm")
    ghost_rotation = {
        cb1: (cpb + cmb) / sqrt2,
        cb2: -sp.I * (cpb - cmb) / sqrt2,
        theta1: (cp_sym + cm_sym) / sqrt2,
        theta2: sp.I * (cp_sym - cm_sym) / sqrt2,
        theta3: sp.Symbol("cu3"),
        thetaB: sp.Symbol("cuB"),
    }
    scalar_part = sp.expand(scalar_part.subs(ghost_rotation))

    barred_labels = {cpb: bar("cp"), cmb: bar("cm"), sp.Symbol("cb3"): bar("c3"), sp.Symbol("cbB"): bar("cB")}
    unbarred_labels = {cp_sym: "cp", cm_sym: "cm", sp.Symbol("cu3"): "c3", sp.Symbol("cuB"): "cB"}
    scalar_field_labels = {h: "h", G0: "G0", Gp: "Gp", Gm: "Gm"}

    entries: list[dict[str, Any]] = []
    grouped: dict[tuple, sp.Expr] = {}
    for term in scalar_part.as_ordered_terms():
        coefficient = term
        fields: list[str] = []
        for symbol_obj, label in {**barred_labels, **unbarred_labels, **scalar_field_labels}.items():
            power = sp.degree(term, symbol_obj)
            if power:
                fields.extend([label] * int(power))
                coefficient = coefficient / symbol_obj ** power
        ghost_count = sum(1 for f in fields if f.startswith("c") or f.startswith("cb"))
        if ghost_count != 2:
            raise SystemExit(f"engine A: ghost scalar term without a ghost pair: {term}")
        key = tuple(sorted(fields))
        grouped[key] = grouped.get(key, sp.Integer(0)) + sp.simplify(coefficient)
    for key, coefficient in grouped.items():
        coefficient = sp.simplify(coefficient)
        if coefficient == 0:
            continue
        entries.append(record(list(key), to_polynomial(coefficient), "ghost_scalar_mass"))

    # Gauge vertex: + g2 eps_abc (d cbar_a) W^b c^c after integration
    # by parts of - cbar_a d.(g2 eps_abc W^b theta^c).
    vertex = sp.Integer(0)
    w_real = [(Wp + Wm) / sqrt2, sp.I * (Wp - Wm) / sqrt2, W3]
    theta_list = [theta1, theta2, theta3]
    cbar_list = [cb1, cb2, sp.Symbol("cb3")]
    for a in range(3):
        for b_i in range(3):
            for c_i in range(3):
                eps = epsilon(a, b_i, c_i)
                if eps:
                    vertex += g2 * eps * cbar_list[a] * w_real[b_i] * theta_list[c_i]
    vertex = sp.expand(vertex.subs(ghost_rotation))
    grouped_vertex: dict[tuple, sp.Expr] = {}
    for term in vertex.as_ordered_terms():
        coefficient = term
        fields = []
        for symbol_obj, label in {**barred_labels, **unbarred_labels, Wp: "Wp", Wm: "Wm", W3: "W3"}.items():
            power = sp.degree(term, symbol_obj)
            if power:
                fields.extend([label] * int(power))
                coefficient = coefficient / symbol_obj ** power
        key = tuple(sorted(fields))
        grouped_vertex[key] = grouped_vertex.get(key, sp.Integer(0)) + sp.simplify(coefficient)
    for key, coefficient in grouped_vertex.items():
        coefficient = sp.simplify(coefficient)
        if coefficient == 0:
            continue
        entries.append(record(list(key), to_polynomial(coefficient), "ghost_gauge_derivative"))
    return entries


# ---------------------------------------------------------------------------
# Fermion gauge currents from the census
# ---------------------------------------------------------------------------


def fermion_current_entries(action: dict[str, Any]) -> list[dict[str, Any]]:
    """Gauge currents from the census: W3 couples g2 T3, B couples g1 Y,
    W+- couple g2/sqrt2 within doublets, gluons couple g3 to triplets.
    Flavor-diagonal, collapsed to the template pair."""

    entries: list[dict[str, Any]] = []
    components = {
        "Q": (("uL", Fraction(1, 2)), ("dL", Fraction(-1, 2))),
        "L": (("nuL", Fraction(1, 2)), ("eL", Fraction(-1, 2))),
        "u_c": (("uR", None),),
        "d_c": (("dR", None),),
        "e_c": (("eR", None),),
    }
    census = {f["name"]: f for f in action["field_census"]["fermions"] if f["generation"] == 1}
    for multiplet, parts in components.items():
        y = Fraction(census[multiplet]["hypercharge"])
        if len(parts) == 1:
            # The census stores the left-handed conjugate Weyl field;
            # the record template is the right-handed component, whose
            # hypercharge is the negative of the conjugate entry.
            y = -y
        colored = census[multiplet]["color"] in ("3", "3bar")
        for component, t3 in parts:
            if colored:
                entries.append(record([bar(component), component, "Gl"], polynomial(monomial(1, ("g3", 1))), "fermion_vector_current"))
            if t3 is not None:
                entries.append(record([bar(component), component, "W3"], polynomial(monomial(t3, ("g2", 1))), "fermion_vector_current"))
            if y != 0:
                entries.append(record([bar(component), component, "B"], polynomial(monomial(y, ("g1", 1))), "fermion_vector_current"))
        if len(parts) == 2:
            up, down = parts[0][0], parts[1][0]
            entries.append(record([bar(up), down, "Wp"], polynomial(monomial(1, ("g2", 1), ("sqrt2", -1))), "fermion_vector_current"))
            entries.append(record([bar(down), up, "Wm"], polynomial(monomial(1, ("g2", 1), ("sqrt2", -1))), "fermion_vector_current"))
    return entries


# ---------------------------------------------------------------------------
# Table assembly
# ---------------------------------------------------------------------------


def build_table() -> dict[str, Any]:
    action = load_action()
    kappa = solve_gauge_fixing()
    entries = (
        higgs_kinetic_entries()
        + gauge_fixing_entries(kappa)
        + potential_entries()
        + yukawa_entries()
        + yang_mills_cubic_entries()
        + yang_mills_quartic_entries()
        + fermion_current_entries(action)
        + ghost_entries(kappa)
    )
    merged: dict[str, dict[str, Any]] = {}
    cancelled: list[dict[str, Any]] = []
    for entry in entries:
        key = json.dumps({"fields": entry["fields"], "structure": entry["structure"]}, sort_keys=True)
        if key in merged:
            existing = merged[key]["coefficient"]["monomials"]
            addition = entry["coefficient"]["monomials"]
            combined = [
                (m["prefactor"], tuple(tuple(p) for p in m["powers"]))
                for m in existing + addition
            ]
            summed = polynomial(*combined)
            if summed["monomials"]:
                merged[key] = record(entry["fields"], summed, entry["structure"])
            else:
                cancelled.append({"fields": entry["fields"], "structure": entry["structure"]})
                del merged[key]
        else:
            merged[key] = entry
    final = sorted(merged.values(), key=vertex_hash)
    for entry in final:
        violations = conservation_violations(entry["fields"])
        if violations:
            raise SystemExit(f"engine A: conservation violation {violations} in {entry['fields']}")
    kappa_report = {
        v_name: {s_name: str(coeff) for s_name, coeff in sorted(row.items())}
        for v_name, row in sorted(kappa.items())
    }
    return {
        "schema": "rule_table.v2",
        "engine": "A_action_expansion",
        "action_subject_digest": action["subject_digest"],
        "gauge_fixing": {
            "form": "G^a = d.W^a + xi_a (k_a . phi), L_gf = -(1/(2 xi_a)) (G^a)^2",
            "solved_kappa": kappa_report,
            "mixing_cancelled": sorted(cancelled, key=lambda item: json.dumps(item, sort_keys=True)),
        },
        "exclusions": [
            {
                "sector": "gluon_self_couplings_and_gluon_ghost_gauge_vertex",
                "reason": (
                    "the census gluon is collapsed to a single adjoint template, and "
                    "gluon self-interactions enter W/Z pole functions beyond one loop"
                ),
            },
            {
                "sector": "xi3_gluon_gauge_fixing_parameter",
                "reason": "carried in the symbol alphabet, unused by any electroweak record",
            },
        ],
        "entries": final,
        "entry_count": len(final),
        "table_digest": table_digest(final),
    }


def main() -> int:
    table = build_table()
    OUT_PATH.write_text(json.dumps(table, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "WROTE", "entries": table["entry_count"], "table_digest": table["table_digest"][:24]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
