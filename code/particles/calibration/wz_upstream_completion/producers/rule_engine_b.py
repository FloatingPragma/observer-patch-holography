#!/usr/bin/env python3
"""Workstream C, engine B: vertex records by generator-algebra assignment.

Independent second derivation of the canonical rule table.  Engine B
shares with engine A only the canonicalization contract of
vertex_format and the pinned action census; the derivation route and
the arithmetic substrate are disjoint:

* No computer-algebra system: every coefficient is assembled in exact
  complex-rational arithmetic (pairs of Fractions) over an in-module
  polynomial container keyed by field multisets and coupling powers.
* The scalar-gauge sector is computed from the real SO(4) generator
  action on the four real scalar components, with the generator
  matrices themselves computed from the complex sigma-matrix action
  by the module's own arithmetic, and per-structure closed forms:
  vector masses (1/2) g g' (T Phi0).(T' Phi0), gauge seagulls from
  (T Phi).(T' Phi), derivative couplings from dPhi.(T Phi), the R_xi
  functions k_a = - g_a T^a Phi0, and the ghost blocks from
  delta Phi = - g_b theta_b T^b Phi.
* The Yang-Mills cubic and quartic are assembled from the structure
  constants with the module's own charged-basis rotation and Lorentz
  pairing bookkeeping.
* The Yukawa sector is assembled from the conjugate doublet built by
  explicit epsilon contraction of the conjugated components.
* Fermion gauge currents follow the closed rules g2 T3, g1 Y, g2/sqrt2
  and g3 applied to the census multiplets.

The output must reproduce the engine A table digest exactly; the
equivalence checker enforces that identity.
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

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
OUT_PATH = ROOT / "outputs" / "rule_table_engine_b.json"

# ---------------------------------------------------------------------------
# Exact complex-rational polynomial container
# ---------------------------------------------------------------------------
# A CPoly maps (fields, couplings) -> (re, im) with fields a sorted
# tuple of component labels ("d:x" marks a derivative factor) and
# couplings a sorted tuple of (symbol, power) pairs; sqrt2 is carried
# as an explicit power and folded by vertex_format at emission.

ZERO = Fraction(0)


def cnum(re: Fraction | int, im: Fraction | int = 0) -> tuple[Fraction, Fraction]:
    return (Fraction(re), Fraction(im))


def cadd(a: tuple, b: tuple) -> tuple:
    return (a[0] + b[0], a[1] + b[1])


def cmul(a: tuple, b: tuple) -> tuple:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cconj(a: tuple) -> tuple:
    return (a[0], -a[1])


def poly(fields: tuple = (), couplings: tuple = (), value: tuple = (Fraction(1), Fraction(0))) -> dict:
    if value == (ZERO, ZERO):
        return {}
    return {(tuple(sorted(fields)), tuple(sorted(couplings))): value}


def padd(*polys: dict) -> dict:
    out: dict = {}
    for p in polys:
        for key, value in p.items():
            total = cadd(out.get(key, (ZERO, ZERO)), value)
            if total == (ZERO, ZERO):
                out.pop(key, None)
            else:
                out[key] = total
    return out


def pscale(p: dict, value: tuple) -> dict:
    out: dict = {}
    if value == (ZERO, ZERO):
        return out
    for key, existing in p.items():
        product = cmul(existing, value)
        if product != (ZERO, ZERO):
            out[key] = product
    return out


def merge_couplings(a: tuple, b: tuple) -> tuple:
    combined: dict[str, int] = {}
    for name, power in list(a) + list(b):
        combined[name] = combined.get(name, 0) + power
    return tuple(sorted((n, p) for n, p in combined.items() if p))


def pmul(a: dict, b: dict) -> dict:
    out: dict = {}
    for (fa, ca), va in a.items():
        for (fb, cb), vb in b.items():
            key = (tuple(sorted(fa + fb)), merge_couplings(ca, cb))
            total = cadd(out.get(key, (ZERO, ZERO)), cmul(va, vb))
            if total == (ZERO, ZERO):
                out.pop(key, None)
            else:
                out[key] = total
    return out


def ppow(p: dict, exponent: int) -> dict:
    out = poly()
    for _ in range(exponent):
        out = pmul(out, p)
    return out


def psubst(p: dict, field: str, replacement: dict) -> dict:
    out: dict = {}
    for (fields, couplings), value in p.items():
        count = sum(1 for f in fields if f == field)
        base = poly(tuple(f for f in fields if f != field), couplings, value)
        term = base
        for _ in range(count):
            term = pmul(term, replacement)
        out = padd(out, term)
    return out


SQRT2_INV = poly(couplings=(("sqrt2", -1),))
SQRT2_POS = poly(couplings=(("sqrt2", 1),))


def coupling(name: str, power: int = 1) -> dict:
    return poly(couplings=((name, power),))


def field(name: str) -> dict:
    return poly(fields=(name,))


# ---------------------------------------------------------------------------
# Complex doublet components and the real generator algebra
# ---------------------------------------------------------------------------
# Complex basis: H = (H1, H2) with H1 = (r1 - I r2)/sqrt2 and
# H2 = (vh + I G0)/sqrt2, vh = v + h.  Real basis Phi = (r1, r2, vh, G0)
# with Phi0 = (0, 0, v, 0).

REAL_BASIS = ("r1", "r2", "hshift", "G0")

# Complex components of H over the real basis: H_i = sum_j C[i][j] phi_j
H_COMPONENTS = (
    {"r1": cnum(1), "r2": cnum(0, -1)},
    {"hshift": cnum(1), "G0": cnum(0, 1)},
)

# sigma-action factors: delta H = -(I/2) K H per unit coupling with
# K in {sigma1, sigma2, sigma3, identity} for (W1, W2, W3, B).
SIGMA = {
    "W1": ((cnum(0), cnum(1)), (cnum(1), cnum(0))),
    "W2": ((cnum(0), cnum(0, -1)), (cnum(0, 1), cnum(0))),
    "W3": ((cnum(1), cnum(0)), (cnum(0), cnum(-1))),
    "B": ((cnum(1), cnum(0)), (cnum(0), cnum(1))),
}
GAUGE_COUPLING = {"W1": "g2", "W2": "g2", "W3": "g2", "B": "g1"}
XI_OF = {"W1": "xi2", "W2": "xi2", "W3": "xi2", "B": "xi1"}


def real_generator(vector: str) -> list[list[Fraction]]:
    """Real 4x4 matrix T with (M Phi)_j = sum_k T[j][k] phi_k for the
    action M = -(I/2) K on the complex components, computed by the
    module's own complex arithmetic.

    The complex image of basis vector e_k under -(I/2) K determines one
    column: expand (M H)_i = sum_j C[i][j] (T phi)_j and match real and
    imaginary parts against the complex component map.
    """

    K = SIGMA[vector]
    minus_i_half = cnum(0, Fraction(-1, 2))
    columns: dict[str, dict[str, Fraction]] = {}
    for k, source in enumerate(REAL_BASIS):
        image: dict[str, Fraction] = {name: ZERO for name in REAL_BASIS}
        for i in (0, 1):
            amplitude = H_COMPONENTS[i].get(source)
            if amplitude is None:
                continue
            for target_i in (0, 1):
                factor = cmul(minus_i_half, cmul(K[target_i][i], amplitude))
                # Decompose the image over the target complex component:
                # H_target = sum_j C[target][j] phi_j with |C| = 1
                # entries, so phi_j receives Re(factor / C[target][j])
                # via the orthonormal real decomposition.
                for j_name, c_val in H_COMPONENTS[target_i].items():
                    # coefficient of phi_j: Re(conj(C) * factor) since
                    # the complex components are unitary in each slot.
                    contribution = cmul(cconj(c_val), factor)
                    image[j_name] += contribution[0]
        columns[source] = image
    return [[columns[src][dst] for src in REAL_BASIS] for dst in REAL_BASIS]


GENERATORS = {name: real_generator(name) for name in SIGMA}

PHI_POLYS = {
    "r1": field("r1"),
    "r2": field("r2"),
    "hshift": padd(coupling("v"), field("h")),
    "G0": field("G0"),
}
PHI0 = {"r1": poly() and {}, "r2": {}, "hshift": coupling("v"), "G0": {}}
D_PHI = {"r1": field("d:r1"), "r2": field("d:r2"), "hshift": field("d:h"), "G0": field("d:G0")}


def apply_generator(vector: str, components: dict[str, dict]) -> dict[str, dict]:
    T = GENERATORS[vector]
    out: dict[str, dict] = {}
    for j, target in enumerate(REAL_BASIS):
        acc: dict = {}
        for k, source in enumerate(REAL_BASIS):
            if T[j][k] != 0:
                acc = padd(acc, pscale(components[source], cnum(T[j][k])))
        out[target] = acc
    return out


def dot(a: dict[str, dict], b: dict[str, dict]) -> dict:
    out: dict = {}
    for name in REAL_BASIS:
        out = padd(out, pmul(a[name], b[name]))
    return out


# Charged-basis substitutions applied at emission time.
def to_charged(p: dict) -> dict:
    p = psubst(p, "r1", pmul(padd(field("Gp"), field("Gm")), SQRT2_INV))
    p = psubst(p, "r2", pscale(pmul(padd(field("Gp"), pscale(field("Gm"), cnum(-1))), SQRT2_INV), cnum(0, 1)))
    p = psubst(p, "d:r1", pmul(padd(field("d:Gp"), field("d:Gm")), SQRT2_INV))
    p = psubst(p, "d:r2", pscale(pmul(padd(field("d:Gp"), pscale(field("d:Gm"), cnum(-1))), SQRT2_INV), cnum(0, 1)))
    p = psubst(p, "W1", pmul(padd(field("Wp"), field("Wm")), SQRT2_INV))
    p = psubst(p, "W2", pscale(pmul(padd(field("Wp"), pscale(field("Wm"), cnum(-1))), SQRT2_INV), cnum(0, 1)))
    return p


# ---------------------------------------------------------------------------
# Emission: canonical records from a CPoly Lagrangian piece
# ---------------------------------------------------------------------------

VECTOR_LABELS = ("Wp", "Wm", "W3", "B", "Gl")
SCALAR_LABELS = ("h", "G0", "Gp", "Gm")


def coupling_monomials(couplings: tuple, value: tuple) -> list:
    out = []
    if value[0] != 0:
        out.append(monomial(value[0], *couplings))
    if value[1] != 0:
        out.append(monomial(value[1], *couplings, ("I", 1)))
    return out


def emit_scalar_sector(p: dict, entries: list[dict[str, Any]]) -> None:
    """Classify a charged-basis CPoly into records with the same
    conventions declared in vertex_format, implemented independently."""

    grouped: dict[tuple, list] = {}
    for (fields, couplings), value in p.items():
        grouped.setdefault(fields, []).append((couplings, value))

    antisymmetric: dict[tuple, dict[str, list]] = {}
    for fields, coefficient in grouped.items():
        derivative_fields = [f[2:] for f in fields if f.startswith("d:")]
        plain_fields = [f for f in fields if not f.startswith("d:")]
        vectors = [f for f in plain_fields if f in VECTOR_LABELS]
        scalars = [f for f in plain_fields if f in SCALAR_LABELS]
        n_vec, n_scal, n_der = len(vectors), len(scalars), len(derivative_fields)
        if n_der == 2 and n_vec == 0:
            continue
        if n_der == 1 and n_vec == 1 and n_scal == 1:
            pair = tuple(sorted([derivative_fields[0], scalars[0]]))
            if derivative_fields[0] == scalars[0]:
                raise SystemExit(f"engine B: symmetric derivative pair {fields}")
            slot = "on_first" if derivative_fields[0] == pair[0] else "on_second"
            antisymmetric.setdefault((pair, vectors[0]), {}).setdefault(slot, []).extend(coefficient)
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
            raise SystemExit(f"engine B: unclassified term {fields}")
        monomials = []
        for couplings, value in coefficient:
            monomials.extend(coupling_monomials(couplings, value))
        combined = polynomial(*monomials)
        if combined["monomials"]:
            entries.append(record(plain_fields + derivative_fields, combined, structure))

    for (pair, vector), buckets in antisymmetric.items():
        first = polynomial(*[m for c, val in buckets.get("on_first", []) for m in coupling_monomials(c, val)])
        second = polynomial(*[m for c, val in buckets.get("on_second", []) for m in coupling_monomials(c, val)])
        negated_second = polynomial(*[
            (str(-Fraction(m["prefactor"])), tuple(tuple(p) for p in m["powers"]))
            for m in second["monomials"]
        ])
        if first != negated_second:
            raise SystemExit(f"engine B: derivative coupling {pair} {vector} is not antisymmetric")
        if first["monomials"]:
            entries.append(record(list(pair) + [vector], first, "scalar_scalar_gauge"))


# ---------------------------------------------------------------------------
# Sector assignments
# ---------------------------------------------------------------------------


def scalar_gauge_entries() -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Vector masses, seagulls, derivative couplings, solved gauge
    fixing, and the xi-weighted Goldstone bilinears from the closed
    generator forms; returns the entries, the solved k vectors, and
    the derived cancelled-mixing cross-check list."""

    entries: list[dict[str, Any]] = []
    phi = {name: PHI_POLYS[name] for name in REAL_BASIS}
    gauged = {
        name: {
            target: pscale(component, cnum(1))
            for target, component in apply_generator(name, phi).items()
        }
        for name in SIGMA
    }
    # (1/2) |sum_a g_a W^a T^a Phi|^2: masses and seagulls.
    seagull: dict = {}
    for a_name, b_name in itertools.product(SIGMA, repeat=2):
        block = dot(gauged[a_name], gauged[b_name])
        block = pmul(block, coupling(GAUGE_COUPLING[a_name]))
        block = pmul(block, coupling(GAUGE_COUPLING[b_name]))
        block = pmul(block, pmul(field(a_name), field(b_name)))
        seagull = padd(seagull, pscale(block, cnum(Fraction(1, 2))))
    # dPhi . (g_a W^a T^a Phi): derivative couplings and mixing.
    cross: dict = {}
    for a_name in SIGMA:
        block = dot(D_PHI, gauged[a_name])
        block = pmul(block, coupling(GAUGE_COUPLING[a_name]))
        cross = padd(cross, pmul(block, field(a_name)))
    # Solved gauge fixing: k_a = - g_a T^a Phi0; the mixing counterterm
    # is + (k_a . dPhi) W^a and the scalar block is -(xi_a/2)(k_a.Phi)^2.
    kappa: dict[str, dict[str, Any]] = {}
    counterterm: dict = {}
    xi_bilinears: dict = {}
    for a_name in SIGMA:
        k_vector = {
            target: pscale(pmul(component, coupling(GAUGE_COUPLING[a_name])), cnum(-1))
            for target, component in apply_generator(a_name, {n: PHI0[n] for n in REAL_BASIS}).items()
        }
        kappa[a_name] = k_vector
        counterterm = padd(counterterm, pmul(dot(k_vector, D_PHI), field(a_name)))
        k_dot_phi = dot(k_vector, phi)
        # Subtract the vacuum constant before squaring so only field
        # monomials survive; k.Phi is linear in fields with no constant
        # because T^a Phi0 has no hshift component.
        xi_bilinears = padd(
            xi_bilinears,
            pscale(pmul(coupling(XI_OF[a_name]), ppow(k_dot_phi, 2)), cnum(Fraction(-1, 2))),
        )
    # Audit trail: the mixing bilinears present in the cross term must
    # each be cancelled exactly by the counterterm; the surviving total
    # must emit no vector_scalar_mixing record.
    cancelled: list[dict[str, Any]] = []
    charged_cross = to_charged(cross)
    for (fields, _couplings), _value in sorted(charged_cross.items()):
        derivative_fields = [f[2:] for f in fields if f.startswith("d:")]
        plain = [f for f in fields if not f.startswith("d:")]
        if len(fields) == 2 and len(derivative_fields) == 1 and plain and plain[0] in VECTOR_LABELS:
            item = {"fields": sorted([derivative_fields[0], plain[0]]), "structure": "vector_scalar_mixing"}
            if item not in cancelled:
                cancelled.append(item)
    cancelled.sort(key=lambda item: json.dumps(item, sort_keys=True))
    total = padd(seagull, cross, counterterm, xi_bilinears)
    emit_scalar_sector(to_charged(total), entries)
    if any(entry["structure"] == "vector_scalar_mixing" for entry in entries):
        raise SystemExit("engine B: mixing bilinear survives the solved gauge fixing")
    return entries, kappa, cancelled


def potential_entries() -> list[dict[str, Any]]:
    """L = +(mu2/2)|Phi|^2 - (lam/4)|Phi|^4 on the shifted real basis,
    constant part subtracted, minimum not imposed."""

    phi_sq: dict = {}
    for name in REAL_BASIS:
        phi_sq = padd(phi_sq, ppow(PHI_POLYS[name], 2))
    lagrangian = padd(
        pscale(pmul(coupling("mu2"), phi_sq), cnum(Fraction(1, 2))),
        pscale(pmul(coupling("lam"), ppow(phi_sq, 2)), cnum(Fraction(-1, 4))),
    )
    lagrangian = {
        key: value
        for key, value in lagrangian.items()
        if key[0]  # drop the field-free vacuum constant
    }
    entries: list[dict[str, Any]] = []
    emit_scalar_sector(to_charged(lagrangian), entries)
    return entries


def yukawa_entries() -> list[dict[str, Any]]:
    """Records from -Qbar Y Hc fR with Hc the epsilon-contracted
    conjugate doublet for the up sector and H itself for the down and
    lepton sectors, assembled in complex-rational arithmetic."""

    h1 = pmul(padd(field("Gp"), poly()), poly())  # H1 = Gp
    h1 = field("Gp")
    h2 = pmul(padd(coupling("v"), field("h"), pscale(field("G0"), cnum(0, 1))), SQRT2_INV)
    h1_conj = field("Gm")
    h2_conj = pmul(padd(coupling("v"), field("h"), pscale(field("G0"), cnum(0, -1))), SQRT2_INV)
    # Htilde = i sigma2 H* = (H2*, -H1*)
    htilde = (h2_conj, pscale(h1_conj, cnum(-1)))
    hdown = (h1, h2)

    entries: list[dict[str, Any]] = []

    def emit(sector: str, doublet: tuple, left_pair: tuple[str, str], right: str) -> None:
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                y_sym = f"{sector}[{i}][{j}]"
                for slot, left in enumerate(left_pair):
                    operator = pscale(pmul(coupling(y_sym), doublet[slot]), cnum(-1))
                    for (fields, couplings), value in operator.items():
                        scalars = [f for f in fields]
                        structure = "fermion_bilinear_mass" if not scalars else "fermion_scalar_yukawa"
                        combined = polynomial(*coupling_monomials(couplings, value))
                        entries.append(record([bar(left), right] + scalars, combined, structure))

    emit("Yu", htilde, ("uL", "dL"), "uR")
    emit("Yd", hdown, ("uL", "dL"), "dR")
    emit("Ye", hdown, ("nuL", "eL"), "eR")

    conjugates: list[dict[str, Any]] = []
    for entry in entries:
        fields = entry["fields"]
        bar_fields = [f for f in fields if f.endswith("_bar")]
        plain_fermions = [f for f in fields if not f.endswith("_bar") and f not in SCALAR_LABELS]
        scalars = [f for f in fields if f in SCALAR_LABELS]
        swapped = [{"Gp": "Gm", "Gm": "Gp"}.get(s, s) for s in scalars]
        new_fields = [bar(f) for f in plain_fermions] + [f[:-4] for f in bar_fields] + swapped
        conjugated = []
        for m in entry["coefficient"]["monomials"]:
            prefactor = Fraction(m["prefactor"])
            powers = []
            for symbol_name, power in m["powers"]:
                if symbol_name == "I":
                    prefactor = -prefactor
                    powers.append((symbol_name, power))
                elif symbol_name[:2] in ("Yu", "Yd", "Ye") and "[" in symbol_name:
                    plain = symbol_name.split("[")[0]
                    i_idx, j_idx = symbol_name[len(plain):].strip("[]").split("][")
                    powers.append((f"{plain}d[{j_idx}][{i_idx}]", power))
                else:
                    powers.append((symbol_name, power))
            conjugated.append((str(prefactor), tuple(powers)))
        conjugates.append(record(new_fields, polynomial(*conjugated), entry["structure"]))
    return entries + conjugates


def structure_constant(a: int, b: int, c: int) -> int:
    return int((b - a) * (c - b) * (c - a) / 2)


CHARGED_VECTOR = {
    # W^1 = (Wp + Wm)/sqrt2, W^2 = I (Wp - Wm)/sqrt2, W^3 = W3
    0: ((("Wp",), cnum(1)), (("Wm",), cnum(1))),
    1: ((("Wp",), cnum(0, 1)), (("Wm",), cnum(0, -1))),
    2: ((("W3",), cnum(1)), ),
}


def rotate_index(index: int) -> list[tuple[str, tuple, tuple]]:
    """Charged-basis expansion of real index: (label, couplings, value)."""

    if index == 2:
        return [("W3", (), cnum(1))]
    out = []
    for label, value in CHARGED_VECTOR[index]:
        out.append((label[0], (("sqrt2", -1),), value))
    return out


def yang_mills_cubic_entries() -> list[dict[str, Any]]:
    """-g2 eps_abc A_a B_b C_c with role order (dW, W^mu, W^nu); the
    record carries the sorted role assignment and total antisymmetry
    is checked over all permutations."""

    tensor: dict[tuple[str, str, str], tuple] = {}
    for a, b, c in itertools.product(range(3), repeat=3):
        eps = structure_constant(a, b, c)
        if not eps:
            continue
        for la, ca, va in rotate_index(a):
            for lb, cb, vb in rotate_index(b):
                for lc, cc, vc in rotate_index(c):
                    key = (la, lb, lc)
                    couplings = merge_couplings(merge_couplings(ca, cb), cc)
                    value = cmul(cmul(va, vb), cmul(vc, cnum(-eps)))
                    existing = tensor.get(key, ((), (ZERO, ZERO)))
                    if existing[0] and existing[0] != couplings:
                        raise SystemExit("engine B: cubic coupling-power mismatch")
                    tensor[key] = (couplings, cadd(existing[1], value))
    entries: list[dict[str, Any]] = []
    seen: set[tuple[str, ...]] = set()
    for key, (couplings, value) in list(tensor.items()):
        multiset = tuple(sorted(key))
        if multiset in seen:
            continue
        seen.add(multiset)
        if len(set(multiset)) != 3:
            for assignment, (_, other) in tensor.items():
                if tuple(sorted(assignment)) == multiset and other != (ZERO, ZERO):
                    raise SystemExit(f"engine B: repeated-field cubic assignment {assignment}")
            continue
        ref_couplings, ref_value = tensor.get(multiset, ((), (ZERO, ZERO)))
        for permutation in itertools.permutations(range(3)):
            assignment = tuple(multiset[p] for p in permutation)
            sign = 1
            for x in range(3):
                for y in range(x + 1, 3):
                    if permutation[x] > permutation[y]:
                        sign = -sign
            _, value_at = tensor.get(assignment, ((), (ZERO, ZERO)))
            expected = (ref_value[0] * sign, ref_value[1] * sign)
            if value_at != expected:
                raise SystemExit(f"engine B: cubic tensor not antisymmetric at {assignment}")
        if ref_value != (ZERO, ZERO):
            monomials = coupling_monomials(merge_couplings(ref_couplings, (("g2", 1),)), ref_value)
            entries.append(record(list(multiset), polynomial(*monomials), "yang_mills_three_point"))
    return entries


def yang_mills_quartic_entries() -> list[dict[str, Any]]:
    """-(g2^2/4) sum_a (eps_abc P_b Q_c)(eps_ade P_d Q_e) with the P/Q
    pairing routed to its Lorentz pairing label."""

    accumulator: dict[tuple, tuple] = {}
    for a in range(3):
        for b, c, d, e in itertools.product(range(3), repeat=4):
            eps1 = structure_constant(a, b, c)
            eps2 = structure_constant(a, d, e)
            if not eps1 or not eps2:
                continue
            for lb, cb, vb in rotate_index(b):
                for lc, cc, vc in rotate_index(c):
                    for ld, cd, vd in rotate_index(d):
                        for le, ce, ve in rotate_index(e):
                            p_fields = tuple(sorted((lb, ld)))
                            q_fields = tuple(sorted((lc, le)))
                            multiset = tuple(sorted(p_fields + q_fields))
                            partition = tuple(sorted([p_fields, q_fields]))
                            couplings = merge_couplings(
                                merge_couplings(cb, cc), merge_couplings(cd, ce)
                            )
                            value = cmul(cmul(vb, vc), cmul(vd, ve))
                            value = cmul(value, cnum(Fraction(-eps1 * eps2, 4)))
                            key = (multiset, partition, couplings)
                            accumulator[key] = cadd(accumulator.get(key, (ZERO, ZERO)), value)
    pairing_labels = (
        ("yang_mills_four_point_12_34", ((0, 1), (2, 3))),
        ("yang_mills_four_point_13_24", ((0, 2), (1, 3))),
        ("yang_mills_four_point_14_23", ((0, 3), (1, 2))),
    )
    merged: dict[tuple, list] = {}
    for (multiset, partition, couplings), value in accumulator.items():
        if value == (ZERO, ZERO):
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
            raise SystemExit(f"engine B: quartic partition {partition} has no label")
        merged.setdefault((multiset, chosen), []).extend(
            coupling_monomials(merge_couplings(couplings, (("g2", 2),)), value)
        )
    entries: list[dict[str, Any]] = []
    for (multiset, label), monomials in sorted(merged.items()):
        combined = polynomial(*monomials)
        if combined["monomials"]:
            entries.append(record(list(multiset), combined, label))
    return entries


def fermion_current_entries(action: dict[str, Any]) -> list[dict[str, Any]]:
    """Closed current rules over the census multiplets."""

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


def ghost_entries(kappa: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """Ghost blocks from delta Phi = - g_b theta_b T^b Phi and the
    solved gauge functions: scalar block
    - cbar_a xi_a (k_a . delta Phi), gauge block
    + g2 eps_abc (d cbar_a) W^b c^c, everything rotated to the charged
    ghost basis cbar1 = (cp_bar + cm_bar)/sqrt2 (conjugate rotation)
    and c1 = (cp + cm)/sqrt2."""

    entries: list[dict[str, Any]] = []
    ghost_names = {"W1": "t1", "W2": "t2", "W3": "t3", "B": "tB"}

    # delta Phi per ghost direction: -g_b T^b Phi.
    delta_phi = {
        b_name: {
            target: pscale(pmul(component, coupling(GAUGE_COUPLING[b_name])), cnum(-1))
            for target, component in apply_generator(b_name, {n: PHI_POLYS[n] for n in REAL_BASIS}).items()
        }
        for b_name in SIGMA
    }
    scalar_block: dict = {}
    for a_name in SIGMA:
        for b_name in SIGMA:
            block = dot(kappa[a_name], delta_phi[b_name])
            block = pmul(block, pmul(coupling(XI_OF[a_name]), pmul(field("BAR" + a_name), field(ghost_names[b_name]))))
            scalar_block = padd(scalar_block, pscale(block, cnum(-1)))

    # Rotate ghosts to the charged basis.
    scalar_block = to_charged(scalar_block)
    scalar_block = psubst(scalar_block, "BARW1", pmul(padd(field(bar("cp")), field(bar("cm"))), SQRT2_INV))
    scalar_block = psubst(scalar_block, "BARW2", pscale(pmul(padd(field(bar("cp")), pscale(field(bar("cm")), cnum(-1))), SQRT2_INV), cnum(0, -1)))
    scalar_block = psubst(scalar_block, "BARW3", field(bar("c3")))
    scalar_block = psubst(scalar_block, "BARB", field(bar("cB")))
    scalar_block = psubst(scalar_block, "t1", pmul(padd(field("cp"), field("cm")), SQRT2_INV))
    scalar_block = psubst(scalar_block, "t2", pscale(pmul(padd(field("cp"), pscale(field("cm"), cnum(-1))), SQRT2_INV), cnum(0, 1)))
    scalar_block = psubst(scalar_block, "t3", field("c3"))
    scalar_block = psubst(scalar_block, "tB", field("cB"))
    for (fields, couplings), value in sorted(scalar_block.items()):
        monomials = coupling_monomials(couplings, value)
        combined = polynomial(*monomials)
        if combined["monomials"]:
            entries.append(record(list(fields), combined, "ghost_scalar_mass"))

    # Gauge block.
    gauge_block: dict = {}
    for a, b, c in itertools.product(range(3), repeat=3):
        eps = structure_constant(a, b, c)
        if not eps:
            continue
        cbar = ("BARW1", "BARW2", "BARW3")[a]
        wfield = ("W1", "W2", "W3")[b]
        theta = ("t1", "t2", "t3")[c]
        gauge_block = padd(
            gauge_block,
            pscale(pmul(coupling("g2"), pmul(field(cbar), pmul(field(wfield), field(theta)))), cnum(eps)),
        )
    gauge_block = to_charged(gauge_block)
    gauge_block = psubst(gauge_block, "BARW1", pmul(padd(field(bar("cp")), field(bar("cm"))), SQRT2_INV))
    gauge_block = psubst(gauge_block, "BARW2", pscale(pmul(padd(field(bar("cp")), pscale(field(bar("cm")), cnum(-1))), SQRT2_INV), cnum(0, -1)))
    gauge_block = psubst(gauge_block, "BARW3", field(bar("c3")))
    gauge_block = psubst(gauge_block, "t1", pmul(padd(field("cp"), field("cm")), SQRT2_INV))
    gauge_block = psubst(gauge_block, "t2", pscale(pmul(padd(field("cp"), pscale(field("cm"), cnum(-1))), SQRT2_INV), cnum(0, 1)))
    gauge_block = psubst(gauge_block, "t3", field("c3"))
    for (fields, couplings), value in sorted(gauge_block.items()):
        monomials = coupling_monomials(couplings, value)
        combined = polynomial(*monomials)
        if combined["monomials"]:
            entries.append(record(list(fields), combined, "ghost_gauge_derivative"))
    return entries


# ---------------------------------------------------------------------------
# Table assembly
# ---------------------------------------------------------------------------


def build_table() -> dict[str, Any]:
    action = json.loads(ACTION_PATH.read_text(encoding="utf-8"))
    scalar_entries, kappa, cancelled_mixing = scalar_gauge_entries()
    entries = (
        scalar_entries
        + potential_entries()
        + yukawa_entries()
        + yang_mills_cubic_entries()
        + yang_mills_quartic_entries()
        + fermion_current_entries(action)
        + ghost_entries(kappa)
    )
    merged: dict[str, dict[str, Any]] = {}
    cancelled = list(cancelled_mixing)
    for entry in entries:
        key = json.dumps({"fields": entry["fields"], "structure": entry["structure"]}, sort_keys=True)
        if key in merged:
            combined = [
                (m["prefactor"], tuple(tuple(p) for p in m["powers"]))
                for m in merged[key]["coefficient"]["monomials"] + entry["coefficient"]["monomials"]
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
            raise SystemExit(f"engine B: conservation violation {violations} in {entry['fields']}")
    kappa_report = {}
    for v_name, k_vector in sorted(kappa.items()):
        row = {}
        for target, component in k_vector.items():
            for (fields, couplings), value in component.items():
                if fields:
                    raise SystemExit("engine B: field-dependent gauge function")
                name = {"r1": "r1", "r2": "r2", "G0": "G0", "hshift": "h"}[target]
                text_parts = []
                if value[0] != 0:
                    text_parts.append(str(value[0]))
                if value[1] != 0:
                    text_parts.append(f"{value[1]}*I")
                row[name] = {"couplings": list(couplings), "value": "+".join(text_parts)}
        kappa_report[{"W1": "W1", "W2": "W2", "W3": "W3", "B": "B"}[v_name]] = row
    return {
        "schema": "rule_table.v2",
        "engine": "B_generator_assignment",
        "action_subject_digest": action["subject_digest"],
        "gauge_fixing": {
            "form": "G^a = d.W^a + xi_a (k_a . phi), L_gf = -(1/(2 xi_a)) (G^a)^2",
            "solved_kappa": kappa_report,
            "mixing_cancelled": sorted(cancelled, key=lambda item: json.dumps(item, sort_keys=True)),
        },
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
