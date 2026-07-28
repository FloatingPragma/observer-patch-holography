#!/usr/bin/env python3
"""Canonical vertex-record format and vocabulary for Workstream C.

Canonicalization contract only: this module fixes the record container,
the component-field vocabulary, the symbol alphabet, the polynomial
coefficient encoding, the structure vocabulary, and the record hash.  It
contains no derivation logic and no coupling assignment, so both rule
engines import it without sharing rule content.

Conventions fixed here, binding on both engines:

* Fields are the broken-phase components listed in FIELD_TABLE with
  their exact quantum numbers; conjugates carry the ``_bar`` suffix.
* Fermion vertices are flavor-collapsed: a flavor-diagonal coupling is
  one record over the template pair; a matrix coupling carries the
  explicit generation indices in its symbols.
* A coefficient is a polynomial: a sorted tuple of monomials, each an
  exact Fraction prefactor with sorted (symbol, power) pairs over the
  symbol alphabet; ``sqrt2`` and the imaginary unit ``I`` are symbols
  held in canonical residue form: even powers fold into the prefactor,
  so the stored exponent of each is 0 or 1 and equal-valued monomials
  have equal encodings.
* Kinematic normalization is fixed per structure label on the sorted
  field list (f1 <= f2 <= ...):
  - ``vector_scalar_mixing``: coefficient of V^mu d_mu phi.
  - ``scalar_scalar_gauge``: coefficient c of (d_mu phi_a) phi_b V^mu
    where (phi_a, phi_b) is the sorted scalar pair; the reversed
    assignment must carry exactly -c (the engine proves antisymmetry).
  - ``yang_mills_three_point``: with roles A = d_mu V_nu, B = V^mu,
    C = V^nu, the coefficient of the role assignment (A, B, C) =
    (f1, f2, f3); total antisymmetry over role permutations is proved
    by the emitting engine.
  - ``yang_mills_four_point_12_34`` / ``_13_24`` / ``_14_23``: the
    Lorentz scalar (f_i . f_j)(f_k . f_l) with the index pairing taken
    from the label positions in the sorted field list.
  - ``ghost_gauge_derivative``: coefficient of (d_mu cbar) V^mu c with
    the roles fixed by the barred ghost, the vector, and the unbarred
    ghost in the record.
* Canonically normalized kinetic bilinears (scalar d phi d phi, vector
  field-strength quadratic, ghost cbar d^2 c, fermion kinetic terms) and
  the pure gauge-fixing terms -(1/(2 xi)) (d . V)^2 are not records; all
  mass-type and mixing bilinears of the gauge-fixed action are.
* Flavor-matrix symbols bind positionally: in any record, the row index
  of ``Y[i][j]`` binds to the barred fermion template and the column
  index to the unbarred one.  ``Yud``/``Ydd``/``Yed`` name the dagger
  matrices, ``Yud[a][b] = conj(Yu[b][a])``, with the same positional
  binding, so hermitian-conjugate records are expressible exactly.
* The vacuum minimum condition is NOT imposed: tadpole records are
  retained, per the direct FJ contract.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any, Sequence

# Component-field vocabulary: name -> (electric charge as a Fraction
# string in units of e with Q = T3 + Y, color triality, ghost number,
# fermion number, spin label).
FIELD_TABLE: dict[str, dict[str, Any]] = {
    "Gl": {"charge": "0", "triality": 0, "ghost": 0, "fermion": 0, "spin": "vector"},
    "Wp": {"charge": "1", "triality": 0, "ghost": 0, "fermion": 0, "spin": "vector"},
    "Wm": {"charge": "-1", "triality": 0, "ghost": 0, "fermion": 0, "spin": "vector"},
    "W3": {"charge": "0", "triality": 0, "ghost": 0, "fermion": 0, "spin": "vector"},
    "B": {"charge": "0", "triality": 0, "ghost": 0, "fermion": 0, "spin": "vector"},
    "h": {"charge": "0", "triality": 0, "ghost": 0, "fermion": 0, "spin": "scalar"},
    "G0": {"charge": "0", "triality": 0, "ghost": 0, "fermion": 0, "spin": "scalar"},
    "Gp": {"charge": "1", "triality": 0, "ghost": 0, "fermion": 0, "spin": "scalar"},
    "Gm": {"charge": "-1", "triality": 0, "ghost": 0, "fermion": 0, "spin": "scalar"},
    "uL": {"charge": "2/3", "triality": 1, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "dL": {"charge": "-1/3", "triality": 1, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "uR": {"charge": "2/3", "triality": 1, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "dR": {"charge": "-1/3", "triality": 1, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "nuL": {"charge": "0", "triality": 0, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "eL": {"charge": "-1", "triality": 0, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "eR": {"charge": "-1", "triality": 0, "ghost": 0, "fermion": 1, "spin": "weyl"},
    "cp": {"charge": "1", "triality": 0, "ghost": 1, "fermion": 0, "spin": "ghost"},
    "cm": {"charge": "-1", "triality": 0, "ghost": 1, "fermion": 0, "spin": "ghost"},
    "c3": {"charge": "0", "triality": 0, "ghost": 1, "fermion": 0, "spin": "ghost"},
    "cB": {"charge": "0", "triality": 0, "ghost": 1, "fermion": 0, "spin": "ghost"},
    "cG": {"charge": "0", "triality": 0, "ghost": 1, "fermion": 0, "spin": "ghost"},
}

SYMBOL_ALPHABET = (
    "g1", "g2", "g3", "v", "lam", "mu2", "xi1", "xi2", "xi3",
    "sqrt2", "I",
) + tuple(
    f"{sector}[{i}][{j}]"
    for sector in ("Yu", "Yd", "Ye", "Yud", "Ydd", "Yed")
    for i in (1, 2, 3)
    for j in (1, 2, 3)
)

STRUCTURE_VOCABULARY = (
    "yang_mills_three_point",       # role-ordered cubic, see module docstring
    "yang_mills_four_point_12_34",  # (f1.f2)(f3.f4) Lorentz pairing
    "yang_mills_four_point_13_24",  # (f1.f3)(f2.f4) Lorentz pairing
    "yang_mills_four_point_14_23",  # (f1.f4)(f2.f3) Lorentz pairing
    "vector_bilinear_mass",         # g_{mu nu} bilinear
    "vector_scalar_mixing",         # V^mu d_mu phi bilinear
    "scalar_gauge_gauge",           # g_{mu nu} S V V
    "scalar_scalar_gauge",          # antisymmetric derivative coupling
    "scalar_scalar_gauge_gauge",    # g_{mu nu} S S V V seagull
    "scalar_potential",             # no tensor structure
    "scalar_bilinear_mass",         # scalar mass-type bilinear
    "scalar_tadpole",               # single-scalar vertex, retained
    "fermion_vector_current",       # gamma^mu chiral current, unit norm
    "fermion_scalar_yukawa",        # chiral scalar coupling
    "fermion_bilinear_mass",        # v-induced fermion bilinear
    "ghost_gauge_derivative",       # (d cbar) V c coupling
    "ghost_scalar_mass",            # R_xi ghost mass and ghost-scalar term
)


def bar(field: str) -> str:
    return f"{field}_bar"


def field_quantum_numbers(label: str) -> dict[str, Any]:
    base = label[:-4] if label.endswith("_bar") else label
    if base not in FIELD_TABLE:
        raise ValueError(f"unknown component field {base}")
    entry = dict(FIELD_TABLE[base])
    if label.endswith("_bar"):
        entry["charge"] = str(-Fraction(entry["charge"]))
        entry["triality"] = -entry["triality"]
        entry["ghost"] = -entry["ghost"]
        entry["fermion"] = -entry["fermion"]
    return entry


def monomial(prefactor: Fraction | int, *powers: tuple[str, int]) -> tuple:
    combined: dict[str, int] = {}
    for symbol_name, power in powers:
        if symbol_name not in SYMBOL_ALPHABET:
            raise ValueError(f"symbol {symbol_name} is not in the alphabet")
        combined[symbol_name] = combined.get(symbol_name, 0) + int(power)
    value = Fraction(prefactor)
    for symbol_name, unit_square in (("sqrt2", Fraction(2)), ("I", Fraction(-1))):
        power = combined.get(symbol_name, 0)
        residue = power % 2
        value *= unit_square ** ((power - residue) // 2)
        combined[symbol_name] = residue
    ordered = tuple(sorted((s, p) for s, p in combined.items() if p != 0))
    return (str(value), ordered)


def polynomial(*monomials: tuple) -> dict[str, Any]:
    combined: dict[tuple, Fraction] = {}
    for raw_prefactor, raw_powers in monomials:
        prefactor, powers = monomial(Fraction(raw_prefactor), *raw_powers)
        combined[powers] = combined.get(powers, Fraction(0)) + Fraction(prefactor)
    cleaned = sorted(
        (powers, value) for powers, value in combined.items() if value != 0
    )
    return {
        "monomials": [
            {"prefactor": str(value), "powers": [list(p) for p in powers]}
            for powers, value in cleaned
        ]
    }


def record(fields: Sequence[str], coefficient: dict[str, Any], structure: str) -> dict[str, Any]:
    if structure not in STRUCTURE_VOCABULARY:
        raise ValueError(f"unknown structure label {structure}")
    for label in fields:
        field_quantum_numbers(label)
    if not coefficient.get("monomials"):
        raise ValueError("a record needs a nonzero coefficient polynomial")
    return {
        "fields": sorted(fields),
        "coefficient": coefficient,
        "structure": structure,
    }


def conservation_violations(fields: Sequence[str]) -> list[str]:
    charge = Fraction(0)
    triality = 0
    ghost = 0
    violations = []
    for label in fields:
        entry = field_quantum_numbers(label)
        charge += Fraction(entry["charge"])
        triality += entry["triality"]
        ghost += entry["ghost"]
    if charge != 0:
        violations.append(f"electric charge {charge}")
    if triality % 3 != 0:
        violations.append(f"color triality {triality}")
    if ghost != 0:
        violations.append(f"ghost number {ghost}")
    return violations


def vertex_hash(entry: dict[str, Any]) -> str:
    payload = json.dumps(entry, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def table_digest(entries: Sequence[dict[str, Any]]) -> str:
    hashes = sorted(vertex_hash(entry) for entry in entries)
    payload = json.dumps(hashes, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
