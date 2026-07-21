#!/usr/bin/env python3
"""Exact finite-A5 and compact coefficient-algebra certificate.

The certificate proves representation arithmetic and enumerates compact-Lie
alternatives under explicit hypotheses.  It does not identify record
projectors, port fluctuations, or their commutators with physical currents.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
SQRT5 = sp.sqrt(5)
PHI = (1 + SQRT5) / 2
PHIBAR = (1 - SQRT5) / 2
CLASS_SIZES = [1, 15, 20, 12, 12]
IRREPS = {
    "1": {"dim": 1, "chi": [1, 1, 1, 1, 1], "rational_orbit": "1"},
    "3": {"dim": 3, "chi": [3, -1, 0, PHI, PHIBAR], "rational_orbit": "3+3prime"},
    "3prime": {"dim": 3, "chi": [3, -1, 0, PHIBAR, PHI], "rational_orbit": "3+3prime"},
    "4": {"dim": 4, "chi": [4, 0, 1, -1, -1], "rational_orbit": "4"},
    "5": {"dim": 5, "chi": [5, 1, -1, 0, 0], "rational_orbit": "5"},
}


def inner(chi_a: list[sp.Expr], chi_b: list[sp.Expr]) -> sp.Expr:
    return sp.simplify(sp.Rational(1, 60) * sum(s*a*b for s, a, b in zip(CLASS_SIZES, chi_a, chi_b, strict=True)))


def decompose(chi: list[sp.Expr]) -> dict[str, int]:
    result: dict[str, int] = {}
    for name, row in IRREPS.items():
        multiplicity = sp.simplify(inner(chi, row["chi"]))
        if multiplicity != 0:
            if not multiplicity.is_Integer:
                raise AssertionError(f"nonintegral multiplicity {name}: {multiplicity}")
            result[name] = int(multiplicity)
    return result


# A5 actions on the icosahedral cell orbits.
ORBIT_CHARS = {
    "vertices_A5_over_C5": [12, 0, 0, 2, 2],
    "edges_A5_over_C2": [30, 2, 0, 0, 0],
    "faces_A5_over_C3": [20, 0, 2, 0, 0],
    "flags_regular_A5": [60, 0, 0, 0, 0],
}
ORBIT_DECOMP = {name: decompose(chi) for name, chi in ORBIT_CHARS.items()}
assert ORBIT_DECOMP["vertices_A5_over_C5"] == {"1": 1, "3": 1, "3prime": 1, "5": 1}

# Gauge adjoint character under opposite-triplet embeddings.
chi1 = IRREPS["1"]["chi"]
chi3 = IRREPS["3"]["chi"]
chi3p = IRREPS["3prime"]["chi"]
chi5 = IRREPS["5"]["chi"]
su3_from_3p = [sp.simplify(x*x-y) for x, y in zip(chi3p, chi1, strict=True)]
assert decompose(su3_from_3p) == {"3prime": 1, "5": 1}
sm_char = [sp.simplify(a+b+c) for a, b, c in zip(su3_from_3p, chi3, chi1, strict=True)]
assert sm_char == ORBIT_CHARS["vertices_A5_over_C5"]

# Restriction multiplicities to face C3 and edge C2 stabilizers.
restrictions: dict[str, dict[str, int]] = {}
for name, row in IRREPS.items():
    d = row["dim"]
    chi = row["chi"]
    restrictions[name] = {
        "C2_trivial": int(sp.Rational(d + chi[1], 2)),
        "C2_sign": int(sp.Rational(d - chi[1], 2)),
        "C3_trivial": int(sp.Rational(d + 2*chi[2], 3)),
        "C3_omega": int(sp.Rational(d - chi[2], 3)),
        "C3_omega2": int(sp.Rational(d - chi[2], 3)),
        "C5_trivial": int(sp.simplify(sp.Rational(1, 5)*(d + 2*chi[3] + 2*chi[4]))),
    }
minimal_c3_dim = min(IRREPS[n]["dim"] for n in IRREPS if restrictions[n]["C3_omega"] > 0)
minimal_c3_irreps = [n for n in IRREPS if restrictions[n]["C3_omega"] > 0 and IRREPS[n]["dim"] == minimal_c3_dim]
assert minimal_c3_dim == 3 and set(minimal_c3_irreps) == {"3", "3prime"}

# Compact-Lie trichotomy.  A group-level A5 action on a compact center torus
# preserves its integral lattice; hence its real center representation is
# rational and must contain 3 and 3' with equal multiplicity.
PORT_IRREPS = ["1", "3", "3prime", "5"]
center_candidates = []
for mask in itertools.product([0, 1], repeat=len(PORT_IRREPS)):
    center = [name for name, keep in zip(PORT_IRREPS, mask, strict=True) if keep]
    if ("3" in center) != ("3prime" in center):
        continue
    complement = [name for name in PORT_IRREPS if name not in center]
    center_candidates.append({
        "center_irreps": center,
        "center_dim": sum(IRREPS[n]["dim"] for n in center),
        "semisimple_irreps": complement,
        "semisimple_dim": sum(IRREPS[n]["dim"] for n in complement),
    })

trichotomy = []
for candidate in center_candidates:
    center = set(candidate["center_irreps"])
    ss = set(candidate["semisimple_irreps"])
    sdim = candidate["semisimple_dim"]
    lie_type = None
    reason = None
    if sdim == 0:
        lie_type = "u(1)^12"
        reason = "all of P12 is central"
    elif sdim == 11 and ss == {"3", "3prime", "5"}:
        lie_type = "su(3)+su(2)+u(1)"
        reason = "the only compact semisimple dimension-11 type is A2+A1"
    elif sdim == 6 and ss == {"3", "3prime"} and center == {"1", "5"}:
        lie_type = "su(2)+su(2)+u(1)^6"
        reason = "the only compact semisimple dimension-6 type is A1+A1"
    elif sdim == 12:
        reason = "A1^4 is the only semisimple type, but A5 cannot permute four ideals and stable A1 ideals cannot supply a 5-irrep"
    elif sdim in {1, 2, 4, 5, 7}:
        reason = "no compact semisimple Lie algebra has this dimension"
    elif sdim == 6:
        reason = "A1+A1 cannot carry the required 1+5 A5 module because A5 cannot permute the two ideals"
    elif sdim == 10:
        reason = "B2/C2 exists dimensionally, but P12 has no compatible rational two-dimensional center and no required A5 module split"
    else:
        reason = "excluded by low-dimensional compact semisimple classification and A5-module structure"
    row = dict(candidate)
    row.update({"viable": lie_type is not None, "lie_type": lie_type, "reason": reason})
    if lie_type is not None:
        trichotomy.append(row)

assert {row["lie_type"] for row in trichotomy} == {
    "u(1)^12",
    "su(2)+su(2)+u(1)^6",
    "su(3)+su(2)+u(1)",
}

# Degree <=2 adjacency-local kinetic relation.  Let the color triplet share a
# kinetic coefficient with the rank-5 color band.  For f(A)=a+bA+cA^2 this
# implies k1 = 3 k2 - 2 k3, independent of which triplet is weak/color.
a, b, c = sp.symbols("a b c", real=True)
r = SQRT5
f = lambda x: a + b*x + c*x*x
kinetic_relations = []
for color_eigenvalue, weak_eigenvalue in [(-r, r), (r, -r)]:
    b_solution = sp.solve(sp.Eq(f(color_eigenvalue), f(-1)), b, dict=True)[0]
    k1 = sp.simplify(f(5).subs(b_solution))
    k2 = sp.simplify(f(weak_eigenvalue).subs(b_solution))
    k3 = sp.simplify(f(color_eigenvalue).subs(b_solution))
    assert sp.simplify(k1 - 3*k2 + 2*k3) == 0
    kinetic_relations.append({
        "color_eigenvalue": str(color_eigenvalue),
        "weak_eigenvalue": str(weak_eigenvalue),
        "color_coherence_condition": str(b_solution),
        "relation": "k1 = 3*k2 - 2*k3",
    })

# A tempting but invalid Higgs identification: the missing A5 irrep 4 is real
# type, with scalar commutant.  Therefore it has no A5-invariant complex
# structure J^2=-1 and cannot be a hypercharge-carrying complex doublet while
# commuting with the A5 action.
higgs_dimension_match_status = {
    "dimension_identity": "P12 + 4 has dimensions 1+3+3+4+5 = 16",
    "physical_identification": "REJECTED_WITHOUT_NONCOMMUTING_EXTENSION",
    "reason": "the real A5 irrep 4 has End_A5(4)=R, so no A5-invariant complex structure exists for a commuting U(1) hypercharge action",
}

payload = {
    "schema": "OPH A5 screen/SM closure certificate v1",
    "vertex_module": {
        "character": [str(x) for x in ORBIT_CHARS["vertices_A5_over_C5"]],
        "decomposition": ORBIT_DECOMP["vertices_A5_over_C5"],
        "sm_adjoint_character": [str(x) for x in sm_char],
        "exact_match": sm_char == ORBIT_CHARS["vertices_A5_over_C5"],
    },
    "cell_orbit_modules": ORBIT_DECOMP,
    "stabilizer_restrictions": restrictions,
    "generation_bridge": {
        "theorem": "A nontrivial face-C3 phase extended to a dimension-minimal irreducible linear A5 carrier forces dimension 3 (3 or 3').",
        "minimal_dimension": minimal_c3_dim,
        "minimal_irreps": minimal_c3_irreps,
        "CKM_phase_count_at_Ng_3": 1,
        "physical_gate": "derive a quotient-visible face-phase-to-family attachment",
    },
    "compact_lie_trichotomy": trichotomy,
    "coefficient_rank5_noncentral_corollary": {
        "theorem": "For a group-level A5 action on the compact coefficient algebra, noncentrality of its W5 summand selects su(3)+su(2)+u(1).",
        "proof_summary": "the other two viable compact outcomes place the W5 coefficient summand in the center",
        "physical_gate": "A5-COMMON-ACTION plus a source-derived W5 current commutator; PORT-CURRENT-INNER is the stronger route",
    },
    "degree2_kinetic_relations": kinetic_relations,
    "higgs_missing_four_audit": higgs_dimension_match_status,
    "theorem_scope": {
        "exact_here": "conditional finite A5 vertex module and compact coefficient-algebra classification",
        "exact_companion_constructions": [
            "u(3)+so(3) coefficient bracket",
            "trace-balanced block integration under coefficient 1/2",
            "six-axis lattice quotient isomorphic to Z6",
        ],
        "physical_gates": [
            "the echosahedral source-selector certificate closes UD12 and RP-A5 only on the declared quotient-visible carrier lineage",
            "PORT-CURRENT-INNER is closed conditionally on the declared response branch by port_current_inner_certificate.py: a full-rank compact skew-adjoint commutator-closed current algebra with inner A5 action and a register-relabeling no-go; the response packet (four band scales, one common odd sign) is a typed branch premise whose physical measurement is the remaining producer",
            "BLOCK-DETERMINANT-BALANCE, PORT-SPIN-LIFT, and AXIS-CENTER-DESCENT identify the physical group and quotient",
            "MAR supplies the matter, hypercharge, and faithful tensor-kernel premises",
        ],
    },
    "claim_boundary": {
        "proved_here": "finite A5 arithmetic and conditional coefficient-algebra classification",
        "work_in_progress": [
            "physical current realization and refinement intertwining",
            "physical block/axis descent and MAR matter realization",
            "face-phase family attachment and gauge kinetic operator",
        ],
    },
}

out = HERE / "a5_screen_sm_closure.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
