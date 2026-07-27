#!/usr/bin/env python3
"""Exact certificate for GitHub issue #611: RESPONSE-GRAMMAR COMPLETENESS.

The port space is R^12 over the icosahedral port adjacency of the certified
carrier, with every arithmetic decision exact in Q(sqrt5).  The commutant of
the sixty orientation-preserving carrier rotations is four-dimensional; the
machine-checked statement is Lean/Screen/A5Commutant.lean
(`commutant_decomposition`, `orbitals_independent`, `equivariant_iff_commutes`):
every matrix equivariant under the listed sixty rotations is a unique rational
combination of the identity, adjacency, distance-two, and antipode orbitals,
and entry equivariance is equivalent to commutation with the linear action.
This certificate derives, rather than assumes:

* the four isotypic projectors P_1, P_3, P_3', P_5 of the A5 port
  decomposition 1 + 3 + 3' + 5, constructed by Lagrange interpolation in the
  adjacency matrix over its exact eigenvalues 5, sqrt5, -sqrt5, -1, with
  idempotence, mutual orthogonality, partition of unity, the eigenvalue
  relations, and ranks 1, 3, 3, 5 verified by exact traces;
* the executable four-orbit transport fact behind the commutant dimension:
  the sixty reconstructed rotations act on ordered port pairs with exactly
  four orbits (diagonal, adjacent, distance-two, antipodal), and the orbital
  basis is exactly the projector basis under an invertible change of basis;
* the complete quotient-visible response grammar: an equivariant involution
  lies in the four-dimensional commutant, the projectors are an exact basis
  of orthogonal idempotents summing to the identity, so every equivariant
  involution is a sign pattern sum(eps_s P_s) with eps_s = +-1 on each
  isotypic sector (on each real-irreducible block the commutant is the real
  scalars by Schur, and a scalar squaring to one is +-1); the mechanism
  space is the sixteen sign patterns, of which fourteen are nontrivial;
* the operational selection: among the fourteen admissible responses exactly
  two are realized as signed graph automorphisms, namely +J and -J for the
  unique nonidentity central involution
  J = (A^3 - 4 A^2 - 5 A + 10 I)/10, which is the antipode permutation; the
  hash-pinned #566 inverse-port producer selects R = -J by the
  impulse/readback sign convention;
* the mandatory countermodel pair: two inequivalent admissible responses
  survive inside the equivariance + involution grammar without the
  signed-automorphism clause, so incidence with equivariance and involution
  alone does not force +-J and the operational readback clause is
  load-bearing;
* presentation invariance: conjugation by every one of the sixty carrier
  rotations fixes each admissible response, and the rotation list is closed
  under composition.

Coverage is proved inside the declared operational grammar, unique up to the
conventional sign.  Forcing the grammar itself from the three axioms stays
open, so the bounded exit is `independence_limited`.  No floating point
appears in a proof step.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import sys
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import port_current_inner_certificate as p566  # noqa: E402

SCHEMA = "oph.response_grammar_completeness_certificate.v1"
RECEIPT_SCHEMA = "oph.response_grammar_completeness_receipt.v1"
NEGATIVE_SCHEMA = "oph.response_grammar_completeness_negative_controls.v1"
PRODUCER_SCHEMA = "oph.port_current_response_manifest.v5"
LEAN_COMMUTANT_MODULE = "Lean/Screen/A5Commutant.lean"

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

F5 = p566.F5
ZERO = p566.ZERO
ONE = p566.ONE
MINUS_ONE = -p566.ONE
SQRT5 = F5(0, 1)

RMat = p566.RMat

# Isotypic sectors of the twelve-port module under A5, with the exact
# adjacency eigenvalue and rank of each sector.
SECTORS = ("1", "3", "3p", "5")
EIGENVALUES: dict[str, F5] = {
    "1": F5(5),
    "3": SQRT5,
    "3p": -SQRT5,
    "5": MINUS_ONE,
}
SECTOR_RANKS = {"1": 1, "3": 3, "3p": 3, "5": 5}

SignPattern = tuple[int, int, int, int]
ALL_PATTERNS: tuple[SignPattern, ...] = tuple(
    itertools.product((1, -1), repeat=4)
)
PATTERN_J: SignPattern = (1, -1, -1, 1)
PATTERN_MINUS_J: SignPattern = (-1, 1, 1, -1)
PATTERN_FIVE_FLIP: SignPattern = (1, 1, 1, -1)
GRAMMAR_CLAUSES = ("equivariant", "involutive", "signed", "nontrivial")
OPERATIONAL_CLAUSE = "signed_graph_automorphism_readback"


# ---------------------------------------------------------------------------
# Small exact matrix helpers over Q(sqrt5)
# ---------------------------------------------------------------------------


def rident(n: int) -> RMat:
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def madd(x: RMat, y: RMat) -> RMat:
    return [[x[i][j] + y[i][j] for j in range(len(x[0]))] for i in range(len(x))]


def msub(x: RMat, y: RMat) -> RMat:
    return [[x[i][j] - y[i][j] for j in range(len(x[0]))] for i in range(len(x))]


def smul(c: F5, x: RMat) -> RMat:
    return [[c * entry for entry in row] for row in x]


def mat_eq(x: RMat, y: RMat) -> bool:
    return all(
        x[i][j] == y[i][j] for i in range(len(x)) for j in range(len(x[0]))
    )


def rtrace(x: RMat) -> F5:
    total = ZERO
    for i in range(len(x)):
        total = total + x[i][i]
    return total


def cross(u: p566.Vec3, v: p566.Vec3) -> p566.Vec3:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def columns3(u: p566.Vec3, v: p566.Vec3, w: p566.Vec3) -> RMat:
    return [[u[i], v[i], w[i]] for i in range(3)]


def inv3(m: RMat) -> RMat:
    det = p566.det3(m)
    require(not det.is_zero(), "ROTATION_GROUP", "singular 3x3 basis matrix")
    scale = det.inv()
    cof = [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
        ],
        [
            m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
        ],
        [
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ]
    return [[cof[i][j] * scale for j in range(3)] for i in range(3)]


def f5_sort_key(value: F5) -> tuple[Any, Any]:
    return (value.a, value.b)


# ---------------------------------------------------------------------------
# The port model: twelve vertices, adjacency, distances, antipode
# ---------------------------------------------------------------------------


def port_model() -> tuple[list[p566.Vec3], list[list[int]], list[list[int]], list[int]]:
    """The exact icosahedral port model from the Q(sqrt5) vertex coordinates."""

    verts = p566.standard_vertices()
    require(len(verts) == 12, "PORT_MODEL", "the port model must carry twelve vertices")
    distance = [
        [p566.vertex_distance(verts[i], verts[j]) for j in range(12)]
        for i in range(12)
    ]
    adjacency = [[1 if distance[i][j] == 1 else 0 for j in range(12)] for i in range(12)]
    antipode: list[int] = []
    for i in range(12):
        require(
            sum(adjacency[i]) == 5,
            "PORT_MODEL",
            f"vertex {i} does not have exactly five neighbours",
        )
        require(
            distance[i].count(3) == 1,
            "PORT_MODEL",
            f"vertex {i} does not have exactly one antipode",
        )
        antipode.append(distance[i].index(3))
    for i in range(12):
        require(
            antipode[antipode[i]] == i and antipode[i] != i,
            "PORT_MODEL",
            "the antipode map is not a fixed-point-free involution",
        )
        for j in range(12):
            require(
                adjacency[antipode[i]][antipode[j]] == adjacency[i][j],
                "PORT_MODEL",
                "the antipode map does not preserve port adjacency",
            )
    return verts, adjacency, distance, antipode


# ---------------------------------------------------------------------------
# The sixty orientation-preserving rotations, reconstructed exactly
# ---------------------------------------------------------------------------


def rotation_permutations(
    verts: Sequence[p566.Vec3], adjacency: Sequence[Sequence[int]]
) -> list[tuple[int, ...]]:
    """All sixty rotations as port permutations, from exact linear maps.

    A rotation is determined by the image of one directed edge together
    with orientation, so the sixty directed edges enumerate the group.
    Each candidate linear map is built exactly in Q(sqrt5), checked to
    have determinant one, and checked to permute the vertex set.
    """

    index = {verts[i]: i for i in range(12)}
    base = 0
    neighbour = next(j for j in range(12) if adjacency[base][j])
    source = columns3(
        verts[base], verts[neighbour], cross(verts[base], verts[neighbour])
    )
    source_inverse = inv3(source)
    perms: list[tuple[int, ...]] = []
    for a in range(12):
        for b in range(12):
            if not adjacency[a][b]:
                continue
            target = columns3(verts[a], verts[b], cross(verts[a], verts[b]))
            rotation = p566.rmul(target, source_inverse)
            require(
                p566.det3(rotation) == ONE,
                "ROTATION_GROUP",
                "a directed-edge candidate is not orientation preserving",
            )
            images: list[int] = []
            for vertex in verts:
                image = p566.apply3(rotation, vertex)
                require(
                    image in index,
                    "ROTATION_GROUP",
                    "a candidate rotation does not preserve the vertex set",
                )
                images.append(index[image])
            perms.append(tuple(images))
    require(
        len(perms) == 60 and len(set(perms)) == 60,
        "ROTATION_GROUP",
        f"expected sixty distinct rotations, got {len(set(perms))}",
    )
    require(
        tuple(range(12)) in perms,
        "ROTATION_GROUP",
        "the identity rotation is missing",
    )
    perm_set = set(perms)
    for g in perms:
        for h in perms:
            composed = tuple(g[h[k]] for k in range(12))
            require(
                composed in perm_set,
                "ROTATION_GROUP",
                "the rotation list is not closed under composition",
            )
    for g in perms:
        for i in range(12):
            for j in range(12):
                require(
                    adjacency[g[i]][g[j]] == adjacency[i][j],
                    "ROTATION_GROUP",
                    "a rotation does not preserve port adjacency",
                )
    return perms


def pair_orbits(
    perms: Sequence[tuple[int, ...]], distance: Sequence[Sequence[int]]
) -> dict[str, list[int]]:
    """Orbits of ordered port pairs under the sixty rotations.

    Exactly four orbits, coinciding with the four distance classes; this is
    the executable counterpart of the machine-checked `pair_transport` and
    `commutant_decomposition` in Lean/Screen/A5Commutant.lean, hence of the
    commutant dimension four.
    """

    unseen = {(i, j) for i in range(12) for j in range(12)}
    orbits: list[set[tuple[int, int]]] = []
    while unseen:
        seed = min(unseen)
        orbit = {seed}
        frontier = [seed]
        while frontier:
            i, j = frontier.pop()
            for g in perms:
                pair = (g[i], g[j])
                if pair not in orbit:
                    orbit.add(pair)
                    frontier.append(pair)
        unseen -= orbit
        orbits.append(orbit)
    require(
        len(orbits) == 4,
        "PAIR_ORBITS",
        f"expected four ordered-pair orbits, got {len(orbits)}",
    )
    names = {0: "diagonal", 1: "adjacent", 2: "distance_two", 3: "antipodal"}
    representatives: dict[str, list[int]] = {}
    for orbit in orbits:
        classes = {distance[i][j] for i, j in orbit}
        require(
            len(classes) == 1,
            "PAIR_ORBITS",
            "an ordered-pair orbit mixes distance classes",
        )
        representatives[names[classes.pop()]] = list(min(orbit))
    require(
        set(representatives) == set(names.values()),
        "PAIR_ORBITS",
        "the four ordered-pair orbits do not match the four distance classes",
    )
    return representatives


# ---------------------------------------------------------------------------
# Isotypic projectors by exact Lagrange interpolation in the adjacency
# ---------------------------------------------------------------------------


def adjacency_powers(adjacency: Sequence[Sequence[int]]) -> tuple[RMat, RMat, RMat, RMat]:
    a_mat = [[F5(adjacency[i][j]) for j in range(12)] for i in range(12)]
    a2 = p566.rmul(a_mat, a_mat)
    a3 = p566.rmul(a2, a_mat)
    return rident(12), a_mat, a2, a3


def isotypic_projectors(adjacency: Sequence[Sequence[int]]) -> dict[str, RMat]:
    """The four isotypic projectors as exact cubic polynomials in A.

    P_s = prod_{t != s} (A - e_t I) / (e_s - e_t) over the four exact
    eigenvalues.  Idempotence, mutual orthogonality, partition of unity,
    the eigenvalue relations, and the ranks are all verified exactly; a
    tampered adjacency breaks the quartic minimal polynomial and fails
    closed here.
    """

    identity, a_mat, a2, a3 = adjacency_powers(adjacency)
    projectors: dict[str, RMat] = {}
    for sector in SECTORS:
        own = EIGENVALUES[sector]
        others = [EIGENVALUES[t] for t in SECTORS if t != sector]
        t1, t2, t3 = others
        c2 = -(t1 + t2 + t3)
        c1 = t1 * t2 + t1 * t3 + t2 * t3
        c0 = -(t1 * t2 * t3)
        denominator = (own - t1) * (own - t2) * (own - t3)
        numerator = madd(madd(a3, smul(c2, a2)), madd(smul(c1, a_mat), smul(c0, identity)))
        projectors[sector] = smul(denominator.inv(), numerator)
    for sector in SECTORS:
        p_s = projectors[sector]
        require(
            mat_eq(p566.rmul(a_mat, p_s), smul(EIGENVALUES[sector], p_s)),
            "PROJECTOR_SYSTEM",
            f"A does not act as the declared eigenvalue on sector {sector}",
        )
        require(
            mat_eq(p566.rmul(p_s, p_s), p_s),
            "PROJECTOR_SYSTEM",
            f"projector {sector} is not idempotent",
        )
        require(
            rtrace(p_s) == F5(SECTOR_RANKS[sector]),
            "PROJECTOR_SYSTEM",
            f"projector {sector} does not have rank {SECTOR_RANKS[sector]}",
        )
    for s, t in itertools.combinations(SECTORS, 2):
        require(
            all(
                entry.is_zero()
                for row in p566.rmul(projectors[s], projectors[t])
                for entry in row
            ),
            "PROJECTOR_SYSTEM",
            f"projectors {s} and {t} are not orthogonal",
        )
    total = rident(12)
    partition = projectors["1"]
    for sector in SECTORS[1:]:
        partition = madd(partition, projectors[sector])
    require(
        mat_eq(partition, total),
        "PROJECTOR_SYSTEM",
        "the projectors do not sum to the identity",
    )
    for s, t in itertools.product(SECTORS, repeat=2):
        expected = F5(SECTOR_RANKS[s]) if s == t else ZERO
        require(
            rtrace(p566.rmul(projectors[s], projectors[t])) == expected,
            "PROJECTOR_SYSTEM",
            "the projector Gram matrix is not diagonal with the sector ranks",
        )
    return projectors


def frame_band_identification(
    verts: Sequence[p566.Vec3], adjacency: Sequence[Sequence[int]]
) -> None:
    """Sector 3 carries the vertex frame: neighbour sums equal sqrt5 times
    the vertex coordinates, exactly."""

    for i in range(12):
        total = (ZERO, ZERO, ZERO)
        for j in range(12):
            if adjacency[i][j]:
                total = (
                    total[0] + verts[j][0],
                    total[1] + verts[j][1],
                    total[2] + verts[j][2],
                )
        expected = (SQRT5 * verts[i][0], SQRT5 * verts[i][1], SQRT5 * verts[i][2])
        require(
            total == expected,
            "SECTOR_IDENTIFICATION",
            "the vertex frame is not the sqrt5 adjacency eigenband",
        )


# ---------------------------------------------------------------------------
# The central involution J and the orbital/projector change of basis
# ---------------------------------------------------------------------------


def permutation_matrix(perm: Sequence[int]) -> RMat:
    return [
        [ONE if j == perm[i] else ZERO for j in range(12)]
        for i in range(12)
    ]


def central_involution(
    adjacency: Sequence[Sequence[int]],
    antipode: Sequence[int],
    projectors: Mapping[str, RMat],
) -> RMat:
    """J = (A^3 - 4 A^2 - 5 A + 10 I)/10 equals the antipode permutation and
    the sign pattern (+, -, -, +) on (1, 3, 3', 5), exactly."""

    identity, a_mat, a2, a3 = adjacency_powers(adjacency)
    polynomial = smul(
        F5(1) / F5(10),
        madd(
            msub(a3, smul(F5(4), a2)),
            msub(smul(F5(10), identity), smul(F5(5), a_mat)),
        ),
    )
    j_perm = permutation_matrix(antipode)
    require(
        mat_eq(polynomial, j_perm),
        "CENTRAL_INVOLUTION",
        "the cubic polynomial identity for J does not reproduce the antipode permutation",
    )
    j_signed = sign_pattern_matrix(projectors, PATTERN_J)
    require(
        mat_eq(j_perm, j_signed),
        "CENTRAL_INVOLUTION",
        "J is not the (+,-,-,+) isotypic sign pattern",
    )
    require(
        mat_eq(p566.rmul(j_perm, j_perm), rident(12)) and not mat_eq(j_perm, rident(12)),
        "CENTRAL_INVOLUTION",
        "J is not a nonidentity involution",
    )
    return j_perm


def orbital_basis_identities(
    adjacency: Sequence[Sequence[int]],
    distance: Sequence[Sequence[int]],
    j_mat: RMat,
    projectors: Mapping[str, RMat],
) -> dict[str, dict[str, str]]:
    """The orbital commutant basis of Lean/Screen/A5Commutant.lean equals the
    projector basis under an exact invertible change of basis."""

    identity, a_mat, a2, _ = adjacency_powers(adjacency)
    n_mat = [
        [ONE if distance[i][j] == 2 else ZERO for j in range(12)]
        for i in range(12)
    ]
    half = F5(1) / F5(2)
    require(
        mat_eq(n_mat, smul(half, msub(msub(a2, smul(F5(5), identity)), smul(F5(2), a_mat)))),
        "COMMUTANT_BASIS",
        "the distance-two orbital is not (A^2 - 5I - 2A)/2",
    )
    coefficients: dict[str, dict[str, F5]] = {
        "identity": {"1": ONE, "3": ONE, "3p": ONE, "5": ONE},
        "adjacency": {"1": F5(5), "3": SQRT5, "3p": -SQRT5, "5": MINUS_ONE},
        "distance_two": {"1": F5(5), "3": -SQRT5, "3p": SQRT5, "5": MINUS_ONE},
        "antipode": {"1": ONE, "3": MINUS_ONE, "3p": MINUS_ONE, "5": ONE},
    }
    orbitals = {
        "identity": identity,
        "adjacency": a_mat,
        "distance_two": n_mat,
        "antipode": j_mat,
    }
    for name, orbital in orbitals.items():
        combined = p566.rzeros(12, 12)
        for sector in SECTORS:
            combined = madd(combined, smul(coefficients[name][sector], projectors[sector]))
        require(
            mat_eq(orbital, combined),
            "COMMUTANT_BASIS",
            f"the {name} orbital is not the declared projector combination",
        )
    all_ones = [[ONE for _ in range(12)] for _ in range(12)]
    require(
        mat_eq(madd(madd(identity, a_mat), madd(n_mat, j_mat)), all_ones)
        and mat_eq(all_ones, smul(F5(12), projectors["1"])),
        "COMMUTANT_BASIS",
        "the four orbitals do not resolve the all-ones matrix as 12 P_1",
    )
    return {
        name: {sector: coefficients[name][sector].text() for sector in SECTORS}
        for name in orbitals
    }


# ---------------------------------------------------------------------------
# The sign-pattern grammar and the signed-automorphism selection
# ---------------------------------------------------------------------------


def pattern_name(pattern: SignPattern) -> str:
    return "".join("+" if sign == 1 else "-" for sign in pattern)


def sign_pattern_matrix(projectors: Mapping[str, RMat], pattern: SignPattern) -> RMat:
    result = p566.rzeros(12, 12)
    for sector, sign in zip(SECTORS, pattern):
        result = madd(result, smul(F5(sign), projectors[sector]))
    return result


def validate_admissible_set(claimed: Sequence[SignPattern]) -> list[SignPattern]:
    """The admissible set is exactly the fourteen nontrivial sign patterns;
    any claimed extra or trivial pattern is rejected."""

    rows = [tuple(pattern) for pattern in claimed]
    require(
        len(rows) == len(set(rows)),
        "ADMISSIBLE_SET",
        "the claimed admissible set repeats a sign pattern",
    )
    for pattern in rows:
        require(
            pattern in ALL_PATTERNS,
            "ADMISSIBLE_SET",
            f"{pattern} is not a +-1 sign pattern on the four sectors",
        )
        require(
            len(set(pattern)) == 2,
            "ADMISSIBLE_SET",
            f"{pattern_name(pattern)} is +-identity and is not admissible",
        )
    require(
        len(rows) == 14,
        "ADMISSIBLE_SET",
        f"the admissible set has exactly fourteen patterns, got {len(rows)}",
    )
    return rows


def signed_automorphism_report(
    response: RMat, adjacency: Sequence[Sequence[int]]
) -> dict[str, Any]:
    """Whether a response matrix is a signed graph automorphism: entries in
    {0, +1, -1}, exactly one nonzero per row and column, and the underlying
    permutation an incidence automorphism."""

    for i in range(12):
        for j in range(12):
            entry = response[i][j]
            if not (entry == ZERO or entry == ONE or entry == MINUS_ONE):
                return {
                    "realized": False,
                    "offending_entry": {"row": i, "col": j, "value": entry.text()},
                }
    permutation: list[int] = []
    signs: list[int] = []
    for i in range(12):
        support = [j for j in range(12) if not response[i][j].is_zero()]
        if len(support) != 1:
            return {
                "realized": False,
                "offending_entry": {
                    "row": i,
                    "col": -1,
                    "value": f"row support size {len(support)}",
                },
            }
        permutation.append(support[0])
        signs.append(1 if response[i][support[0]] == ONE else -1)
    if sorted(permutation) != list(range(12)):
        return {
            "realized": False,
            "offending_entry": {"row": -1, "col": -1, "value": "columns not a permutation"},
        }
    for i in range(12):
        for j in range(12):
            if adjacency[permutation[i]][permutation[j]] != adjacency[i][j]:
                return {
                    "realized": False,
                    "offending_entry": {
                        "row": i,
                        "col": j,
                        "value": "permutation breaks incidence",
                    },
                }
    return {"realized": True, "permutation": permutation, "signs": signs}


def require_signed_automorphism(
    response: RMat, adjacency: Sequence[Sequence[int]]
) -> dict[str, Any]:
    report = signed_automorphism_report(response, adjacency)
    require(
        report["realized"],
        "SIGNED_AUTOMORPHISM",
        "the response is not a signed graph automorphism: "
        + str(report.get("offending_entry")),
    )
    return report


def check_equivariance(
    response: RMat, perms: Sequence[tuple[int, ...]]
) -> None:
    for g in perms:
        for i in range(12):
            for j in range(12):
                require(
                    response[g[i]][g[j]] == response[i][j],
                    "PRESENTATION_INVARIANCE",
                    "conjugation by a carrier rotation moves an admissible response",
                )


def entry_value_texts(response: RMat) -> list[str]:
    values = {entry for row in response for entry in row}
    return [value.text() for value in sorted(values, key=f5_sort_key)]


# ---------------------------------------------------------------------------
# Manifest validation and the pinned producer reference
# ---------------------------------------------------------------------------


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    e565.enforce_source_firewall(manifest)
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    for key in (
        "mass_target",
        "measured_coupling",
        "particle_target",
        "coupling_target",
    ):
        require(key not in manifest, "FORBIDDEN_DEPENDENCY", f"forbidden manifest key {key}")
    require(
        tuple(manifest.get("grammar_clauses", ())) == GRAMMAR_CLAUSES,
        "GRAMMAR_DECLARATION",
        "the manifest must declare exactly the equivariant, involutive, signed, "
        "nontrivial grammar clauses in order",
    )
    require(
        manifest.get("operational_clause") == OPERATIONAL_CLAUSE,
        "GRAMMAR_DECLARATION",
        "the manifest must declare the signed-graph-automorphism readback clause",
    )
    require(
        manifest.get("declared_selected_response") == "minus_central_involution",
        "PRODUCER_SELECTION",
        "the declared selected response must be the producer convention R = -J",
    )
    require(
        manifest.get("lean_commutant_module") == LEAN_COMMUTANT_MODULE,
        "LEAN_REFERENCE",
        f"the manifest must cite {LEAN_COMMUTANT_MODULE} for the commutant dimension",
    )
    require(
        (MODULE_DIR.parents[1] / LEAN_COMMUTANT_MODULE).is_file(),
        "LEAN_REFERENCE",
        "the cited Lean commutant module is missing from the repository",
    )


def load_producer_manifest(
    manifest: Mapping[str, Any], base_dir: Path
) -> dict[str, Any]:
    path_raw = manifest.get("producer_manifest_path")
    require(
        isinstance(path_raw, str),
        "PRODUCER_REFERENCE",
        "producer_manifest_path is missing",
    )
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    producer = load_json(path)
    require(
        manifest.get("producer_manifest_sha256") == sha256_json(producer),
        "UPSTREAM_HASH",
        "the pinned #566 producer manifest hash does not match",
    )
    require(
        producer.get("schema") == PRODUCER_SCHEMA,
        "PRODUCER_REFERENCE",
        "the pinned manifest is not the #566 inverse-port response manifest",
    )
    blob = json.dumps(producer, sort_keys=True)
    require(
        "negative-antipode response" in blob,
        "PRODUCER_REFERENCE",
        "the pinned producer manifest does not carry the negative-antipode "
        "response convention",
    )
    return producer


# ---------------------------------------------------------------------------
# Certificate payload
# ---------------------------------------------------------------------------


def certificate_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    producer = load_producer_manifest(manifest, base)

    verts, adjacency, distance, antipode = port_model()
    perms = rotation_permutations(verts, adjacency)
    orbit_representatives = pair_orbits(perms, distance)
    projectors = isotypic_projectors(adjacency)
    frame_band_identification(verts, adjacency)
    j_mat = central_involution(adjacency, antipode, projectors)
    change_of_basis = orbital_basis_identities(adjacency, distance, j_mat, projectors)

    # --- Grammar enumeration: sixteen patterns, fourteen admissible ---------
    trivial = [pattern for pattern in ALL_PATTERNS if len(set(pattern)) == 1]
    admissible = validate_admissible_set(
        [pattern for pattern in ALL_PATTERNS if len(set(pattern)) == 2]
    )
    identity = rident(12)
    pattern_matrices: dict[SignPattern, RMat] = {}
    for pattern in ALL_PATTERNS:
        matrix = sign_pattern_matrix(projectors, pattern)
        require(
            mat_eq(p566.rmul(matrix, matrix), identity),
            "INVOLUTION_CHECK",
            f"sign pattern {pattern_name(pattern)} is not an involution",
        )
        pattern_matrices[pattern] = matrix
    require(
        mat_eq(pattern_matrices[(1, 1, 1, 1)], identity)
        and mat_eq(pattern_matrices[(-1, -1, -1, -1)], smul(MINUS_ONE, identity)),
        "INVOLUTION_CHECK",
        "the two trivial patterns are not +-identity",
    )

    # --- Presentation invariance for every admissible response --------------
    for pattern in admissible:
        check_equivariance(pattern_matrices[pattern], perms)

    # --- Operational selection: signed graph automorphisms ------------------
    realized: list[SignPattern] = []
    non_realized: list[dict[str, Any]] = []
    for pattern in admissible:
        report = signed_automorphism_report(pattern_matrices[pattern], adjacency)
        if report["realized"]:
            realized.append(pattern)
        else:
            non_realized.append(
                {
                    "pattern": pattern_name(pattern),
                    "offending_entry": report["offending_entry"],
                }
            )
    require(
        realized == [PATTERN_J, PATTERN_MINUS_J],
        "SIGNED_AUTOMORPHISM_SET",
        f"expected exactly +-J as signed automorphisms, got {realized}",
    )
    require(
        len(non_realized) == 12,
        "SIGNED_AUTOMORPHISM_SET",
        "twelve admissible patterns must fail the signed-automorphism clause",
    )
    j_report = require_signed_automorphism(pattern_matrices[PATTERN_J], adjacency)
    minus_j_report = require_signed_automorphism(
        pattern_matrices[PATTERN_MINUS_J], adjacency
    )
    require(
        j_report["permutation"] == list(antipode)
        and set(j_report["signs"]) == {1}
        and minus_j_report["permutation"] == list(antipode)
        and set(minus_j_report["signs"]) == {-1},
        "SIGNED_AUTOMORPHISM_SET",
        "the realized responses are not +-1 times the antipode permutation",
    )

    # --- Producer selection --------------------------------------------------
    selected = pattern_matrices[PATTERN_MINUS_J]
    require(
        mat_eq(selected, smul(MINUS_ONE, j_mat)),
        "PRODUCER_SELECTION",
        "the selected producer response is not -J",
    )

    # --- Countermodel pair: the operational clause is load-bearing ----------
    model_a = pattern_matrices[PATTERN_J]
    model_b = pattern_matrices[PATTERN_FIVE_FLIP]
    require(
        not mat_eq(model_a, model_b)
        and not mat_eq(model_a, smul(MINUS_ONE, model_b)),
        "COUNTERMODEL",
        "the countermodel pair is equivalent up to sign",
    )
    b_report = signed_automorphism_report(model_b, adjacency)
    require(
        b_report["realized"] is False,
        "COUNTERMODEL",
        "the second countermodel unexpectedly passes the operational clause",
    )
    countermodel_pair = {
        "clause_dropped": OPERATIONAL_CLAUSE,
        "model_a": {
            "pattern": pattern_name(PATTERN_J),
            "sector_signs": dict(zip(SECTORS, PATTERN_J)),
            "trace": rtrace(model_a).text(),
            "entry_values": entry_value_texts(model_a),
            "signed_graph_automorphism": True,
            "description": "the central involution J, the antipode permutation",
        },
        "model_b": {
            "pattern": pattern_name(PATTERN_FIVE_FLIP),
            "sector_signs": dict(zip(SECTORS, PATTERN_FIVE_FLIP)),
            "trace": rtrace(model_b).text(),
            "entry_values": entry_value_texts(model_b),
            "signed_graph_automorphism": False,
            "offending_entry": b_report["offending_entry"],
            "description": "I - 2 P_5, the reflection through the quintet sector",
        },
        "both_satisfy_equivariance_involution_signedness_nontriviality": True,
        "inequivalent": True,
        "conclusion": (
            "incidence with equivariance and involution alone retains fourteen "
            "admissible responses, two exhibited here as an inequivalent pair, "
            "so the operational readback clause is load-bearing for the "
            "uniqueness of the inverse-port law"
        ),
    }

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 611,
        "manifest_sha256": sha256_json(manifest),
        "producer_reference": {
            "manifest_path": str(manifest.get("producer_manifest_path")),
            "manifest_sha256": sha256_json(producer),
            "schema": PRODUCER_SCHEMA,
            "convention": (
                "the incidence-derived negative-antipode response of the #566 "
                "inverse-port producer; the common sign is a charge-conjugation "
                "convention"
            ),
        },
        "port_space": {
            "dimension": 12,
            "arithmetic_field": "Q(sqrt5), no floating point in any proof decision",
            "vertex_model": "cyclic permutations of (0, +-1, +-phi)",
            "adjacency_eigenvalues": {
                sector: EIGENVALUES[sector].text() for sector in SECTORS
            },
            "isotypic_decomposition": "1 + 3 + 3' + 5",
            "frame_band_identification": (
                "neighbour sums equal sqrt5 times the vertex coordinates, so "
                "sector 3 carries the vertex frame"
            ),
        },
        "commutant": {
            "dimension": 4,
            "pair_orbit_count": 4,
            "pair_orbit_representatives": orbit_representatives,
            "lean_reference": LEAN_COMMUTANT_MODULE,
            "lean_theorems": [
                "commutant_decomposition",
                "orbitals_independent",
                "equivariant_iff_commutes",
            ],
            "lean_statement": (
                "every matrix equivariant under the listed sixty rotations is a "
                "unique rational combination of the identity, adjacency, "
                "distance-two, and antipode orbitals; entry equivariance is "
                "equivalent to commutation with the linear action; dimension "
                "four matches the rational isotypic structure 1 + 5 + (3 + 3') "
                "with endomorphism fields Q, Q, Q(sqrt5)"
            ),
            "orbital_in_projector_basis": change_of_basis,
            "basis_note": (
                "the orbital-to-projector coefficient matrix is a Vandermonde "
                "system in the four distinct eigenvalues, hence invertible: the "
                "projector quadruple spans exactly the machine-checked "
                "four-dimensional commutant"
            ),
        },
        "isotypic_projectors": {
            "construction": (
                "Lagrange interpolation in the adjacency matrix over the four "
                "exact eigenvalues"
            ),
            "idempotent": True,
            "mutually_orthogonal": True,
            "partition_of_unity": True,
            "eigenvalue_relations_verified": True,
            "ranks_by_exact_trace": {
                sector: SECTOR_RANKS[sector] for sector in SECTORS
            },
            "diagonal_entries": {
                sector: projectors[sector][0][0].text() for sector in SECTORS
            },
        },
        "involution_completeness": {
            "argument": (
                "an equivariant involution lies in the four-dimensional "
                "commutant; the projectors are an exact basis of orthogonal "
                "idempotents summing to the identity, so R = sum(c_s P_s) has "
                "R^2 = sum(c_s^2 P_s) and R^2 = I forces c_s^2 = 1 on every "
                "sector; on each real-irreducible isotypic block the "
                "equivariant commutant is the real scalars by Schur, so every "
                "coefficient is +-1"
            ),
            "pattern_count": 16,
            "trivial_patterns": [pattern_name(pattern) for pattern in trivial],
            "admissible_count": 14,
            "admissible_patterns": [
                pattern_name(pattern) for pattern in admissible
            ],
            "sector_order": list(SECTORS),
        },
        "signed_automorphism_selection": {
            "realized_patterns": [pattern_name(p) for p in realized],
            "realized_count": 2,
            "central_involution_identity": "J = (A^3 - 4*A^2 - 5*A + 10*I)/10",
            "identity_verified": True,
            "J_sign_pattern": list(PATTERN_J),
            "J_is_antipode_permutation": True,
            "non_realized": non_realized,
        },
        "selected_producer_response": {
            "pattern": list(PATTERN_MINUS_J),
            "pattern_name": pattern_name(PATTERN_MINUS_J),
            "response": "R = -J, minus the antipode permutation",
            "convention": (
                "impulse/readback sign convention of the pinned #566 "
                "inverse-port producer"
            ),
            "realized_as_signed_automorphism": True,
        },
        "countermodel_pair": countermodel_pair,
        "presentation_invariance": {
            "rotations_checked": 60,
            "group_closure_verified": True,
            "conjugation_fixes_every_admissible_response": True,
            "refinement_note": (
                "each admissible response is a polynomial in the adjacency "
                "matrix, so it is fixed by conjugation with every incidence "
                "automorphism; the sixty rotations are checked exhaustively"
            ),
            "verdict": "presentation_invariant",
        },
        "verdict": {
            "grammar_internal_uniqueness": "up_to_sign",
            "axiom_forcing": "independence_limited",
            "operational_clause": "load_bearing",
            "bounded_exit": "independence_limited",
        },
        "claim_boundary": {
            "proves": (
                "the complete mechanism space of the declared quotient-visible "
                "response grammar (sixteen sign patterns, fourteen admissible), "
                "the exactly-two signed-automorphism realizations +-J, the "
                "producer selection R = -J unique up to the conventional sign "
                "inside the operational grammar, the load-bearing status of the "
                "operational readback clause via the retained countermodel "
                "pair, and presentation invariance under the sixty carrier "
                "rotations"
            ),
            "status": "exact_named_realization",
            "does_not_close": [
                "forcing the equivariance, involution, and signed-readback "
                "grammar itself from the three axioms",
                "laboratory attachment of the response",
                "continuum dynamics",
            ],
        },
        "verifier_command": (
            "python3 code/a5_closure/response_grammar_completeness_certificate.py "
            "verify --manifest "
            "code/a5_closure/manifests/response_grammar_completeness_reference.json "
            "--receipt "
            "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json"
        ),
    }


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def tampered_adjacency_control() -> None:
    """Breaking one incidence edge must break the projector system."""

    _, adjacency, _, _ = port_model()
    broken = [row[:] for row in adjacency]
    neighbour = next(j for j in range(12) if broken[0][j])
    broken[0][neighbour] = 0
    broken[neighbour][0] = 0
    isotypic_projectors(broken)


def fifteenth_pattern_control() -> None:
    """A claimed fifteenth admissible pattern must be rejected."""

    claimed = [pattern for pattern in ALL_PATTERNS if len(set(pattern)) == 2]
    claimed.append((1, 1, 1, 1))
    validate_admissible_set(claimed)


def non_automorphism_promotion_control() -> None:
    """Promoting the quintet-flip response to an automorphism must fail the
    entry check."""

    _, adjacency, _, _ = port_model()
    projectors = isotypic_projectors(adjacency)
    require_signed_automorphism(
        sign_pattern_matrix(projectors, PATTERN_FIVE_FLIP), adjacency
    )


def negative_control_cases(
    manifest: Mapping[str, Any]
) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    wrong_schema = copy.deepcopy(dict(manifest))
    wrong_schema["schema"] = "oph.response_grammar_completeness_certificate.v0"
    cases.append(("wrong_schema", wrong_schema, "SCHEMA"))

    clause_dropped = copy.deepcopy(dict(manifest))
    clause_dropped["grammar_clauses"] = ["equivariant", "signed", "nontrivial"]
    cases.append(("involution_clause_dropped", clause_dropped, "GRAMMAR_DECLARATION"))

    operational_dropped = copy.deepcopy(dict(manifest))
    operational_dropped["operational_clause"] = "none"
    cases.append(("operational_clause_dropped", operational_dropped, "GRAMMAR_DECLARATION"))

    plus_j = copy.deepcopy(dict(manifest))
    plus_j["declared_selected_response"] = "plus_central_involution"
    cases.append(("plus_J_promotion", plus_j, "PRODUCER_SELECTION"))

    pin_drift = copy.deepcopy(dict(manifest))
    pin_drift["producer_manifest_sha256"] = "0" * 64
    cases.append(("producer_pin_drift", pin_drift, "UPSTREAM_HASH"))

    target_injection = copy.deepcopy(dict(manifest))
    target_injection["mass_target"] = {"target": "lepton"}
    cases.append(("mass_target_injection", target_injection, "FORBIDDEN_DEPENDENCY"))

    firewall = copy.deepcopy(dict(manifest))
    firewall["description"] = "hypercharge assignment for the response"
    cases.append(("source_firewall_token", firewall, "FORBIDDEN_DEPENDENCY"))

    return cases


def negative_control_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    results = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir)
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append(
            {
                "name": name,
                "expected_error": expected_code,
                "actual_error": actual_code,
                "passed": True,
            }
        )
    tamper_cases: list[tuple[str, Callable[[], None], str]] = [
        ("tampered_adjacency_breaks_projectors", tampered_adjacency_control, "PROJECTOR_SYSTEM"),
        ("fifteenth_admissible_pattern_rejected", fifteenth_pattern_control, "ADMISSIBLE_SET"),
        ("non_automorphism_promotion_rejected", non_automorphism_promotion_control, "SIGNED_AUTOMORPHISM"),
    ]
    for name, control, expected_code in tamper_cases:
        actual_code = "ACCEPTED"
        try:
            control()
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append(
            {
                "name": name,
                "expected_error": expected_code,
                "actual_error": actual_code,
                "passed": True,
            }
        )
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 611,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "grammar_without_operational_clause": (
                "fourteen admissible responses survive; the retained pair "
                "(+,-,-,+) and (+,+,+,-) is inequivalent, so the readback "
                "clause carries the uniqueness"
            ),
            "tampered_adjacency": (
                "one broken incidence edge destroys the quartic minimal "
                "polynomial and the projector system fails closed"
            ),
            "fifteenth_pattern": (
                "the +-identity patterns are the only two excluded patterns; a "
                "claimed fifteenth admissible pattern is rejected"
            ),
            "non_automorphism_promotion": (
                "the quintet-flip response carries the dense entry 1/6 on the "
                "diagonal and fails the signed-automorphism entry check"
            ),
        },
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def verify_receipt(
    manifest: Mapping[str, Any],
    receipt: Mapping[str, Any],
    base_dir: Path | None = None,
) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "response_grammar_completeness_reference.json",
        MODULE_DIR / "receipts" / "response_grammar_completeness_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_611_negative_controls.json",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all")
    all_cmd.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS"}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "receipt": str(default_receipt),
                    "negative_controls": str(default_negative),
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
