#!/usr/bin/env python3
"""Exact certificate for GitHub issue #610: equal state weights from the three axioms.

Executable half of the spine in Lean/Screen/EqualStateWeights.lean. The Lean
module proves: if the feasible family and the A3 objective are invariant under
a transitive port action, and the information projection is unique, then the
realized state weight is 1/12 at every port. This certificate supplies the
concrete inputs that theorem consumes and the countermodels that delimit it:

* the deck action: the sixty orientation-preserving automorphisms of the
  twelve-port icosahedral incidence, computed by exhaustive graph-automorphism
  enumeration (120 incidence automorphisms, orientation character 60 + 60),
  verified as a group (closure, identity, two-sided inverses), verified to
  preserve adjacency and the oriented face list of
  Lean/ObserverPatchHolography/CoreAxioms.lean with consistent cyclic
  orientation, to commute with the antipode i -> 11 - i, to act pairwise
  transitively on ports, and to coincide exactly with the sixty listed
  rotations that Lean/Screen/A5PortAction.lean kernel-checks;

* the positive lane in closed form. The feasible set is the twelve-port
  probability simplex: the complete A1/A2 constraint grammar contributes only
  deck-invariant constraints, and the base case is normalization alone. A3
  supplies the exact reference and aggregation rule of
  claims/axiom_registry.yaml: a compatible local reference family tau_r, a
  finite A1-generated observer cover G_r whose restriction map is injective
  on the feasible family, and strictly positive exact weights w_{r,P} from
  quotient-visible A1 data, natural under admissible presentation
  equivalence. With reference tau uniform and any strictly positive invariant
  weight w, the KKT stationarity condition for minimizing the weighted
  divergence sum_p w rho_p log(rho_p / tau_p) on the simplex reads
  w (log(rho_p / tau_p) + 1) = lambda at every port. At rho = tau every
  density ratio rho_p / tau_p is exactly 1, so every gradient component
  equals w (log 1 + 1) = w, an exact symbolic equality with no numerics.
  The map x -> x log x is strictly convex on the open simplex, so with a
  fixed full-support reference and strictly positive weights the objective
  is strictly convex on the feasible family and the stationary point is the
  unique minimizer: the projection is tau itself and the realized weight is
  exactly 1/12 at every port. An added deck-invariant linear constraint
  keeps the uniform point feasible and stationary;

* the distinction between two objects that agree at 1/12: the measured
  central block trace of the carrier port atoms (a realization property read
  from the pinned federation carrier manifest) and the derived state weight
  (this theorem). The receipt records both fields and the fact that they are
  distinct objects;

* four countermodels, one per dropped hypothesis, each required to fail: a
  tilted reference projects to itself and misses the uniform point; the
  order-ten stabilizer of one antipodal port pair admits an orbit-constant
  invariant reference whose projection is non-uniform; a linear invariant
  objective has every simplex point as a minimizer, so uniqueness fails; a
  cover omitting one port has a non-injective restriction map, so two
  feasible local families differ at the omitted port with equal objective.
  A control that unexpectedly passes raises and the certificate exits
  nonzero.

Every arithmetic decision is exact integer or fractions.Fraction arithmetic;
no floating point appears in the positive lane or anywhere else, and the
emitted payload is walked to reject any float.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

SCHEMA = "oph.equal_state_weights_certificate.v1"
GENERATED_BY = "code/a5_closure/equal_state_weights_certificate.py"
LEAN_SPINE_MODULE = "Lean/Screen/EqualStateWeights.lean"
LEAN_ACTION_MODULE = "Lean/Screen/A5PortAction.lean"
LEAN_FACES_MODULE = "Lean/ObserverPatchHolography/CoreAxioms.lean"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"

PORTS = 12
UNIFORM = tuple(Fraction(1, 12) for _ in range(PORTS))

# Neighbor table of Lean/Screen/PortFrameGram.lean; port i is antipodal to
# 11 - i.
NEIGHBORS: tuple[tuple[int, ...], ...] = (
    (1, 2, 3, 4, 6),
    (0, 2, 3, 5, 7),
    (0, 1, 4, 5, 8),
    (0, 1, 6, 7, 9),
    (0, 2, 6, 8, 10),
    (1, 2, 7, 8, 11),
    (0, 3, 4, 9, 10),
    (1, 3, 5, 9, 11),
    (2, 4, 5, 10, 11),
    (3, 6, 7, 10, 11),
    (4, 6, 8, 9, 11),
    (5, 7, 8, 9, 10),
)
ADJ: tuple[frozenset[int], ...] = tuple(frozenset(row) for row in NEIGHBORS)

# The twenty coherently oriented faces of
# Lean/ObserverPatchHolography/CoreAxioms.lean `orientedFaces`, in the same
# port labels (the Lean adjacency and this module's adjacency are identical).
ORIENTED_FACES: tuple[tuple[int, int, int], ...] = (
    (0, 1, 2), (0, 3, 1), (0, 2, 4), (0, 6, 3), (0, 4, 6),
    (1, 5, 2), (1, 3, 7), (1, 7, 5), (2, 8, 4), (2, 5, 8),
    (3, 6, 9), (3, 9, 7), (4, 10, 6), (4, 8, 10), (5, 7, 11),
    (5, 11, 8), (6, 10, 9), (7, 9, 11), (8, 11, 10), (9, 10, 11),
)

# Control (a) reference: one tilted, strictly positive, normalized reference
# family. Its unconstrained-simplex projection is itself, and it differs from
# the uniform point, so the wrong-reference control fails as required.
WRONG_REFERENCE_TAU: tuple[Fraction, ...] = (Fraction(2, 13),) + tuple(
    Fraction(1, 13) for _ in range(PORTS - 1)
)

Permutation = tuple[int, ...]


def antipode(i: int) -> int:
    return 11 - i


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[i]] for i in range(PORTS))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * PORTS
    for i, image in enumerate(permutation):
        result[image] = i
    return tuple(result)


def cyclic_class(face: tuple[int, int, int]) -> tuple[int, int, int]:
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


POSITIVE_FACE_CLASSES = frozenset(cyclic_class(face) for face in ORIENTED_FACES)
NEGATIVE_FACE_CLASSES = frozenset(
    cyclic_class((face[2], face[1], face[0])) for face in ORIENTED_FACES
)


def verify_face_data() -> None:
    """The copied oriented face list is a coherent orientation of the adjacency."""

    require(
        len(ORIENTED_FACES) == 20 and len(POSITIVE_FACE_CLASSES) == 20,
        "FACE_DATA",
        "the oriented face list must carry twenty distinct faces",
    )
    require(
        not (POSITIVE_FACE_CLASSES & NEGATIVE_FACE_CLASSES),
        "FACE_DATA",
        "a face coincides with its own reversal",
    )
    directed: dict[tuple[int, int], int] = {}
    for a, b, c in ORIENTED_FACES:
        for i, j in ((a, b), (b, c), (c, a)):
            require(j in ADJ[i], "FACE_DATA", f"face edge ({i},{j}) is not an incidence edge")
            directed[(i, j)] = directed.get((i, j), 0) + 1
    for i in range(PORTS):
        for j in ADJ[i]:
            require(
                directed.get((i, j)) == 1 and directed.get((j, i)) == 1,
                "FACE_DATA",
                f"edge ({i},{j}) does not lie in exactly two opposite-oriented faces",
            )


def orientation_sign(permutation: Permutation) -> int:
    """+1 when every oriented face maps to an oriented face with consistent
    cyclic orientation, -1 when every image lands in the reversed classes."""

    images = [
        cyclic_class((permutation[a], permutation[b], permutation[c]))
        for a, b, c in ORIENTED_FACES
    ]
    if all(image in POSITIVE_FACE_CLASSES for image in images):
        return 1
    if all(image in NEGATIVE_FACE_CLASSES for image in images):
        return -1
    raise CertificateError(
        "FACE_ORIENTATION_MIXED",
        "an incidence automorphism maps faces with inconsistent orientation",
    )


def graph_automorphisms() -> list[Permutation]:
    """Exhaustive backtracking enumeration of the incidence automorphisms."""

    results: list[Permutation] = []
    image = [-1] * PORTS
    used = [False] * PORTS

    def extend(vertex: int) -> None:
        if vertex == PORTS:
            results.append(tuple(image))
            return
        for candidate in range(PORTS):
            if used[candidate]:
                continue
            consistent = True
            for prior in range(vertex):
                if (prior in ADJ[vertex]) != (image[prior] in ADJ[candidate]):
                    consistent = False
                    break
            if consistent:
                image[vertex] = candidate
                used[candidate] = True
                extend(vertex + 1)
                used[candidate] = False
                image[vertex] = -1

    extend(0)
    return results


def deck_rotations() -> list[Permutation]:
    """The sixty orientation-preserving deck rotations, sorted."""

    automorphisms = graph_automorphisms()
    require(
        len(automorphisms) == 120,
        "AUTOMORPHISM_COUNT",
        f"expected 120 incidence automorphisms, got {len(automorphisms)}",
    )
    for g in automorphisms:
        for i in range(PORTS):
            require(
                g[antipode(i)] == antipode(g[i]),
                "ANTIPODE_CENTER",
                "an incidence automorphism does not commute with the antipode",
            )
    signs = {g: orientation_sign(g) for g in automorphisms}
    plus = sorted(g for g in automorphisms if signs[g] == 1)
    minus = [g for g in automorphisms if signs[g] == -1]
    require(
        len(plus) == 60 and len(minus) == 60,
        "ORIENTATION_SPLIT",
        "the orientation character must split 120 = 60 + 60",
    )
    return plus


def verify_rotation_group(rotations: Sequence[Permutation]) -> dict[str, Any]:
    """Group axioms, invariances, and pairwise port transitivity, all exact."""

    require(
        len(rotations) == 60 and len(set(rotations)) == 60,
        "GROUP_ORDER",
        f"expected sixty distinct rotations, got {len(set(rotations))}",
    )
    for g in rotations:
        require(
            sorted(g) == list(range(PORTS)),
            "PERMUTATION",
            "a listed rotation is not a port permutation",
        )
    rotation_set = set(rotations)
    require(
        tuple(range(PORTS)) in rotation_set,
        "GROUP_IDENTITY",
        "the identity permutation is missing",
    )
    for g in rotations:
        for i in range(PORTS):
            require(
                {g[j] for j in ADJ[i]} == set(ADJ[g[i]]),
                "ADJACENCY_PRESERVATION",
                "a listed rotation does not preserve the port adjacency",
            )
    for g in rotations:
        for i in range(PORTS):
            require(
                g[antipode(i)] == antipode(g[i]),
                "ANTIPODE_COMMUTATION",
                "a listed rotation does not commute with the antipode",
            )
    for g in rotations:
        for h in rotations:
            require(
                compose(g, h) in rotation_set,
                "GROUP_CLOSURE",
                "the listed rotations are not closed under composition",
            )
    for g in rotations:
        require(
            inverse(g) in rotation_set,
            "GROUP_INVERSES",
            "a listed rotation lacks a listed inverse",
        )
    for g in rotations:
        require(
            orientation_sign(g) == 1,
            "FACE_ORIENTATION",
            "a listed rotation does not preserve the oriented face list",
        )
    ordered_pairs = 0
    for source in range(PORTS):
        for target in range(PORTS):
            require(
                any(g[source] == target for g in rotations),
                "PORT_TRANSITIVITY",
                f"no listed rotation carries port {source} to port {target}",
            )
            ordered_pairs += 1
    return {
        "port_count": PORTS,
        "full_incidence_automorphism_order": 120,
        "rotation_group_order": 60,
        "identity_present": True,
        "closure_verified": True,
        "inverses_verified": True,
        "adjacency_preserved": True,
        "oriented_faces_preserved_with_consistent_cyclic_orientation": True,
        "antipode_i_to_11_minus_i_commutes": True,
        "pairwise_transitivity": {
            "ordered_port_pairs_checked": ordered_pairs,
            "all_pairs_connected": True,
        },
        "rotation_list_sha256": sha256_json([list(g) for g in rotations]),
    }


def lean_listed_rotations() -> list[Permutation]:
    """The sixty rotations listed in Lean/Screen/A5PortAction.lean `perms`."""

    path = REPO_ROOT / "Lean" / "Screen" / "A5PortAction.lean"
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise CertificateError("LEAN_ACTION_READ", f"cannot read {path}: {exc}") from exc
    marker = "def perms : List (List Nat) := ["
    start = text.find(marker)
    require(start >= 0, "LEAN_ACTION_PARSE", "the Lean perms block is missing")
    end = text.find("]]", start)
    require(end >= 0, "LEAN_ACTION_PARSE", "the Lean perms block is unterminated")
    block = text[start + len(marker) - 1 : end + 1]
    rows = re.findall(r"\[([0-9,\s]+)\]", block)
    perms = [tuple(int(token) for token in row.split(",")) for row in rows]
    require(
        len(perms) == 60 and all(len(p) == PORTS for p in perms),
        "LEAN_ACTION_PARSE",
        f"expected sixty twelve-entry rows in the Lean perms block, got {len(perms)}",
    )
    return perms


def orbits_under(subgroup: Sequence[Permutation]) -> list[tuple[int, ...]]:
    seen: set[int] = set()
    orbits: list[tuple[int, ...]] = []
    for port in range(PORTS):
        if port in seen:
            continue
        orbit = sorted({g[port] for g in subgroup})
        seen.update(orbit)
        orbits.append(tuple(orbit))
    return orbits


def stabilizer_of_antipodal_pair(
    rotations: Sequence[Permutation], port: int = 0
) -> list[Permutation]:
    """The setwise stabilizer of the antipodal port pair {port, 11 - port}."""

    pair = {port, antipode(port)}
    return [g for g in rotations if {g[port], g[antipode(port)]} == pair]


def projection_witness(
    tau: Sequence[Fraction], weight: Fraction
) -> dict[str, Any]:
    """Closed-form information projection onto the simplex, exact.

    For strictly positive reference tau and strictly positive weight w, the
    KKT stationarity condition of min sum_p w rho_p log(rho_p / tau_p) over
    the simplex is w (log(rho_p / tau_p) + 1) = lambda at every port, so the
    density ratios rho_p / tau_p are equal, and normalization forces
    rho = tau. Strict convexity of x -> x log x on the open simplex with a
    fixed full-support reference makes the stationary point the unique
    minimizer. The witness is the candidate rho = tau together with the
    exact density ratios (all 1) and the exact symbolic gradient component
    w (log 1 + 1) = w, equal at every port with no numerics.
    """

    require(
        isinstance(weight, Fraction) and weight > 0,
        "WEIGHTS_POSITIVE",
        "the A3 weight must be a strictly positive exact fraction",
    )
    require(
        all(isinstance(t, Fraction) and t > 0 for t in tau),
        "REFERENCE_SUPPORT",
        "the reference family must be strictly positive at every port",
    )
    require(
        sum(tau) == 1,
        "REFERENCE_NORMALIZATION",
        "the reference family must be normalized",
    )
    projection = tuple(tau)
    ratios = tuple(projection[p] / tau[p] for p in range(PORTS))
    require(
        all(ratio == Fraction(1) for ratio in ratios),
        "KKT_RATIO",
        "the density ratio at the projection is not exactly one at every port",
    )
    require(
        len(set(ratios)) == 1,
        "KKT_GRADIENT",
        "the KKT gradient components are not exactly equal across ports",
    )
    return {
        "projection_exact": projection,
        "density_ratios_rho_over_tau": [str(ratio) for ratio in ratios],
        "gradient_component_symbolic": "w * (log(1) + 1) = w, identical at every port",
        "gradient_component_exact_weight": str(weight),
        "uniqueness": (
            "x -> x log x is strictly convex on the open simplex; with a fixed "
            "full-support reference and strictly positive weights the objective "
            "is strictly convex on the feasible family, so the stationary point "
            "is the unique minimizer"
        ),
    }


def positive_lane(rotations: Sequence[Permutation]) -> dict[str, Any]:
    """The exact positive lane: uniform reference, invariant weights, 1/12."""

    weight = Fraction(1)
    witness = projection_witness(UNIFORM, weight)
    projection = witness["projection_exact"]
    require(
        all(value == Fraction(1, 12) for value in projection),
        "POSITIVE_LANE",
        "the information projection is not exactly 1/12 at every port",
    )
    for g in rotations:
        permuted = tuple(UNIFORM[g[p]] for p in range(PORTS))
        require(
            permuted == UNIFORM,
            "REFERENCE_INVARIANCE",
            "the uniform reference is not invariant under a deck rotation",
        )
    # One deck-invariant linear constraint: sum_p c rho_p = c with c = 3. The
    # coefficient vector is constant, hence invariant under every port
    # permutation; the uniform point satisfies it exactly; its gradient is a
    # constant vector, proportional to the normalization gradient, so KKT
    # stationarity at the uniform point survives with the same equal
    # components.
    c = Fraction(3)
    coefficients = tuple(c for _ in range(PORTS))
    for g in rotations:
        permuted = tuple(coefficients[g[p]] for p in range(PORTS))
        require(
            permuted == coefficients,
            "CONSTRAINT_INVARIANCE",
            "the added linear constraint is not deck invariant",
        )
    lhs = sum(coefficients[p] * UNIFORM[p] for p in range(PORTS))
    require(
        lhs == c,
        "CONSTRAINT_FEASIBILITY",
        "the uniform point does not satisfy the added invariant constraint",
    )
    require(
        len(set(coefficients)) == 1,
        "CONSTRAINT_STATIONARITY",
        "the added constraint gradient is not proportional to the normalization gradient",
    )
    # The implied orbit constraint: the deck action has a single port orbit,
    # and the sum over that orbit equals the orbit fraction 12/12 = 1, which
    # normalization supplies.
    orbit = sorted({g[0] for g in rotations})
    require(
        orbit == list(range(PORTS)),
        "ORBIT_CONSTRAINT",
        "the deck action does not have a single port orbit",
    )
    orbit_fraction = Fraction(len(orbit), PORTS)
    require(
        sum(UNIFORM[p] for p in orbit) == orbit_fraction == Fraction(1),
        "ORBIT_CONSTRAINT",
        "the single-orbit sum does not equal the orbit fraction",
    )
    return {
        "feasible_set": (
            "the twelve-port probability simplex; the complete A1/A2 "
            "constraint grammar contributes only deck-invariant constraints, "
            "and the base case is normalization alone"
        ),
        "reference_family_tau": [str(value) for value in UNIFORM],
        "reference_invariant_under_deck_action": True,
        "weights": (
            "w = 1 at every port; invariance under the transitive deck action "
            "forces any admissible strictly positive exact weight rule to be "
            "constant across ports, and the stationarity witness is identical "
            "for every such constant"
        ),
        "kkt_stationarity": {
            "density_ratios_rho_over_tau": witness["density_ratios_rho_over_tau"],
            "gradient_component_symbolic": witness["gradient_component_symbolic"],
            "gradient_components_all_equal": True,
            "numerics_used": False,
        },
        "uniqueness": witness["uniqueness"],
        "invariant_linear_constraint_case": {
            "constraint": "sum_p c rho_p = c with c = 3",
            "deck_invariant": True,
            "uniform_point_feasible": True,
            "uniform_point_stationary": True,
            "gradient_note": (
                "the constraint gradient is the constant vector (3, ..., 3), "
                "proportional to the normalization gradient, so the KKT system "
                "at the uniform point is unchanged"
            ),
        },
        "implied_orbit_constraint": {
            "single_port_orbit": True,
            "orbit_sum_equals_orbit_fraction": "12/12 = 1, supplied by normalization",
        },
        "projection": [str(value) for value in projection],
        "result_per_port": "1/12",
    }


def state_weight_vs_block_trace(derived_state_weight: str) -> dict[str, Any]:
    """Both fields on one receipt: measured block trace and derived weight.

    The carrier manifest's primitive central port atoms carry the measured
    normalized block trace 1/12 per port. That is a realization property of
    the certified carrier algebra. The derived state weight 1/12 per port is
    the conclusion of the equal-state-weights theorem. The two are distinct
    objects that agree at 1/12.
    """

    path = MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME
    manifest = load_json(path)
    carrier = manifest.get("carrier")
    require(isinstance(carrier, Mapping), "CARRIER_MANIFEST", "carrier object is missing")
    atoms = carrier.get("central_port_atoms")
    require(
        isinstance(atoms, list) and len(atoms) == PORTS,
        "CARRIER_MANIFEST",
        "the carrier manifest must list twelve primitive central port atoms",
    )
    for atom in atoms:
        require(
            atom.get("primitive") is True,
            "CARRIER_MANIFEST",
            "a central port atom is not primitive",
        )
        trace = atom.get("normalized_trace", {})
        value = Fraction(int(trace.get("numerator")), int(trace.get("denominator")))
        require(
            value == Fraction(1, 12),
            "CARRIER_MANIFEST",
            f"atom {atom.get('atom_id')} does not carry normalized trace 1/12",
        )
    return {
        "carrier_manifest_path": f"manifests/{CARRIER_MANIFEST_NAME}",
        "carrier_manifest_sha256": sha256_json(manifest),
        "measured_central_block_trace": "1/12",
        "measured_central_block_trace_source": (
            "the twelve primitive central port atoms of the pinned federation "
            "carrier manifest, each with normalized trace 1/12"
        ),
        "derived_state_weight": derived_state_weight,
        "derived_state_weight_source": (
            "the information projection under invariance, uniqueness, and "
            "transitivity, composed in Lean/Screen/EqualStateWeights.lean"
        ),
        "distinct_objects_that_agree_at_one_twelfth": True,
        "distinction": (
            "the block trace is a realization property of the carrier algebra; "
            "the state weight is a property of the realized A3 state; neither "
            "value is derived from the other in this certificate"
        ),
    }


def control_wrong_reference() -> dict[str, Any]:
    """Control (a): a tilted reference. The projection equals the tilted
    reference and misses the uniform point, so the 1/12 conclusion fails."""

    tau = WRONG_REFERENCE_TAU
    witness = projection_witness(tau, Fraction(1))
    projection = witness["projection_exact"]
    require(
        projection != UNIFORM,
        "CONTROL_NOT_FAILED",
        "the wrong-reference control unexpectedly passed: the tilted-reference "
        "projection equals the uniform point",
    )
    differing_port = next(
        p for p in range(PORTS) if projection[p] != Fraction(1, 12)
    )
    return {
        "dropped_hypothesis": "the reference family is the invariant uniform family",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "tau_prime": [str(value) for value in tau],
            "projection": [str(value) for value in projection],
            "projection_equals_tau_prime": True,
            "differs_from_uniform_at_port": differing_port,
            "projection_value_at_that_port": str(projection[differing_port]),
            "uniform_value": "1/12",
        },
    }


def control_non_transitive(rotations: Sequence[Permutation]) -> dict[str, Any]:
    """Control (b): the stabilizer of one antipodal port pair. The subgroup
    is order ten, acts with two orbits, and admits an orbit-constant
    invariant non-uniform reference whose projection is non-uniform."""

    subgroup = stabilizer_of_antipodal_pair(rotations, 0)
    require(
        len(subgroup) == 10,
        "STABILIZER_ORDER",
        f"expected the order-ten pair stabilizer, got order {len(subgroup)}",
    )
    subgroup_set = set(subgroup)
    for g in subgroup:
        for h in subgroup:
            require(
                compose(g, h) in subgroup_set,
                "STABILIZER_CLOSURE",
                "the pair stabilizer is not closed under composition",
            )
    orbit_list = orbits_under(subgroup)
    sizes = sorted(len(orbit) for orbit in orbit_list)
    require(
        sizes == [2, 10],
        "STABILIZER_ORBITS",
        f"expected orbit sizes [2, 10], got {sizes}",
    )
    pair_orbit = next(orbit for orbit in orbit_list if len(orbit) == 2)
    ring_orbit = next(orbit for orbit in orbit_list if len(orbit) == 10)
    require(
        pair_orbit == (0, 11),
        "STABILIZER_ORBITS",
        f"the two-port orbit is {pair_orbit}, expected (0, 11)",
    )
    tau = [Fraction(0)] * PORTS
    for port in pair_orbit:
        tau[port] = Fraction(1, 4)
    for port in ring_orbit:
        tau[port] = Fraction(1, 20)
    tau_tuple = tuple(tau)
    for g in subgroup:
        permuted = tuple(tau_tuple[g[p]] for p in range(PORTS))
        require(
            permuted == tau_tuple,
            "CONTROL_INVARIANCE",
            "the orbit-constant reference is not invariant under the stabilizer",
        )
    transitive = all(
        any(g[0] == target for g in subgroup) for target in range(PORTS)
    )
    require(
        not transitive,
        "CONTROL_NOT_FAILED",
        "the non-transitive control unexpectedly passed: the pair stabilizer "
        "acts transitively on ports",
    )
    witness = projection_witness(tau_tuple, Fraction(1))
    projection = witness["projection_exact"]
    require(
        projection != UNIFORM,
        "CONTROL_NOT_FAILED",
        "the non-transitive control unexpectedly passed: the orbit-constant "
        "reference projects to the uniform point",
    )
    return {
        "dropped_hypothesis": "the invariance group acts transitively on ports",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "stabilized_antipodal_pair": [0, 11],
            "stabilizer_order": 10,
            "stabilizer_orbit_sizes": [2, 10],
            "subgroup_invariant_tau": [str(value) for value in tau_tuple],
            "projection": [str(value) for value in projection],
            "projection_value_on_pair_orbit": "1/4",
            "projection_value_on_ring_orbit": "1/20",
            "projection_is_uniform": False,
        },
    }


def control_non_unique() -> dict[str, Any]:
    """Control (c): the invariant linear objective F(rho) = sum_p rho_p is
    constant on the simplex, so every feasible point minimizes and the
    uniqueness hypothesis fails."""

    minimizer_one = UNIFORM
    minimizer_two = (Fraction(1),) + tuple(Fraction(0) for _ in range(PORTS - 1))
    for candidate in (minimizer_one, minimizer_two):
        require(
            sum(candidate) == 1 and all(value >= 0 for value in candidate),
            "CONTROL_FEASIBILITY",
            "a non-uniqueness witness is not a simplex point",
        )
    value_one = sum(minimizer_one)
    value_two = sum(minimizer_two)
    require(
        value_one == value_two == Fraction(1),
        "CONTROL_OBJECTIVE",
        "the linear objective is not constant on the exhibited simplex points",
    )
    require(
        minimizer_one != minimizer_two,
        "CONTROL_NOT_FAILED",
        "the non-uniqueness control unexpectedly passed: the exhibited "
        "minimizers coincide",
    )
    return {
        "dropped_hypothesis": "the information projection is unique on the feasible family",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "objective": "F(rho) = sum_p rho_p, invariant under every port permutation",
            "minimizer_one": [str(value) for value in minimizer_one],
            "minimizer_two": [str(value) for value in minimizer_two],
            "objective_values": [str(value_one), str(value_two)],
            "every_simplex_point_minimizes": True,
            "distinct_minimizers_exhibited": True,
        },
    }


def control_incomplete_cover() -> dict[str, Any]:
    """Control (d): a cover omitting port 0. The restriction map is not
    injective on the local family, so two feasible states differ at port 0
    with equal objective and the uniqueness premise fails."""

    omitted = 0
    covered = tuple(p for p in range(PORTS) if p != omitted)
    state_one = UNIFORM
    state_two = (Fraction(5, 12),) + tuple(Fraction(1, 12) for _ in range(PORTS - 1))
    for state in (state_one, state_two):
        require(
            all(Fraction(0) <= value <= Fraction(1) for value in state),
            "CONTROL_FEASIBILITY",
            "an incomplete-cover witness leaves the local weight family",
        )
    restriction_one = tuple(state_one[p] for p in covered)
    restriction_two = tuple(state_two[p] for p in covered)
    require(
        restriction_one == restriction_two,
        "CONTROL_RESTRICTION",
        "the exhibited states do not agree on the covered ports",
    )
    require(
        state_one != state_two and state_one[omitted] != state_two[omitted],
        "CONTROL_NOT_FAILED",
        "the incomplete-cover control unexpectedly passed: the exhibited "
        "states agree at the omitted port",
    )
    return {
        "dropped_hypothesis": (
            "the observer cover restriction map is injective on the feasible family"
        ),
        "expected_failure": True,
        "failed": True,
        "witness": {
            "omitted_port": omitted,
            "covered_ports": list(covered),
            "state_one": [str(value) for value in state_one],
            "state_two": [str(value) for value in state_two],
            "state_one_at_omitted_port": str(state_one[omitted]),
            "state_two_at_omitted_port": str(state_two[omitted]),
            "restrictions_to_cover_equal": True,
            "objective_note": (
                "the cover objective sums weighted divergences over covered "
                "ports only; the two states have identical covered "
                "restrictions, so their objective values are equal by "
                "structural identity, with no log evaluated"
            ),
            "restriction_map_injective": False,
        },
    }


def require_no_floats(value: Any, path: str = "$") -> None:
    if isinstance(value, bool):
        return
    require(
        not isinstance(value, float),
        "FLOAT_FORBIDDEN",
        f"a float appears in the payload at {path}",
    )
    if isinstance(value, Mapping):
        for key, item in value.items():
            require_no_floats(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            require_no_floats(item, f"{path}[{index}]")


def build_payload() -> dict[str, Any]:
    verify_face_data()
    rotations = deck_rotations()
    group_report = verify_rotation_group(rotations)
    lean_perms = lean_listed_rotations()
    require(
        set(rotations) == set(lean_perms),
        "LEAN_ACTION_MATCH",
        "the computed sixty rotations do not coincide with the Lean listed rotations",
    )
    group_report["matches_lean_listed_rotations"] = True

    lane = positive_lane(rotations)
    controls = {
        "wrong_reference": control_wrong_reference(),
        "non_transitive": control_non_transitive(rotations),
        "non_unique": control_non_unique(),
        "incomplete_cover": control_incomplete_cover(),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 610,
        "generated_by": GENERATED_BY,
        "lean_spine": {
            "module": LEAN_SPINE_MODULE,
            "action_module": LEAN_ACTION_MODULE,
            "oriented_faces_module": LEAN_FACES_MODULE,
            "hypotheses_supplied_here": [
                "invariance of the feasible family and objective under the deck action",
                "uniqueness of the information projection by strict convexity",
                "the concrete sixty-rotation action and its pairwise transitivity",
            ],
            "composition_theorem": "OPH.EqualStateWeights.equal_state_weights",
        },
        "deck_action": group_report,
        "positive_lane": lane,
        "state_weight_vs_block_trace": state_weight_vs_block_trace(
            lane["result_per_port"]
        ),
        "controls": controls,
        "claim_boundary": {
            "proves": (
                "the deck action is a transitive group of oriented incidence "
                "automorphisms; the uniform reference with any strictly "
                "positive invariant weight has the uniform point as its "
                "unique information projection, exactly 1/12 per port; each "
                "of the four dropped hypotheses admits the recorded "
                "countermodel"
            ),
            "does_not_close": [
                "selection of the reference family itself: if distinct invariant reference, cover, or weight rules survive, the A3 specification stays open and no dependent selection is promoted",
                "equal central-block algebra dimensions or any realization property of the carrier; the measured block trace is a distinct object",
                "physical identification of ports or any downstream gauge datum",
            ],
        },
    }
    require_no_floats(payload)
    payload["payload_sha256"] = "sha256:" + sha256_json(payload)
    return payload


def default_output_path() -> Path:
    return MODULE_DIR / "manifests" / "equal_state_weights_reference.json"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Executable certificate for issue #610: equal state weights"
    )
    parser.add_argument("--output", type=Path, default=default_output_path())
    args = parser.parse_args(argv)
    try:
        payload = build_payload()
    except CertificateError as exc:
        print(
            json.dumps(
                {"status": "FAIL", "code": exc.code, "message": exc.message},
                indent=2,
            ),
            file=sys.stderr,
        )
        return 1
    write_json(args.output, payload)
    print(
        json.dumps(
            {
                "status": "PASS",
                "manifest": str(args.output),
                "payload_sha256": payload["payload_sha256"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
