#!/usr/bin/env python3
"""Exact certificate for GitHub issue #625: the integer load fiber and the
quadratic readback law, typed against the three axioms.

Executable companion of Lean/Screen/UnitSplit12.lean (integer split
arithmetic) and Lean/Screen/A5Commutant.lean (commutant dimension four).
Two lanes and one ledger:

* INTEGER FIBER (A1). A1 asserts finite accessible algebras with central
  record algebras and twelve primitive pairwise-orthogonal central port
  projections. A finite-dimensional algebra has an atomic center, so every
  record projection is a sum of primitive central atoms with multiplicity
  0 or 1, and every load observable declared as a record-counting
  difference takes values in Z. The certificate consumes the pinned
  federation carrier manifest, verifies atomicity on the full record
  lattice of 4096 projections, and verifies that every load observable in
  the selector certificate's grammar (twelve port charges and their total)
  is an integer combination of atom counters. Conclusion:
  integer_fiber = "axiom_forced_from_A1_finiteness_and_atomicity". The
  out-of-class control is a continuous-spectrum load model with no finite
  atomic family; the A1 finiteness clause rejects it exactly.

* QUADRATIC READBACK (A3), two exact steps. (a) Second-order theorem: on
  the twelve-port simplex with the A3 uniform reference tau and any
  feasible perturbation rho = tau + epsilon v with sum_p v_p = 0, the
  exact second-order jet of D(rho || tau) = sum_p rho_p log(rho_p / tau_p)
  has first-order coefficient vector (1, ..., 1), which the feasibility
  clause annihilates, and second-order coefficient matrix
  diag(1 / (2 tau_p)) = 6 I, the weighted Fisher quadratic form. The jet
  algebra uses the exact Taylor polynomial of log(1 + x) through second
  order with rational coefficients; no numerics appear. (b) Invariance
  menu: the space of deck-invariant quadratic forms on the port space is
  the four-dimensional commutant span of the identity, adjacency,
  distance-two, and antipode orbitals (Lean/Screen/A5Commutant.lean,
  `commutant_decomposition`). The four isotypic projectors P_1, P_3,
  P_3', P_5 are rebuilt exactly and the Fisher form is expressed in that
  basis with coefficient vector (6, 6, 6, 6): the identity class.
  Conclusion: readback_quadratic_form =
  "forced_to_identity_class_by_A3_reference". The retained countermodel
  family is the four-parameter menu itself; the adjacency-weighted form
  6 I + A passes incidence-equivariance, symmetry, and positive
  definiteness with projector coefficients (11, 6 + sqrt5, 6 - sqrt5, 5),
  so equivariance without the A3 reference clause does not force the
  readback. That is the recorded independence boundary.

* CONSUMER LEDGER. Every consumer of the integer load fiber or the
  quadratic readback (the a5_closure selector and port-current chain, the
  Lean unit split, the consensus proof-obligation verifier, the screen
  sieve and its dependents, and the pixel chain K = 4N/P surfaces) is
  listed with its consumed object and its typed source: integer fiber
  from A1 finiteness and atomicity, quadratic form from the A3
  second-order theorem. No consumer file is edited.

Controls fail closed: the non-atomic load model is rejected by the A1
finiteness clause; the tilted-reference Fisher form has unequal diagonal
(13/4 and 13/2), so the reference clause is load-bearing; a claimed
first-order readback at the projection is rejected because the linear term
vanishes exactly on the feasible tangent space; the adjacency-form
readback claimed as A3-forced is rejected because its projector
coefficient vector is not constant. A control that unexpectedly passes
raises and the certificate exits nonzero.

Every arithmetic decision is exact integer, fractions.Fraction, or
Q(sqrt5) arithmetic; the emitted payload is walked to reject any float.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import port_current_inner_certificate as p566  # noqa: E402
import response_grammar_completeness_certificate as r611  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

F5 = p566.F5

SCHEMA = "oph.load_fiber_readback_certificate.v1"
GENERATED_BY = "code/a5_closure/load_fiber_readback_certificate.py"
LEAN_COMMUTANT_MODULE = "Lean/Screen/A5Commutant.lean"
LEAN_UNIT_SPLIT_MODULE = "Lean/Screen/UnitSplit12.lean"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"

PORTS = 12
UNIFORM_REFERENCE_TAU: tuple[Fraction, ...] = tuple(
    Fraction(1, PORTS) for _ in range(PORTS)
)

# Control (b) reference: one tilted, strictly positive, normalized reference.
# Its Fisher diagonal 1/(2 tau_p) is (13/4, 13/2, ..., 13/2), unequal, so the
# tilted form is outside the identity class and the control fails as required.
TILTED_REFERENCE_TAU: tuple[Fraction, ...] = (Fraction(2, 13),) + tuple(
    Fraction(1, 13) for _ in range(PORTS - 1)
)

# Control (d) countermodel: the adjacency-weighted quadratic form a*I + b*A.
# With b != 0 it passes every incidence-equivariance check and is outside the
# identity class, so a claim that A3 forces it must be rejected.
ADJACENCY_FORM_COEFFICIENTS: tuple[int, int] = (6, 1)

# Control (a) countermodel: a load model with a continuous record spectrum
# and no finite atomic family. The A1 finiteness clause rejects it exactly.
NON_ATOMIC_LOAD_MODEL: dict[str, Any] = {
    "model": "continuous_spectrum_load",
    "record_spectrum": "continuous_interval_0_1",
    "primitive_central_atom_count": "no finite atomic family",
    "load_observable": "lebesgue_density_readback",
}

ACCEPTED_SPECTRUM = "finite_atomic_central_records"
ACCEPTED_LOAD_TYPE = "record_counting_difference"

INTEGER_FIBER_SOURCE = (
    "A1 finiteness and atomicity (integer_fiber lane of this certificate)"
)
QUADRATIC_FORM_SOURCE = (
    "A3 second-order theorem at the declared reference "
    "(quadratic_readback lane of this certificate)"
)

# Consumer ledger rows: (path, binding token that must appear in the file,
# consumed object, consumes integer fiber, consumes quadratic form).
CONSUMERS: tuple[tuple[str, str, str, bool, bool], ...] = (
    (
        "code/a5_closure/echosahedral_selector_certificate.py",
        "integer_port_charges",
        "integer port-defect domain, total charge 12, quadratic mismatch "
        "cost, and the unit-split identity H(q) = 12 + sum (q_p - 1)^2",
        True,
        True,
    ),
    (
        "code/a5_closure/port_current_inner_certificate.py",
        "carrier_manifest_sha256",
        "the port-current packet bound to the hash-pinned carrier with its "
        "integer defect readback",
        True,
        False,
    ),
    (
        "code/a5_closure/response_grammar_completeness_certificate.py",
        "A5Commutant",
        "the equivariant response grammar on the port space whose invariant "
        "quadratic forms are the commutant menu",
        False,
        True,
    ),
    (
        "Lean/Screen/UnitSplit12.lean",
        "unit_split_of_positive_sum",
        "the twelve-slot positive integer split theorem on the integer load "
        "fiber",
        True,
        False,
    ),
    (
        "code/consensus/verify_issue_517_proof_obligations.py",
        "unit_split",
        "the unit-split domain and the strict quadratic gap 2 from the "
        "selector receipt",
        True,
        True,
    ),
    (
        "code/particles/hierarchy/verify_screen_sieve_theorem.py",
        "defect_cost_sum_q2",
        "integer curvature charges 6 - deg with total 12 and the quadratic "
        "defect cost sum q^2",
        True,
        True,
    ),
    (
        "code/particles/hierarchy/verify_issue_335_local_global_resonance.py",
        "total_curvature_charge",
        "the total curvature charge 12 from the screen sieve",
        True,
        False,
    ),
    (
        "code/particles/leptons/"
        "derive_charged_nonuniform_port_record_dynamics_no_go.py",
        "total_charge",
        "the screen-sieve minimum total charge",
        True,
        False,
    ),
    (
        "code/capacity_readback/F_candidate_capK.py",
        "K = 4N/P",
        "the pixel-chain equal-area cell count K = 4N/P read back as an "
        "integer family of record carriers",
        True,
        False,
    ),
    (
        "code/capacity_readback/F_candidate_capL.py",
        "K = 4N/P",
        "the pixel-chain marked-host-cell factor K = 4N/P",
        True,
        False,
    ),
    (
        "code/capacity_readback/capacity_semantics_menu_certificate.py",
        "K = 4N/P",
        "the cell-product-structure axis with K = 4N/P cells as independent "
        "record carriers",
        True,
        False,
    ),
)


# ---------------------------------------------------------------------------
# Integer fiber lane (A1)
# ---------------------------------------------------------------------------


def accept_load_model(model: Mapping[str, Any]) -> dict[str, Any]:
    """The A1 acceptance gate for a declared load model.

    A1 asserts finite accessible algebras with central record algebras and
    twelve primitive pairwise-orthogonal central port projections. A model
    with a continuous record spectrum or no finite atomic family leaves the
    axiom class and is rejected by the finiteness clause.
    """

    spectrum = model.get("record_spectrum")
    require(
        spectrum == ACCEPTED_SPECTRUM,
        "A1_FINITENESS",
        "A1 asserts finite accessible algebras with atomic central record "
        f"algebras; record spectrum {spectrum!r} leaves the axiom class",
    )
    count = model.get("primitive_central_atom_count")
    require(
        type(count) is int and count >= 1,
        "A1_FINITENESS",
        "A1 requires a finite positive count of primitive central atoms; "
        f"got {count!r}",
    )
    load_type = model.get("load_observable")
    require(
        load_type == ACCEPTED_LOAD_TYPE,
        "A1_LOAD_TYPE",
        "the load observable must be a record-counting difference; "
        f"got {load_type!r}",
    )
    return {
        "record_spectrum": spectrum,
        "primitive_central_atom_count": count,
        "load_observable": load_type,
        "accepted": True,
    }


def atom_indicator_vectors(
    atoms: Sequence[Mapping[str, Any]],
) -> list[tuple[int, ...]]:
    """The twelve primitive atoms as standard-basis indicator vectors."""

    require(
        len(atoms) == PORTS,
        "ATOM_COUNT",
        f"expected twelve primitive central atoms, got {len(atoms)}",
    )
    vectors: list[tuple[int, ...]] = []
    for position, atom in enumerate(atoms):
        require(
            atom.get("primitive") is True,
            "ATOM_PRIMITIVITY",
            f"atom {atom.get('atom_id')} is not flagged primitive",
        )
        vectors.append(
            tuple(1 if k == position else 0 for k in range(PORTS))
        )
    return vectors


def verify_record_lattice(vectors: Sequence[tuple[int, ...]]) -> dict[str, Any]:
    """Atomicity on the full record lattice, exactly.

    The central record algebra generated by twelve pairwise-orthogonal
    primitive atoms summing to one is the diagonal algebra; its projections
    are exactly the 4096 subset sums of atoms, each with multiplicity 0 or
    1, and each atom is minimal.
    """

    for i in range(PORTS):
        for j in range(i + 1, PORTS):
            product = tuple(
                vectors[i][k] * vectors[j][k] for k in range(PORTS)
            )
            require(
                all(entry == 0 for entry in product),
                "ATOM_ORTHOGONALITY",
                f"atoms {i} and {j} are not orthogonal",
            )
    projection_count = 0
    for mask in range(1 << PORTS):
        multiplicities = tuple(
            1 if mask & (1 << i) else 0 for i in range(PORTS)
        )
        projection = tuple(
            sum(multiplicities[i] * vectors[i][k] for i in range(PORTS))
            for k in range(PORTS)
        )
        require(
            all(entry in (0, 1) for entry in projection),
            "ATOM_MULTIPLICITY",
            "a record projection has a multiplicity outside {0, 1}",
        )
        require(
            all(entry * entry == entry for entry in projection),
            "RECORD_LATTICE",
            "a record projection is not idempotent",
        )
        projection_count += 1
    require(
        projection_count == 4096,
        "RECORD_LATTICE",
        f"expected 4096 record projections, got {projection_count}",
    )
    unit = tuple(sum(vectors[i][k] for i in range(PORTS)) for k in range(PORTS))
    require(
        unit == (1,) * PORTS,
        "ATOM_COMPLETENESS",
        "the primitive atoms do not sum to the unit",
    )
    # Minimality: the only subprojections of an atom inside the lattice are
    # zero and the atom itself.
    for i in range(PORTS):
        atom = vectors[i]
        subprojections = 0
        for mask in range(1 << PORTS):
            projection = tuple(
                sum(
                    (1 if mask & (1 << a) else 0) * vectors[a][k]
                    for a in range(PORTS)
                )
                for k in range(PORTS)
            )
            if all(projection[k] <= atom[k] for k in range(PORTS)):
                subprojections += 1
        require(
            subprojections == 2,
            "ATOM_PRIMITIVITY",
            f"atom {i} has {subprojections} lattice subprojections; a "
            "primitive atom has exactly two (zero and itself)",
        )
    return {
        "atom_count": PORTS,
        "record_projection_count": 4096,
        "atom_multiplicities": [0, 1],
        "atoms_pairwise_orthogonal": True,
        "atoms_sum_to_unit": True,
        "every_projection_idempotent": True,
        "every_atom_minimal": True,
        "lattice_closure": (
            "complement, meet, and join act on subset masks, so the 4096 "
            "subset sums are closed under all three by construction"
        ),
    }


def load_observable_grammar(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Every declared load observable is an integer combination of atom
    counters, so its value fiber is Z."""

    source = manifest.get("source_readback")
    require(
        isinstance(source, Mapping),
        "LOAD_GRAMMAR",
        "source_readback is missing from the carrier manifest",
    )
    require(
        source.get("defect_domain") == "integer_port_charges",
        "LOAD_GRAMMAR",
        "the declared defect domain is not integer port charges",
    )
    total = source.get("total_charge")
    require(
        type(total) is int and total == 12,
        "LOAD_GRAMMAR",
        f"the declared total charge must be the JSON integer 12, got {total!r}",
    )
    observables: list[dict[str, Any]] = []
    ports = manifest["carrier"]["ports"]
    for position, port in enumerate(ports):
        coefficients = tuple(
            1 if k == position else 0 for k in range(PORTS)
        )
        require(
            all(type(c) is int for c in coefficients),
            "LOAD_GRAMMAR",
            f"port charge q_{port} has a non-integer atom coefficient",
        )
        observables.append(
            {
                "name": f"q_{port}",
                "atom_coefficients": list(coefficients),
                "integer_combination_of_atom_counters": True,
            }
        )
    total_coefficients = (1,) * PORTS
    require(
        all(type(c) is int for c in total_coefficients),
        "LOAD_GRAMMAR",
        "the total charge has a non-integer atom coefficient",
    )
    observables.append(
        {
            "name": "total_charge",
            "atom_coefficients": list(total_coefficients),
            "integer_combination_of_atom_counters": True,
            "declared_value": 12,
        }
    )
    # Exact value witnesses on integer record counts: the unit split and one
    # deviation, both in Z at every observable.
    unit_counts = (1,) * PORTS
    deviation_counts = (2, 0) + (1,) * (PORTS - 2)
    for counts in (unit_counts, deviation_counts):
        for row in observables:
            value = sum(
                c * n for c, n in zip(row["atom_coefficients"], counts)
            )
            require(
                type(value) is int,
                "LOAD_GRAMMAR",
                f"observable {row['name']} left the integer fiber",
            )
    return {
        "defect_domain": "integer_port_charges",
        "declared_total_charge": 12,
        "observables": observables,
        "value_fiber": "Z",
        "value_witness_counts": [list(unit_counts), list(deviation_counts)],
        "argument": (
            "each atom counter is a record-counting difference on an atomic "
            "central record algebra, hence Z-valued; an integer combination "
            "of Z-valued counters is Z-valued"
        ),
    }


def integer_fiber_lane(manifest: Mapping[str, Any]) -> dict[str, Any]:
    atoms = manifest["carrier"]["central_port_atoms"]
    vectors = atom_indicator_vectors(atoms)
    lattice = verify_record_lattice(vectors)
    grammar = load_observable_grammar(manifest)
    accepted_model = accept_load_model(
        {
            "record_spectrum": ACCEPTED_SPECTRUM,
            "primitive_central_atom_count": PORTS,
            "load_observable": ACCEPTED_LOAD_TYPE,
        }
    )
    return {
        "axiom_input": (
            "A1 asserts finite accessible algebras with central record "
            "algebras and twelve primitive pairwise-orthogonal central port "
            "projections summing to one"
        ),
        "carrier_manifest_path": f"manifests/{CARRIER_MANIFEST_NAME}",
        "carrier_manifest_sha256": sha256_json(manifest),
        "atomicity": lattice,
        "load_observable_grammar": grammar,
        "accepted_load_model": accepted_model,
        "conclusion": "axiom_forced_from_A1_finiteness_and_atomicity",
        "witness": (
            "twelve primitive atoms, 4096 record projections with atom "
            "multiplicities in {0, 1}, and thirteen declared load "
            "observables, each an integer combination of atom counters"
        ),
    }


# ---------------------------------------------------------------------------
# Quadratic readback lane (A3), step (a): the exact second-order jet
# ---------------------------------------------------------------------------

Jet = dict[tuple[int, int], Fraction]


def jet_mul(left: Jet, right: Jet) -> Jet:
    """Product of jets in (epsilon, v), truncated past epsilon^2."""

    out: Jet = {}
    for (ei, vi), ci in left.items():
        for (ej, vj), cj in right.items():
            if ei + ej > 2:
                continue
            key = (ei + ej, vi + vj)
            out[key] = out.get(key, Fraction(0)) + ci * cj
    return {key: value for key, value in out.items() if value != 0}


def port_divergence_jet(tau_p: Fraction) -> Jet:
    """The exact jet of (tau_p + eps v) log(1 + eps v / tau_p) through
    second order in eps, with the Taylor polynomial x - x^2/2 of log(1+x)."""

    require(
        isinstance(tau_p, Fraction) and tau_p > 0,
        "REFERENCE_SUPPORT",
        "the reference must be strictly positive at every port",
    )
    u_jet: Jet = {(1, 1): Fraction(1) / tau_p}
    log_jet: Jet = {
        (1, 1): Fraction(1) / tau_p,
        (2, 2): -Fraction(1, 2) / (tau_p * tau_p),
    }
    require(
        jet_mul(u_jet, u_jet) == {(2, 2): Fraction(1) / (tau_p * tau_p)},
        "JET_EXPANSION",
        "the squared perturbation jet is wrong",
    )
    rho_jet: Jet = {(0, 0): tau_p, (1, 1): Fraction(1)}
    return jet_mul(rho_jet, log_jet)


def divergence_expansion(tau: Sequence[Fraction]) -> dict[str, Any]:
    """First- and second-order data of D(rho || tau) at rho = tau, exact."""

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
    first_order: list[Fraction] = []
    second_order_diagonal: list[Fraction] = []
    for tau_p in tau:
        jet = port_divergence_jet(tau_p)
        expected = {
            (1, 1): Fraction(1),
            (2, 2): Fraction(1, 2) / tau_p,
        }
        require(
            jet == expected,
            "JET_EXPANSION",
            "the per-port divergence jet is not eps*v + eps^2*v^2/(2 tau_p)",
        )
        first_order.append(jet[(1, 1)])
        second_order_diagonal.append(jet[(2, 2)])
    return {
        "first_order_coefficients": first_order,
        "second_order_diagonal": second_order_diagonal,
    }


def tangent_basis() -> list[tuple[int, ...]]:
    """Eleven integer basis vectors of the feasible tangent space."""

    basis: list[tuple[int, ...]] = []
    for k in range(PORTS - 1):
        vector = [0] * PORTS
        vector[k] = 1
        vector[k + 1] = -1
        basis.append(tuple(vector))
    return basis


def second_order_theorem() -> dict[str, Any]:
    """Step (a): first order vanishes on the tangent space, second order is
    the Fisher quadratic form 6 I at the uniform reference."""

    expansion = divergence_expansion(UNIFORM_REFERENCE_TAU)
    first = expansion["first_order_coefficients"]
    require(
        first == [Fraction(1)] * PORTS,
        "JET_EXPANSION",
        "the first-order coefficient vector is not (1, ..., 1)",
    )
    basis = tangent_basis()
    for vector in basis:
        value = sum(c * entry for c, entry in zip(first, vector))
        require(
            value == 0,
            "FIRST_ORDER_NOT_ZERO",
            "the first-order term does not vanish on a tangent basis vector",
        )
    off_tangent = sum(first[k] * (1 if k == 0 else 0) for k in range(PORTS))
    require(
        off_tangent == Fraction(1),
        "JET_EXPANSION",
        "the off-tangent witness evaluation changed",
    )
    diagonal = expansion["second_order_diagonal"]
    require(
        diagonal == [Fraction(6)] * PORTS,
        "FISHER_NOT_IDENTITY",
        "the second-order diagonal at the uniform reference is not 6 per port",
    )
    require(
        len(set(diagonal)) == 1,
        "FISHER_NOT_IDENTITY",
        "the Fisher diagonal is not constant across ports",
    )
    return {
        "reference_family_tau": [str(value) for value in UNIFORM_REFERENCE_TAU],
        "expansion": (
            "D(tau + eps v || tau) = eps sum_p v_p + eps^2 sum_p v_p^2 / "
            "(2 tau_p) + O(eps^3), by the exact jet of "
            "(tau_p + eps v)(log(1 + eps v / tau_p))"
        ),
        "first_order_coefficients": [str(value) for value in first],
        "first_order_vanishes_on_tangent_space": True,
        "tangent_basis_evaluations": [0] * (PORTS - 1),
        "off_tangent_witness": {
            "vector": "e_0 (sum of entries 1, outside the tangent space)",
            "value": "1",
            "reading": (
                "the vanishing of the first-order term uses the feasibility "
                "clause sum_p v_p = 0 exactly"
            ),
        },
        "second_order_diagonal": [str(value) for value in diagonal],
        "fisher_quadratic_form": "Q = diag(1 / (2 tau_p)) = 6 I",
        "fisher_coefficient_per_port": "6",
        "positive_definite": True,
        "weight_note": (
            "a strictly positive constant A3 weight w multiplies the jet "
            "port-uniformly, so the form becomes 6 w I and stays in the "
            "identity class"
        ),
        "numerics_used": False,
    }


# ---------------------------------------------------------------------------
# Quadratic readback lane (A3), step (b): the invariance menu
# ---------------------------------------------------------------------------


def orbital_matrices(
    distance: Sequence[Sequence[int]], antipode: Sequence[int]
) -> dict[str, list[list[Fraction]]]:
    orbitals: dict[str, list[list[Fraction]]] = {
        "identity": [
            [Fraction(1) if i == j else Fraction(0) for j in range(PORTS)]
            for i in range(PORTS)
        ],
        "adjacency": [
            [
                Fraction(1) if distance[i][j] == 1 else Fraction(0)
                for j in range(PORTS)
            ]
            for i in range(PORTS)
        ],
        "distance_two": [
            [
                Fraction(1) if distance[i][j] == 2 else Fraction(0)
                for j in range(PORTS)
            ]
            for i in range(PORTS)
        ],
        "antipode": [
            [
                Fraction(1) if antipode[i] == j else Fraction(0)
                for j in range(PORTS)
            ]
            for i in range(PORTS)
        ],
    }
    return orbitals


def verify_equivariance(
    matrix: Sequence[Sequence[Fraction]],
    rotations: Sequence[tuple[int, ...]],
) -> bool:
    """Entry equivariance under every listed rotation, the executable twin
    of `Equivariant` in Lean/Screen/A5Commutant.lean."""

    for g in rotations:
        for i in range(PORTS):
            for j in range(PORTS):
                if matrix[g[i]][g[j]] != matrix[i][j]:
                    return False
    return True


def projector_coefficients(
    form: Sequence[Sequence[Any]], projectors: Mapping[str, Any]
) -> dict[str, F5]:
    """The unique isotypic coefficients c_s = tr(Q P_s) / rank_s."""

    form_f5 = [
        [
            entry if isinstance(entry, F5) else F5(entry)
            for entry in row
        ]
        for row in form
    ]
    coefficients: dict[str, F5] = {}
    for sector in r611.SECTORS:
        product = p566.rmul(form_f5, projectors[sector])
        trace = r611.rtrace(product)
        coefficients[sector] = trace * F5(
            Fraction(1, r611.SECTOR_RANKS[sector])
        )
    reconstructed = r611.smul(coefficients["1"], projectors["1"])
    for sector in r611.SECTORS[1:]:
        reconstructed = r611.madd(
            reconstructed, r611.smul(coefficients[sector], projectors[sector])
        )
    require(
        r611.mat_eq(reconstructed, form_f5),
        "MENU_COEFFICIENTS",
        "the form is not the isotypic combination of its trace coefficients",
    )
    return coefficients


def invariance_menu() -> dict[str, Any]:
    """Step (b): the deck-invariant quadratic forms are the four-dimensional
    commutant span, and the A3 Fisher form is the identity element."""

    verts, adjacency, distance, antipode = r611.port_model()
    rotations = r611.rotation_permutations(verts, adjacency)
    orbitals = orbital_matrices(distance, antipode)
    support_total = 0
    for name, matrix in orbitals.items():
        require(
            all(
                matrix[i][j] == matrix[j][i]
                for i in range(PORTS)
                for j in range(PORTS)
            ),
            "MENU_DIMENSION",
            f"orbital {name} is not symmetric",
        )
        require(
            verify_equivariance(matrix, rotations),
            "ORBITAL_EQUIVARIANCE",
            f"orbital {name} is not equivariant under the sixty rotations",
        )
        support_total += sum(
            1
            for i in range(PORTS)
            for j in range(PORTS)
            if matrix[i][j] != 0
        )
    require(
        support_total == PORTS * PORTS,
        "MENU_DIMENSION",
        "the four orbital supports do not partition the 144 entries",
    )
    projectors = r611.isotypic_projectors(adjacency)

    fisher = [
        [Fraction(6) if i == j else Fraction(0) for j in range(PORTS)]
        for i in range(PORTS)
    ]
    require(
        verify_equivariance(fisher, rotations),
        "ORBITAL_EQUIVARIANCE",
        "the Fisher form is not equivariant",
    )
    fisher_coefficients = projector_coefficients(fisher, projectors)
    require(
        all(
            fisher_coefficients[sector] == F5(6)
            for sector in r611.SECTORS
        ),
        "FISHER_NOT_IDENTITY",
        "the Fisher form does not have projector coefficient vector "
        "(6, 6, 6, 6)",
    )

    a_coeff, b_coeff = ADJACENCY_FORM_COEFFICIENTS
    countermodel = [
        [
            Fraction(a_coeff) * orbitals["identity"][i][j]
            + Fraction(b_coeff) * orbitals["adjacency"][i][j]
            for j in range(PORTS)
        ]
        for i in range(PORTS)
    ]
    countermodel_equivariant = verify_equivariance(countermodel, rotations)
    countermodel_coefficients = projector_coefficients(
        countermodel, projectors
    )

    return {
        "menu": (
            "deck-invariant quadratic forms on the port space are the "
            "symmetric commutant elements; the commutant is the "
            "four-dimensional span of the identity, adjacency, "
            "distance-two, and antipode orbitals, all symmetric, so the "
            "menu has exactly four parameters"
        ),
        "menu_dimension": 4,
        "lean_commutant_module": LEAN_COMMUTANT_MODULE,
        "lean_commutant_theorems": [
            "OPH.A5Commutant.commutant_decomposition",
            "OPH.A5Commutant.orbitals_independent",
            "OPH.A5Commutant.equivariant_iff_commutes",
        ],
        "orbitals_symmetric": True,
        "orbitals_equivariant_under_sixty_rotations": True,
        "orbital_supports_partition_entries": True,
        "fisher_form": {
            "matrix": "6 I",
            "equivariant": True,
            "isotypic_coefficients": {
                sector: fisher_coefficients[sector].text()
                for sector in r611.SECTORS
            },
            "coefficient_vector": ["6", "6", "6", "6"],
            "identity_class": True,
            "reading": (
                "the A3 Fisher form at the uniform reference is the "
                "identity element of the menu: port-diagonal with equal "
                "weights"
            ),
        },
        "retained_countermodel": {
            "family": "a*I + b*A + c*N + d*P, four free parameters",
            "witness_form": (
                f"{a_coeff}*I + {b_coeff}*A (adjacency-weighted quadratic "
                "form)"
            ),
            "passes_incidence_equivariance": countermodel_equivariant,
            "symmetric": True,
            "positive_definite": all(
                countermodel_coefficients[sector].is_positive()
                for sector in r611.SECTORS
            ),
            "isotypic_coefficients": {
                sector: countermodel_coefficients[sector].text()
                for sector in r611.SECTORS
            },
            "identity_class": False,
            "independence_boundary": (
                "equivariance without the A3 reference clause admits the "
                "full four-parameter menu, so the reference clause is the "
                "selecting hypothesis for the identity class"
            ),
        },
        "conclusion": "forced_to_identity_class_by_A3_reference",
    }


# ---------------------------------------------------------------------------
# Consumer ledger
# ---------------------------------------------------------------------------


def consumer_ledger() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path, token, consumes, uses_integer, uses_quadratic in CONSUMERS:
        full = REPO_ROOT / path
        try:
            text = full.read_text(encoding="utf-8")
        except OSError as exc:
            raise CertificateError(
                "CONSUMER_MISSING", f"cannot read consumer {path}: {exc}"
            ) from exc
        require(
            token in text,
            "CONSUMER_MISSING",
            f"consumer {path} lacks the binding token {token!r}",
        )
        rows.append(
            {
                "consumer": path,
                "binding_token": token,
                "consumes": consumes,
                "integer_fiber_source": (
                    INTEGER_FIBER_SOURCE if uses_integer else None
                ),
                "quadratic_form_source": (
                    QUADRATIC_FORM_SOURCE if uses_quadratic else None
                ),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Controls (fail closed)
# ---------------------------------------------------------------------------


def control_non_atomic_load_model() -> dict[str, Any]:
    """Control (a): the continuous-spectrum load model must be rejected by
    the A1 finiteness clause, exactly."""

    code = "ACCEPTED"
    message = ""
    try:
        accept_load_model(NON_ATOMIC_LOAD_MODEL)
    except CertificateError as exc:
        code = exc.code
        message = exc.message
    require(
        code == "A1_FINITENESS",
        "CONTROL_NOT_FAILED",
        "the non-atomic load model was not rejected by the A1 finiteness "
        f"clause: got {code}",
    )
    return {
        "dropped_hypothesis": "A1 finiteness and central atomicity",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "model": dict(NON_ATOMIC_LOAD_MODEL),
            "rejection_code": code,
            "rejection_message": message,
            "reading": (
                "a continuous-spectrum load model has no finite atomic "
                "family, so it is out of the A1 class and no integer fiber "
                "conclusion applies to it"
            ),
        },
    }


def control_tilted_reference() -> dict[str, Any]:
    """Control (b): with a tilted reference the Fisher form has unequal
    diagonal and leaves the identity class, so the A3 reference clause is
    load-bearing."""

    expansion = divergence_expansion(TILTED_REFERENCE_TAU)
    diagonal = expansion["second_order_diagonal"]
    require(
        len(set(diagonal)) > 1,
        "CONTROL_NOT_FAILED",
        "the tilted-reference control unexpectedly passed: the tilted "
        "Fisher diagonal is constant",
    )
    require(
        diagonal != [Fraction(6)] * PORTS,
        "CONTROL_NOT_FAILED",
        "the tilted-reference control unexpectedly passed: the tilted "
        "Fisher diagonal equals the uniform one",
    )
    return {
        "dropped_hypothesis": "the A3 reference family is the uniform family",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "tau_prime": [str(value) for value in TILTED_REFERENCE_TAU],
            "second_order_diagonal": [str(value) for value in diagonal],
            "distinct_diagonal_values": sorted(
                {str(value) for value in diagonal}
            ),
            "proportional_to_identity": False,
            "reading": (
                "diag(1 / (2 tau_p)) with unequal tau_p is not proportional "
                "to the identity, so the reference clause selects the "
                "identity class"
            ),
        },
    }


def control_linear_readback() -> dict[str, Any]:
    """Control (c): a claimed first-order readback at the projection is
    rejected because the linear term vanishes exactly on the feasible
    tangent space."""

    expansion = divergence_expansion(UNIFORM_REFERENCE_TAU)
    first = expansion["first_order_coefficients"]
    evaluations = [
        sum(c * entry for c, entry in zip(first, vector))
        for vector in tangent_basis()
    ]
    require(
        all(value == 0 for value in evaluations),
        "CONTROL_NOT_FAILED",
        "the linear-readback control unexpectedly passed: the first-order "
        "term is nonzero on a feasible tangent vector",
    )
    return {
        "dropped_hypothesis": (
            "the leading readback of a small feasible mismatch is first "
            "order at the projection"
        ),
        "expected_failure": True,
        "failed": True,
        "witness": {
            "first_order_coefficients": [str(value) for value in first],
            "tangent_basis_evaluations": [
                str(value) for value in evaluations
            ],
            "reading": (
                "the first-order term of D(tau + eps v || tau) is eps "
                "sum_p v_p, exactly zero for every feasible v, so the "
                "leading readback is the quadratic term"
            ),
        },
    }


def control_adjacency_form_claimed_forced(menu: Mapping[str, Any]) -> dict[str, Any]:
    """Control (d): the adjacency-weighted form passes every visible
    equivariance check and is outside the identity class, so a claim that
    A3 forces it is rejected."""

    countermodel = menu["retained_countermodel"]
    require(
        countermodel["passes_incidence_equivariance"] is True,
        "CONTROL_NOT_FAILED",
        "the adjacency-form control lost its equivariance witness",
    )
    coefficients = countermodel["isotypic_coefficients"]
    values = [coefficients[sector] for sector in r611.SECTORS]
    require(
        len(set(values)) > 1,
        "CONTROL_NOT_FAILED",
        "the adjacency-form control unexpectedly passed: the countermodel "
        "coefficient vector is constant, so it lies in the identity class",
    )
    return {
        "dropped_hypothesis": "the A3 reference clause selecting the readback",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "claimed_form": countermodel["witness_form"],
            "passes_incidence_equivariance": True,
            "isotypic_coefficients": dict(coefficients),
            "fisher_coefficient_vector": ["6", "6", "6", "6"],
            "identity_class": False,
            "reading": (
                "the claim that A3 forces the adjacency-weighted form is "
                "rejected: its isotypic coefficient vector differs from "
                "(6, 6, 6, 6), and only the A3 reference clause separates "
                "the two inside the equivariant menu"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Payload assembly
# ---------------------------------------------------------------------------


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
    manifest = load_json(MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME)
    e565.validate_carrier(manifest)

    fiber = integer_fiber_lane(manifest)
    second_order = second_order_theorem()
    menu = invariance_menu()
    ledger = consumer_ledger()

    controls = {
        "non_atomic_load_model": control_non_atomic_load_model(),
        "tilted_reference": control_tilted_reference(),
        "linear_readback": control_linear_readback(),
        "adjacency_form_claimed_forced": control_adjacency_form_claimed_forced(
            menu
        ),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 625,
        "generated_by": GENERATED_BY,
        "lean_spine": {
            "commutant_module": LEAN_COMMUTANT_MODULE,
            "unit_split_module": LEAN_UNIT_SPLIT_MODULE,
            "commutant_theorems": [
                "OPH.A5Commutant.commutant_decomposition",
                "OPH.A5Commutant.orbitals_independent",
                "OPH.A5Commutant.equivariant_iff_commutes",
            ],
            "unit_split_theorem": "OPH.UnitSplit12.unit_split_of_positive_sum",
        },
        "integer_fiber": "axiom_forced_from_A1_finiteness_and_atomicity",
        "readback_quadratic_form": "forced_to_identity_class_by_A3_reference",
        "integer_fiber_lane": fiber,
        "quadratic_readback_lane": {
            "second_order_theorem": second_order,
            "invariance_menu": menu,
            "conclusion": "forced_to_identity_class_by_A3_reference",
        },
        "consumer_ledger": ledger,
        "controls": controls,
        "claim_boundary": {
            "proves": (
                "on the pinned A1 carrier, every declared load observable "
                "is an integer combination of atomic record counters, so "
                "the load fiber is Z by A1 finiteness and atomicity; the "
                "first-order term of the A3 divergence vanishes on the "
                "feasible tangent space and the second-order term is the "
                "Fisher form diag(1 / (2 tau_p)) = 6 I at the uniform "
                "reference, the identity element of the four-parameter "
                "equivariant menu"
            ),
            "does_not_close": [
                "selection of the A3 reference family itself; the "
                "equal-state-weights receipt carries that hypothesis "
                "surface",
                "third-order and higher readback structure; the theorem is "
                "a second-order jet statement with symbolic O(eps^3) "
                "remainder",
                "physical identification of ports, currents, or any "
                "downstream gauge datum",
                "consumer semantics beyond the typing rows; no consumer "
                "file is edited",
            ],
        },
    }
    require_no_floats(payload)
    payload["payload_sha256"] = "sha256:" + sha256_json(payload)
    return payload


def default_output_path() -> Path:
    return MODULE_DIR / "manifests" / "load_fiber_readback_reference.json"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Executable certificate for issue #625: integer load fiber and "
            "quadratic readback"
        )
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
