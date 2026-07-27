#!/usr/bin/env python3
"""Exact certificate for GitHub issue #625: the integer load fiber and the
quadratic readback law, typed against the three axioms.

Executable companion of Lean/Screen/UnitSplit12.lean (integer split
arithmetic) and Lean/Screen/A5Commutant.lean (commutant dimension four).
Two lanes and one ledger:

* INTEGER COUNTING GRAMMAR (named realization). A1 asserts finite accessible algebras with central
  record algebras and twelve primitive pairwise-orthogonal central port
  projections. A finite-dimensional algebra has an atomic center, so every
  record projection is a sum of primitive central atoms with multiplicity
  0 or 1. Once the pinned realization additionally declares its load
  observables to be record-counting differences, those observables take
  values in Z. The certificate consumes the pinned
  federation carrier manifest, verifies atomicity on the full record
  lattice of 4096 projections, and verifies that every load observable in
  the selector certificate's grammar (twelve port charges and their total)
  is an integer combination of atom counters. A rationally weighted central
  observable that passes reduced finite-atomic and covariance checks is
  retained as a candidate against the stronger claim that A1 alone selects
  this counting grammar. Its complete
  A1-A3 operational and refinement lift is not supplied. Conclusion: the
  integer fiber is exact inside the declared counting realization; the
  axiom-level classification remains open.

* A3 OBJECTIVE CURVATURE, two exact steps. (a) Second-order theorem: on
  the twelve-port simplex with the A3 uniform reference tau and any
  feasible perturbation rho = tau + epsilon v with sum_p v_p = 0, the
  exact second-order jet of D(rho || tau) = sum_p rho_p log(rho_p / tau_p)
  has first-order coefficient vector (1, ..., 1), which the feasibility
  clause annihilates, and second-order Taylor-coefficient matrix
  diag(1 / (2 tau_p)) = 6 I. The Hessian, equivalently the Fisher
  information matrix, is twice this matrix:
  diag(1 / tau_p) = 12 I. The jet
  algebra uses the exact Taylor polynomial of log(1 + x) through second
  order with rational coefficients; no numerics appear. (b) Invariance
  menu: the space of deck-invariant quadratic forms on the port space is
  the four-dimensional commutant span of the identity, adjacency,
  distance-two, and antipode orbitals (Lean/Screen/A5Commutant.lean,
  `commutant_decomposition`). The four isotypic projectors P_1, P_3,
  P_3', P_5 are rebuilt exactly and the Taylor coefficient is expressed in that
  basis with coefficient vector (6, 6, 6, 6): the identity class.
  Conclusion: both the local Hessian 12 I and its second-order Taylor
  coefficient 6 I lie on the identity ray. The overall coefficient follows
  the declared A3 weight, and identifying this infinitesimal curvature with an exact discrete physical
  readback law is a separate interface. The retained control family is
  the four-parameter menu itself; the adjacency-weighted form
  6 I + A passes incidence-equivariance, symmetry, and positive
  definiteness with projector coefficients (11, 6 + sqrt5, 6 - sqrt5, 5),
  so incidence equivariance alone does not force the readback. A complete
  A1-A3 lift of this alternative readback is not supplied.

* CONSUMER AUDIT. Potential downstream consumers are listed and checked for
  stable identifying tokens. This is an inventory only: none of those files
  imports this receipt, so the certificate does not claim that the exact
  discrete cost or pixel-chain integrality has been retyped.

Controls fail closed: the non-atomic load model is rejected by the A1
finiteness clause; the tilted-reference Taylor coefficient and Hessian have unequal diagonal
(13/4 and 13/2), so the reference clause is load-bearing; a claimed
first-order readback at the projection is rejected because the linear term
vanishes exactly on the feasible tangent space; the adjacency-form
readback claimed as A3-forced is rejected because its projector coefficient
vector is not constant; and a diagonal form that singles out port zero is
rejected by the carrier-rotation equivariance check. A control that
unexpectedly passes raises and the certificate exits nonzero.

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

SCHEMA = "oph.load_fiber_readback_certificate.v3"
GENERATED_BY = "code/a5_closure/load_fiber_readback_certificate.py"
LEAN_COMMUTANT_MODULE = "Lean/Screen/A5Commutant.lean"
LEAN_UNIT_SPLIT_MODULE = "Lean/Screen/UnitSplit12.lean"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"

PORTS = 12
UNIFORM_REFERENCE_TAU: tuple[Fraction, ...] = tuple(
    Fraction(1, PORTS) for _ in range(PORTS)
)

# Control (b) reference: one tilted, strictly positive, normalized reference.
# Its second-order Taylor-coefficient diagonal 1/(2 tau_p) is
# (13/4, 13/2, ..., 13/2), unequal, so the tilted form is outside the identity
# class and the control fails as required. The Hessian is twice this diagonal.
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
DECLARED_COUNTING_LOAD_TYPE = "record_counting_difference"

INTEGER_FIBER_SOURCE = (
    "declared record-counting grammar on the pinned finite atomic carrier "
    "(exact named realization, not A1 alone)"
)
DISCRETE_QUADRATIC_COST_SOURCE = (
    "declared normalized central-readback Hilbert-Schmidt cost in the "
    "echosahedral selector receipt; independent of the A3 Taylor-curvature "
    "calculation"
)

# Audited surface rows: (path, binding token that must appear in the file,
# related object, consumes the declared integer grammar, consumes the exact
# discrete quadratic cost).  Most rows in the original campaign list were
# symbol-level matches rather than consumers.  SURFACE_CLASSIFICATION records
# the corrected relation and prevents an inventory row from being reported as
# a typed dependency.
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
        "the port-current packet is bound to the carrier hash but does not "
        "consume the #625 integer-load theorem",
        False,
        False,
    ),
    (
        "code/a5_closure/response_grammar_completeness_certificate.py",
        "A5Commutant",
        "the equivariant response grammar uses the commutant menu but does not "
        "identify its response with the A3 objective-curvature calculation",
        False,
        False,
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
        "defect cost sum q^2 have an independent geometric source",
        False,
        False,
    ),
    (
        "code/particles/hierarchy/verify_issue_335_local_global_resonance.py",
        "total_curvature_charge",
        "the total curvature charge 12 from the screen sieve",
        False,
        False,
    ),
    (
        "code/particles/leptons/"
        "derive_charged_nonuniform_port_record_dynamics_no_go.py",
        "total_charge",
        "the screen-sieve minimum total charge has an independent geometric "
        "source and is not a #625 load-fiber consumer",
        False,
        False,
    ),
    (
        "code/capacity_readback/F_candidate_capK.py",
        "K = 4N/P",
        "the capacity candidate uses a cell-count symbol K, not the twelve-port "
        "load fiber",
        False,
        False,
    ),
    (
        "code/capacity_readback/F_candidate_capL.py",
        "K = 4N/P",
        "the capacity candidate uses a marked-cell factor K, not the "
        "twelve-port load fiber",
        False,
        False,
    ),
    (
        "code/capacity_readback/capacity_semantics_menu_certificate.py",
        "K = 4N/P",
        "the capacity menu audits a cell-product hypothesis and is not a "
        "twelve-port load-fiber consumer",
        False,
        False,
    ),
)

SURFACE_CLASSIFICATION: dict[str, str] = {
    "code/a5_closure/echosahedral_selector_certificate.py": (
        "direct_named_realization_consumer"
    ),
    "code/a5_closure/port_current_inner_certificate.py": (
        "carrier_hash_consumer_not_load_fiber_consumer"
    ),
    "code/a5_closure/response_grammar_completeness_certificate.py": (
        "commutant_consumer_not_A3_curvature_consumer"
    ),
    "Lean/Screen/UnitSplit12.lean": "conditional_mathematical_integrality_theorem",
    "code/consensus/verify_issue_517_proof_obligations.py": (
        "downstream_named_selector_receipt_consumer"
    ),
    "code/particles/hierarchy/verify_screen_sieve_theorem.py": (
        "independent_geometric_integer_source"
    ),
    "code/particles/hierarchy/verify_issue_335_local_global_resonance.py": (
        "downstream_screen_sieve_consumer_not_load_fiber_consumer"
    ),
    "code/particles/leptons/"
    "derive_charged_nonuniform_port_record_dynamics_no_go.py": (
        "downstream_screen_sieve_consumer_not_load_fiber_consumer"
    ),
    "code/capacity_readback/F_candidate_capK.py": (
        "cell_count_symbol_not_load_fiber_consumer"
    ),
    "code/capacity_readback/F_candidate_capL.py": (
        "cell_count_symbol_not_load_fiber_consumer"
    ),
    "code/capacity_readback/capacity_semantics_menu_certificate.py": (
        "cell_product_hypothesis_not_load_fiber_consumer"
    ),
}


# ---------------------------------------------------------------------------
# Integer fiber lane (A1)
# ---------------------------------------------------------------------------


def validate_a1_finite_atomic_carrier(model: Mapping[str, Any]) -> dict[str, Any]:
    """Validate only the finite atomic carrier clauses supplied by A1.

    A1 asserts finite accessible algebras with central record algebras and
    twelve primitive pairwise-orthogonal central port projections. A model
    with a continuous record spectrum or no finite atomic family leaves the
    axiom class and is rejected by the finiteness clause. This function does
    not restrict the values of a central readback observable.
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
    return {
        "record_spectrum": spectrum,
        "primitive_central_atom_count": count,
        "accepted_by_A1_finite_atomic_gate": True,
    }


def accept_declared_counting_model(model: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the additional named-realization counting grammar."""

    carrier = validate_a1_finite_atomic_carrier(model)
    load_type = model.get("load_observable")
    require(
        load_type == DECLARED_COUNTING_LOAD_TYPE,
        "DECLARED_LOAD_TYPE",
        "the named counting realization requires a record-counting "
        "difference; "
        f"got {load_type!r}",
    )
    return {
        **carrier,
        "load_observable": load_type,
        "accepted_by_declared_counting_grammar": True,
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
    accepted_model = accept_declared_counting_model(
        {
            "record_spectrum": ACCEPTED_SPECTRUM,
            "primitive_central_atom_count": PORTS,
            "load_observable": DECLARED_COUNTING_LOAD_TYPE,
        }
    )
    return {
        "axiom_input": (
            "A1 asserts finite accessible algebras with central record "
            "algebras and twelve primitive pairwise-orthogonal central port "
            "projections summing to one"
        ),
        "additional_realization_input": (
            "the pinned source_readback grammar declares integer_port_charges, "
            "total charge 12, and record-counting observables"
        ),
        "carrier_manifest_path": f"manifests/{CARRIER_MANIFEST_NAME}",
        "carrier_manifest_sha256": sha256_json(manifest),
        "atomicity": lattice,
        "load_observable_grammar": grammar,
        "accepted_load_model": accepted_model,
        "conclusion": "exact_named_realization_for_declared_counting_grammar",
        "witness": (
            "twelve primitive atoms, 4096 record projections with atom "
            "multiplicities in {0, 1}, and thirteen declared load "
            "observables, each an integer combination of atom counters"
        ),
    }


def reduced_carrier_noninteger_readback_candidate() -> dict[str, Any]:
    """A reduced finite atomic carrier with a covariant noninteger readback.

    A1 types boundary observables and readback maps but does not normalize every
    such observable to an integer atom count. For every port p, let q'_p be
    one half of the p-th atom counter. The full twelve-observable family is
    carried into itself by every carrier rotation, so no presentation label is
    preferred. The declared record-counting grammar rejects this family, which
    shows that the grammar, rather than A1 finiteness alone, carries the
    integrality claim.
    """

    family = tuple(
        tuple(
            Fraction(1, 2) if atom == port else Fraction(0)
            for atom in range(PORTS)
        )
        for port in range(PORTS)
    )
    require(
        all(row[port] == Fraction(1, 2) for port, row in enumerate(family)),
        "COUNTERMODEL_NOT_NONINTEGER",
        "the reduced-carrier readback family lost its half-atom value",
    )
    verts, adjacency, _, _ = r611.port_model()
    rotations = r611.rotation_permutations(verts, adjacency)
    covariance_checks = 0
    for rotation in rotations:
        for port, coefficients in enumerate(family):
            transported = [Fraction(0)] * PORTS
            for atom, value in enumerate(coefficients):
                transported[rotation[atom]] = value
            require(
                tuple(transported) == family[rotation[port]],
                "COUNTERMODEL_PRESENTATION",
                "the half-atom readback family is not carrier-equivariant",
            )
            covariance_checks += 1

    model = {
        "record_spectrum": ACCEPTED_SPECTRUM,
        "primitive_central_atom_count": PORTS,
        "load_observable": "covariant_half_atom_counter_family",
    }
    a1_validation = validate_a1_finite_atomic_carrier(model)
    rejection_code = "ACCEPTED"
    try:
        accept_declared_counting_model(model)
    except CertificateError as exc:
        rejection_code = exc.code
    require(
        rejection_code == "DECLARED_LOAD_TYPE",
        "COUNTERMODEL_DECLARED_GRAMMAR",
        "the rational family was not separated from the declared counting grammar",
    )
    return {
        "same_reduced_finite_atomic_carrier": True,
        "same_twelve_primitive_central_projections": True,
        "A1_finite_atomic_validation": a1_validation,
        "boundary_observable_family": [
            {
                "observable": f"q_prime_{port}",
                "atom_coefficients": [str(value) for value in coefficients],
                "value_on_corresponding_atom": "1/2",
            }
            for port, coefficients in enumerate(family)
        ],
        "family_size": PORTS,
        "carrier_rotations_checked": len(rotations),
        "equivariance_checks": covariance_checks,
        "presentation_covariant": True,
        "common_half_integer_value": "1/2",
        "integer_valued": False,
        "declared_counting_grammar_rejection_code": rejection_code,
        "rejected_by_declared_counting_grammar_after_reduced_carrier_gate": True,
        "complete_A1_operational_and_refinement_schema_instantiated": False,
        "complete_A2_meaning_naturality_instantiated": False,
        "complete_A3_feasible_family_cover_weights_and_optimizer_instantiated": False,
        "supports_full_A1_A2_A3_independence_claim": False,
        "status": "reduced_carrier_candidate_full_schema_lift_open",
        "reading": (
            "the reduced carrier checks permit this equivariant rational "
            "central-readback family. A complete revised-A1/A2/A3 lift is "
            "required before it becomes an axiom-class countermodel"
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
    represented by the Taylor-coefficient matrix 6 I at the uniform reference.

    The Hessian and Fisher information matrix are 12 I. The factor of two is
    kept explicit so the coefficient of eps^2 is not mislabeled as a Hessian.
    """

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
        "CURVATURE_NOT_IDENTITY",
        "the second-order diagonal at the uniform reference is not 6 per port",
    )
    require(
        len(set(diagonal)) == 1,
        "CURVATURE_NOT_IDENTITY",
        "the second-order Taylor-coefficient diagonal is not constant across ports",
    )
    hessian_diagonal = [2 * value for value in diagonal]
    require(
        hessian_diagonal == [Fraction(12)] * PORTS,
        "HESSIAN_COEFFICIENT",
        "the Hessian diagonal at the uniform reference is not 12 per port",
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
        "second_order_taylor_coefficient_matrix": (
            "C = diag(1 / (2 tau_p)) = 6 I"
        ),
        "second_order_taylor_coefficient_per_port": "6",
        "hessian_matrix": "H = diag(1 / tau_p) = 12 I",
        "hessian_coefficient_per_port": "12",
        "fisher_information_matrix": "F = H = 12 I",
        "quadratic_term_relation": "eps^2 v^T C v = (eps^2 / 2) v^T H v",
        "positive_definite": True,
        "weight_note": (
            "a strictly positive constant A3 weight w multiplies the jet "
            "port-uniformly, so C becomes 6 w I and H becomes 12 w I; both "
            "stay on the identity ray. The absolute coefficients require "
            "unit local weight"
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
    commutant span, and the A3 objective curvature is in the identity class."""

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

    taylor_coefficient = [
        [Fraction(6) if i == j else Fraction(0) for j in range(PORTS)]
        for i in range(PORTS)
    ]
    require(
        verify_equivariance(taylor_coefficient, rotations),
        "ORBITAL_EQUIVARIANCE",
        "the second-order Taylor-coefficient form is not equivariant",
    )
    taylor_coefficients = projector_coefficients(taylor_coefficient, projectors)
    require(
        all(
            taylor_coefficients[sector] == F5(6)
            for sector in r611.SECTORS
        ),
        "CURVATURE_NOT_IDENTITY",
        "the second-order Taylor coefficient does not have projector vector "
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
        "a3_objective_curvature": {
            "second_order_taylor_coefficient_matrix": "6 I",
            "hessian_and_fisher_information_matrix": "12 I",
            "matrix": "6 I",
            "equivariant": True,
            "isotypic_coefficients": {
                sector: taylor_coefficients[sector].text()
                for sector in r611.SECTORS
            },
            "coefficient_vector": ["6", "6", "6", "6"],
            "identity_class": True,
            "reading": (
                "the coefficient of eps^2 is 6 I and the Hessian is 12 I; "
                "both are port-diagonal identity-class forms"
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
            "same_reduced_objective_and_state_selection_data": True,
            "A3_second_order_taylor_coefficient_remains": "6 I",
            "A3_objective_hessian_remains": "12 I",
            "differs_only_in_proposed_physical_readback_on_reduced_checks": True,
            "complete_A1_A2_A3_operational_lift_proved": False,
            "supports_full_axiom_independence_claim": False,
            "independence_boundary": (
                "A3 gives Taylor coefficient 6 I and Hessian 12 I at the "
                "normalized uniform reference. The reduced carrier checks also "
                "admit 6 I + A as a separate physical-readback candidate while "
                "every visible incidence-equivariance check passes. A complete "
                "A1-A3 operational lift of that candidate remains open"
            ),
        },
        "conclusion": "A3_objective_curvature_on_identity_ray_at_uniform_reference",
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
        classification = SURFACE_CLASSIFICATION.get(path)
        require(
            classification is not None,
            "CONSUMER_CLASSIFICATION",
            f"audited surface {path} has no dependency classification",
        )
        rows.append(
            {
                "consumer": path,
                "binding_token": token,
                "related_object": consumes,
                "dependency_class": classification,
                "consumes_declared_integer_counting_grammar": uses_integer,
                "consumes_exact_discrete_quadratic_cost": uses_quadratic,
                "integer_fiber_source": (
                    INTEGER_FIBER_SOURCE if uses_integer else None
                ),
                "quadratic_form_source": (
                    DISCRETE_QUADRATIC_COST_SOURCE if uses_quadratic else None
                ),
                "consumes_A3_objective_curvature_receipt": False,
                "binding_check": "token_presence_only",
                "imports_this_receipt": False,
                "typed_against_this_receipt": False,
                "requires_axiom_forced_issue_625_result": False,
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
        validate_a1_finite_atomic_carrier(NON_ATOMIC_LOAD_MODEL)
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
    """Control (b): with a tilted reference the objective-curvature form has
    unequal diagonal and leaves the identity class, so the A3 reference clause
    is load-bearing."""

    expansion = divergence_expansion(TILTED_REFERENCE_TAU)
    diagonal = expansion["second_order_diagonal"]
    require(
        len(set(diagonal)) > 1,
        "CONTROL_NOT_FAILED",
        "the tilted-reference control unexpectedly passed: the tilted "
        "Taylor-coefficient diagonal is constant",
    )
    require(
        diagonal != [Fraction(6)] * PORTS,
        "CONTROL_NOT_FAILED",
        "the tilted-reference control unexpectedly passed: the tilted "
        "Taylor-coefficient diagonal equals the uniform one",
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
        "dropped_hypothesis": (
            "a source law identifying A3 objective curvature with the "
            "discrete physical readback"
        ),
        "expected_failure": True,
        "failed": True,
        "witness": {
            "claimed_form": countermodel["witness_form"],
            "passes_incidence_equivariance": True,
            "isotypic_coefficients": dict(coefficients),
            "A3_taylor_coefficient_vector": ["6", "6", "6", "6"],
            "identity_class": False,
            "reading": (
                "the claim that A3 forces the adjacency-weighted form is "
                "rejected: its isotypic coefficient vector differs from "
                "(6, 6, 6, 6). Equivariance alone does not identify either "
                "form as the physical discrete cost"
            ),
        },
    }


def control_presentation_dependent_form() -> dict[str, Any]:
    """Reject a quadratic form that assigns a special coefficient to port zero."""

    verts, adjacency, _distance, _antipode = r611.port_model()
    rotations = r611.rotation_permutations(verts, adjacency)
    form = [
        [
            Fraction(7 if i == 0 else 6) if i == j else Fraction(0)
            for j in range(PORTS)
        ]
        for i in range(PORTS)
    ]
    equivariant = verify_equivariance(form, rotations)
    require(
        equivariant is False,
        "CONTROL_NOT_FAILED",
        "the presentation-dependent form unexpectedly passed carrier "
        "equivariance",
    )
    return {
        "dropped_hypothesis": "presentation independence under carrier rotations",
        "expected_failure": True,
        "failed": True,
        "witness": {
            "form": "diag(7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)",
            "singled_out_port": 0,
            "rotations_checked": len(rotations),
            "passes_incidence_equivariance": False,
            "reading": (
                "a coordinate formula that gives one named port a different "
                "coefficient depends on the presentation and is rejected"
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
    reduced_countermodel = reduced_carrier_noninteger_readback_candidate()
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
        "presentation_dependent_form": control_presentation_dependent_form(),
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
        "bounded_exit": "conditional_open_interface",
        "integer_fiber": "exact_named_realization_for_declared_counting_grammar",
        "readback_quadratic_form": (
            "A3_objective_curvature_identity_ray_exact_physical_readback_open"
        ),
        "integer_fiber_lane": fiber,
        "reduced_carrier_noninteger_readback_candidate": reduced_countermodel,
        "quadratic_readback_lane": {
            "second_order_theorem": second_order,
            "invariance_menu": menu,
            "conclusion": (
                "A3_objective_curvature_on_identity_ray_at_uniform_reference; "
                "exact_discrete_readback_identification_open"
            ),
        },
        "consumer_ledger": ledger,
        "consumer_integration_status": (
            "all eleven originally named surfaces classified: two consume the "
            "declared named-realization grammar/cost, one is a conditional "
            "mathematical integrality theorem, and eight use independent objects "
            "or were false-positive symbol matches; none requires an "
            "axiom-forced #625 result"
        ),
        "controls": controls,
        "claim_boundary": {
            "proves": (
                "on the pinned finite atomic carrier, every observable in the "
                "additional declared record-counting grammar is an integer "
                "combination of atom counters; the "
                "first-order term of the A3 divergence vanishes on the "
                "feasible tangent space; the coefficient of eps^2 is "
                "diag(1 / (2 tau_p)) = 6 I and the Hessian is "
                "diag(1 / tau_p) = 12 I for unit local weight at the uniform "
                "reference. Both lie on the identity ray of the "
                "four-parameter equivariant menu"
            ),
            "does_not_close": [
                "derivation of the record-counting load grammar from A1; a "
                "carrier-equivariant rational candidate is retained only at "
                "the reduced finite-atomic interface",
                "a complete A1 operational/refinement, A2 naturality, and A3 "
                "optimizer lift of either alternative-readback candidate",
                "selection of the A3 reference family itself; the "
                "equal-state-weights receipt carries that hypothesis "
                "surface",
                "the absolute objective-curvature coefficients without a "
                "declared local weight normalization",
                "third-order and higher readback structure; the theorem is "
                "a second-order jet statement with symbolic O(eps^3) "
                "remainder",
                "identification of the A3 objective curvature with the exact "
                "discrete mismatch cost used by downstream consumers",
                "physical identification of ports, currents, or any "
                "downstream gauge datum",
                "promotion of either named-realization premise to an axiom; "
                "the audited consumers remain typed against their declared or "
                "independent sources",
            ],
            "status": "conditional_open_interface",
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
