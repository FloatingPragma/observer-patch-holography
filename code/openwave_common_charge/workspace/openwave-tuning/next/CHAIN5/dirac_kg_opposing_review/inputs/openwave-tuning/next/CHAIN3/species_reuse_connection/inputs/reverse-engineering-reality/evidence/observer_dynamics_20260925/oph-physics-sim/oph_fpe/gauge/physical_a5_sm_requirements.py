"""Fail-closed physical A5-to-Standard-Model requirements contract.

This module turns ``survival-proof/SIMULATOR_RECEIPT_REQUIREMENTS.md`` into a
stable simulator-side stage DAG.  It deliberately does much less than a
physics producer: it can replay a ROOT *inventory*, inventory candidate
evidence, and evaluate the dependency contract.  It cannot turn a JSON boolean
or a theorem-side checklist into physical evidence.  A hash/role inventory is
not the normative ROOT receipt: typed role semantics, an actual code/build
binding, and a pre-outcome commitment are still required.

Only verifiers registered in :data:`REGISTERED_PHYSICAL_PRODUCERS` may admit a
receipt.  The registry is initially empty.  Consequently a well-formed ROOT
inventory remains ``OPEN`` and every physical emergence lane remains false.
That is the intended fail-closed state until narrow source-derived producers
and independent replay verifiers are implemented.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from oph_fpe.evidence.production_envelope import (
    COMMON_STAGE,
    PRODUCTION_BUNDLE_REPORT_ARTIFACT_TYPE,
    PRODUCTION_BUNDLE_REPORT_SCHEMA,
    verify_production_bundle_manifest,
)


ROOT_MANIFEST_SCHEMA = "oph.physical-a5-sm.root-manifest/1.0.0"
ROOT_REPORT_SCHEMA = "oph.physical-a5-sm.root-verification/1.0.0"
ROOT_REPORT_ARTIFACT_TYPE = "OPH_PHYSICAL_A5_SM_ROOT_VERIFICATION"
REQUIREMENTS_REPORT_SCHEMA = "oph.physical-a5-sm.requirements-audit/1.0.0"
REQUIREMENTS_REPORT_ARTIFACT_TYPE = "OPH_PHYSICAL_A5_SM_REQUIREMENTS_AUDIT"
MAX_MANIFEST_BYTES = 4 * 1024 * 1024
MAX_ARTIFACT_BYTES = 256 * 1024 * 1024
MAX_JSON_ARTIFACT_BYTES = 16 * 1024 * 1024
MAX_ARTIFACTS = 4096
MAX_ANCESTRY_EDGES = 32768

_SHA256_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
_COMMIT_RE = re.compile(r"^[0-9a-f]{40}(?:[0-9a-f]{24})?$")
_IDENTIFIER_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.:+/-]{0,127}$")
class TerminalStatus(str, Enum):
    """Only terminal classes licensed by the simulator receipt runbook."""

    PASS = "PASS"
    OPEN = "OPEN"
    UNRESOLVED = "UNRESOLVED"
    FAIL = "FAIL"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class ClaimScope(str, Enum):
    """Which optional downstream claim tiers are being evaluated."""

    STRUCTURAL = "structural"
    FULL_INTERACTING = "full_interacting"
    CONTINUUM = "continuum"


EXACT_SMALL_REQUIRED_CHECKS = (
    "FULL_STATE_SPACE_ENUMERATION",
    "COMPLETE_SOURCE_REGISTRY_ENUMERATION",
    "EXACT_HAMILTONIAN_AND_GAUSS_KERNEL",
    "EXACT_RIESZ_PROJECTOR_REPLAY",
    "PRIMITIVE_RESIDUE_REPLAY",
    "SCALAR_COMPETITOR_REPLAY",
    "FULL_SCREEN_J_ALL_REPLAY",
    "REFINED_COMPLEMENT_REPLAY",
    "INJECTED_NEGATIVE_CONTROL_SUITE",
)

SCALE_CAMPAIGN_REQUIRED_CHECKS = (
    "EXACT_SMALL_ORACLE_PASS",
    "SAME_VERIFIER_IDENTITY",
    "SAME_ROOT_PROVENANCE",
    "FROZEN_THRESHOLD_POLICY",
    "FROZEN_NEGATIVE_CONTROL_SUITE",
    "FROZEN_SEED_SCHEDULE",
    "FROZEN_SCALE_GRID",
)

TARGET_SELECTION_FIELDS = (
    "standard_model_multiplicities_used_to_select_rank",
    "measured_standard_model_couplings_used_to_select_candidate",
    "desired_family_count_used_to_define_source_registry",
    "desired_spectrum_labels_used_to_define_primitive_composite_status",
    "desired_global_quotient_used_to_define_public_category",
)

REQUIRED_ROOT_ROLES = (
    "simulator_code",
    "dirty_tree_manifest",
    "complete_configuration",
    "regulator_manifest",
    "state_space_or_quotient_manifest",
    "dynamics_generator",
    "source_operator_registry",
    "boundary_and_sector_manifest",
    "gauge_representation_manifest",
    "spin_representation_manifest",
    "a5_representation_manifest",
    "refinement_representation_manifest",
    "reproducibility_manifest",
    "numerical_backend_manifest",
    "candidate_domain",
    "selection_law",
)

ROOT_RECEIPTS = (
    "ROOT_MANIFEST_SCHEMA_RECEIPT",
    "ROOT_ARTIFACT_CONTAINMENT_RECEIPT",
    "ROOT_ARTIFACT_HASH_REPLAY_RECEIPT",
    "ROOT_REQUIRED_OBJECT_CENSUS_RECEIPT",
    "ROOT_REPRODUCIBILITY_POLICY_RECEIPT",
    "ROOT_SOURCE_ANCESTRY_DAG_RECEIPT",
    "ROOT_CANDIDATE_SELECTION_COMMITMENT_RECEIPT",
    "ROOT_TARGET_SELECTION_DECLARATIONS_RECEIPT",
    "ROOT_INVENTORY_REPLAY_RECEIPT",
    "ROOT_TYPED_ROLE_SEMANTICS_RECEIPT",
    "ROOT_CODE_BUILD_BINDING_RECEIPT",
    "ROOT_PREOUTCOME_COMMITMENT_RECEIPT",
    "ROOT_IMMUTABLE_PACKET_REPLAY_RECEIPT",
)

ROOT_INVENTORY_RECEIPTS = (
    "ROOT_MANIFEST_SCHEMA_RECEIPT",
    "ROOT_ARTIFACT_CONTAINMENT_RECEIPT",
    "ROOT_ARTIFACT_HASH_REPLAY_RECEIPT",
    "ROOT_REQUIRED_OBJECT_CENSUS_RECEIPT",
    "ROOT_REPRODUCIBILITY_POLICY_RECEIPT",
    "ROOT_SOURCE_ANCESTRY_DAG_RECEIPT",
    "ROOT_CANDIDATE_SELECTION_COMMITMENT_RECEIPT",
    "ROOT_TARGET_SELECTION_DECLARATIONS_RECEIPT",
)

GEOMETRY_SOURCE_RECEIPTS = (
    "GEOMETRY_ORIENTED_CLOSED_S2_SUPPORT_RECEIPT",
    "GEOMETRY_CALIBRATED_CURVATURE_RESPONSE_RECEIPT",
    "GEOMETRY_EXACT_CURVATURE_RISK_LEDGER_RECEIPT",
    "GEOMETRY_COMPLETE_SETTLEMENT_RECEIPT",
    "GEOMETRY_ATOMIC_DEFECT_EXPOSURE_RECEIPT",
    "GEOMETRY_PAIR_FISHER_INVERSE_SQUARE_RECEIPT",
    "GEOMETRY_ORIENTED_TWELVE_PORT_A5_FRAME_RECEIPT",
    "GEOMETRY_REGULAR_EDGEWISE_COFINAL_REFINEMENT_RECEIPT",
    "GEOMETRY_16_16_12_TRAP_NEGATIVE_CONTROL_RECEIPT",
    "GEOMETRY_READBACK_A5_EQUIVARIANCE_RECEIPT",
    "GEOMETRY_565_SOURCE_DERIVED_RECEIPT",
)

GEOMETRY_HARDWARE_RECEIPTS = (
    "GEOMETRY_HARDWARE_TWELVE_ORTHOGONAL_PORT_PROJECTIONS_RECEIPT",
    "GEOMETRY_HARDWARE_FIXED_POINT_FREE_ANTIPODE_RECEIPT",
    "GEOMETRY_HARDWARE_ORIENTED_REGULAR_ICOSAHEDRAL_DIRECTIONS_RECEIPT",
    "GEOMETRY_HARDWARE_FAITHFUL_A5_READBACK_ACTION_RECEIPT",
    "GEOMETRY_HARDWARE_COFINAL_PERSISTENCE_RECEIPT",
    "GEOMETRY_565_OPERATIONAL_HARDWARE_RECEIPT",
)

CURRENT_RECEIPTS = (
    "CURRENT_RAW_TWO_SIDED_RESPONSE_HISTORY_RECEIPT",
    "CURRENT_CLOCK_CALIBRATION_AND_NOISE_RECEIPT",
    "CURRENT_LOSSLESS_RECIPROCAL_RESPONSE_RECEIPT",
    "CURRENT_GENERATOR_MAP_RANK12_RECEIPT",
    "CURRENT_ORDERED_LOOP_COMMUTATOR_RECONSTRUCTION_RECEIPT",
    "CURRENT_COMMUTATOR_CLOSED_IN_PORT_TANGENT_RECEIPT",
    "CURRENT_COMPACT_SKEW_ADJOINT_RESPONSE_RECEIPT",
    "CURRENT_COMPLETE_A5_COVARIANCE_RECEIPT",
    "CURRENT_1_5_3_3PRIME_DECOMPOSITION_RECEIPT",
    "CURRENT_A5_OVERLAP_HOLONOMY_FULLNESS_RECEIPT",
    "CURRENT_SAME_SOURCE_INTERNAL_IMPLEMENTER_RECEIPT",
    "CURRENT_NO_HIDDEN_PUBLIC_TANGENT_RECEIPT",
    "CURRENT_REFINEMENT_INTERTWINER_RECEIPT",
    "CURRENT_566_PHYSICAL_RECEIPT",
)

GLOBAL_FORM_RECEIPTS = (
    "GLOBAL_FORM_FROZEN_CARRIER_CANDIDATE_DOMAIN_RECEIPT",
    "GLOBAL_FORM_MINIMAL_PUBLIC_CARRIER_SELECTION_RECEIPT",
    "GLOBAL_FORM_ORIENTED_PRIMITIVE_VOLUME_CLOCK_RECEIPT",
    "GLOBAL_FORM_PRIMITIVE_3P_PLUS_2Q_BALANCE_RECEIPT",
    "GLOBAL_FORM_CURRENT_CLOCK_INTERTWINER_RECEIPT",
    "GLOBAL_FORM_PUBLIC_TANNAKA_COMPLETENESS_RECEIPT",
    "GLOBAL_FORM_CARRIER_TENSOR_GENERATION_RECEIPT",
    "GLOBAL_FORM_PHYSICAL_PORT_LOOP_COCYCLE_RECEIPT",
    "GLOBAL_FORM_DECK_KERNEL_EXHAUSTION_RECEIPT",
    "GLOBAL_FORM_DERIVED_Z6_LINE_LATTICE_RECEIPT",
    "GLOBAL_FORM_UV_POLARIZATION_EIGHT_CLASS_SELECTOR_RECEIPT",
    "GLOBAL_FORM_567_PHYSICAL_RECEIPT",
)

SPIN_EXCHANGE_RECEIPTS = (
    "SPIN_ORIENTED_4D_COMPLEX_AND_FRAME_BUNDLE_RECEIPT",
    "SPIN_W2_ZERO_CERTIFICATE_RECEIPT",
    "SPIN_PLUS_31_COCYCLE_AND_LIFT_RECEIPT",
    "SPIN_CENTRAL_MINUS_ONE_ACTION_RECEIPT",
    "SPIN_WEYL_MATRICES_AND_CHIRALITY_SELECTOR_RECEIPT",
    "SPIN_EXTERNAL_GRADING_IDENTIFICATION_RECEIPT",
    "SPIN_OPERATOR_CAR_RECEIPT",
    "SPIN_GRADED_LOCALITY_RECEIPT",
    "SPIN_EXCHANGE_VEC_VS_SVEC_RECEIPT",
    "SPIN_G6_CHARGE_CONJUGATION_COMPATIBILITY_RECEIPT",
    "SPIN_LOCAL_AND_GLOBAL_ANOMALY_RECEIPT",
    "SPIN_COFINAL_INTERTWINER_RECEIPT",
    "SPIN_RANK15_SOURCE_DERIVED_MODULE_PROJECTOR_RECEIPT",
    "SPIN_RANK15_NONZERO_RESIDUE_REFINEMENT_RECEIPT",
    "SPIN_EXCHANGE_314_PHYSICAL_RECEIPT",
)

SOURCE_REGISTRY_RECEIPTS = (
    "SOURCE_COMPLETE_REGISTRY_ENUMERATION_RECEIPT",
    "SOURCE_OPERATOR_HASH_AND_LOCALITY_METADATA_RECEIPT",
    "SOURCE_DERIVED_ACTION_METADATA_RECEIPT",
    "SOURCE_ANCESTRY_AND_REFINEMENT_METADATA_RECEIPT",
    "SOURCE_HIDDEN_ORTHOGONAL_CONTROL_RECEIPT",
    "SOURCE_EXACT_HAMILTONIAN_OR_TRANSFER_RECEIPT",
    "SOURCE_VACUUM_GROUND_STATE_RECEIPT",
    "SOURCE_MOMENT_IDENTIFIABILITY_RECEIPT",
    "SOURCE_PRIMITIVE_QUOTIENT_OR_FILTRATION_DESCENT_RECEIPT",
    "SOURCE_PRIMITIVE_RESIDUE_EXHAUSTION_RECEIPT",
    "SOURCE_POSITIVE_AND_CHARGED_LEDGER_SPLIT_RECEIPT",
    "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES_PHYSICAL_RECEIPT",
)

SCALAR_RECEIPTS = (
    "SCALAR_COMPLETE_GAUGE_COVARIANT_PRIMITIVE_SEARCH_RECEIPT",
    "SCALAR_MANDATORY_COMPETITOR_CONTROL_CENSUS_RECEIPT",
    "SCALAR_UNIQUE_W3_RESIDUE_CHANNEL_RECEIPT",
    "SCALAR_COMPETITOR_ZERO_RESIDUE_OR_POSITIVE_GAP_RECEIPT",
    "SCALAR_MULTIPLICITY_ONE_RECEIPT",
    "SCALAR_COMPLEMENT_COMPLETE_REFINEMENT_RECEIPT",
    "SCALAR_CHANNEL_PHYSICAL_RECEIPT",
)

FAMILY_RECEIPTS = (
    "FAMILY_FULL_REAL_SCREEN_PROJECTOR_DECOMPOSITION_RECEIPT",
    "FAMILY_FROZEN_FIRST_POSITIVE_BAND_RECEIPT",
    "FAMILY_COMPLEX_RANK3_BAND_RECEIPT",
    "FAMILY_SOURCE_DERIVED_RANK15_DOMAIN_RECEIPT",
    "FAMILY_RAW_FULL_SCREEN_MIXED_RESPONSE_RECEIPT",
    "FAMILY_J_DERIVED_FROM_J_ALL_RECEIPT",
    "FAMILY_COMPLEX_RANK45_AND_IMAGE_EQUALITY_RECEIPT",
    "FAMILY_A5_G6_SPIN_EQUIVARIANCE_RECEIPT",
    "FAMILY_EXCLUDED_BAND_NULL_RESPONSE_RECEIPT",
    "FAMILY_ALL45_CROSS_CORRELATION_COVERAGE_RECEIPT",
    "FAMILY_COFINAL_REFINEMENT_AND_COMPLEMENT_RECEIPT",
    "FAMILY_569_NEGATIVE_CONTROL_SUITE_RECEIPT",
    "FAMILY_ATTACHMENT_569_PHYSICAL_RECEIPT",
)

Q1_LOCAL_ACTION_RECEIPTS = (
    "Q1_COMPLETE_FINITE_REGULATOR_DYNAMICS_RECEIPT",
    "Q1_GAUGE_ACTION_RECEIPT",
    "Q1_SCALAR_KINETIC_AND_POTENTIAL_RECEIPT",
    "Q1_FERMION_KINETIC_ACTION_RECEIPT",
    "Q1_ALLOWED_INTERACTION_BASIS_RECEIPT",
    "Q1_GAUGE_FIXING_OR_CONSTRAINT_TREATMENT_RECEIPT",
    "Q1_INTEGRATION_OR_PHASE_SPACE_DOMAIN_RECEIPT",
    "Q1_LOCALITY_DECOMPOSITION_RECEIPT",
    "Q1_BOUNDEDNESS_AND_STABILITY_RECEIPT",
    "Q1_G6_SPIN_A5_REFINEMENT_COVARIANCE_RECEIPT",
    "Q1_IDENTICAL_ACTION_BOUND_IN_Q2_RECEIPT",
    "Q1_LOCAL_ACTION_PHYSICAL_RECEIPT",
)

COMPLETE_COUPLED_DYNAMICS_RECEIPTS = (
    "COUPLED_IDENTICAL_Q1_Q2_ACTION_BINDING_RECEIPT",
    "COUPLED_GAUGE_SCALAR_FERMION_MEASURE_CENSUS_RECEIPT",
    "COUPLED_COMPLETE_PARTITION_OR_EVOLUTION_FUNCTIONAL_RECEIPT",
    "COUPLED_INTERACTION_SECTOR_NONTRIVIALITY_RECEIPT",
    "COUPLED_FINITENESS_AND_STABILITY_RECEIPT",
    "COUPLED_GAUGE_SPIN_A5_WARD_IDENTITY_RECEIPT",
    "COUPLED_SOURCE_VERTEX_DYNAMICS_BINDING_RECEIPT",
    "COUPLED_REFINEMENT_AND_COMPLEMENT_CONTROL_RECEIPT",
    "COMPLETE_COUPLED_DYNAMICS_PHYSICAL_RECEIPT",
)

FAMILY_EXACT_SYMMETRY_COMPATIBILITY_RECEIPTS = (
    "FAMILY_EXACT_ACTION_COMMUTANT_RECEIPT",
    "FAMILY_NONDEGENERATE_MATRIX_COMPATIBILITY_RECEIPT",
    "FAMILY_EXACT_CUBIC_SLOT_COUPLING_COMPATIBILITY_RECEIPT",
    "FAMILY_EXACT_REFINEMENT_STABILITY_RECEIPT",
    "FAMILY_EXACT_SYMMETRY_COMPATIBILITY_PHYSICAL_RECEIPT",
)

FAMILY_BREAKING_OR_DESCENT_RECEIPTS = (
    "FAMILY_EXACT_ACTION_COMMUTANT_RECEIPT",
    "FAMILY_BREAKING_CANDIDATE_DOMAIN_PREOUTCOME_FREEZE_RECEIPT",
    "FAMILY_BREAKING_RAW_SOURCE_HISTORY_RECEIPT",
    "FAMILY_BREAKING_TARGET_BLINDNESS_RECEIPT",
    "FAMILY_BREAKING_ORBIT_AND_STABILIZER_RECEIPT",
    "FAMILY_BREAKING_CUBIC_SLOT_COUPLING_RECEIPT",
    "FAMILY_BREAKING_REFINEMENT_STABILITY_RECEIPT",
    "FAMILY_BREAKING_OR_DESCENT_PHYSICAL_RECEIPT",
)

VERTEX_RECEIPTS = (
    "VERTEX_BOSONIC_AND_GRASSMANN_SOURCE_REGISTRY_RECEIPT",
    "VERTEX_CONNECTED_GENERATING_FUNCTIONAL_RECEIPT",
    "VERTEX_GRASSMANN_LEGENDRE_EFFECTIVE_ACTION_RECEIPT",
    "VERTEX_FULL_1PI_AMPUTATION_RECEIPT",
    "VERTEX_REDUCIBLE_CONTACT_MIXING_SUBTRACTION_RECEIPT",
    "VERTEX_FROZEN_RENORMALIZATION_POINT_RECEIPT",
    "VERTEX_THREE_INVARIANT_TENSOR_PROJECTION_RECEIPT",
    "VERTEX_COMPLEX_3X3_FAMILY_MATRIX_INTERVAL_RECEIPT",
    "VERTEX_FULL_COVARIANCE_AND_REFINEMENT_RECEIPT",
    "VERTEX_BASIS_QUOTIENTED_OBSERVABLE_RECEIPT",
    "VERTEX_GAUSSIAN_NONLINEAR_INTERPOLATOR_NEGATIVE_CONTROL_RECEIPT",
    "VERTEX_1PI_PHYSICAL_RECEIPT",
)

Q2_H_RECEIPTS = (
    "Q2_H_FINITE_KINEMATIC_HILBERT_RECEIPT",
    "Q2_H_SELF_ADJOINT_LOCAL_HAMILTONIAN_RECEIPT",
    "Q2_H_IDENTICAL_Q1_HAMILTONIAN_BINDING_RECEIPT",
    "Q2_H_MANDATORY_UPSTREAM_OBJECT_BINDING_RECEIPT",
    "Q2_H_NONZERO_COMPLETE_GAUSS_GENERATORS_RECEIPT",
    "Q2_H_GAUSS_COMMUTATOR_AND_GLOBAL_QUOTIENT_RECEIPT",
    "Q2_H_FAITHFUL_NONTRIVIAL_GAUGE_ACTION_RECEIPT",
    "Q2_H_DERIVED_NONZERO_PHYSICAL_PROJECTOR_RECEIPT",
    "Q2_H_HAMILTONIAN_GAUSS_COMPATIBILITY_RECEIPT",
    "Q2_H_PHYSICAL_GROUND_AND_OBSERVABLE_ALGEBRA_RECEIPT",
    "Q2_H_CHIRAL_INDEX_AND_ANOMALY_RECEIPT",
    "Q2_H_MIRROR_EXTRA_SECTOR_EXHAUSTION_AND_GAP_RECEIPT",
    "Q2_H_UNIFORM_LOCALITY_RECEIPT",
    "Q2_H_COMPLEMENT_COMPLETE_REFINEMENT_RECEIPT",
    "Q2_H_PHYSICAL_RECEIPT",
)

Q2_E_RECEIPTS = (
    "Q2_E_CELL_REGULATOR_AND_ADMISSIBLE_GAUGE_DOMAIN_RECEIPT",
    "Q2_E_SMOOTH_DIRAC_OPERATOR_RECEIPT",
    "Q2_E_EXACT_GINSPARG_WILSON_RECEIPT",
    "Q2_E_GAMMA5_HERMITICITY_RECEIPT",
    "Q2_E_CONSTANT_RANK_CHIRAL_PROJECTOR_RECEIPT",
    "Q2_E_UNIFORM_LOCALITY_RECEIPT",
    "Q2_E_DETERMINANT_LINE_GLOBAL_GROUPOID_RECEIPT",
    "Q2_E_MEASURE_CURRENT_ANOMALY_CURVATURE_RECEIPT",
    "Q2_E_EXHAUSTIVE_HOLONOMY_TRIVIALIZATION_RECEIPT",
    "Q2_E_EXACT_GAUGE_INVARIANCE_AND_WARD_RECEIPT",
    "Q2_E_REFINEMENT_TRANSPORT_RECEIPT",
    "Q2_E_IDENTICAL_Q1_ACTION_BINDING_RECEIPT",
    "Q2_E_MANDATORY_UPSTREAM_OBJECT_BINDING_RECEIPT",
    "Q2_E_COMPLETE_BOSONIC_GAUGE_FERMION_MEASURE_RECEIPT",
    "Q2_E_FULL_COUPLED_PARTITION_FUNCTIONAL_RECEIPT",
    "Q2_E_COUPLED_FINITENESS_AND_STABILITY_RECEIPT",
    "Q2_E_FINITE_EUCLIDEAN_STRUCTURAL_RECEIPT",
)

POSITIVITY_REFLECTION_RECEIPTS = (
    "POSITIVITY_COMPLETE_COUPLED_ALGEBRA_BINDING_RECEIPT",
    "Q2_E_REFLECTION_POSITIVITY_RECEIPT",
    "POSITIVITY_OR_POSITIVE_TRANSFER_PHYSICAL_RECEIPT",
)

POSITIVITY_TRANSFER_RECEIPTS = (
    "POSITIVITY_COMPLETE_COUPLED_ALGEBRA_BINDING_RECEIPT",
    "Q2_E_EQUIVALENT_POSITIVE_TRANSFER_HAMILTONIAN_RECEIPT",
    "POSITIVITY_OR_POSITIVE_TRANSFER_PHYSICAL_RECEIPT",
)

REFINEMENT_RECEIPTS = (
    "REFINEMENT_COMPLETE_GENERATING_ARROW_TRANSPORT_RECEIPT",
    "REFINEMENT_ISOMETRY_OR_TYPED_CONTROL_RECEIPT",
    "REFINEMENT_COMPOSITION_AND_DIAMOND_COHERENCE_RECEIPT",
    "REFINEMENT_UNIFORM_OPERATOR_NORM_BOUND_RECEIPT",
    "REFINEMENT_SUMMABLE_DEFECT_BUDGET_RECEIPT",
    "REFINEMENT_MATCHED_RIESZ_AND_RESIDUE_RANK_RECEIPT",
    "REFINEMENT_COMPLEMENT_HIDDEN_SECTOR_CONTROL_RECEIPT",
    "REFINEMENT_COFINAL_BRANCH_PATH_INDEPENDENCE_RECEIPT",
    "REFINEMENT_HELD_OUT_OBSERVABLE_ERROR_BUDGET_RECEIPT",
    "REFINEMENT_H_DIRECT_SUM_ZERO_NEGATIVE_CONTROL_RECEIPT",
    "REFINEMENT_COMPLETENESS_PHYSICAL_RECEIPT",
)

PHYSICAL_IDENTIFICATION_COMMON_RECEIPTS = (
    "PHYSICAL_ID_RAW_APPARATUS_HISTORY_RECEIPT",
    "PHYSICAL_ID_CALIBRATED_UNITS_AND_INTERVALS_RECEIPT",
    "PHYSICAL_ID_INFORMATIONAL_COMPLETENESS_RECEIPT",
    "PHYSICAL_ID_HIDDEN_SECTOR_NEGATIVE_CONTROL_RECEIPT",
    "PHYSICAL_ID_ALTERNATIVE_GEOMETRY_NEGATIVE_CONTROL_RECEIPT",
    "PHYSICAL_ID_HELD_OUT_PREDICTION_RECEIPT",
    "PHYSICAL_ID_RECORD_TO_SIMULATOR_HASH_ANCESTRY_RECEIPT",
)

PHYSICAL_IDENTIFICATION_FINITE_CLASS_RECEIPTS = (
    *PHYSICAL_IDENTIFICATION_COMMON_RECEIPTS,
    "PHYSICAL_ID_FROZEN_FINITE_ALTERNATIVE_CLASS_RECEIPT",
    "PHYSICAL_IDENTIFICATION_RECEIPT",
)

PHYSICAL_IDENTIFICATION_REDUCTION_THEOREM_RECEIPTS = (
    *PHYSICAL_IDENTIFICATION_COMMON_RECEIPTS,
    "PHYSICAL_ID_INFINITE_CLASS_FINITE_TEST_REDUCTION_THEOREM_RECEIPT",
    "PHYSICAL_IDENTIFICATION_RECEIPT",
)

Q4_OS_RECEIPTS = (
    "Q4_OS_COMPLETE_SCHWINGER_N_POINT_FAMILY_RECEIPT",
    "Q4_OS_DISTRIBUTIONAL_CONVERGENCE_RECEIPT",
    "Q4_OS_EUCLIDEAN_COVARIANCE_RECEIPT",
    "Q4_OS_GRADED_SYMMETRY_RECEIPT",
    "Q4_OS_REFLECTION_POSITIVITY_RECEIPT",
    "Q4_OS_REGULARITY_AND_GROWTH_RECEIPT",
    "Q4_OS_LOCALITY_AND_CLUSTERING_RECEIPT",
    "Q4_OS_TIGHTNESS_AND_NONCOLLAPSE_RECEIPT",
    "Q4_OS_NONPERTURBATIVE_CONTINUUM_RECEIPT",
)


@dataclass(frozen=True)
class ReceiptRoute:
    """One complete alternative receipt route for a stage."""

    route_id: str
    receipts: tuple[str, ...]


@dataclass(frozen=True)
class StageSpec:
    """Static, topologically ordered physical-requirement stage."""

    stage_id: str
    all_dependencies: tuple[str, ...]
    any_dependency_groups: tuple[tuple[str, ...], ...]
    routes: tuple[ReceiptRoute, ...]
    claim_boundary: str


def _route(route_id: str, receipts: Sequence[str]) -> ReceiptRoute:
    return ReceiptRoute(route_id, tuple(receipts))


STAGE_SPECS: tuple[StageSpec, ...] = (
    StageSpec(
        "ROOT",
        (),
        (),
        (_route("immutable_root_packet", ROOT_RECEIPTS),),
        "A replayed immutable packet is a provenance root, not a physical emergence result.",
    ),
    StageSpec(
        "GEOMETRY_565",
        ("ROOT",),
        (),
        (
            _route("source_derived_geometry", GEOMETRY_SOURCE_RECEIPTS),
            _route("operational_hardware_geometry", GEOMETRY_HARDWARE_RECEIPTS),
        ),
        "The hardware route assumes the carrier geometry and does not derive issue #565 from legacy source axioms.",
    ),
    StageSpec(
        "CURRENT_566",
        ("ROOT", "GEOMETRY_565"),
        (),
        (_route("two_sided_current_tomography", CURRENT_RECEIPTS),),
        "A5 module arithmetic without ordered response tomography and same-current internal overlap holonomy does not pass this lane.",
    ),
    StageSpec(
        "GLOBAL_FORM_567",
        ("ROOT", "CURRENT_566"),
        (),
        (_route("physical_global_form", GLOBAL_FORM_RECEIPTS),),
        "The Z6 quotient and UV polarization must be derived from physical loops and readbacks.",
    ),
    StageSpec(
        "SPIN_EXCHANGE_314",
        ("ROOT", "GEOMETRY_565", "GLOBAL_FORM_567"),
        (),
        (_route("spin_exchange", SPIN_EXCHANGE_RECEIPTS),),
        "Exterior-algebra dimensions do not select the physical rank-15 module.",
    ),
    StageSpec(
        "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
        ("ROOT", "GLOBAL_FORM_567", "SPIN_EXCHANGE_314"),
        (),
        (_route("complete_source_registry", SOURCE_REGISTRY_RECEIPTS),),
        "Completeness is only inside the frozen finite source algebra.",
    ),
    StageSpec(
        "SCALAR_CHANNEL",
        ("ROOT", "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES"),
        (),
        (_route("primitive_scalar_search", SCALAR_RECEIPTS),),
        "This selects Higgs quantum numbers, not a Higgs phase or condensate.",
    ),
    StageSpec(
        "FAMILY_ATTACHMENT_569",
        (
            "ROOT",
            "CURRENT_566",
            "GLOBAL_FORM_567",
            "SPIN_EXCHANGE_314",
            "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
        ),
        (),
        (_route("full_screen_family_attachment", FAMILY_RECEIPTS),),
        "The rank-45 map must be derived from J_all; a chosen 15x3 reshape is inadmissible.",
    ),
    StageSpec(
        "Q1_LOCAL_ACTION",
        (
            "ROOT",
            "GLOBAL_FORM_567",
            "SPIN_EXCHANGE_314",
            "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
            "SCALAR_CHANNEL",
            "FAMILY_ATTACHMENT_569",
        ),
        (),
        (_route("finite_local_coupled_action", Q1_LOCAL_ACTION_RECEIPTS),),
        "The identical target-blind local action must be consumed by Q2; a schematic Lagrangian is insufficient.",
    ),
    StageSpec(
        "Q2_H",
        (
            "ROOT",
            "GLOBAL_FORM_567",
            "SPIN_EXCHANGE_314",
            "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
            "SCALAR_CHANNEL",
            "FAMILY_ATTACHMENT_569",
            "Q1_LOCAL_ACTION",
        ),
        (),
        (_route("finite_hamiltonian", Q2_H_RECEIPTS),),
        "A trivial gauge action, zero Gauss generators, or empty physical subspace fails.",
    ),
    StageSpec(
        "Q2_E",
        (
            "ROOT",
            "GLOBAL_FORM_567",
            "SPIN_EXCHANGE_314",
            "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
            "SCALAR_CHANNEL",
            "FAMILY_ATTACHMENT_569",
            "Q1_LOCAL_ACTION",
        ),
        (),
        (_route("finite_euclidean", Q2_E_RECEIPTS),),
        "This lane alone is finite Euclidean structure, not a positive/unitary quantum theory.",
    ),
    StageSpec(
        "POSITIVITY_OR_POSITIVE_TRANSFER",
        ("ROOT", "Q2_E"),
        (),
        (
            _route("reflection_positivity", POSITIVITY_REFLECTION_RECEIPTS),
            _route("positive_transfer_hamiltonian", POSITIVITY_TRANSFER_RECEIPTS),
        ),
        "One of reflection positivity or an equivalent positive transfer Hamiltonian is required.",
    ),
    StageSpec(
        "REFINEMENT_COMPLETENESS",
        (
            "ROOT",
            "GEOMETRY_565",
            "CURRENT_566",
            "GLOBAL_FORM_567",
            "SPIN_EXCHANGE_314",
            "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
            "SCALAR_CHANNEL",
            "FAMILY_ATTACHMENT_569",
            "Q1_LOCAL_ACTION",
        ),
        (("Q2_H",), ("Q2_E", "POSITIVITY_OR_POSITIVE_TRANSFER")),
        (_route("complement_complete_refinement", REFINEMENT_RECEIPTS),),
        "Old-sector convergence alone does not exclude new zero-energy sectors.",
    ),
    StageSpec(
        "PHYSICAL_IDENTIFICATION",
        ("ROOT", "REFINEMENT_COMPLETENESS"),
        (("Q2_H",), ("Q2_E", "POSITIVITY_OR_POSITIVE_TRANSFER")),
        (
            _route(
                "frozen_finite_alternative_class",
                PHYSICAL_IDENTIFICATION_FINITE_CLASS_RECEIPTS,
            ),
            _route(
                "infinite_class_reduced_by_checked_theorem",
                PHYSICAL_IDENTIFICATION_REDUCTION_THEOREM_RECEIPTS,
            ),
        ),
        "Without independent apparatus evidence the result is only an internally complete finite model.",
    ),
    StageSpec(
        "COMPLETE_COUPLED_DYNAMICS",
        (
            "ROOT",
            "Q1_LOCAL_ACTION",
            "REFINEMENT_COMPLETENESS",
            "PHYSICAL_IDENTIFICATION",
        ),
        (("Q2_H",), ("Q2_E", "POSITIVITY_OR_POSITIVE_TRANSFER")),
        (_route("complete_coupled_dynamics", COMPLETE_COUPLED_DYNAMICS_RECEIPTS),),
        "This binds one nontrivial interacting action, measure/evolution law, and refinement family across Q1 and Q2.",
    ),
    StageSpec(
        "FAMILY_BREAKING_OR_DESCENT",
        (
            "ROOT",
            "FAMILY_ATTACHMENT_569",
            "Q1_LOCAL_ACTION",
            "REFINEMENT_COMPLETENESS",
            "PHYSICAL_IDENTIFICATION",
        ),
        (),
        (
            _route(
                "exact_surviving_family_symmetry_compatible",
                FAMILY_EXACT_SYMMETRY_COMPATIBILITY_RECEIPTS,
            ),
            _route(
                "target_blind_family_breaking_or_descent",
                FAMILY_BREAKING_OR_DESCENT_RECEIPTS,
            ),
        ),
        "Exact surviving family symmetry may pass only if it permits the emitted matrices; otherwise target-blind breaking/descent is required.",
    ),
    StageSpec(
        "VERTEX_1PI",
        (
            "ROOT",
            "SCALAR_CHANNEL",
            "FAMILY_ATTACHMENT_569",
            "Q1_LOCAL_ACTION",
            "COMPLETE_COUPLED_DYNAMICS",
            "FAMILY_BREAKING_OR_DESCENT",
        ),
        (),
        (_route("amputated_1pi_vertex", VERTEX_RECEIPTS),),
        "Required for numerical Yukawa matrices and the full interacting claim.",
    ),
    StageSpec(
        "Q4_OS",
        (
            "ROOT",
            "Q2_E",
            "POSITIVITY_OR_POSITIVE_TRANSFER",
            "REFINEMENT_COMPLETENESS",
            "PHYSICAL_IDENTIFICATION",
            "COMPLETE_COUPLED_DYNAMICS",
            "FAMILY_BREAKING_OR_DESCENT",
            "VERTEX_1PI",
        ),
        (),
        (_route("osterwalder_schrader_continuum", Q4_OS_RECEIPTS),),
        "Optional and required only for a nonperturbative continuum Wightman claim.",
    ),
)

STAGE_IDS = tuple(spec.stage_id for spec in STAGE_SPECS)
STAGE_DAG_EDGES = tuple(
    (dependency, spec.stage_id)
    for spec in STAGE_SPECS
    for dependency in spec.all_dependencies
) + tuple(
    (dependency, spec.stage_id)
    for spec in STAGE_SPECS
    for group in spec.any_dependency_groups
    for dependency in group
)

BASE_GLOBAL_PASS_STAGES = (
    "ROOT",
    "GEOMETRY_565",
    "CURRENT_566",
    "GLOBAL_FORM_567",
    "SPIN_EXCHANGE_314",
    "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
    "SCALAR_CHANNEL",
    "FAMILY_ATTACHMENT_569",
    "Q1_LOCAL_ACTION",
    "REFINEMENT_COMPLETENESS",
    "PHYSICAL_IDENTIFICATION",
)

FULL_INTERACTING_PASS_STAGES = (
    "COMPLETE_COUPLED_DYNAMICS",
    "FAMILY_BREAKING_OR_DESCENT",
    "VERTEX_1PI",
)

ISSUE_CLOSURE_STAGES = {
    "565": ("GEOMETRY_565",),
    "566": ("GEOMETRY_565", "CURRENT_566"),
    "567": ("CURRENT_566", "GLOBAL_FORM_567", "SPIN_EXCHANGE_314"),
    "569": (
        "GEOMETRY_565",
        "CURRENT_566",
        "GLOBAL_FORM_567",
        "SPIN_EXCHANGE_314",
        "SOURCE_REGISTRY_AND_PRIMITIVE_RESIDUES",
        "FAMILY_ATTACHMENT_569",
        "PHYSICAL_IDENTIFICATION",
    ),
}

ISSUE_DELIMITATION_STATUS = {
    "590": {
        "status": "DELIMITATION_CLOSED",
        "physical_pass": None,
        "claim_boundary": (
            "Issue #590 is closed only by sharp delimitation. A stronger structural "
            "or interacting claim still requires the corresponding physical stages."
        ),
    }
}

REGISTERED_INVENTORY_VERIFIERS = (
    {
        "verifier_id": "physical_a5_sm_root_inventory_replayer_v1",
        "artifact_type": ROOT_REPORT_ARTIFACT_TYPE,
        "schema": ROOT_REPORT_SCHEMA,
        "receipt_keys": ROOT_INVENTORY_RECEIPTS,
    },
)

REGISTERED_PHYSICAL_PRODUCERS: tuple[Mapping[str, Any], ...] = ()

_KNOWN_DIAGNOSTIC_SCHEMAS = frozenset(
    {
        "PHYSICAL-FAMILY-POLE-RECEIPT-v1",
        "oph.efsm-conditional-theorem/v2",
        "oph_a5_twelve_port_structural_certificate_v2",
        "oph.a5-sm-theorem-evidence/1.0.0",
        "oph.a5-sm-theorem-application/2.0.0",
        "oph.source-to-eft/1.0.0",
    }
)

_NONSCIENTIFIC_LABEL_TOKENS = (
    "ASSUMPTION",
    "DEMO",
    "FORCED",
    "FROZEN_TARGET",
    "MOCK",
    "PLACEHOLDER",
    "SYNTHETIC",
    "TARGET_VALUE",
    "VISUALIZATION",
)


class RootManifestError(ValueError):
    """Strict ROOT manifest parsing or replay failure."""


@dataclass(frozen=True)
class RootArtifact:
    artifact_id: str
    path: str
    sha256: str
    media_type: str
    semantic_role: str
    provenance_kind: str


@dataclass(frozen=True)
class RootAncestryEdge:
    parent_artifact_id: str
    child_artifact_id: str
    operation_artifact_id: str


@dataclass(frozen=True)
class ParsedRootManifest:
    path: Path
    base_dir: Path
    raw_sha256: str
    packet_id: str
    code_commit: str
    dirty_tree_digest: str
    regulator_id: str
    boundary_condition_id: str
    superselection_sector_id: str
    random_seeds: tuple[int, ...]
    reproducibility: Mapping[str, Any]
    numerical_policy: Mapping[str, Any]
    artifacts: tuple[RootArtifact, ...]
    ancestry_edges: tuple[RootAncestryEdge, ...]
    candidate_domain_artifact_id: str
    candidate_domain_hash: str
    selection_law_artifact_id: str
    selection_law_hash: str
    target_selection_dependencies: Mapping[str, Any]


def verify_physical_a5_sm_stage_envelope(
    production_bundle_manifest_path: str | Path,
    *,
    expected_root_hash: str,
    expected_stage_id: str,
    expected_claim_scope: ClaimScope | str,
    expected_branch_id: str | None = None,
) -> dict[str, Any]:
    """Replay generic P0 identity for one stage without admitting science.

    The generic production envelope is deliberately only a provenance
    inventory.  Matching the root, lane, scope, and optional branch is a
    necessary admission condition, never a sufficient physical receipt.
    """

    if _SHA256_RE.fullmatch(expected_root_hash) is None:
        raise ValueError("expected_root_hash_must_be_sha256")
    if expected_stage_id not in STAGE_IDS:
        raise ValueError(f"unknown_expected_stage_id:{expected_stage_id}")
    try:
        scope = ClaimScope(expected_claim_scope)
    except ValueError as exc:
        raise ValueError(f"unsupported_claim_scope:{expected_claim_scope}") from exc
    if expected_branch_id is not None and (
        _IDENTIFIER_RE.fullmatch(expected_branch_id) is None
    ):
        raise ValueError("expected_branch_id_must_be_bounded_identifier")

    report = verify_production_bundle_manifest(production_bundle_manifest_path)
    return _stage_envelope_admission_from_report(
        report,
        expected_root_hash=expected_root_hash,
        expected_stage_id=expected_stage_id,
        expected_claim_scope=scope,
        expected_branch_id=expected_branch_id,
    )


def _stage_envelope_admission_from_report(
    report: Mapping[str, Any],
    *,
    expected_root_hash: str,
    expected_stage_id: str,
    expected_claim_scope: ClaimScope,
    expected_branch_id: str | None,
) -> dict[str, Any]:
    markers_match = bool(
        report.get("schema") == PRODUCTION_BUNDLE_REPORT_SCHEMA
        and report.get("artifact_type") == PRODUCTION_BUNDLE_REPORT_ARTIFACT_TYPE
    )
    bundle_inventory_passed = report.get("inventory_replay_passed") is True
    envelope_mapping = report.get("envelopes")
    if not isinstance(envelope_mapping, Mapping):
        envelope_mapping = {}
    common_rows = [
        (artifact_id, row)
        for artifact_id, row in envelope_mapping.items()
        if isinstance(artifact_id, str)
        and isinstance(row, Mapping)
        and row.get("profile") == COMMON_STAGE
    ]
    stage_rows = [
        (artifact_id, row)
        for artifact_id, row in common_rows
        if row.get("claim_lane") == expected_stage_id
    ]
    bundle_scope_coherent = bool(common_rows) and all(
        row.get("claim_scope") == expected_claim_scope.value
        for _, row in common_rows
    )
    bundle_stage_identity_coherent = bool(common_rows) and all(
        row.get("stage_id") in STAGE_IDS
        and row.get("claim_lane") == row.get("stage_id")
        for _, row in common_rows
    )
    unique_stage_row = len(stage_rows) == 1
    artifact_id: str | None = None
    row: Mapping[str, Any] = {}
    if unique_stage_row:
        artifact_id, row = stage_rows[0]
    row_inventory_passed = row.get("inventory_replay_passed") is True
    stage_identity_passed = row.get("stage_id") == expected_stage_id
    root_identity_passed = row.get("source_root_hash") == expected_root_hash
    scope_identity_passed = row.get("claim_scope") == expected_claim_scope.value
    branch_identity_passed = bool(
        expected_branch_id is None or row.get("branch_id") == expected_branch_id
    )
    generic_promotion_claimed = bool(
        report.get("scientific_replay_passed") is True
        or report.get("promotion_allowed") is True
        or any(
            row.get("scientific_replay_passed") is True
            or row.get("promotion_allowed") is True
            or row.get("producer_status_trusted") is True
            for _, row in common_rows
        )
    )
    nonscientific_labels = [
        value.upper()
        for candidate in (report, *(row for _, row in common_rows))
        for key in (
            "mode",
            "run_mode",
            "status",
            "epistemic_status",
            "receipt_type",
            "evidence_class",
        )
        if isinstance((value := candidate.get(key)), str)
    ]
    visualization_only = bool(
        report.get("visualization_only") is True
        or any(candidate.get("visualization_only") is True for _, candidate in common_rows)
        or any(
            token in label
            for label in nonscientific_labels
            for token in _NONSCIENTIFIC_LABEL_TOKENS
        )
    )
    identity_passed = bool(
        markers_match
        and bundle_inventory_passed
        and bundle_scope_coherent
        and bundle_stage_identity_coherent
        and unique_stage_row
        and row_inventory_passed
        and stage_identity_passed
        and root_identity_passed
        and scope_identity_passed
        and branch_identity_passed
        and not generic_promotion_claimed
        and not visualization_only
    )
    blockers: list[str] = []
    if not markers_match:
        blockers.append("production_bundle_report_marker_mismatch")
    if not bundle_inventory_passed:
        blockers.append("production_bundle_inventory_replay_not_passed")
    if not bundle_scope_coherent:
        blockers.append("production_bundle_mixed_or_wrong_claim_scope")
    if not bundle_stage_identity_coherent:
        blockers.append("production_bundle_stage_lane_identity_mismatch")
    if not unique_stage_row:
        blockers.append(f"common_stage_envelope_row_count:{len(stage_rows)}")
    if unique_stage_row and not row_inventory_passed:
        blockers.append("stage_envelope_inventory_replay_not_passed")
    if unique_stage_row and not stage_identity_passed:
        blockers.append("stage_envelope_stage_id_mismatch")
    if unique_stage_row and not root_identity_passed:
        blockers.append("stage_envelope_source_root_mismatch")
    if unique_stage_row and not scope_identity_passed:
        blockers.append("stage_envelope_claim_scope_mismatch")
    if unique_stage_row and not branch_identity_passed:
        blockers.append("stage_envelope_branch_mismatch")
    if generic_promotion_claimed:
        blockers.append("generic_inventory_claimed_unlicensed_scientific_promotion")
    if visualization_only:
        blockers.append("demo_assumption_is_visualization_only")
    blockers.append("no_registered_stage_scientific_producer")
    integrity_failed = bool(
        not markers_match
        or not bundle_inventory_passed
        or not bundle_scope_coherent
        or not bundle_stage_identity_coherent
        or len(stage_rows) > 1
        or (
            unique_stage_row
            and not (
                row_inventory_passed
                and stage_identity_passed
                and root_identity_passed
                and scope_identity_passed
                and branch_identity_passed
            )
        )
        or generic_promotion_claimed
    )
    status = (
        TerminalStatus.FAIL
        if integrity_failed
        else TerminalStatus.UNRESOLVED
        if visualization_only
        else TerminalStatus.OPEN
    )
    return {
        "status": status.value,
        "passed": status is TerminalStatus.PASS,
        "inventory_identity_passed": identity_passed,
        "scientific_admission_passed": False,
        "promotion_allowed": False,
        "visualization_only": visualization_only,
        "expected_root_hash": expected_root_hash,
        "expected_stage_id": expected_stage_id,
        "expected_claim_scope": expected_claim_scope.value,
        "expected_branch_id": expected_branch_id,
        "matched_envelope_artifact_id": artifact_id,
        "matched_envelope": dict(row),
        "generic_bundle_report": dict(report),
        "blockers": blockers,
        "claim_boundary": (
            "Generic P0 replay proves only same-root production identity. It cannot "
            "admit a physical stage without a registered stage-specific scientific "
            "verifier; DEMO_ASSUMPTION evidence remains visualization-only."
        ),
    }


def _missing_stage_envelope_admission(
    *,
    expected_root_hash: str | None,
    expected_stage_id: str,
    expected_claim_scope: ClaimScope,
) -> dict[str, Any]:
    return {
        "status": TerminalStatus.OPEN.value,
        "passed": False,
        "inventory_identity_passed": False,
        "scientific_admission_passed": False,
        "promotion_allowed": False,
        "visualization_only": False,
        "expected_root_hash": expected_root_hash,
        "expected_stage_id": expected_stage_id,
        "expected_claim_scope": expected_claim_scope.value,
        "expected_branch_id": None,
        "matched_envelope_artifact_id": None,
        "matched_envelope": {},
        "generic_bundle_report": None,
        "blockers": [
            "production_bundle_manifest_not_supplied",
            "no_registered_stage_scientific_producer",
        ],
        "claim_boundary": (
            "No generic production envelope was supplied; the stage remains OPEN."
        ),
    }


def verify_physical_a5_sm_root_manifest(manifest_path: str | Path) -> dict[str, Any]:
    """Replay one on-disk ROOT inventory without promoting it to physical ROOT."""

    try:
        parsed = _parse_root_manifest(manifest_path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, RootManifestError, TypeError) as exc:
        return _incomplete_root_report(
            f"root_manifest_parse_failed:{type(exc).__name__}:{exc}",
            str(manifest_path) if isinstance(manifest_path, (str, Path)) else None,
        )

    artifact_rows, artifact_values = _replay_root_artifacts(parsed)
    all_artifacts_valid = bool(artifact_rows) and all(
        row["passed"] is True for row in artifact_rows.values()
    )
    containment = bool(artifact_rows) and all(
        row["containment_passed"] is True for row in artifact_rows.values()
    )
    hash_replay = bool(artifact_rows) and all(
        row["hash_replay_passed"] is True for row in artifact_rows.values()
    )
    role_counts = {
        role: sum(artifact.semantic_role == role for artifact in parsed.artifacts)
        for role in REQUIRED_ROOT_ROLES
    }
    root_object_bindings = _verify_root_object_bindings(parsed, artifact_rows)
    role_census = all(count == 1 for count in role_counts.values()) and root_object_bindings[
        "passed"
    ]
    metadata = _verify_root_metadata(parsed)
    ancestry = _verify_source_ancestry(parsed)
    commitments = _verify_candidate_selection_commitments(parsed, artifact_rows)
    target_declarations = _verify_target_selection_declarations(parsed)

    receipts: dict[str, bool] = {
        "ROOT_MANIFEST_SCHEMA_RECEIPT": True,
        "ROOT_ARTIFACT_CONTAINMENT_RECEIPT": containment,
        "ROOT_ARTIFACT_HASH_REPLAY_RECEIPT": hash_replay and all_artifacts_valid,
        "ROOT_REQUIRED_OBJECT_CENSUS_RECEIPT": role_census,
        "ROOT_REPRODUCIBILITY_POLICY_RECEIPT": metadata["passed"],
        "ROOT_SOURCE_ANCESTRY_DAG_RECEIPT": ancestry["passed"],
        "ROOT_CANDIDATE_SELECTION_COMMITMENT_RECEIPT": commitments["passed"],
        "ROOT_TARGET_SELECTION_DECLARATIONS_RECEIPT": target_declarations["passed"],
    }
    inventory_passed = all(receipts[key] for key in ROOT_INVENTORY_RECEIPTS)
    receipts.update(
        {
            "ROOT_INVENTORY_REPLAY_RECEIPT": inventory_passed,
            # These are normative semantic/chronology checks.  The inventory
            # schema deliberately does not contain enough evidence to admit
            # any of them.
            "ROOT_TYPED_ROLE_SEMANTICS_RECEIPT": False,
            "ROOT_CODE_BUILD_BINDING_RECEIPT": False,
            "ROOT_PREOUTCOME_COMMITMENT_RECEIPT": False,
            "ROOT_IMMUTABLE_PACKET_REPLAY_RECEIPT": False,
        }
    )

    blockers: list[str] = []
    for artifact_id, row in artifact_rows.items():
        blockers.extend(f"artifact:{artifact_id}:{item}" for item in row["blockers"])
    blockers.extend(
        f"required_role_count:{role}:{count}"
        for role, count in role_counts.items()
        if count != 1
    )
    blockers.extend(metadata["blockers"])
    blockers.extend(root_object_bindings["blockers"])
    blockers.extend(ancestry["blockers"])
    blockers.extend(commitments["blockers"])
    blockers.extend(target_declarations["blockers"])

    open_requirements = [
        "no_registered_typed_role_semantics_verifier",
        "no_registered_code_build_binding_verifier",
        "no_registered_preoutcome_commitment_verifier",
    ]
    status = TerminalStatus.OPEN if inventory_passed else TerminalStatus.FAIL
    return {
        "schema": ROOT_REPORT_SCHEMA,
        "artifact_type": ROOT_REPORT_ARTIFACT_TYPE,
        "manifest_path": parsed.path.name,
        "manifest_sha256": parsed.raw_sha256,
        "packet_id": parsed.packet_id,
        "code_commit": parsed.code_commit,
        "dirty_tree_digest": parsed.dirty_tree_digest,
        "regulator_id": parsed.regulator_id,
        "artifact_rows": artifact_rows,
        "decoded_json_artifact_ids": sorted(artifact_values),
        "required_role_counts": role_counts,
        "root_object_bindings": root_object_bindings,
        "metadata_verification": metadata,
        "source_ancestry": ancestry,
        "candidate_selection_commitments": commitments,
        "target_selection_declarations": target_declarations,
        "receipts": receipts,
        "inventory_passed": inventory_passed,
        "status": status.value,
        "passed": status is TerminalStatus.PASS,
        "inventory_blockers": sorted(set(blockers)),
        "blockers": sorted(set(blockers + open_requirements)),
        "open_requirements": open_requirements,
        "claim_boundary": (
            "This is an inventory replay only. Hashes, declared target-selection fields, "
            "and a typed ancestry graph do not prove per-role semantics, bind the packet "
            "to the executable build, establish pre-outcome chronology, or prove any "
            "downstream physical lane."
        ),
    }


def write_physical_a5_sm_root_report(
    manifest_path: str | Path,
    output_path: str | Path | None = None,
) -> Path:
    """Write a deterministic ROOT verification report beside its manifest."""

    manifest = _resolve_regular_nonsymlink_file(Path(manifest_path))
    report = verify_physical_a5_sm_root_manifest(manifest)
    target = (
        Path(output_path)
        if output_path is not None
        else manifest.with_name("physical_a5_sm_root_verification.json")
    )
    bundle_root = manifest.parent.resolve(strict=True)
    try:
        target_parent = target.parent.resolve(strict=True)
    except OSError as exc:
        raise RootManifestError("root_report_output_parent_missing") from exc
    if target_parent != bundle_root:
        raise RootManifestError("root_report_output_must_be_beside_manifest")
    target = target_parent / target.name
    if target == manifest:
        raise RootManifestError("root_report_output_cannot_overwrite_manifest")
    if target.is_symlink():
        raise RootManifestError("root_report_output_is_symlink")
    if target.exists():
        raise RootManifestError("root_report_output_already_exists")
    try:
        with target.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    except FileExistsError as exc:
        raise RootManifestError("root_report_output_already_exists") from exc
    return target


def verify_physical_a5_sm_root_report_file(report_path: str | Path) -> dict[str, Any]:
    """Recompute a stored ROOT report and require exact report equality."""

    try:
        path = _resolve_regular_nonsymlink_file(Path(report_path))
        if path.stat().st_size > MAX_MANIFEST_BYTES:
            raise RootManifestError("root_report_exceeds_size_limit")
        candidate = _strict_json_loads(path.read_bytes())
        if not isinstance(candidate, Mapping):
            raise RootManifestError("root_report_must_be_object")
        if candidate.get("schema") != ROOT_REPORT_SCHEMA:
            raise RootManifestError("root_report_schema_mismatch")
        if candidate.get("artifact_type") != ROOT_REPORT_ARTIFACT_TYPE:
            raise RootManifestError("root_report_artifact_type_mismatch")
        manifest_name = candidate.get("manifest_path")
        if not isinstance(manifest_name, str):
            raise RootManifestError("root_report_manifest_path_missing")
        relative = Path(manifest_name)
        if relative.is_absolute() or ".." in relative.parts or len(relative.parts) != 1:
            raise RootManifestError("root_report_manifest_path_unsafe")
        recomputed = verify_physical_a5_sm_root_manifest(path.parent / relative)
        exact = candidate == recomputed
        inventory_replay_passed = exact and recomputed.get("inventory_passed") is True
        status = (
            TerminalStatus(recomputed["status"])
            if exact
            else TerminalStatus.FAIL
        )
        return {
            "status": status.value,
            "passed": exact and status is TerminalStatus.PASS,
            "inventory_replay_passed": inventory_replay_passed,
            "exact_report_replay": exact,
            "recomputed": recomputed,
            "blockers": [] if exact else ["stored_root_report_differs_from_replay"],
        }
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, RootManifestError, TypeError) as exc:
        return {
            "status": TerminalStatus.FAIL.value,
            "passed": False,
            "inventory_replay_passed": False,
            "exact_report_replay": False,
            "recomputed": None,
            "blockers": [f"root_report_replay_failed:{type(exc).__name__}:{exc}"],
        }


def verify_physical_a5_sm_requirements(
    root_manifest_path: str | Path,
    *,
    claim_scope: ClaimScope | str = ClaimScope.STRUCTURAL,
    production_bundle_manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """Evaluate the full physical requirements DAG from independently replayed evidence.

    The input is a ROOT manifest path, not a mapping of lane booleans.  JSON
    artifacts inside the root are inventoried, but no downstream receipt is
    admitted unless an in-code producer verifier is registered for its exact
    schema and artifact type.
    """

    try:
        scope = ClaimScope(claim_scope)
    except ValueError as exc:
        raise ValueError(f"unsupported_claim_scope:{claim_scope}") from exc

    root = verify_physical_a5_sm_root_manifest(root_manifest_path)
    admitted: dict[str, dict[str, str]] = {}
    evidence_inventory = _inventory_embedded_json_evidence(root_manifest_path, root)
    production_bundle_report = (
        verify_production_bundle_manifest(production_bundle_manifest_path)
        if production_bundle_manifest_path is not None
        else None
    )
    optional_stages = {
        ClaimScope.STRUCTURAL: set(FULL_INTERACTING_PASS_STAGES) | {"Q4_OS"},
        ClaimScope.FULL_INTERACTING: {"Q4_OS"},
        ClaimScope.CONTINUUM: set(),
    }[scope]
    stages: dict[str, dict[str, Any]] = {}
    for spec in STAGE_SPECS:
        if spec.stage_id == "ROOT":
            stage_envelope_admission = _missing_stage_envelope_admission(
                expected_root_hash=root.get("manifest_sha256"),
                expected_stage_id=spec.stage_id,
                expected_claim_scope=scope,
            )
            stage_envelope_admission["blockers"] = [
                "generic_common_stage_envelope_cannot_admit_normative_root",
                "no_registered_root_scientific_producer",
            ]
            stage_envelope_admission["claim_boundary"] = (
                "Generic P0 inventory cannot discharge typed ROOT semantics, build "
                "identity, or pre-outcome chronology."
            )
        elif production_bundle_report is None:
            stage_envelope_admission = _missing_stage_envelope_admission(
                expected_root_hash=root.get("manifest_sha256"),
                expected_stage_id=spec.stage_id,
                expected_claim_scope=scope,
            )
        else:
            stage_envelope_admission = _stage_envelope_admission_from_report(
                production_bundle_report,
                expected_root_hash=root.get("manifest_sha256"),
                expected_stage_id=spec.stage_id,
                expected_claim_scope=scope,
                expected_branch_id=None,
            )
        all_dependencies_passed = all(
            stages[dependency]["passed"] is True
            for dependency in spec.all_dependencies
        )
        if spec.any_dependency_groups:
            any_dependency_group_passed = any(
                all(stages[dependency]["passed"] is True for dependency in group)
                for group in spec.any_dependency_groups
            )
        else:
            any_dependency_group_passed = True

        route_rows: dict[str, dict[str, Any]] = {}
        for route in spec.routes:
            receipt_rows = {
                key: {
                    "admitted": key in admitted,
                    "producer_id": admitted.get(key, {}).get("producer_id"),
                }
                for key in route.receipts
            }
            route_passed = all(
                row["admitted"] is True for row in receipt_rows.values()
            )
            route_status = (
                TerminalStatus.NOT_APPLICABLE
                if spec.stage_id in optional_stages
                else TerminalStatus.PASS
                if route_passed
                else TerminalStatus.OPEN
            )
            route_rows[route.route_id] = {
                "receipts": receipt_rows,
                "status": route_status.value,
                "passed": route_status is TerminalStatus.PASS,
                "missing_receipts": [
                    key for key, row in receipt_rows.items() if row["admitted"] is not True
                ],
            }
        evidence_passed = any(row["passed"] is True for row in route_rows.values())
        production_identity_passed = bool(
            stage_envelope_admission.get("inventory_identity_passed") is True
        )
        blockers: list[str] = []
        blockers.extend(
            f"dependency_not_passed:{dependency}"
            for dependency in spec.all_dependencies
            if stages[dependency]["passed"] is not True
        )
        if spec.any_dependency_groups and not any_dependency_group_passed:
            blockers.append(
                "no_alternative_dependency_group_passed:"
                + "|".join("+".join(group) for group in spec.any_dependency_groups)
            )
        if not evidence_passed:
            missing = sorted(
                {key for route in route_rows.values() for key in route["missing_receipts"]}
            )
            blockers.extend(f"no_registered_physical_producer:{key}" for key in missing)
        if (
            spec.stage_id != "ROOT"
            and spec.stage_id not in optional_stages
            and not production_identity_passed
        ):
            blockers.append("same_root_production_identity_not_passed")
        if spec.stage_id == "ROOT":
            status = TerminalStatus(root.get("status", TerminalStatus.FAIL.value))
            blockers = list(root.get("blockers", []))
        elif spec.stage_id in optional_stages:
            status = TerminalStatus.NOT_APPLICABLE
            blockers = [f"not_applicable_to_claim_scope:{scope.value}"]
        elif stage_envelope_admission["status"] == TerminalStatus.FAIL.value:
            status = TerminalStatus.FAIL
            blockers.extend(
                f"production_envelope:{item}"
                for item in stage_envelope_admission["blockers"]
            )
        elif stage_envelope_admission["status"] == TerminalStatus.UNRESOLVED.value:
            status = TerminalStatus.UNRESOLVED
            blockers.extend(
                f"production_envelope:{item}"
                for item in stage_envelope_admission["blockers"]
            )
        else:
            dependency_statuses = [
                TerminalStatus(stages[dependency]["status"])
                for dependency in spec.all_dependencies
            ]
            if any(item is TerminalStatus.FAIL for item in dependency_statuses):
                status = TerminalStatus.FAIL
            elif any(
                item is TerminalStatus.UNRESOLVED for item in dependency_statuses
            ):
                status = TerminalStatus.UNRESOLVED
            elif _physical_stage_pass_authorized(
                all_dependencies_passed=all_dependencies_passed,
                any_dependency_group_passed=any_dependency_group_passed,
                evidence_passed=evidence_passed,
                production_identity_passed=production_identity_passed,
            ):
                status = TerminalStatus.PASS
            else:
                status = TerminalStatus.OPEN
        stages[spec.stage_id] = {
            "stage_id": spec.stage_id,
            "all_dependencies": list(spec.all_dependencies),
            "any_dependency_groups": [list(group) for group in spec.any_dependency_groups],
            "dependency_gate_passed": all_dependencies_passed and any_dependency_group_passed,
            "routes": route_rows,
            "production_envelope_admission": stage_envelope_admission,
            "production_identity_gate_passed": production_identity_passed,
            "evidence_gate_passed": evidence_passed,
            "status": status.value,
            "passed": status is TerminalStatus.PASS,
            "blockers": blockers,
            "claim_boundary": spec.claim_boundary,
        }

    claim_tiers = _compute_claim_tiers(stages)
    q2_physical_branch = claim_tiers["Q2_PHYSICAL_BRANCH_PASS"]
    structural_global_pass = claim_tiers["PHYSICAL_A5_SM_STRUCTURAL_GLOBAL_PASS"]
    full_interacting_pass = claim_tiers["PHYSICAL_A5_SM_FULL_INTERACTING_PASS"]
    continuum_pass = claim_tiers[
        "NONPERTURBATIVE_CONTINUUM_WIGHTMAN_PHYSICAL_PASS"
    ]
    issue_closure: dict[str, dict[str, Any]] = {
        issue: {
            "required_stages": list(required),
            "status": _aggregate_stage_statuses(
                [TerminalStatus(stages[stage_id]["status"]) for stage_id in required]
            ).value,
        }
        for issue, required in ISSUE_CLOSURE_STAGES.items()
    }
    for row in issue_closure.values():
        row["passed"] = row["status"] == TerminalStatus.PASS.value
    issue_closure.update({issue: dict(row) for issue, row in ISSUE_DELIMITATION_STATUS.items()})

    exact_small_status = TerminalStatus.OPEN
    scale_status = TerminalStatus.OPEN
    exact_small_checks = {
        check_id: {
            "status": TerminalStatus.OPEN.value,
            "passed": False,
            "producer_id": None,
            "verified_envelope_sha256": None,
        }
        for check_id in EXACT_SMALL_REQUIRED_CHECKS
    }
    scale_checks = {
        check_id: {
            "status": TerminalStatus.OPEN.value,
            "passed": False,
            "producer_id": None,
            "verified_envelope_sha256": None,
        }
        for check_id in SCALE_CAMPAIGN_REQUIRED_CHECKS
    }
    campaign_gates = {
        "EXACT_SMALL_ORACLE": {
            "status": exact_small_status.value,
            "passed": exact_small_status is TerminalStatus.PASS,
            "receipt": False,
            "required_check_order": list(EXACT_SMALL_REQUIRED_CHECKS),
            "required_checks": exact_small_checks,
            "all_required_checks_passed": False,
            "blockers": [
                f"required_exact_small_check_not_passed:{check_id}"
                for check_id in EXACT_SMALL_REQUIRED_CHECKS
            ],
        },
        "SCALE_CAMPAIGN_ALLOWED": {
            "status": scale_status.value,
            "passed": scale_status is TerminalStatus.PASS,
            "receipt": False,
            "required_check_order": list(SCALE_CAMPAIGN_REQUIRED_CHECKS),
            "required_checks": scale_checks,
            "all_required_checks_passed": False,
            "blockers": [
                f"required_scale_check_not_passed:{check_id}"
                for check_id in SCALE_CAMPAIGN_REQUIRED_CHECKS
            ],
        },
    }
    q2_status = _aggregate_alternative_statuses(
        (
            TerminalStatus(stages["Q2_H"]["status"]),
            _aggregate_stage_statuses(
                (
                    TerminalStatus(stages["Q2_E"]["status"]),
                    TerminalStatus(
                        stages["POSITIVITY_OR_POSITIVE_TRANSFER"]["status"]
                    ),
                )
            ),
        )
    )
    scope_required_stage_ids = _scope_required_stage_ids(scope)
    report_status = _report_status_for_scope(stages, scope=scope, q2_status=q2_status)

    return {
        "schema": REQUIREMENTS_REPORT_SCHEMA,
        "artifact_type": REQUIREMENTS_REPORT_ARTIFACT_TYPE,
        "root_verification": root,
        "claim_scope": scope.value,
        "scope_required_stage_ids": list(scope_required_stage_ids),
        "status": report_status.value,
        "registered_inventory_verifiers": [
            dict(row) for row in REGISTERED_INVENTORY_VERIFIERS
        ],
        "registered_physical_producers": [dict(row) for row in REGISTERED_PHYSICAL_PRODUCERS],
        "registered_physical_producer_count": len(REGISTERED_PHYSICAL_PRODUCERS),
        "registered_downstream_physical_producer_count": sum(
            row.get("stage_id") != "ROOT" for row in REGISTERED_PHYSICAL_PRODUCERS
        ),
        "stage_order": list(STAGE_IDS),
        "stage_dag_edges": [list(edge) for edge in STAGE_DAG_EDGES],
        "production_bundle_verification": production_bundle_report,
        "stages": stages,
        "embedded_evidence_inventory": evidence_inventory,
        "issue_closure": issue_closure,
        "nonphysical_campaign_gates": campaign_gates,
        "receipts": {
            "FINITE_EUCLIDEAN_STRUCTURAL_PASS": stages["Q2_E"]["passed"],
            "Q2_PHYSICAL_BRANCH_PASS": q2_physical_branch,
            "PHYSICAL_A5_SM_STRUCTURAL_GLOBAL_PASS": structural_global_pass,
            "PHYSICAL_A5_SM_GLOBAL_PASS": structural_global_pass,
            "PHYSICAL_A5_SM_FULL_INTERACTING_PASS": full_interacting_pass,
            "NUMERICAL_YUKAWA_MATRICES_PHYSICAL_PASS": full_interacting_pass,
            "NONPERTURBATIVE_CONTINUUM_WIGHTMAN_PHYSICAL_PASS": continuum_pass,
            "EXACT_SMALL_ORACLE_RECEIPT": False,
            "SCALE_CAMPAIGN_ALLOWED": False,
        },
        "claim_tiers": claim_tiers,
        "passed": report_status is TerminalStatus.PASS,
        "blockers": [
            f"required_stage_not_passed:{stage_id}"
            for stage_id in BASE_GLOBAL_PASS_STAGES
            if stages[stage_id]["passed"] is not True
        ]
        + ([] if q2_physical_branch else ["no_physical_Q2_branch_passed"])
        + [
            f"claim_scope_stage_not_passed:{stage_id}"
            for stage_id in scope_required_stage_ids
            if stages[stage_id]["passed"] is not True
        ],
        "full_interacting_blockers": [
            f"interacting_stage_not_passed:{stage_id}"
            for stage_id in FULL_INTERACTING_PASS_STAGES
            if stages[stage_id]["passed"] is not True
        ],
        "claim_boundary": (
            "OPEN, synthetic, theorem-only, structural, family/pole checklist, and arbitrary "
            "JSON receipts are diagnostic-only. Exact-small and scale admission are "
            "nonphysical campaign gates and are not part of the scientific conjunction. "
            "No stage can pass without a same-root verified production envelope from a "
            "registered producer."
        ),
    }


def _aggregate_stage_statuses(statuses: Sequence[TerminalStatus]) -> TerminalStatus:
    """Aggregate conjunctive statuses without inventing a success state."""

    if statuses and all(status is TerminalStatus.PASS for status in statuses):
        return TerminalStatus.PASS
    if any(status is TerminalStatus.FAIL for status in statuses):
        return TerminalStatus.FAIL
    if any(status is TerminalStatus.UNRESOLVED for status in statuses):
        return TerminalStatus.UNRESOLVED
    if statuses and all(status is TerminalStatus.NOT_APPLICABLE for status in statuses):
        return TerminalStatus.NOT_APPLICABLE
    return TerminalStatus.OPEN


def _physical_stage_pass_authorized(
    *,
    all_dependencies_passed: bool,
    any_dependency_group_passed: bool,
    evidence_passed: bool,
    production_identity_passed: bool,
) -> bool:
    """Require all four independently recomputed gates for a stage PASS."""

    return all(
        value is True
        for value in (
            all_dependencies_passed,
            any_dependency_group_passed,
            evidence_passed,
            production_identity_passed,
        )
    )


def _scope_required_stage_ids(scope: ClaimScope) -> tuple[str, ...]:
    if scope is ClaimScope.STRUCTURAL:
        return ()
    if scope is ClaimScope.FULL_INTERACTING:
        return FULL_INTERACTING_PASS_STAGES
    return (*FULL_INTERACTING_PASS_STAGES, "Q4_OS")


def _report_status_for_scope(
    stages: Mapping[str, Mapping[str, Any]],
    *,
    scope: ClaimScope,
    q2_status: TerminalStatus,
) -> TerminalStatus:
    """Aggregate exactly the claim tier selected by ``scope``."""

    statuses = [
        TerminalStatus(stages[stage_id]["status"])
        for stage_id in (
            *BASE_GLOBAL_PASS_STAGES,
            *_scope_required_stage_ids(scope),
        )
    ]
    return _aggregate_stage_statuses((*statuses, q2_status))


def _aggregate_alternative_statuses(statuses: Sequence[TerminalStatus]) -> TerminalStatus:
    """Aggregate disjunctive branches, preserving OPEN over premature failure."""

    if any(status is TerminalStatus.PASS for status in statuses):
        return TerminalStatus.PASS
    if any(status is TerminalStatus.OPEN for status in statuses):
        return TerminalStatus.OPEN
    if any(status is TerminalStatus.UNRESOLVED for status in statuses):
        return TerminalStatus.UNRESOLVED
    if statuses and all(status is TerminalStatus.NOT_APPLICABLE for status in statuses):
        return TerminalStatus.NOT_APPLICABLE
    return TerminalStatus.FAIL


def _compute_claim_tiers(stages: Mapping[str, Mapping[str, Any]]) -> dict[str, bool]:
    """Compute claim tiers only from already verifier-admitted stage statuses."""

    q2_physical_branch = bool(
        stages["Q2_H"]["passed"]
        or (
            stages["Q2_E"]["passed"]
            and stages["POSITIVITY_OR_POSITIVE_TRANSFER"]["passed"]
        )
    )
    structural_global_pass = bool(
        q2_physical_branch
        and all(stages[stage_id]["passed"] for stage_id in BASE_GLOBAL_PASS_STAGES)
    )
    full_interacting_pass = bool(
        structural_global_pass
        and all(stages[stage_id]["passed"] for stage_id in FULL_INTERACTING_PASS_STAGES)
    )
    return {
        "Q2_PHYSICAL_BRANCH_PASS": q2_physical_branch,
        "PHYSICAL_A5_SM_STRUCTURAL_GLOBAL_PASS": structural_global_pass,
        "PHYSICAL_A5_SM_FULL_INTERACTING_PASS": full_interacting_pass,
        "NONPERTURBATIVE_CONTINUUM_WIGHTMAN_PHYSICAL_PASS": bool(
            full_interacting_pass and stages["Q4_OS"]["passed"]
        ),
    }


def _parse_root_manifest(manifest_path: str | Path) -> ParsedRootManifest:
    if isinstance(manifest_path, Mapping):
        raise RootManifestError("root_manifest_must_be_an_on_disk_file")
    path = _resolve_regular_nonsymlink_file(Path(manifest_path))
    if path.stat().st_size > MAX_MANIFEST_BYTES:
        raise RootManifestError("root_manifest_exceeds_size_limit")
    raw = path.read_bytes()
    payload = _strict_json_loads(raw)
    if not isinstance(payload, Mapping):
        raise RootManifestError("root_manifest_must_be_object")
    _strict_keys(
        payload,
        required={
            "schema",
            "packet_id",
            "code_commit",
            "dirty_tree_digest",
            "regulator_id",
            "boundary_condition_id",
            "superselection_sector_id",
            "random_seeds",
            "reproducibility",
            "numerical_policy",
            "artifacts",
            "source_ancestry_edges",
            "candidate_domain_artifact_id",
            "candidate_domain_hash",
            "selection_law_artifact_id",
            "selection_law_hash",
            "target_selection_dependencies",
        },
        object_name="root manifest",
    )
    if payload["schema"] != ROOT_MANIFEST_SCHEMA:
        raise RootManifestError("root_manifest_schema_mismatch")

    artifacts_raw = payload["artifacts"]
    if not isinstance(artifacts_raw, list) or not (1 <= len(artifacts_raw) <= MAX_ARTIFACTS):
        raise RootManifestError("artifacts_must_be_a_bounded_nonempty_list")
    artifacts: list[RootArtifact] = []
    artifact_ids: set[str] = set()
    artifact_paths: set[str] = set()
    for index, row in enumerate(artifacts_raw):
        if not isinstance(row, Mapping):
            raise RootManifestError(f"artifact_{index}_must_be_object")
        _strict_keys(
            row,
            required={
                "artifact_id",
                "path",
                "sha256",
                "media_type",
                "semantic_role",
                "provenance_kind",
            },
            object_name=f"artifact[{index}]",
        )
        artifact = RootArtifact(
            artifact_id=_identifier(row["artifact_id"], f"artifact[{index}].artifact_id"),
            path=_relative_path_string(row["path"], f"artifact[{index}].path"),
            sha256=_hash(row["sha256"], f"artifact[{index}].sha256"),
            media_type=_media_type(row["media_type"], f"artifact[{index}].media_type"),
            semantic_role=_identifier(row["semantic_role"], f"artifact[{index}].semantic_role"),
            provenance_kind=_provenance_kind(
                row["provenance_kind"], f"artifact[{index}].provenance_kind"
            ),
        )
        if artifact.artifact_id in artifact_ids:
            raise RootManifestError(f"duplicate_artifact_id:{artifact.artifact_id}")
        if artifact.path in artifact_paths:
            raise RootManifestError(f"duplicate_artifact_path:{artifact.path}")
        artifact_ids.add(artifact.artifact_id)
        artifact_paths.add(artifact.path)
        artifacts.append(artifact)

    edges_raw = payload["source_ancestry_edges"]
    if not isinstance(edges_raw, list) or len(edges_raw) > MAX_ANCESTRY_EDGES:
        raise RootManifestError("source_ancestry_edges_must_be_a_bounded_list")
    edges: list[RootAncestryEdge] = []
    edge_keys: set[tuple[str, str, str]] = set()
    for index, row in enumerate(edges_raw):
        if not isinstance(row, Mapping):
            raise RootManifestError(f"source_ancestry_edge_{index}_must_be_object")
        _strict_keys(
            row,
            required={
                "parent_artifact_id",
                "child_artifact_id",
                "operation_artifact_id",
            },
            object_name=f"source_ancestry_edges[{index}]",
        )
        edge = RootAncestryEdge(
            _identifier(row["parent_artifact_id"], "parent_artifact_id"),
            _identifier(row["child_artifact_id"], "child_artifact_id"),
            _identifier(row["operation_artifact_id"], "operation_artifact_id"),
        )
        edge_key = (
            edge.parent_artifact_id,
            edge.child_artifact_id,
            edge.operation_artifact_id,
        )
        if edge_key in edge_keys:
            raise RootManifestError(f"duplicate_source_ancestry_edge:{edge_key}")
        edge_keys.add(edge_key)
        edges.append(edge)

    seeds_raw = payload["random_seeds"]
    if not isinstance(seeds_raw, list) or not seeds_raw or len(seeds_raw) > 4096:
        raise RootManifestError("random_seeds_must_be_a_bounded_nonempty_list")
    seeds: list[int] = []
    for index, value in enumerate(seeds_raw):
        if isinstance(value, bool) or not isinstance(value, int) or not (0 <= value < 2**63):
            raise RootManifestError(f"random_seeds[{index}]_invalid")
        seeds.append(value)
    if len(set(seeds)) != len(seeds):
        raise RootManifestError("random_seeds_must_be_unique")

    reproducibility = payload["reproducibility"]
    numerical_policy = payload["numerical_policy"]
    target_dependencies = payload["target_selection_dependencies"]
    if not isinstance(reproducibility, Mapping):
        raise RootManifestError("reproducibility_must_be_object")
    if not isinstance(numerical_policy, Mapping):
        raise RootManifestError("numerical_policy_must_be_object")
    if not isinstance(target_dependencies, Mapping):
        raise RootManifestError("target_selection_dependencies_must_be_object")

    code_commit = payload["code_commit"]
    if not isinstance(code_commit, str) or _COMMIT_RE.fullmatch(code_commit) is None:
        raise RootManifestError("code_commit_must_be_40_or_64_lowercase_hex")

    return ParsedRootManifest(
        path=path,
        base_dir=path.parent.resolve(strict=True),
        raw_sha256="sha256:" + hashlib.sha256(raw).hexdigest(),
        packet_id=_identifier(payload["packet_id"], "packet_id"),
        code_commit=code_commit,
        dirty_tree_digest=_hash(payload["dirty_tree_digest"], "dirty_tree_digest"),
        regulator_id=_identifier(payload["regulator_id"], "regulator_id"),
        boundary_condition_id=_identifier(
            payload["boundary_condition_id"], "boundary_condition_id"
        ),
        superselection_sector_id=_identifier(
            payload["superselection_sector_id"], "superselection_sector_id"
        ),
        random_seeds=tuple(seeds),
        reproducibility=dict(reproducibility),
        numerical_policy=dict(numerical_policy),
        artifacts=tuple(artifacts),
        ancestry_edges=tuple(edges),
        candidate_domain_artifact_id=_identifier(
            payload["candidate_domain_artifact_id"], "candidate_domain_artifact_id"
        ),
        candidate_domain_hash=_hash(payload["candidate_domain_hash"], "candidate_domain_hash"),
        selection_law_artifact_id=_identifier(
            payload["selection_law_artifact_id"], "selection_law_artifact_id"
        ),
        selection_law_hash=_hash(payload["selection_law_hash"], "selection_law_hash"),
        target_selection_dependencies=dict(target_dependencies),
    )


def _replay_root_artifacts(
    parsed: ParsedRootManifest,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    decoded_json: dict[str, Any] = {}
    resolved_paths: set[Path] = set()
    for artifact in parsed.artifacts:
        blockers: list[str] = []
        containment_passed = False
        hash_replay_passed = False
        actual_sha256: str | None = None
        byte_length: int | None = None
        resolved_relative_path: str | None = None
        try:
            path = _resolve_contained_artifact(parsed.base_dir, artifact.path)
            containment_passed = True
            if path in resolved_paths:
                raise RootManifestError("duplicate_resolved_artifact_path")
            resolved_paths.add(path)
            byte_length = path.stat().st_size
            if byte_length > MAX_ARTIFACT_BYTES:
                raise RootManifestError("artifact_exceeds_size_limit")
            raw = path.read_bytes()
            actual_sha256 = "sha256:" + hashlib.sha256(raw).hexdigest()
            if actual_sha256 != artifact.sha256:
                raise RootManifestError("declared_sha256_mismatch")
            hash_replay_passed = True
            resolved_relative_path = path.relative_to(parsed.base_dir).as_posix()
            if artifact.media_type == "application/json":
                if byte_length > MAX_JSON_ARTIFACT_BYTES:
                    raise RootManifestError("json_artifact_exceeds_decode_size_limit")
                decoded_json[artifact.artifact_id] = _strict_json_loads(raw)
        except (OSError, UnicodeDecodeError, json.JSONDecodeError, RootManifestError) as exc:
            blockers.append(str(exc))
        rows[artifact.artifact_id] = {
            "artifact_id": artifact.artifact_id,
            "declared_path": artifact.path,
            "declared_sha256": artifact.sha256,
            "actual_sha256": actual_sha256,
            "byte_length": byte_length,
            "media_type": artifact.media_type,
            "semantic_role": artifact.semantic_role,
            "provenance_kind": artifact.provenance_kind,
            "resolved_relative_path": resolved_relative_path,
            "containment_passed": containment_passed,
            "hash_replay_passed": hash_replay_passed,
            "passed": containment_passed and hash_replay_passed and not blockers,
            "blockers": blockers,
        }
    return rows, decoded_json


def _verify_root_object_bindings(
    parsed: ParsedRootManifest,
    artifact_rows: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Bind scalar ROOT commitments to their uniquely typed artifact rows."""

    blockers: list[str] = []
    by_role: dict[str, list[RootArtifact]] = {}
    for artifact in parsed.artifacts:
        by_role.setdefault(artifact.semantic_role, []).append(artifact)
    dirty_rows = by_role.get("dirty_tree_manifest", [])
    dirty_binding_passed = bool(
        len(dirty_rows) == 1
        and dirty_rows[0].sha256 == parsed.dirty_tree_digest
        and artifact_rows.get(dirty_rows[0].artifact_id, {}).get("actual_sha256")
        == parsed.dirty_tree_digest
    )
    if not dirty_binding_passed:
        blockers.append("dirty_tree_digest_not_bound_to_dirty_tree_manifest")
    provenance_rows: dict[str, bool] = {}
    for role in REQUIRED_ROOT_ROLES:
        rows = by_role.get(role, [])
        expected_kind = (
            "source_operation" if role == "dynamics_generator" else "source_primitive"
        )
        passed = len(rows) == 1 and rows[0].provenance_kind == expected_kind
        provenance_rows[role] = passed
        if not passed:
            blockers.append(f"required_role_provenance_mismatch:{role}:{expected_kind}")
    return {
        "passed": not blockers,
        "dirty_tree_digest_binding": dirty_binding_passed,
        "required_role_provenance_typing": provenance_rows,
        "blockers": blockers,
    }


def _verify_root_metadata(parsed: ParsedRootManifest) -> dict[str, Any]:
    blockers: list[str] = []
    expected_reproducibility = {
        "mode",
        "seed_schedule",
        "parallel_reduction",
        "blas_threads_per_worker",
    }
    if set(parsed.reproducibility) != expected_reproducibility:
        blockers.append("reproducibility_schema_mismatch")
    else:
        if parsed.reproducibility["mode"] not in {
            "deterministic",
            "seeded_stochastic",
        }:
            blockers.append("reproducibility_mode_invalid")
        if parsed.reproducibility["seed_schedule"] != "explicit_manifest_order":
            blockers.append("seed_schedule_not_explicit_manifest_order")
        if parsed.reproducibility["parallel_reduction"] not in {
            "deterministic",
            "not_applicable",
        }:
            blockers.append("parallel_reduction_policy_invalid")
        threads = parsed.reproducibility["blas_threads_per_worker"]
        if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
            blockers.append("blas_threads_per_worker_invalid")

    expected_numerical = {
        "backend",
        "backend_version",
        "precision_bits",
        "interval_method",
        "rank_certification_method",
    }
    if set(parsed.numerical_policy) != expected_numerical:
        blockers.append("numerical_policy_schema_mismatch")
    else:
        for field_name in ("backend", "backend_version"):
            value = parsed.numerical_policy[field_name]
            if not isinstance(value, str) or not value or len(value) > 128:
                blockers.append(f"numerical_policy_{field_name}_invalid")
        precision = parsed.numerical_policy["precision_bits"]
        if isinstance(precision, bool) or not isinstance(precision, int) or precision < 1:
            blockers.append("numerical_policy_precision_bits_invalid")
        for field_name in ("interval_method", "rank_certification_method"):
            value = parsed.numerical_policy[field_name]
            if not isinstance(value, str) or _IDENTIFIER_RE.fullmatch(value) is None:
                blockers.append(f"numerical_policy_{field_name}_invalid")
    return {"passed": not blockers, "blockers": blockers}


def _verify_source_ancestry(parsed: ParsedRootManifest) -> dict[str, Any]:
    by_id = {artifact.artifact_id: artifact for artifact in parsed.artifacts}
    blockers: list[str] = []
    parents: dict[str, set[str]] = {artifact_id: set() for artifact_id in by_id}
    children: dict[str, set[str]] = {artifact_id: set() for artifact_id in by_id}
    serialized_edges: list[dict[str, str]] = []
    for edge in parsed.ancestry_edges:
        missing = [
            artifact_id
            for artifact_id in (
                edge.parent_artifact_id,
                edge.child_artifact_id,
                edge.operation_artifact_id,
            )
            if artifact_id not in by_id
        ]
        if missing:
            blockers.append("ancestry_edge_missing_artifact:" + ",".join(sorted(set(missing))))
            continue
        if edge.parent_artifact_id == edge.child_artifact_id:
            blockers.append(f"ancestry_self_loop:{edge.child_artifact_id}")
        if by_id[edge.parent_artifact_id].provenance_kind == "source_operation":
            blockers.append(
                f"ancestry_parent_is_source_operation:{edge.parent_artifact_id}"
            )
        if by_id[edge.child_artifact_id].provenance_kind == "source_operation":
            blockers.append(
                f"ancestry_child_is_source_operation:{edge.child_artifact_id}"
            )
        operation = by_id[edge.operation_artifact_id]
        if operation.provenance_kind != "source_operation":
            blockers.append(
                f"ancestry_operation_not_source_operation:{edge.operation_artifact_id}"
            )
        parents[edge.child_artifact_id].add(edge.parent_artifact_id)
        children[edge.parent_artifact_id].add(edge.child_artifact_id)
        serialized_edges.append(
            {
                "parent_artifact_id": edge.parent_artifact_id,
                "child_artifact_id": edge.child_artifact_id,
                "operation_artifact_id": edge.operation_artifact_id,
            }
        )

    for artifact in parsed.artifacts:
        incoming = parents[artifact.artifact_id]
        if artifact.provenance_kind in {"source_primitive", "source_operation"} and incoming:
            blockers.append(f"{artifact.provenance_kind}_has_ancestry:{artifact.artifact_id}")
        if artifact.provenance_kind == "derived_array" and not incoming:
            blockers.append(f"derived_array_missing_ancestry:{artifact.artifact_id}")

    operation_ids = {
        artifact.artifact_id
        for artifact in parsed.artifacts
        if artifact.provenance_kind == "source_operation"
    }
    referenced_operation_ids = {
        edge.operation_artifact_id
        for edge in parsed.ancestry_edges
        if edge.operation_artifact_id in by_id
    }
    for operation_id in sorted(operation_ids - referenced_operation_ids):
        blockers.append(f"source_operation_not_referenced:{operation_id}")

    indegree = {artifact_id: len(parent_ids) for artifact_id, parent_ids in parents.items()}
    frontier = sorted(artifact_id for artifact_id, degree in indegree.items() if degree == 0)
    visited: list[str] = []
    while frontier:
        node = frontier.pop(0)
        visited.append(node)
        for child in sorted(children[node]):
            indegree[child] -= 1
            if indegree[child] == 0:
                frontier.append(child)
                frontier.sort()
    acyclic = len(visited) == len(by_id)
    if not acyclic:
        blockers.append("source_ancestry_cycle_detected")

    primitive_ids = {
        artifact.artifact_id
        for artifact in parsed.artifacts
        if artifact.provenance_kind == "source_primitive"
    }
    reachable = set(primitive_ids)
    for node in visited:
        if node in reachable:
            reachable.update(children[node])
    for artifact in parsed.artifacts:
        if artifact.provenance_kind == "derived_array" and artifact.artifact_id not in reachable:
            blockers.append(f"derived_array_not_source_reachable:{artifact.artifact_id}")

    return {
        "passed": not blockers,
        "acyclic": acyclic,
        "source_primitive_ids": sorted(primitive_ids),
        "source_operation_ids": sorted(operation_ids),
        "topological_order": visited if acyclic else None,
        "edges": serialized_edges,
        "blockers": sorted(set(blockers)),
    }


def _verify_candidate_selection_commitments(
    parsed: ParsedRootManifest,
    artifact_rows: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    blockers: list[str] = []
    by_id = {artifact.artifact_id: artifact for artifact in parsed.artifacts}
    checks = (
        (
            "candidate_domain",
            parsed.candidate_domain_artifact_id,
            parsed.candidate_domain_hash,
        ),
        ("selection_law", parsed.selection_law_artifact_id, parsed.selection_law_hash),
    )
    rows: dict[str, Any] = {}
    for expected_role, artifact_id, declared_hash in checks:
        artifact = by_id.get(artifact_id)
        artifact_row = artifact_rows.get(artifact_id, {})
        passed = bool(
            artifact is not None
            and artifact.semantic_role == expected_role
            and artifact.provenance_kind == "source_primitive"
            and artifact.sha256 == declared_hash
            and artifact_row.get("actual_sha256") == declared_hash
            and artifact_row.get("passed") is True
        )
        if not passed:
            blockers.append(f"{expected_role}_commitment_mismatch")
        rows[expected_role] = {
            "artifact_id": artifact_id,
            "declared_hash": declared_hash,
            "passed": passed,
        }
    if parsed.candidate_domain_artifact_id == parsed.selection_law_artifact_id:
        blockers.append("candidate_domain_and_selection_law_must_be_distinct")
    return {"passed": not blockers, "rows": rows, "blockers": blockers}


def _verify_target_selection_declarations(parsed: ParsedRootManifest) -> dict[str, Any]:
    blockers: list[str] = []
    if set(parsed.target_selection_dependencies) != set(TARGET_SELECTION_FIELDS):
        blockers.append("target_selection_dependency_schema_mismatch")
    rows: dict[str, bool] = {}
    for field_name in TARGET_SELECTION_FIELDS:
        value = parsed.target_selection_dependencies.get(field_name)
        passed = value is False
        rows[field_name] = passed
        if not passed:
            blockers.append(f"target_selection_dependency_not_false:{field_name}")
    return {"passed": not blockers, "rows": rows, "blockers": blockers}


def _inventory_embedded_json_evidence(
    root_manifest_path: str | Path,
    root_report: Mapping[str, Any],
) -> list[dict[str, Any]]:
    if root_report.get("inventory_passed") is not True:
        return []
    try:
        parsed = _parse_root_manifest(root_manifest_path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, RootManifestError, TypeError):
        return []
    inventory: list[dict[str, Any]] = []
    for artifact in parsed.artifacts:
        if artifact.media_type != "application/json":
            continue
        row = root_report.get("artifact_rows", {}).get(artifact.artifact_id, {})
        if row.get("passed") is not True:
            continue
        try:
            path = _resolve_contained_artifact(parsed.base_dir, artifact.path)
            payload = _strict_json_loads(path.read_bytes())
            classification, reason = _classify_json_evidence(payload)
        except (OSError, UnicodeDecodeError, json.JSONDecodeError, RootManifestError) as exc:
            classification = "unadmitted"
            reason = f"json_inventory_failed:{type(exc).__name__}:{exc}"
        inventory.append(
            {
                "artifact_id": artifact.artifact_id,
                "semantic_role": artifact.semantic_role,
                "classification": classification,
                "physical_receipt_admission": False,
                "reason": reason,
            }
        )
    return inventory


def _classify_json_evidence(payload: Any) -> tuple[str, str]:
    if not isinstance(payload, Mapping):
        return "unadmitted", "json_evidence_is_not_an_object"
    labels = [
        value.upper()
        for key in (
            "mode",
            "run_mode",
            "status",
            "epistemic_status",
            "receipt_type",
            "evidence_class",
        )
        if isinstance((value := payload.get(key)), str)
    ]
    if payload.get("visualization_only") is True or any(
        token in label
        for label in labels
        for token in _NONSCIENTIFIC_LABEL_TOKENS
    ):
        return (
            "diagnostic_only",
            "demo_assumption_is_visualization_only_and_cannot_promote_science",
        )
    schema = payload.get("schema") or payload.get("schema_version")
    if schema in _KNOWN_DIAGNOSTIC_SCHEMAS:
        if schema == "PHYSICAL-FAMILY-POLE-RECEIPT-v1":
            return (
                "diagnostic_only",
                "survival_proof_family_pole_contract_is_not_a_registered_physical_producer",
            )
        return "diagnostic_only", "known_theorem_or_structural_schema_is_nonpromoting"
    open_markers = (
        payload.get("receipt_verdict"),
        payload.get("status"),
        payload.get("verdict"),
    )
    if any(isinstance(value, str) and value.upper() == "OPEN" for value in open_markers):
        return "diagnostic_only", "open_receipt_is_not_evidence"
    if payload.get("synthetic") is True or payload.get("fixture") is True:
        return "diagnostic_only", "synthetic_fixture_is_not_evidence"
    boolean_receipt_keys = sorted(
        key
        for key, value in payload.items()
        if isinstance(key, str)
        and (key.endswith("_RECEIPT") or key.endswith("_PASS") or key == "passed")
        and isinstance(value, bool)
    )
    if boolean_receipt_keys:
        return (
            "unadmitted",
            "caller_authored_receipt_booleans_have_no_registered_replay_verifier",
        )
    return "unadmitted", "artifact_schema_has_no_registered_physical_producer"


def _incomplete_root_report(blocker: str, manifest_path: str | None) -> dict[str, Any]:
    receipts = {key: False for key in ROOT_RECEIPTS}
    return {
        "schema": ROOT_REPORT_SCHEMA,
        "artifact_type": ROOT_REPORT_ARTIFACT_TYPE,
        "manifest_path": manifest_path,
        "manifest_sha256": None,
        "packet_id": None,
        "artifact_rows": {},
        "decoded_json_artifact_ids": [],
        "required_role_counts": {role: 0 for role in REQUIRED_ROOT_ROLES},
        "root_object_bindings": {"passed": False, "blockers": [blocker]},
        "metadata_verification": {"passed": False, "blockers": [blocker]},
        "source_ancestry": {"passed": False, "blockers": [blocker]},
        "candidate_selection_commitments": {"passed": False, "blockers": [blocker]},
        "target_selection_declarations": {"passed": False, "blockers": [blocker]},
        "receipts": receipts,
        "inventory_passed": False,
        "status": TerminalStatus.FAIL.value,
        "passed": False,
        "inventory_blockers": [blocker],
        "blockers": [blocker],
        "open_requirements": [],
        "claim_boundary": "No ROOT receipt is emitted from malformed or off-disk evidence.",
    }


def _resolve_regular_nonsymlink_file(path: Path) -> Path:
    if path.is_symlink():
        raise RootManifestError("path_is_symlink")
    resolved = path.resolve(strict=True)
    if not resolved.is_file():
        raise RootManifestError("path_is_not_regular_file")
    return resolved


def _resolve_contained_artifact(base_dir: Path, declared_path: str) -> Path:
    relative = Path(declared_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise RootManifestError("artifact_path_is_absolute_or_contains_parent_traversal")
    cursor = base_dir
    for part in relative.parts:
        cursor = cursor / part
        if cursor.is_symlink():
            raise RootManifestError("artifact_path_contains_symlink")
    resolved = cursor.resolve(strict=True)
    try:
        resolved.relative_to(base_dir)
    except ValueError as exc:
        raise RootManifestError("artifact_path_escapes_bundle_root") from exc
    if not resolved.is_file():
        raise RootManifestError("artifact_path_is_not_regular_file")
    return resolved


def _strict_json_loads(raw: bytes) -> Any:
    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise RootManifestError(f"duplicate_json_key:{key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> Any:
        raise RootManifestError(f"nonfinite_json_constant:{value}")

    value = json.loads(
        raw.decode("utf-8"),
        object_pairs_hook=pairs_hook,
        parse_constant=reject_constant,
    )
    _validate_finite_json(value)
    return value


def _validate_finite_json(value: Any) -> None:
    """Reject overflow-to-infinity floats as well as explicit NaN tokens."""

    if isinstance(value, float) and not math.isfinite(value):
        raise RootManifestError("nonfinite_json_number")
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str):
                raise RootManifestError("json_object_key_not_string")
            _validate_finite_json(child)
    elif isinstance(value, list):
        for child in value:
            _validate_finite_json(child)


def _strict_keys(
    payload: Mapping[str, Any],
    *,
    required: set[str],
    object_name: str,
) -> None:
    actual = set(payload)
    missing = sorted(required - actual)
    unexpected = sorted(actual - required)
    if missing:
        raise RootManifestError(f"{object_name}_missing_fields:{missing}")
    if unexpected:
        raise RootManifestError(f"{object_name}_unexpected_fields:{unexpected}")


def _identifier(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or _IDENTIFIER_RE.fullmatch(value) is None:
        raise RootManifestError(f"{field_name}_must_be_bounded_identifier")
    return value


def _hash(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise RootManifestError(f"{field_name}_must_be_sha256")
    return value


def _relative_path_string(value: Any, field_name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 512
        or "\\" in value
    ):
        raise RootManifestError(f"{field_name}_must_be_bounded_relative_path")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value:
        raise RootManifestError(f"{field_name}_must_be_normalized_relative_path")
    return value


def _media_type(value: Any, field_name: str) -> str:
    allowed = {
        "application/json",
        "application/octet-stream",
        "application/x-npy",
        "application/x-npz",
        "text/plain",
    }
    if value not in allowed:
        raise RootManifestError(f"{field_name}_unsupported")
    return str(value)


def _provenance_kind(value: Any, field_name: str) -> str:
    if value not in {"source_primitive", "source_operation", "derived_array"}:
        raise RootManifestError(f"{field_name}_invalid")
    return str(value)


__all__ = [
    "BASE_GLOBAL_PASS_STAGES",
    "ClaimScope",
    "EXACT_SMALL_REQUIRED_CHECKS",
    "FAMILY_BREAKING_OR_DESCENT_RECEIPTS",
    "FAMILY_EXACT_SYMMETRY_COMPATIBILITY_RECEIPTS",
    "FULL_INTERACTING_PASS_STAGES",
    "ISSUE_CLOSURE_STAGES",
    "ISSUE_DELIMITATION_STATUS",
    "Q2_E_RECEIPTS",
    "Q2_H_RECEIPTS",
    "REQUIRED_ROOT_ROLES",
    "REGISTERED_INVENTORY_VERIFIERS",
    "REGISTERED_PHYSICAL_PRODUCERS",
    "REQUIREMENTS_REPORT_ARTIFACT_TYPE",
    "REQUIREMENTS_REPORT_SCHEMA",
    "ROOT_MANIFEST_SCHEMA",
    "ROOT_INVENTORY_RECEIPTS",
    "ROOT_RECEIPTS",
    "ROOT_REPORT_ARTIFACT_TYPE",
    "ROOT_REPORT_SCHEMA",
    "SCALE_CAMPAIGN_REQUIRED_CHECKS",
    "STAGE_DAG_EDGES",
    "STAGE_IDS",
    "STAGE_SPECS",
    "TARGET_SELECTION_FIELDS",
    "TerminalStatus",
    "RootManifestError",
    "verify_physical_a5_sm_requirements",
    "verify_physical_a5_sm_root_manifest",
    "verify_physical_a5_sm_root_report_file",
    "verify_physical_a5_sm_stage_envelope",
    "write_physical_a5_sm_root_report",
]
