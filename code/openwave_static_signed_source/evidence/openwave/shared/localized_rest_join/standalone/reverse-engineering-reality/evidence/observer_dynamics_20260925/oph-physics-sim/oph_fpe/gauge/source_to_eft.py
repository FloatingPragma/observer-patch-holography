"""Fail-closed source-to-EFT contracts for the OPH simulator.

This module deliberately stops before FJ, BRST, and pole production.  It
implements the exact finite mathematics that can be checked on the source
side without importing Standard Model target values:

* typed quotient-visible moment packets;
* an exact rational affine identifiability theorem, certified by a left
  inverse in the infinity norm;
* deterministic/stochastic source-law separation;
* coherent provenance commitments and perturbative masks; and
* exact affine matching-tower Jacobian and error composition.

The mathematical receipts below are not physical receipts.  In particular,
no caller-provided parameter vector, hash, covariance, or boolean can promote
``MATCH1``, ``FJ1``, ``COV1``, ``BRST1``, or a physical pole.  Those receipts
remain false until independent primitive producers and replay verifiers are
implemented.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import re
from dataclasses import dataclass, field
from enum import Enum
from fractions import Fraction
from typing import Any, Mapping, Sequence


SOURCE_TO_EFT_SCHEMA_VERSION = "oph.source-to-eft/1.0.0"
QUOTIENT_MOMENT_SCHEMA_VERSION = "oph.quotient-moment-packet/1.0.0"
SOURCE_LAW_SCHEMA_VERSION = "oph.finite-source-law/1.0.0"

CONDITIONAL_POLE_STATUS = "CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL"

_SHA256_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
_IDENTIFIER_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.:/+-]{0,127}$")

COMMITMENT_FIELDS = (
    "source_root_hash",
    "branch_hash",
    "field_census_hash",
    "scheme_hash",
    "threshold_hash",
    "fj_convention_hash",
    "perturbative_mask_hash",
    "analytic_sheet_hash",
    "units_clock_hash",
)

REQUIRED_SOURCE_COMPONENTS = (
    "oph_source",
    "clock",
    "d10",
    "d11",
    "matching",
    "fj",
    "source_law",
    "brst",
    "analytic_sheet",
)

_PHYSICAL_RECEIPT_KEYS = (
    "COHERENT_PHYSICAL_SOURCE_PACKET_RECEIPT",
    "SOURCE_PARAMETER_CANDIDATE_PRODUCER_RECEIPT",
    "SOURCE_LAW_PRIMITIVE_PRODUCER_RECEIPT",
    "OPH_SOURCE_PRIMITIVE_PRODUCER_RECEIPT",
    "OPH_QUOTIENT_MOMENT_PHYSICAL_RECEIPT",
    "OPH_NATIVE_EFT_REALIZATION_RECEIPT",
    "MATCH1_RECEIPT",
    "FJ1_RECEIPT",
    "COV1_RECEIPT",
    "BRST1_RECEIPT",
    "PHYSICAL_POLE_RECEIPT",
    "PHYSICAL_COHOMOLOGY_POLE_RECEIPT",
    "OPH_NATIVE_STRICT_1L_PHYSICAL_RECEIPT",
)


class MomentSemantics(str, Enum):
    """How quotient moments were defined, not how they are interpreted later."""

    DETERMINISTIC_SETTLED_STATE = "deterministic_settled_state"
    STOCHASTIC_SOURCE_LAW = "stochastic_source_law"


class SourceLawSemantics(str, Enum):
    """Mutually exclusive probability semantics for a finite source packet."""

    DETERMINISTIC = "deterministic"
    STOCHASTIC = "stochastic"


class CalculationLane(str, Enum):
    """The imported validation lane and OPH-native prediction lane never mix."""

    IMPORTED_SM_STRICT_1L_VALIDATION = "IMPORTED_SM_STRICT_1L_VALIDATION"
    OPH_NATIVE_STRICT_1L = "OPH_NATIVE_STRICT_1L"


def _require_identifier(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or _IDENTIFIER_RE.fullmatch(value) is None:
        raise ValueError(f"{field_name} must be a bounded identifier")
    return value


def _require_hash(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise ValueError(f"{field_name} must be sha256:<64 lowercase hex>")
    return value


def _fraction(value: Any, field_name: str) -> Fraction:
    """Parse exact rational input while rejecting binary floating-point values."""

    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{field_name} must be an exact integer, fraction string, or Fraction")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        try:
            parsed = Fraction(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{field_name} is not an exact rational") from exc
        return parsed
    raise TypeError(f"{field_name} must be an exact integer, fraction string, or Fraction")


def _fraction_tuple(values: Sequence[Any], field_name: str) -> tuple[Fraction, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{field_name} must be a sequence of exact rationals")
    return tuple(_fraction(value, f"{field_name}[{index}]") for index, value in enumerate(values))


def _qstr(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _canonical_hash(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _strict_keys(
    payload: Mapping[str, Any],
    *,
    required: set[str],
    optional: set[str] | None = None,
    object_name: str,
) -> None:
    optional = optional or set()
    actual = set(payload)
    missing = sorted(required - actual)
    unexpected = sorted(actual - required - optional)
    if missing:
        raise ValueError(f"{object_name} missing required fields: {missing}")
    if unexpected:
        raise ValueError(f"{object_name} has unexpected fields: {unexpected}")


@dataclass(frozen=True)
class CoherentSourceCommitments:
    """The nine commitments shared by every coherent source component."""

    source_root_hash: str
    branch_hash: str
    field_census_hash: str
    scheme_hash: str
    threshold_hash: str
    fj_convention_hash: str
    perturbative_mask_hash: str
    analytic_sheet_hash: str
    units_clock_hash: str

    def __post_init__(self) -> None:
        for name in COMMITMENT_FIELDS:
            _require_hash(getattr(self, name), name)

    @classmethod
    def from_artifact(cls, payload: Mapping[str, Any]) -> "CoherentSourceCommitments":
        if not isinstance(payload, Mapping):
            raise TypeError("commitments must be a mapping")
        _strict_keys(
            payload,
            required=set(COMMITMENT_FIELDS),
            object_name="coherent source commitments",
        )
        return cls(**{name: payload[name] for name in COMMITMENT_FIELDS})

    def to_artifact(self) -> dict[str, str]:
        return {name: getattr(self, name) for name in COMMITMENT_FIELDS}

    @property
    def bundle_hash(self) -> str:
        return _canonical_hash(self.to_artifact())


@dataclass(frozen=True)
class SourceComponentEnvelope:
    """A component identity bound to the complete coherent commitment tuple."""

    component_kind: str
    artifact_hash: str
    commitments: CoherentSourceCommitments

    def __post_init__(self) -> None:
        if self.component_kind not in REQUIRED_SOURCE_COMPONENTS:
            raise ValueError(f"unknown coherent source component: {self.component_kind!r}")
        _require_hash(self.artifact_hash, "artifact_hash")


def verify_coherent_source_components(
    components: Sequence[SourceComponentEnvelope],
) -> dict[str, Any]:
    """Verify component census and equality of all nine commitments.

    Equality is structural only.  Hash strings are identifiers, not evidence
    that any component was physically produced.
    """

    by_kind: dict[str, SourceComponentEnvelope] = {}
    duplicates: list[str] = []
    for component in components:
        if not isinstance(component, SourceComponentEnvelope):
            raise TypeError("components must contain SourceComponentEnvelope values")
        if component.component_kind in by_kind:
            duplicates.append(component.component_kind)
        else:
            by_kind[component.component_kind] = component

    missing = sorted(set(REQUIRED_SOURCE_COMPONENTS) - set(by_kind))
    mismatches: dict[str, list[str]] = {}
    baseline = by_kind.get("oph_source")
    if baseline is None and by_kind:
        baseline = next(iter(by_kind.values()))
    if baseline is not None:
        for kind, component in by_kind.items():
            fields = [
                name
                for name in COMMITMENT_FIELDS
                if getattr(component.commitments, name)
                != getattr(baseline.commitments, name)
            ]
            if fields:
                mismatches[kind] = fields

    valid = not missing and not duplicates and not mismatches
    blockers = [*(f"missing_component:{item}" for item in missing)]
    blockers.extend(f"duplicate_component:{item}" for item in sorted(set(duplicates)))
    blockers.extend(
        f"commitment_mismatch:{kind}:{field_name}"
        for kind, fields in sorted(mismatches.items())
        for field_name in fields
    )
    return {
        "schema_version": SOURCE_TO_EFT_SCHEMA_VERSION,
        "required_components": list(REQUIRED_SOURCE_COMPONENTS),
        "present_components": sorted(by_kind),
        "missing_components": missing,
        "duplicate_components": sorted(set(duplicates)),
        "commitment_mismatches": mismatches,
        "commitment_bundle_hash": baseline.commitments.bundle_hash if baseline else None,
        "receipts": {
            "COHERENT_SOURCE_COMPONENT_CENSUS_RECEIPT": not missing and not duplicates,
            "COHERENT_SOURCE_COMMITMENT_EQUALITY_RECEIPT": valid,
            "COHERENT_PHYSICAL_SOURCE_PACKET_RECEIPT": False,
        },
        "blockers": blockers
        + ["primitive_component_producers_and_replay_verifiers_missing"],
    }


@dataclass(frozen=True)
class PerturbativeMask:
    """A finite, exact, downward-closed set of coupling/loop multiindices."""

    axes: tuple[str, ...]
    monomials: tuple[tuple[int, ...], ...]

    def __post_init__(self) -> None:
        axes = tuple(self.axes)
        monomials = tuple(tuple(row) for row in self.monomials)
        if not axes or len(axes) > 16 or len(set(axes)) != len(axes):
            raise ValueError("mask axes must contain 1..16 unique identifiers")
        for index, axis in enumerate(axes):
            _require_identifier(axis, f"axes[{index}]")
        if not monomials or len(monomials) > 10_000:
            raise ValueError("mask must contain 1..10000 monomials")
        for row_index, row in enumerate(monomials):
            if len(row) != len(axes):
                raise ValueError(f"monomials[{row_index}] has the wrong dimension")
            for exponent in row:
                if isinstance(exponent, bool) or not isinstance(exponent, int):
                    raise TypeError("mask exponents must be integers")
                if exponent < 0 or exponent > 16:
                    raise ValueError("mask exponents must lie in [0, 16]")
        if len(set(monomials)) != len(monomials):
            raise ValueError("mask monomials must be unique")
        present = set(monomials)
        for row in monomials:
            for lower in itertools.product(*(range(value + 1) for value in row)):
                if lower not in present:
                    raise ValueError(
                        f"perturbative mask is not downward closed: {lower} missing below {row}"
                    )
        object.__setattr__(self, "axes", axes)
        object.__setattr__(self, "monomials", tuple(sorted(monomials)))

    @property
    def commitment(self) -> str:
        return _canonical_hash(
            {"axes": list(self.axes), "monomials": [list(row) for row in self.monomials]}
        )

    def to_artifact(self) -> dict[str, Any]:
        return {
            "axes": list(self.axes),
            "monomials": [list(row) for row in self.monomials],
            "perturbative_mask_hash": self.commitment,
        }


def verify_common_perturbative_mask(
    direct: PerturbativeMask,
    converted: PerturbativeMask,
) -> dict[str, Any]:
    equal = direct == converted and direct.commitment == converted.commitment
    return {
        "direct_mask_hash": direct.commitment,
        "converted_mask_hash": converted.commitment,
        "receipts": {
            "DIRECT_MASK_DOWNWARD_CLOSED_RECEIPT": True,
            "CONVERTED_MASK_DOWNWARD_CLOSED_RECEIPT": True,
            "COMMON_PERTURBATIVE_MASK_EQUALITY_RECEIPT": equal,
            "FJ1_RECEIPT": False,
            "BRST1_RECEIPT": False,
        },
        "blockers": [] if equal else ["perturbative_monomial_mask_mismatch"],
        "physical_boundary": (
            "Mask equality is necessary but does not produce a source Higgs action, "
            "an FJ conversion, or a BRST-complete kernel."
        ),
    }


@dataclass(frozen=True)
class QuotientMomentPacket:
    """Exact quotient-visible observables, with no EFT parameter slot."""

    regulator_id: str
    quotient_artifact_hash: str
    observable_map_hash: str
    observable_names: tuple[str, ...]
    moments: tuple[Fraction, ...]
    semantics: MomentSemantics
    admissible_branch_hashes: tuple[str, ...]
    commitments: CoherentSourceCommitments
    settled_state_hash: str | None = None
    source_law_hash: str | None = None
    refinement_defects: tuple[Fraction, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        _require_identifier(self.regulator_id, "regulator_id")
        _require_hash(self.quotient_artifact_hash, "quotient_artifact_hash")
        _require_hash(self.observable_map_hash, "observable_map_hash")
        names = tuple(self.observable_names)
        if not names or len(names) > 256 or len(set(names)) != len(names):
            raise ValueError("observable_names must contain 1..256 unique identifiers")
        for index, name in enumerate(names):
            _require_identifier(name, f"observable_names[{index}]")
        moments = _fraction_tuple(self.moments, "moments")
        if len(moments) != len(names):
            raise ValueError("moments must have one value per observable")
        if not isinstance(self.semantics, MomentSemantics):
            raise TypeError("semantics must be a MomentSemantics value")
        branches = tuple(self.admissible_branch_hashes)
        if not branches or len(branches) > 4096 or len(set(branches)) != len(branches):
            raise ValueError("admissible_branch_hashes must be nonempty and unique")
        for index, branch_hash in enumerate(branches):
            _require_hash(branch_hash, f"admissible_branch_hashes[{index}]")
        if self.commitments.branch_hash not in branches:
            raise ValueError("the committed branch must occur in admissible_branch_hashes")
        if self.semantics is MomentSemantics.DETERMINISTIC_SETTLED_STATE:
            _require_hash(self.settled_state_hash, "settled_state_hash")
            if self.source_law_hash is not None:
                raise ValueError("deterministic moments cannot carry a source_law_hash")
        else:
            _require_hash(self.source_law_hash, "source_law_hash")
            if self.settled_state_hash is not None:
                raise ValueError("stochastic moments cannot carry a settled_state_hash")
        defects = _fraction_tuple(self.refinement_defects, "refinement_defects")
        if any(value < 0 for value in defects):
            raise ValueError("refinement_defects must be nonnegative")
        object.__setattr__(self, "observable_names", names)
        object.__setattr__(self, "moments", moments)
        object.__setattr__(self, "admissible_branch_hashes", branches)
        object.__setattr__(self, "refinement_defects", defects)

    @classmethod
    def from_artifact(cls, payload: Mapping[str, Any]) -> "QuotientMomentPacket":
        """Parse a strict artifact; an injected ``theta`` field is rejected."""

        if not isinstance(payload, Mapping):
            raise TypeError("quotient moment packet must be a mapping")
        required = {
            "schema_version",
            "regulator_id",
            "quotient_artifact_hash",
            "observable_map_hash",
            "observable_names",
            "moments",
            "semantics",
            "admissible_branch_hashes",
            "commitments",
            "refinement_defects",
        }
        optional = {"settled_state_hash", "source_law_hash"}
        _strict_keys(
            payload,
            required=required,
            optional=optional,
            object_name="quotient moment packet",
        )
        if payload["schema_version"] != QUOTIENT_MOMENT_SCHEMA_VERSION:
            raise ValueError("unsupported quotient moment schema_version")
        try:
            semantics = MomentSemantics(payload["semantics"])
        except (TypeError, ValueError) as exc:
            raise ValueError("unknown quotient moment semantics") from exc
        return cls(
            regulator_id=payload["regulator_id"],
            quotient_artifact_hash=payload["quotient_artifact_hash"],
            observable_map_hash=payload["observable_map_hash"],
            observable_names=tuple(payload["observable_names"]),
            moments=_fraction_tuple(payload["moments"], "moments"),
            semantics=semantics,
            admissible_branch_hashes=tuple(payload["admissible_branch_hashes"]),
            commitments=CoherentSourceCommitments.from_artifact(payload["commitments"]),
            settled_state_hash=payload.get("settled_state_hash"),
            source_law_hash=payload.get("source_law_hash"),
            refinement_defects=_fraction_tuple(
                payload["refinement_defects"], "refinement_defects"
            ),
        )

    def to_artifact(self) -> dict[str, Any]:
        artifact: dict[str, Any] = {
            "schema_version": QUOTIENT_MOMENT_SCHEMA_VERSION,
            "regulator_id": self.regulator_id,
            "quotient_artifact_hash": self.quotient_artifact_hash,
            "observable_map_hash": self.observable_map_hash,
            "observable_names": list(self.observable_names),
            "moments": [_qstr(value) for value in self.moments],
            "semantics": self.semantics.value,
            "admissible_branch_hashes": list(self.admissible_branch_hashes),
            "commitments": self.commitments.to_artifact(),
            "refinement_defects": [_qstr(value) for value in self.refinement_defects],
        }
        if self.settled_state_hash is not None:
            artifact["settled_state_hash"] = self.settled_state_hash
        if self.source_law_hash is not None:
            artifact["source_law_hash"] = self.source_law_hash
        return artifact


@dataclass(frozen=True)
class ExactAffineMomentMap:
    """A frozen exact affine map ``F(theta) = A theta + b``."""

    parameter_names: tuple[str, ...]
    observable_names: tuple[str, ...]
    matrix: tuple[tuple[Fraction, ...], ...]
    offset: tuple[Fraction, ...]
    map_artifact_hash: str

    def __post_init__(self) -> None:
        parameters = tuple(self.parameter_names)
        observables = tuple(self.observable_names)
        if not parameters or len(parameters) > 256 or len(set(parameters)) != len(parameters):
            raise ValueError("parameter_names must contain 1..256 unique identifiers")
        if not observables or len(observables) > 256 or len(set(observables)) != len(observables):
            raise ValueError("observable_names must contain 1..256 unique identifiers")
        for index, name in enumerate(parameters):
            _require_identifier(name, f"parameter_names[{index}]")
        for index, name in enumerate(observables):
            _require_identifier(name, f"observable_names[{index}]")
        matrix = _matrix(self.matrix, len(observables), len(parameters), "matrix")
        offset = _fraction_tuple(self.offset, "offset")
        if len(offset) != len(observables):
            raise ValueError("offset has the wrong dimension")
        _require_hash(self.map_artifact_hash, "map_artifact_hash")
        object.__setattr__(self, "parameter_names", parameters)
        object.__setattr__(self, "observable_names", observables)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "offset", offset)

    def evaluate(self, theta: Sequence[Any]) -> tuple[Fraction, ...]:
        vector = _fraction_tuple(theta, "theta")
        if len(vector) != len(self.parameter_names):
            raise ValueError("theta has the wrong dimension")
        return tuple(
            self.offset[row_index]
            + sum(coefficient * value for coefficient, value in zip(row, vector, strict=True))
            for row_index, row in enumerate(self.matrix)
        )


@dataclass(frozen=True)
class UntrustedParameterCandidate:
    """A candidate used to check a theorem; never evidence of source production."""

    values: tuple[Fraction, ...]
    candidate_artifact_hash: str

    def __post_init__(self) -> None:
        values = _fraction_tuple(self.values, "candidate values")
        if not values:
            raise ValueError("candidate values cannot be empty")
        _require_hash(self.candidate_artifact_hash, "candidate_artifact_hash")
        object.__setattr__(self, "values", values)


@dataclass(frozen=True)
class ExactLeftInverseCertificate:
    """Exact rational matrix claimed to be a left inverse of the moment map."""

    matrix: tuple[tuple[Fraction, ...], ...]
    certificate_artifact_hash: str

    def __post_init__(self) -> None:
        rows = tuple(tuple(row) for row in self.matrix)
        if not rows or not rows[0]:
            raise ValueError("left inverse matrix cannot be empty")
        width = len(rows[0])
        if any(len(row) != width for row in rows):
            raise ValueError("left inverse matrix must be rectangular")
        parsed = _matrix(rows, len(rows), width, "left_inverse")
        _require_hash(self.certificate_artifact_hash, "certificate_artifact_hash")
        object.__setattr__(self, "matrix", parsed)


@dataclass(frozen=True)
class ExactIdentifiabilityCertificate:
    """Replay result for the affine source-parameter identifiability theorem."""

    left_inverse_exact: bool
    lower_inverse_modulus: Fraction
    actual_residual: Fraction
    declared_residual_bound: Fraction
    ambiguity_diameter_bound: Fraction | None
    observable_binding_valid: bool
    branch_unique: bool
    blockers: tuple[str, ...]

    @property
    def mathematical_receipt(self) -> bool:
        return bool(
            self.left_inverse_exact
            and self.observable_binding_valid
            and self.lower_inverse_modulus > 0
            and self.actual_residual <= self.declared_residual_bound
        )

    def to_report(self) -> dict[str, Any]:
        return {
            "norm": "linf",
            "left_inverse_exact": self.left_inverse_exact,
            "lower_inverse_modulus_sigma": _qstr(self.lower_inverse_modulus),
            "actual_residual": _qstr(self.actual_residual),
            "declared_residual_bound": _qstr(self.declared_residual_bound),
            "two_candidate_ambiguity_diameter_bound": (
                _qstr(self.ambiguity_diameter_bound)
                if self.ambiguity_diameter_bound is not None
                else None
            ),
            "observable_binding_valid": self.observable_binding_valid,
            "branch_unique": self.branch_unique,
            "receipts": {
                "EXACT_AFFINE_LEFT_INVERSE_RECEIPT": self.left_inverse_exact,
                "EXACT_SOURCE_PARAMETER_IDENTIFIABILITY_MATH_RECEIPT": (
                    self.mathematical_receipt
                ),
                "SOURCE_PARAMETER_CANDIDATE_PRODUCER_RECEIPT": False,
                "OPH_NATIVE_EFT_REALIZATION_RECEIPT": False,
            },
            "blockers": list(self.blockers),
            "physical_boundary": (
                "The candidate is intentionally untrusted. Exact residual and inverse-modulus "
                "checks prove an implication about this map; they do not prove that the "
                "candidate or map came from the OPH source."
            ),
        }


def verify_exact_affine_identifiability(
    packet: QuotientMomentPacket,
    moment_map: ExactAffineMomentMap,
    candidate: UntrustedParameterCandidate,
    left_inverse: ExactLeftInverseCertificate,
    *,
    residual_bound: Any,
) -> ExactIdentifiabilityCertificate:
    """Replay the exact affine identifiability theorem in the infinity norm.

    If ``B A = I``, then ``||A delta||_inf >= sigma ||delta||_inf`` with
    ``sigma = 1 / ||B||_inf``.  All arithmetic is rational and exact.
    """

    epsilon = _fraction(residual_bound, "residual_bound")
    if epsilon < 0:
        raise ValueError("residual_bound must be nonnegative")
    if len(candidate.values) != len(moment_map.parameter_names):
        raise ValueError("candidate dimension does not match the moment map")
    if len(left_inverse.matrix) != len(moment_map.parameter_names):
        raise ValueError("left inverse row count must equal the parameter dimension")
    if any(len(row) != len(moment_map.observable_names) for row in left_inverse.matrix):
        raise ValueError("left inverse column count must equal the observable dimension")

    product = _matmul(left_inverse.matrix, moment_map.matrix)
    identity = _identity(len(moment_map.parameter_names))
    exact = product == identity
    operator_norm = max(sum(abs(value) for value in row) for row in left_inverse.matrix)
    sigma = Fraction(0) if operator_norm == 0 else Fraction(1, 1) / operator_norm
    evaluated = moment_map.evaluate(candidate.values)
    if len(packet.moments) != len(evaluated):
        raise ValueError("moment packet dimension does not match the moment map")
    actual_residual = max(
        (abs(actual - expected) for actual, expected in zip(evaluated, packet.moments, strict=True)),
        default=Fraction(0),
    )
    observable_binding = bool(
        packet.observable_names == moment_map.observable_names
        and packet.observable_map_hash == moment_map.map_artifact_hash
    )
    blockers: list[str] = []
    if not exact:
        blockers.append("left_inverse_identity_failed")
    if sigma <= 0:
        blockers.append("positive_lower_inverse_modulus_missing")
    if not observable_binding:
        blockers.append("quotient_observable_map_binding_mismatch")
    if actual_residual > epsilon:
        blockers.append("candidate_residual_exceeds_declared_bound")
    branch_unique = len(packet.admissible_branch_hashes) == 1
    if not branch_unique:
        blockers.append("source_branch_not_unique")
    ambiguity_bound = Fraction(2) * epsilon / sigma if exact and sigma > 0 else None
    blockers.extend(
        (
            "source_parameter_candidate_producer_verifier_missing",
            "frozen_moment_map_physical_realization_verifier_missing",
        )
    )
    return ExactIdentifiabilityCertificate(
        left_inverse_exact=exact,
        lower_inverse_modulus=sigma,
        actual_residual=actual_residual,
        declared_residual_bound=epsilon,
        ambiguity_diameter_bound=ambiguity_bound,
        observable_binding_valid=observable_binding,
        branch_unique=branch_unique,
        blockers=tuple(blockers),
    )


def _matrix(
    values: Sequence[Sequence[Any]],
    rows: int,
    columns: int,
    field_name: str,
) -> tuple[tuple[Fraction, ...], ...]:
    if isinstance(values, (str, bytes)) or len(values) != rows:
        raise ValueError(f"{field_name} must have {rows} rows")
    result: list[tuple[Fraction, ...]] = []
    for row_index, row in enumerate(values):
        parsed = _fraction_tuple(row, f"{field_name}[{row_index}]")
        if len(parsed) != columns:
            raise ValueError(f"{field_name}[{row_index}] must have {columns} columns")
        result.append(parsed)
    return tuple(result)


def _identity(size: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(int(row == column)) for column in range(size))
        for row in range(size)
    )


def _matmul(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    if not left or not right or not right[0]:
        raise ValueError("cannot multiply empty matrices")
    inner = len(left[0])
    if any(len(row) != inner for row in left):
        raise ValueError("left matrix is ragged")
    if len(right) != inner:
        raise ValueError("matrix dimensions do not compose")
    width = len(right[0])
    if any(len(row) != width for row in right):
        raise ValueError("right matrix is ragged")
    return tuple(
        tuple(
            sum(
                (left[row][index] * right[index][column] for index in range(inner)),
                Fraction(0),
            )
            for column in range(width)
        )
        for row in range(len(left))
    )


def _matvec(
    matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    if not matrix or any(len(row) != len(vector) for row in matrix):
        raise ValueError("matrix and vector dimensions do not compose")
    return tuple(
        sum((coefficient * value for coefficient, value in zip(row, vector, strict=True)), Fraction(0))
        for row in matrix
    )


@dataclass(frozen=True)
class FiniteSourceLawPacket:
    """An exact finite source-law representation with explicit semantics.

    A deterministic source point and a stochastic finite ensemble are
    different types of evidence even though both have computable moments.
    The packet never accepts a caller-authored covariance matrix.
    """

    semantics: SourceLawSemantics
    coordinate_names: tuple[str, ...]
    admissible_branch_hashes: tuple[str, ...]
    source_law_hash: str
    commitments: CoherentSourceCommitments
    deterministic_point: tuple[Fraction, ...] = field(default_factory=tuple)
    support_points: tuple[tuple[Fraction, ...], ...] = field(default_factory=tuple)
    weights: tuple[Fraction, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not isinstance(self.semantics, SourceLawSemantics):
            raise TypeError("semantics must be a SourceLawSemantics value")
        names = tuple(self.coordinate_names)
        if not names or len(names) > 256 or len(set(names)) != len(names):
            raise ValueError("coordinate_names must contain 1..256 unique identifiers")
        for index, name in enumerate(names):
            _require_identifier(name, f"coordinate_names[{index}]")
        branches = tuple(self.admissible_branch_hashes)
        if not branches or len(branches) > 4096 or len(set(branches)) != len(branches):
            raise ValueError("admissible_branch_hashes must be nonempty and unique")
        for index, branch_hash in enumerate(branches):
            _require_hash(branch_hash, f"admissible_branch_hashes[{index}]")
        if self.commitments.branch_hash not in branches:
            raise ValueError("the committed branch must occur in admissible_branch_hashes")
        _require_hash(self.source_law_hash, "source_law_hash")

        point = _fraction_tuple(self.deterministic_point, "deterministic_point")
        support = tuple(
            _fraction_tuple(row, f"support_points[{index}]")
            for index, row in enumerate(self.support_points)
        )
        weights = _fraction_tuple(self.weights, "weights")
        if self.semantics is SourceLawSemantics.DETERMINISTIC:
            if len(point) != len(names):
                raise ValueError("deterministic_point has the wrong dimension")
            if support or weights:
                raise ValueError("a deterministic source cannot carry ensemble points or weights")
        else:
            if point:
                raise ValueError("a stochastic source cannot carry a deterministic point")
            if not support or len(support) != len(weights):
                raise ValueError("a stochastic source requires one weight per support point")
            if len(support) > 100_000:
                raise ValueError("stochastic support is too large for this finite contract")
            if any(len(row) != len(names) for row in support):
                raise ValueError("every stochastic support point must have the declared dimension")
            if any(weight < 0 for weight in weights):
                raise ValueError("stochastic weights must be nonnegative")
            if sum(weights, Fraction(0)) != 1:
                raise ValueError("stochastic weights must sum exactly to one")

        object.__setattr__(self, "coordinate_names", names)
        object.__setattr__(self, "admissible_branch_hashes", branches)
        object.__setattr__(self, "deterministic_point", point)
        object.__setattr__(self, "support_points", support)
        object.__setattr__(self, "weights", weights)

    @classmethod
    def from_artifact(cls, payload: Mapping[str, Any]) -> "FiniteSourceLawPacket":
        if not isinstance(payload, Mapping):
            raise TypeError("source law packet must be a mapping")
        common = {
            "schema_version",
            "semantics",
            "coordinate_names",
            "admissible_branch_hashes",
            "source_law_hash",
            "commitments",
        }
        try:
            semantics = SourceLawSemantics(payload.get("semantics"))
        except (TypeError, ValueError) as exc:
            raise ValueError("unknown source-law semantics") from exc
        if semantics is SourceLawSemantics.DETERMINISTIC:
            required = common | {"deterministic_point"}
        else:
            required = common | {"support_points", "weights"}
        _strict_keys(payload, required=required, object_name="finite source law packet")
        if payload["schema_version"] != SOURCE_LAW_SCHEMA_VERSION:
            raise ValueError("unsupported source law schema_version")
        return cls(
            semantics=semantics,
            coordinate_names=tuple(payload["coordinate_names"]),
            admissible_branch_hashes=tuple(payload["admissible_branch_hashes"]),
            source_law_hash=payload["source_law_hash"],
            commitments=CoherentSourceCommitments.from_artifact(payload["commitments"]),
            deterministic_point=_fraction_tuple(
                payload.get("deterministic_point", ()), "deterministic_point"
            ),
            support_points=tuple(
                _fraction_tuple(row, f"support_points[{index}]")
                for index, row in enumerate(payload.get("support_points", ()))
            ),
            weights=_fraction_tuple(payload.get("weights", ()), "weights"),
        )

    def mean(self) -> tuple[Fraction, ...]:
        if self.semantics is SourceLawSemantics.DETERMINISTIC:
            return self.deterministic_point
        return tuple(
            sum(
                (weight * point[column] for point, weight in zip(self.support_points, self.weights, strict=True)),
                Fraction(0),
            )
            for column in range(len(self.coordinate_names))
        )

    def covariance(self) -> tuple[tuple[Fraction, ...], ...] | None:
        """Return exact source covariance only when source semantics license it.

        A branch-ambiguous deterministic packet has no selected law, so it has
        no covariance.  A unique deterministic source has exactly zero source
        covariance; deterministic enclosures belong in ``DeterministicErrorBudget``.
        """

        dimension = len(self.coordinate_names)
        if self.semantics is SourceLawSemantics.DETERMINISTIC:
            if len(self.admissible_branch_hashes) != 1:
                return None
            return tuple(tuple(Fraction(0) for _ in range(dimension)) for _ in range(dimension))
        mean = self.mean()
        return tuple(
            tuple(
                sum(
                    (
                        weight
                        * (point[row] - mean[row])
                        * (point[column] - mean[column])
                        for point, weight in zip(self.support_points, self.weights, strict=True)
                    ),
                    Fraction(0),
                )
                for column in range(dimension)
            )
            for row in range(dimension)
        )


@dataclass(frozen=True)
class DeterministicErrorBudget:
    """Non-probabilistic errors that must never be relabelled as covariance."""

    root_enclosure: Fraction = Fraction(0)
    matching_remainder: Fraction = Fraction(0)
    pole_truncation: Fraction = Fraction(0)
    scale_variation: Fraction = Fraction(0)
    rounding_error: Fraction = Fraction(0)
    refinement_error: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        for name in (
            "root_enclosure",
            "matching_remainder",
            "pole_truncation",
            "scale_variation",
            "rounding_error",
            "refinement_error",
        ):
            parsed = _fraction(getattr(self, name), name)
            if parsed < 0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, parsed)

    def to_artifact(self) -> dict[str, str]:
        return {
            name: _qstr(getattr(self, name))
            for name in (
                "root_enclosure",
                "matching_remainder",
                "pole_truncation",
                "scale_variation",
                "rounding_error",
                "refinement_error",
            )
        }


def finite_source_law_report(packet: FiniteSourceLawPacket) -> dict[str, Any]:
    branch_unique = len(packet.admissible_branch_hashes) == 1
    deterministic_valid = bool(
        packet.semantics is SourceLawSemantics.DETERMINISTIC and branch_unique
    )
    stochastic_valid = packet.semantics is SourceLawSemantics.STOCHASTIC
    covariance = packet.covariance()
    blockers: list[str] = []
    if packet.semantics is SourceLawSemantics.DETERMINISTIC and not branch_unique:
        blockers.append("deterministic_source_forbidden_by_branch_ambiguity")
    blockers.append("source_law_primitive_producer_or_theorem_verifier_missing")
    return {
        "schema_version": SOURCE_TO_EFT_SCHEMA_VERSION,
        "semantics": packet.semantics.value,
        "branch_count": len(packet.admissible_branch_hashes),
        "source_mean": [_qstr(value) for value in packet.mean()],
        "source_covariance": (
            [[_qstr(value) for value in row] for row in covariance]
            if covariance is not None
            else None
        ),
        "deterministic_error_channels": [
            "root_enclosure",
            "matching_remainder",
            "pole_truncation",
            "scale_variation",
            "rounding_error",
            "refinement_error",
            "discrete_branch_set",
        ],
        "receipts": {
            "FINITE_SOURCE_LAW_STRUCTURAL_RECEIPT": deterministic_valid or stochastic_valid,
            "DETERMINISTIC_SOURCE_SEMANTICS_RECEIPT": deterministic_valid,
            "STOCHASTIC_SOURCE_SEMANTICS_RECEIPT": stochastic_valid,
            "SOURCE_LAW_PRIMITIVE_PRODUCER_RECEIPT": False,
            "COV1_RECEIPT": False,
        },
        "blockers": blockers,
    }


def verify_finite_source_law_artifact(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Fail-closed artifact verifier, including interval-as-covariance mutations."""

    try:
        packet = FiniteSourceLawPacket.from_artifact(payload)
    except (KeyError, TypeError, ValueError) as exc:
        return {
            "schema_version": SOURCE_TO_EFT_SCHEMA_VERSION,
            "receipts": {
                "FINITE_SOURCE_LAW_STRUCTURAL_RECEIPT": False,
                "SOURCE_LAW_PRIMITIVE_PRODUCER_RECEIPT": False,
                "COV1_RECEIPT": False,
            },
            "blockers": [f"invalid_source_law_artifact:{exc}"],
        }
    return finite_source_law_report(packet)


@dataclass(frozen=True)
class EFTIntervalCensus:
    """Complete field census and scheme for one open EFT scale interval."""

    interval_id: str
    field_names: tuple[str, ...]
    renormalization_scheme: str

    def __post_init__(self) -> None:
        _require_identifier(self.interval_id, "interval_id")
        fields = tuple(self.field_names)
        if not fields or len(fields) > 4096 or len(set(fields)) != len(fields):
            raise ValueError("field_names must contain a complete nonempty unique census")
        for index, name in enumerate(fields):
            _require_identifier(name, f"field_names[{index}]")
        _require_identifier(self.renormalization_scheme, "renormalization_scheme")
        object.__setattr__(self, "field_names", fields)

    @property
    def census_hash(self) -> str:
        return _canonical_hash(
            {"interval_id": self.interval_id, "field_names": list(self.field_names)}
        )


@dataclass(frozen=True)
class ExactAffineTransformation:
    """An exact rational affine transformation with an exact L-infinity norm."""

    input_coordinates: tuple[str, ...]
    output_coordinates: tuple[str, ...]
    matrix: tuple[tuple[Fraction, ...], ...]
    offset: tuple[Fraction, ...]
    artifact_hash: str

    def __post_init__(self) -> None:
        inputs = tuple(self.input_coordinates)
        outputs = tuple(self.output_coordinates)
        if not inputs or not outputs:
            raise ValueError("affine transformations require nonempty coordinate sets")
        if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
            raise ValueError("affine transformation coordinates must be unique")
        for index, name in enumerate(inputs):
            _require_identifier(name, f"input_coordinates[{index}]")
        for index, name in enumerate(outputs):
            _require_identifier(name, f"output_coordinates[{index}]")
        matrix = _matrix(self.matrix, len(outputs), len(inputs), "matrix")
        offset = _fraction_tuple(self.offset, "offset")
        if len(offset) != len(outputs):
            raise ValueError("offset has the wrong dimension")
        _require_hash(self.artifact_hash, "artifact_hash")
        object.__setattr__(self, "input_coordinates", inputs)
        object.__setattr__(self, "output_coordinates", outputs)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "offset", offset)

    @property
    def linf_lipschitz(self) -> Fraction:
        return max(sum(abs(value) for value in row) for row in self.matrix)

    def evaluate(self, value: Sequence[Any]) -> tuple[Fraction, ...]:
        vector = _fraction_tuple(value, "affine input")
        if len(vector) != len(self.input_coordinates):
            raise ValueError("affine input has the wrong dimension")
        linear = _matvec(self.matrix, vector)
        return tuple(a + b for a, b in zip(linear, self.offset, strict=True))


@dataclass(frozen=True)
class ExactMatchingStep:
    """One RG-flow then finite matching/decoupling step."""

    step_id: str
    source_census: EFTIntervalCensus
    target_census: EFTIntervalCensus
    flow: ExactAffineTransformation
    matching: ExactAffineTransformation
    implementation_remainder_bound: Fraction
    perturbative_mask: PerturbativeMask

    def __post_init__(self) -> None:
        _require_identifier(self.step_id, "step_id")
        if not isinstance(self.source_census, EFTIntervalCensus):
            raise TypeError("source_census must be an EFTIntervalCensus")
        if not isinstance(self.target_census, EFTIntervalCensus):
            raise TypeError("target_census must be an EFTIntervalCensus")
        if not isinstance(self.flow, ExactAffineTransformation):
            raise TypeError("flow must be an ExactAffineTransformation")
        if not isinstance(self.matching, ExactAffineTransformation):
            raise TypeError("matching must be an ExactAffineTransformation")
        if self.flow.output_coordinates != self.matching.input_coordinates:
            raise ValueError("flow output coordinates do not feed the matching map")
        remainder = _fraction(
            self.implementation_remainder_bound, "implementation_remainder_bound"
        )
        if remainder < 0:
            raise ValueError("implementation_remainder_bound must be nonnegative")
        if not isinstance(self.perturbative_mask, PerturbativeMask):
            raise TypeError("perturbative_mask must be a PerturbativeMask")
        object.__setattr__(self, "implementation_remainder_bound", remainder)

    @property
    def combined_matrix(self) -> tuple[tuple[Fraction, ...], ...]:
        return _matmul(self.matching.matrix, self.flow.matrix)

    @property
    def combined_offset(self) -> tuple[Fraction, ...]:
        transformed = _matvec(self.matching.matrix, self.flow.offset)
        return tuple(
            left + right for left, right in zip(transformed, self.matching.offset, strict=True)
        )

    @property
    def lipschitz_bound(self) -> Fraction:
        return self.matching.linf_lipschitz * self.flow.linf_lipschitz


def matching_field_census_hash(steps: Sequence[ExactMatchingStep]) -> str:
    if not steps:
        raise ValueError("matching tower cannot be empty")
    intervals = [steps[0].source_census, *(step.target_census for step in steps)]
    return _canonical_hash(
        [
            {
                "interval_id": interval.interval_id,
                "field_names": list(interval.field_names),
            }
            for interval in intervals
        ]
    )


def matching_scheme_hash(steps: Sequence[ExactMatchingStep]) -> str:
    if not steps:
        raise ValueError("matching tower cannot be empty")
    intervals = [steps[0].source_census, *(step.target_census for step in steps)]
    return _canonical_hash(
        [
            {
                "interval_id": interval.interval_id,
                "renormalization_scheme": interval.renormalization_scheme,
            }
            for interval in intervals
        ]
    )


def matching_threshold_hash(steps: Sequence[ExactMatchingStep]) -> str:
    if not steps:
        raise ValueError("matching tower cannot be empty")
    return _canonical_hash(
        [
            {
                "step_id": step.step_id,
                "source_interval": step.source_census.interval_id,
                "target_interval": step.target_census.interval_id,
                "flow_artifact_hash": step.flow.artifact_hash,
                "matching_artifact_hash": step.matching.artifact_hash,
            }
            for step in steps
        ]
    )


@dataclass(frozen=True)
class ExactMatchingTower:
    """A field-censused exact affine RG/matching tower."""

    steps: tuple[ExactMatchingStep, ...]
    commitments: CoherentSourceCommitments

    def __post_init__(self) -> None:
        steps = tuple(self.steps)
        if not steps or len(steps) > 1024:
            raise ValueError("matching tower must contain 1..1024 steps")
        for index, step in enumerate(steps):
            if not isinstance(step, ExactMatchingStep):
                raise TypeError(f"steps[{index}] must be an ExactMatchingStep")
        for previous, current in zip(steps, steps[1:]):
            if previous.target_census != current.source_census:
                raise ValueError("matching tower field-census chain is discontinuous")
            if previous.matching.output_coordinates != current.flow.input_coordinates:
                raise ValueError("matching tower coordinate chain is discontinuous")
        masks = {step.perturbative_mask.commitment for step in steps}
        if len(masks) != 1:
            raise ValueError("matching tower steps use different perturbative masks")
        mask_hash = next(iter(masks))
        expected = {
            "field_census_hash": matching_field_census_hash(steps),
            "scheme_hash": matching_scheme_hash(steps),
            "threshold_hash": matching_threshold_hash(steps),
            "perturbative_mask_hash": mask_hash,
        }
        for field_name, value in expected.items():
            if getattr(self.commitments, field_name) != value:
                raise ValueError(f"matching tower {field_name} commitment mismatch")
        object.__setattr__(self, "steps", steps)


@dataclass(frozen=True)
class ExactMatchingCompositionCertificate:
    endpoint_error_bound: Fraction
    endpoint_jacobian: tuple[tuple[Fraction, ...], ...]
    total_matrix: tuple[tuple[Fraction, ...], ...]
    total_offset: tuple[Fraction, ...]
    step_error_bounds: tuple[Fraction, ...]
    step_lipschitz_bounds: tuple[Fraction, ...]

    def to_report(self) -> dict[str, Any]:
        return {
            "arithmetic": "exact_rational",
            "norm": "linf",
            "endpoint_error_bound": _qstr(self.endpoint_error_bound),
            "step_error_bounds": [_qstr(value) for value in self.step_error_bounds],
            "step_lipschitz_bounds": [
                _qstr(value) for value in self.step_lipschitz_bounds
            ],
            "endpoint_jacobian": [
                [_qstr(value) for value in row] for row in self.endpoint_jacobian
            ],
            "total_affine_map": {
                "matrix": [[_qstr(value) for value in row] for row in self.total_matrix],
                "offset": [_qstr(value) for value in self.total_offset],
            },
            "receipts": {
                "MATCHING_TOWER_FIELD_CENSUS_BINDING_RECEIPT": True,
                "MATCHING_TOWER_COMMON_MASK_RECEIPT": True,
                "MATCHING_TOWER_EXACT_ERROR_COMPOSITION_RECEIPT": True,
                "MATCHING_TOWER_EXACT_JACOBIAN_COMPOSITION_RECEIPT": True,
                "MATCH1_RECEIPT": False,
            },
            "blockers": [
                "matching_maps_have_no_independent_primitive_producer_or_replay_verifier"
            ],
        }


def compose_exact_matching_tower(
    tower: ExactMatchingTower, *, initial_error_bound: Any
) -> ExactMatchingCompositionCertificate:
    """Compose the exact Jacobian and deterministic error recurrence."""

    error = _fraction(initial_error_bound, "initial_error_bound")
    if error < 0:
        raise ValueError("initial_error_bound must be nonnegative")
    input_size = len(tower.steps[0].flow.input_coordinates)
    total_matrix = _identity(input_size)
    total_offset = tuple(Fraction(0) for _ in range(input_size))
    errors: list[Fraction] = []
    lipschitz: list[Fraction] = []
    for step in tower.steps:
        step_matrix = step.combined_matrix
        step_offset = step.combined_offset
        total_offset = tuple(
            left + right
            for left, right in zip(
                _matvec(step_matrix, total_offset), step_offset, strict=True
            )
        )
        total_matrix = _matmul(step_matrix, total_matrix)
        step_lipschitz = step.lipschitz_bound
        error = step_lipschitz * error + step.implementation_remainder_bound
        lipschitz.append(step_lipschitz)
        errors.append(error)
    return ExactMatchingCompositionCertificate(
        endpoint_error_bound=error,
        endpoint_jacobian=total_matrix,
        total_matrix=total_matrix,
        total_offset=total_offset,
        step_error_bounds=tuple(errors),
        step_lipschitz_bounds=tuple(lipschitz),
    )


def source_to_eft_contract_report(
    lane: CalculationLane,
    *,
    moment_packet: QuotientMomentPacket | None = None,
    moment_map: ExactAffineMomentMap | None = None,
    candidate: UntrustedParameterCandidate | None = None,
    left_inverse: ExactLeftInverseCertificate | None = None,
    residual_bound: Any = Fraction(0),
    source_law: FiniteSourceLawPacket | None = None,
    source_components: Sequence[SourceComponentEnvelope] = (),
    matching_tower: ExactMatchingTower | None = None,
    initial_matching_error: Any = Fraction(0),
) -> dict[str, Any]:
    """Aggregate typed mathematical checks while keeping physics fail closed.

    There is intentionally no argument for caller-authored primitive receipt
    booleans.  The module has no primitive OPH, FJ, MATCH, COV, BRST, or pole
    producer/verifier, so every such receipt is unconditionally false.
    """

    if not isinstance(lane, CalculationLane):
        raise TypeError("lane must be a CalculationLane")
    blockers: list[str] = []
    receipts: dict[str, bool] = {
        "IMPORTED_SM_VALIDATION_LANE_DECLARATION_RECEIPT": (
            lane is CalculationLane.IMPORTED_SM_STRICT_1L_VALIDATION
        ),
        "OPH_NATIVE_LANE_DECLARATION_RECEIPT": (
            lane is CalculationLane.OPH_NATIVE_STRICT_1L
        ),
        "IMPORTED_SM_STRICT_1L_VALIDATION_RECEIPT": False,
    }
    receipts.update({key: False for key in _PHYSICAL_RECEIPT_KEYS})

    coherent_report = verify_coherent_source_components(source_components)
    receipts.update(coherent_report["receipts"])
    blockers.extend(coherent_report["blockers"])

    identifiability_report: dict[str, Any] | None = None
    identification_inputs = (moment_packet, moment_map, candidate, left_inverse)
    if all(value is not None for value in identification_inputs):
        assert moment_packet is not None
        assert moment_map is not None
        assert candidate is not None
        assert left_inverse is not None
        certificate = verify_exact_affine_identifiability(
            moment_packet,
            moment_map,
            candidate,
            left_inverse,
            residual_bound=residual_bound,
        )
        identifiability_report = certificate.to_report()
        receipts.update(identifiability_report["receipts"])
        blockers.extend(identifiability_report["blockers"])
    else:
        receipts["EXACT_SOURCE_PARAMETER_IDENTIFIABILITY_MATH_RECEIPT"] = False
        blockers.append("complete_identifiability_input_tuple_missing")

    source_law_result: dict[str, Any] | None = None
    if source_law is not None:
        source_law_result = finite_source_law_report(source_law)
        receipts.update(source_law_result["receipts"])
        blockers.extend(source_law_result["blockers"])
    else:
        receipts["FINITE_SOURCE_LAW_STRUCTURAL_RECEIPT"] = False
        receipts["COV1_RECEIPT"] = False
        blockers.append("finite_source_law_packet_missing")

    matching_report: dict[str, Any] | None = None
    if matching_tower is not None:
        matching_report = compose_exact_matching_tower(
            matching_tower, initial_error_bound=initial_matching_error
        ).to_report()
        receipts.update(matching_report["receipts"])
        blockers.extend(matching_report["blockers"])
    else:
        receipts["MATCHING_TOWER_EXACT_ERROR_COMPOSITION_RECEIPT"] = False
        receipts["MATCHING_TOWER_EXACT_JACOBIAN_COMPOSITION_RECEIPT"] = False
        receipts["MATCH1_RECEIPT"] = False
        blockers.append("exact_matching_tower_missing")

    cross_commitment_valid = True
    commitment_values: list[CoherentSourceCommitments] = []
    if moment_packet is not None:
        commitment_values.append(moment_packet.commitments)
    if source_law is not None:
        commitment_values.append(source_law.commitments)
    if matching_tower is not None:
        commitment_values.append(matching_tower.commitments)
    if source_components:
        commitment_values.extend(component.commitments for component in source_components)
    if commitment_values:
        baseline = commitment_values[0]
        cross_commitment_valid = all(value == baseline for value in commitment_values[1:])
    else:
        cross_commitment_valid = False
    receipts["SOURCE_TO_EFT_CROSS_ARTIFACT_COMMITMENT_RECEIPT"] = cross_commitment_valid
    if not cross_commitment_valid:
        blockers.append("source_to_eft_cross_artifact_commitment_mismatch")

    coordinate_binding_valid = bool(
        moment_map is not None
        and source_law is not None
        and matching_tower is not None
        and moment_map.parameter_names == source_law.coordinate_names
        and moment_map.parameter_names
        == matching_tower.steps[0].flow.input_coordinates
    )
    receipts["SOURCE_TO_MATCHING_COORDINATE_BINDING_RECEIPT"] = coordinate_binding_valid
    if not coordinate_binding_valid:
        blockers.append("source_parameter_coordinate_binding_mismatch_or_missing")

    source_semantics_binding_valid = False
    if moment_packet is not None and source_law is not None:
        same_branches = (
            moment_packet.admissible_branch_hashes == source_law.admissible_branch_hashes
        )
        if moment_packet.semantics is MomentSemantics.DETERMINISTIC_SETTLED_STATE:
            same_semantics = source_law.semantics is SourceLawSemantics.DETERMINISTIC
            same_law_reference = moment_packet.source_law_hash is None
        else:
            same_semantics = source_law.semantics is SourceLawSemantics.STOCHASTIC
            same_law_reference = moment_packet.source_law_hash == source_law.source_law_hash
        source_semantics_binding_valid = bool(
            same_branches and same_semantics and same_law_reference
        )
    receipts["QUOTIENT_MOMENT_SOURCE_LAW_BINDING_RECEIPT"] = (
        source_semantics_binding_valid
    )
    if not source_semantics_binding_valid:
        blockers.append("quotient_moment_source_law_binding_mismatch_or_missing")

    # Reassert the physical boundary after merging every subreport.  This is
    # intentionally redundant: no future structural report can accidentally
    # promote a physical receipt through a key collision.
    receipts.update({key: False for key in _PHYSICAL_RECEIPT_KEYS})
    if lane is CalculationLane.IMPORTED_SM_STRICT_1L_VALIDATION:
        blockers.append("imported_sm_validation_engine_not_implemented_in_this_contract")
        lane_boundary = (
            "Imported, convention-labelled inputs may validate a separate field-theory "
            "engine, but they are never OPH predictions."
        )
    else:
        blockers.extend(
            (
                "oph_native_source_primitive_producers_missing",
                "fj1_match1_cov1_brst1_and_physical_pole_verifiers_missing",
            )
        )
        lane_boundary = (
            "The OPH-native lane consumes only replay-verified primitive source artifacts. "
            "No such complete producer chain exists in this module."
        )

    return {
        "schema_version": SOURCE_TO_EFT_SCHEMA_VERSION,
        "lane": lane.value,
        "status": CONDITIONAL_POLE_STATUS,
        "lane_boundary": lane_boundary,
        "identifiability": identifiability_report,
        "source_law": source_law_result,
        "coherent_source": coherent_report,
        "matching_tower": matching_report,
        "receipts": receipts,
        "blockers": sorted(set(blockers)),
    }


__all__ = [
    "COMMITMENT_FIELDS",
    "CONDITIONAL_POLE_STATUS",
    "CalculationLane",
    "CoherentSourceCommitments",
    "DeterministicErrorBudget",
    "EFTIntervalCensus",
    "ExactAffineMomentMap",
    "ExactAffineTransformation",
    "ExactIdentifiabilityCertificate",
    "ExactLeftInverseCertificate",
    "ExactMatchingCompositionCertificate",
    "ExactMatchingStep",
    "ExactMatchingTower",
    "FiniteSourceLawPacket",
    "MomentSemantics",
    "PerturbativeMask",
    "QUOTIENT_MOMENT_SCHEMA_VERSION",
    "QuotientMomentPacket",
    "REQUIRED_SOURCE_COMPONENTS",
    "SOURCE_LAW_SCHEMA_VERSION",
    "SOURCE_TO_EFT_SCHEMA_VERSION",
    "SourceComponentEnvelope",
    "SourceLawSemantics",
    "UntrustedParameterCandidate",
    "compose_exact_matching_tower",
    "finite_source_law_report",
    "matching_field_census_hash",
    "matching_scheme_hash",
    "matching_threshold_hash",
    "source_to_eft_contract_report",
    "verify_coherent_source_components",
    "verify_common_perturbative_mask",
    "verify_exact_affine_identifiability",
    "verify_finite_source_law_artifact",
]
