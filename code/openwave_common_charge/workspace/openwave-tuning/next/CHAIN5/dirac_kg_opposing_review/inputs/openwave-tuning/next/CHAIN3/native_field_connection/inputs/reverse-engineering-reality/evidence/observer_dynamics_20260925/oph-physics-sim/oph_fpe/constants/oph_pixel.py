from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Any

import numpy as np


ALPHA_INV_SOURCE_MAP = 136.99483517741294
ALPHA_INV_GAUGE_WIDTH_MAP = 137.03566013694658
ALPHA_INV_MEASURED_ENDPOINT_COMPARISON = 137.035999177
P_SOURCE_MAP = (1.0 + math.sqrt(5.0)) / 2.0 + math.sqrt(math.pi) / ALPHA_INV_SOURCE_MAP
P_GAUGE_WIDTH_MAP = (1.0 + math.sqrt(5.0)) / 2.0 + math.sqrt(math.pi) / ALPHA_INV_GAUGE_WIDTH_MAP
P_MEASURED_ENDPOINT_COMPARISON = 1.6309682094039593

# Compatibility aliases used by older run configurations.  ``P_STAR`` is a
# reproducibility value, not an epistemically valid source default.
ALPHA_INV_SOURCE_CANDIDATE = ALPHA_INV_SOURCE_MAP
ALPHA_INV_ENDPOINT_CALIBRATED = ALPHA_INV_MEASURED_ENDPOINT_COMPARISON
P_SOURCE_CANDIDATE = P_SOURCE_MAP
P_ENDPOINT_CALIBRATED = P_MEASURED_ENDPOINT_COMPARISON
P_STAR = P_MEASURED_ENDPOINT_COMPARISON

P_EMPIRICAL_HADRON_CLOSURE_CENTRAL = 1.6310415203592243
P_EMPIRICAL_HADRON_CLOSURE_INTERVAL = (1.6310314631270647, 1.6310515775913838)

PIXEL_MODE_SOURCE_CANDIDATE = "source_candidate"
PIXEL_MODE_ENDPOINT_CALIBRATED = "endpoint_calibrated"
PIXEL_MODE_SOURCE_MAP = "source_map"
PIXEL_MODE_GAUGE_WIDTH_MAP = "gauge_width_map"
PIXEL_MODE_MEASURED_ENDPOINT_COMPARISON = "measured_endpoint_comparison"


class PixelParameterProfile(str, Enum):
    SOURCE_MAP = "source_map"
    GAUGE_WIDTH_MAP = "gauge_width_map"
    SOURCE_CANDIDATE = "source_candidate"
    HIERARCHY_PUBLIC = "hierarchy_public"
    EMPIRICAL_HADRON_CLOSURE = "empirical_hadron_closure"
    MEASURED_COMPARISON = "measured_comparison"


@dataclass(frozen=True)
class PixelProfileSpec:
    profile: PixelParameterProfile
    P: float | None
    interval: tuple[float, float] | None
    epistemic_class: str
    simulation_role: str
    recovered_core_allowed: bool = False
    public_prediction_allowed: bool = False
    generative_receipt_allowed: bool = False

    def as_jsonable(self) -> dict[str, Any]:
        return {
            "profile": self.profile.value,
            "P": self.P,
            "interval": list(self.interval) if self.interval is not None else None,
            "epistemic_class": self.epistemic_class,
            "simulation_role": self.simulation_role,
            "recovered_core_allowed": self.recovered_core_allowed,
            "public_prediction_allowed": self.public_prediction_allowed,
            "generative_receipt_allowed": self.generative_receipt_allowed,
        }


PIXEL_PARAMETER_PROFILES: dict[PixelParameterProfile, PixelProfileSpec] = {
    PixelParameterProfile.SOURCE_MAP: PixelProfileSpec(
        PixelParameterProfile.SOURCE_MAP,
        P_SOURCE_MAP,
        None,
        "source_map_fixed_point",
        "source-side mathematical map; physical Thomson promotion remains blocked",
        generative_receipt_allowed=True,
    ),
    PixelParameterProfile.GAUGE_WIDTH_MAP: PixelProfileSpec(
        PixelParameterProfile.GAUGE_WIDTH_MAP,
        P_GAUGE_WIDTH_MAP,
        None,
        "gauge_width_map_fixed_point",
        "independent mathematical map and comparison branch; no measured endpoint promotion",
        generative_receipt_allowed=True,
    ),
    PixelParameterProfile.SOURCE_CANDIDATE: PixelProfileSpec(
        PixelParameterProfile.SOURCE_CANDIDATE,
        P_SOURCE_MAP,
        None,
        "legacy_source_map_alias",
        "deprecated alias for source_map",
        generative_receipt_allowed=True,
    ),
    PixelParameterProfile.HIERARCHY_PUBLIC: PixelProfileSpec(
        PixelParameterProfile.HIERARCHY_PUBLIC,
        P_ENDPOINT_CALIBRATED,
        None,
        "endpoint_calibrated_hierarchy_surface",
        "legacy public hierarchy comparison",
    ),
    PixelParameterProfile.EMPIRICAL_HADRON_CLOSURE: PixelProfileSpec(
        PixelParameterProfile.EMPIRICAL_HADRON_CLOSURE,
        P_EMPIRICAL_HADRON_CLOSURE_CENTRAL,
        P_EMPIRICAL_HADRON_CLOSURE_INTERVAL,
        "external_data_hadron_closure",
        "empirical uncertainty propagation only",
    ),
    PixelParameterProfile.MEASURED_COMPARISON: PixelProfileSpec(
        PixelParameterProfile.MEASURED_COMPARISON,
        P_ENDPOINT_CALIBRATED,
        None,
        "measured_endpoint_comparison",
        "comparison only; forbidden in generative receipts",
    ),
}


@dataclass(frozen=True)
class OPHPixelConstants:
    """OPH local pixel/cell constants.

    P is the local screen-cell area in Planck-area units. It calibrates
    area/capacity weights and later quantitative readouts. It is not used to
    change the BW normalization or force a 3D dimension estimate.
    """

    P: float = P_STAR
    source: str = "endpoint_public"
    pixel_mode: str = PIXEL_MODE_ENDPOINT_CALIBRATED
    epistemic_profile: str = PixelParameterProfile.HIERARCHY_PUBLIC.value
    phi: float = (1.0 + math.sqrt(5.0)) / 2.0
    sqrt_pi: float = math.sqrt(math.pi)

    @property
    def alpha_from_P(self) -> float:
        return (self.P - self.phi) / self.sqrt_pi

    @property
    def alpha_inverse_from_P(self) -> float:
        alpha = self.alpha_from_P
        if abs(alpha) < 1e-30:
            return float("inf")
        return 1.0 / alpha

    @property
    def cell_area_planck(self) -> float:
        return self.P

    @property
    def cell_entropy_capacity(self) -> float:
        return self.P / 4.0

    def as_jsonable(self) -> dict[str, Any]:
        return {
            "P": self.P,
            "P_source": self.source,
            "pixel_mode": self.pixel_mode,
            "epistemic_profile": self.epistemic_profile,
            "profile": pixel_parameter_profile(self.epistemic_profile).as_jsonable(),
            "phi": self.phi,
            "sqrt_pi": self.sqrt_pi,
            "alpha_from_P": self.alpha_from_P,
            "alpha_inverse_from_P": self.alpha_inverse_from_P,
            "source_candidate": {
                "alpha_inverse": ALPHA_INV_SOURCE_MAP,
                "P": P_SOURCE_MAP,
            },
            "gauge_width_map": {
                "alpha_inverse": ALPHA_INV_GAUGE_WIDTH_MAP,
                "P": P_GAUGE_WIDTH_MAP,
            },
            "endpoint_calibrated": {
                "alpha_inverse": ALPHA_INV_MEASURED_ENDPOINT_COMPARISON,
                "P": P_MEASURED_ENDPOINT_COMPARISON,
            },
            "role": (
                "local screen-cell area and entropy/capacity normalization; "
                "not a BW normalization factor and not a dimension-forcing input. "
                "endpoint_calibrated values must stay out of recovered-core receipts"
            ),
        }


def pixel_constants_for_mode(pixel_mode: str) -> OPHPixelConstants:
    mode = str(pixel_mode).lower().replace("-", "_")
    if mode == PIXEL_MODE_SOURCE_CANDIDATE:
        return OPHPixelConstants(
            P=P_SOURCE_MAP,
            source="source_candidate",
            pixel_mode=PIXEL_MODE_SOURCE_CANDIDATE,
            epistemic_profile=PixelParameterProfile.SOURCE_CANDIDATE.value,
        )
    if mode == PIXEL_MODE_SOURCE_MAP:
        return OPHPixelConstants(
            P=P_SOURCE_MAP,
            source="source_map_fixed_point",
            pixel_mode=PIXEL_MODE_SOURCE_MAP,
            epistemic_profile=PixelParameterProfile.SOURCE_MAP.value,
        )
    if mode == PIXEL_MODE_GAUGE_WIDTH_MAP:
        return OPHPixelConstants(
            P=P_GAUGE_WIDTH_MAP,
            source="gauge_width_map_fixed_point",
            pixel_mode=PIXEL_MODE_GAUGE_WIDTH_MAP,
            epistemic_profile=PixelParameterProfile.GAUGE_WIDTH_MAP.value,
        )
    if mode in {PIXEL_MODE_ENDPOINT_CALIBRATED, "endpoint_public"}:
        return OPHPixelConstants(
            P=P_ENDPOINT_CALIBRATED,
            source="endpoint_public",
            pixel_mode=PIXEL_MODE_ENDPOINT_CALIBRATED,
            epistemic_profile=PixelParameterProfile.HIERARCHY_PUBLIC.value,
        )
    raise ValueError(f"unknown OPH pixel mode: {pixel_mode}")


def pixel_parameter_profile(profile: str | PixelParameterProfile) -> PixelProfileSpec:
    try:
        key = profile if isinstance(profile, PixelParameterProfile) else PixelParameterProfile(str(profile))
    except ValueError as exc:
        raise ValueError(f"unknown OPH pixel parameter profile: {profile}") from exc
    return PIXEL_PARAMETER_PROFILES[key]


def pixel_constants_for_profile(profile: str | PixelParameterProfile) -> OPHPixelConstants:
    spec = pixel_parameter_profile(profile)
    if spec.profile is PixelParameterProfile.EMPIRICAL_HADRON_CLOSURE:
        raise ValueError(
            "empirical_hadron_closure is an interval-valued comparison profile; "
            "select a sampled value explicitly and preserve its interval provenance"
        )
    if spec.P is None:
        raise ValueError(f"pixel profile {spec.profile.value} has no scalar P")
    return OPHPixelConstants(
        P=spec.P,
        source=spec.epistemic_class,
        pixel_mode=spec.profile.value,
        epistemic_profile=spec.profile.value,
    )


def require_generative_pixel_profile(profile: str | PixelParameterProfile) -> PixelProfileSpec:
    """Reject measured/empirical P profiles at a source-producing boundary."""
    spec = pixel_parameter_profile(profile)
    if not spec.generative_receipt_allowed:
        raise ValueError(
            f"pixel profile {spec.profile.value} is comparison-only and cannot enter a generative receipt"
        )
    return spec


def pixel_value_provenance(
    P: float,
    declared_profile: str | PixelParameterProfile | None,
    *,
    tolerance: float = 1.0e-15,
) -> dict[str, Any]:
    """Audit whether a numeric ``P`` is bound to a generative profile.

    Numeric equality alone is not provenance.  In particular, the historical
    ``P_STAR`` value remains usable for replay/comparison, but it cannot enter a
    source-producing receipt unless a declared generative profile both permits
    that role and canonically owns the supplied value.
    """

    value = float(P)
    tol = float(tolerance)
    if not math.isfinite(value):
        raise ValueError("P must be finite")
    if not math.isfinite(tol) or tol < 0.0:
        raise ValueError("tolerance must be finite and nonnegative")
    if declared_profile is None:
        return {
            "P": value,
            "declared_profile": None,
            "profile_known": False,
            "canonical_P": None,
            "value_matches_profile": False,
            "generative_receipt_allowed": False,
            "GENERATIVE_PIXEL_PROFILE_RECEIPT": False,
            "blockers": ["pixel_profile_not_declared"],
        }
    try:
        spec = pixel_parameter_profile(declared_profile)
    except ValueError:
        return {
            "P": value,
            "declared_profile": str(declared_profile),
            "profile_known": False,
            "canonical_P": None,
            "value_matches_profile": False,
            "generative_receipt_allowed": False,
            "GENERATIVE_PIXEL_PROFILE_RECEIPT": False,
            "blockers": ["pixel_profile_unknown"],
        }
    canonical = spec.P
    value_matches = bool(
        canonical is not None
        and math.isclose(value, float(canonical), rel_tol=0.0, abs_tol=tol)
    )
    allowed = bool(spec.generative_receipt_allowed and value_matches)
    blockers = []
    if canonical is None:
        blockers.append("pixel_profile_has_no_scalar_value")
    elif not value_matches:
        blockers.append("pixel_value_does_not_match_declared_profile")
    if not spec.generative_receipt_allowed:
        blockers.append("pixel_profile_is_comparison_only")
    return {
        "P": value,
        "declared_profile": spec.profile.value,
        "profile_known": True,
        "canonical_P": canonical,
        "epistemic_class": spec.epistemic_class,
        "simulation_role": spec.simulation_role,
        "value_matches_profile": value_matches,
        "generative_receipt_allowed": spec.generative_receipt_allowed,
        "GENERATIVE_PIXEL_PROFILE_RECEIPT": allowed,
        "blockers": blockers,
    }


def equal_cell_area_planck(patch_count: int, P: float = P_STAR) -> np.ndarray:
    return np.full(int(patch_count), float(P), dtype=np.float64)


def equal_cell_entropy(patch_count: int, P: float = P_STAR) -> np.ndarray:
    return np.full(int(patch_count), float(P) / 4.0, dtype=np.float64)


def total_area_planck(patch_count: int, P: float = P_STAR) -> float:
    return float(int(patch_count) * float(P))


def total_entropy_capacity(patch_count: int, P: float = P_STAR) -> float:
    return float(int(patch_count) * float(P) / 4.0)


def screen_radius_planck(patch_count: int, P: float = P_STAR) -> float:
    return math.sqrt((int(patch_count) * float(P)) / (4.0 * math.pi))


def scale_factor_from_patch_count(N_t: int, N_0: int) -> float:
    if int(N_0) <= 0:
        raise ValueError("N_0 must be positive")
    return math.sqrt(int(N_t) / int(N_0))


def cap_area_planck(cap_weights: np.ndarray, cell_area_planck: float | np.ndarray = P_STAR) -> float:
    weights = np.asarray(cap_weights, dtype=np.float64)
    area = np.asarray(cell_area_planck, dtype=np.float64)
    return float(np.sum(area * weights))


def cap_entropy_capacity(cap_weights: np.ndarray, cell_entropy: float | np.ndarray = P_STAR / 4.0) -> float:
    weights = np.asarray(cap_weights, dtype=np.float64)
    entropy = np.asarray(cell_entropy, dtype=np.float64)
    return float(np.sum(entropy * weights))
