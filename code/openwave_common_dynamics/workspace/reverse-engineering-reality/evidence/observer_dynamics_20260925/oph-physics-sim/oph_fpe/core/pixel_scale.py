from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Any

from oph_fpe.constants.oph_pixel import (
    OPHPixelConstants,
    PIXEL_MODE_ENDPOINT_CALIBRATED,
    PIXEL_MODE_SOURCE_CANDIDATE,
    P_STAR,
    PixelParameterProfile,
    pixel_constants_for_mode,
    pixel_parameter_profile,
)


DEFAULT_PIXEL_RATIO = P_STAR


@dataclass(frozen=True)
class PixelScale:
    """OPH local pixel-ratio scale used by the quantitative closure branch.

    In the compact paper, P belongs to the D10/Phase-II quantitative closure
    surface. It supplies the local ruler and coupling readout scale. It should
    be recorded for OPH-FPE runs, but it is not by itself a BW/Lorentz
    dimension-forcing mechanism.
    """

    ratio_p: float = DEFAULT_PIXEL_RATIO
    planck_area: float = 1.0
    c: float = 1.0
    hbar: float = 1.0
    k_b: float = 1.0
    source: str = "endpoint_public"
    pixel_mode: str = PIXEL_MODE_ENDPOINT_CALIBRATED
    epistemic_profile: str = PixelParameterProfile.HIERARCHY_PUBLIC.value

    @property
    def constants(self) -> OPHPixelConstants:
        return OPHPixelConstants(
            P=self.ratio_p,
            source=self.source,
            pixel_mode=self.pixel_mode,
            epistemic_profile=self.epistemic_profile,
        )

    @property
    def a_cell(self) -> float:
        return self.ratio_p * self.planck_area

    @property
    def cell_area_planck(self) -> float:
        return self.ratio_p

    @property
    def cell_entropy_capacity(self) -> float:
        return self.ratio_p / 4.0

    @property
    def ellbar_shared(self) -> float:
        return self.ratio_p / 4.0

    @property
    def g_natural(self) -> float:
        return self.a_cell / (4.0 * self.ellbar_shared)

    @property
    def length_unit(self) -> float:
        return sqrt(self.a_cell)

    @property
    def time_unit(self) -> float:
        return self.length_unit / self.c

    @property
    def energy_unit(self) -> float:
        return self.hbar * self.c / self.length_unit

    @property
    def temperature_unit(self) -> float:
        return self.hbar * self.c / (self.k_b * self.length_unit)

    def as_jsonable(self) -> dict[str, Any]:
        constants = self.constants
        return {
            "P": self.ratio_p,
            "ratio_P": self.ratio_p,
            "P_source": self.source,
            "source": self.source,
            "pixel_mode": self.pixel_mode,
            "epistemic_profile": self.epistemic_profile,
            "profile": pixel_parameter_profile(self.epistemic_profile).as_jsonable(),
            "phi": constants.phi,
            "sqrt_pi": constants.sqrt_pi,
            "alpha_from_P": constants.alpha_from_P,
            "alpha_inverse_from_P": constants.alpha_inverse_from_P,
            "role": (
                "local screen-cell area and entropy/capacity normalization for the "
                "D10/Phase-II quantitative layer; not a BW normalization factor and "
                "not a dimension-forcing input. endpoint_calibrated values are excluded "
                "from recovered-core receipt claims"
            ),
            "planck_area": self.planck_area,
            "cell_area_planck": self.cell_area_planck,
            "cell_entropy_capacity": self.cell_entropy_capacity,
            "a_cell": self.a_cell,
            "ellbar_shared": self.ellbar_shared,
            "G_natural": self.g_natural,
            "unit_readout": {
                "length_unit": self.length_unit,
                "time_unit": self.time_unit,
                "energy_unit": self.energy_unit,
                "temperature_unit": self.temperature_unit,
                "c": self.c,
                "hbar": self.hbar,
                "k_B": self.k_b,
            },
        }


def pixel_scale_from_config(config: dict[str, Any]) -> PixelScale:
    oph_constants = config.get("oph_constants", {}) or {}
    pixel = config.get("pixel", {}) or {}
    pixel_mode = str(
        oph_constants.get("pixel_mode", pixel.get("pixel_mode", PIXEL_MODE_ENDPOINT_CALIBRATED))
    )
    default_profile = (
        PixelParameterProfile.SOURCE_CANDIDATE.value
        if pixel_mode == PIXEL_MODE_SOURCE_CANDIDATE
        else PixelParameterProfile.HIERARCHY_PUBLIC.value
    )
    epistemic_profile = str(
        oph_constants.get("epistemic_profile", pixel.get("epistemic_profile", default_profile))
    )
    profile_spec = pixel_parameter_profile(epistemic_profile)
    if profile_spec.profile is PixelParameterProfile.EMPIRICAL_HADRON_CLOSURE:
        raise ValueError(
            "interval-valued empirical_hadron_closure cannot be selected implicitly as the global pixel scale"
        )
    explicit_p = oph_constants.get("P", oph_constants.get("ratio_P", pixel.get("ratio_P", pixel.get("P"))))
    if explicit_p is None:
        mode_constants = pixel_constants_for_mode(pixel_mode)
        ratio_p = mode_constants.P
        default_source = mode_constants.source
        pixel_mode = mode_constants.pixel_mode
    else:
        ratio_p = explicit_p
        default_source = "endpoint_public"
    source = oph_constants.get(
        "P_source",
        oph_constants.get("source", pixel.get("source", default_source)),
    )
    return PixelScale(
        ratio_p=float(ratio_p),
        planck_area=float(pixel.get("planck_area", pixel.get("ell_P_squared", 1.0))),
        c=float(pixel.get("c", 1.0)),
        hbar=float(pixel.get("hbar", 1.0)),
        k_b=float(pixel.get("k_B", pixel.get("k_b", 1.0))),
        source=str(source),
        pixel_mode=pixel_mode,
        epistemic_profile=epistemic_profile,
    )
