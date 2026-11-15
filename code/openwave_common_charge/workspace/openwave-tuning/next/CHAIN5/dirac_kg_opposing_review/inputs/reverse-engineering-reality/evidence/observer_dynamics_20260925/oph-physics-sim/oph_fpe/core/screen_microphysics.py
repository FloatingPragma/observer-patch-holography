from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from typing import Any

from oph_fpe.constants.oph_pixel import (
    screen_radius_planck,
    total_area_planck,
    total_entropy_capacity,
)
from oph_fpe.core.pixel_scale import PixelScale, pixel_scale_from_config
from oph_fpe.core.screen_ports import (
    echosahedral_patch_architecture_report,
    echosahedral_port_names,
)


@dataclass(frozen=True)
class ScreenMicrophysics:
    """Fixed-cutoff OPH screen architecture for a run.

    This is the simulator's finite carrier description: a cellulated S2 chart
    with echosahedral twelve-port patches, explicit cell area from P, and
    finite edge-sector bookkeeping. It is still a regulator chart, not a
    microscopic bulk lattice.
    """

    pixel_scale: PixelScale
    patch_count: int
    edge_count: int
    ports_per_patch: int = 12
    chart: str = "support_visible_s2_cellulation"
    carrier: str = "federated_echosahedral_patch"
    cap_family: str = "round_caps_on_s2"
    edge_sector_law: str = "fixed_cutoff_heat_kernel_casimir_surrogate"
    screen_units_mode: str = "numerical_regulator"

    @property
    def cell_area(self) -> float:
        return self.pixel_scale.a_cell

    @property
    def screen_area(self) -> float:
        return self.patch_count * self.cell_area

    @property
    def effective_radius(self) -> float:
        return sqrt(self.screen_area / (4.0 * pi))

    @property
    def physical_radius_planck(self) -> float | None:
        if self.screen_units_mode != "physical_cell_toy_universe":
            return None
        return screen_radius_planck(self.patch_count, self.pixel_scale.cell_area_planck)

    @property
    def angular_cell_area(self) -> float:
        if self.patch_count <= 0:
            return 0.0
        return 4.0 * pi / self.patch_count

    @property
    def port_budget(self) -> int:
        return self.patch_count * self.ports_per_patch

    @property
    def routed_port_fraction(self) -> float:
        budget = max(1, self.port_budget)
        return (2.0 * self.edge_count) / budget

    def as_jsonable(self) -> dict[str, Any]:
        total_area = total_area_planck(self.patch_count, self.pixel_scale.cell_area_planck)
        total_entropy = total_entropy_capacity(self.patch_count, self.pixel_scale.cell_area_planck)
        physical_mode = self.screen_units_mode == "physical_cell_toy_universe"
        local_architecture = echosahedral_patch_architecture_report(
            self.ports_per_patch
        )
        return {
            "chart": self.chart,
            "carrier": self.carrier,
            "cap_family": self.cap_family,
            "edge_sector_law": self.edge_sector_law,
            "patch_count": self.patch_count,
            "edge_count": self.edge_count,
            "ports_per_patch": self.ports_per_patch,
            "port_names": echosahedral_port_names(self.ports_per_patch),
            "port_budget": self.port_budget,
            "federation_architecture": {
                "cell_count": self.patch_count,
                "local_body_per_cell": "echosahedral_patch",
                "local_template_reused_by_all_cells": True,
                "local_patch_architecture": local_architecture,
                "ECHOSAHEDRAL_FEDERATION_ARCHITECTURE_RECEIPT": bool(
                    self.patch_count > 0
                    and local_architecture[
                        "ECHOSAHEDRAL_LOCAL_PATCH_ARCHITECTURE_RECEIPT"
                    ]
                    is True
                ),
            },
            "routed_port_fraction": self.routed_port_fraction,
            "pixel_scale": self.pixel_scale.as_jsonable(),
            "cell_area": self.cell_area,
            "cell_area_planck": self.pixel_scale.cell_area_planck,
            "cell_entropy_capacity": self.pixel_scale.cell_entropy_capacity,
            "screen_area": self.screen_area,
            "effective_radius": self.effective_radius,
            "angular_cell_area": self.angular_cell_area,
            "screen_units": {
                "mode": self.screen_units_mode,
                "cell_area_planck": self.pixel_scale.cell_area_planck,
                "cell_entropy_capacity": self.pixel_scale.cell_entropy_capacity,
                "total_area_planck": total_area if physical_mode else None,
                "total_entropy_capacity": total_entropy if physical_mode else None,
                "physical_radius_planck": self.physical_radius_planck,
                "regulator_area_weight_sum": total_area,
                "regulator_entropy_weight_sum": total_entropy,
                "use_P_for": [
                    "area_weight_metadata",
                    "entropy_weight_metadata",
                    "cap_capacity",
                    "normalized_residual_weights",
                ],
            },
            "claim_boundary": (
                "fixed-cutoff screen microphysics architecture: cellulated S2 chart, "
                "echosahedral port budget, records, and edge-sector bookkeeping. "
                "P normalizes local area/capacity; BW remains lambda_C(2*pi*t)."
            ),
        }


def ports_per_patch_from_config(config: dict[str, Any]) -> int:
    screen = config.get("screen", {}) or {}
    return int(screen.get("ports_per_patch", screen.get("echosahedral_ports", 12)))


def screen_microphysics_from_config(
    config: dict[str, Any], patch_count: int, edge_count: int
) -> ScreenMicrophysics:
    screen = config.get("screen", {}) or {}
    screen_units = config.get("screen_units", {}) or {}
    return ScreenMicrophysics(
        pixel_scale=pixel_scale_from_config(config),
        patch_count=int(patch_count),
        edge_count=int(edge_count),
        ports_per_patch=ports_per_patch_from_config(config),
        chart=str(screen.get("chart", "support_visible_s2_cellulation")),
        carrier=str(screen.get("carrier", "federated_echosahedral_patch")),
        cap_family=str(screen.get("cap_family", "round_caps_on_s2")),
        edge_sector_law=str(screen.get("edge_sector_law", "fixed_cutoff_heat_kernel_casimir_surrogate")),
        screen_units_mode=str(screen_units.get("mode", "numerical_regulator")),
    )
