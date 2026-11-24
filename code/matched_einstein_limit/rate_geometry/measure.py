"""Static isotropic clock and signal comparison with an independent rim clock."""
from dataclasses import dataclass
from fractions import Fraction
from math import ceil, log2, pi, sqrt

import numpy as np

from ..seam_overlap.protected import audit


@dataclass(frozen=True)
class Source:
    """Clock rate, reciprocal signal speed, and central vacuum normalization."""

    energy_density: float
    spatial_source_density: float

    def rates(self, r):
        clock = 1 + self.energy_density * r**2 / 8
        speed = clock / (1 - self.spatial_source_density * r**2 / 8)
        return clock, speed


@dataclass(frozen=True)
class Tapes:
    radial_echoes: np.ndarray
    interior_laps: np.ndarray
    rim_lap: int
    central_ticks_per_emitter_unit: int


def expose(source, aperture, shells, tick):
    """A shell has constant rates at its midpoint and two interior ring guides."""
    width = aperture / shells
    centers = (np.arange(shells) + 0.5) * width
    ring_radii = centers[:, None] + np.array([-1, 1]) * width / sqrt(12)
    clock, speed = source.rates(centers)
    rim_clock, rim_speed = source.rates(aperture)
    radial = np.rint((2 * width / speed) * clock / tick).astype(np.int64)
    angular = np.rint((2 * pi * ring_radii / speed[:, None]) * clock[:, None] / tick)
    rim_lap = round((2 * pi * aperture / rim_speed) * rim_clock / tick)
    central_ticks = round(1 / (rim_clock * tick))
    return Tapes(radial, angular.astype(np.int64), rim_lap, central_ticks)


def interior_volume(tapes, tick):
    lengths = tapes.radial_echoes * (tick / 2)
    areas = (tapes.interior_laps * tick)**2 / pi
    return lengths @ areas.mean(axis=1)


def match_volume(source, ell, shells, tick):
    """Choose the aperture with interior records alone, then read its rim."""
    target = (4 * pi / 3) * ell**3
    lo, hi = 0.8 * ell, 1.2 * ell

    def survey(aperture):
        return expose(source, aperture, shells, tick)

    if not interior_volume(survey(lo), tick) < target < interior_volume(survey(hi), tick):
        raise ValueError("Volume outside aperture bracket")
    for _ in range(ceil(log2((hi - lo) / tick))):
        mid = (lo + hi) / 2
        if interior_volume(survey(mid), tick) < target:
            lo = mid
        else:
            hi = mid
    tapes = survey((lo + hi) / 2)
    area = (tapes.rim_lap * tick)**2 / pi
    return interior_volume(tapes, tick), area, 1 / (tapes.central_ticks_per_emitter_unit * tick)


def run():
    stored = audit()["protected_logical_registers"]["energy_contrast"]
    held = float(Fraction(stored))
    rows = []
    for level, ell in enumerate((0.4, 0.2, 0.1, 0.05)):
        shells = 32 * 2**level
        tick = ell**5 / 2**24
        _, vacuum_area, _ = match_volume(Source(0, 0), ell, shells, tick)
        factor = 4 * pi * ell**4 / 15
        for name, density, spatial in (("ordinary", 1, 1), ("protected", held, held),
                                       ("protected rate only", held, 0), ("ordinary rate only", 1, 0)):
            volume, area, rate = match_volume(Source(density, spatial), ell, shells, tick)
            time_curvature = (vacuum_area - area) / factor
            space_curvature = 4 * (rate - 1) / ell**2 - time_curvature / 3
            rows.append(dict(ell=ell, arm=name, rim_clock=rate, area_shift=area-vacuum_area,
                             volume_error=volume-(4*pi/3)*ell**3,
                             inferred_coupling=time_curvature/density,
                             spatial_component_per_density=space_curvature/density))
    return rows


if __name__ == "__main__":
    for row in run():
        print(("{ell:.2f}  {arm:20s}  {rim_clock:.9f}  {area_shift:+.10g}  "
               "{inferred_coupling:.7f}  {spatial_component_per_density:+.7f}").format(**row))
