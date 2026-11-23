"""Static relativistic halo observables and a finite-edge test on public shear stacks."""

import argparse
import hashlib
import json
import math
import tarfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from .. import lensing_profile


HERE = Path(__file__).resolve().parent
MAX_SHEAR_FACTOR = 4 / math.pi


@dataclass(frozen=True)
class Halo:
    """Mass lengths in geometric units, on the exterior of the compact baryons."""

    core: float
    slope: float
    edge: float

    def __post_init__(self):
        if not all(math.isfinite(x) for x in (self.core, self.slope, self.edge)):
            raise ValueError("The finite metric needs finite parameters")
        if self.core < 0 or self.slope < 0 or self.edge <= 0:
            raise ValueError("Mass terms must be nonnegative and radius positive")
        if self.edge <= 2 * self.mass(self.edge):
            raise ValueError("The outer match must lie outside the horizon")

    def mass(self, radius):
        return self.core + self.slope * min(radius, self.edge)

    def compactness(self, radius):
        if radius <= 0 or radius <= 2 * self.mass(radius):
            raise ValueError("Static instrument inside excluded region")
        return self.mass(radius) / radius

    def orbital_beta2(self, radius):
        u = self.compactness(radius)
        x = u / (1 - 2 * u)
        if x >= 1:
            raise ValueError("No timelike circular orbit at that radius")
        return x

    def log_lapse_squared(self, radius):
        u = self.compactness(radius)
        if radius >= self.edge:
            return math.log1p(-2 * u)
        factor = 1 - 2 * self.slope
        core_radius = 2 * self.core / factor
        return (
            math.log1p(-2 * self.mass(self.edge) / self.edge)
            + (
                2 * self.slope * math.log(radius / self.edge)
                + math.log1p(-core_radius / radius)
                - math.log1p(-core_radius / self.edge)
            )
            / factor
        )

    def deflection(self, turn):
        """Bending at infinity, parameterized by closest areal radius."""
        self.orbital_beta2(turn)
        at_turn = self.log_lapse_squared(turn)

        def departure(angle):
            sine = math.sin(angle)
            r = turn / sine
            denominator = (
                1
                + math.expm1(at_turn - self.log_lapse_squared(r)) / math.cos(angle) ** 2
            )
            if denominator <= 0:
                raise ValueError("Turning ray is outside the direct scattering branch")
            # Evaluate the excess over the flat polar angle without cancellation.
            return math.expm1(
                -0.5 * (math.log1p(-2 * self.compactness(r)) + math.log(denominator))
            )

        cut = [math.asin(turn / self.edge)] if turn < self.edge else []
        return (
            2
            * quad(departure, 0, math.pi / 2, points=cut, epsrel=2e-9, epsabs=3e-14)[0]
        )


def shear_factor(y):
    if y < 0 or not math.isfinite(y):
        raise ValueError("Expected finite nonnegative radius ratio")
    if y >= 1:
        return MAX_SHEAR_FACTOR / y
    return (4 * y / (1 + math.sqrt(1 - y * y)) + 2 * math.acos(y)) / math.pi


def shear_stack(profile, a0, upper_a0):
    lower = np.asarray(profile["lower_bin_edge_m_s2"])
    keep = lower >= 1e-13
    center = np.asarray(profile["g_bar_m_s2"])[keep]
    upper = 2 * center - lower[keep]
    signal = np.asarray(profile["tangential_esd"])[keep]
    cross = np.asarray(profile["cross_esd"])[keep]
    cov = np.asarray(profile["covariance"])[np.ix_(keep, keep)]
    weights = 1 / cov.diagonal()
    weights /= weights.sum()
    sigma = math.sqrt(weights @ cov @ weights)
    limits = {}
    for name, b, scale in (
        ("centres", center, a0),
        ("bin_edge", upper, a0),
        ("mass_and_bin_edge", 10**0.2 * upper, upper_a0),
    ):
        allowed = (
            4 / math.pi * b + MAX_SHEAR_FACTOR * np.sqrt(b * scale)
        ) / lensing_profile.FOUR_G
        gap = float(weights @ (signal - allowed))
        limits[name] = dict(
            upper_by_bin=allowed.tolist(),
            weighted_upper=float(weights @ allowed),
            observed_minus_upper=gap,
            gap_over_statistical_se=gap / sigma,
        )
    # Every positive mixture is bounded by its most favorable sharp edge.
    # Convexity between consecutive knots makes this a finite maximum.
    core = 4 / math.pi * center / lensing_profile.FOUR_G
    amplitude = np.sqrt(a0 * center) / lensing_profile.FOUR_G
    bounds = []
    for cutoff in (0.0, *center):
        ratios = [shear_factor(math.sqrt(cutoff / b)) for b in center]
        bounds.append(float(weights @ (core + amplitude * ratios)))
    common_bin_limits = {}
    for label, mass, scale in (
        ("nominal", 1.0, a0),
        ("high_mass_calibration", 10**0.2, upper_a0),
    ):
        low = mass * lower[keep]
        high = mass * upper
        relaxed, sharp = [], []
        for edge in (0.0, *np.unique(np.r_[low, high])):
            at_edge = np.clip(edge, low, high)

            def factors(b):
                return np.array([shear_factor(math.sqrt(edge / x)) for x in b])

            prediction = (
                4 / math.pi * high + np.sqrt(high * scale) * factors(at_edge)
            ) / lensing_profile.FOUR_G
            relaxed.append(float(weights @ prediction))
            profiles = [
                (4 / math.pi * b + np.sqrt(b * scale) * factors(b))
                / lensing_profile.FOUR_G
                for b in (low, high, at_edge)
            ]
            sharp.append(float(weights @ np.max(profiles, axis=0)))
        upper_bound = max(sharp)
        common_bin_limits[label] = dict(
            weighted_upper=upper_bound,
            relaxed_upper=max(relaxed),
            gap_over_statistical_se=float((weights @ signal - upper_bound) / sigma),
        )
    index = int(np.argmax(bounds))
    common_upper = bounds[index]
    return dict(
        common_edge_arbitrary_bin_weights=common_bin_limits,
        common_cutoff=dict(
            max_weighted_prediction=common_upper,
            maximizing_edge_acceleration=(0.0, *center)[index],
            gap_over_statistical_se=float((weights @ signal - common_upper) / sigma),
        ),
        bins=len(signal),
        centres=center.tolist(),
        positive_weights=weights.tolist(),
        weighted_signal=float(weights @ signal),
        statistical_se=sigma,
        cross_over_statistical_se=float(weights @ cross / sigma),
        envelopes=limits,
    )


def run(archive_path):
    receipt = json.loads(lensing_profile.RECEIPT.read_text())
    if (
        hashlib.sha256(archive_path.read_bytes()).hexdigest()
        != lensing_profile.ARCHIVE_SHA256
    ):
        raise ValueError("Unknown KiDS release")
    surveys = {}
    with tarfile.open(archive_path) as archive:
        for name in ("KiDS-bright", "GAMA"):
            profile = lensing_profile.load_profile(
                archive, lensing_profile.PROFILES[name]
            )
            surveys[name] = shear_stack(
                profile,
                receipt["fitted_sparc_a0_m_s2"],
                receipt["sparc_cluster_bootstrap_95pct_m_s2"][1],
            )
    rows = []
    for edge in (2.0, 10.0, 100.0, 1000.0):
        halo = Halo(0.0, 1e-6, edge)
        rows.append(
            dict(
                outer_over_turn=edge,
                orbital_beta2=halo.orbital_beta2(1),
                outer_over_inner_clock_rate=math.exp(
                    (halo.log_lapse_squared(2) - halo.log_lapse_squared(1)) / 2
                ),
                deflection_rad=halo.deflection(1),
                shear_factor=shear_factor(1 / edge),
            )
        )
    return dict(
        archive_sha256=lensing_profile.ARCHIVE_SHA256,
        sparc_a0_m_s2=receipt["fitted_sparc_a0_m_s2"],
        weak_shear_ceiling=MAX_SHEAR_FACTOR,
        public_data=lensing_profile.ARCHIVE,
        controlled_metric=rows,
        surveys=surveys,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.archive)
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for key, survey in result["surveys"].items():
        print(
            key,
            survey["weighted_signal"],
            survey["statistical_se"],
            {k: v["gap_over_statistical_se"] for k, v in survey["envelopes"].items()},
        )
