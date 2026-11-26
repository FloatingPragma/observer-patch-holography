"""Independent rods, string action and volume-matched balls in the product chart."""
import argparse
from dataclasses import dataclass
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

from scipy.optimize import brentq

from .archive import source_densities

# Share exactly the ten-frame convention and exponential-slice quadrature of
# the matched Einstein theorem. The normal model CLI does not set PYTHONPATH.
CODE = Path(__file__).resolve().parents[4] / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
from matched_einstein_limit.finite_entropy_region.refinement import product_measures, quadrature  # noqa: E402
from matched_einstein_limit.frame import FRAME  # noqa: E402

KAPPA = 2 * math.pi
ETA = 1.0
RADII = (1/32, 1/64, 1/128, 1/256)
HERE = Path(__file__).resolve().parent


def _check_densities(ordinary, collar):
    if not (0 <= ordinary <= 1 and 0 <= collar <= 1):
        raise ValueError("each homogeneous string density must be between zero and one")


def curvature(ordinary, collar, *, collar_strength=1):
    """Transverse angle transported per proper area, in oriented radians."""
    _check_densities(ordinary, collar)
    if not 0 <= collar_strength <= 1:
        raise ValueError("comparison collar strengths range from zero to one")
    return KAPPA * (ordinary + collar_strength * collar)


def energy_tensor(ordinary, collar):
    """Covariant orthonormal components ordered (t,x,y,z), from tension one."""
    _check_densities(ordinary, collar)
    rho = ordinary + collar
    return ((rho, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, -rho))


def sheet_action(sheet_number, metric):
    """Action on a unit coordinate (t,z) patch at fixed material label measure."""
    determinant = metric[0][0]*metric[3][3] - metric[0][3]*metric[3][0]
    if sheet_number < 0 or determinant >= 0:
        raise ValueError("sheet number is nonnegative and its worldsheet timelike")
    return -sheet_number * math.sqrt(-determinant)


def einstein_tensor(K):
    """Contract the transverse surface Ricci curvature in the product metric."""
    ricci_diagonal = (0, K, K, 0)
    metric = (-1, 1, 1, 1)
    scalar = sum(r*g for r, g in zip(ricci_diagonal, metric, strict=True))
    return tuple(tuple(r - g*scalar/2 if i == j else 0
                       for j in range(4))
                 for i, (r, g) in enumerate(zip(ricci_diagonal, metric, strict=True)))


def project(tensor, u):
    return sum(tensor[i][j]*u[i]*u[j] for i in range(4) for j in range(4))


def spatial_stretch(u):
    return 1 + u[1]**2 + u[2]**2


def transverse_survey(K, rod_length=1/8):
    """Circle circumference, surface area and signed holonomy of a transverse disk."""
    if K < 0 or rod_length <= 0 or math.sqrt(K)*rod_length >= math.pi:
        raise ValueError("needs a positive rod and a transverse disk within its injectivity radius")
    if K == 0:
        return dict(circumference=2*math.pi*rod_length, disk_area=math.pi*rod_length**2,
                    holonomy=0.0)
    angle = math.sqrt(K)*rod_length
    # Transport in the radial frame gives 2π(1−S'_K) independently of the disk integral.
    holonomy = 4*math.pi*math.sin(angle/2)**2
    return dict(circumference=2*math.pi*math.sin(angle)/math.sqrt(K),
                disk_area=4*math.pi*(math.sin(angle/2)/math.sqrt(K))**2,
                holonomy=holonomy)


@dataclass(frozen=True)
class Ball:
    radius: float
    volume: float
    area: float


@lru_cache(maxsize=256)
def _matched(K, ell, b2, order):
    if not (0 <= K <= 2*KAPPA and 0 < ell <= 1/32 and 1 <= b2 <= 9):
        raise ValueError("survey outside the enrolled curvature, radius or frame range")
    target = 4*math.pi*ell**3/3
    if K == 0:
        return Ball(ell, target, 4*math.pi*ell**2)
    def scaled_volume(t):
        return product_measures(K, ell*t, b2, order)[0]/target - 1
    t = brentq(scaled_volume, 0.9, 1.1, xtol=1e-15)
    volume, area = product_measures(K, ell*t, b2, order)
    return Ball(ell*t, volume, area)


def measured_ball(K, ell, u, order=16):
    """Independent geometric survey at fixed proper volume for an enrolled boost."""
    return _matched(float(K), float(ell), float(spatial_stretch(u)), order)


@lru_cache(maxsize=256)
def _weighted(K, r, b2, order):
    if not (0 <= K <= 2*KAPPA and 0 < r <= 1/16 and 1 <= b2 <= 9):
        raise ValueError("weighted survey outside the enrolled normal neighborhood")
    if K == 0:
        return 4*math.pi*r**4/15
    nodes, weights = quadrature(b2, order)[-2:]
    return r * sum(w * s * product_measures(K, r*s, b2, order)[0]
                   for s, w in zip(nodes, weights, strict=True))


def weighted_volume(K, radius, u, order=16):
    """Integrate (r²−s²)/(2r) over the true proper volume, by parts."""
    return float(_weighted(float(K), float(radius), float(spatial_stretch(u)), order))


def survey_report():
    densities = source_densities()
    ordinary, collar = densities["held"]
    circles = {name: transverse_survey(curvature(*pair)) for name, pair in densities.items()}
    calibration = circles["ordinary"]
    gain = calibration["holonomy"] / calibration["disk_area"] / energy_tensor(ordinary, 0)[0][0]
    records = []
    for ell in RADII:
        for j, u in enumerate(FRAME):
            flat = measured_ball(0, ell, u)
            base = measured_ball(curvature(ordinary, 0), ell, u)
            both = measured_ball(curvature(ordinary, collar), ell, u)
            weakened = measured_ball(curvature(ordinary, collar, collar_strength=0.5), ell, u)
            missing = measured_ball(curvature(ordinary, collar, collar_strength=0), ell, u)
            theoretical = -4*math.pi/15*gain*float(project(energy_tensor(0, collar), u))
            records.append(dict(ell=ell, frame=j, b2=float(spatial_stretch(u)),
                target_volume=flat.volume, ordinary_radius=base.radius, held_radius=both.radius,
                ordinary_area=base.area, held_area=both.area,
                ordinary_weighted_volume=weighted_volume(curvature(ordinary, 0), base.radius, u),
                held_weighted_volume=weighted_volume(curvature(ordinary, collar), both.radius, u),
                ordinary_slope=(base.area-flat.area)/ell**4,
                held_slope=(both.area-base.area)/ell**4,
                half_strength_slope=(weakened.area-base.area)/ell**4,
                zero_strength_slope=(missing.area-base.area)/ell**4,
                stress_forecast=theoretical))
    return dict(schema="oph.partition_repair.parallel_sheets.v1", units="proper rod and clock",
                theta_per_sheet=KAPPA, sheet_tension=1, flat_longitudinal_axis="z",
                measured_ordinary_gain=gain, densities={key: list(pair) for key, pair in densities.items()},
                records_are_model_integrals=True, circle_surveys=circles, balls=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = survey_report()
    if args.write:
        (HERE/"geometry.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        # Floats reflect the platform's quadrature, so check physical tolerances.
        from math import isclose
        archive = json.loads((HERE/"geometry.json").read_text())
        for key, value in data.items():
            if key not in ("balls", "circle_surveys", "measured_ordinary_gain"):
                if archive[key] != value:
                    raise ValueError(f"geometry replay differs at {key}")
        if not isclose(data["measured_ordinary_gain"], archive["measured_ordinary_gain"], rel_tol=1e-13):
            raise ValueError("ordinary gain changed")
        for key, readings in data["circle_surveys"].items():
            for name, val in readings.items():
                if not isclose(val, archive["circle_surveys"][key][name], rel_tol=1e-13):
                    raise ValueError(f"circle replay differs at {key}, {name}")
        for before, after in zip(archive["balls"], data["balls"], strict=True):
            for key in after:
                if not isclose(before[key], after[key], rel_tol=2e-9, abs_tol=1e-14):
                    raise ValueError(f"ball replay differs at {key}")
        print("ordinary curvature calibration and 40 held moving-ball measurements agree")


if __name__ == "__main__":
    main()
