"""Signed radial-pressure contrast from two rings and measured stellar populations."""

import argparse
from itertools import product
import json
from math import ceil, floor, inf, nextafter, pi, sqrt
from pathlib import Path

from mpmath import iv
import numpy as np
from scipy.integrate import cumulative_simpson, quad, simpson

from ..aperture_virial.measure import G, KPC_KM, Tracer
from ..stellar_lensing import ARCSEC_PER_RADIAN, C_KM_S

HERE = Path(__file__).parent
Z = (0.222, 0.609, 2.035)
THETA = (1.43, 2.07)
CORE = np.array([0.066, 0.082])
EDGE = np.array([0.50, 6.05])
LOG_MASSES = {"Chabrier": (10.85, 11.27), "Salpeter": (11.13, 11.52)}


def acceleration():
    receipt = HERE.parent / "receipts/sparc_deep_regime_diagnostic.json"
    values = json.loads(receipt.read_text())["deep_radial_acceleration"]
    return values["a0_fixed_exponent_m_s2"], values[
        "a0_galaxy_cluster_bootstrap_95pct_m_s2"
    ]


def fractions(radius):
    projected = 1 - (np.hypot(radius, EDGE) - np.hypot(radius, CORE)) / (EDGE - CORE)
    spherical = (
        2
        / pi
        * (EDGE * np.arctan(radius / EDGE) - CORE * np.arctan(radius / CORE))
        / (EDGE - CORE)
    )
    return spherical, projected


def geometry():
    chi = [
        C_KM_S / 70 * 1000 * quad(lambda y: 1 / sqrt(0.3 * (1 + y) ** 3 + 0.7), 0, z)[0]
        for z in Z
    ]
    dl = chi[0] / (1 + Z[0])
    critical = [
        C_KM_S**2 / (4 * pi * G * dl * (1 - chi[0] / ch)) / 1e6 for ch in chi[1:]
    ]
    return dl / ARCSEC_PER_RADIAN, critical, 1 - chi[1] / chi[2]


def fibre(log_masses, a0, scale):
    tracer = Tracer(size=1201)
    angular = tracer.radius[:, None]
    light = 10 ** (-0.4 * np.array([18.75, 17.15]))
    # Each spherical brightness component integrates to one before weighting.
    # Radii are expressed in arcseconds. A shared normalization cancels.
    tracer.nu = (
        (1 / (angular**2 + CORE**2) - 1 / (angular**2 + EDGE**2)) / (EDGE - CORE)
    ) @ (light / sum(light))
    aperture = tracer.aperture(1.0)
    masses = 10 ** np.array(log_masses)
    enclosed = fractions(angular)[0] @ masses
    vc2 = G * enclosed / (scale * tracer.radius) + np.sqrt(
        G * enclosed * a0 * KPC_KM / 1000
    )
    # Direct isotropic Jeans integral j_r(r) = integral_r^infty nu*vc² dlog(r).
    pressure = cumulative_simpson(
        (tracer.nu * vc2)[::-1], x=-tracer.logr[::-1], initial=0
    )[::-1]
    weights = tracer.radius**3 * (aperture[0] + aperture[1])
    numerator = simpson(weights * pressure, x=tracer.logr)
    light_in_fibre = simpson(weights * tracer.nu, x=tracer.logr)
    return dict(
        isotropic_km_s=sqrt(numerator / light_in_fibre),
        anisotropy_bounds_km_s=np.sqrt(tracer.moment(aperture, vc2)).tolist(),
    )


def nominal(log_masses, a0):
    scale, crit, efficiency = geometry()
    masses = 10 ** np.array(log_masses)
    force = a0 * (KPC_KM / 1000)
    s1, s2 = THETA

    def enclosed(theta):
        baryons = fractions(theta)[0] @ masses
        return baryons + theta * scale * sqrt(force * baryons / G)

    def projected(theta):
        baryons = fractions(theta)[1] @ masses
        radial = quad(
            lambda ang: sqrt(fractions(theta / np.cos(ang))[0] @ masses), 0, pi / 2
        )[0]
        return (baryons + theta * scale * sqrt(force / G) * radial) / (
            pi * (1000 * theta * scale) ** 2
        )

    orbital = [projected(t) for t in THETA]
    contrast = orbital[0] - orbital[1] - (crit[0] - crit[1])
    s = sqrt(1 - (s1 / s2) ** 2)
    truncated = (
        s
        / (pi * (1000 * scale) ** 2)
        * (
            quad(lambda x: enclosed(s1 / sqrt(1 - (s * x) ** 2)), 0, 1)[0] / s1**2
            - enclosed(s2) / (s2**2 * (1 + s))
        )
    )
    # SIS-equivalent scale at this ray, relevant also for a local deflection
    # of another shape. An external primary sheet is eliminated in the difference.
    b = s2 * contrast / crit[1]
    sigma = (
        C_KM_S * sqrt(b / ARCSEC_PER_RADIAN / (4 * pi * efficiency)) if b >= 0 else None
    )
    return dict(
        stellar_fibre=fibre(log_masses, a0, scale),
        unlimited_orbital_means=orbital,
        upper_pressure_contrast_unlimited=contrast,
        upper_pressure_contrast_any_positive_exterior=truncated - (crit[0] - crit[1]),
        zero_pressure_primary_sheet_kappa1=1 - orbital[0] / crit[0],
        zero_pressure_secondary_bend_arcsec=b,
        zero_pressure_intermediate_ray_arcsec=s2
        * (1 - (orbital[1] + crit[0] - orbital[0]) / crit[0]),
        zero_pressure_secondary_equivalent_km_s=sigma,
    )


def decimal(x):
    return iv.mpf(str(x))


def distance_bounds():
    # Trapezoidal error ≤ z^3 max|f''|/(12 n²). On 0≤z≤2.035 the direct
    # derivative of (0.3(1+z)^3+0.7)^(-1/2) has |f''|<60.
    n = 256
    distances = []
    for z in map(decimal, Z):
        h = z / n
        values = [
            1 / iv.sqrt(decimal(".3") * (1 + i * h) ** 3 + decimal(".7"))
            for i in range(n + 1)
        ]
        value = h * (sum(values[1:-1]) + (values[0] + values[-1]) / 2)
        error = 60 * z**3 / (12 * n * n)
        distances.append(
            (value + iv.mpf([-error.b, error.b])) * decimal(C_KM_S) / 70 * 1000
        )
    dl = distances[0] / (1 + decimal(Z[0]))
    scale = dl * iv.pi / (180 * 3600)
    critical = [
        decimal(C_KM_S) ** 2
        / (4 * iv.pi * decimal(G) * dl * (1 - distances[0] / x))
        / 10**6
        for x in distances[1:]
    ]
    return scale, critical[0] - critical[1]


def frac(t, core, edge):
    return (
        2
        / iv.pi
        * (edge * iv.atan2(t, edge) - core * iv.atan2(t, core))
        / (edge - core)
    )


def partition(lo, hi, count):
    return [
        (lo + (hi - lo) * i / count, lo + (hi - lo) * (i + 1) / count)
        for i in range(count)
    ]


def upper_orbital(scale, a0_upper, ranges, steps=192):
    # Monotonicity of the truncated projected contrast makes the smallest
    # inner and largest outer radius extremal for every allowed mass profile.
    r1, r2 = decimal(ranges[0][0]), decimal(ranges[1][1])
    s = iv.sqrt(1 - (r1 / r2) ** 2)
    radii = [r1 / iv.sqrt(1 - (s * i / steps) ** 2) for i in range(1, steps + 1)]
    core_low = tuple(map(decimal, (".056", ".072")))
    core_high = tuple(map(decimal, (".076", ".092")))
    edges = [
        partition(decimal(lo), decimal(hi), 2)
        for lo, hi in ((".45", ".55"), ("5.95", "6.15"))
    ]
    # The positivity tests below put both stellar masses at their Chabrier
    # upper limits, even when the two components vary independently.
    masses = [decimal(10) ** decimal(x) for x in ("10.94", "11.32")]
    gain = iv.sqrt(decimal(a0_upper) * decimal(KPC_KM) / 1000 / decimal(G))
    upper = iv.mpf(0)
    for (t0, T0), (t1, T1) in product(*edges):
        compact = list(zip(core_low, (t0, t1)))
        diffuse = list(zip(core_high, (T0, T1)))
        end = [frac(r2, *v) for v in diffuse]
        for a, b in zip(compact, diffuse):
            assert frac(r1, *b) / frac(r2, *a) > r1 / (r2 * (1 + s))
        samples = [[frac(r, *shape) for shape in compact] for r in radii]
        # All integrands increase with radius. A right-endpoint Riemann sum
        # over x and a lower enclosure of the outer boundary term give an
        # upper bound on the complete continuous parameter box.
        star = sum(
            (sum(v[i] for v in samples) / steps / r1**2 - end[i] / (r2**2 * (1 + s)))
            * masses[i]
            for i in (0, 1)
        )
        pos_root = (
            sum(
                r * iv.sqrt(sum(f * m for f, m in zip(v, masses)))
                for r, v in zip(radii, samples)
            )
            / steps
        )
        end_root = iv.sqrt(sum(f * m for f, m in zip(end, masses)))
        excess = pos_root / r1**2 - end_root / (r2 * (1 + s))
        assert excess > 0
        candidate = s / (iv.pi * scale**2 * 10**6) * (star + gain * scale * excess)
        upper = max(upper, candidate.b)
    return upper


def certified():
    scale, critical = distance_bounds()
    _, a0_box = acceleration()
    report = {}
    catalogue = json.loads((HERE.parent / "aperture_virial/sample.json").read_text())
    inner = next(row for row in catalogue if row["name"] == "J0946+1006")
    panels = {
        "ring_errors_1": ((1.42, 1.44), (2.05, 2.09)),
        "ring_errors_2": ((1.41, 1.45), (2.03, 2.11)),
        "extended_to_bolton_radii": (
            (min(inner["bSIE"], inner["bLTM"]), 1.45),
            (2.03, 2.11),
        ),
    }
    for name, ranges in panels.items():
        upper = upper_orbital(scale, a0_box[1], ranges)
        lower = critical - upper
        amin = decimal(ranges[0][0])
        bmax = decimal(ranges[1][1])
        weight = (1 / amin**2 - 1 / bmax**2) / (iv.pi * scale**2 * 10**6)
        report[name] = dict(
            ring_ranges_arcsec=ranges,
            orbital_contrast_upper=ceil(nextafter(float(upper.b), inf)),
            critical_contrast_lower=floor(nextafter(float(critical.a), -inf)),
            pressure_contrast_upper=ceil(nextafter(float((-lower).b), inf)),
            additional_interior_mass_lower_msun=max(
                0, 10**8 * floor(nextafter(float((lower / weight / 10**8).a), -inf))
            ),
        )
    return report


def run():
    scale, critical, efficiency = geometry()
    a0, a0_box = acceleration()
    return dict(
        radii_kpc=[scale * t for t in THETA],
        critical_msun_pc2=critical,
        critical_contrast_msun_pc2=critical[0] - critical[1],
        secondary_efficiency=efficiency,
        primary_a0_m_s2=a0,
        a0_box_m_s2=a0_box,
        nominal={name: nominal(v, a0) for name, v in LOG_MASSES.items()},
        all_positive_exteriors=certified(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = HERE / "results.json"
    data = run()
    if args.write:
        record.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(record.read_text())

        def check(actual, expected):
            if isinstance(actual, dict):
                assert actual.keys() == expected.keys()
                for key in actual:
                    check(actual[key], expected[key])
            elif actual is None:
                assert expected is None
            else:
                np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=1e-11)

        check(data, saved)
        print("Replayed the two-ring pressure contrast and the finite-radius bound.")
