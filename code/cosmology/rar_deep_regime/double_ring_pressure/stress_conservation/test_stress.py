import json
from math import exp, inf, pi, sqrt
from pathlib import Path
import sys

from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from cosmology.rar_deep_regime.double_ring_pressure.stress_conservation.measure import (
    HERE,
    run,
    sightline_weight,
)


def test_projection_and_conservation_with_an_outward_rising_column():
    inner, outer, length = 1.0, 1.75, 4.0

    def radial(r):
        return r**2 * exp(-((r / length) ** 2))

    def derivative(r):
        return 2 * r * (1 - (r / length) ** 2) * exp(-((r / length) ** 2))

    def enthalpy_term(r):
        return (2 * radial(r) + exp(-((r / length) ** 2))) * r / (1 + r**2) ** 1.5

    def integrate(function):
        return sum(quad(function, a, b)[0] for a, b in [(inner, outer), (outer, inf)])

    # Independent sightline integration. Both local radial pressure and both columns are positive.
    def column(r):
        return 2 * quad(lambda z: radial(sqrt(r * r + z * z)), 0, inf)[0]

    difference = column(outer) - column(inner)
    assert difference > 0

    def weight(r):
        return sightline_weight(r, inner, outer)

    assert abs(integrate(lambda r: weight(r) / r) - pi * (outer - inner)) < 1e-9

    def tangential_excess(r):
        return r * (derivative(r) + enthalpy_term(r)) / 2

    weighted_excess = integrate(lambda r: weight(r) * tangential_excess(r) / r)
    gravity = integrate(lambda r: weight(r) * enthalpy_term(r))
    assert abs(difference - (2 * weighted_excess - gravity)) < 1e-9
    bound = 3 * sqrt(3) * pi * inner**2 * difference / 2
    stress_integral = 4 * pi * integrate(lambda r: r**2 * max(0, tangential_excess(r)))
    assert gravity > 0 and stress_integral > bound > 0
    optimum = sqrt(1.5) * inner
    assert abs(weight(optimum) / optimum**3 - 4 / (3 * sqrt(3) * inner**2)) < 1e-14


def test_archive_and_heavy_stellar_control():
    result = run()
    assert json.loads((HERE / "results.json").read_text()) == result
    assert result["nominal_salpeter_zero_pressure_secondary_bend_arcsec"] > 0
    light = result["spherical_chabrier_force_and_aligned_arcs"]["ring_errors_1"]
    assert light["integrated_positive_anisotropy_lower_msun"] == 24_000_000_000


def test_exterior_shell_can_only_lower_slow_kinetic_contrast():
    inner, outer = 0.9, 1.5
    for r in (outer, outer * 1.1, outer * 3, outer * 25):
        # Independent solid-angle integration gives both projected shell fractions.
        cuts = [sqrt(1 - (R / r) ** 2) for R in (inner, outer)]
        densities = [
            quad(lambda mu: float(mu >= c), 0, 1, points=[c])[0] / (pi * R**2)
            for R, c in zip((inner, outer), cuts, strict=True)
        ]
        weight = densities[0] - densities[1]
        kinetic = sightline_weight(r, inner, outer) / (4 * pi * r**3)
        assert abs(kinetic / -weight - (1 + cuts[0]) * (1 + cuts[1]) / 2) < 1e-8
        assert weight + kinetic / 2 < 0
    branch = run()["spherical_chabrier_force_and_aligned_arcs"]["ring_errors_1"]
    assert branch["necessary_kinetic_speed_ceiling_in_thousandths_of_c"] == 343
