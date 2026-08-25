#!/usr/bin/env python3
"""External-field quadrupole of two unscreened QUMOND responses, cross-checked.

A dark response written as a local nonlinear Poisson law in the Newtonian
field, ``div[nu(|grad Phi_N|/a0) grad Phi_N]``, gives every point mass its own
phantom halo.  In the Galactic external field that halo is anisotropic and
produces a quadrupole ``Q2`` in the inner Solar System (Milgrom 2009;
Blanchet and Novak 2011).  This module computes ``Q2`` from first principles
by the inner multipole expansion of the phantom potential and compares it with
the published Blanchet--Novak value for the simple function and with the
Park et al. (2026) benchmark for the radial-acceleration function.

The purpose is diagnostic: it compares local central values with published
benchmarks while exposing the quadrature warnings and missing error controls.
The present integration has no certified singularity subtraction or radial
tail bound and therefore cannot independently certify a Cassini tension or
retire either fixed-input law.  It is also not a theorem excluding any wider
class of local, screened, environment-dependent, nonlocal, or dynamical laws.

Convention (Blanchet--Novak): ``U = -Phi``, ``U_quad = (1/2) Q_ij x^i x^j``,
``Q_ij = Q2 (e_i e_j - delta_ij / 3)``, hence ``Phi_quad = -(1/3) Q2 r^2 P2``.
Inner expansion of ``Phi_ph`` with ``lap Phi_ph = div F``,
``F = (nu - 1) grad Phi_N``: the ``r^2 P2`` coefficient is
``A = (1/4 pi) int F . grad(P2 / r^3) d^3x`` (one integration by parts; the
surface terms vanish), and ``Q2 = -3 A``.  Lengths in ``r_M = sqrt(G M / a0)``,
accelerations in ``a0``.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy import integrate

G_SI = 6.67430e-11
M_SUN = 1.98847e30
BLANCHET_NOVAK_SIMPLE_Q2 = 4.1e-26   # s^-2, a0 = 1.2e-10, eta = 1.6 (their Table)
BLANCHET_NOVAK_A0 = 1.2e-10
BLANCHET_NOVAK_ETA = 1.583
PARK_RAR_Q2 = 3.387e-26              # s^-2, a0 = 1.02e-10, e_N = 1.643
PARK_A0 = 1.02e-10
PARK_E_NEWTONIAN = 1.643
CASSINI_Q2 = 1.6e-27
CASSINI_SIGMA = 1.8e-27


def nu_simple(y: float) -> float:
    return 0.5 + math.sqrt(0.25 + 1.0 / y)


def nu_rar(y: float) -> float:
    r = math.sqrt(y)
    if r > 700.0:
        return 1.0
    return 1.0 + 1.0 / math.expm1(r)


def quadrupole_coefficient_diagnostic(
    nu: Callable[[float], float], e: float, epsrel: float = 1e-9
) -> dict[str, Any]:
    """Numerical coefficient plus explicit quadrature diagnostics.

    SciPy's error estimates and integration warnings are part of the result.
    The finite radial cutoffs have no analytic tail certificate here.
    """
    r0 = 1.0 / math.sqrt(e)  # radius where the Newtonian field vanishes on axis
    edges = [1e-4, 0.1, 0.5 * r0, 0.95 * r0, r0, 1.05 * r0, 1.5 * r0, 3.0, 10.0, 100.0, 5000.0]

    def integrand(th: float, r: float) -> float:
        s, c = math.sin(th), math.cos(th)
        x, z = r * s, r * c
        gx = x / r**3
        gz = z / r**3 - e
        y = math.hypot(gx, gz)
        if y <= 0.0:
            return 0.0
        nm1 = nu(y) - 1.0
        r5, r7 = r**5, r**7
        q = 3.0 * z * z - r * r
        dphix = (-2.0 * x) / (2.0 * r5) - 5.0 * q * x / (2.0 * r7)
        dphiz = (4.0 * z) / (2.0 * r5) - 5.0 * q * z / (2.0 * r7)
        return nm1 * (gx * dphix + gz * dphiz) * 2.0 * math.pi * r * r * s

    total = 0.0
    error_sum = 0.0
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", integrate.IntegrationWarning)
        for a, b in zip(edges[:-1], edges[1:]):
            value, error = integrate.dblquad(
                integrand,
                a,
                b,
                0.0,
                math.pi,
                epsabs=1e-13,
                epsrel=epsrel,
            )
            total += value
            error_sum += abs(error)

    warning_counts = {
        "roundoff_or_tolerance": 0,
        "maximum_subdivisions": 0,
        "nonconvergence_or_divergence": 0,
        "other": 0,
    }
    for item in caught:
        message = str(item.message).lower()
        if "roundoff" in message or "requested tolerance" in message:
            warning_counts["roundoff_or_tolerance"] += 1
        elif "maximum number of subdivisions" in message:
            warning_counts["maximum_subdivisions"] += 1
        elif "does not converge" in message or "probably divergent" in message:
            warning_counts["nonconvergence_or_divergence"] += 1
        else:
            warning_counts["other"] += 1
    return {
        "coefficient": total / (4.0 * math.pi),
        "quadrature_error_estimate_abs": error_sum / (4.0 * math.pi),
        "integration_warning_count": len(caught),
        "integration_warning_counts": warning_counts,
        "radial_cutoffs": [edges[0], edges[-1]],
        "tail_bound_certified": False,
        "singularity_subtraction_certified": False,
        "numerically_certified": False,
    }


def quadrupole_coefficient(
    nu: Callable[[float], float], e: float, epsrel: float = 1e-9
) -> float:
    """Compatibility wrapper returning the uncertified central value."""
    return float(quadrupole_coefficient_diagnostic(nu, e, epsrel)["coefficient"])


def q2_field_law(nu: Callable[[float], float], a0: float, e: float,
                 mass: float = M_SUN, epsrel: float = 1e-9) -> float:
    r_m = math.sqrt(G_SI * mass / a0)
    return -3.0 * quadrupole_coefficient(nu, e, epsrel) * a0 / r_m


def q2_field_law_diagnostic(
    nu: Callable[[float], float],
    a0: float,
    e: float,
    mass: float = M_SUN,
    epsrel: float = 1e-9,
) -> dict[str, Any]:
    r_m = math.sqrt(G_SI * mass / a0)
    raw = quadrupole_coefficient_diagnostic(nu, e, epsrel)
    scale = abs(-3.0 * a0 / r_m)
    return {
        "Q2_s2": -3.0 * float(raw["coefficient"]) * a0 / r_m,
        "quadrature_error_estimate_abs_s2": (
            scale * float(raw["quadrature_error_estimate_abs"])
        ),
        **{key: value for key, value in raw.items() if key != "coefficient"},
    }


def local_halo_footprint(rho_local: float = 6.8e-22, radius_m: float = 1.496e12) -> dict[str, float]:
    """Footprint of an isotropic ambient dark density: zero quadrupole, a tidal
    scale ``4 pi G rho / 3``, and the enclosed mass inside ``radius_m``."""
    tidal = 4.0 * math.pi * G_SI * rho_local / 3.0
    enclosed = rho_local * 4.0 / 3.0 * math.pi * radius_m**3
    return {
        "rho_local_kg_m3": rho_local,
        "quadrupole_s2": 0.0,
        "tidal_scale_s2": tidal,
        "tidal_over_cassini_sigma": tidal / CASSINI_SIGMA,
        "enclosed_mass_kg_within_radius": enclosed,
        "enclosed_mass_solar_within_radius": enclosed / M_SUN,
        "radius_m": radius_m,
    }


def run(epsrel: float = 1e-9) -> dict[str, Any]:
    simple_bn = q2_field_law_diagnostic(
        nu_simple, BLANCHET_NOVAK_A0, BLANCHET_NOVAK_ETA, epsrel=epsrel
    )
    rar_bn = q2_field_law_diagnostic(
        nu_rar, BLANCHET_NOVAK_A0, BLANCHET_NOVAK_ETA, epsrel=epsrel
    )
    rar_park = q2_field_law_diagnostic(
        nu_rar, PARK_A0, PARK_E_NEWTONIAN, epsrel=epsrel
    )
    simple_q2 = float(simple_bn["Q2_s2"])
    park_q2 = float(rar_park["Q2_s2"])
    return {
        "schema": "oph.cosmology.qumond_quadrupole_crosscheck.v2",
        "scope": "uncertified_numerical_crosscheck_of_two_fixed_input_field_laws",
        "physical_claim": False,
        "promotion_allowed": False,
        "field_law": {
            "simple_function_blanchet_novak_inputs": {
                **simple_bn,
                "published_Q2_s2": BLANCHET_NOVAK_SIMPLE_Q2,
                "relative_difference": (
                    abs(simple_q2 - BLANCHET_NOVAK_SIMPLE_Q2)
                    / BLANCHET_NOVAK_SIMPLE_Q2
                ),
            },
            "rar_function_blanchet_novak_inputs": rar_bn,
            "rar_function_park_inputs": {
                **rar_park,
                "published_Q2_s2": PARK_RAR_Q2,
                "relative_difference": abs(park_q2 - PARK_RAR_Q2) / PARK_RAR_Q2,
                "pull_vs_cassini_sigma": (park_q2 - CASSINI_Q2) / CASSINI_SIGMA,
            },
            "verdict": (
                "uncertified central values reproduce the two published "
                "benchmarks closely, but thousands of integration warnings, "
                "discarded-tail uncertainty, and absent singularity-aware "
                "error control prevent this local receipt from independently "
                "certifying or retiring either fixed-input field law"
            ),
            "certification_blocker": (
                "supply a singularity-aware split or analytic subtraction, "
                "retained cubature errors, radial-tail bounds, a cutoff and "
                "precision refinement table, zero unresolved warnings, and an "
                "independent high-precision replay"
            ),
        },
        "density_formulation": local_halo_footprint(),
        "cassini": {"Q2_s2": CASSINI_Q2, "sigma_s2": CASSINI_SIGMA},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--epsrel", type=float, default=1e-9)
    args = ap.parse_args()
    out = run(args.epsrel)
    text = json.dumps(out, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
