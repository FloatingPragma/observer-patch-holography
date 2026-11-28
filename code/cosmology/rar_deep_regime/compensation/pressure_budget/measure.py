"""Static spherical pressure support and the observed lensing pseudo-mass budget."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from .. import measure as aperture

HERE = Path(__file__).parent
C = 299792.458  # km/s


def layers(implied, cap):
    return np.ceil(implied / cap).astype(int).tolist()


def compare(case, radius):
    total = np.array(case["total_enclosed_mass_msun"])
    ordinary = np.array(case["subtracted_baryon_cap_msun"])
    epsilon = eta = 0.001
    # Treat the extracted lens mass as the linear, spherical lensing proxy.
    # Give the tiny potential-gradient and compactness corrections their
    # largest admitted signs before subtracting the allowed baryon mass.
    factor = (1 + epsilon) / (1 - eta)
    ratio = total / (ordinary * factor) - 1
    caps = {
        "arbitrary_matching": 1,
        "native_external_coherent": 3 / 7,
        "native_external_sharp_quantum": 17 / 80,
        "sharp_classical": 17 / 82,
    }
    excess = total / factor - ordinary * (1 + caps["native_external_coherent"])
    # At R=r an imposed OUTER traction term 2πr R² p_r(R)/c²
    # has units of mass; quote its equivalent mean pressure only at R=r.
    masskg = (1 + epsilon) * excess * aperture.MSUN_KG
    radius_m = radius * 3.085677581491367e19
    return dict(
        supported_lens_factor=factor,
        additional_mass_over_baryons=ratio.tolist(),
        closed_stock_layer_equivalents={k: layers(ratio, v) for k, v in caps.items()},
        minimum_additional_mass_msun=excess.tolist(),
        missing_mass_equivalent_pa_at_R_equal_r=(
            masskg * (1000 * C) ** 2 / (2 * np.pi * radius_m**3)
        ).tolist(),
        extra_25_percent_downward_optical_mass_dark_ratio=(
            0.75 * total / (ordinary * factor) - 1
        ).tolist(),
    )


def profile():
    # Finite-mass anisotropic profile with flat-space conservation. G=c=1.
    # The common amplitude cancels from all ratios and can be arbitrarily small.
    def rho(r):
        return 1 / (1 + r * r) ** 2

    mass_inf = np.pi**2
    out = []
    for r in (0.5, 1.0, 2.0, 4.0, 10.0):
        m = 2 * np.pi * (np.arctan(r) - r / (1 + r * r))
        for w in (0.0, 0.5, 1.0):
            pressure = w * rho(r)

            def tangential(s):
                return w * rho(s) * (1 - s * s) / (1 + s * s)

            tail = (
                -4
                * np.pi
                * r
                * quad(lambda s: s * tangential(s), r, np.inf, epsabs=1e-12)[0]
            )
            support = 2 * np.pi * r**3 * pressure
            np.testing.assert_allclose(support, tail, atol=3e-12)
            assert support <= mass_inf - m + 1e-12
            out.append(
                dict(
                    r=r,
                    w=w,
                    mass_inside=m / mass_inf,
                    lens_inside=(m + support) / mass_inf,
                    traction_mass=support / mass_inf,
                    exterior_mass=1 - m / mass_inf,
                )
            )
    return out


def outer_wall():
    # An imposed wall at R supports a uniform isotropic pressure up to the edge.
    # At leading weak field, r=R, so the tail is empty and only traction remains.
    r = 1.0
    rho = 1.0
    mass = 4 * np.pi * r**3 * rho / 3
    rows = []
    for w in (0.0, 0.5, 1.0):
        traction = 2 * np.pi * r**3 * w * rho
        lens = mass + traction
        rows.append(
            dict(
                pressure_fraction=w,
                lens_over_mass=lens / mass,
                boundary_over_mass=traction / mass,
            )
        )
    return rows


def paid_trails(cases, radius_kpc):
    """Transfer the three actual conserved board excesses, using the same ratio per stack."""
    from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock import (
        trail,
    )

    boards = json.loads((trail.HERE / "trail.json").read_text())["boards"]
    result = []
    for board in boards:
        initial = Fraction(board["initial_energy"])
        paid = Fraction(board["net_work"])
        assert Fraction(board["terminal_matched_blank"]) == initial
        ratio = paid / initial
        rows = {}
        for label in (
            "nominal",
            "two_systematic_double_baryons",
            "bonferroni_four_normal_95_double_baryons",
        ):
            case = cases[label]
            baryons = np.array(case["subtracted_baryon_cap_msun"])
            observed = np.array(case["total_enclosed_mass_msun"])
            maximum = (1.001 / 0.999) * (1 + float(ratio)) * baryons
            to_speed = aperture.G / (radius_kpc / 1000)
            rows[label] = dict(
                maximum_equivalent_speed_kms=np.sqrt(maximum * to_speed).tolist(),
                compared_speed_kms=np.sqrt(observed * to_speed).tolist(),
                can_supply_optical_lower_mass=(maximum >= observed).tolist(),
            )
        result.append(
            dict(
                live_seams=board["live_seams"],
                net_work=str(paid),
                incident=str(initial),
                conserved_ratio=str(ratio),
                cases=rows,
            )
        )
    return result


def run():
    inputs = json.loads((aperture.HERE / "inputs.json").read_text())
    row = inputs["mistele_all"]["outer"]
    cases = aperture.mass_budget(row)["cases"]
    return dict(
        radius_kpc=row["radius_kpc"],
        mass_allowance="twice each catalog bin ceiling",
        potential_gradient_bound=0.001,
        compactness_bound=0.001,
        caps=[str(Fraction(3, 7)), str(Fraction(17, 80)), str(Fraction(17, 82))],
        cases={
            k: compare(cases[k], row["radius_kpc"])
            for k in (
                "two_systematic_double_baryons",
                "bonferroni_four_normal_95_double_baryons",
            )
        },
        commissioned_trails=paid_trails(cases, row["radius_kpc"]),
        smooth_stress_checks=profile(),
        externally_confined_checks=outer_wall(),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "results.json"
    if p.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.measure import (
            check,
        )

        check(json.loads(path.read_text()), result)
