"""Convert a two-sightline pressure difference into a conserved stress bound."""

import argparse
import hashlib
import json
from math import floor, inf, nextafter, sqrt
from pathlib import Path

from mpmath import iv

from cosmology.rar_deep_regime.double_ring_pressure.measure import (
    C_KM_S,
    G,
    KPC_KM,
    critical_density_bounds,
    decimal,
)

HERE = Path(__file__).resolve().parent


def sightline_weight(radius, inner, outer):
    """Length of the line-of-sight interval for which the two radii straddle r."""
    if radius <= inner:
        return 0.0
    first = sqrt((radius - inner) * (radius + inner))
    if radius <= outer:
        return 2 * first
    second = sqrt((radius - outer) * (radius + outer))
    return 2 * (outer - inner) * (outer + inner) / (first + second)


def lower_to_quantum(value, quantum):
    return floor(nextafter(float((value / decimal(quantum)).a), -inf)) * quantum


def bounds(margin, ranges, scale, acceleration_upper):
    inner_lower = 1000 * scale * decimal(ranges[0][0])
    separation_upper = 1000 * scale * (decimal(ranges[1][1]) - decimal(ranges[0][0]))
    anisotropy = margin / (2 * iv.pi * separation_upper)
    integral = 3 * iv.sqrt(3) * iv.pi / 2 * inner_lower**2 * margin
    # Enclose the orbital mass by placing both permitted stellar masses inside R2.
    baryons = decimal(10) ** decimal("10.94") + decimal(10) ** decimal("11.32")
    ring_kpc = scale * decimal(ranges[1][1])
    acceleration = acceleration_upper * decimal(KPC_KM) / 1000
    orbital_upper = baryons + ring_kpc * iv.sqrt(baryons * acceleration / decimal(G))
    speed_lower = iv.sqrt(2 * integral / orbital_upper)
    assert speed_lower < 1 / iv.sqrt(2)
    # A collisionless, spherically symmetric source obeys 2 p_t <= rho v_max^2.
    mass = integral * 2 * (decimal(C_KM_S) / 1000) ** 2
    return {
        "positive_pressure_column_difference_msun_pc2": margin,
        "weighted_anisotropy_density_lower_millionths_msun_pc3": lower_to_quantum(
            anisotropy * 10**6, 1
        ),
        "integrated_positive_anisotropy_lower_msun": lower_to_quantum(integral, 10**7),
        "orbital_mass_inside_second_ring_upper_msun": -lower_to_quantum(
            -orbital_upper, 10**8
        ),
        "necessary_kinetic_speed_ceiling_in_thousandths_of_c": lower_to_quantum(
            speed_lower * 1000, 1
        ),
        "kinetic_mass_lower_msun_if_every_speed_below_1000_km_s": lower_to_quantum(
            mass, 10**12
        ),
    }


def run():
    source_files = [HERE.parent / name for name in ("results.json", "elliptical.json")]
    source, optical = [json.loads(path.read_text()) for path in source_files]
    scale, _ = critical_density_bounds()
    acceleration_upper = decimal(source["a0_box_m_s2"][1])
    aligned, fitted = {}, {}
    for label, report in source["all_positive_exteriors"].items():
        aligned[label] = bounds(
            -report["pressure_contrast_upper"],
            report["ring_ranges_arcsec"],
            scale,
            acceleration_upper,
        )
        pieces = optical["panels"][label]["radius_partitions"]
        partitions = [
            bounds(
                piece["optical_minus_orbital_lower"],
                piece["ring_ranges_arcsec"],
                scale,
                acceleration_upper,
            )
            for piece in pieces
        ]
        # The true angular radii belong to one of these partition boxes.
        fitted[label] = {
            k: (max if k == "orbital_mass_inside_second_ring_upper_msun" else min)(
                row[k] for row in partitions
            )
            for k in partitions[0]
        }
    salpeter = source["nominal"]["Salpeter"]
    return {
        "source_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files
        },
        "spherical_chabrier_force_and_aligned_arcs": aligned,
        "spherical_chabrier_force_and_primary_elliptical_fit": fitted,
        "nominal_salpeter_zero_pressure_secondary_bend_arcsec": salpeter[
            "zero_pressure_secondary_bend_arcsec"
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    record = HERE / "results.json"
    data = run()
    if parser.parse_args().write:
        record.write_text(json.dumps(data, indent=2) + "\n")
    elif json.loads(record.read_text()) != data:
        raise SystemExit("Run does not reproduce the archived stress bounds")
    print("Replayed the stress-conservation bounds for both optical reductions.")
