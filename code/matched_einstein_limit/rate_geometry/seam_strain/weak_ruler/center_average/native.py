"""Small-ball area averaged over observer centers on the measured five-station star."""

import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from . import cone
from .. import aperture, hinge
from ... import live_bounds, live_modular, live_radar, live_source
from ...source import certificate as cert, measure as source

HERE = Path(__file__).parent
ARMS = ("ordinary", "decorrelated", "left", "right")


def area_coefficient(alpha):
    """3F and its alpha derivative for the mean-volume protocol."""
    x = iv.pi - alpha / 2
    cot = iv.cos(x) / iv.sin(x)
    csc2 = 1 / iv.sin(x) ** 2
    f = x * csc2 - cot
    # df/dalpha = csc²(x) (x cot(x) - 1).
    return (2 * alpha * f / 5, 2 * (f + alpha * csc2 * (x * cot - 1)) / 5)


@lru_cache(maxsize=1)
def regional_energy():
    """Independent native H_R=(PH+HP)/2 and number on the twelve-port cut."""
    board = source.formation.federation(1)
    ports = aperture.enroll(board, 1)[0]
    inside = live_radar.station_rows(board, ports).ravel()
    hi, lo, colerr, _ = live_source.pump_center(True)
    packed = hi.astype(object) * (1 << 36) + lo
    re, im = np.hsplit(packed, 2)
    h = source.system(board, source.formation.region(board), -1)[0]
    part = h[inside].tocoo()
    hr = np.zeros((len(inside), 8), dtype=object)
    hi_part = hr.copy()
    np.add.at(hr, part.row, part.data[:, None] * re[part.col])
    np.add.at(hi_part, part.row, part.data[:, None] * im[part.col])

    def pair(br, bi):
        return (
            sp.Matrix(re[inside].T @ br + im[inside].T @ bi)
            + sp.I * sp.Matrix(re[inside].T @ bi - im[inside].T @ br)
        ) / (1 << 164)

    number = pair(re[inside], im[inside])
    ph = pair(hr, hi_part)
    energy = (ph + ph.H) / 80  # h is the integer 40H.
    norm = sp.Rational(int(abs(h).sum(axis=1).max()), 40)
    (ref, prod), ferr, _ = live_source.formation_center()
    rotations = live_modular.rotation_centers(ref, ferr)
    preparations = dict(
        reference=((2, 1), ref, ferr, 1),
        ordinary=((-1, 1), ref, ferr, 1),
        decorrelated=((2, 1), prod - ref, 3 * ferr, 2),
    )
    for name, (rho, err) in rotations.items():
        preparations[name] = ((2, 1), rho, err, 1)
    out = {}
    for key, operator, bound in (("number", number, 1), ("energy", energy, norm)):
        row = {}
        for arm, (weights, rho, prep_error, trace_bound) in preparations.items():
            src = sp.kronecker_product(sp.diag(*weights) / 3, rho)
            center = sp.re(sp.trace(operator * src))
            err = bound * (
                prep_error + (trace_bound + prep_error) * (6 * colerr + 9 * colerr**2)
            )
            row[arm] = cert.enclosed(center, err)
        out[key] = row
    return dict(ports=ports, **out)


def match_two(f, g, quantity):
    """At fixed first two scalar readings, per unit partial re-pairing.

    Correct the ordinary rail and right marginal-axis rotation. The reference
    is interior in both controls. Positive small re-pairing fractions are physical.
    """
    det = f["ordinary"] * g["right"] - f["right"] * g["ordinary"]
    if det.a <= 0 <= det.b:
        raise ValueError("the two controls do not resolve the constraints")
    rail = (f["right"] * g["decorrelated"] - f["decorrelated"] * g["right"]) / det
    turn = (f["decorrelated"] * g["ordinary"] - f["ordinary"] * g["decorrelated"]) / det
    return dict(
        control_determinant=live_bounds.endpoints(det),
        ordinary_per_repairing=live_bounds.endpoints(rail),
        right_angle_per_repairing=live_bounds.endpoints(turn),
        matched_slopes={
            name: live_bounds.endpoints(
                row["decorrelated"] + rail * row["ordinary"] + turn * row["right"]
            )
            for name, row in quantity.items()
        },
    )


def ordinary_prediction(area, regional):
    """Calibrate the area-loss coefficient on the ordinary rail, then re-pair."""
    ordinary = regional["ordinary"]
    if ordinary.a <= 0 <= ordinary.b:
        raise ValueError("ordinary source contrast does not calibrate this observable")
    gain = area["ordinary"] / ordinary
    predicted = gain * regional["decorrelated"]
    return dict(
        gain=live_bounds.endpoints(gain),
        predicted_repairing=live_bounds.endpoints(predicted),
        measured_minus_predicted=live_bounds.endpoints(area["decorrelated"] - predicted),
    )


def measurements():
    edges, den, err = hinge.exact_ranges()
    inputs = hinge.paired_inputs(
        edges, den, hinge.edge_errors(edges, den, err), hinge.sources()
    )
    lengths = inputs.pop("reference")
    radial = hinge.radial_jets(lengths, inputs)
    _, volume, _ = live_bounds.geometry_jets(lengths, inputs)
    value = iv.mpf(0)
    slopes = {k: iv.mpf(0) for k in ARMS}
    individual_value = iv.mpf(0)
    individual = {k: iv.mpf(0) for k in ARMS}
    rows = []
    for length, angle in radial:
        if not (angle.v.a > iv.pi.b and angle.v.b < (2 * iv.pi).a):
            raise ValueError(
                "averaging protocol requires all enrolled angles between pi and 2pi"
            )
        coef, dcoef = area_coefficient(angle.v)
        endpoints = live_bounds.endpoints(angle.v)
        coeff_i = iv.mpf(list(cone.individual_interval(*endpoints, bins=1024)))
        slope_i = iv.mpf(list(cone.individual_slope_interval(*endpoints)))
        # Enroll the middle third of each radial edge, with equal proper
        # weight per unit observer volume. The total center volume is fixed.
        value += length.v * coef / 3
        individual_value += length.v * coeff_i / 3
        portion = {}
        axis = {}
        for tag in ARMS:
            i = live_bounds.KEYS.index(tag)
            slopes[tag] += (length.g[i] * coef + length.v * dcoef * angle.g[i]) / 3
            individual[tag] += (
                length.g[i] * coeff_i + length.v * slope_i * angle.g[i]
            ) / 3
            portion[tag] = live_bounds.endpoints(angle.g[i])
            # Area response on the axis at the same flat-volume radius ell,
            # divided by ell^2; the parent star and its aperture are distinct.
            axis[tag] = live_bounds.endpoints(
                iv.mpf(2) / 3 * (angle.v / (2 * iv.pi)) ** (-iv.mpf(2) / 3) * angle.g[i]
            )
        rows.append(
            dict(
                length=live_bounds.endpoints(length.v),
                angle=live_bounds.endpoints(angle.v),
                angle_slopes=portion,
                axis_area_slopes_per_radius_squared=axis,
                individually_fixed_volume_area_coefficient=live_bounds.endpoints(
                    coeff_i
                ),
                individually_fixed_volume_angle_derivative=live_bounds.endpoints(
                    slope_i
                ),
                mean_fixed_volume_area_deficit_coefficient=live_bounds.endpoints(coef),
            )
        )
    (region,) = (
        row
        for row in json.loads((hinge.HERE / "regional.json").read_text())["regions"]
        if row["radius"] == 1
    )
    material = regional_energy()
    assert material["ports"] == region["ports"]

    def arms(table):
        return {k: iv.mpf(table[k]) for k in ARMS}

    E = arms(material["energy"])
    N = arms(material["number"])
    K = arms(region["regional_entropy_slopes"])
    V = {k: volume.g[live_bounds.KEYS.index(k)] for k in ARMS}
    # Flat-volume radius ell is common to all preparations. With the center
    # sampling volume W fixed, W times mean area loss is value*ell^4+O(ell^6).
    return dict(
        radial_edges=rows,
        total_mean_volume_area_coefficient=live_bounds.endpoints(value),
        total_mean_volume_area_slopes={
            k: live_bounds.endpoints(x) for k, x in slopes.items()
        },
        total_individually_fixed_volume_area_coefficient=live_bounds.endpoints(
            individual_value
        ),
        individual_area_slopes={
            k: live_bounds.endpoints(x) for k, x in individual.items()
        },
        equal_individual_area_coefficient_and_regional_number=match_two(
            individual, N, dict(energy=E, regional_entropy=K, star_volume=V)
        ),
        ordinary_calibration={
            protocol: {
                name: ordinary_prediction(area, table)
                for name, table in (("regional_entropy", K), ("regional_energy", E))
            }
            for protocol, area in (("individual", individual), ("common", slopes))
        },
        reference_volume=live_bounds.endpoints(volume.v),
        regional=material,
        equal_mean_area_coefficient_and_regional_number=match_two(
            slopes, N, dict(energy=E, regional_entropy=K, star_volume=V)
        ),
        equal_regional_energy_and_number=match_two(
            E,
            N,
            dict(
                mean_area_coefficient=slopes,
                individual_area_coefficient=individual,
                regional_entropy=K,
            ),
        ),
    )


def run():
    saved = iv.dps
    try:
        iv.dps = 35
        return measurements()
    finally:
        iv.dps = saved


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    msg = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "native.json").write_text(msg)
    else:
        print(msg, end="")
