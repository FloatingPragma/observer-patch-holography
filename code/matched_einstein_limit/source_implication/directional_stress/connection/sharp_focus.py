"""Ordinary-calibrated transfer from regional marked force to the original sharp images."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response.four_ticks import parse_interval, raw_effect
from ...jacobi_screen import certificate as optical
from ...jacobi_screen import source_response as screen
from .. import bridge
from ..traction.certificate import four_tick_columns
from ..translation import certificate as material
from . import measure

HERE = measure.HERE
SOURCE = HERE.parent.parent
RATIO_CENTER = "3.366"


def archived(name):
    return json.loads((SOURCE / name).read_text())


def focus_tangents(moments, variations, volume_ratio):
    """Differentiate the actual image after combining source tangents, before time-jet algebra."""
    image = screen.measure.matrix_jets(moments)

    def read(delta):
        change = screen.measure.response_jets(moments, delta)
        db, dc = screen.measure.deformation_response(image, change)
        return np.array([db.trace(), -dc.trace()], dtype=object)

    rot = read(variations["rotation"])
    raw = {k: read(variations[k]) for k in ("ordinary", "unpair")}

    # Common volume rotation cancels even when a source changes the rotation angle.
    def scale(panel, scalar):
        return np.vectorize(lambda value: value * scalar, otypes=[object])(panel)

    diff = read(variations["ordinary"] - scale(variations["unpair"], volume_ratio))
    return {
        "ordinary": raw["ordinary"] + scale(rot, screen.intervals()["ordinary"]),
        "pairing": raw["unpair"] + scale(rot, screen.intervals()["unpair"]),
        "difference": diff,
        "difference_at_center": read(
            variations["ordinary"] - scale(variations["unpair"], iv.mpf(RATIO_CENTER))
        ),
        "pairing_before_rotation": raw["unpair"],
        "rotation": rot,
    }


def material_difference(effects, source, rho):
    """Trace-norm bound on d_o F - rho d_p F, with the shared rotation removed."""
    # The centered operator gives an equal expectation on any traceless tangent.
    ident = sum(x.trace() for x in effects) / 8
    centered = [x - ident * sp.eye(4) for x in effects]
    return raw_effect(
        centered, source[1]["ordinary"], sp.Rational("12e-7"), sp.Rational("4e-15")
    ) - rho * raw_effect(
        centered, source[1]["unpair"], sp.Rational("18e-7"), sp.Rational("4e-15")
    )


def ratio(read, rho):
    # d_o = rho d_p + difference on the identical stationary source.
    return 1 / (rho + read["difference"] / read["pairing"])


def encode(read):
    return {k: str(v) for k, v in read.items()}


def common_angle_transfer(geometry, material, rho):
    """Enclose the held prediction using identical angle and angle ratio on both instruments."""

    # All rectangles enclose their endpoints even when the dyadic split rounds.
    def boxes(x, n):
        for i in range(n):
            lo = x.a + (x.b - x.a) * i / n
            hi = x.a + (x.b - x.a) * (i + 1) / n
            yield iv.mpf([lo.a, hi.b])

    values = []
    center = iv.mpf(RATIO_CENTER)
    for alpha in boxes(screen.intervals()["unpair"], 16):
        gp = geometry["pairing_before_rotation"] + alpha * geometry["rotation"]
        mp = material["pairing_before_rotation"] + alpha * material["rotation"]
        for r in boxes(rho, 16):
            gd = (
                geometry["difference_at_center"]
                + (center - r) * geometry["pairing_before_rotation"]
            )
            md = (
                material["difference_at_center"]
                + (center - r) * material["pairing_before_rotation"]
            )
            values.append(1 + (md / mp - gd / gp) / (r + gd / gp))
    return iv.mpf([min(v.a for v in values), max(v.b for v in values)])


def joint_prediction(ordinary, meter_ratios, geometry):
    """Fit two coefficients on ordinary data for two charges, then predict both collar reads."""
    (j0, m0), (j1, m1) = ordinary
    g0, g1 = (g["ordinary"] for g in geometry)
    minor = j0 * m1 - j1 * m0
    assert 0 not in minor
    a = (g0 * m1 - g1 * m0) / minor
    b = (j0 * g1 - j1 * g0) / minor
    predicted = []
    for (j, _), (rj, rm), g in zip(ordinary, meter_ratios, geometry):
        # Ordinary fit gives a*j/g + b*m/g = 1 exactly. Using this identity
        # keeps the large common source rotation out of the held residual.
        fraction = a * j / g["ordinary"]
        held = rm + fraction * (rj - rm)
        residual = g["ratio"] - held
        predicted.append(
            dict(
                current_or_work_fraction=str(fraction),
                predicted_pairing_over_ordinary=str(held),
                observed_minus_predicted_ratio=str(residual),
            )
        )
    return dict(
        ordinary_minor=str(minor),
        current_or_work_gain=str(a),
        marked_force_gain=str(b),
        held=predicted,
    )


def clock_comparison(rows, direction):
    """Derive a common clock acceleration from the two ordinary reads only."""
    force_ratio = rows["-1"][0]["ordinary"] / rows["1"][0]["ordinary"]
    negative, positive = [rows[q][1][direction]["ordinary"] for q in ("-1", "1")]
    den = negative[0] - force_ratio * positive[0]
    if 0 in den:
        return {"ordinary_minor_contains_zero": True}
    u = -(negative[1] - force_ratio * positive[1]) / den
    c = (positive[1] + u * positive[0]) / rows["1"][0]["ordinary"]
    constraints = {}
    for q in ("-1", "1"):
        m, rays = rows[q]
        read = rays[direction]
        correction = m["difference"] / m["pairing"]
        t, g = (
            da - correction * dp for da, dp in zip(read["difference"], read["pairing"])
        )
        entry = {"ordinary_calibrated_compatibility_residual": str(g + u * t)}
        if 0 not in t:
            entry["required_clock_acceleration_with_free_species_gain"] = str(-g / t)
        constraints[q] = entry
    return dict(
        ordinary_nonaffinity=str(u), ordinary_gain=str(c), held_constraints=constraints
    )


def run():
    iv.dps = 32
    atlas = measure.translation.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    src = bridge.action.front.inputs()
    patch = measure.translation.chosen(atlas)
    weights = bridge.native.prior.weights()[[7]]
    rho = parse_interval(
        archived("action_response/four_ticks.json")["ratio_of_volume_adjustment_angles"]
    )
    rays = archived("jacobi_screen/source_response.json")
    assert rays["dwell"] == 4 and (
        rays["launches"],
        rays["opposite_launches"],
    ) == tuple(rings)
    all_read = archived("directional_stress/translation/certificate.json")
    late = archived("directional_stress/connection/late_read.json")
    current = archived("directional_stress/certificate.json")
    work = archived("directional_stress/traction/certificate.json")["charges"]
    charges = {}
    rows = {}
    for q in (-1, 1):
        key = str(q)
        blocks = []
        waves = {"forward": [], "reverse": []}
        for odd in (False, True):
            h, ops = measure.translation.op_table(board, route, patch, q, odd)
            force, unit = ops["marked_force"]
            assert unit == 4
            ids = np.unique(force.nonzero()[0])
            w = four_tick_columns(
                h, ids, bridge.beam.indices(rings[0] + rings[1]), weights[0]
            )
            blocks.append(material.compress(w, force[ids, :][:, ids], unit))
            amp = np.array(optical.amplitude_jets(h, *rings, weights), dtype=object)
            for direction, wave in (
                ("forward", amp),
                ("reverse", amp.transpose(0, 1, 2, 4, 3, 6, 5)),
            ):
                waves[direction].append(optical.probability_jets(wave)[:, :, 0])
        x = {k: parse_interval(v) for k, v in all_read[key]["marked_force"].items()}
        for name, val in late["charges"][key]["reconstructed_full_slopes"].items():
            assert val in x[name]
        x["difference"] = material_difference(blocks, src, rho)
        x["difference_at_center"] = material_difference(
            blocks, src, iv.mpf(RATIO_CENTER)
        )
        pair, turn = bridge.action.front.variable_coefficients(blocks, src)
        error = iv.mpf(["-4e-15", "4e-15"])
        x["pairing_before_rotation"] = pair + error
        x["rotation"] = turn + error
        x["ratio"] = ratio(x, rho)
        optical_reads = {}
        for direction, obs in waves.items():
            m, raw = optical.enclose_moments(
                np.array(obs, dtype=object), src[0], sp.Rational("6e-7"), src
            )
            optical_reads[direction] = focus_tangents(m, raw, rho)
        fwd, rev = optical_reads["forward"], optical_reads["reverse"]
        for direction, sign in (("odd", -1), ("even", 1)):
            optical_reads[direction] = {
                name: fwd[name] + sign * rev[name] for name in fwd
            }
        display = {}
        for direction, read in optical_reads.items():
            g = {k: v[1] for k, v in read.items()}
            g["ratio"] = ratio(g, rho)
            # Compare the centered slopes before adding back the common rho.
            # The terms below equal the normalized ratio minus one.
            force_ratio_minus_one = (
                x["difference"] / x["pairing"] - g["difference"] / g["pairing"]
            ) * g["ratio"]
            direct_ratio = g["ratio"] / x["ratio"]
            force_ratio = 1 + force_ratio_minus_one
            together = common_angle_transfer(g, x, rho)
            lo = max(force_ratio.a, direct_ratio.a, together.a)
            hi = min(force_ratio.b, direct_ratio.b, together.b)
            assert lo <= hi
            force_ratio = iv.mpf([lo, hi])
            g.update(
                ordinary_gain=g["ordinary"] / x["ordinary"],
                force_predicts_pairing=g["ordinary"] * x["ratio"],
                observed_over_force_prediction=force_ratio,
            )
            display[direction] = encode(g)
            read["focus"] = g
        charges[key] = dict(material=encode(x), sharp_focusing=display)
        rows[key] = (x, optical_reads)
    additive = {}
    for direction, kind in (("odd", "energy_current"), ("even", "seam_work")):
        ordinary, ratios, geometric = [], [], []
        for q in ("-1", "1"):
            x, optics = rows[q]
            if direction == "odd":
                primary = parse_interval(current[q]["ordinary"])
                r = parse_interval(current[q]["pairing_over_ordinary"])
            else:
                primary = parse_interval(
                    work[q]["observables"]["full_seam"]["ordinary"]
                )
                r = 1 / parse_interval(work[q]["whole_seam_ordinary_over_pairing"])
            ordinary.append((primary, x["ordinary"]))
            ratios.append((r, x["ratio"]))
            geometric.append(optics[direction]["focus"])
        additive[kind] = joint_prediction(ordinary, ratios, geometric)
    return dict(
        dwell=4,
        route=route,
        rings=rings,
        force_chart_center=int(patch["images"][patch["cells"].index(0)]),
        ordinary_to_pairing_volume_rotation_ratio=str(rho),
        charges=charges,
        two_meter_ordinary_calibrations=additive,
        clock_calibrations={
            k: clock_comparison(rows, k) for k in ("forward", "reverse", "odd", "even")
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "sharp_focus.json").write_text(data)
    else:
        print(data, end="")
