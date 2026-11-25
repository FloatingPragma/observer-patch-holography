"""Direction-resolved focusing on the stationary survey at one common clock appointment."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ..common_front import certificate as front
from ..action_response import four_ticks
from . import certificate, measure

HERE = measure.HERE


def intervals():
    angles = json.loads((HERE.parent / "null_optics/certificate.json").read_text())[
        "matched_pairing"
    ]
    return {
        name: iv.mpf([int(np.floor(a)), int(np.ceil(b))])
        for name, (a, b) in [
            ("ordinary", angles["ordinary_angle_interval"]),
            ("unpair", angles["angle_interval"]),
        ]
    }


def decode(operators, source):
    m, variations = certificate.enclose_moments(
        operators, source[0], sp.Rational("6e-7"), source
    )
    jac = measure.matrix_jets(m)
    angles = intervals()
    unmixed = {}
    for name, changes in variations.items():
        db, dc = measure.deformation_response(jac, measure.response_jets(m, changes))
        unmixed[name] = np.array([db.trace(), -dc.trace()], dtype=object)

    def matched(name):
        return unmixed[name] + unmixed["rotation"] * angles[name]

    _, b, c = measure.optical_matrices(jac)
    return (
        {s: matched(s) for s in angles},
        unmixed,
        np.array([b.trace(), -c.trace()], dtype=object),
    )


def named(x):
    return {"expansion": str(x[0]), "clock_focusing": str(x[1])}


def reference_calibration(panels, currents):
    """Allow a common source-dependent nonaffinity and test the material conversion independently."""
    parse = four_ticks.parse_interval
    q_minus, q_plus = "-1", "1"
    rays = {q: data["forward_minus_reverse"] for q, data in panels.items()}
    ref = {
        q: (
            parse(row["reference"]["expansion"]),
            parse(row["reference"]["clock_focusing"]),
            parse(currents[q]["reference"]),
        )
        for q, row in rays.items()
    }
    neg, pos = ref[q_minus], ref[q_plus]
    reference_ratio = neg[2] / pos[2]
    u = -(neg[1] - reference_ratio * pos[1]) / (neg[0] - reference_ratio * pos[0])
    c = (pos[1] + u * pos[0]) / pos[2]
    assert c.b < 0
    minor = -neg[0] * pos[2] + neg[2] * pos[0]
    assert 0 not in minor
    tests = {}
    station = json.loads((HERE.parent / "full_algebra/certificate.json").read_text())
    eta = parse(station["calibrated_gain"])
    for name, meter in [("ordinary", "ordinary"), ("unpair", "pairing")]:
        err = {}
        du = {}
        for q, row in rays.items():
            err[q] = (
                c * parse(currents[q][meter])
                - parse(row[name]["clock_focusing"])
                - u * parse(row[name]["expansion"])
            )
            du[q] = err[q] / ref[q][0]
        assert du[q_minus].b < du[q_plus].a
        dc = (neg[0] * err[q_plus] - pos[0] * err[q_minus]) / minor
        area_change = -parse(station["signed_transfers"][meter + "_entropy"]) / eta
        assert area_change.b < 0 and dc.a > 0
        tests[name] = dict(
            fixed_volume_survey_area_derivative=str(area_change),
            required_du_for_fixed_conversion={q: str(v) for q, v in du.items()},
            required_conversion_derivative_if_fitted=str(dc),
        )
    return dict(
        reference_nonaffinity=str(u),
        reference_signed_current_factor=str(c),
        clock_and_material_conversion_minor=str(minor),
        source_tangents=tests,
    )


def combined_responses(panels, work):
    """Absolute reference test and scale-independent bound on four held responses."""
    parse = four_ticks.parse_interval
    data, lower_bounds = {}, {}
    odd, even = {}, {}
    for q, rows in panels.items():
        ref = rows["forward_plus_reverse"]["reference"]
        data[q] = (
            parse(ref["expansion"]),
            parse(ref["clock_focusing"]),
            parse(work[q]["observables"]["full_seam"]["reference"]),
        )
        odd[q] = parse(
            rows["forward_minus_reverse"]["held_focusing_over_current_prediction"]
        )
        even[q] = parse(
            rows["forward_plus_reverse"][
                "held_focusing_over_whole_seam_work_prediction"
            ]
        )

    def minimax_lower(high, low):
        assert high.a > 1 and low.b < 1
        return str(((high.a - low.b) / (high.a + low.b)).a)

    for q in panels:
        lower_bounds[q] = minimax_lower(odd[q], even[q])
    highest = max(odd.values(), key=lambda v: float(v.a))
    lowest = min(even.values(), key=lambda v: float(v.b))
    lower_bounds["both_charges"] = minimax_lower(highest, lowest)
    s0, f0, w0 = data["-1"]
    s1, f1, w1 = data["1"]
    u = -(f0 - f1 * w0 / w1) / (s0 - s1 * w0 / w1)
    gain = (f1 + u * s1) / w1
    assert gain.b < 0 and w0.a > 0 and w1.a > 0 and s0.b < 0 < s1.a
    return dict(
        common_collar_rescaling_minimax_relative_error_lower=lower_bounds,
        absolute_even_if_whole_seam_is_entire_source=dict(
            common_nonaffinity=str(u),
            signed_gain=str(gain),
            positive_gain_thresholds={q: str(-f / s) for q, (s, f, w) in data.items()},
        ),
    )


def run():
    iv.dps = 32
    f = certificate.beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = certificate.beam.screens(board, route)
    coeff = certificate.prior.weights()[[7]]
    source = front.inputs()
    out = {}
    flux_inputs = {}
    energy = json.loads((HERE.parent / "action_response/four_ticks.json").read_text())[
        "charges"
    ]
    work = json.loads(
        (HERE.parent / "directional_stress/traction/certificate.json").read_text()
    )["charges"]
    currents = json.loads(
        (HERE.parent / "directional_stress/certificate.json").read_text()
    )
    rho = four_ticks.parse_interval(
        json.loads((HERE.parent / "action_response/four_ticks.json").read_text())[
            "ratio_of_volume_adjustment_angles"
        ]
    )
    for q in (-1, 1):
        effects = {"forward": [], "reverse": []}
        for odd in (False, True):
            h = certificate.native.generator(board, route, q, odd)
            assert (h - h.T).count_nonzero() == 0
            amps = np.array(certificate.amplitude_jets(h, *rings, coeff), dtype=object)
            for name, wave in [
                ("forward", amps),
                ("reverse", amps.transpose(0, 1, 2, 4, 3, 6, 5)),
            ]:
                effects[name].append(certificate.probability_jets(wave)[:, :, 0])
        differences = {}
        reference = {}
        panels = {}
        for name, tables in effects.items():
            slopes, raw, reference[name] = decode(
                np.array(tables, dtype=object), source
            )
            panels[name] = {s: named(x) for s, x in slopes.items()}
            differences[name] = raw
            if name == "forward":
                e = four_ticks.parse_interval(
                    energy[str(q)]["surveys"]["equal_beam_launches"][
                        "energy_ordinary_over_pairing"
                    ]
                )
                panels[name]["held_focusing_over_seam_energy_prediction"] = str(
                    slopes["unpair"][1] / slopes["ordinary"][1] * e
                )
        angles = intervals()
        for label, sign in [("forward_minus_reverse", -1), ("forward_plus_reverse", 1)]:
            raw = {
                name: val + sign * differences["reverse"][name]
                for name, val in differences["forward"].items()
            }
            matched = {
                name: raw[name] + raw["rotation"] * angles[name] for name in angles
            }
            base = reference["forward"] + sign * reference["reverse"]
            combination = {name: named(value) for name, value in matched.items()}
            combination["reference"] = named(base)
            # Both source tangents use the same surveyed-volume rotation.
            # Subtract rho times the unpairing direction before bounding the ratio.
            correction = raw["ordinary"][1] - rho * raw["unpair"][1]
            ratio = 1 / (rho + correction / matched["unpair"][1])
            combination["pairing_over_ordinary_focusing"] = str(ratio)
            panels[label] = combination
            if sign < 0:
                flux_inputs[q] = matched
                odd_ratio = ratio
            else:
                work_ratio = four_ticks.parse_interval(
                    work[str(q)]["whole_seam_ordinary_over_pairing"]
                )
                combination["held_focusing_over_whole_seam_work_prediction"] = str(
                    ratio * work_ratio
                )
        current_ratio = four_ticks.parse_interval(
            currents[str(q)]["pairing_over_ordinary"]
        )
        panels["forward_minus_reverse"]["held_focusing_over_current_prediction"] = str(
            odd_ratio / current_ratio
        )
        out[str(q)] = panels
    # A source-independent common reparametrization gives F + u theta = c J.
    # Solve the two signed ordinary slopes for u and c, with no collar fit.
    current = {
        q: four_ticks.parse_interval(currents[str(q)]["ordinary"]) for q in (-1, 1)
    }
    current_ratio = current[-1] / current[1]
    minus, plus = [flux_inputs[q]["ordinary"] for q in (-1, 1)]
    clock_u = -(minus[1] - current_ratio * plus[1]) / (
        minus[0] - current_ratio * plus[0]
    )
    clock_c = (plus[1] + clock_u * plus[0]) / current[1]
    assert clock_c.a > 0
    return dict(
        dwell=4,
        launches=rings[0],
        opposite_launches=rings[1],
        species=out,
        ordinary_calibration=dict(
            nonaffinity=str(clock_u), signed_current_factor=str(clock_c)
        ),
        source_dependent_calibration=reference_calibration(out, currents),
        odd_even_comparison=combined_responses(out, work),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "source_response.json").write_text(output)
    else:
        print(output, end="")
