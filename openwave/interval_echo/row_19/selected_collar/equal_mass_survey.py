"""Test the smallest far survey at fixed enclosed ordinary mass and mean volume."""

from math import log, sqrt

from mpmath import iv
import numpy as np

from . import stationary_work as work

SEED = 28122927


def neutral_mix(kinetic, localized, mass_key, area_key):
    gain, loss = (iv.mpf(v[mass_key]) for v in (kinetic, localized))
    assert gain.a > 0 > loss.b
    fraction = gain / (gain - loss)
    area = (1 - fraction) * iv.mpf(kinetic[area_key]) + fraction * iv.mpf(
        localized[area_key]
    )
    return fraction, area


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    data = work.local.inputs("stationary_work")
    board = work.local.survey.bulk.federation(0)
    setup, _ = work.steady.symmetric_instrument(board)
    ports, gates, edges = setup[2], setup[3], setup[6]
    far = gates[1]
    neighbors = {a + b - far for a, b in edges if far in (a, b)}
    region = next(
        r
        for r in work.local.inputs("restored_balls")["surveys"]
        if r["center"] == ports[far] and r["hop_radii"] == [0, 1]
    )
    sites, volumes, areas = (region[k] for k in ("ball_sites", "volumes", "cuts"))
    clock, rate = work.coherent.balls.proper_clock(rng)
    half = (
        sqrt(log(2 * 100000 / 0.01) / (2 * work.N))
        + data["finite_wait_trace_error_upper"]
        + 1e-13
    )
    records = {}
    for reader, direction, states in work.preparations(board, setup):
        old = next(
            r
            for r in data["surveys"]
            if r["instrument"] == reader
            and r["direction"] == direction
            and r["center"] == ports[far]
            and r["radii"] == [0, 1]
        )
        k, d = (old["changes"][key] for key in ("work", "ordinary"))
        tangent, area = neutral_mix(k, d, "ordinary_mass_tangent", "area_tangent")
        entropy = (1 - tangent) * iv.mpf(k["entropy_tangent"]) + tangent * iv.mpf(
            d["entropy_tangent"]
        )
        finite, change = neutral_mix(
            k, d, "held_equal_volume_mass_change", "held_equal_volume_area_change"
        )
        # Choose from the independent calibration. A fresh realization of each
        # protocol includes its own public radius at the registered probability.
        f = float(finite.mid)
        pos = states["fresh"].diagonal()
        baseline = work.survey(pos, sites, neighbors, 0.5, rng)
        branches = []
        for name in ("work", "ordinary"):
            p = 0.75 * pos + 0.25 * states[name].diagonal()
            q = old["changes"][name]["public_probability"]
            branches.append(np.array(work.survey_probabilities(p, sites, neighbors, q)))
        # Record the branch label as well as every radius, clock and load read.
        joint = work.coherent.draws(
            np.concatenate(((1 - f) * branches[0], f * branches[1])), rng
        )
        counts = [a + b for a, b in zip(joint[:8], joint[8:], strict=True)]
        dv, da, dm = work.held_differences(counts, baseline, volumes, areas, rate, half)
        records[reader + "_" + direction] = {
            "center": ports[far],
            "radii": [0, 1],
            "mass_neutral_tangent_localization_fraction": work.bounds(tangent),
            "mass_neutral_area_tangent": work.bounds(area),
            "mass_neutral_entropy_tangent": work.bounds(entropy),
            "matched_finite_localization_fraction": work.bounds(finite),
            "set_fraction": f,
            "predicted_matched_quarter_area_change": work.bounds(change),
            "fresh_counts": baseline,
            "held_branch_counts": joint,
            "held_volume_change": work.bounds(dv),
            "held_mass_change": work.bounds(dm),
            "held_area_change": work.bounds(da),
        }
    return {
        "seed": SEED,
        "starts": work.N,
        "far_site": ports[far],
        "clocks": clock,
        "clock_ratio": work.bounds(rate),
        "records": records,
        "inputs": work.steady.exp.save_inputs(
            "EQUAL_MASS_SURVEY.md", "stationary_work.json"
        ),
    }


if __name__ == "__main__":
    result = run()
    work.steady.exp.emit("equal_mass_survey", result)
    for arm, row in result["records"].items():
        print(
            arm,
            row["held_volume_change"],
            row["held_mass_change"],
            row["held_area_change"],
        )
