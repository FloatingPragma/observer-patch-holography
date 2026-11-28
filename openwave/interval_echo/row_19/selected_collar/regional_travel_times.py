"""Predict a regional travel-time variance from word clocks and the measured modular source."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.executable_model import cell, federation
from model.experiments.interval_echo_einstein_station import shared_variance, volume
from . import word_routes as routes

HERE = routes.HERE
SEED = 28120841


def costs(code):
    times = []
    for matching in code:
        graph = cell().adjacency().toarray().astype(bool)
        for a, b in matching:
            graph[a, b] = graph[b, a] = False
        steps = shortest_path(graph, unweighted=True, indices=0)
        times.append(int(2 * steps.sum()))
    return np.array(times)


def law(total, admitted, base, first, second, weights):
    # Count the clock times of all starts directly, including one or zero
    # admissions. No admission label is needed to read the final tape.
    p2 = F(admitted * (admitted - 1), total * (total - 1))
    p1 = F(admitted, total) - p2  # each particular slot alone
    result = {0: 1 - 2 * p1 - p2}

    def put(t, p):
        t = int(t)
        result[t] = result.get(t, F(0)) + p

    for i, j in np.ndindex(9, 9):
        w = F(int(weights[i, j]), 2592)
        put(first[i] - base, p1 * w)
        put(second[j] - base, p1 * w)
        put(first[i] + second[j] - 2 * base, p2 * w)
    assert sum(result.values()) == 1
    return result


def measured(tapes, bins, power):
    # The two clock histograms have known equal population means.
    score = np.array(bins) ** power
    val, err = routes.contrast(tapes)
    return routes.bounds(routes.read_score((val, err), score))


def count_passed(starts, probability, rng):
    batches, tail = divmod(starts, routes.lab.BATCH)
    return int(rng.binomial(routes.lab.BATCH, probability, size=batches).sum()) + int(
        rng.binomial(tail, probability)
    )


def survival_band(failed):
    f = failed / routes.lab.SHOTS
    err = routes.reg.src.probability_error(np.asarray(f), routes.lab.SHOTS)
    return 1 - routes.reg.band(f - err, f + err)


def deadline(tape, bins, rate, source_change, rng, calibration_rng):
    # Draw a deadline for each actual timing outcome. Split large cohorts
    # before the binomial draws to keep their sizes exactly representable.
    readings, observations = [], []
    for arm in tape:
        failed = 0
        for n, t in zip(arm, bins, strict=True):
            p = float(1 - rate ** (t // 2))
            failed += count_passed(n, p, rng)
        readings.append(survival_band(failed))
        observations.append([failed, routes.lab.SHOTS - failed])
    scale = 2 / routes.iv.log(routes.shells.interval(rate))
    elapsed = scale * routes.iv.log(readings[1] / readings[0])
    baseline = readings[0] + routes.iv.mpf(source_change)
    passed = count_passed(observations[1][1], float(rate), calibration_rng)
    extra = survival_band(routes.lab.SHOTS - passed)
    corrected = scale * routes.iv.log(readings[1] / baseline)
    control = scale * routes.iv.log(extra / baseline)
    return {
        "failed_passed": observations,
        "extra_round_trip_failed_passed": [routes.lab.SHOTS - passed, passed],
        "source_corrected_duration": routes.bounds(corrected),
        "source_corrected_extra_round_trip": routes.bounds(control),
        "survival_change": routes.bounds(readings[1] - readings[0]),
        "equivalent_duration_change": routes.bounds(elapsed),
    }


def run():
    iv = routes.iv
    iv.dps = 35
    src = volume.radar_source_and_hand(federation(0))[-1]
    _, fresh, shared = np.rint(2592 * shared_variance.source_words(src)[0]).astype(int)
    d = np.array([0, 1, -1, -1, 0, 1, 1, -1, 0])
    np.testing.assert_array_equal(shared - fresh, 2 * np.outer(d, d))
    code, common = routes.codebook()
    calibration_port = int(np.flatnonzero(common[0])[0])
    assert all((0, calibration_port) not in pattern for pattern in code)
    native, _ = routes.calibrate(code)
    assert native["1"] == native["2"]
    clock = costs(code)
    first = clock[1:]
    full = law(2, 2, clock[0], first, first, fresh)
    quarter_floor = sum(weight * F(1, 4) ** (t // 2) for t, weight in full.items())
    assert quarter_floor == F(3025, 9216) and quarter_floor > F(9, 32)
    rng = np.random.default_rng(SEED)
    timer_rng = np.random.default_rng(SEED + 1)
    calibration_rng = np.random.default_rng(SEED + 2)
    prior = json.loads((HERE / "resolved_shells.json").read_bytes())
    kappa = iv.mpf(prior["conditional_entropy_slope"])
    rows = []
    for z in prior["levels"]:
        N, r = z["cells"], z["radius_hops"]
        n = sum(z["shell_populations"][: r + 1])
        p2 = F(n * (n - 1), N * (N - 1))
        assert p2 == F(z["eligible_pair_probability"])
        for layout, second in (
            ("same", first),
            ("turned", first.reshape(3, 3).T.flatten()),
        ):
            gain = F(2 * int(first @ d) * int(second @ d), 1296)
            exact = p2 * gain
            arms = [law(N, n, clock[0], first, second, w) for w in (fresh, shared)]
            bins = sorted(set(arms[0]) | set(arms[1]))
            delta = [
                hi - lo
                for lo, hi in zip(
                    *[[a.get(t, F(0)) for t in bins] for a in arms], strict=True
                )
            ]
            assert sum(t * x for t, x in zip(bins, delta, strict=True)) == 0
            assert sum(t * t * x for t, x in zip(bins, delta, strict=True)) == exact
            forecast = {}
            for rate in (F(1, 2), F(1, 4)):
                sign = 1 if layout == "same" else -1
                cond = sign * (rate * (1 - rate)) ** 2 / 1296
                assert (
                    sum(x * rate ** (t // 2) for t, x in zip(bins, delta, strict=True))
                    == p2 * cond
                )
                forecast[str(rate)] = {
                    "population": str(p2 * cond),
                    "from_modular_read": routes.bounds(
                        iv.mpf(z["direct_entropy_tangent"])
                        * routes.shells.interval(cond)
                        / kappa
                    ),
                }
            tags = {}
            for q in (1, 2):
                tape = [
                    routes.lab.counts([float(a.get(t, 0)) for t in bins], rng)
                    for a in arms
                ]
                tags[str(q)] = {
                    "counts": tape,
                    "mean_change": measured(tape, bins, 1),
                    "variance_change": measured(tape, bins, 2),
                    "deadlines": {
                        str(rate): deadline(
                            tape,
                            bins,
                            rate,
                            forecast[str(rate)]["from_modular_read"],
                            timer_rng,
                            calibration_rng,
                        )
                        for rate in (F(1, 2), F(1, 4))
                    },
                }
            rows.append(
                {
                    "level": z["level"],
                    "layout": layout,
                    "eligible": n,
                    "bins": list(map(int, bins)),
                    "conditional_variance_gain": str(gain),
                    "deadline_prediction": forecast,
                    "variance_change": str(exact),
                    "from_modular_read": routes.bounds(
                        iv.mpf(z["direct_entropy_tangent"])
                        * routes.shells.interval(gain)
                        / kappa
                    ),
                    "probes": tags,
                }
            )
    return {
        "seed": SEED,
        "quarter_survival_lower_bound_same_code": str(quarter_floor),
        # log(1/4) = 2 log(1/2), and both fresh survival fractions tend to one.
        "limiting_shortening_ratio_half_over_quarter": str(
            2 * (F(1, 2) * (1 - F(1, 2)) / (F(1, 4) * (1 - F(1, 4)))) ** 2
        ),
        "deadline_seed": SEED + 1,
        "round_trip_seed": SEED + 2,
        "starts": routes.lab.SHOTS,
        "calibration_path_in_center": [0, calibration_port, 0],
        "sharp_return_appointments": clock.tolist(),
        "directed_word_score": d.tolist(),
        "levels": rows,
        "inputs": {
            k: sha256((HERE / k).read_bytes()).hexdigest()
            for k in (
                "REGIONAL_TRAVEL_TIMES.md",
                "word_routes.json",
                "resolved_shells.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "regional_travel_times.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for z in data["levels"]:
        print(
            z["level"],
            z["layout"],
            z["from_modular_read"],
            z["probes"]["1"]["variance_change"],
        )
