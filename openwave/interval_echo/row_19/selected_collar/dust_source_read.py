"""Compare the pressureless tower's binary coincidence source with the full collar read."""

from collections import Counter
from fractions import Fraction
import json
from math import log, sqrt
from os.path import relpath

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

from model.candidate_laws.interval_echo_partition_repair.published_strings import (
    archive,
)
from . import reference_arrival as reference
from model.candidate_laws.interval_echo_published_records.spherical_dust import (
    exterior_reader,
)


def run():
    iv = MPIntervalContext()
    iv.dps = 50
    rng = np.random.default_rng(28122929)
    previous = reference.local.inputs("reference_arrival")
    sharp, fresh, shared = reference.local.survey.standards(
        reference.local.survey.bulk.federation(0)
    )[0][:3]
    arms = reference.programs(sharp, fresh, shared)
    groups = [np.flatnonzero(sharp[0] == x) for x in (32, 48, 16)]
    reps = [g[0] for g in groups]
    classes = np.empty(9, int)
    for i, g in enumerate(groups):
        classes[g] = i
    # Check the grouping on the entire calibrated 81-word support.
    block = sharp[np.ix_(reps, reps)]
    assert np.array_equal(sharp, block[np.ix_(classes, classes)])
    bounds = reference.bounds

    def interactions(p):
        return [
            [iv.log(p[i][0] * p[0][j] / (p[i][j] * p[0][0])) for j in range(3)]
            for i in range(3)
        ]

    def minor(d):
        return d[1][1] * d[2][2] - d[1][2] * d[2][1]

    exact = interactions([[iv.mpf(int(a)) for a in row] for row in block])
    n = previous["batched_starts"]
    error = sqrt(log(2 * 81 / 0.01) / (2 * n)) + 5e-14
    counts = np.array(previous["reference_word_counts"], dtype=object).reshape(9, 9)
    measured = interactions(
        [
            [iv.mpf([counts[i, j] / n - error, counts[i, j] / n + error]) for j in reps]
            for i in reps
        ]
    )
    # At fixed marginal probabilities, constants added to a row make zero
    # contribution to a contrast. Encode each of the two remaining rows by
    # a local randomized bit, with its scale fixed solely on the reference.
    tickets = 2**32
    settings = []
    gains = []
    for row in exact[1:]:
        lo = min(row, key=lambda z: float(z.mid))
        hi = max(row, key=lambda z: float(z.mid))
        gain = hi - lo
        thresholds = []
        for z in row:
            if z is lo:
                threshold = 0
            elif z is hi:
                threshold = tickets
            else:
                pos = (z - lo) * tickets / gain
                threshold = int(pos.a)
                assert threshold == int(pos.b)
            thresholds.append(threshold)
        settings.append(thresholds)
        gains.append(gain)
    bits = [
        x[1]
        for x in reference.local.inputs("pair_clock_information")["local_reader"][
            "calibrations"
        ]["fixed_contact"]["sharp_records_in_calibration_order"]
    ]
    shots = 2**44
    half = sqrt(log(2 * 16 / 0.01) / (2 * shots)) + 1e-14
    results = {}
    for name, table in arms.items():
        assert np.array_equal(table.sum(axis=0), fresh.sum(axis=0))
        assert np.array_equal(table.sum(axis=1), fresh.sum(axis=1))
        sums = [
            sum(
                int(table[i, j])
                for i in range(9)
                for j in range(9)
                if (bits[i], bits[j]) == key
            )
            for key in archive.PAIRS
        ]
        k = []
        for channel, thresholds in enumerate(settings, 1):
            num = sum(
                int(table[i, j]) * thresholds[classes[j]]
                for i in groups[channel]
                for j in range(9)
            )
            k.append(
                int(
                    rng.binomial(
                        shots, float(Fraction(num, int(table.sum()) * tickets))
                    )
                )
            )
        results[name] = {"arrival_pair_counts_per_2592": sums, "coincidence_counts": k}
    for name, row in results.items():
        delta = sum(
            g * (k - f) / shots
            for g, k, f in zip(
                gains,
                row["coincidence_counts"],
                results["fresh"]["coincidence_counts"],
                strict=True,
            )
        )
        # Each rounded conditional probability lies below the exact one by
        # at most one ticket, so the contrast of two preparations inherits
        # the same bound per channel. Two independent count means enter.
        err = sum(gains) * (2 * half + 1 / tickets)
        row["charge_from_two_counts"] = bounds(delta + iv.mpf([-1, 1]) * err)
        row["full_charge"] = previous["arms"][name]["measured_modular_contrast"]
    words = archive.posted_words()
    q = Counter(words)
    scale = iv.log(iv.mpf(q[0, 1]) * q[1, 0] / (q[0, 0] * q[1, 1]))
    changes = []
    for name in ("ordinary", "held"):
        weights = archive.address_lotteries(words)[name]
        changes.append(
            sum(
                -iv.mpf(w.numerator) / w.denominator * iv.log(iv.mpf(q[x]) / len(words))
                for w, x in zip(weights, words, strict=True)
            )
        )
    charge = changes[1] - changes[0]
    # Use the travelling observer's single frequency as ordinary rest mass.
    # Infer the added stock from the source words before comparing the rods.
    single = Fraction(sum(t for t, _ in words), len(words))
    mass = iv.mpf(single.numerator) / single.denominator + charge / scale
    tower_path = exterior_reader.HERE / "results.json"
    tower = json.loads(tower_path.read_bytes())
    projections = []
    for row in tower["ball_surveys"]["held"]:
        den = row["inverse_radius"]
        gain = tower["ordinary_calibration"][str(den)]["gain"]
        expected = [6 * float(mass.mid) * gain] + [0.0] * 9
        errors = [
            abs(a - b)
            for a, b in zip(row["rod_einstein_tensor"], expected, strict=True)
        ]
        projections.append(
            {
                "inverse_radius": den,
                "ordinary_gain": gain,
                "predicted_rest_component": expected[0],
                "maximum_component_difference": max(errors),
            }
        )
    paths = [archive.PUBLISHED, tower_path]
    return {
        "source_classes": [g.tolist() for g in groups],
        "sharp_counts": block.tolist(),
        "sharp_minor": bounds(minor(exact)),
        "observed_minor": bounds(minor(measured)),
        "tickets": tickets,
        "thresholds": settings,
        "gains": [bounds(g) for g in gains],
        "shots": shots,
        "seed": 28122929,
        "arms": results,
        "dust_address_reference_counts": [q[x] for x in archive.PAIRS],
        "dust_log_odds_coefficient": bounds(scale),
        "dust_charge_change": bounds(charge),
        "dust_inferred_change_in_counted_mass": bounds(charge / scale),
        "dust_held_rod_comparison": projections,
        "inputs": reference.coherent.steady.exp.save_inputs(
            "reference_arrival.json",
            "pair_clock_information.json",
            *(relpath(p, reference.local.HERE) for p in paths),
        ),
    }


if __name__ == "__main__":
    result = run()
    reference.coherent.steady.exp.emit("dust_source_read", result)
    print("Measured interaction determinant:", result["observed_minor"])
    for name, arm in result["arms"].items():
        print(name, arm["arrival_pair_counts_per_2592"], arm["charge_from_two_counts"])
