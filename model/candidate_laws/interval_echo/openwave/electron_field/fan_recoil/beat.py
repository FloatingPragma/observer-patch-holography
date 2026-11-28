"""Read signed terminal quadratures on the same two energy-admitted fan lines."""

import argparse
from fractions import Fraction
from math import acos, asin, cos, hypot, pi
from pathlib import Path
import json

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.stats import beta

from ..fan_path.channel_bound import dyadic_vectors
from ..fan_path.fixture import FanPath, archived_fixture
from ..fan_path.spectral import separated_reference, verified_inertia
from .transfer_series import integer_generator, integer_port_columns, model_fingerprint

HERE = Path(__file__).resolve().parent
TIMES = (0, 1, 2, 3, 27.0625, 30.8125)
AXES = ("X", "Y", "Z")
INPUTS = ("X+", "X−", "Y+", "Y−", "Z+", "Z−")
PORTS = (
    "filter_refused",
    "admitted_charge_outside_terminal_fan",
    "terminal_fan_charge_off_gateway",
    "gateway_charge_orbital_code_refused",
    "gateway_charge_and_code_spectator_vacuum_refused",
    "all_accepted_port_plus",
    "all_accepted_port_minus",
)
SHOTS = 2**28
ALPHA = 0.01
MARGINALS = 64
SEED = 983241


def states():
    waves = np.array([[1, 1], [1, -1], [1, 1j], [1, -1j], [1, 0], [0, 1]])
    return waves / np.linalg.norm(waves, axis=1)[:, None]


def resolution(k, n, model_pad=0):
    """Clopper–Pearson marginal with the joint-comparison tail allocation."""
    tail = ALPHA / (2 * MARGINALS)
    lo = beta.ppf(tail, k, n - k + 1) if k else 0
    hi = beta.ppf(1 - tail, k + 1, n - k) if k < n else 1
    return [max(0.0, float(lo) - model_pad), min(1.0, float(hi) + model_pad)]


def subtract(a, b):
    return [a[0] - b[1], a[1] - b[0]]


def divide(a, b):
    assert b[0] * b[1] > 0
    values = [v / w for v in a for w in b]
    return [min(values), max(values)]


def multiply(a, b):
    values = [v * w for v in a for w in b]
    return [min(values), max(values)]


def combine_intervals(a, b):
    return [a[0] + b[0], a[1] + b[1]]


def overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def terminal_effects(path, columns, probes):
    """Measure gateway, orbital and every spectator's vacuum before the Pauli read."""
    receiver = path.size - 1
    cell = columns[path.owner == receiver]
    hub = columns[path.blocks[receiver, 0].rows]
    family = path.decode(columns)
    joined = family.sum(axis=0) / np.sqrt(len(family))
    np.testing.assert_allclose(
        joined, path.source(receiver).T @ columns, rtol=0, atol=3e-16
    )
    local_code = np.einsum("exa,exb->ab", family.conj(), family)
    joint_code = joined.conj().T @ joined
    in_fan = cell.T @ cell
    at_hub = hub.T @ hub
    accepted_rows = probes.conj() @ joined
    pair_port_effects = np.einsum("xa,xb->xab", accepted_rows.conj(), accepted_rows)
    undecoded = np.array(
        [
            np.eye(2) - in_fan,
            in_fan - at_hub,
            at_hub - local_code,
            local_code - joint_code,
        ]
    )
    return np.stack(
        [
            np.concatenate((undecoded, terminal))
            for terminal in pair_port_effects.reshape(3, 2, 2, 2)
        ]
    )


def neutral_clock(path, rng):
    neutral = path.field.ring_laplacian(path.field.words[path.heads])
    e, v = eigh(neutral)
    zero, odd = path.vacuum, path.code[:, 0]
    start = (zero + odd) / np.sqrt(2)
    frequencies = []
    for t in (0, 1):
        state = v @ (np.exp(-1j * t * e) * (v.T @ start))
        contrast = 2 * np.conj(zero @ state) * (odd @ state)
        outcomes = [contrast.real, contrast.imag]
        counts = rng.binomial(
            SHOTS, np.clip((1 + np.asarray(outcomes)) / 2, 0, 1)
        ).tolist()
        estimate = (2 * counts[0] / SHOTS - 1) + 1j * (2 * counts[1] / SHOTS - 1)
        intervals = [[2 * v - 1 for v in resolution(k, SHOTS, 1e-12)] for k in counts]
        radius = hypot(
            max(abs(estimate.real - bound) for bound in intervals[0]),
            max(abs(estimate.imag - bound) for bound in intervals[1]),
        )
        phase_error = asin(radius / abs(estimate))
        frequencies.append((estimate, phase_error, counts))
    initial, last = frequencies
    angle = -np.angle(last[0] * initial[0].conjugate())
    error = initial[1] + last[1]
    # This neutral fan's exact vacuum-to-odd clock has frequencies in [0,2].
    assert angle + error < pi
    return {
        "plus_counts_at_dwell_zero_and_one_for_x_and_y": [
            data[2] for data in frequencies
        ],
        "observed_angular_frequency": float(angle),
        "frequency_interval": [float(angle - error), float(angle + error)],
    }


def acquire_band(path, e, v, rng, pad):
    probes = states()
    source = v.T @ path.source()
    admission = (source @ probes.T).T
    pass_prob = np.sum(abs(admission) ** 2, axis=1)
    effects = terminal_effects(path, v, probes)
    observable = {}
    for coherence in ("admitted_coherent", "admitted_phase_erased"):
        rows = []
        for time in TIMES:
            relative = np.exp(-1j * time * (e - e[0]))
            propagated = admission * relative
            if coherence == "admitted_coherent":
                p = np.einsum("sa,rcab,sb->rsc", propagated.conj(), effects, propagated)
            else:
                # Equivalent to measuring and then erasing the fine line label.
                p = np.einsum("sa,rcaa->rsc", abs(admission) ** 2, effects)
            assert np.max(abs(p.imag)) < 1e-14
            q = np.concatenate(
                (np.repeat((1 - pass_prob)[None, :, None], len(AXES), axis=0), p.real),
                axis=2,
            )
            assert np.min(q) > -1e-12
            np.testing.assert_allclose(q.sum(axis=2), 1, rtol=0, atol=1e-12)
            counts = np.empty(q.shape, dtype=np.int64)
            for axis, law in enumerate(q):
                draws = np.maximum(law.flatten(), 0)
                counts[axis] = rng.multinomial(SHOTS, draws / draws.sum()).reshape(
                    len(probes), len(PORTS)
                )
            total_accepted = int(counts[:, :, -2:].sum())
            total_trials = SHOTS * len(AXES)
            rows.append(
                {
                    "dwell": time,
                    "categorical_counts": counts.tolist(),
                    "joint_accepted_clicks_pooled": total_accepted,
                    "joint_accepted_trials_pooled": total_trials,
                    "joint_probability_interval": resolution(
                        total_accepted, total_trials, pad
                    ),
                }
            )
        observable[coherence] = rows
    return observable


def signed_difference(bank, time, pad):
    data = next(row for row in bank if row["dwell"] == time)
    y_last_input = data["categorical_counts"][AXES.index("Y")][INPUTS.index("Z−")]
    reflected_plus, reflected_minus = y_last_input[-2:]
    launched = sum(y_last_input)
    return {
        "quadrature": "terminal Y+ minus Y−, source Z−, same joint acceptance",
        "counts": [reflected_plus, reflected_minus],
        "trials_at_this_input": launched,
        "probability_difference_interval": subtract(
            resolution(reflected_plus, launched, pad),
            resolution(reflected_minus, launched, pad),
        ),
    }


def infer(records, clock, pad):
    coherent = {row["dwell"]: row for row in records["admitted_coherent"]}
    erased = {row["dwell"]: row for row in records["admitted_phase_erased"]}
    erased_zero_three = [erased[t] for t in (0, 3)]
    baseline = resolution(
        sum(row["joint_accepted_clicks_pooled"] for row in erased_zero_three),
        sum(row["joint_accepted_trials_pooled"] for row in erased_zero_three),
        pad,
    )
    assert all(
        overlap(baseline, row["joint_probability_interval"]) for row in erased.values()
    )
    zeroth, third = [coherent[t]["joint_probability_interval"] for t in (0, 3)]
    a0, a3 = subtract(zeroth, baseline), subtract(third, baseline)
    assert a0[0] * a0[1] > 0
    cosine_interval = divide(a3, a0)
    assert cos(3) < cosine_interval[0] < cosine_interval[1] < 1
    gap = [acos(cosine_interval[1]) / 3, acos(cosine_interval[0]) / 3]
    p0, p3 = [
        row["joint_accepted_clicks_pooled"] / row["joint_accepted_trials_pooled"]
        for row in (coherent[0], coherent[3])
    ]
    erased_rate = sum(
        row["joint_accepted_clicks_pooled"] for row in erased_zero_three
    ) / sum(row["joint_accepted_trials_pooled"] for row in erased_zero_three)
    center = acos((p3 - erased_rate) / (p0 - erased_rate)) / 3
    held = []
    for t in (1, 2):
        prediction = combine_intervals(
            baseline, multiply(a0, [cos(gap[1] * t), cos(gap[0] * t)])
        )
        measured = coherent[t]["joint_probability_interval"]
        assert overlap(prediction, measured)
        held.append(
            {
                "dwell": t,
                "predicted_probability_interval": prediction,
                "observed_probability_interval": measured,
            }
        )
    return {
        "contrast_at_zero_interval": a0,
        "contrast_at_three_interval": a3,
        "measured_gap_in_common_clock_units": center
        / clock["observed_angular_frequency"],
        "gap_interval_in_common_clock_units": divide(gap, clock["frequency_interval"]),
        "clock_phase_interval_uses_dwell_one": True,
        "erased_joint_probability_interval": baseline,
        "held_phase_probability_comparisons": held,
        "held_signed_y_at_dwell_one": {
            "coherent": signed_difference(records["admitted_coherent"], 1, pad),
            "erased": signed_difference(records["admitted_phase_erased"], 1, pad),
        },
    }


def trial(path, rng):
    _, bracket, _ = separated_reference(path)
    certificates = [verified_inertia(path.h, cut) for cut in bracket]
    assert [cert["count_below"] for cert in certificates] == [2, 2]
    e, v = eigsh(
        path.h, k=4, which="SA", v0=rng.standard_normal(path.h.shape[0]), tol=1e-13
    )
    e, v = e[:2], v[:, :2]
    _, _, centers, radii, error = dyadic_vectors(path.h, e, v, bracket[0])
    largest_rounding = Fraction(4 * len(v), 2**40)
    pad = (
        float(8 * (error + largest_rounding) + 4 * Fraction(max(TIMES)) * max(radii))
        + 1e-12
    )
    clock = neutral_clock(path, rng)
    readings = acquire_band(path, e, v, rng, pad)
    comparison = infer(readings, clock, pad)
    # The eigensolve generates the forward process and is excluded from infer().
    exact_gap = centers[-1] - centers[0]
    exact_radius = sum(radii)
    measured = comparison["gap_interval_in_common_clock_units"]
    assert (
        measured[0]
        < float(exact_gap - exact_radius)
        < float(exact_gap + exact_radius)
        < measured[1]
    )
    generator, _ = integer_generator(path)
    source, _ = integer_port_columns(path, 0)
    receiver, _ = integer_port_columns(path, path.size - 1)
    return {
        "fan_count": path.size,
        "charge_sign": path.sign,
        "source_and_receiving_hubs": [path.cells[0][0], path.cells[-1][0]],
        "native_bridge_addresses": path.bridges,
        "configuration_matrix_and_read_columns_sha256": model_fingerprint(
            generator, source, receiver
        ),
        "exact_filter_inertia_at_threshold_bracket_ends": certificates,
        "energy_for_forward_model_only": {
            "centres": list(map(str, centers)),
            "error_radii": list(map(str, radii)),
        },
        "added_model_probability_pad": pad,
        "independent_neutral_clock": clock,
        "count_only_comparison": comparison,
        "count_banks": readings,
    }


def run():
    rng = np.random.default_rng(SEED)
    with archived_fixture() as archive:
        records = [
            trial(FanPath(archive, size, -1, asymmetric=True), rng) for size in (3, 4)
        ]
    first, last = [
        row["count_only_comparison"]["gap_interval_in_common_clock_units"]
        for row in records
    ]
    difference = subtract(first, last)
    assert difference[0] > 0
    return {
        "trials_per_receiver_axis_and_clock_quadrature": SHOTS,
        "source_label_order": INPUTS,
        "receiver_axis_order": AXES,
        "outcome_order": PORTS,
        "launch_and_quadrature_shot_seed": SEED,
        "marginal_interval_budget": MARGINALS,
        "joint_confidence_lower": 1 - ALPHA,
        "three_fan_minus_four_fan_beat_ratio_interval": difference,
        "paths": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as file:
            file.write(content)
    else:
        from ..replay import compare

        compare(
            json.loads(content), json.loads((HERE / "beat.json").read_text()), "beat"
        )
        print("Reproduced the coherent, dephased and clocked beat records.")


if __name__ == "__main__":
    main()
