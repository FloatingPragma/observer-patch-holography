"""Start a local read from sharp stock, and calibrate the same lines separately."""

import argparse
from fractions import Fraction
from hashlib import sha256
from math import acos, asin, ceil, floor, hypot, pi, sqrt
from pathlib import Path
from random import Random
import json

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.stats import beta, binom

from ...fan_path.channel_bound import dyadic_vectors
from ...fan_path.fixture import FanPath, archived_fixture
from ...fan_path.spectral import separated_reference
from ..transfer_series import integer_generator, integer_port_columns, model_fingerprint

HERE = Path(__file__).resolve().parent
DWELLS = (0, 1, 2, 3, 4, 6, 9)
FLIGHT = 27.0625
SHOTS = 2**34
RAMSEY_SHOTS = 2**30
TAIL = 0.01 / 4096
SEED = 832018
PRE = ["first_energy_refused", "local_charge_herald_refused", "second_energy_refused"]
PROBABILITY_TICKS = 2**40
DRAW_PAD = 1e-11


def bound(k, n, pad=0):
    lo = float(beta.ppf(TAIL, k, n - k + 1)) if k else 0
    hi = float(beta.ppf(1 - TAIL, k + 1, n - k)) if k < n else 1
    return [max(0, lo - pad - DRAW_PAD), min(1, hi + pad + DRAW_PAD)]


def minus(a, b):
    return [a[0] - b[1], a[1] - b[0]]


def divide(a, b):
    assert b[0] * b[1] > 0
    endpoints = [x / y for x in a for y in b]
    return [min(endpoints), max(endpoints)]


def draw(probability, shots, rng):
    probability = np.array(probability, dtype=float)
    assert probability.min() > -2e-12
    assert abs(probability.sum() - 1) < 2e-12
    probability = np.maximum(probability, 0)
    probability /= probability.sum()
    ticks = np.rint(PROBABILITY_TICKS * probability).astype(np.int64).tolist()
    largest = max(range(len(ticks)), key=ticks.__getitem__)
    ticks[largest] += PROBABILITY_TICKS - sum(ticks)
    assert min(ticks) >= 0
    assert max(abs(np.asarray(ticks) / PROBABILITY_TICKS - probability)) < DRAW_PAD
    # One deterministic midpoint uniform per conditional bin, including empty
    # bins. A roundoff-scale zero therefore cannot shift later instruments'
    # pseudorandom streams. The last count is fixed by the number enrolled.
    remaining, denominator = shots, PROBABILITY_TICKS
    counts = []
    for mass in ticks[:-1]:
        u = (rng.getrandbits(52) + 0.5) / 2**52
        if mass == 0 or remaining == 0:
            k = 0
        elif mass == denominator:
            k = remaining
        else:
            k = int(binom.ppf(u, remaining, mass / denominator))
        counts.append(k)
        remaining -= k
        denominator -= mass
    return [*counts, remaining]


def reader(columns, receiver, before):
    at = columns[receiver]
    left = columns[before]
    right_port = (at + 1j * left) / sqrt(2)
    left_port = (left + 1j * at) / sqrt(2)
    return (
        at.conj().T @ at,
        right_port.conj().T @ right_port,
        left_port.conj().T @ left_port,
    )


def read(used, available, effects, rng, pad):
    norm = np.trace(used).real
    assert abs(norm - (available[-1])) < 1e-11
    a, b, c = [np.trace(used @ effect).real for effect in effects]
    failures = minus_stages(available)
    charge = draw([*failures, norm - a, a], SHOTS, rng)
    current = draw([*failures, norm - b - c, b, c], SHOTS, rng)
    return {
        "charge_counts": charge,
        "current_counts": current,
        "simulated_charge_and_port_probabilities": [float(x) for x in (a, b, c)],
        "terminal_charge_probability_interval": bound(charge[-1], SHOTS, pad),
        "signed_current_port_difference_interval": minus(
            bound(current[-2], SHOTS, pad), bound(current[-1], SHOTS, pad)
        ),
    }


def minus_stages(available):
    chain = [1, *available]
    return [a - b for a, b in zip(chain, chain[1:])]


def bank(rho, e, available, effects, rng, pad, slow=1):
    records = []
    for t in DWELLS:
        phase = np.exp(-1j * slow * t * e)
        propagated = phase[:, None] * rho * phase.conj()
        records.append({"dwell": t, **read(propagated, available, effects, rng, pad)})
    return records


def contrast_gap(coherent, erased, clock, delay):
    a = {row["dwell"]: row for row in coherent}
    b = {row["dwell"]: row for row in erased}
    baseline = bound(
        b[0]["charge_counts"][-1] + b[delay]["charge_counts"][-1], 2 * SHOTS
    )
    # Use the separately padded original intervals for the reference envelope.
    baseline = [
        min(
            baseline[0],
            *(b[t]["terminal_charge_probability_interval"][0] for t in (0, delay)),
        ),
        max(
            baseline[1],
            *(b[t]["terminal_charge_probability_interval"][1] for t in (0, delay)),
        ),
    ]
    zero, late = [
        minus(a[t]["terminal_charge_probability_interval"], baseline)
        for t in (0, delay)
    ]
    assert zero[0] > 0
    c = divide(late, zero)
    assert -1 < c[0] < c[1] < 1
    return divide([acos(c[1]) / delay, acos(c[0]) / delay], clock)


def quadratures(fringe, weight, rng):
    real_imag = [fringe.real, fringe.imag]
    return [
        draw(
            [
                (1 + weight + 2 * r) / 4,
                (1 + weight - 2 * r) / 4,
                (1 - weight) / 4,
                (1 - weight) / 4,
            ],
            RAMSEY_SHOTS,
            rng,
        )
        for r in real_imag
    ]


def phase_from_counts(records, dwell, pad):
    sign = []
    for counts in records:
        plus = [sum(row[::2]) for row in counts]
        means = 2 * np.array(plus) / RAMSEY_SHOTS - 1
        intervals = [[2 * y - 1 for y in bound(k, RAMSEY_SHOTS, pad)] for k in plus]
        radius = hypot(
            *(
                max(abs(x - y) for y in band)
                for x, band in zip(means, intervals, strict=True)
            )
        )
        center = complex(*means)
        sign.append((center, asin(radius / abs(center))))
    first, last = sign
    phase = -np.angle(last[0] * first[0].conjugate())
    error = last[1] + first[1]
    assert 0 < phase - error < phase + error < pi
    return [(phase - error) / dwell, (phase + error) / dwell]


def formation(path, v, e, stock_rows, rng, vector_pad):
    local = path.field.ring_laplacian(path.field.words[path.heads])
    levels, modes = eigh(local)
    start = (path.vacuum + path.code[:, 0]) / sqrt(2)
    clock_banks = []
    for t in (0, 1):
        evolved = modes @ (np.exp(-1j * t * levels) * (modes.T @ start))
        contrast = 2 * np.conj(path.vacuum @ evolved) * (path.code[:, 0] @ evolved)
        clock_banks.append(quadratures(contrast, 1, rng))
    clock = phase_from_counts(clock_banks, 1, 1e-12)
    # All-vacuum source, with an extra addressed neutral-to-pair feeder off
    # during waiting. A fine filter is applied only in the charged arm.
    weights = (
        abs(v[stock_rows].sum(axis=0) / sqrt(stock_rows.stop - stock_rows.start)) ** 2
    )
    found = []
    for mu in (0, 0.5):
        for j, (energy, weight) in enumerate(zip(e, weights, strict=True)):
            arrays = [
                quadratures(weight * np.exp(-1j * t * (energy + 2 * mu)), weight, rng)
                for t in (0, 1 / 4)
            ]
            observed = phase_from_counts(arrays, 1 / 4, vector_pad)
            lower, upper = divide(observed, clock)
            # The inverse beta tails enter a small-angle division. Report
            # outward to seven decimals, with a numerical endpoint allowance.
            formation_interval = [
                floor((lower - 1e-9) * 1e7) / 1e7,
                ceil((upper + 1e-9) * 1e7) / 1e7,
            ]
            found.append(
                {
                    "per_defect_offset": mu,
                    "line": j,
                    "coherent_arm_weight": float(weight),
                    "quadrature_trial_counts": arrays,
                    "formation_frequency_in_neutral_clock_units": formation_interval,
                }
            )
    return {
        "ring_clock_counts": clock_banks,
        "ring_clock_angular_frequency_interval": clock,
        "formation_pointer_outcome_order": [
            "filter_pass_plus",
            "filter_pass_minus",
            "filter_refused_plus",
            "filter_refused_minus",
        ],
        "formation_lines": found,
    }


def run_path(path, prior, rng):
    h = path.h
    a, _ = integer_generator(path)
    launch, _ = integer_port_columns(path, 0)
    receive, _ = integer_port_columns(path, path.size - 1)
    digest = model_fingerprint(a, launch, receive)
    assert digest == prior["configuration_matrix_and_read_columns_sha256"]
    assert [
        row["count_below"]
        for row in prior["exact_filter_inertia_at_threshold_bracket_ends"]
    ] == [2, 2]
    _, bracket, _ = separated_reference(path)
    spectral_rng = np.random.default_rng(SEED + path.size)
    e, v = eigsh(
        h, k=4, which="SA", v0=spectral_rng.standard_normal(h.shape[0]), tol=1e-13
    )
    e, v = e[:2], v[:, :2]
    _, _, fine_centers, raw_radii, frame = dyadic_vectors(h, e, v, bracket[0], bits=50)
    # Different BLAS libraries can round a trial vector across a dyadic bin.
    # Verify a common outward certificate instead of recording implementation-
    # dependent raw residuals as if they were observable readings.
    centers = [Fraction(round(c * 2**40), 2**40) for c in fine_centers]
    radii = [Fraction(1, 2**36)] * len(centers)
    assert all(
        abs(c - f) + r <= upper
        for c, f, r, upper in zip(centers, fine_centers, raw_radii, radii, strict=True)
    )
    frame_upper = Fraction(1, 2**34)
    assert frame <= frame_upper
    vector_pad = float(frame_upper) + 4 * len(v) / 2**40
    phase_pad = float(max(radii)) * max(DWELLS)
    stock_rows = path.blocks[0, 0].rows
    terminal = path.blocks[path.size - 1, 0].rows
    receiver = np.arange(terminal.start, terminal.stop)
    before = np.array(
        [path.ids[path.bridges[-1][0], path.basis[row][1]] for row in receiver]
    )
    np.testing.assert_array_equal(path.bridge_term[before, receiver].A.ravel(), -4)
    # Calibrate the terminal read as the quarter occupation exchange on the
    # isolated pair: U = exp(-i π/4)(I + i X)/sqrt(2).
    isolated = np.array([[4, -4], [-4, 4]])
    ev, vv = eigh(isolated)
    actual = (vv * np.exp(-1j * pi * ev / 16)) @ vv.T
    np.testing.assert_allclose(
        actual,
        np.exp(-1j * pi / 4) * np.array([[1, 1j], [1j, 1]]) / sqrt(2),
        atol=1e-14,
    )
    effects = reader(v, receiver, before)
    D = effects[0]
    initial_weight = abs(v[stock_rows.start]) ** 2
    initial_total = initial_weight.sum()
    pad = (
        float(
            2 * sqrt(initial_total) * vector_pad
            + vector_pad**2
            + initial_total * (12 * vector_pad + 2 * phase_pad)
        )
        + 2e-12
    )
    first = np.diag(initial_weight)
    post_herald = (D * initial_weight) @ D
    if D[0, 1] < 0:
        assert post_herald[0, 1] < 0
    else:
        assert post_herald[0, 1] > 0
    herald_seen = np.trace(D @ first)
    available = [initial_total, herald_seen, np.trace(post_herald)]
    clock = formation(path, v, e, stock_rows, rng, 4 * vector_pad + 2e-12)
    source = np.zeros(len(v))
    source[stock_rows.start] = 1
    exact_command = -1j * FLIGHT * h
    flight = expm_multiply(
        exact_command, source, traceA=complex(exact_command.diagonal().sum())
    )
    arrived = np.zeros_like(flight)
    arrived[receiver] = flight[receiver]
    sharp_outcome = float(np.vdot(arrived, arrived).real)
    admitted = v.T.conj() @ arrived
    sharp_rho = np.outer(admitted, admitted.conj())
    sharp_available = [1, sharp_outcome, np.trace(sharp_rho).real]
    raw = []
    for t in DWELLS:
        if t:
            command = -1j * t * h
            free = expm_multiply(
                command, arrived, traceA=complex(command.diagonal().sum())
            )
        else:
            free = arrived
        direct_effects = reader(free[:, None], receiver, before)
        raw.append(
            {
                "dwell": t,
                **read(
                    np.eye(1) * sharp_outcome,
                    [1, sharp_outcome, sharp_outcome],
                    [effect / sharp_outcome for effect in direct_effects],
                    rng,
                    1e-10,
                ),
            }
        )
    banks = {
        "sharp_arrival_native_unfiltered": raw,
        "sharp_arrival_same_two_line_filter": bank(
            sharp_rho, e, sharp_available, effects, rng, 8 * vector_pad
        ),
        "stationary_two_lines_herald_omitted": bank(
            first, e, [initial_total] * 3, effects, rng, pad
        ),
        "stationary_two_lines_local_charge_herald": bank(
            post_herald, e, available, effects, rng, pad
        ),
        "herald_then_fine_phase_erasure": bank(
            np.diag(np.diag(post_herald)), e, available, effects, rng, pad
        ),
        "herald_then_half_flight_generator": bank(
            post_herald, e, available, effects, rng, pad, slow=0.5
        ),
    }
    silent = banks["stationary_two_lines_herald_omitted"]
    awake = banks["stationary_two_lines_local_charge_herald"]
    erased = banks["herald_then_fine_phase_erasure"]
    assert all(
        row["signed_current_port_difference_interval"][0]
        <= 0
        <= row["signed_current_port_difference_interval"][1]
        for row in [*silent, *erased]
    )
    measured_gap = contrast_gap(
        awake, erased, clock["ring_clock_angular_frequency_interval"], 3
    )
    measured_slow = contrast_gap(
        banks["herald_then_half_flight_generator"],
        erased,
        clock["ring_clock_angular_frequency_interval"],
        6,
    )
    assert measured_slow[1] < measured_gap[0]
    assert max(measured_slow[0], measured_gap[0] / 2) < min(
        measured_slow[1], measured_gap[1] / 2
    )
    signed = next(
        row["signed_current_port_difference_interval"]
        for row in awake
        if row["dwell"] == 3
    )
    assert signed[1] < 0
    native_formation = [
        row["formation_frequency_in_neutral_clock_units"]
        for row in clock["formation_lines"]
        if row["per_defect_offset"] == 0
    ]
    formation_gap = minus(native_formation[1], native_formation[0])
    assert max(formation_gap[0], measured_gap[0]) < min(
        formation_gap[1], measured_gap[1]
    )
    shifted_formation = [
        row["formation_frequency_in_neutral_clock_units"]
        for row in clock["formation_lines"]
        if row["per_defect_offset"] != 0
    ]
    measured_shifts = [
        minus(after, before)
        for before, after in zip(native_formation, shifted_formation, strict=True)
    ]
    assert all(low > 0 and low <= 1 <= high for low, high in measured_shifts)
    assert measured_gap[1] < min(
        row["formation_frequency_in_neutral_clock_units"][0]
        for row in clock["formation_lines"]
    )
    return {
        "path_fans": path.size,
        "reproduced_filter_and_path_sha256": digest,
        "source_rim_word": path.basis[stock_rows.start],
        "source_and_receiving_hubs": [path.cells[0][0], path.cells[-1][0]],
        "reference_energy_line_centers": [str(x) for x in centers],
        "reference_energy_line_radii": [str(x) for x in radii],
        "stationary_stock_before_herald_admission": float(initial_total),
        "stationary_stock_after_herald_after_filter_probability": float(available[-1]),
        "full_no_phase_source_gateway_arrival_before_local_wait": sharp_outcome,
        "stationary_stock_probability_margin": float(pad),
        "terminal_charge_clock_frequency_interval": measured_gap,
        "half_rate_control_frequency_interval": measured_slow,
        "independent_formation_line_separation_interval": formation_gap,
        "per_defect_control_formation_frequency_shift_intervals": measured_shifts,
        "same_defect_flight_unchanged_by_uniform_defect_energy_offset": True,
        "raw_count_banks": banks,
        "independent_neutral_and_formation_counts": clock,
    }


def run():
    prior = HERE.parent / "beat.json"
    earlier = json.loads(prior.read_text())["paths"]
    rng = Random(SEED)
    with archived_fixture() as archived:
        result = [
            run_path(FanPath(archived, n, -1, asymmetric=True), earlier[n - 3], rng)
            for n in (3, 4)
        ]
    return {
        "source_filter_certificate_sha256": sha256(prior.read_bytes()).hexdigest(),
        "joint_confidence_lower": 0.99,
        "receiver_shots_per_read": SHOTS,
        "quadrature_shots_per_read": RAMSEY_SHOTS,
        "seed": SEED,
        "probability_ticks": PROBABILITY_TICKS,
        "count_probability_rounding_allowance": DRAW_PAD,
        "original_acquisition_sha256": sha256(
            (HERE / "records.json").read_bytes()
        ).hexdigest(),
        "charge_outcome_order": [*PRE, "no_gateway_charge", "terminal_gateway_charge"],
        "current_outcome_order": [
            *PRE,
            "outside_isolated_bridge",
            "gateway_side",
            "entry_side",
        ],
        "paths": result,
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
        from ...replay import compare

        repeat = json.loads(content)
        compare(repeat, json.loads((HERE / "replicate.json").read_text()), "self_start")
        original = json.loads((HERE / "records.json").read_text())
        for first, second in zip(original["paths"], repeat["paths"], strict=True):
            for field in (
                "terminal_charge_clock_frequency_interval",
                "independent_formation_line_separation_interval",
                "half_rate_control_frequency_interval",
            ):
                lo, hi = first[field]
                other_lo, other_hi = second[field]
                assert max(lo, other_lo) < min(hi, other_hi)
            for name, old_records in first["raw_count_banks"].items():
                for a, b in zip(
                    old_records, second["raw_count_banks"][name], strict=True
                ):
                    np.testing.assert_allclose(
                        a["simulated_charge_and_port_probabilities"],
                        b["simulated_charge_and_port_probabilities"],
                        rtol=0,
                        atol=2e-13,
                    )
        print(
            "Reproduced the portable count replicate and matched the original detector laws and clock intervals."
        )


if __name__ == "__main__":
    main()
