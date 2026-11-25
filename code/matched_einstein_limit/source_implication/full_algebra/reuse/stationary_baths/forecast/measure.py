"""Calibrate a future radar survey entirely from standard and held observer counts."""

import argparse
import json
from pathlib import Path
from itertools import combinations

import numpy as np

from .. import certificate
from .. import measure as baths
from . import apparatus as app

HERE = Path(__file__).resolve().parent
TRIALS = 2**44
UNIT = 2**34
PAULI = 2 * baths.BASIS
IDENTITY = np.eye(4)
ODDS = baths.ODDS
# Includes all 62*5*12 standards, 15 head settings per source, every target
# range at every delay, the changed flights, and the erasure controls.
TAIL = np.sqrt(np.log(2 * 20000 / 0.01) / (2 * TRIALS))
RANGE_ERROR = 2.5 * (TAIL + 24 / UNIT) + 1e-10
HEAD_ERROR = 2 * (TAIL + 8 / UNIT) + 1e-10


def draw(rng, probabilities):
    p = np.asarray(probabilities)
    if p.min() < -1e-12 or abs(p.sum() - 1) > 2e-11:
        raise ValueError(
            f"terminal probabilities do not form an instrument: min={p.min()}, sum={p.sum()}"
        )
    # Dyadic panel with an explicit total variation allowance in each read.
    tickets = np.rint(np.maximum(p, 0) * UNIT).astype(np.int64)
    tickets[np.argmax(p)] += UNIT - int(tickets.sum())
    return rng.multinomial(TRIALS, tickets / UNIT)


def radar(rng, state, clocks):
    return np.stack([draw(rng, row) for row in app.census(state, clocks)])


def parked(rng, state):
    moments = np.einsum("aij,bji->ab", PAULI, state).real
    weights = np.trace(state, axis1=1, axis2=2).real
    probs = np.stack([(weights + moments) / 2, (weights - moments) / 2], axis=-1)
    counts = np.stack([draw(rng, row.ravel()).reshape(2, 2) for row in probs])
    return counts


def standards():
    out = np.zeros((2, 31, 2, 4, 4), complex)
    for b in range(2):
        out[b, 0, b] = IDENTITY / 4
        for a, p in enumerate(PAULI):
            out[b, 1 + 2 * a, b] = (IDENTITY + p) / 4
            out[b, 2 + 2 * a, b] = (IDENTITY - p) / 4
    return out


def response_counts(rng, channel, clocks):
    panels = []
    for n in app.DELAYS:
        panels.append(
            [
                [radar(rng, app.evolve(channel, state, n), clocks) for state in sector]
                for sector in standards()
            ]
        )
    return np.array(panels)


def decode_calibration(counts):
    ranges = counts @ app.SCORES / TRIALS
    return ranges[:, :, 0], (ranges[:, :, 1::2] - ranges[:, :, 2::2]) / 2


def predict(head, offset, slopes):
    # Range expectation in sector b is its identity offset times its trace,
    # plus fifteen measured signed expectations times calibrated slopes.
    x = (head[:, :, 0] - head[:, :, 1]).T / TRIALS
    # The parity weights are the same fixed preparation controls; obtain them
    # independently from all retained tag outcomes in the first head setting.
    p = head[0].sum(axis=1) / TRIALS
    result = np.einsum("nbr,b->nr", offset, p) + np.einsum("nbar,ba->nr", slopes, x)
    # An opposite probability shift in the two parity masses changes the
    # intercept only by their *difference*. Conditional moments are unneeded.
    radius = RANGE_ERROR * (1 + abs(x).sum()) + (TAIL + 8 / UNIT) * abs(
        offset[:, 1] - offset[:, 0]
    )
    radius += HEAD_ERROR * (abs(slopes) + RANGE_ERROR).sum(axis=(1, 2))
    return result, radius


def length_powers():
    # Each squared edge is a monomial in the twelve positive directed ranges.
    log_ranges = np.zeros((4, 4, 12))
    for i, (a, b) in enumerate(app.ROUTES):
        log_ranges[a, b, i] = 1
    gains = (log_ranges - log_ranges.transpose(1, 0, 2)).mean(axis=1)
    return np.array(
        [
            log_ranges[a, b] + log_ranges[b, a] - gains[a] - gains[b]
            for a, b in combinations(range(4), 2)
        ]
    )


POWERS = length_powers()


def determinant(ranges):
    iv = certificate.c.iv
    log_ranges = [iv.ln(r) for r in ranges]
    squared = [
        iv.exp(sum(float(x) * r for x, r in zip(row, log_ranges, strict=True) if x))
        for row in POWERS
    ]
    a, b, c, e, f, g = squared
    p, q, r = (a + b - e) / 2, (a + c - f) / 2, (b + c - g) / 2
    det = a * b * c + 2 * p * q * r - a * r * r - b * q * q - c * p * p
    dp, dq, dr = 2 * (q * r - c * p), 2 * (p * r - b * q), 2 * (p * q - a * r)
    jac = [
        b * c - r * r + (dp + dq) / 2,
        a * c - q * q + (dp + dr) / 2,
        a * b - p * p + (dq + dr) / 2,
        -dp / 2,
        -dq / 2,
        -dr / 2,
    ]
    return det, [j * d for j, d in zip(jac, squared, strict=True)]


def volume_band(means, radii):
    # Enclose all twelve first derivatives on the range box, then apply the
    # mean-value bound to the nonlinear detector with its gains refitted.
    iv = certificate.c.iv
    radii = np.broadcast_to(radii, means.shape)
    mids = [iv.mpf(float(m)) for m in means]
    widths = [iv.mpf(float(e)) for e in radii]
    ranges = [m + iv.mpf([-1, 1]) * e for m, e in zip(mids, widths, strict=True)]
    det, weights = determinant(ranges)
    if det.a <= 0:
        raise ArithmeticError("unresolved tetrahedron")
    width = iv.mpf(0)
    for r, e, column in zip(ranges, widths, POWERS.T, strict=True):
        deriv = sum(float(k) * w for k, w in zip(column, weights, strict=True) if k)
        width += abs(deriv / r) * e
    center = iv.sqrt(determinant(mids)[0]) / 6
    box = center + iv.mpf([-1, 1]) * width / (12 * iv.sqrt(det))
    return [
        float(np.nextafter(float(box.a), -np.inf)),
        float(np.nextafter(float(box.b), np.inf)),
    ]


def schedules(channels):
    pure = ODDS[:, None, None] * np.diag([1, 0, 0, 0])[None]
    out = {}
    for marked in range(-1, 12):
        weights = np.full(12, 1 / 12 if marked == -1 else 15 / (32 * 11))
        if marked >= 0:
            weights[marked] = 17 / 32
        phi = np.einsum("d,bdij->bij", weights, channels)
        state = app.evolve(phi, pure, 320)
        traces = np.trace(state, axis1=1, axis2=2)
        np.testing.assert_allclose(traces, ODDS, atol=5e-10, rtol=0)
        state = (
            (state + state.conj().transpose(0, 2, 1))
            / 2
            * (ODDS / traces.real)[:, None, None]
        )
        fixed, contraction = baths.prepare(channels, weights)
        label = "uniform" if marked == -1 else f"route_{marked:02d}"
        out[label] = (state, weights, np.linalg.norm(fixed - state), contraction)
    return out


def control_sources(slopes):
    difference = slopes[3, 0] - slopes[3, 1]
    axis, route = np.unravel_index(np.argmax(abs(difference)), difference.shape)
    delta = ODDS[1] * PAULI[axis] / 4
    origin = ODDS[:, None, None] * IDENTITY / 4
    states = {
        f"balanced_{name}": (
            origin + sign * np.stack([delta, -delta]),
            None,
            None,
            None,
        )
        for name, sign in (("plus", 1), ("minus", -1))
    }
    return states, {"axis": baths.AXES[axis], "directed_route_index": int(route)}


def mean_matched(off, slopes):
    frame = slopes[:, :, :, :].reshape(len(app.DELAYS), 30, 12)
    at_four = ODDS @ off[3]
    gradient = baths.station.range_gradient(at_four)[1] @ frame[3].T
    # Traceless Pauli coordinates are the subnormalized moments recorded by
    # the two tag outcomes. Project away the complete initial range frame.
    _, singular, vh = np.linalg.svd(frame[0].T, full_matrices=True)
    if singular[-1] <= singular[0] * 1e-10:
        raise ArithmeticError("calibration range frame unresolved")
    kernel = vh[12:].T
    null = kernel @ (kernel.T @ gradient)
    perturbation = np.einsum("ba,aij->bij", null.reshape(2, 15), PAULI) / 4
    denom = max(
        np.linalg.norm(d, 2) / p for d, p in zip(perturbation, ODDS, strict=True)
    )
    if denom <= 1e-12:
        raise ArithmeticError("calibrated volume gradient unresolved")
    delta = perturbation / (8 * denom)
    base = ODDS[:, None, None] * IDENTITY / 4
    return {
        f"mean_matched_{name}": (base + sign * delta, None, None, None)
        for name, sign in (("plus", 1), ("minus", -1))
    }


def study(channels, clocks, changed):
    rng = np.random.default_rng(1350817)
    uniform = channels.mean(axis=1)
    standards_bank = response_counts(rng, uniform, clocks)
    off, slope = decode_calibration(standards_bank)
    sources = schedules(channels)
    controls, chosen = control_sources(slope)
    sources.update(controls)
    sources.update(mean_matched(off, slope))
    banks = {"standard_return_ticks": standards_bank}
    rows = {}
    origin = ODDS[:, None, None] * IDENTITY / 4
    for label, (state, weights, residual, contraction) in sources.items():
        head = parked(rng, state)
        predicted, error = predict(head, off, slope)
        target = np.stack(
            [radar(rng, app.evolve(uniform, state, n), clocks) for n in app.DELAYS]
        )
        changed_target = radar(rng, app.evolve(changed, state, 4), clocks)
        reset_target = radar(rng, origin, clocks)
        # The erasure reads the same two fixed sector weights; all coin moments vanish.
        reset_prediction = np.einsum("br,b->r", off[0], head[0].sum(axis=1) / TRIALS)
        reset_error = RANGE_ERROR + (TAIL + 8 / UNIT) * abs(off[0, 1] - off[0, 0])
        for key, raw in [
            ("parked", head),
            ("release", target),
            ("contrast_four", changed_target),
            ("reset", reset_target),
        ]:
            banks[label + "_" + key] = raw
        means = target @ app.SCORES / TRIALS
        held_radius = error + RANGE_ERROR
        discrepancy = np.abs(predicted - means)
        if np.any(discrepancy > held_radius):
            raise ArithmeticError(f"uncovered held range: {label}")
        reset_observed = reset_target @ app.SCORES / TRIALS
        if np.any(abs(reset_prediction - reset_observed) > reset_error + RANGE_ERROR):
            raise ArithmeticError("reset forecast failed")
        bands = [volume_band(m, e) for m, e in zip(predicted, error, strict=True)]
        held = [volume_band(m, RANGE_ERROR) for m in means]
        rows[label] = dict(
            route_weights=None if weights is None else weights.tolist(),
            fixed_point_distance=residual,
            conditional_contraction=contraction,
            max_standardized_range_discrepancy=float(np.max(discrepancy / held_radius)),
            largest_directed_range_error=float(discrepancy.max()),
            forecast_volume=[float(baths.station.geometry(m)[1]) for m in predicted],
            forecast_volume_band=bands,
            independent_volume_band=held,
            frozen_volume_excluded=[
                b[1] < held[0][0] or held[0][1] < b[0] for b in held
            ],
            frozen_ranges_excluded_per_delay=(abs(means - means[0]) > 2 * RANGE_ERROR)
            .sum(axis=1)
            .tolist(),
            contrast_four_ranges_excluding_native=int(
                np.sum(
                    abs(changed_target @ app.SCORES / TRIALS - predicted[3])
                    > held_radius[3]
                )
            ),
            reset_control_ranges_excluding_native=int(
                np.sum(abs(reset_observed - predicted[3]) > held_radius[3])
            ),
        )
    return dict(
        delays=app.DELAYS,
        trials_per_setting=TRIALS,
        one_range_half_width=RANGE_ERROR,
        signed_head_half_width=HEAD_ERROR,
        capped_bin_scores=app.SCORES.tolist(),
        balanced_control=chosen,
        held=rows,
    ), banks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    ch, fx = app.instruments(8)
    changed, _ = app.instruments(4)
    out, banks = study(ch, fx, changed.mean(axis=1))
    out["history"] = history(ch, fx)
    if args.output:
        args.output.mkdir(exist_ok=True)
        (args.output / "results.json").write_text(
            json.dumps(out, indent=2, allow_nan=False) + "\n"
        )
        np.savez_compressed(args.output / "counts.npz", **banks)
    else:
        saved = json.loads((HERE / "results.json").read_text())
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(json.loads(json.dumps(out)), saved, "radar forecast")
        with np.load(HERE / "counts.npz") as data:
            if set(data.files) != set(banks):
                raise ValueError("acquisition bank mismatch")
            for k, v in banks.items():
                np.testing.assert_array_equal(data[k], v)
    print(
        "held sources:",
        len(out["held"]),
        "; worst error / joint range allowance:",
        max(v["max_standardized_range_discrepancy"] for v in out["held"].values()),
    )


def history(channels, effects):
    uniform = channels.mean(axis=1)
    initials = standards()
    operators = []
    for n in app.DELAYS:
        mean = np.array(
            [
                [
                    app.census(app.evolve(uniform, s, n), effects) @ app.SCORES
                    for s in sector
                ]
                for sector in initials
            ]
        )
        operators.append(((mean[:, 1::2] - mean[:, 2::2]) / 2).reshape(30, 12).T)
    ans = []
    for k in (1, 2, 3):
        frame = np.vstack(operators[:k])
        sv = np.linalg.svd(frame, compute_uv=False)
        goal = operators[3]
        weights = goal @ np.linalg.pinv(frame, rcond=1e-10)
        leftover = goal - weights @ frame
        ans.append(
            dict(
                delays=app.DELAYS[:k],
                numerical_rank=int(np.sum(sv > sv[0] * 1e-10)),
                extreme_singular_values=[float(sv[0]), float(sv[-1])],
                worst_unexplained_four_cycle_slope=float(
                    np.linalg.norm(leftover, axis=1).max()
                ),
                worst_score_one_norm=float(abs(weights).sum(axis=1).max()),
                worst_count_noise_from_range_inputs=float(
                    RANGE_ERROR * abs(weights).sum(axis=1).max()
                ),
            )
        )
    return ans


if __name__ == "__main__":
    main()
