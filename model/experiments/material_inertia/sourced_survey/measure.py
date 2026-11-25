"""Compare a free incident packet's measured exchange energy with 30 independent clocks."""

import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

from mpmath import iv
import numpy as np

from matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen.range_contraction import (
    held_stations as clocks,
    measure as stationary,
)
from . import source

RECORD = Path(__file__).with_name("results.json")
STARTS = 2**46
SEED = 261129
# 14 clock panels * 30 clocks, 11 incoming panels * 3 probabilities, 14 parked-energy readings, and 3 early position readings.
FAMILY = 500
TAIL = 0.01


def stocks():
    x = stationary.preparations()["reference"]
    right = np.trace(x.sum(axis=0).reshape(2, 2, 2, 2), axis1=0, axis2=2)
    coin = np.kron(np.eye(2) / 2, right)
    blank = np.zeros_like(coin)
    return {
        "reference": x,
        "even": np.stack([coin, blank]),
        "odd": np.stack([blank, coin]),
    }


def acquire_clock(bins, x, rng):
    probs = np.einsum("rtpij,pji->rt", bins, x).real
    assert probs.min() > -1e-15 and probs.sum(axis=1).max() < 1
    probs = np.maximum(probs, 0)
    rates = np.column_stack([probs, 1 - probs.sum(axis=1)])
    return [rng.multinomial(STARTS, r).tolist() for r in rates]


def numbers(x):
    return clocks.geometry.interval_numbers(x)


def volume_ranges(means):
    table = np.full((6, 6), iv.mpf(0), dtype=object)
    for route, value in zip(clocks.ROUTES, means, strict=True):
        table[route] = value
    squared = clocks.geometry.observed_squared(table)
    out = []
    for i, *others in combinations(range(6), 4):
        gram = np.array(
            [[clocks.geometry.dot(squared, i, j, k) for k in others] for j in others],
            dtype=object,
        )
        det = clocks.geometry.curved.det_interval(gram)
        assert det.a > 0
        out.append(iv.sqrt(det) / 6)
    return out


def decode(raw):
    """Uses integer records alone after detector calibration; inclusive for misses and timeouts."""
    iv.dps = 35
    relative = iv.sqrt(iv.log(2 * FAMILY / iv.mpf(".01")) / (2 * STARTS))
    # Fixed-point formation error, timing effect, and binary64 probability contraction.
    score_radius = 2.5 * relative + iv.mpf("1e-9")
    prob_radius = relative + iv.mpf("3e-10")
    assert float(clocks.history.timing_error()[1]) < 1e-12
    cone = iv.mpf([-1, 1])
    out = {}
    raw_means, raw_energies, vols = {}, {}, {}
    for name, bank in raw.items():
        m = []
        for count in bank["round_trips"]:
            assert len(count) == 11 and sum(count) == STARTS and min(count) >= 0
            mid = iv.mpf(sum(i * k for i, k in enumerate(count, 2))) / (4 * STARTS)
            m.append(mid + score_radius * cone)
        raw_means[name] = m
        vols[name] = volume_ranges(m)
        n = bank["parked_odd"]
        assert 0 <= n <= STARTS
        raw_energies[name] = iv.mpf(n) / STARTS + cone * prob_radius
        out[name] = {
            "volumes": [numbers(v) for v in vols[name]],
            "parked_energy": numbers(raw_energies[name]),
        }

    for name, bank in raw.items():
        if "incoming_orientation" not in bank:
            continue
        missed, even, odd = bank["incoming_orientation"]
        assert missed + even + odd == STARTS and min(missed, even, odd) >= 0
        probabilities = [
            iv.mpf(n) / STARTS + prob_radius * cone for n in (missed, even, odd)
        ]
        baseline = raw_energies["reference"]
        added = probabilities[2] - (1 - probabilities[0]) * baseline
        measured = raw_energies[name] - baseline
        assert (measured - added).a <= 0 <= (measured - added).b
        # Coefficient constraints sum to one, so the reference term is exact.
        forecast = [
            r + probabilities[1] * (s - r) + probabilities[2] * (a - r)
            for r, s, a in zip(
                raw_means["reference"], raw_means["even"], raw_means["odd"], strict=True
            )
        ]
        for predicted, observed in zip(forecast, raw_means[name], strict=True):
            residue = predicted - observed
            assert residue.a <= 0 <= residue.b
        delta = [
            numbers(v - b) for v, b in zip(vols[name], vols["reference"], strict=True)
        ]
        out[name].update(
            incident_exchange=numbers(probabilities[2]),
            net_energy_prediction=numbers(added),
            net_energy_observed=numbers(measured),
            max_abs_central_clock_prediction_error=max(
                abs(float((p.mid - q.mid).mid))
                for p, q in zip(forecast, raw_means[name], strict=True)
            ),
            volume_changes=delta,
            contractions=sum(b < 0 for a, b in delta),
            expansions=sum(a > 0 for a, b in delta),
        )
    # The exact weight is determined by native count probabilities, independently of all clocks.
    # Flight energy and its complete distribution are constant throughout the two exposures.
    forward = raw_energies["forward_1"] - raw_energies["forward_4"]
    reverse = raw_energies["reverse_1"] - raw_energies["reverse_4"]
    assert forward.b < 0 < reverse.a
    weight = reverse / (reverse - forward)
    exact = []
    for dwell in [1, 4]:
        m = [
            r + weight * (f - r)
            for f, r in zip(
                raw_means[f"forward_{dwell}"],
                raw_means[f"reverse_{dwell}"],
                strict=True,
            )
        ]
        exact.append(volume_ranges(m))
    signs = [numbers(a - b) for a, b in zip(*exact, strict=True)]
    held = [
        numbers(a - b)
        for a, b in zip(vols["matched_1"], vols["matched_4"], strict=True)
    ]
    out["equal_energy_time_control"] = dict(
        exact_forward_fraction=numbers(weight),
        exact_weight_volume_differences=signs,
        observed_held_volume_differences=held,
        observed_held_energy_difference=numbers(
            raw_energies["matched_1"] - raw_energies["matched_4"]
        ),
    )
    return out


def mixer(raw):
    # Only the incoming energy assay and a blank control enter this hardware coin.
    baseline = Fraction(raw["reference"]["parked_odd"], STARTS)

    def gained(name, dwell):
        miss, _, odd = raw[f"{name}_{dwell}"]["incoming_orientation"]
        return Fraction(odd, STARTS) - (1 - Fraction(miss, STARTS)) * baseline

    diff = {k: gained(k, 1) - gained(k, 4) for k in ("forward", "reverse")}
    assert diff["forward"] < 0 < diff["reverse"]
    weight = diff["reverse"] / (diff["reverse"] - diff["forward"])
    return Fraction(round(weight * 2**24), 2**24)


def native_match(standard, launches, bins):
    """Directed deterministic transfer at exactly matched parked energy."""
    iv.dps = 35
    cone = iv.mpf([-1, 1])
    p = iv.mpf(float(standard["reference"][1].trace().real))
    errors = iv.mpf("1e-10") * cone

    def incoming(name, dwell):
        panel = launches[f"{name}_{dwell}"]
        return [iv.mpf(panel[key]) + errors for key in ("capture", "exchange")]

    def net(name, dwell):
        capture, e = incoming(name, dwell)
        return e - capture * p

    df = net("forward", 1) - net("forward", 4)
    dr = net("reverse", 1) - net("reverse", 4)
    assert df.b < 0 < dr.a
    w = dr / (dr - df)
    # A 4x4 block contributes 16 real products per contraction. Clock bin
    # effect entries are at most 1, state entries are dyadic and at most 1.
    # Binary64 conversion/accumulation over ten bins plus native kernel error
    # is below 1e-10 in each standard capped mean (see proof).
    clock_means = []
    for x in standard.values():
        hits = np.einsum("rtpij,pji->rt", bins, x).real
        mean = 3 + hits @ ((np.arange(2, 12) - 12) / 4)
        clock_means.append([iv.mpf(float(v)) + errors for v in mean])
    grid = []
    for dwell in (1, 4):
        forward, reverse = (incoming(k, dwell) for k in ("forward", "reverse"))
        c, e = [r + w * (f - r) for f, r in zip(forward, reverse, strict=True)]
        m = [
            a + (c - e) * (b - a) + e * (d - a)
            for a, b, d in zip(*clock_means, strict=True)
        ]
        grid.append(volume_ranges(m))
    deltas = [numbers(a - b) for a, b in zip(*grid, strict=True)]
    assert deltas[0][1] < -0.0002
    return dict(forward_weight=numbers(w), volume_differences=deltas)


def run():
    source.local_proof()
    launches, bound = source.flights()
    ports, bins = clocks.return_bins()
    standard = stocks()
    all_states = dict(standard)
    rng = np.random.default_rng(SEED)
    raw = {}
    # Freeze the analyser standards and the independently measured incoming records first.
    for name, x in standard.items():
        raw[name] = {"parked_odd": int(rng.binomial(STARTS, np.trace(x[1]).real))}

    def acquire_source(name, law):
        c, e = law["capture"], law["exchange"]
        weights = np.array([1 - c, c - e, e])
        assert weights.min() >= 0
        raw[name] = {"incoming_orientation": rng.multinomial(STARTS, weights).tolist()}
        all_states[name] = sum(
            w * x for w, x in zip(weights, standard.values(), strict=True)
        )

    for name, law in launches.items():
        acquire_source(name, law)
    coin = mixer(raw)
    for dwell in (1, 4):
        forward, reverse = (launches[f"{k}_{dwell}"] for k in ("forward", "reverse"))
        law = {
            k: float(coin) * v + (1 - float(coin)) * reverse[k]
            for k, v in forward.items()
        }
        name = f"matched_{dwell}"
        launches[name] = law
        acquire_source(name, law)
    for name, x in all_states.items():
        bank = raw[name]
        if "parked_odd" not in bank:
            bank["parked_odd"] = int(rng.binomial(STARTS, np.trace(x[1]).real))
        bank["round_trips"] = acquire_clock(bins, x, rng)
    # The direction meter is a separate, unconditional single-port read.
    early = {
        name: int(rng.binomial(STARTS, p))
        for name, p in bound["port_11_at_one_over_32"].items()
    }
    assert early["forward"] / STARTS > 0.51 and early["reverse"] / STARTS < 0.49
    return {
        "early_port_11_counts": early,
        "exact_equal_energy": native_match(standard, launches, bins),
        "source_ports": [5, 11],
        "capture_ports": [0, 1],
        "clock_ports": list(ports),
        "starts_per_setting": STARTS,
        "risk": TAIL,
        "seed": SEED,
        "quartets_in_station_indices": list(map(list, combinations(range(6), 4))),
        "directed_station_indices": list(map(list, clocks.ROUTES)),
        "return_or_timeout_ticks": list(range(2, 13)),
        "incoming_bins": ["miss", "even", "odd"],
        "source_probabilities": launches,
        "native_enclosure": bound,
        "measured_forward_mixer": [coin.numerator, coin.denominator],
        "counts": raw,
        "inference": decode(raw),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        RECORD.write_text(
            json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        assert result == json.loads(RECORD.read_text())
        print(
            "Replayed incident orientations, parked parity and all 420 return-clock banks."
        )
