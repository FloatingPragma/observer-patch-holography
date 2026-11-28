"""Read the native repair timetable and test its equal-exposure alternatives."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components

from .. import measure as stock
from . import poisson_precision as precision

weak = precision.experiment
lab = precision.lab
HERE = precision.HERE
TICKS = 512
SEED = 918472
ARMS = ("recorded", "even", "grouped", "static_mixture")


def acquire(placement, seed, ticks, groups, edges, cumulative):
    rng = np.random.default_rng(seed)
    w = np.array(list(map(int, placement)) + [0, 0, 0])
    readings = []
    for _ in range(ticks):
        for stage in groups:
            for edge in stage:
                a, b = edges[edge]
                code = 9 * w[4 + edge] + 3 * w[a] + w[b]
                for table in cumulative:
                    code = np.searchsorted(table[code], rng.integers(64), side="right")
                w[a], w[b], w[4 + edge] = code % 9 // 3, code % 3, code // 9
        readings.append(int(stock.MAG[w[5]]))
    return readings


def source_clock():
    W, V, _, _ = stock.kernels()
    _, _, groups, edges, powers, ix, words, codes = stock.carrier()
    graph = stock.graph(W, edges, powers, ix, codes) + stock.graph(
        V, edges, powers, ix, codes
    )
    _, labels = connected_components(graph, directed=False)
    starts = ix[(words[:, 4:] == 0).all(axis=1)]
    classes = {}
    for label in np.unique(labels[starts]):
        members = labels == label
        flags = stock.MAG[words[members, 4:]] == 2
        placements = starts[labels[starts] == label]
        key = "".join(map(str, words[placements[0], :4]))
        classes[key] = dict(
            placements=["".join(map(str, words[x, :4])) for x in placements],
            words=int(members.sum()),
            seam_fractions=[
                str(Fraction(int(k), len(flags))) for k in flags.sum(axis=0)
            ],
        )
    assert sum(len(x["placements"]) for x in classes.values()) == 81
    cumulative = [np.cumsum(m, axis=1) for m in (W, V)]
    return classes, {
        k: acquire(k, SEED, TICKS, groups, edges, cumulative) for k in ("2200", "1100")
    }


def time_record(unitaries, start, bits):
    wave = start.copy()
    ages = []
    for bit in bits:
        wave = unitaries[int(bit)] @ wave
        ages.append(wave)
    waves = np.column_stack(ages)
    rho = waves @ waves.conj().T / len(bits)
    np.testing.assert_allclose(np.trace(rho), 1, atol=2e-11, rtol=0)
    return rho


def measure(rho, meters, rng):
    p = rho.diagonal().real
    return dict(
        sites=rng.multinomial(weak.READS, p / p.sum()).tolist(),
        fringes=[lab.fringe_counts(c, rho, rng) for _, _, c in meters],
    )


def decode(banks, meters, tape, arms=ARMS):
    out = {}
    gain = lab.ordinary_gain(tape)
    for name in arms:
        data = {}
        for k, (r, w, c) in enumerate(meters):
            norm = np.trace(abs(c)) + np.triu(abs(c), 1).sum()
            spans = []
            for tag in ("blank", name):
                pos = banks[tag]["sites"]
                fr = banks[tag]["fringes"][k]
                assert min(pos) >= 0 and len(pos) == len(w) and sum(pos) == weak.READS
                assert len(fr) == 3 and min(fr) >= 0 and sum(fr) == lab.READS
                spans.append((w @ pos / weak.READS, (fr[0] - fr[1]) / lab.READS))
            a, b = np.subtract(spans[1], spans[0])
            # Each bounded read uses the common family-wise risk allocation.
            signs = np.array([-1, 1])
            site = a + signs * np.ptp(w) * np.sqrt(2 * lab.lab.RISK / weak.READS)
            fringe = norm * (b + signs * np.sqrt(8 * lab.lab.RISK / lab.READS))
            data[str(r)] = (
                (lab.lab.signed_product(gain, site) + fringe) / (len(w) - 1)
            ).tolist()
        out[name] = data
    return out


def release_contrast(rho, blank, meters):
    """Exact native initial acceleration difference at each enrolled boundary."""
    delta = (rho - blank).real
    return {
        str(radius): float(
            (4 * lab.lab.K * (weights @ delta.diagonal()) + (fringe * delta).sum())
            / (len(weights) - 1)
        )
        for radius, weights, fringe in meters
    }


def apparatus(tape):
    b, keep, marks, y, h0, _, hw, _, _ = weak.fixture(tape)
    H, i, j = precision.flight.detector(b, keep)
    meters = precision.meters(b, keep, marks, y, h0, H, i, j)["written"]
    assert [row[0] for row in meters] == [2, 3, 4, 6]
    return h0, hw, meters


def run(tape):
    classes, readings = source_clock()
    assert 2 not in readings["1100"]
    flags = np.array(readings["2200"]) == 2
    total = int(flags.sum())
    ages = np.arange(TICKS + 1)
    sequences = dict(
        recorded=flags,
        even=np.diff(ages * total // TICKS),
        grouped=np.arange(TICKS) < total,
    )
    h0, hw, meters = apparatus(tape)
    psi = weak.filtered(*eigh(h0))[0]
    unit = [weak.unitary(h, weak.STEP) for h in (h0, hw)]
    blank = time_record(unit, psi, np.zeros(TICKS, dtype=int))
    p = total / TICKS
    density = {"blank": blank}
    timing = {}
    for name, bits in sequences.items():
        assert int(sum(bits)) == total
        delta = weak.STEP * (np.r_[0, np.cumsum(bits)] - p * ages)
        density[name] = time_record(unit, psi, bits)
        timing[name] = float(np.max(abs(delta)))
    density["static_mixture"] = (1 - p) * blank + p * time_record(
        unit, psi, np.ones(TICKS, dtype=int)
    )
    rng = np.random.default_rng(SEED + 1)
    banks = {tag: measure(rho, meters, rng) for tag, rho in density.items()}
    decoded = decode(banks, meters, tape)
    native = {tag: release_contrast(density[tag], blank, meters) for tag in ARMS}
    for tag, cuts in native.items():
        for r, value in cuts.items():
            lo, hi = decoded[tag][r]
            assert lo <= value <= hi
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CONTACT_TIMING.md",
                HERE / "mediated_release.json",
                stock.HERE / "results.json",
            )
        },
        classes=classes,
        seed=SEED,
        readings=readings,
        ticks=TICKS,
        tick_width=weak.STEP,
        active_ticks=total,
        maximum_time_imbalance=timing,
        source_banks=banks,
        decoded=decoded,
        native_acceleration_change=native,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    modes = cli.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "contact_timing.json"
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    if args.counts_only:
        result = json.loads(path.read_bytes())
        assert (
            decode(result["source_banks"], apparatus(tape)[2], tape)
            == result["decoded"]
        )
    else:
        result = run(tape)
        if args.write:
            path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), result)
    print(result["active_ticks"], result["maximum_time_imbalance"], result["decoded"])
