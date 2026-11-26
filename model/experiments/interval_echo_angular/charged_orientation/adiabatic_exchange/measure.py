"""Record guiding failures, winding and exchange on graded native vacancy sources."""

import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from . import certificate, exchange

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "results.json.gz"
SHOTS = 2**22
GRID = 2**22
SEED = 70261
EPS = float(np.sqrt(np.log(2 * 1_000_000 / 0.01) / (2 * SHOTS)) + 56 / GRID)


def bank(prob, rng):
    p = np.asarray(prob, dtype=float)
    if abs(p.sum() - 1) > 5e-9 or p.min() < -1e-14:
        raise ValueError("Invalid outcome law")
    p = np.maximum(p, 0) / p.sum()
    base = np.rint(p.ravel() * GRID).astype(np.int64)
    base[np.argmax(base)] += GRID - base.sum()
    if np.any(base < 0):
        raise ValueError("Invalid finite sampler grid")
    counts = rng.multinomial(SHOTS, base / GRID).reshape(p.shape)
    return dict(probability=p.tolist(), counts=counts.tolist())


def panel(n, kind, duration, h, rng):
    words, snapshots, exits, incoherent = exchange.trace(n, kind, duration, h)
    timing = []
    for leg, half, probs in snapshots:
        arms = []
        for axis, p in enumerate(probs):
            sign = 1 if axis == 0 else -1
            if leg:
                a, b, park = (exchange.ONE if n == 1 else exchange.PAIR)[leg - 1]
                tracked = (a, b) if half == 1 else (b,)
                choices = [
                    tuple(
                        sorted(
                            [(sign * x) % 8]
                            + ([] if park is None else [(sign * park) % 8])
                        )
                    )
                    for x in tracked
                ]
            else:
                choices = [exchange.source(n)]
            which = [words.index(word) for word in choices]
            raw = bank(p, rng)
            arms.append(
                dict(
                    **raw,
                    tracked_words=[list(x) for x in choices],
                    localized_probability=float(np.sum(p[which])),
                    localized_counts=sum(raw["counts"][i] for i in which),
                )
            )
        timing.append(dict(leg=leg, half=half, arms=arms))
    fringe = bank(exits, rng)
    return dict(
        holes=n,
        law=kind,
        step_duration=duration,
        appointment=h,
        stages=timing,
        fringe=fringe,
        incoherent=bank(incoherent, rng),
        dark_probability=float(exits[1].sum()),
        dark_counts=sum(fringe["counts"][1]),
    )


def run():
    rng = np.random.default_rng(SEED)
    h = exchange.LAW.seam_term(8)
    from_col, to_col = [exchange.LAW.SRC.index(w) for w in ((0, 1), (1, 0))]
    clocks = []
    for t in (0.25, 0.125):
        p = abs(expm(-1j * t * h)[to_col, from_col]) ** 2
        raw = bank([1 - p, p], rng)
        q = raw["counts"][1] / SHOTS
        clocks.append(
            dict(
                **raw,
                setting=t,
                reading=float(2 * np.arcsin(np.sqrt(q))),
                interval=[
                    float(2 * np.arcsin(np.sqrt(val)))
                    for val in (max(0, q - EPS), min(1, q + EPS))
                ],
            )
        )
    panels = []
    for dur, step in [(256, 0.25), (1024, 0.25), (4096, 0.25), (4096, 0.125)]:
        for n in (1, 2):
            for kind in (
                "native",
                "exchange_minus",
                "flux_pi",
                "zero_contrast",
                "parked",
            ):
                panels.append(panel(n, kind, dur, step, rng))
    return dict(
        shots=SHOTS,
        seed=SEED,
        probability_grid=GRID,
        source="PCG64 multinomial aggregate draws from the calculated finite pulse law; independent cohorts at each snapshot and interferometer setting",
        radians="native exchange",
        simultaneous_probability_radius=EPS,
        geometry=exchange.geometry(),
        stock_paths=exchange.initial_stocks(),
        certificate=certificate.local_certificate(),
        reflection=certificate.reflection_sectors(),
        clock_calibration=clocks,
        word_orders={
            str(n): [list(w) for w in exchange.generator(n, 8, "native")[0]]
            for n in (1, 2)
        },
        panels=panels,
    )


def verify(a, b):
    if isinstance(b, dict):
        if a.keys() != b.keys():
            raise ValueError("Different receipt fields")
        for k in b:
            verify(a[k], b[k])
    elif isinstance(b, list):
        if len(a) != len(b):
            raise ValueError("Different receipt sizes")
        if b and isinstance(b[0], float):
            if not np.allclose(a, b, atol=8e-9, rtol=1e-8):
                raise ValueError("Different calculated probability")
        else:
            for x, y in zip(a, b, strict=True):
                verify(x, y)
    elif isinstance(b, float):
        if not np.isclose(a, b, atol=8e-9, rtol=1e-8):
            raise ValueError("Different observable")
    elif a != b:
        raise ValueError("Different setting or count")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    target = ap.parse_args().output
    data = run()
    if target:
        with target.open("xb") as f:
            f.write(
                gzip.compress(
                    json.dumps(data, separators=(",", ":"), allow_nan=False).encode(),
                    mtime=0,
                )
            )
    else:
        verify(data, json.loads(gzip.decompress(RECEIPT.read_bytes())))
        print("Adiabatic exchange and winding receipts verified")
