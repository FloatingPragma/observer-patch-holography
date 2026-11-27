"""Single complete Fourier clock on an unknown sharp three-vacancy stock."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh
from scipy.stats import beta, binom

from model.candidate_laws.interval_echo.openwave.row_16_binding.signed_pair.material import (
    Assembly,
)

HERE = Path(__file__).resolve().parent
DELAYS = 65536
PERIOD = 512
WINDOW = 0.25
SHOTS = 2**26
SEED = 26112697
CONTRASTS = (8, 4, 0)
RESOLUTION = PERIOD / DELAYS
TAIL = 0.01 / (4 * DELAYS)


def population(energies):
    """Sum normalized Fourier clock laws for the equal sharp-stock ensemble."""
    bins = np.arange(DELAYS, dtype=np.int64)
    probabilities = np.zeros(DELAYS)
    for energy in energies:
        centre = float(energy) / RESOLUTION
        nearest = round(centre)
        fraction = centre - nearest
        offset = (nearest - bins + DELAYS // 2) % DELAYS - DELAYS // 2
        distance = offset + fraction
        # At each integer offset the squared geometric-series numerator is the same.
        largest = nearest % DELAYS
        distance[largest] = 1.0
        line = (
            np.sin(np.pi * fraction) / (DELAYS * np.sin(np.pi * distance / DELAYS))
        ) ** 2
        line[largest] = (np.sinc(fraction) / np.sinc(fraction / DELAYS)) ** 2
        probabilities += line
    probabilities /= len(energies)
    np.testing.assert_allclose(probabilities.sum(), 1, rtol=0, atol=3e-13)
    return probabilities / probabilities.sum()


def acquire():
    rng = np.random.default_rng(SEED)
    rows = []
    for k in CONTRASTS:
        source = Assembly(0, 3, k).h.toarray()
        if source.shape != (220, 220):
            raise ValueError("trimer stock count changed")
        probs = population(eigvalsh(source))
        counts = rng.multinomial(SHOTS, probs).tolist()
        rows.append(dict(k=k, counts=counts))
    host = rng.multinomial(SHOTS, population([0])).tolist()
    return dict(
        seed=SEED, shots=SHOTS, delays=DELAYS, period=PERIOD, rows=rows, host=host
    )


def load(path=HERE / "levels_counts.json.gz"):
    return json.loads(gzip.decompress(path.read_bytes()))


def save(path, data):
    text = json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n"
    path.write_bytes(gzip.compress(text.encode("ascii"), compresslevel=6, mtime=0))


def survey_bounds(survey, k):
    rows = {r["n"]: r for r in survey["rows"] if r["sign"] == -1 and r["k"] == k}
    if any(
        sum(c for _, _, c in rows[n]["census"]) != size
        for n, size in ((1, 12), (2, 66), (3, 220))
    ):
        raise ValueError("incomplete vacancy census")
    floors = {n: k * min(a for a, b, c in rows[n]["census"]) for n in (1, 2, 3)}
    ceiling = 30 + k * max(a for a, b, c in rows[3]["census"])
    return floors[3], ceiling, min(floors[1] + floors[2], 3 * floors[1])


def decode(data, survey):
    if (
        data["shots"] != SHOTS
        or data["delays"] != DELAYS
        or data["period"] != PERIOD
        or [r["k"] for r in data["rows"]] != list(CONTRASTS)
        or data["host"] != [SHOTS] + [0] * (DELAYS - 1)
    ):
        raise ValueError("incomplete trimer stocks or failed clock standard")
    leak = 1 / (DELAYS * np.sin(np.pi * WINDOW / PERIOD)) ** 2
    # The exact one-sided lower binomial limit is above leak precisely when
    # the corresponding right-tail test rejects at that same alpha.
    critical = int(binom.isf(TAIL, SHOTS, leak)) + 1
    assert (
        binom.sf(critical - 1, SHOTS, leak) < TAIL < binom.sf(critical - 2, SHOTS, leak)
    )
    rows = []
    for r in data["rows"]:
        counts = r["counts"]
        k = r["k"]
        if (
            len(counts) != DELAYS
            or sum(counts) != SHOTS
            or any(type(c) is not int or c < 0 for c in counts)
        ):
            raise ValueError("missing a clock bin or original launch")
        lo, hi, threshold = survey_bounds(survey, k)
        if not (lo >= 0 and hi + WINDOW < PERIOD and threshold + WINDOW < PERIOD):
            raise ValueError("possible clock alias in the claimed threshold interval")
        # Only ordinary windows intersecting the supported energy range have content.
        bins = [
            i
            for i, c in enumerate(counts)
            if c >= critical
            and i * RESOLUTION + WINDOW >= lo
            and i * RESOLUTION - WINDOW <= hi
            and i * RESOLUTION + WINDOW < threshold
        ]
        bins.sort(key=lambda i: (-counts[i], i))
        selected = []
        for i in bins:
            if all(abs(i - j) * RESOLUTION > 2 * WINDOW for j in selected):
                selected.append(i)
        evidence = []
        for i in sorted(selected):
            low = float(beta.ppf(TAIL, counts[i], SHOTS - counts[i] + 1))
            assert low > leak
            evidence.append(
                dict(
                    centre=i * RESOLUTION,
                    window=[i * RESOLUTION - WINDOW, i * RESOLUTION + WINDOW],
                    bright=counts[i],
                    probability_lower=low,
                )
            )
        rows.append(
            dict(
                k=k,
                energy_support=[lo, hi],
                breakup_lower=threshold,
                above_cut=len(bins),
                windows=evidence,
            )
        )
    return dict(leak=leak, min_counts=critical, rows=rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        save(args.acquire, acquire())
    else:
        print(
            json.dumps(
                decode(load(), json.loads((HERE / "counts.json").read_text())), indent=2
            )
        )
