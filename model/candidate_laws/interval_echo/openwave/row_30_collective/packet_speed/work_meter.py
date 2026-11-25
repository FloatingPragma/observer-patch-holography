"""A separate local energy instrument for the initial velocity of a bare wall."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import measure
from .. import packets


def run():
    observed = json.loads((measure.HERE / "results.json").read_text())
    rng = np.random.default_rng(8376245)
    launches = 2**32
    raw, result = {}, []
    for row in observed:
        if row["contrast"] != 8:
            continue
        n, idx = row["dimers"], row["kick_index"]
        x = np.arange(n + 1)
        shape = np.exp(-((x - n / 2) ** 2) / (4 * min(4.0, (n + 1) / 7) ** 2))
        counts = []
        for shift in (np.pi / 2, -np.pi / 2):
            f = packets.synthesize(
                shape
                * np.exp(1j * x * (measure.MOMENTA[idx] + shift))
                / np.linalg.norm(shape)
            )
            mean_x = 2 * (f[:-1].conj() * f[1:]).real
            counts.append(rng.binomial(launches, (1 + mean_x) / 2))
        up, down = counts
        raw[f"n{n}_p{idx}"] = np.array(counts)
        # E=-h sum X_j + constant, h=1/2 in dimer-exchange-gap units.
        derivative = (down - up).sum() / (2 * launches)
        # Independent setting/site Bernoulli counts. Hoeffding for their weighted sum.
        radius = np.sqrt(n * np.log(2 * 14 / 0.01) / (4 * launches))
        band = [float(derivative - radius), float(derivative + radius)]
        flights = [r for r in observed if r["dimers"] == n and r["kick_index"] == idx]
        result.append(
            {
                "dimers": n,
                "kick_index": idx,
                "work_derivative": float(derivative),
                "work_band": band,
                "disjoint_from_flight": {
                    str(int(r["contrast"])): bool(
                        band[1] < r["velocity_band"][0]
                        or r["velocity_band"][1] < band[0]
                    )
                    for r in flights
                },
            }
        )
    return result, raw


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    out = parser.parse_args().output
    result, raw = run()
    if out:
        out.mkdir(exist_ok=False)
        np.savez_compressed(out / "work_counts.npz", **raw)
        (out / "work.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        expected = json.loads((measure.HERE / "work.json").read_text())
        assert len(result) == len(expected)
        for a, b in zip(result, expected, strict=True):
            np.testing.assert_allclose(a["work_band"], b["work_band"], atol=1e-14)
            assert a["disjoint_from_flight"] == b["disjoint_from_flight"]
        with np.load(measure.HERE / "work_counts.npz") as archive:
            assert set(raw) == set(archive.files)
            for key in raw:
                np.testing.assert_array_equal(archive[key], raw[key])
    for row in result:
        print(
            row["dimers"],
            row["kick_index"],
            row["work_band"],
            row["disjoint_from_flight"],
        )
