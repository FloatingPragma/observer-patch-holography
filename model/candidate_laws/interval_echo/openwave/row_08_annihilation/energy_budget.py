"""Count the regional energy throughout free evolution, before observing the source."""

import argparse
import json
from itertools import product

import numpy as np
from scipy import sparse

from . import bare_receivers, experiment as ex, walls
from .verify import compare

TIMES = (0, 2, 8, 32)
SHOTS = 2**27
SEED = 81137


def terms(k, cut):
    bits = np.asarray(list(product((0, 1), repeat=10)))
    words = np.empty((len(bits), 20), int)
    words[:, ::2] = 2 * bits
    words[:, 1::2] = 2 - 2 * bits
    _, anchors, _ = walls.geometry()
    layout = []
    for j in range(10):
        group = "ring" if j < 4 else ("left" if j < 7 else "right")
        layout.append((group, 2 * j, 2 * j + 1, True))
    layout += [("ring", j, (j + 1) % 8, False) for j in (1, 3, 5, 7)]
    layout += [
        (group, j, j + 1, False)
        for group, j in (("left", 9), ("left", 11), ("right", 15), ("right", 17))
    ]
    if not cut:
        layout += [
            ("attachments", anchors[0], 8, False),
            ("attachments", anchors[1], 14, False),
        ]
    indices = np.arange(len(bits))
    out = []
    for group, a, b, inner in layout:
        if inner:
            out.append((group, a, b, 4 * k, 1, indices ^ (1 << (9 - a // 2)), True))
        else:
            out.append((group, a, b, 0, 4 * k + 0.5, words[:, a] != words[:, b], False))
    return out


def acquire(rng, k, cut):
    basis = list(product((0, 1), repeat=4))
    ref = walls.reference(0)
    _, ports, _ = walls.geometry()
    q = np.asarray(walls.loads(ref))[ports]
    ground = np.outer(*(walls.wing_ref(k, int(s), cut, 3)[2] for s in q))
    starts = np.zeros((16, 64, 2))
    for wall in (0, 1):
        selected = tuple(z ^ int(wall and j >= 2) for j, z in enumerate(ref))
        starts[basis.index(selected), :, wall] = ground.ravel()
    starts = starts.reshape(-1, 2)
    h = sparse.csr_matrix(walls.operators(k, cut, 3))
    evolved, _ = bare_receivers.flights(h, starts, TIMES)
    crew = terms(k, cut)
    readings = []
    for time, y in zip(TIMES, evolved, strict=True):
        frame = []
        exact = np.zeros(2)
        for group, a, b, low, gap, read, inner in crew:
            if inner:
                bright = np.sum(abs(y - y[read]) ** 2, axis=0) / 4
            else:
                bright = np.sum(abs(y[read]) ** 2, axis=0)
            assert np.all(bright >= 0) and np.all(bright <= 1 + 1e-10)
            draws = rng.binomial(SHOTS, np.minimum(bright.round(12), 1)).tolist()
            frame.append(
                dict(
                    region=group,
                    ports=[a, b],
                    low=low,
                    span=gap,
                    high_counts=draws,
                    probabilities=bright.round(12).tolist(),
                )
            )
            exact += low + gap * bright
        np.testing.assert_allclose(
            exact, np.sum(y.conj() * (h @ y), axis=0), rtol=0, atol=2e-8
        )
        readings.append(dict(time=time, terms=frame))
    return dict(contrast=k, cut=cut, readings=readings)


def summarize(case):
    records = case["readings"]
    m = 2 * sum(len(r["terms"]) for r in records)
    eps = np.sqrt(np.log(2 * m / 0.01) / (2 * SHOTS))
    result = []
    for r in records:
        regions = {}
        for group in ("ring", "left", "right", "attachments", "total"):
            terms = [v for v in r["terms"] if group == "total" or group == v["region"]]
            # Within this setup the initial wing and contact distributions are identical for both sources.
            diff = sum(
                v["span"] * (v["high_counts"][1] - v["high_counts"][0]) / SHOTS
                for v in terms
            )
            err = 2 * eps * sum(v["span"] for v in terms)
            exact = sum(
                v["span"] * (v["probabilities"][1] - v["probabilities"][0])
                for v in terms
            )
            regions[group] = dict(
                measured_surplus=diff,
                simultaneous_halfwidth=err,
                pooled_99_halfwidth=float(
                    np.sqrt(
                        np.log(2 * 5 * len(records) / 0.01)
                        * sum(v["span"] ** 2 for v in terms)
                        / SHOTS
                    )
                ),
                expected_surplus=exact,
            )
        result.append(dict(time=r["time"], regions=regions))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    rng = np.random.default_rng(SEED)
    cases = [acquire(rng, k, cut) for k, cut in ((8, False), (8, True), (0, False))]
    for c in cases:
        c["balance"] = summarize(c)
    result = dict(shots=SHOTS, seed=SEED, cases=cases)
    path = ex.HERE / "energy_budget.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        compare(result, json.loads(path.read_text()), "energy_budget")
        print(
            "Verified the full-ensemble regional collision energy accounts and all controls."
        )


if __name__ == "__main__":
    main()
