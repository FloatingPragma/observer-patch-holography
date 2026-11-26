"""Read the prepared retained source with independent parked Pauli and radar banks."""

import argparse
from itertools import combinations
import json
import platform
from math import log, sqrt
from pathlib import Path

import mpmath as mp
import numpy as np

from .. import control, held
from ...stationary_baths.forecast import apparatus as app, measure as forecast
from .... import station

HERE = Path(__file__).resolve().parent
STEP = 1e-5
POW2 = 66
N = 1 << POW2
Q = 1 << 36
BATCH = 1 << 51
SEED = 362741903
# 2 endpoints × (30 signed block moments + 2 block masses + 12 clock means).
TAIL = sqrt(log(2 * 88 / 0.01) / (2 * N))
SETTLE = 2 * (1 - (control.TWIRL - STEP)) ** 256
NUM = 1e-10
HEAD_ERROR = 2 * (TAIL + 4 / Q + SETTLE) + NUM
CLOCK_ERROR = 2.5 * (TAIL + 12 / Q + SETTLE) + NUM
EDGES = list(combinations(range(4), 2))


def tickets(p):
    p = np.asarray(p)
    assert p.min() >= -2e-13 and abs(p.sum() - 1) < 1e-11
    num = np.rint(np.maximum(p, 0) * Q).astype(np.int64)
    num[np.argmax(num)] += Q - int(num.sum())
    return num.tolist()


def prepare():
    x, ch, _ = control.compare()
    source_at, clocks, meta = held.panel(ch, x)
    assert meta["phase_axes"] == [2, 0]
    out = []
    for t in (-STEP, STEP):
        source = source_at(t)
        assert (
            np.linalg.eigvalsh(source).min() > 0
            and abs(np.trace(source, axis1=1, axis2=2).sum() - 1) < 1e-11
        )
        moments = np.einsum("aij,bji->ab", forecast.PAULI, source).real
        weights = source.trace(axis1=1, axis2=2).real
        head = np.stack(
            [(weights + moments) / 2, (weights - moments) / 2], axis=-1
        ).reshape(15, 4)
        out.append(
            dict(
                pauli=[tickets(row) for row in head],
                clock=[tickets(row) for row in app.census(source, clocks)],
            )
        )
    return out


def counts(tables):
    rng = np.random.Generator(np.random.PCG64(SEED))

    def draw(row):
        # Every binomial call sees an integer <=2^51. Sum in 128-batch chunks to
        # keep subtotals representable as int64 and the final counts as Python ints.
        p = np.asarray(row, dtype=float) / Q
        out = [0] * len(row)
        for _ in range(N // (128 * BATCH)):
            batch = rng.multinomial(BATCH, p, size=128).sum(axis=0)
            out = [a + int(b) for a, b in zip(out, batch, strict=True)]
        assert sum(out) == N
        return [str(v) for v in out]

    return [
        {key: [draw(row) for row in arm[key]] for key in ("pauli", "clock")}
        for arm in tables
    ]


def area_band(means):
    iv = mp.iv
    mid = [iv.mpf(float(v)) for v in means]
    e = iv.mpf(float(CLOCK_ERROR))
    whole = [r + iv.mpf([-1, 1]) * e for r in mid]

    def faces(ranges):
        lg = [iv.ln(r) for r in ranges]
        sq = [
            iv.exp(sum(float(k) * x for k, x in zip(row, lg, strict=True) if k))
            for row in forecast.POWERS
        ]
        grad = [iv.mpf(0) for _ in sq]
        area = iv.mpf(0)
        for vertices in combinations(range(4), 3):
            ind = [EDGES.index(edge) for edge in combinations(vertices, 2)]
            a, b, c = [sq[i] for i in ind]
            root = iv.sqrt(2 * (a * b + b * c + c * a) - a * a - b * b - c * c)
            area += root / 4
            for i, a, b, c in zip(ind, (a, b, c), (b, c, a), (c, a, b), strict=True):
                grad[i] += a * (b + c - a) / (4 * root)
        return area, grad

    center = faces(mid)[0]
    grad = faces(whole)[1]
    rad = sum(
        abs(sum(float(k) * g for k, g in zip(col, grad, strict=True) if k) / r) * e
        for col, r in zip(forecast.POWERS.T, whole, strict=True)
    )
    return bounds(center + iv.mpf([-1, 1]) * rad)


def bounds(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def summarize(arm):
    h = np.array([[int(v) / N for v in row] for row in arm["pauli"]]).reshape(15, 2, 2)
    w = h[0].sum(axis=1)
    blocks = (
        w[:, None, None] * np.eye(4)
        + np.einsum("ab,aij->bij", h[:, :, 0] - h[:, :, 1], forecast.PAULI)
    ) / 4
    eigen = np.linalg.eigvalsh(blocks).ravel()
    # Pauli orthogonality gives Frobenius <=2*HEAD_ERROR on each block.
    lo = eigen - 2 * HEAD_ERROR
    hi = eigen + 2 * HEAD_ERROR
    assert np.min(lo) > 0
    entropy = -sum(eigen * np.log(eigen))
    slope = np.maximum(abs(1 + np.log(lo)), abs(1 + np.log(hi)))
    delta = float(slope.sum() * (2 * HEAD_ERROR) + 1e-13)
    c = np.array([[int(v) / N for v in row] for row in arm["clock"]])
    means = c @ app.SCORES
    area, volume = station.geometry(means)
    return dict(
        entropy=dict(value=float(entropy), band=[entropy - delta, entropy + delta]),
        area=dict(value=float(area), band=area_band(means)),
        volume=dict(value=float(volume), band=forecast.volume_band(means, CLOCK_ERROR)),
    )


def decode(samples):
    minus, plus = [summarize(a) for a in samples]
    results = {}
    iv = mp.iv

    def interval(x):
        return iv.mpf(x)

    for key in ("entropy", "area", "volume"):
        m, p = minus[key], plus[key]
        results[key] = dict(
            value=(p["value"] - m["value"]) / (2 * STEP),
            band=bounds(
                (interval(p["band"]) - interval(m["band"])) / iv.mpf("0.00002")
            ),
        )
    for name, scale in [("ordinary_residual", 1), ("double_gain_control", 2)]:
        # Published ordinary-gain enclosure, also used for its preselected doubled control.
        eta = iv.mpf(["23.10423617469", "23.10423621780"]) * scale
        results[name] = dict(
            value=results["entropy"]["value"]
            + held.ORDINARY_GAIN * scale * results["area"]["value"],
            band=bounds(
                interval(results["entropy"]["band"])
                + eta * interval(results["area"]["band"])
            ),
        )
    return dict(endpoint_summaries=[minus, plus], secant=results)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--replay-seed", action="store_true")
    args = p.parse_args()
    path = HERE / "records.json"
    if args.replay_seed:
        reference = json.loads((HERE / "certificate.json").read_text())["sampler"]
        versions = dict(
            python=platform.python_version(),
            numpy=np.__version__,
            system=platform.system(),
            machine=platform.machine(),
            bit_generator="PCG64",
        )
        if versions != reference:
            p.exit(
                2,
                "The registered sampler differs on this platform. Use audit.py for platform-independent verification.\n",
            )
        archive = json.loads(path.read_text())
        assert counts(archive["tickets"]) == archive["counts"]
        print("Replayed all integer counts on the registered sampler.")
        return
    model = prepare()
    if args.write:
        samples = counts(model)
        result = dict(
            trials_per_setting=str(N),
            trial_chunks=str(N // BATCH),
            tickets_per_setting=Q,
            seed=SEED,
            signed_twirl_offsets=[-STEP, STEP],
            tickets=model,
            counts=samples,
            **decode(samples),
        )
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        archive = json.loads(path.read_text())
        for arm, ref in zip(model, archive["tickets"], strict=True):
            for key in ("pauli", "clock"):
                assert np.max(np.abs(np.array(arm[key]) - ref[key])) <= 1
        for key, section in decode(archive["counts"])["secant"].items():
            np.testing.assert_allclose(
                section["band"], archive["secant"][key]["band"], atol=2e-7
            )
        print("Reproduced the independent paired head and arrival records.")


if __name__ == "__main__":
    main()
