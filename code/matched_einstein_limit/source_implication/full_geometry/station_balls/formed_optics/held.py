"""Fresh high-count response at both endpoints of the reachable formation segment."""

import json
from math import log
import numpy as np
from mpmath import iv

from . import bounds, certificate, measure

BATCH = 2**48
BATCHES = 4096
STARTS = BATCH * BATCHES
TAIL = log(2 * (2 * 5 * 2 * 5 * 25) / 0.01)


def bands(bank):
    iv.dps = 34
    raw = np.asarray(bank)[:, :, :5].transpose(0, 2, 1)
    value = raw.astype(float)
    width = np.sqrt(2 * value * TAIL) + 2 * TAIL
    lo = np.maximum(0.0, (value - width) / STARTS - 2e-17)
    hi = (value + width) / STARTS + 2e-17
    return np.array(
        [iv.mpf([float(a), float(b)]) for a, b in zip(lo.flat, hi.flat, strict=True)],
        dtype=object,
    ).reshape(raw.shape)


def report():
    data = json.loads((measure.HERE / "certificate.json").read_text())
    rng = np.random.default_rng(56118026)
    result = {}
    for q, views in data["charges"].items():
        rows = {}
        for edge, source in views.items():
            p = np.asarray(source["probabilities_times_1e18"]) / 1e18
            obs = {}
            for key, panels in (("forward", p), ("reverse", p.transpose(0, 2, 1))):
                bank = []
                for table in panels:
                    counts = [
                        rng.multinomial(BATCH, np.r_[a, 1 - a.sum()], size=BATCHES)
                        .sum(axis=0)
                        .tolist()
                        for a in table.T
                    ]
                    bank.append(counts)
                t = bands(bank)
                exact = certificate.decode(source["probabilities_times_1e18"])
                if key == "reverse":
                    exact = exact.transpose(0, 2, 1)
                for z, b in zip(exact.flat, t.flat, strict=True):
                    assert b.a <= z.a and b.b >= z.b
                obs[key] = dict(
                    counts=bank,
                    at_independent=bounds.analyze(t),
                    at_formed=bounds.analyze(t, formed=True),
                )
            rows[edge] = obs
        result[q] = rows
    return dict(
        starts_per_emitter=STARTS,
        batch_size=BATCH,
        batches=BATCHES,
        stocks=[*measure.STOCKS, "formed_ordinary"],
        simultaneous_probability=0.99,
        charges=result,
    )


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
