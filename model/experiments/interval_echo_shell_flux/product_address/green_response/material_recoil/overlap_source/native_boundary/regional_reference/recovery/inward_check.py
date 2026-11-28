"""Held overlap preparation with an inward recorded shift on two apparatuses."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import joint_records as joint

source = joint.source
HERE = joint.HERE
SHOTS = 2**30
SEED = 116807


def run():
    rng = np.random.default_rng(SEED)
    states, b, d, other = source.preparation()
    target = 0.68
    w = (4 * target - 1) / 3
    v = np.sqrt(w) * d + 1j * np.sqrt(1 - w) * other
    held = np.outer(b, b) / 4 + 0.75 * np.outer(v, v.conj())
    u = float((held[1, 1] + held[2, 2]).real)
    count = int(rng.binomial(SHOTS, u))
    # Union of 3*7 equipped expectations, two weighted bare-radius differences,
    # and one binary source read. Simulated counts on independent preparations.
    eps = np.sqrt(np.log(2 * 24 / 0.01) / (2 * SHOTS))
    _, dist, ports, blocks = source.matched.calibrated()
    bare_counts = {}
    bare_means = {}
    for k, rho in (("paired", states["P"]), ("held", held)):
        bare_counts[k] = {}
        radius = {}
        for p in (None, *ports):
            probs = source.matched.trace(rho, blocks[1, p]["number"][2])
            counts = rng.binomial(SHOTS, np.clip(probs, 0, 1))
            radius[p] = counts @ dist / SHOTS
            bare_counts[k][str(p)] = counts.tolist()
        bare_means[k] = np.array([radius[p] - radius[None] for p in ports])
    bare = bare_means["held"] - bare_means["paired"]
    bare_width = 2 * np.linalg.norm(dist) * eps
    stocks, _ = joint.inputs()
    empty = source.timing.get_arms()[0]["empty"]
    stocks = {k: stocks[k] for k in ("baseline", "unit")}
    stocks["held"] = [empty[0] / 2, (empty[1] + held) / 2]
    laws = joint.joint_laws(8, stocks)
    counts = {}
    observed = {}
    for k, law in laws.items():
        banks = [
            rng.multinomial(
                SHOTS, np.maximum(h, 0).ravel() / np.maximum(h, 0).sum()
            ).reshape(h.shape)
            for h in law
        ]
        observed[k] = joint.decode(banks)
        counts[k] = [h.tolist() for h in banks]
    data = {}
    for k in ("unit", "held"):
        delta = observed[k] - observed["baseline"]
        score = (delta[0] + 10 * delta[2]) / 11
        data[k] = dict(
            witness=float(score),
            witness_half_width=float(42 * eps / 11),
            connected_radius_shifts_hops=(delta[1:3] * (2 * max(dist))).tolist(),
            radius_half_width_hops=float(4 * eps * 2 * max(dist)),
        )
    return dict(
        protocol_sha256=sha256((HERE / "INWARD_CHECK.md").read_bytes()).hexdigest(),
        seed=SEED,
        trials=SHOTS,
        source_agreement_count=count,
        source_agreement_interval=[count / SHOTS - eps, count / SHOTS + eps],
        bare_connected_radius_hops=bare.tolist(),
        bare_half_width=float(bare_width),
        equipped=data,
        bare_counts=bare_counts,
        equipped_counts=counts,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "inward_check.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        source.region.check(json.loads(path.read_text()), data)
    print(
        data["source_agreement_interval"],
        data["bare_connected_radius_hops"],
        data["bare_half_width"],
        data["equipped"],
    )
