"""Use ordinary interventions and one local source read to forecast two instruments."""

import argparse
import json
import numpy as np
from . import measure as source

HERE = source.HERE


def run():
    stock = source.matched.stocks()
    baseline = stock["target"][0][1]
    standards = [a[1] for key in ("exchange_1_1", "balanced_sharp") for a in stock[key]]
    eye = np.eye(4)
    bright = (eye[1] + eye[2]) / np.sqrt(2)
    dark = (eye[1] - eye[2]) / np.sqrt(2)
    other = (eye[0] - eye[3]) / np.sqrt(2)
    states = [baseline, *standards]
    for theta in (np.pi / 6, np.pi / 3):
        vector = np.cos(theta) * dark + 1j * np.sin(theta) * other
        states.append(
            0.25 * np.outer(bright, bright) + 0.75 * np.outer(vector, vector.conj())
        )
    board, dist, origins, blocks = source.matched.calibrated()
    order = [
        "baseline",
        "orbital_low",
        "orbital_high",
        "endpoint_low",
        "endpoint_high",
        "pi_6",
        "pi_3",
    ]
    keep = np.equal.outer([1, 0, 0, -1], [1, 0, 0, -1])
    physical = {
        name: [np.zeros((2, 2)), (r * keep).real]
        for name, r in zip(source.timing.NAMES, states, strict=True)
    }
    banks, _ = source.timing.bank(physical)
    rng = np.random.default_rng(191228)
    shots = 2**24
    times = np.array([0.5, 1, 2, 4, 6, 8, 9])
    n = len(states) * (2 * len(times) + len(blocks[1, None]["number"])) + 2
    eps = np.sqrt(np.log(2 * n / 0.01) / (2 * shots))
    delta = np.r_[
        np.repeat(np.sqrt(1.5 * (dist @ dist)) * eps, len(blocks[1, None]["number"])),
        np.full(2 * len(times), eps),
    ]
    raw = []
    measured = []
    theory = []
    for j, r in enumerate(states):
        port_probs = {
            str(p): source.matched.trace(r, blocks[1, p]["number"])
            for p in (None, *origins)
        }
        counts = {k: rng.binomial(shots, v.clip(0, 1)) for k, v in port_probs.items()}

        def radius(arr):
            return sum((arr[str(p)] - arr["None"]) @ dist for p in origins) / len(
                origins
            )

        clock_prob = [
            (1 + np.cos(times[:, None] * freq) @ weights[:, j]) / 2
            for freq, weights, _ in banks
        ]
        clock_counts = [rng.binomial(shots, v.clip(0, 1)) for v in clock_prob]
        measured.append(np.concatenate([radius(counts), *clock_counts]) / shots)
        theory.append(np.concatenate([radius(port_probs), *clock_prob]))
        raw.append(
            dict(
                name=order[j],
                ports={k: v.tolist() for k, v in counts.items()},
                clocks=[z.tolist() for z in clock_counts],
            )
        )
    observed = np.array(measured)
    theory = np.array(theory)
    assert np.all(abs(observed - theory) < delta)
    ds, dc = observed[2] - observed[1], observed[4] - observed[3]
    ds0, dc0 = theory[2] - theory[1], theory[4] - theory[3]
    held = []
    for j in (5, 6):
        u = float((states[j][1, 1] + states[j][2, 2]).real)
        alpha = 2 * u - 1
        # This is an exact operator identity after pinching by the conserved
        # quantity. It is enough to determine ALL the enrolled instruments.
        difference = (
            states[j]
            - baseline
            - alpha * (standards[1] - standards[0])
            - (1 - alpha) * (standards[3] - standards[2])
        ) * keep
        assert np.max(abs(difference)) < 1e-14
        hits = int(rng.binomial(shots, u))
        a = 2 * np.clip(hits / shots, 0.25, 1) - 1
        prediction = observed[0] + a * ds + (1 - a) * dc
        limit = 6 * delta + 2 * eps * abs(ds - dc)
        np.testing.assert_allclose(
            theory[0] + alpha * ds0 + (1 - alpha) * dc0, theory[j], atol=4e-12
        )
        assert np.all(abs(observed[j] - prediction) < limit)
        held.append(
            dict(
                name=order[j],
                source_agreements=hits,
                ordinary_orbital_weight=a,
                predicted=prediction.tolist(),
                direct=observed[j].tolist(),
                simultaneous_half_width=limit.tolist(),
            )
        )
    return dict(
        trials=shots,
        clock_times=times.tolist(),
        flight_times=[0, 1, 4, 16, 64, "infinite_mean"],
        labels=order,
        held=held,
        counts=raw,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    r = run()
    path = HERE / "coincidence.json"
    if args.write:
        path.write_text(json.dumps(r, indent=2) + "\n")
    else:
        source.region.check(json.loads(path.read_text()), r)
