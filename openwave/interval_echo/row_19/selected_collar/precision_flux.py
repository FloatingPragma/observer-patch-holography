"""Independent high-count phase meter on the frozen held optical records."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import measure, observer_records

HERE = Path(__file__).resolve().parent
SHOTS = 2**78
BATCH = 2**62


def attach(held_path, extra):
    if sha256(held_path.read_bytes()).hexdigest() != extra["optical_record_sha256"]:
        raise ValueError("Held clock or reception record differs")
    bank = json.loads(held_path.read_text())
    bank["shots_per_flux_setting"] = SHOTS
    for name, settings in extra["flux_plus_counts"].items():
        bank["stocks"][name]["flux_plus_counts"] = settings
    return observer_records.decode(bank)


def acquire(path, time, seed, states, board, route, rings):
    _, masks = observer_records.phase_read.sums_for_charge(
        board, route, rings, -1, time
    )
    rng = np.random.default_rng(seed)
    counts = {}
    for name, x in states.items():
        energy = np.einsum("btij,bji->t", np.array(masks["active"]), x).real
        probs = (1 + energy / 5376) / 2
        counts[name] = [
            sum(map(int, rng.binomial(BATCH, float(p), size=SHOTS // BATCH)))
            for p in probs
        ]
    return dict(
        shots_per_setting=SHOTS,
        seed=seed,
        optical_record_sha256=sha256(path.read_bytes()).hexdigest(),
        flux_plus_counts=counts,
    )


def run(counts_only, write):
    if not counts_only:
        search = json.loads((HERE / "search.json").read_text())
        api = measure.jacobi.beam
        prior = measure.bridge.tangents()
        ref = prior["reference"]
        basis = measure.joint_basis()
        direction = np.einsum(
            "i,ibjk->bjk", search["joint_reading_coefficients"], basis
        )
        states = {
            k + "_" + tag: ref + sign * h * d
            for k, h, d in [
                ("ordinary", 1 / 2048, prior["ordinary"]),
                ("collar", search["step"], direction),
            ]
            for tag, sign in [("minus", -1), ("plus", 1)]
        }
        formation = api.front.optics.pairing.formation
        board = formation.federation(1)
        route = formation.region(board)
        rings = api.screens(board, route)
    for t, seed in [(3.5, 190817), (4.5, 190818)]:
        held = HERE / f"held_{t}.json"
        path = HERE / f"flux_{t}.json"
        bank = (
            json.loads(path.read_text())
            if counts_only
            else acquire(held, t, seed, states, board, route, rings)
        )
        decoded = attach(held, bank)
        if write:
            bank["decoded"] = decoded
            path.write_text(json.dumps(bank, indent=2, allow_nan=False) + "\n")
        else:
            assert decoded == json.loads(path.read_text())["decoded"]
        print(
            t,
            decoded["predicted_odd_focusing"],
            decoded["measured_odd_focusing"],
            flush=True,
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts-only", action="store_true")
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    run(args.counts_only, args.write)
