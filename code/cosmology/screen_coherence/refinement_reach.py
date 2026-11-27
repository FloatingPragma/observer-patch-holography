"""Compare clock receipts across levels and test settled one-tick gateway replies."""

import json
from collections import Counter
from pathlib import Path

from strict_history.local_replay import LocalReplay

HERE = Path(__file__).resolve().parent


def gateway_trials(model, owner):
    trials = []
    for p in model.ports(owner):
        u = 12 * owner + p
        ((v, closure),) = [(v, t) for v, t in model.links(u) if v // 12 != owner]
        peer, q = divmod(v, 12)
        if (u, closure) not in model.links(v):
            raise ValueError("Remote endpoint disagrees on material transmission")
        own_initial = (int.from_bytes(model.row(owner)[:2], "little") >> p) & 1
        peer_initial = (int.from_bytes(model.row(peer)[:2], "little") >> q) & 1
        trials.append(
            dict(
                port=p,
                peer=peer,
                peer_port=q,
                initial_loads=[2 * own_initial, 2 * peer_initial],
                settled_direct_echo_wait=None if closure else 2,
            )
        )
        if bool(closure) != (own_initial != peer_initial):
            raise ValueError("Gateway control disagrees with initial endpoint readings")
    return trials


def depth_bounds(depths):
    # q is survival to depth r-1 on a binary tree. The root has three offers.
    q = 1.0
    bounds = {}
    for r in range(1, max(depths) + 1):
        if r in depths:
            bounds[r] = 1 - (1 - q / 2) ** 3
        q *= 1 - q / 4
    return bounds


def compare(level):
    clocks = json.loads((HERE / f"level{level}_clocks.json").read_bytes())
    model = LocalReplay(level)
    stations = clocks["stations"]
    trials = {s["observer"]: gateway_trials(model, s["observer"]) for s in stations}
    hist = Counter(
        sum(r["settled_direct_echo_wait"] is not None for r in rows)
        for rows in trials.values()
    )
    return dict(
        direct_receipts=trials,
        responding_direct_neighbors=hist,
        paired_pulse_shifts={
            start: Counter(
                z for s in stations for z in s["native_summary"][start]["pulse_shifts"]
            )
            for start in ("0", "16", "96")
        },
    )


if __name__ == "__main__":
    report = dict(
        levels={level: compare(level) for level in (12, 13)},
        clock_seed=20261126,
        ensemble_conditions="Independent fair initial port loads, eventual service of each eligible external contact, three distinct gateway ports per observer. Fixed-seed samples are reported separately.",
        settled_direct_neighbor_count_probabilities=["1/8", "3/8", "3/8", "1/8"],
        settled_material_minimum_handoff_depth_probability_upper_bounds=depth_bounds(
            (1, 10, 100, 1000)
        ),
        level202_observer_count=str(20 * 4**202),
    )
    print(json.dumps(report, indent=1, allow_nan=False))
