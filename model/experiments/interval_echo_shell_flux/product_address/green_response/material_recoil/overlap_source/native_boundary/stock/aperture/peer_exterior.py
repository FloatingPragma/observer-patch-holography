"""Read one ordinary bearing at finite deadlines and at first contact with the source."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import landmark_arrival as arrival

HERE = arrival.HERE


def run():
    board, marks, tables, mass = arrival.pulse_inputs()
    adj, slots = arrival.passage.contacts(board, marks)
    pins = [25, 36, *marks[-1]]
    F, den = arrival.blank_exits(adj, pins)
    fields = arrival.eventual(tables[2], pins, F, den)
    # An independently prepared ordinary double defines the position-only
    # bearing at any station by its probability of reaching 36 before 25.
    # Both preparations and both released loads use this single fixed chart.
    coord = mass["original"] @ fields
    n = board.n
    kernels = tables[2].reshape(11, n, 24)[:2] % n
    for t, table in tables.items():
        for m, row in enumerate(table.reshape(11, n, 24)):
            if t == 1 and m in arrival.source.PARTIAL:
                continue
            assert np.array_equal(row % n, kernels[int(m != 0)])
            assert np.all((row // n == 0) == (m == 0))
    live = np.ones(n, bool)
    live[pins] = False
    kernels[:, pins[:2]] = np.array(pins[:2])[None, :, None]
    assert np.array_equal(kernels[0, live], kernels[1, live])
    drifts = []
    for table in kernels:
        diff = coord[table].sum(1) - 24 * coord
        assert not diff[live].any() and not diff[pins[:2]].any()
        drifts.append([str(diff[i] / 24) for i in pins[2:]])
    stopped = kernels.copy()
    stopped[:, pins[2:]] = np.array(pins[2:])[None, :, None]
    assert np.array_equal(stopped[0], stopped[1])
    assert np.all(coord[stopped[0]].sum(1) == 24 * coord)
    measured = json.loads((HERE / "peer_register_lab.json").read_bytes())
    da = Q(measured["readings"]["blank"]["exact_change"])
    x = 104
    distances = shortest_path(adj, unweighted=True, indices=pins[2:]).min(0)
    start_distance = int(distances[x])
    grid = 2**47
    lower = np.tile([math.floor(v * grid) for v in coord], (2, 1))
    upper = np.tile([math.ceil(v * grid) for v in coord], (2, 1))
    steps = {
        0,
        1,
        2,
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
        4096,
        start_distance,
        start_distance + 1,
    }
    panel = {}
    for k in range(max(steps) + 1):
        if k in steps:
            margin = [
                Q(int(lower[0, x] - upper[1, x]), grid),
                Q(int(upper[0, x] - lower[1, x]), grid),
            ]
            panel[str(k)] = list(map(float, [da * z for z in margin]))
        if k <= start_distance:
            assert lower[0, x] == lower[1, x] and upper[0, x] == upper[1, x]
        lower = np.stack(
            [v[t].sum(1) // 24 for v, t in zip(lower, kernels, strict=True)]
        )
        upper = np.stack(
            [(v[t].sum(1) + 23) // 24 for v, t in zip(upper, kernels, strict=True)]
        )
    # Independent launch tickets, common random numbers across the two
    # strata. Up to the first contact their paths coincide ticket for ticket.
    rng = np.random.default_rng(519073)
    count = 2**17
    sites = np.full((2, count), x)
    escaped = np.zeros(count, bool)
    fixed = np.array(coord, float)
    held = {}
    deadlines = (16, 64, 256)
    for step in range(1, max(deadlines) + 1):
        escaped |= np.isin(sites[0], pins[2:])
        tickets = rng.integers(24, size=count)
        for j in (0, 1):
            sites[j] = kernels[j, sites[j], tickets]
        assert np.array_equal(sites[0, ~escaped], sites[1, ~escaped])
        if step in deadlines:
            d = fixed[sites[0]] - fixed[sites[1]]
            mean = float(d.mean()) * float(da)
            width = float(da) * math.sqrt(
                2 * math.log(2 * len(deadlines) / 0.01) / count
            )
            truth = panel[str(step)]
            assert mean - width <= truth[0] <= truth[1] <= mean + width
            held[str(step)] = dict(
                seen_source_before_step=int(escaped.sum()),
                mean=mean,
                bounds=[mean - width, mean + width],
            )
    eventual = da * (fields[0, x] - fields[1, x])
    assert eventual == Q(measured["exact_arrival_change"])
    return dict(
        sha256={
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("peer_register_lab.json", "LANDMARK_ARRIVAL.md")
        },
        launch=x,
        beacons=pins[:2],
        source_contacts=pins[2:],
        distance_to_source=start_distance,
        ordinary_bearing_at_launch=str(coord[x]),
        open_closed_local_bearing_steps=drifts,
        first_contact_stopped_change_all_stations="0",
        fixed_deadline_changes=panel,
        eventual_change=str(eventual),
        held_trials_per_pair=count,
        joint_coverage=0.99,
        observed=held,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = HERE / "peer_exterior.json"
    if args.write:
        path.write_text(raw)
    else:
        assert raw == path.read_text()
    print(
        out["distance_to_source"],
        out["fixed_deadline_changes"]["64"],
        float(Q(out["eventual_change"])),
    )
