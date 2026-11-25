"""Compare the compact radar replay with the original native two-leg path reader."""

import numpy as np
from observer_zero_radar import echo_table, local_graph
from model.candidate_laws.interval_echo.openwave.row_28_acceleration import (
    experiment,
    radar,
)


def test_native_round_trips_and_misses():
    n, edges, ranks, *_ = experiment.layout(0)
    rng = np.random.default_rng(28)
    initial = 2 * rng.integers(0, 2, n, dtype=np.int8)
    closures, _, settled = experiment.history(n, edges, ranks, initial, 0.25, rng)
    adjacency = [[] for _ in range(n)]
    for (a, b), closure in zip(edges, closures):
        adjacency[a].append((b, closure))
        adjacency[b].append((a, closure))
    cap = 18
    ports, depth, ball = local_graph(lambda u: adjacency[u], 0, cap)
    targets = [i for i in range(1, len(ports), 5) if ports[i] // 12 != 0]
    emissions = [0, 2, 8, int(settled), int(settled) + 4]
    rows = echo_table(ball, targets, emissions, cap, settled)
    control = echo_table(
        [{q: 0 for q in row} for row in ball], targets, emissions, cap, 0
    )
    assert all(row == [2 * depth[i] for i in targets] for row in control)
    for t, row in zip(emissions, rows):
        receipts = radar.echoes(0, closures, 0, [ports[i] for i in targets], t, False)
        assert row == [
            x - t if x is not None and x - t <= cap else None for x in receipts
        ]
    assert any(v is None for row in rows for v in row)
