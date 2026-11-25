"""Cross-check the complete echo survey with separate forward single-source flights."""

from collections import deque
from fractions import Fraction
import json
from pathlib import Path
from random import Random

from observer_zero_reach import echo_waits, return_delays

HERE = Path(__file__).resolve().parent


def forward(links, source, send):
    times = {source: send}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v, closure in links(u):
            if v not in times and (not closure or times[u] < closure):
                times[v] = times[u] + 1
                queue.append(v)
    return times


def test_complete_probe_without_cutoff():
    for seed in range(7):
        rng = Random(seed)
        edges = [[] for _ in range(40)]
        for u in range(len(edges)):
            for v in range(u):
                if rng.random() < 0.09 or v == u - 1:
                    closure = rng.choice((0, 0, 1, 3, 7, 10))
                    edges[u].append((v, closure))
                    edges[v].append((u, closure))
        links = edges.__getitem__
        back = return_delays(links, 10)
        for send in (0, 1, 3, 6, 10, 10**6):
            expected = {}
            for station, arrived in forward(links, 0, send).items():
                if station // 12:
                    received = forward(links, station, arrived).get(0)
                    if received is not None:
                        expected[station] = received - send
            assert echo_waits(links, back, send) == expected


def test_archived_complete_and_aperture_receipts():
    reach = json.loads((HERE / "observer_zero_reach.json").read_bytes())
    old = json.loads((HERE / "observer_zero_radar.json").read_bytes())
    fixed = json.loads((HERE / "observer_zero_capacity.json").read_bytes())
    assert (
        reach["manifest_sha256"] == old["manifest_sha256"] == fixed["manifest_sha256"]
    )
    labels = list(map(tuple, reach["source_labels"]))
    tail = list(map(tuple, reach["fixed_cohort_labels"]))
    assert set(tail) == set(map(tuple, fixed["reflectors"]))
    intervals = reach["send_intervals"]
    assert (len(labels), len(tail)) == (174, 156)
    assert [
        i
        for r in intervals
        for i in range(r["first_send_tick"], r["last_send_tick"] + 1)
    ] == list(range(97))
    for send in fixed["origin_emission_ticks"]:
        row = next(
            r
            for r in intervals
            if r["first_send_tick"] <= min(send, 96) <= r["last_send_tick"]
        )
        waits = dict(zip(labels, row["round_trip_ticks"], strict=True))
        assert [waits[tuple(p)] for p in fixed["reflectors"]] == fixed["echo_waits"][
            send
        ]
        assert max(t for t in waits.values() if t is not None) == 32
        assert Fraction(row["fixed_cohort_mean_range"]) == Fraction(
            sum(waits[p] for p in tail), 2 * len(tail)
        )
        if send in old["emission_ticks"]:
            assert [waits.get((a, b)) for a, b, _ in old["reflectors"]] == old[
                "echo_waits"
            ][old["emission_ticks"].index(send)]
