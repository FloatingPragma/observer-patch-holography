"""Every remote material echo accessible to observer zero during repair and settlement."""

import argparse
import json
from collections import deque
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from observer_zero_capacity import component
from observer_zero_radar import ArchiveRadar, MANIFEST_SHA


def return_delays(links, settled):
    """At each departure tick, find every port with a finite flight time home."""
    delays = [None] * (settled + 1)
    delays[-1] = component(links, 0)
    for t in reversed(range(settled)):
        current = {0: 0}
        for v, remaining in delays[t + 1].items():
            for u, closure in links(v):
                if closure and t >= closure:
                    continue
                if u not in current or remaining + 1 < current[u]:
                    current[u] = remaining + 1
        delays[t] = current
    return delays


def echo_waits(links, delays, send):
    """Flood forward only through sites that can still send a reply home."""
    queue = deque([(0, send)])
    seen = {0}
    echoes = {}
    while queue:
        u, t = queue.popleft()
        remaining = delays[min(t, len(delays) - 1)].get(u)
        if remaining is None:
            continue
        if u // 12 != 0:
            echoes[u] = t - send + remaining
        for v, closure in links(u):
            if v not in seen and (not closure or t < closure):
                seen.add(v)
                queue.append((v, t + 1))
    return echoes


def measure(archive):
    reader = ArchiveRadar(archive)
    settled = reader.manifest["settlement_round"]

    @lru_cache(maxsize=None)
    def links(port):
        return tuple(reader.links(port))

    delays = return_delays(links, settled)
    if any((u, t) not in links(v) for u in delays[0] for v, t in links(u)):
        raise ValueError("A material contact has mismatched endpoint histories")
    surveys = [echo_waits(links, delays, s) for s in range(settled + 1)]
    labels = sorted(surveys[0])
    fixed = sorted(surveys[-1])
    if any(row.keys() - previous.keys() for previous, row in zip(surveys, surveys[1:])):
        raise ValueError("Lost material routes cannot reopen")
    if surveys[-1] != {p: 2 * d for p, d in delays[-1].items() if p // 12}:
        raise ValueError("Settled echoes disagree with the static return times")
    intervals = []
    for s, row in enumerate(surveys):
        if s and row == surveys[s - 1]:
            intervals[-1]["last_send_tick"] = s
            continue
        intervals.append(
            dict(
                first_send_tick=s,
                last_send_tick=s,
                remote_ports=len(row),
                remote_observers=len({p // 12 for p in row}),
                maximum_radar_range=str(Fraction(max(row.values()), 2)),
                fixed_cohort_mean_range=str(
                    Fraction(sum(row[p] for p in fixed), 2 * len(fixed))
                ),
                round_trip_ticks=[row.get(p) for p in labels],
            )
        )
    return dict(
        archive_manifest="https://dovybzrnzxfo1.cloudfront.net/manifest.json",
        manifest_sha256=MANIFEST_SHA,
        receiver=[0, 0],
        settlement_round=settled,
        apparatus="Matched material relay probes. A transfer takes one tick. Each named remote port reflects immediately. Each entry is the earliest returning echo from that port. No radius or flight-time cutoff.",
        source_labels=[list(divmod(p, 12)) for p in labels],
        fixed_cohort_labels=[list(divmod(p, 12)) for p in fixed],
        send_intervals=intervals,
        stationary_continuation="For any later send, each reachable beacon returns the same round-trip wait as at settlement for as long as the clock and relay service are operated.",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(measure(args.archive), indent=1) + "\n")
