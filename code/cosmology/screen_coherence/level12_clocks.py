"""Reciprocal clock pulses and material radar on a native strict sphere."""

import argparse
from collections import deque
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path

from strict_history import archive, local_replay
from oph_model_sim import bigboard, boards

EPOCHS = (0, 4, 16, 20, 96, 100)
LIMIT = 12
LEVEL = 12


def arrivals(model, origin, send, targets, transparent=False):
    pending = set(targets)
    times = {origin: send}
    queue = deque([origin])
    while queue and pending:
        u = queue.popleft()
        tick = times[u]
        pending.discard(u)
        if tick == send + LIMIT:
            continue
        for v, closure in model.links(u):
            if v not in times and (transparent or not closure or tick < closure):
                times[v] = tick + 1
                queue.append(v)
    return {target: times.get(target) for target in targets}


def survey(model, observer, transparent=False):
    origin = 12 * observer
    peers = [12 * model.neighbor(model.level, observer, slot)[0] for slot in range(3)]
    results = []
    for send in (0, 4) if transparent else EPOCHS:
        outward = arrivals(model, origin, send, peers, transparent)
        for peer in peers:
            remote_hit = outward[peer]
            received = (
                None
                if remote_hit is None
                else arrivals(model, peer, remote_hit, [origin], transparent)[origin]
            )
            pulse = arrivals(model, peer, send, [origin], transparent)[origin]
            results.append(
                dict(
                    send=send,
                    peer=peer // 12,
                    challenge_received=remote_hit,
                    echo_received=received,
                    incoming_pulse=pulse,
                )
            )
    return results


def summarize(results, starts):
    by_time = {(r["send"], r["peer"]): r for r in results}
    peers = sorted({r["peer"] for r in results})
    counts = {}
    for start in starts:
        shifts = []
        ratios = []
        for p in peers:
            first = by_time[start, p]
            last = by_time[start + 4, p]
            a, b = first["incoming_pulse"], last["incoming_pulse"]
            if a is not None and b is not None:
                shifts.append(str(Fraction(b - a, 4) - 1))
            a, b = first["echo_received"], last["echo_received"]
            if a is not None and b is not None:
                ratios.append(str(Fraction(b - start - 4, a - start)))
        counts[start] = dict(
            pulse_shifts=shifts,
            timed_echo_ratios=ratios,
            unequal_echo_ratio_pairs=sum(a != b for a, b in combinations(ratios, 2)),
            eligible_echo_ratio_pairs=len(list(combinations(ratios, 2))),
        )
    return counts


def measure(level=LEVEL):
    if not 0 <= level <= 13:
        raise ValueError(
            "This survey uses the native 32-bit observer address space (levels 0 through 13)"
        )
    model = local_replay.LocalReplay(level)
    origins = list(
        dict.fromkeys([0] + [(face + 1) * 4**level - 1 for face in range(20)])
    )
    stations = []
    for observer in origins:
        real = survey(model, observer)
        clear = survey(model, observer, True)
        stations.append(
            dict(
                observer=observer,
                native=real,
                transparent=clear,
                native_summary=summarize(real, (0, 16, 96)),
                transparent_summary=summarize(clear, (0,)),
            )
        )
    sources = [Path(m.__file__) for m in (local_replay, bigboard, boards, archive)] + [
        archive.HERE / "repair.cu",
        archive.HERE / "seams.json",
    ]
    return dict(
        level=level,
        seed=archive.SEED,
        board_size=20 * 4**level,
        service="native sharp strict repair",
        enrollment="Observer zero plus the central descendant of each base face. Each station and its three directly connected observers use port zero. Addresses select the instrument placements.",
        clock="Equal paced local clocks, one tick per serviced material transfer. Both sides schedule labelled pulses. Differences of local arrival times measure received periods. Challenges trigger immediate replies. Every direction has a 12-tick flight window.",
        survey="The three peer emitters and the two-way challenges run on separate matched histories. Late time 96 is a survey setting. Transparent control passes all contacts. Ratios compare round trips with common local send times.",
        input_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        stations=stations,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, choices=range(14), default=LEVEL)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(
        json.dumps(measure(args.level), indent=1, allow_nan=False) + "\n"
    )
