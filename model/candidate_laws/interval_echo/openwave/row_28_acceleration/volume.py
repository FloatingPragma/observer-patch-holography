"""Count distinct reflected marker IDs on one local clock during native repair."""

import argparse
import json
from tempfile import NamedTemporaryFile
from pathlib import Path
from functools import cache
from itertools import combinations
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from .experiment import HERE, layout, calibrate
from .signals import adjacency

EMISSIONS = (0, 4, 8, 12, 16, 24, 32, 48, 52, 56)
WINDOW = 12


def source_clock(length):
    law = candidate.native_law
    kinetic, contrast = law.local_term_parts()
    pulse = expm(-1j * np.pi * (kinetic + candidate.KAPPA * contrast))
    start = law.SRC.index((0, 1))
    other = law.SRC.index((1, 0))
    rng = np.random.default_rng(280335)
    bit = 0
    log = [bit]
    error = 0
    for _ in range(length):
        pos = start if bit == 0 else other
        dest = other if bit == 0 else start
        error = max(error, abs(1 - abs(pulse[dest, pos]) ** 2))
        if error > 1e-12:
            raise ArithmeticError("clock transfer failed")
        born_one = abs(pulse[other, pos]) ** 2
        bit = int(rng.random() < min(1, born_one))
        log.append(bit)
    return log


def replies(graphs, emission, timeout):
    # Each column is a different reflector on an independent matched copy.
    # Pending mail may be retained. This implementation permits waiting.
    n = graphs[0].shape[0]
    outgoing = np.zeros(n, dtype=np.int8)
    outgoing[0] = 1
    returning = np.zeros((n, n), dtype=np.int8)
    returning[0, 0] = 1
    result = np.full(n, -1, dtype=int)
    result[0] = 0
    for delay in range(1, timeout + 1):
        graph = graphs[min(emission + delay - 1, len(graphs) - 1)]
        next_out = np.maximum(outgoing, graph @ outgoing > 0)
        returning = np.maximum(returning, graph @ returning > 0)
        new = np.flatnonzero(next_out > outgoing)
        returning[new, new] = 1
        first = (result < 0) & (returning[0] > 0)
        result[first] = delay
        outgoing = next_out
    return result


def stationary_replies(graph):
    # Once settled both legs use the same stationary undirected graph.
    n = graph.shape[0]
    distances = np.full(n, -1, dtype=int)
    distances[0] = 0
    active = np.zeros(n, dtype=np.int8)
    active[0] = 1
    tick = 0
    while active.any():
        tick += 1
        active = ((graph @ active > 0) & (distances < 0)).astype(np.int8)
        distances[active > 0] = tick
    return 2 * distances.clip(min=0) - (distances < 0).astype(int)


def summary(reads):
    largest = int(reads.max())
    cutoffs = list(range(0, largest + 3, 2))
    counts = [int(np.sum((reads >= 0) & (reads <= window))) for window in cutoffs]
    return dict(
        cutoffs=cutoffs,
        marker_counts=counts,
        cube_root_scales=[round(v ** (1 / 3), 10) for v in counts],
    )


def tetrahedra(graphs):
    # A message carries its ordered destinations. An observer needs only
    # their own stamps; the propagation internals below generate those stamps.
    @cache
    def local_hits(origin, t):
        seen = np.full(graphs[0].shape[0], -1, dtype=int)
        seen[origin] = t
        front = np.zeros(len(seen), np.int8)
        front[origin] = 1
        while front.any() and np.any(seen[:12] < 0):
            graph = graphs[min(t, len(graphs) - 1)]
            t += 1
            front = ((graph @ front > 0) & (seen < 0)).astype(np.int8)
            seen[front > 0] = t
        return seen[:12]

    targets = range(1, 12)
    records = []
    for t in EMISSIONS:

        def tour(waypoints):
            at = 0
            now = t
            for dest in waypoints:
                now = int(local_hits(at, now)[dest])
                if now < 0:
                    return None
                at = dest
            return now - t

        echoes = {a: tour([a, 0]) for a in targets}
        trip = {(a, b): tour([a, b, 0]) for a in targets for b in targets if a != b}
        volumes = {}
        for triple in combinations(targets, 3):
            ids = (0,) + triple
            doubled = np.zeros((4, 4), dtype=int)
            okay = True
            for i, j in combinations(range(4), 2):
                a, b = ids[i], ids[j]
                if a == 0:
                    value = echoes[b]
                elif (
                    trip[a, b] is None
                    or trip[b, a] != trip[a, b]
                    or echoes[a] is None
                    or echoes[b] is None
                ):
                    value = None
                else:
                    value = 2 * trip[a, b] - echoes[a] - echoes[b]
                if value is None or value <= 0:
                    okay = False
                    break
                doubled[i, j] = doubled[j, i] = value
            if not okay:
                outcome = "unresolved"
            else:
                d = doubled**2
                gram8 = d[0, 1:, None] + d[0, None, 1:] - d[1:, 1:]
                det = int(gram8[0] @ np.cross(gram8[1], gram8[2]))
                minors = [
                    int(gram8[i, i] * gram8[j, j] - gram8[i, j] ** 2)
                    for i, j in combinations(range(3), 2)
                ]
                outcome = det if det >= 0 and min(minors) >= 0 else "non_euclidean"
            volumes[",".join(map(str, triple))] = outcome
        records.append(
            dict(
                emission=t,
                echo=[echoes[a] for a in targets],
                tours=[[0 if a == b else trip[a, b] for b in targets] for a in targets],
                determinants=volumes,
            )
        )
    return records


def run():
    n, edges, *_ = layout(0)
    passing = adjacency(n, edges)
    control = replies([passing], 0, WINDOW)
    local_transfers = calibrate()
    ref_tetra = tetrahedra([passing])[0]["determinants"]
    label = next(k for k, v in ref_tetra.items() if isinstance(v, int) and v > 0)
    panels = []
    with np.load(HERE / "records.npz") as saved:
        for seed in range(280210, 280218):
            commit = saved[f"epoch_{seed}"]
            last = int(max(commit))
            graphs = [
                adjacency(n, edges[(commit == 0) | (commit > t)])
                for t in range(last + 1)
            ]
            survey = [replies(graphs, t, WINDOW) for t in EMISSIONS]
            settled = stationary_replies(graphs[-1])
            counts = [int((row >= 0).sum()) for row in survey]
            if any(b > a for a, b in zip(counts, counts[1:])):
                raise AssertionError("fixed-window marker count grew")
            panels.append(
                dict(
                    seed=seed,
                    tetrahedra=tetrahedra(graphs),
                    last_write=last,
                    writes=int(np.sum(commit > 0)),
                    counts=counts,
                    return_delays=[r.tolist() for r in survey],
                    settled_return_delays=settled.tolist(),
                    growing_aperture=summary(settled),
                )
            )
    return dict(
        selected_blank_tetra=label,
        local_transfers=local_transfers,
        clock_bits=source_clock(
            max(
                max(EMISSIONS) + WINDOW,
                max(
                    p["last_write"] + max(p["settled_return_delays"]) + 2
                    for p in panels
                ),
            )
        ),
        emissions=EMISSIONS,
        window=WINDOW,
        marker_count=n,
        transmitting_control=control.tolist(),
        panels=panels,
    )


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--write", action="store_true")
    data = run()
    path = HERE / "volume.json"
    as_text = json.dumps(data, separators=(",", ":")) + "\n"
    if args.parse_args().write:
        with NamedTemporaryFile("w", dir=path.parent, delete=False) as new:
            new.write(as_text)
        Path(new.name).replace(path)
    elif json.loads(path.read_text()) != json.loads(as_text):
        raise ValueError("Volume or clock record differs")
    else:
        print(
            "Verified all 19,200 time-limited returns, local clock, and eight settled surveys."
        )


if __name__ == "__main__":
    main()
