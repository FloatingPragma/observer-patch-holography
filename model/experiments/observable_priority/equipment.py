"""Install a candidate in the sharp native service and measure its addressed hits."""

from functools import lru_cache
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import consensus
from model.executable_model import federation
from model.experiments.interval_echo_consensus.__main__ import local_ops
from model.experiments.interval_echo_consensus.asynchronous import SharpService, verify_native_service

from .observed import Screen, need
from .source import bind_source, read


@lru_cache(maxsize=1)
def operators():
    _, _, hit, probes = local_ops()
    return hit, [op for _, op in probes]


class Loads(dict):
    def __missing__(self, _):
        return 1


class Collars(dict):
    def __missing__(self, _):
        return 0


def install(level, candidate):
    bind_source()
    car = federation(level)
    screen = Screen(car.n, car.seams)
    service = SharpService(tuple(range(car.n)), car.seams, reverse=candidate == "reverse")
    if candidate not in ("native", "reverse"):
        # Relabeling breaks ties on disjoint seams and preserves incident comparisons.
        ranks = ([r * (len(service.edges) + 1) - j for j, r in enumerate(service.rank)]
                 if candidate == "relabeled" else read(Path(candidate)))
        need(isinstance(ranks, list) and len(ranks) == len(screen.edges) and
             all(type(i) is int for i in ranks), "invalid candidate ranks")
        higher = [[] for _ in ranks]
        for i, j in screen.wedges:
            need(ranks[i] != ranks[j], "incident priority tie")
            a, b = (i, j) if ranks[i] < ranks[j] else (j, i)
            higher[b].append(a)
        service.rank = ranks
        service.higher = higher
    return screen, service


@lru_cache(maxsize=1)
def addressed_bits():
    # Both signed writes and every blank basis column have deterministic probe bits.
    verify_native_service()
    write, probes = operators()
    occupied = probes[0]
    bits = {}
    for a, b in ((0, 2), (2, 0)):
        state = occupied @ write[:, 3 * a + b]
        p = float(np.vdot(state, state).real)
        need(abs(p - 1) < 1e-11, "native occupied read")
        bits[a - b] = round(p)
    blank = np.linalg.norm(occupied[:, :9], axis=0) ** 2
    need(np.allclose(blank, 0, rtol=0, atol=1e-11), "native blank read")
    bits[0] = 0
    return bits


def apparatus(service):
    # The observer receives only the two occupied/blank outcomes.
    bits = addressed_bits()

    def observe(i, j, preparation):
        center, a, b = preparation
        loads, marks = Loads({center: 2, a: 0, b: 0}), Collars()
        for e in (i, j):
            if service.enabled(e, loads, marks):
                service.commit(e, loads, marks)
        return [bits[marks[i]], bits[marks[j]]]

    return observe


def collect(screen, observe):
    return [[i, j, *observe(i, j, prep)] for (i, j), prep in screen.wedges.items()]


def tensor_wedge(screen, service, instruction):
    """Fresh three-port material wedge, with a separate addressed probe on each collar."""
    i, j = instruction["pair"]
    center = instruction["double"]
    ports = sorted({center, *instruction["vacancies"]})
    loads = Loads({center: 2, **{v: 0 for v in instruction["vacancies"]}})
    marks = Collars()
    ket = np.zeros((11, 11, 3, 3, 3))
    ket[(0, 0, *(loads[v] for v in ports))] = 1
    write, probes = operators()

    def measured(address, state, matrices):
        a, b = screen.edges[address]
        axes = ((i, j).index(address), 2 + ports.index(a), 2 + ports.index(b))
        outcomes = [consensus._apply(state, op, axes) for op in matrices]
        weights = [float(np.vdot(s, s).real) for s in outcomes]
        populated = [k for k, p in enumerate(weights) if p > 1e-11]
        need(len(populated) == 1 and abs(sum(weights) - 1) < 1e-11, "non-sharp wedge")
        k = populated[0]
        return k, outcomes[k] / np.sqrt(weights[k])

    for address in (i, j):
        if service.enabled(address, loads, marks):
            branch, ket = measured(address, ket, [*probes[:2], write])
            need(branch == 2, "native wedge write")
            service.commit(address, loads, marks)
    result = []
    for address in (i, j):
        outcome, after = measured(address, ket, probes)
        need(np.linalg.norm(after - ket) < 1e-11, "disturbed addressed wedge")
        result.append(int(outcome == 0))
    return result
