"""Embed an extended neutral receiver using only real contacts of the federation."""

from collections import deque
from model.executable_model import federation
from ..passive_decay.material import commission, edge


def embed(board, receivers):
    near = tuple(set(words) for words in board.neighbours)
    dimers = tuple(
        tuple(
            (first, other)
            for other in sorted(near[first])
            if board.degree[first] == board.degree[other]
        )
        for first in range(board.n)
    )
    tried = 0

    def attach(port, occupied, triangle):
        for p in sorted(near[port] - occupied):
            for y in dimers[p]:
                if y[1] not in occupied and (not triangle or y[1] in near[port]):
                    yield y

    def wire(core_a, core_b, ordered):
        a, b = core_a, core_b
        stock = (a, b, *ordered)
        used = set(sum(stock, ()))
        heads = {edge(a[0], b[0]): 32}
        pin_requests = [(a[1], 32), (b[1], 33)]
        for x, y in zip((a, *ordered[:-1]), ordered, strict=True):
            heads[edge(x[0], y[0])] = 33
            heads[edge(x[0], y[1])] = 32
            pin_requests.append((x[1], 32))
        for rail, cost in pin_requests:
            choices = sorted(p for p in near[rail] - used if edge(rail, p) not in heads)
            if not choices:
                return None
            heads[edge(rail, choices[0])] = cost
        return stock, heads

    def trail(a, b, ordered, occupied):
        nonlocal tried
        if len(ordered) == receivers:
            tried += 1
            return wire(a, b, ordered)
        port = ordered[-1][0] if ordered else a[0]
        for next_dimer in attach(port, occupied, True):
            hit = trail(a, b, (*ordered, next_dimer), occupied | set(next_dimer))
            if hit:
                return hit
        return None

    for choices in dimers:
        for a in choices:
            for b in attach(a[0], set(a), False):
                hit = trail(a, b, (), set(a + b))
                if hit:
                    return hit, tried
    return None, tried


def route(lengths=(4, 3), crossing=False):
    search = []
    for receivers in lengths:
        for level in range(3):
            board = federation(level)
            hit, tried = cross_owner(board) if crossing else embed(board, receivers)
            search.append(
                dict(
                    receivers=receivers,
                    cells=board.cells,
                    full_routes_visited=tried,
                    found=bool(hit),
                )
            )
            if hit:
                dimers, heads = hit
                prepared = commission(board, dimers, heads)
                sources = set(sum(dimers[:2], ()))
                ends = set(dimers[-1])
                dist = {p: 0 for p in sources}
                queue = deque(sources)
                while queue:
                    p = queue.popleft()
                    for q in board.neighbours[p]:
                        if q not in dist:
                            dist[q] = dist[p] + 1
                            queue.append(q)
                geometry = dict(
                    owner_ids=[[a // 12, b // 12] for a, b in dimers],
                    raw_seams_to_each_end=[dist[q] for q in sorted(ends)],
                    dimer_links_to_end=len(dimers) - 2,
                    stock=prepared,
                )
                return board, dimers, heads, search, geometry
    raise ValueError(search)


def cross_owner(board):
    near = tuple(set(ns) for ns in board.neighbours)
    pairs = [
        (p, q)
        for p in range(board.n)
        for q in sorted(near[p])
        if board.degree[p] == board.degree[q]
    ]
    by_first = {}
    for p, q in pairs:
        by_first.setdefault(p, []).append((p, q))

    def next_dimer(p, unused, tri):
        for rail in sorted(near[p] - unused):
            for x, y in by_first.get(rail, []):
                if y not in unused and (not tri or y in near[p]):
                    yield x, y

    tries = 0
    for c in pairs:
        if c[0] // 12 == c[1] // 12:
            continue
        for a in next_dimer(c[0], set(c), True):
            for b in next_dimer(a[0], set(c + a), False):
                core = set(a + b)
                for d in next_dimer(c[1], set(c + a + b), True):
                    for e in next_dimer(d[0], set(c + a + b + d), True):
                        tries += 1
                        if {p // 12 for p in e} & {p // 12 for p in core}:
                            continue
                        if any(near[q] & core for q in e):
                            continue
                        dimers = (a, b, c, d, e)
                        used = set(sum(dimers, ()))
                        heads = {edge(a[0], b[0]): 32}
                        pins = [(a[1], 32), (b[1], 33)]
                        for x, y in ((c, a), (c[::-1], d), (d, e)):
                            heads[edge(x[0], y[0])] = 33
                            heads[edge(x[0], y[1])] = 32
                            pins.append((x[1], 32))
                        for rail, cost in pins:
                            choices = sorted(
                                p
                                for p in near[rail] - used
                                if edge(rail, p) not in heads
                            )
                            if not choices:
                                break
                            heads[edge(rail, choices[0])] = cost
                        else:
                            return (dimers, heads), tries
    return None, tries
