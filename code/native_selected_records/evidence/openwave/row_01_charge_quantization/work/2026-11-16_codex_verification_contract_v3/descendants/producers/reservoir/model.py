"""Native source family and a controller that sees only returned records."""

# V3 verification contract: these checks are unconditional, also under python -O.
def _require(condition, message="Scientific check failed"):
    if not condition:
        raise AssertionError(message)

from native_function import native
import random


def graph(cell, bound):
    edges = [(12*c+a, 12*c+b) for c in range(bound) for a, b in cell]
    edges += [(12*c+v, 12*(c+1)+v) for c in range(bound-1) for v in range(12)]
    return edges


def tree(n, edges):
    seen, queue, result = {0}, [0], []
    for a in queue:
        for u, v in edges:
            b = v if u == a else u if v == a else None
            if b is not None and b not in seen:
                seen.add(b)
                queue.append(b)
                result.append((a, b))
    _require(len(seen) == n, 'Scientific check at original line 21')
    return result


class Source:
    """Only event() is the controller's source interface."""
    def __init__(self, initial, law):
        self._loads = tuple(initial)
        self._law = law

    def event(self, edge, coin, probe=False):
        if probe or self._law == 'native':
            self._loads, (_, current) = native(self._loads, edge, coin)
        else:
            a, b = edge
            out = list(self._loads)
            current = coin * (out[a] - out[b])
            out[a] -= current
            out[b] += current
            self._loads = tuple(out)
        return [edge[0], edge[1], coin, current]


def history(initial, law, seed, edges, horizon, cap):
    source, rng = Source(initial, law), random.Random(seed)
    body = [source.event(edges[rng.randrange(len(edges))], rng.randrange(2))
            for _ in range(horizon)]
    probe = []
    for edge in tree(len(initial), edges):
        first = None
        for _ in range(cap):
            coin = rng.randrange(2)
            probe.append(source.event(edge, coin, probe=True))
            if first is None:
                first = coin
            elif coin != first:
                break
        else:
            return body, probe, 'inconclusive'
    return body, probe, 'conclusive'
