"""Literal integer source. Only the event closure is given to the controller."""
import random
CALLS = {'native': 0, 'swap_hold': 0}
CELL = ((0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3), (2, 4), (2, 6), (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9), (4, 11), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11))
TRI = ((0, 1), (1, 2), (0, 2))

def native(q, edge, coin, law='native'):
    CALLS[law] += 1
    a, b = edge
    d = q[a] - q[b]
    j = (d - (2 * coin - 1) * (abs(d) % 2)) // 2 if law == 'native' else coin * d
    out = list(q)
    out[a] -= j
    out[b] += j
    return (out, [a, b, coin, j])

def graph(bound, m):
    n = 12 * bound
    core = [(12 * c + a, 12 * c + b) for c in range(bound) for a, b in CELL]
    core += [(12 * c + v, 12 * (c + 1) + v) for c in range(bound - 1) for v in range(12)]
    target = [0, n, n + 1]
    copy = [n + 2, n + 3, n + 4]
    blank = [n + 5, n + 6]
    rings = [list(range(n + 7, n + 7 + 3 * m)), list(range(n + 7 + 3 * m, n + 7 + 6 * m))]
    edges = list(core)
    for vs in (target, copy):
        edges += [(vs[a], vs[b]) for a, b in TRI]
    edges += list(zip(target, copy)) + [tuple(blank)]
    for ring in rings:
        edges += list(zip(ring, ring[1:] + ring[:1]))
    bridges = [(0, blank[0]), (0, rings[0][0]), (0, rings[1][0])]
    edges += bridges
    return {'N': n + 7 + 6 * m, 'core_N': n, 'core': core, 'edges': edges, 'target': target, 'copy': copy, 'blank': blank, 'rings': rings, 'idle_bridges': bridges}

def prepare(first_cell, g):
    q = list(first_cell) + [0] * (g['N'] - 12)
    for ring in g['rings']:
        q[ring[2 * (len(ring) // 3)]] = 1
    return q

def source(initial, seed, control=False, forced_coin=None):
    state = list(initial)
    rng = random.Random(seed)
    tape = []

    def request(edge, stage):
        nonlocal state
        coin = rng.randrange(2) if forced_coin is None else forced_coin
        law = 'swap_hold' if control and stage == 'balance' else 'native'
        state, record = native(state, edge, coin, law)
        tape.append([stage, record])
        return list(record)

    def audit():
        return {'final': list(state), 'events': list(tape)}
    return (request, audit)

def require(condition, message='scientific check failed'):
    if not condition:
        raise AssertionError(message)
