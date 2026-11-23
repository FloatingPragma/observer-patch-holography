"""Coherent, addressed transport of the compact four-port release across a real bridge."""
from collections import Counter
from functools import lru_cache
import json

import numpy as np

from model.executable_model import cell, federation
from .local_release import SUPPORT, all_cells, native_instruments, prepare, record


def walk(neighbors, start, end, avoid, length):
    """An exact-length walk. A predecessor at each layer permits repeated vertices."""
    parents = [{start: start}]
    for _ in range(length):
        step = {}
        for v in parents[-1]:
            for w in neighbors[v]:
                if w not in avoid:
                    step.setdefault(w, v)
        parents.append(step)
    if end not in parents[-1]:
        return None
    path = [end]
    for layer in reversed(parents[1:]):
        path.append(layer[path[-1]])
    return path[::-1]


@lru_cache(maxsize=None)
def plan(left, right):
    tile = cell()
    edges = tile.seams + [(12+a, 12+b) for a, b in tile.seams] + [(left, 12+right)]
    neighbors = [[] for _ in range(24)]
    for a, b in edges:
        neighbors[a].append(b)
        neighbors[b].append(a)
    for length in range(1, 49):
        routes = []
        for j, v in enumerate(SUPPORT):
            parked = set(SUPPORT[j+1:]) | {12+w for w in SUPPORT[:j]}
            route = walk(neighbors, v, 12+v, parked, length)
            if route is None:
                break
            routes.append(route)
        if len(routes) == len(SUPPORT):
            return tuple((a, b) for route in routes for a, b in zip(route, route[1:]))
    raise ArithmeticError("No equal-length route within the registered budget")


def columns(pulses):
    # Each sharp one-load pi contact is SWAP times exp(-i*pi*k) on its two ports.
    # A spectator sees 00 and is fixed. Keep the phase exponent for each input column.
    occupants = list(range(24))
    exponents = [0]*24
    for a, b in pulses:
        exponents[occupants[a]] += 1
        exponents[occupants[b]] += 1
        occupants[a], occupants[b] = occupants[b], occupants[a]
    return [occupants.index(v) for v in SUPPORT], [exponents[v] for v in SUPPORT]


def transport():
    native_instruments()
    signed, _ = prepare()
    all_cells(signed)
    charts = {}
    census = {}
    for level in range(3):
        board = federation(level)
        edges = set(board.seams)
        counts = Counter()
        for a, b in board.external:
            for src, dest in ((a, b), (b, a)):
                left, right = src % 12, dest % 12
                pulses = plan(left, right)
                ports = [12*(src//12)+i for i in range(12)] + [12*(dest//12)+i for i in range(12)]
                if any(tuple(sorted((ports[v], ports[w]))) not in edges for v, w in pulses):
                    raise ArithmeticError("The route leaves its two physical cells")
                counts[len(pulses)] += 1
                if (left, right) in charts:
                    continue
                landing, powers = columns(pulses)
                if landing != [12+v for v in SUPPORT] or len(set(powers)) != 1:
                    raise ArithmeticError("The transfer fails for a material column or a relative phase")
                unlinked, _ = columns([(v, w) for v, w in pulses if v//12 == w//12])
                if any(v >= 12 for v in unlinked):
                    raise ArithmeticError("The receiver fired without its bridge")
                # A round trip from the first destination to an empty neighbor leaves
                # the address intact and marks just one amplitude with two extra phases.
                extra = ((12+4, 12+9), (12+9, 12+4))
                if (4, 9) not in cell().seams:
                    raise ArithmeticError("Missing control seam")
                moved, changed = columns((*pulses, *extra))
                if moved != landing or changed != [powers[0]+2, *powers[1:]]:
                    raise ArithmeticError("The phase error does not isolate the first rail")
                reads = {}
                for strength in (0, .5, 8):
                    state = np.zeros(24, complex)
                    state[moved] = signed[list(SUPPORT)]/2 * np.exp(-1j*np.pi*strength*np.array(changed))
                    reads[str(strength)] = record(state, 1, reverse=True).tolist()
                    expected = 1/16 if strength == .5 else 1/4
                    np.testing.assert_allclose(reads[str(strength)][5], expected, rtol=0, atol=1e-12)
                charts[left, right] = {"pulses": pulses, "hops_per_rail": powers[0],
                                       "round_trip_inverse_ports_and_quiet": reads}
        census[board.cells] = dict(sorted(counts.items()))
    return {"pulse_count_by_directed_interface": census,
            "charts": {f"{a},{b}": v for (a, b), v in sorted(charts.items())}}


if __name__ == "__main__":
    print(json.dumps(transport(), indent=2, allow_nan=False))
