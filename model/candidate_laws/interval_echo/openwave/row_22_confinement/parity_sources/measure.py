"""Two held plaquette parities, with and without conservation of the background."""
import argparse
from collections import defaultdict, deque
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from ..measure import skeleton

HERE = Path(__file__).resolve().parent
SEED = 221125
SHOTS = 2**17
TIMES = np.exp2(np.arange(-14, 1, 2))
PHASE_ERROR = float(np.arcsin(np.sqrt(4 * np.log(2 * 8192 / 1e-6) / SHOTS)))
RELEASE_SHOTS = 2**25


def pair(a, b):
    return (a, b) if a < b else (b, a)


def edges(path):
    return [pair(a, b) for a, b in zip(path, path[1:])]


def survey(board):
    graph, ports = skeleton(board)
    planar, embedding = nx.check_planarity(graph)
    if not planar or nx.node_connectivity(graph) < 3:
        raise ValueError('Need a spherical polyhedron with unique faces')
    visited, rings = set(), []
    for a in sorted(graph):
        for b in sorted(graph[a]):
            if (a, b) not in visited:
                rings.append(embedding.traverse_face(a, b, visited))
    rings.sort(key=lambda ring: sorted(ring))
    sides = defaultdict(list)
    for f, ring in enumerate(rings):
        for e in edges(ring + ring[:1]):
            sides[e].append(f)
    dual, gates = nx.Graph(), {}
    for (a, b), adjacent in sorted(sides.items()):
        f, g = adjacent
        dual.add_edge(f, g)
        gates[pair(f, g)] = (ports[a, b], ports[b, a])
    # A supplied readiness load, with three distinct internal swap partners per owner.
    parked, partner = set(), {}
    for a in sorted(graph):
        gateways = sorted(ports[a, b] for b in graph[a])
        private = sorted(set(range(12*a, 12*a+12)) - set(gateways))
        ready = next(option for option in permutations(private, 3)
                     if all(q in board.neighbours[p] for p, q in zip(gateways, option, strict=True)))
        partner.update(zip(gateways, ready, strict=True))
        parked.update(ready + tuple(p for p in private if p not in ready)[:3])
    return rings, dual, gates, parked, partner


def prepare(route_edges, gates, parked, partner):
    doubles = parked.copy()
    for e in route_edges:
        p = gates[e][0]
        doubles.remove(partner[p])
        doubles.add(p)
    return doubles


def read_ports(doubles, gates):
    mismatch = [e for e, (a, b) in gates.items() if (a in doubles) != (b in doubles)]
    odd = set()
    for a, b in mismatch:
        odd.symmetric_difference_update((a, b))
    return sorted(odd), sorted(mismatch)


def source_block():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    if np.max(abs(h[9:, :9])) > 1e-11:
        raise ValueError('Blank collar failed to close')
    return exchange[:9, :9], contrast[:9, :9]


def fringe(energy, rng):
    phase = TIMES * energy
    x = rng.binomial(SHOTS, np.clip((1 + np.cos(phase)) / 2, 0, 1))
    y = rng.binomial(SHOTS, np.clip((1 + np.sin(phase)) / 2, 0, 1))
    return [x.tolist(), y.tolist()]


def decode(counts):
    x, y = 2*np.asarray(counts) / SHOTS - 1
    estimate = 0.
    for t, angle in zip(TIMES, np.arctan2(y, x), strict=True):
        branch = round((t * estimate - angle) / (2*np.pi))
        estimate = (angle + 2*np.pi * branch) / t
    return float(estimate)


def read_work(doubles, gates, block, rng):
    # Sum actual sharp-word diagonal entries. The inventory guards remove
    # inter-owner hopping. The 00 bridge reference has zero phase.
    energy = sum(block[3*(2*(a in doubles)) + 2*(b in doubles),
                       3*(2*(a in doubles)) + 2*(b in doubles)] for a, b in gates.values())
    counts = fringe(energy, rng)
    return dict(x_y=counts, measured_work=round(decode(counts), 8))


def detour(path, dual):
    for i, (f, g) in enumerate(zip(path, path[1:])):
        other = sorted((set(dual[f]) & set(dual[g])) - set(path))
        if other:
            return path[:i+1] + [other[0]] + path[i+1:]
    raise ValueError('No triangular detour')


def source_only(f, g, dual):
    if g in dual[f]:
        return [pair(f, g)]
    return [pair(f, min(dual[f])), pair(g, min(dual[g]))]


def syndrome_oracle(dual):
    """Unweighted BFS over all reachable parity bit vectors, independent of paths."""
    moves = [(1 << i) | (1 << j) for i, j in dual.edges]
    distance, front = {0: 0}, deque([0])
    while front:
        bits = front.popleft()
        for move in moves:
            nxt = bits ^ move
            if nxt not in distance:
                distance[nxt] = distance[bits] + 1
                front.append(nxt)
    checked = 0
    for f in sorted(dual):
        for g in sorted(x for x in dual if x > f):
            required = (1 << f) | (1 << g)
            exact = distance[required]
            relaxed = min(n for m, n in distance.items() if m & required == required)
            d = nx.shortest_path_length(dual, f, g)
            if exact != d or relaxed != min(d, 2):
                raise ValueError('Syndrome enumeration disagrees')
            checked += 1
    return dict(reachable_syndromes=len(distance), all_pairs=checked)


def release(blocks, gates, parked, partner, path, rng):
    # Choose a mismatch away from both sources; its monitored material read
    # identifies an actual internal port and physical exterior partner.
    cut = edges(path)[1]
    p, b = gates[cut]
    q = partner[p]
    doubles = prepare(edges(path), gates, parked, partner)
    erased = doubles ^ {p, q}
    before, _ = read_ports(doubles, gates)
    after, _ = read_ports(erased, gates)
    if len(before) != 2 or len(after) != 4 or not set(before) <= set(after):
        raise ValueError('Release did not preserve both held sources')
    settings = {}
    # Fixed clock appointments use nominal registered contrast; other laws get
    # the identical appointments, with their own peak time as a cross-check.
    for name, (internal, exterior) in blocks.items():
        # Order of the three ports: p (gateway), q (same owner), b (across).
        # Their idle partners enforce total one double in p,q and b vacant.
        total = np.kron(internal, np.eye(3))
        ext = np.kron(exterior.reshape(3, 3, 3, 3), np.eye(3)).reshape(3, 3, 3, 3, 3, 3)
        total += ext.transpose(0, 4, 1, 2, 5, 3).reshape(27, 27)
        kept = [18, 6]  # ternary 200 and 020
        h = total[np.ix_(kept, kept)]
        splitting = np.ptp(np.linalg.eigvalsh(h))
        peak = np.pi / splitting if splitting else 1.
        times = [0., np.pi/65, 2*np.pi/65, float(peak)]
        probabilities = [float(abs(expm(-1j*t*h)[1, 0])**2) for t in times]
        counts = [int(rng.binomial(RELEASE_SHOTS, min(1., x))) for x in probabilities]
        settings[name] = dict(h=np.round(h, 10).tolist(), times=times,
                              extra_pair_probability=probabilities, extra_pair_counts=counts)
    return dict(ports_p_q_b=[p, q, b], before=before, after=after, shots=RELEASE_SHOTS,
                variants=settings)


def run():
    rng = np.random.default_rng(SEED)
    x, z = source_block()
    native = x + candidate.KAPPA*z
    blocks = dict(registered=(native, native), shuttered=(np.zeros_like(x), native),
                  exchange_only=(x, x), contrast_only=(candidate.KAPPA*z, candidate.KAPPA*z))
    calibration = fringe(float(native[6, 6]-native[0, 0]), rng)
    levels = []
    for level in (0, 1, 2):
        board = carrier.federation(level)
        rings, dual, gates, parked, partner = survey(board)
        distances = nx.single_source_shortest_path_length(dual, 0)
        rows = []
        for d in range(1, max(distances.values())+1):
            g = min(f for f, dist in distances.items() if dist == d)
            path = nx.shortest_path(dual, 0, g)
            by_trial = {}
            for name, moves in [('least', edges(path)), ('detour', edges(detour(path, dual))),
                                ('sources_only', source_only(0, g, dual))]:
                doubles = prepare(moves, gates, parked, partner)
                odd, mismatched = read_ports(doubles, gates)
                if not {0, g} <= set(odd):
                    raise ValueError('Lost endpoint')
                by_trial[name] = dict(odd_faces=odd, mismatch_edges=mismatched,
                                      **read_work(doubles, gates, native, rng))
            rows.append(dict(distance=d, endpoints=[0, g], trials=by_trial))
        levels.append(dict(level=level, carrier_sha256=carrier.seam_digest(board),
                           owners=board.cells, faces=len(rings), rows=rows,
                           release=release(blocks, gates, parked, partner, path, rng)))
        if level == 0:
            levels[-1]['independent_oracle'] = syndrome_oracle(dual)
    return dict(protocol_sha256=sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest(),
                law='interval_echo', kappa=candidate.KAPPA, seed=SEED, phase_shots=SHOTS,
                times=TIMES.tolist(), simultaneous_energy_half_width=PHASE_ERROR,
                calibration=dict(x_y=calibration, measured_work=decode(calibration)), levels=levels)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    report = run()
    path = HERE/'records.json'
    encoded = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.write:
        path.write_text(encoded)
    elif path.read_text() != encoded:
        raise SystemExit('Record replay changed')
    for level in report['levels']:
        far = level['rows'][-1]
        print(level['owners'], far['distance'],
              {n: v['measured_work'] for n, v in far['trials'].items()})


if __name__ == '__main__':
    main()
