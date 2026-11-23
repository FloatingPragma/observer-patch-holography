"""Measure paired, freely released load rays on two equally energetic native sources."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_null_balance.regional import verify
from model.experiments.material_inertia.dark_collar import label


HERE = Path(__file__).resolve().parent
TIMES = (.5, 1., 2., 4., 8.)
SHOTS = 2**29
SEED = 26112359
ARMS = {"sham": (0, 1), "coin_even": (2, 1), "coin_odd": (2, -1)}
TAGS = (1, 2)


def local_calibration():
    law = candidate.native_law
    exchange, cost = law.local_term_parts()
    native = exchange + candidate.KAPPA * cost
    record_exchange, record_cost = law.generator_parts()
    result, blanks = {}, {}
    for q in TAGS:
        indices = [law.SRC.index(word) for word in ((q, 0), (0, q))]
        cols = np.eye(law.DIM)[:, indices]
        reduced = cols.T @ native @ cols
        np.testing.assert_allclose(native @ cols, cols @ reduced, atol=3e-13, rtol=0)
        blanks[q] = reduced

    for arm, (mag, coin) in ARMS.items():
        columns = np.kron(label(mag, 1, coin)[:, None], np.eye(law.N_SRC))
        analyzer_columns = law.echo() @ columns
        readings = {}
        for name, operator in (("odd_reversal", record_exchange), ("contrast_squared", record_cost)):
            means = np.sum(analyzer_columns * (operator @ analyzer_columns), axis=0)
            # Sharp outcomes certify the separated source and contact-energy measurements.
            np.testing.assert_allclose(operator @ analyzer_columns, analyzer_columns * means,
                                       atol=3e-13, rtol=0)
            readings[name] = means.round(12)
        energy = readings["odd_reversal"] + candidate.KAPPA * readings["contrast_squared"]
        # Check all nine load words, including both orientations of the stationary source.
        projected = columns.T @ native @ columns
        np.testing.assert_allclose(native @ columns, columns @ projected, atol=3e-13, rtol=0)
        np.testing.assert_allclose(projected, np.diag(energy), atol=3e-13, rtol=0)
        physical = np.diag(projected)
        words = [law.SRC.index((1, q)) for q in (0, *TAGS)]
        result[arm] = {
            "odd_reversal": readings["odd_reversal"][words].tolist(),
            "contrast_squared": readings["contrast_squared"][words].tolist(),
            "source_energy": float(energy[words[0]]),
            "energy_with_calibration_tag": energy[words[1:]].tolist(),
            "measured_potential": (energy[words[1:]] - energy[words[0]]).tolist(),
            "native_potential": (physical[words[1:]] - physical[words[0]]).round(12).tolist(),
        }
    return blanks, result


def build_station(board, blanks, alternate_home=False):
    protected, parked = set(), []
    for cell in range(board.cells):
        local = set(range(12 * cell, 12 * (cell + 1)))
        gateways = {v for v in local if board.degree[v] == 6}
        p = min(local - gateways, key=lambda v: (-len(gateways.intersection(board.neighbours[v])), v))
        parked.append(p)
        neighbours = set(board.neighbours[p])
        protected.update(tuple(sorted((p, v))) for v in neighbours)
        remaining = local - neighbours - {p}
        available = sorted((u, v) for u in remaining for v in board.neighbours[u] if u < v and v in remaining)
        matchings = (edges for edges in itertools.combinations(available, 3)
                     if {v for edge in edges for v in edge} == remaining)
        matching = next(matchings)
        if cell == 0 and alternate_home:
            matching = next(matchings)
        protected.update(matching)

    mobile_mask = np.ones(board.n, dtype=bool)
    mobile_mask[parked] = False
    enabled = np.flatnonzero(mobile_mask)
    edges = np.asarray([edge for edge in board.seams if edge not in protected])
    ends = edges.ravel()
    degree = np.bincount(ends, minlength=board.n)
    np.testing.assert_array_equal(degree[enabled], board.degree[enabled] - 1)
    np.testing.assert_array_equal(degree[parked], 0)
    assert set(board.external).isdisjoint(protected)
    u, v = edges.T
    adjacency = sparse.csr_matrix((np.ones(ends.size), (np.r_[u, v], np.r_[v, u])),
                                  shape=(board.n, board.n))
    trimmed = adjacency[enabled][:, enabled]
    assert connected_components(trimmed, directed=False, return_labels=False) == 1
    free = (sparse.diags(degree[enabled].astype(float)) - trimmed)/2
    for q, block in blanks.items():
        diagonal = (np.bincount(u, weights=np.full(len(u), block[0, 0]), minlength=board.n)
                    + np.bincount(v, weights=np.full(len(v), block[1, 1]), minlength=board.n)
                    + candidate.KAPPA * (6 - board.degree) * q**2)
        native = sparse.diags(diagonal[enabled] - 5*candidate.KAPPA*q**2) + block[0, 1]*trimmed
        np.testing.assert_allclose((native - free).data, 0., atol=3e-13, rtol=0)

    # The cell distance counts working inter-owner links. Check there is a native route
    # within each cell between every pair of its eleven mobile ports.
    owners = enabled // 12
    for c in range(board.cells):
        local = np.flatnonzero(owners == c)
        assert connected_components(trimmed[local][:, local], directed=False, return_labels=False) == 1
    exterior = np.asarray(board.external) // 12
    a, c = exterior.T
    cell_graph = sparse.csr_matrix((np.ones(2*len(a)), (np.r_[a, c], np.r_[c, a])),
                                   shape=(board.cells, board.cells))
    dist = shortest_path(cell_graph, directed=False, unweighted=True).astype(int)
    return free.tocsr(), parked, enabled, protected, dist


def evolve(free, gate, prepared, potential):
    lens = sparse.csr_matrix(([potential], ([gate], [gate])), shape=free.shape)
    generator = -1j*(free + lens)
    waves = expm_multiply(generator, prepared, start=0., stop=8., num=17,
                          traceA=generator.diagonal().sum())
    result = np.abs(waves[[int(2*t) for t in TIMES]])**2
    np.testing.assert_allclose(result.sum(axis=1), 1., rtol=0, atol=3e-12)
    return result


def acquisition(board, station, meter, rng):
    free, parked, enabled, protected, dist = station
    inverse = dict(zip(enabled.tolist(), range(len(enabled)), strict=True))
    hosts = [0, *np.flatnonzero(dist[0] == 1)[:2].tolist()]
    placements = []
    distance_bins = dist.ravel()
    for host in hosts:
        p = parked[host]
        gate = min(i for i in board.neighbours[p] if board.degree[i] == 6)
        near = next(i for i in board.neighbours[gate] if i//12 != host)
        neighbours = [i for i in board.neighbours[near]
                      if i//12 == near//12 and tuple(sorted((i, near))) not in protected]
        launches = [near, *neighbours[:2]]
        prepared = np.zeros((len(enabled), 3))
        prepared[[inverse[i] for i in launches], np.arange(3)] = 1
        potentials = set(v for arm in meter.values() for v in arm["native_potential"])
        owner_laws = {}
        for v in potentials:
            port_laws = evolve(free, inverse[gate], prepared, v)
            owner_laws[v] = np.stack([[np.bincount(enabled//12, weights=law[:, j], minlength=board.cells)
                                       for j in range(3)] for law in port_laws])
        rows = []
        for arm, calibration in meter.items():
            for q in TAGS:
                # Independent echo energy is the forecast. The trajectories use the native compression.
                laws = owner_laws[calibration["native_potential"][q-1]]
                forecast = owner_laws[calibration["measured_potential"][q-1]]
                np.testing.assert_allclose(laws, forecast, atol=3e-13, rtol=0)
                for t, ports in zip(TIMES, laws, strict=True):
                    entry = float(ports[0, host])
                    rays = []
                    for side in ports[1:]:
                        hist = np.bincount(distance_bins, weights=np.outer(ports[0], side).ravel())
                        hist /= hist.sum()
                        counts = rng.multinomial(SHOTS, hist)
                        rays.append(dict(distance_counts=counts.tolist(),
                                         exact_different_owner=float(1 - hist[0]),
                                         exact_mean_squared_crossings=float(hist @ np.arange(len(hist))**2)))
                    rows.append(dict(arm=arm, tag=q, time=t, entry=int(rng.binomial(SHOTS, entry)),
                                     exact_entry=entry, rays=rays))
        placements.append(dict(cell=int(host), parked_port=int(p), gateway=int(gate),
                               launches=list(map(int, launches)), records=rows))
    return dict(cells=board.cells, mobile_ports=len(enabled), kept_external_links=len(board.external),
                guards=len(protected), placements=placements)


def run():
    rng = np.random.default_rng(SEED)
    blanks, meter = local_calibration()
    boards = []
    for level in (0, 1, 2):
        board = federation(level)
        boards.append(acquisition(board, build_station(board, blanks), meter, rng))
    total = sum(3*len(p["records"]) for b in boards for p in b["placements"])
    radius = float(np.sqrt(np.log(2*total/.01)/(2*SHOTS)))
    for b in boards:
        for p in b["placements"]:
            for r in p["records"]:
                errors = [abs(r["entry"]/SHOTS - r["exact_entry"])]
                errors += [abs(1 - s["distance_counts"][0]/SHOTS - s["exact_different_owner"])
                           for s in r["rays"]]
                if max(errors) > radius:
                    raise ValueError("outside the frozen simultaneous record envelope")
    return dict(law=candidate.NAME, protocol_sha256=hashlib.sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest(),
                seed=SEED, shots=SHOTS, simultaneous_confidence=.99, record_radius=radius,
                primary_frequencies=total, local_source_meter=meter, boards=boards)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")
    else:
        verify(json.loads((HERE/"results.json").read_text()), result)
    for b in result["boards"]:
        p = b["placements"][0]
        chosen = {r["arm"]: r for r in p["records"] if r["time"] == 2 and r["tag"] == 1}
        print(b["cells"], {a: (r["entry"]/SHOTS, 1-r["rays"][0]["distance_counts"][0]/SHOTS)
                           for a, r in chosen.items()})


if __name__ == "__main__":
    main()
