"""Independent reservoir, edge and drain records for native load transport."""

import argparse
from collections import Counter, deque
from fractions import Fraction
from math import gcd
import hashlib
import json
from pathlib import Path

import numpy as np
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ

from model.executable_model import federation
from model.executable_model.carrier import Carrier
from model.experiments.interval_echo_consensus.native_path_reference.occupation_action.measure import (
    native_kernel,
)
from model.experiments.interval_echo_consensus.refinement_cooling.dimer_band.measure import (
    pointer,
)

HERE = Path(__file__).resolve().parent
DEPTH, BURN, SHOTS, GRID = 6, 2048, 2**31, 2**144


def surveyed_board(board, full_transfer):
    # Each count is deterministic on a sharp launch at a full native contact.
    # Retain every offered address via its canonical enrollment and per-pattern
    # counts rather than storing thousands of identical binary rows.
    trials = 32
    accepted = []
    count_groups = Counter()
    for a, b in board.seams:
        reads = [int(trials * full_transfer[old][1 - old]) for old in (0, 1)]
        count_groups[tuple(reads)] += 1
        if reads == [trials, trials]:
            accepted.append((a, b))
    if len(count_groups) != 1:
        raise ArithmeticError("heterogeneous survey requires uncompressed edge records")
    key = json.dumps(board.seams, separators=(",", ":")).encode()
    receipt = dict(
        ordered_enrollment_sha256=hashlib.sha256(key).hexdigest(),
        trials_per_direction=trials,
        per_edge_forward_reverse_hits=list(next(iter(count_groups))),
        seams_offered=len(board.seams),
    )
    return Carrier(board.name, board.n, accepted), receipt


def distances(board, root=0, stop=DEPTH + 1):
    levels = {root: 0}
    pending = deque([root])
    while pending:
        i = pending.popleft()
        if levels[i] < stop:
            for j in board.neighbours[i]:
                if j not in levels:
                    levels[j] = levels[i] + 1
                    pending.append(j)
    return levels


def geometry(board):
    layer = distances(board)
    ports = sorted(layer)
    index = {p: i for i, p in enumerate(ports)}
    neighbors = [
        [index[q] for q in board.neighbours[p]] if layer[p] <= DEPTH else []
        for p in ports
    ]
    return ports, [layer[p] for p in ports], neighbors


def all_root_counts(board):
    counts = Counter()
    for root in range(board.n):
        layer = distances(board, root, 3)
        cuts = [
            sum(
                layer[q] > d
                for p in layer
                if layer[p] == d
                for q in board.neighbours[p]
            )
            for d in (1, 2)
        ]
        counts[tuple(cuts)] += 1
    return [[*pair, count] for pair, count in sorted(counts.items())]


def stationary(neighbors):
    interior = [v for v, ns in enumerate(neighbors) if ns]
    rank = {p: i for i, p in enumerate(interior)}
    lap = {}
    for i, v in enumerate(interior):
        lap[i] = {rank[w]: -1 for w in neighbors[v] if w in rank}
        lap[i][i] = len(neighbors[v])
    size = len(interior)
    mat = DomainMatrix(lap, (size, size), ZZ)
    rhs = DomainMatrix({rank[0]: {0: 1}}, (size, 1), ZZ)
    solution, denom = mat.solve_den(rhs)
    numerators = [int(x) for x in solution.to_Matrix()]
    common = gcd(int(denom), *numerators) * (1 if denom > 0 else -1)
    g = int(denom) // common
    z = [0] * len(neighbors)
    for i, v in enumerate(interior):
        z[v] = numerators[i] // common
    pi = [16 * value for value in z]
    for v in interior:
        for w in neighbors[v]:
            if not neighbors[w]:
                pi[w] += z[v]
        if sum(z[v] - z[w] for w in neighbors[v]) != g * (v == 0):
            raise ArithmeticError("stationary material source mismatch")
    norm = sum(pi)
    if sum(pi[v] for v, ns in enumerate(neighbors) if not ns) != g:
        raise ArithmeticError("drain balance mismatch")
    return dict(central=g, potential=z, stationary_numerators=pi, denominator=norm)


def after_burn(neighbors, mode):
    # Lower enclosures under the exact dyadic kernel. Lost weight is a joint
    # (one-sided) bound, including roundoff from all 2048 appointments.
    p = [GRID] + [0] * (len(neighbors) - 1)
    if mode == "parked":
        return p
    degrees = [len(ns) for ns in neighbors]
    for _ in range(BURN):
        mass = [0] * len(p)
        for v, ns in enumerate(neighbors):
            if ns:
                value = p[v]
                mass[v] += (16 - degrees[v]) * value
                for w in ns:
                    mass[w] += value
            else:
                dst = 0 if mode == "native" else v
                mass[dst] += 16 * p[v]
        p = [n // 16 for n in mass]
    if not 0 <= GRID - sum(p) <= len(p) * BURN:
        raise ArithmeticError("finite-time probability enclosure")
    return p


def arm(rng, weights, neighbors, mode, meter, source_hit, inner_drain=None):
    distribution = [Fraction(w, GRID) for w in weights]
    # Place all unallocated lower-interval mass on the last category; this
    # changes every event probability by at most the certified deficit.
    distribution[-1] += 1 - sum(distribution)
    counts = rng.multinomial(SHOTS, [float(q) for q in distribution])
    source = drain = inner = 0
    ledger = []
    for n, ns in zip(counts, neighbors, strict=True):
        if not ns:
            # Independently observe the requested launch at the central port.
            source += int(n * source_hit)
            ledger.append([int(n)])
            continue
        flip = int(mode != "parked")
        # First spare tags, then rejected and accepted contacts per neighbor.
        bins = [2 * (8 - len(ns))] + [v for _ in ns for v in (2 - flip, flip)]
        row = rng.multinomial(n, np.array(bins) / 16)
        for k, w in zip(row[2::2], ns, strict=True):
            if w == inner_drain:
                inner += int(k)
            elif not neighbors[w]:
                drain += int(k)
        ledger.append(row.tolist())
    if meter == "edge":
        return ledger
    return [{"source": source, "drain": drain, "inner": inner}[meter], SHOTS]


def burn_bound(cert, weights):
    # TV from the stationary law is bounded by the sum of positive deficits
    # relative to a componentwise lower enclosure of the finite-burn law.
    norm = cert["denominator"]
    overhang = sum(
        max(0, a * GRID - p * norm)
        for a, p in zip(cert["stationary_numerators"], weights, strict=True)
    )
    return Fraction(overhang, GRID * norm)


def panel(board, level, rng, full_transfer, parked_transfer):
    board, survey = surveyed_board(board, full_transfer)
    _, blocked_survey = surveyed_board(board, parked_transfer)
    ports, layer, neighbors = geometry(board)
    cert = stationary(neighbors)
    records = {}
    for mode in ("native", "parked", "source_off"):
        weights = after_burn(neighbors, mode)
        central_hit = (full_transfer if mode == "native" else parked_transfer)[0][1]
        record = {
            "burn_roundoff_units": GRID - sum(weights),
            "source_bank": arm(rng, weights, neighbors, mode, "source", central_hit),
            "edge_bank": arm(rng, weights, neighbors, mode, "edge", central_hit),
            "drain_bank": arm(rng, weights, neighbors, mode, "drain", central_hit),
        }
        if mode == "native":
            record["burn_tv_upper"] = str(burn_bound(cert, weights))
        records[mode] = record
    capture = layer.index(3)
    leaking = [ns.copy() for ns in neighbors]
    leaking[capture] = []
    leak_weights = after_burn(leaking, "native")
    leak_rng = np.random.default_rng(330215 + level)
    leak = {"capture_slot": capture}
    for meter in ("source", "edge", "drain", "inner"):
        leak[meter + "_bank"] = arm(
            leak_rng,
            leak_weights,
            leaking,
            "native",
            meter,
            full_transfer[0][1],
            capture,
        )
    return dict(
        level=level,
        cells=board.cells,
        full_contact_survey=survey,
        parked_survey=blocked_survey,
        interior_collector=leak,
        all_origins_first_two_cuts=all_root_counts(board),
        ports=ports,
        layers=layer,
        neighbor_slots=neighbors,
        certificate=cert,
        arms=records,
    )


def isolated_cell(rng):
    host = federation(0)
    edges = [(a, b) for a, b in host.seams if b < 12]
    board, survey = surveyed_board(
        Carrier("shuttered cell", 12, edges), native_kernel([np.pi] * 2)
    )
    layer = distances(board, stop=3)
    receivers = [p for p in range(12) if layer[p] == 3]
    if len(receivers) != 1:
        raise ValueError("antipodal collector is ambiguous")
    sink = receivers[0]
    foil = (0, board.neighbours[0][0])
    banks = []
    for defect in (False, True):
        ns = [
            [q for q in slots if not defect or tuple(sorted((p, q))) != foil]
            if p != sink
            else []
            for p, slots in enumerate(board.neighbours)
        ]
        cert = stationary(ns)
        weights = after_burn(ns, "native")
        rec = {
            meter + "_bank": arm(rng, weights, ns, "native", meter, 1)
            for meter in ("source", "edge", "drain")
        }
        banks.append(
            dict(
                ports=list(range(12)),
                layers=[layer[p] for p in range(12)],
                neighbor_slots=ns,
                certificate=cert,
                cut=list(foil) if defect else [],
                records=rec,
                burn_tv_upper=str(burn_bound(cert, weights)),
            )
        )
    return dict(survey=survey, banks=banks)


def run():
    pointer()
    for dt, exchange, flip in [
        (np.pi, True, 1),
        (np.pi / 2, True, Fraction(1, 2)),
        (np.pi / 2, False, 0),
    ]:
        kernel = native_kernel([dt, dt], exchange)
        if kernel[0][1] != flip or kernel[1][0] != flip:
            raise ArithmeticError("ranging or transport contrast")
    rng = np.random.default_rng(2026112531)
    return dict(
        antipodal=isolated_cell(np.random.default_rng(330219)),
        shots_per_bank=SHOTS,
        burn=BURN,
        depth=DEPTH,
        panels=[
            panel(
                federation(level),
                level,
                rng,
                native_kernel([np.pi] * 2),
                native_kernel([np.pi] * 2, False),
            )
            for level in range(3)
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=1) + "\n"
    path = HERE / "records.json"
    if args.write:
        path.write_text(output)
    elif path.read_text() != output:
        raise ValueError("shell-flux replay differs")
