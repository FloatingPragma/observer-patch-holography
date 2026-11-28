"""Integer action of the pinned binary seam law on a compact pair and a mobile load."""

import argparse
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, diags

from model.candidate_laws.interval_echo.openwave.row_04 import open_pair, open_pair_release
from model.candidate_laws.interval_echo.openwave.row_04.open_pair_check import BOARD, BOARD_HASH


def add(target, key, value):
    target[key] += value
    if not target[key]:
        del target[key]


def twice_action(state, live, references, strength=8):
    out = Counter()
    for word, amplitude in state.items():
        walls = 0
        for port in word:
            for dest in live[port] - word:
                add(out, word ^ {port, dest}, -amplitude)
                walls += 1
        reference = sum(references[i] for i in word)
        add(out, word, ((1+2*strength)*walls + 2*strength*reference)*amplitude)
    return out


PAIR = {frozenset((156, 157)): 1, frozenset((156, 166)): -1}


def insert(state):
    """Attach the compact pair to a one-port amplitude, enforcing exclusion."""
    out = Counter()
    for single, value in state.items():
        for pair, sign in PAIR.items():
            if single.isdisjoint(pair):
                add(out, single | pair, sign*value)
    return out


def transition_probabilities(generator, start, terminals):
    """Rational Taylor coefficients of each transition, with a unitary remainder."""
    twice = generator*2
    if (twice != twice.T).nnz or not np.array_equal(twice.data, twice.data.astype(int)):
        raise ValueError("expected a real symmetric half-integer generator")
    # Subtract a scalar phase to reduce the Taylor remainder.
    twice = (twice-diags(np.full(twice.shape[0], 75.0))).tocsr()
    norm = F(int(abs(twice).sum(axis=1).max()), 2)
    rows = [tuple(zip(row.indices.tolist(), row.data.astype(int).tolist(), strict=True)) for row in twice]
    real = [F(0)]*len(terminals)
    imag = real.copy()
    power = [int(i == start) for i in range(len(rows))]
    denom, order, tail = 1, 0, F(1)
    while True:
        part = (real, imag)[order % 2]
        sign = (1, -1, -1, 1)[order % 4]
        for j, terminal in enumerate(terminals):
            part[j] += F(sign*power[terminal], denom)
        order += 1
        tail *= norm/order
        if tail < F(1, 10**13):
            break
        power = [sum(coef*power[k] for k, coef in row) for row in rows]
        denom *= 2*order
    # At t=1 the integral remainder has norm at most norm**order/order!.
    return [float(a*a+b*b) for a, b in zip(real, imag, strict=True)]


def contact_records(live, references, bridge, distance):
    """Born probabilities for a single sharp terminal question at the bridge."""
    ports = [p for p in range(len(live)) if p not in (156, 157, 166)]
    h = np.zeros((len(live), len(live)))
    for p in range(len(live)):
        for out, amplitude in twice_action({frozenset([p]): 1}, live, references).items():
            dest, = out
            h[dest, p] = amplitude/2
    compressed = csr_matrix(h[np.ix_(ports, ports)])
    contact = np.zeros(len(ports))
    contact[ports.index(bridge)] = -17
    arms = {"vacant_full_board": (csr_matrix(h), list(range(len(live)))),
            "vacant_exterior_compression": (compressed, ports),
            "compact_pair_free_release": (compressed + diags(contact), ports)}
    launches = [min(p for p in ports if distance[p] == d) for d in (0, 1, 6, 18)]
    prob = {}
    for name, (generator, addresses) in arms.items():
        # Real symmetry permits the same launch at the bridge for all four readouts.
        prob[name] = transition_probabilities(generator, addresses.index(bridge),
                                              [addresses.index(p) for p in launches])
    # Preparation is coherent. The trace distance to its compact component is
    # sqrt(1-F), giving an all-effects bound independent of flight dynamics.
    acceptance = 1-F(364, 729)**8
    prepared, setup = open_pair_release.source()
    target = np.array([-13., 0., 14.])
    if (not np.allclose(prepared, float(acceptance)*np.outer(target, target)/365, atol=2e-13, rtol=0)
            or abs(setup["exhaustion"]-float(1-acceptance)) > 2e-13):
        raise ArithmeticError("pinned monitored preparation disagrees with the accepted compact overlap")
    fidelity = F((13+14)**2, 2*365)
    delta = sqrt(1-fidelity)
    rows = []
    for i, port in enumerate(launches):
        p = prob["compact_pair_free_release"][i]
        lower = float(acceptance)*max(0, p-delta)
        upper = float(acceptance)*min(1, p+delta)
        rows.append({"launch": port, "distance": distance[port],
                     "bridge_probability": {k: round(float(v[i]), 10) for k, v in prob.items()},
                     "accepted_and_bridge_click_interval_full_preparation":
                         [max(0, np.floor(lower*1e8-1)*1e-8), min(1, np.ceil(upper*1e8+1)*1e-8)]})
    return rows


def measure():
    if sha256(BOARD.read_bytes()).hexdigest() != BOARD_HASH:
        raise ValueError("unexpected pinned board")
    with np.load(BOARD, allow_pickle=False) as archive:
        edges = [tuple(sorted(map(int, edge))) for edge in archive["edges"]]
    live = [set() for _ in range(240)]
    references = [6] * len(live)
    written = set(open_pair.geometry()[0])
    for a, b in edges:
        references[a] -= 1
        references[b] -= 1
        if (a, b) not in written:
            live[a].add(b)
            live[b].add(a)
    if len(written) != 11 or len(set(edges)) != len(edges):
        raise ValueError("unexpected written seam set or duplicate edge")
    site = 163
    distance = {site: 0}
    queue = deque([site])
    while queue:
        p = queue.popleft()
        for n in live[p] - distance.keys():
            distance[n] = distance[p]+1
            queue.append(n)
    if len(distance) != len(live):
        raise ArithmeticError("the written board is disconnected")
    panel = {}
    for strength in (8, 0):
        twice_energy = 2 + 6*strength
        if twice_action(PAIR, live, references, strength) != {p: twice_energy*a for p, a in PAIR.items()}:
            raise ArithmeticError("pair is not an exact eigenvector on the written board")
        per_shell = Counter()
        anomalies = {}
        for port in range(len(live)):
            if port in (156, 157, 166):
                continue
            single = {frozenset([port]): 1}
            joint = insert(single)
            residual = twice_action(joint, live, references, strength)
            for x, value in joint.items():
                add(residual, x, -twice_energy*value)
            # Compress the empty-board probe to the 237 exterior addresses,
            # retaining the physical diagonal at the bridge.
            free = {x: value for x, value in twice_action(single, live, references, strength).items()
                    if x.isdisjoint((156, 157, 166))}
            for x, value in insert(free).items():
                add(residual, x, -value)
            if residual:
                coefficient = -(2 + 4*strength)
                if residual != {x: coefficient*v for x, v in joint.items()}:
                    raise ArithmeticError("source attachment leaks or produces an off-diagonal interaction")
                anomalies[port] = coefficient
            else:
                per_shell[distance[port]] += 1
        if set(anomalies) != {site}:
            raise ArithmeticError("unexpected source/probe contact support")
        panel[str(strength)] = {"twice_pair_energy": twice_energy,
                                "ports_with_zero_residual_by_live_distance_to_bridge": dict(sorted(per_shell.items())),
                                "twice_contact_shift_at_bridge": anomalies[site]}
    return {"board_sha256": BOARD_HASH, "written_seams": sorted(written),
            "normalized_accepted_weight_in_compact_component": "729/730",
            "strengths": panel, "contact_detector_at_unit_exposure": contact_records(live, references, site, distance)}


if __name__ == "__main__":
    result = json.dumps(measure(), indent=2) + "\n"
    archive = Path(__file__).with_name("written_source.json")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        archive.write_text(result)
    elif result != archive.read_text():
        raise ArithmeticError("integer certificate differs from archive")
    print("Verified source/probe generator on every exterior port and at their common contact.")
