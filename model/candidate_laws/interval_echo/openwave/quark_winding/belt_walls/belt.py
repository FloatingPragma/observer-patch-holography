"""Material transport and winding of period-three defects on the native ten-port belt."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier


TIMES = (0.0, 1 / 64, 1 / 16, 1 / 4, 1.0, 4.0, 16.0)


def belts():
    car = carrier.cell()
    g = nx.Graph(car.seams)
    result = {}
    for a, b in combinations(range(12), 2):
        if nx.shortest_path_length(g, a, b) != 3:
            continue
        interior = g.subgraph(set(g) - {a, b})
        north = set(car.neighbours[a])
        perimeter = nx.Graph(
            (u, v) for u, v in interior.edges if (u in north) != (v in north)
        )
        (cycle,) = nx.cycle_basis(perimeter)
        order = min(
            tuple(seq[k:] + seq[:k])
            for seq in (cycle, cycle[::-1])
            for k in range(len(seq))
        )
        expected = {
            tuple(sorted((order[i], order[(i + d) % 10])))
            for i in range(10)
            for d in (1, 2)
        }
        assert (
            expected == {tuple(sorted(e)) for e in interior.edges} and len(order) == 10
        )
        result[a, b] = order
    assert len(result) == 6
    return result


def ring(length):
    seams = sorted(
        {tuple(sorted((i, (i + d) % length))) for i in range(length) for d in (1, 2)}
    )
    return carrier.Carrier("squared cycle", length, seams)


def signed_distance(x, origin, length):
    return (np.asarray(x) - origin + length / 2) % length - length / 2


def selected(sector):
    return np.flatnonzero(sector.agreed == sector.agreed.min())


def wall(sites, length):
    gaps = [(a, (b - a) % length) for a, b in zip(sites, (*sites[1:], sites[0]))]
    ((a, gap),) = [(a, gap) for a, gap in gaps if gap != 3]
    return a, gap, (a + gap / 2) % length


def static_read(sites, length):
    start, gap, centre = wall(sites, length)
    n = np.isin(np.arange(length), sites).astype(int)

    # Clocks at both sides use the same coordinate unwrapped from the left occupied port.
    def clock(offset):
        x = np.arange(offset, offset + 3)
        return [
            int(sum(n[(start + j) % length] for j in x if j % 3 == i)) for i in range(3)
        ]

    left, right = clock(-2), clock(gap)
    w = np.roll(n, 1) + n + np.roll(n, -1) - 1
    unit = np.zeros(length, int)
    unit[int(np.floor(centre))] = 1
    injection = np.roll(unit, 1) + unit + np.roll(unit, -1)
    cuts = []
    for radius in (2, 3):
        inside = abs(signed_distance(np.arange(length), centre, length)) <= radius
        cuts.append(
            {
                "radius": radius,
                "sum_above_bulk": int(sum(w[inside])),
                "unit_calibration": int(sum(injection[inside])),
            }
        )
    return dict(
        gap=gap, windows=w.tolist(), left_clock=left, right_clock=right, cuts=cuts
    )


def bands(sector, kappa):
    h = sector.hamiltonian(kappa)
    energy, vectors = np.linalg.eigh(h)
    basis = selected(sector)
    d = len(basis)
    top = vectors[:, -d:]
    least = np.linalg.svd(top[basis], compute_uv=False)[-1] ** 2
    return (
        energy,
        vectors,
        dict(
            selected_count=d,
            minimum_occupied_seams=int(sector.agreed.min()),
            gap_below_top_group=energy[-d] - energy[-d - 1],
            group_width=energy[-1] - energy[-d],
            smallest_squared_overlap=least,
            every_time_retention_lower=max(0.0, 1 - 4 * (1 - least)),
            bottom_state_selected_weight=sum(abs(vectors[basis, 0]) ** 2),
        ),
    )


def transport(sector, kappa):
    length = sector.car.n
    energy, vectors, spectral = bands(sector, kappa)
    basis = selected(sector)
    geometry = {i: wall(sector.states[i], length) for i in basis}
    (start,) = [i for i in basis if geometry[i][0] == 0]
    transitions = []
    for i, j in combinations(basis, 2):
        if sector.K[j, i] == 0:
            continue
        a, b = set(sector.states[i]), set(sector.states[j])
        (lost,) = a - b
        (gained,) = b - a
        dw = signed_distance(geometry[j][2], geometry[i][2], length)
        dm = signed_distance(gained, lost, length)
        transitions.append((int(i), int(j), float(dm), float(dw)))
    (dest,) = [
        i
        for i in basis
        if sector.K[i, start] != 0
        and signed_distance(geometry[i][2], geometry[start][2], length) == 3
    ]
    (old,) = set(sector.states[start]) - set(sector.states[dest])
    (new,) = set(sector.states[dest]) - set(sector.states[start])
    dm = float(signed_distance(new, old, length))
    mw = signed_distance(np.arange(length), old + dm / 2, length)
    ww = signed_distance(
        [geometry[i][2] for i in basis], geometry[start][2] + 1.5, length
    )
    occ = np.array([[int(x in s) for x in range(length)] for s in sector.states])
    inputs = np.zeros((len(energy), 2), complex)
    inputs[start, :] = 1 / np.sqrt(2)
    inputs[dest, :] = np.array([1j, -1j]) / np.sqrt(2)
    modes = vectors.T @ inputs
    panels = []
    for t in TIMES:
        prob = abs(vectors @ (np.exp(-1j * t * energy)[:, None] * modes)) ** 2
        ports = occ.T @ prob
        wall_prob = prob[basis]
        material = float(mw @ (ports[:, 0] - ports[:, 1]) / 2)
        travel = float(ww @ (wall_prob[:, 0] - wall_prob[:, 1]) / 2)
        panels.append(
            dict(
                time=t,
                port_probabilities=ports.T.tolist(),
                wall_probabilities=wall_prob.T.tolist(),
                outside=(1 - wall_prob.sum(axis=0)).tolist(),
                signed_material_displacement=material,
                signed_wall_displacement=travel,
                ratio=material / travel if abs(travel) > 1e-9 else None,
            )
        )
    return dict(
        length=length,
        loads=sector.N,
        kappa=kappa,
        spectral=spectral,
        wall_addresses=[list(sector.states[i]) for i in basis],
        centres=[geometry[i][2] for i in basis],
        launch=[list(sector.states[i]) for i in (start, dest)],
        native_hops=transitions,
        static=static_read(sector.states[start], length),
        launch_material_step=dm,
        launch_wall_step=3,
        reads=panels,
    )


def opened(order, guided):
    sector = BinarySector(carrier.cell(), guided["loads"])
    h = sector.hamiltonian(candidate.KAPPA)
    energy, vectors = np.linalg.eigh(h)

    def idx(configuration):
        return sector.index[tuple(sorted(order[x] for x in configuration))]

    allowed = [idx(s) for s in guided["wall_addresses"]]
    launch = [idx(s) for s in guided["launch"]]
    initial = np.zeros(len(energy), complex)
    initial[launch] = [1 / np.sqrt(2), 1j / np.sqrt(2)]
    modes = vectors.T @ initial
    clear = [i for i, s in enumerate(sector.states) if set(s).issubset(order)]
    return [
        dict(time=t, original_wall=sum(p[allowed]), poles_empty=sum(p[clear]))
        for t in TIMES
        for p in [abs(vectors @ (modes * np.exp(-1j * t * energy))) ** 2]
    ]


def run():
    native = belts()
    order = next(iter(native.values()))
    result = {}
    for length in (10, 13):
        for n in (length // 3, length // 3 + 1):
            sector = BinarySector(ring(length), n)
            for kappa in (candidate.KAPPA, 0.0):
                label = f"L{length}_N{n}_k{int(kappa)}"
                result[label] = transport(sector, kappa)
                if length == 10 and kappa:
                    result[label]["opened_cell"] = opened(order, result[label])
    reference = BinarySector(ring(9), 3)
    return currents._portable(
        dict(
            axes=[dict(poles=k, belt=v) for k, v in native.items()],
            experiments=result,
            reference=dict(
                words=[reference.states[i] for i in selected(reference)],
                spectral=bands(reference, candidate.KAPPA)[2],
            ),
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = run()
    with args.output.open("x") as fp:
        json.dump(data, fp, indent=2, allow_nan=False)
        fp.write("\n")
