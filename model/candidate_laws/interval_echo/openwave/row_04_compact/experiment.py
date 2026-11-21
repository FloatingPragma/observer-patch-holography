"""Compact free modes and locally visible beats on the registered Interval Echo sphere."""

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, null_space
from scipy.sparse import eye
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation


HERE = Path(__file__).resolve().parent
TIMES = (0, 1, 4, 16, 64)
ATOL = 1e-9


def geometry():
    car = cell()
    distances = shortest_path(car.adjacency(), directed=False, unweighted=True)
    antipodes = [
        (i, int(np.argmax(row)))
        for i, row in enumerate(distances)
        if i < np.argmax(row)
    ]
    if len(antipodes) != 6 or np.any(np.sum(distances == 3, axis=1) != 1):
        raise ValueError("unexpected antipodal geometry")
    return car, antipodes


def gateway_table(car):
    out = [[] for _ in range(car.cells)]
    for edge in car.external:
        for p in edge:
            out[p // 12].append(p % 12)
    return [sorted(v) for v in out]


def port_reader(local, port):
    """Intersect two native unsigned-contrast reads at the addressed port."""
    law = candidate.native_law
    u = law.echo()
    tagged = np.repeat([False, *(abs(d) == 1 for _, d in law.LABELS)], len(law.SRC))
    columns = [law.SRC.index(pair) for pair in ((0, 0), (1, 0), (0, 1))]
    inputs = np.eye(law.DIM)[:, columns]
    recovered = u.T @ (tagged[:, None] * (u @ inputs))
    if not np.allclose(recovered, inputs * [0, 1, 1], rtol=0, atol=1e-12):
        raise ArithmeticError("an unsigned contrast probe fails the sharp subset read")
    neighbours = local.neighbours[port][:2]
    hits = np.ones(12)
    for other in neighbours:
        allowed = np.zeros(12)
        allowed[[port, other]] = 1
        hits *= allowed
    np.testing.assert_array_equal(hits, np.eye(12)[port])
    return neighbours, hits


def dark_bands(local_h, gateways):
    energies, vectors = np.linalg.eigh(local_h)
    bands = []
    for e in sorted(set(energies.round(9))):
        basis = vectors[:, abs(energies - e) < 1e-8]
        invisible = basis @ null_space(basis[gateways], rcond=1e-10)
        if invisible.shape[1]:
            bands.append(
                (
                    float(energies[abs(energies - e) < 1e-8].mean()),
                    invisible @ invisible.T,
                )
            )
    return bands


def prepare(local, gateways, target):
    """Pulse one sharp input. Every routing swap and partial contact uses the native seam block."""
    ports = list(np.flatnonzero(abs(target) > 1e-10))
    source = ports[0]
    if set(ports) & set(gateways):
        raise ValueError("degree-reference phase requires an interior port")
    term = candidate.native_law.seam_term(candidate.KAPPA)
    block = term[np.ix_([3, 1], [3, 1])]
    state = np.eye(12, dtype=complex)[:, source]
    swaps = beams = phases = 0
    tree = shortest_path(
        local.adjacency(),
        directed=False,
        unweighted=True,
        return_predecessors=True,
        indices=source,
    )[1]

    def pulse(a, b, angle):
        state[[a, b]] = expm(-1j * angle * block) @ state[[a, b]]

    for dest in ports[1:]:
        path = [dest]
        while path[-1] != source:
            path.append(int(tree[path[-1]]))
        # Move the target output next to the source, then restore every routing rail.
        links = list(zip(path[:-2], path[1:-1]))
        for a, b in links:
            pulse(a, b, np.pi)
        ratio = abs(target[dest]) / abs(state[source])
        if ratio > 1 + ATOL:
            raise ArithmeticError("impossible source split")
        pulse(source, path[-2], 2 * np.arcsin(min(ratio, 1)))
        for a, b in reversed(links):
            pulse(a, b, np.pi)
        swaps += 2 * len(links)
        beams += 1
    for p in ports:
        duration = (
            (np.angle(state[p]) - np.angle(target[p])) % (2 * np.pi) / candidate.KAPPA
        )
        state[p] *= np.exp(-1j * candidate.KAPPA * duration)
        phases += 1
    if np.linalg.norm(target - state) > ATOL:
        raise ArithmeticError("native preparation missed its specified target")
    return state, dict(
        initial_sharp_port=int(source),
        routed_pi_contacts=swaps,
        partial_contacts=beams,
        addressed_degree_phases=phases,
    )


def free(h, state, duration):
    if duration == 0:
        return state
    # A common phase cancels from every local occupation. This speeds up sparse evolution.
    offset = float(h.diagonal().mean())
    shifted = h - offset * eye(h.shape[0], format="csr")
    return expm_multiply(-1j * duration * shifted, state, traceA=0.0)


def audit_stencils(car, h, gateways, axes):
    counts = Counter()
    dark_dimensions = Counter()
    charts = {}
    highest_sharp_floor = 0.0
    local_h = candidate.generator(cell()).toarray()
    local_energy = None
    twice_h = (2 * h).copy()
    if np.max(abs(twice_h.data - np.rint(twice_h.data))) > 1e-12:
        raise ArithmeticError(
            "the registered matrix differs from its half-integer form"
        )
    twice_h.data = np.rint(twice_h.data).astype(np.int64)
    for home, tips in enumerate(gateways):
        chart = tuple(tips)
        if chart not in charts:
            charts[chart] = sum(p for _, p in dark_bands(local_h, tips))
        projector = charts[chart]
        dark_dimensions[int(round(np.trace(projector)))] += 1
        highest_sharp_floor = max(
            highest_sharp_floor, float(projector.diagonal().max())
        )
        available = [pair for pair in axes if not set(tips).intersection(pair)]
        counts[len(available) - 1] += 1
        section = twice_h[:, 12 * home : 12 * (home + 1)]
        for positive, negative in combinations(available, 2):
            v = np.zeros(12, dtype=int)
            v[list(positive)] = 1
            v[list(negative)] = -1
            full = np.zeros(car.n, dtype=int)
            full[12 * home : 12 * (home + 1)] = v
            integers = section @ v
            eigenvalue = integers[12 * home + positive[0]]
            if not np.array_equal(integers, eigenvalue * full):
                raise ArithmeticError("four-port stencil leaked or changed shape")
            if local_energy is not None and eigenvalue != local_energy:
                raise ArithmeticError("different compact energies")
            local_energy = int(eigenvalue)
    return dict(
        independent_four_port_modes_per_cell={
            str(k): v for k, v in sorted(counts.items())
        },
        all_dark_modes_per_cell={str(k): v for k, v in sorted(dark_dimensions.items())},
        best_all_time_home_floor_single_sharp_port=highest_sharp_floor,
        certified_cells=car.cells,
        energy=local_energy / 2,
    )


def release(car, local, h, axes, gateways):
    home = 0
    available = [pair for pair in axes if not set(gateways[home]).intersection(pair)]
    plus, minus = available[:2]
    v = np.zeros(12)
    v[list(plus)] = 0.5
    v[list(minus)] = -0.5
    signed, circuit = prepare(local, gateways[home], v)
    aligned, _ = prepare(local, gateways[home], abs(v))
    sharp = np.eye(12)[:, plus[0]]
    disturbed = np.sqrt(0.9) * signed + np.sqrt(0.1) * np.eye(12)[:, gateways[home][0]]
    states = np.zeros((car.n, 4), complex)
    states[:12] = np.stack((signed, aligned, sharp, disturbed), axis=1)
    bands = dark_bands(candidate.generator(local).toarray(), gateways[home])
    protected = sum(p for _, p in bands)
    floors = np.einsum("ij,ij->j", states[:12].conj(), protected @ states[:12]).real
    trials = []
    for t in TIMES:
        evolved = free(h, states, t)
        occupied = abs(evolved[:12]) ** 2
        inside = occupied.sum(axis=0)
        if (
            np.max(abs((abs(evolved) ** 2).sum(axis=0) - 1)) > ATOL
            or np.any(inside < floors - ATOL)
            or np.max(abs(occupied[:, 0] - v * v)) > ATOL
        ):
            raise AssertionError("home reader violated the predicted bound")
        trials.append(
            dict(
                t=t,
                home=inside.tolist(),
                signed_port_probabilities=occupied[:, 0].tolist(),
            )
        )
    return dict(
        home_cell=home,
        gateways=gateways[home],
        plus_pair=plus,
        minus_pair=minus,
        circuit=circuit,
        preparations=["signed", "all_positive", "single_port", "ten_percent_gateway"],
        all_time_home_lower_bounds=floors.tolist(),
        mean_dark_weight_over_sharp_ports=np.trace(protected) / 12,
        readings=trials,
    )


def beats(car, local, h, axes, gateways, rng):
    home = next(
        (
            i
            for i, tips in enumerate(gateways)
            if any(set(pair) <= set(tips) for pair in axes)
        ),
        None,
    )
    if home is None:
        return dict(eligible=False)
    indices = slice(12 * home, 12 * (home + 1))
    bands = dark_bands(candidate.generator(local).toarray(), gateways[home])
    if len(bands) < 2:
        raise AssertionError("two distinct energies required")
    e1, p = bands[0]
    e2, q = bands[-1]
    scores = p.diagonal() * q.diagonal()
    read = int(np.flatnonzero(scores >= scores.max() - 1e-10)[0])
    reader_neighbours, effect = port_reader(local, read)
    a, b = p[:, read] / np.sqrt(p[read, read]), q[:, read] / np.sqrt(q[read, read])
    source, circuit = prepare(local, gateways[home], (a + b) / np.sqrt(2))
    states = np.zeros((car.n, 3), complex)
    states[indices] = np.column_stack((source, a, b))
    period = 2 * np.pi / (e2 - e1)
    rows = []
    for cycles in (0, 0.5, 1, 16.5):
        evolved = free(h, states, cycles * period)
        if np.max(abs(1 - (abs(evolved[indices]) ** 2).sum(axis=0))) > ATOL:
            raise ArithmeticError("time-dependent excitation left its home cell")
        probabilities = effect @ (abs(evolved[indices]) ** 2)
        profile = (
            a[read] ** 2 / 2
            + b[read] ** 2 / 2
            + a[read] * b[read] * np.cos(2 * np.pi * cycles)
        )
        incoherent = probabilities[1:].mean()
        if (
            max(
                abs(probabilities[0] - profile),
                abs(incoherent - (a[read] ** 2 + b[read] ** 2) / 2),
            )
            > ATOL
        ):
            raise ArithmeticError("full free generator failed the local beat forecast")
        rows.append(
            dict(
                beat_periods=cycles,
                active_probability=probabilities[0],
                mixed_probability=incoherent,
                active_counts=int(rng.binomial(8192, np.clip(probabilities[0], 0, 1))),
                mixed_counts=int(rng.binomial(8192, np.clip(incoherent, 0, 1))),
            )
        )
    return dict(
        eligible=True,
        home_cell=home,
        gateways=gateways[home],
        read_port=read,
        incident_reader_neighbours=reader_neighbours,
        circuit=circuit,
        energies=[e1, e2],
        period=period,
        shots_per_setting=8192,
        readings=rows,
    )


def portable(obj):
    if isinstance(obj, dict):
        return {k: portable(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [portable(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return round(float(obj), 9) + 0.0
    return obj


def run():
    local, axes = geometry()
    rng = np.random.default_rng(731590)
    boards = []
    for level in range(3):
        car = federation(level)
        h = candidate.generator(car)
        gateways = gateway_table(car)
        boards.append(
            dict(
                cells=car.cells,
                external_contacts=len(car.external),
                audit=audit_stencils(car, h, gateways, axes),
                stationary=release(car, local, h, axes, gateways),
                beating=beats(car, local, h, axes, gateways, rng),
            )
        )
    return portable(
        dict(
            candidate=candidate.NAME,
            kappa=candidate.KAPPA,
            shot_seed=731590,
            boards=boards,
        )
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a fresh receipt")
    output = cli.parse_args().output
    result = run()
    if output:
        with output.open("x") as f:
            f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads((HERE / "results.json").read_text()):
        raise SystemExit("compact mode results differ from the archive")
    print(
        [
            (b["cells"], b["audit"]["energy"], b["beating"].get("period"))
            for b in result["boards"]
        ]
    )
