"""Two-mode Ramsey clocks running freely outside the two stored heads."""

import argparse
import json
from fractions import Fraction
from functools import lru_cache
from math import lcm
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as sp
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import (
    compile_wave,
)
from model.executable_model import cell, federation
from model.experiments.interval_echo_space.shell_source import prepare
from ..full_algebra import native_precision as native
from .measure import density_center

HERE = Path(__file__).resolve().parent


@lru_cache(None)
def modes(gateways):
    adjacency = sp.Matrix(cell().adjacency().toarray().astype(int))
    guard = sp.eye(12)[list(gateways), :]
    # A rational vector and its neighbour sum must both vanish at every exit.
    equations = (
        (adjacency**2 - 5 * sp.eye(12)).col_join(guard).col_join(guard * adjacency)
    )
    basis = equations.nullspace()
    if not basis:
        return None
    vector = basis[0] * lcm(*(int(x.q) for x in basis[0]))
    adjacent = adjacency * vector
    assert (
        adjacency * adjacent == 5 * vector
        and guard * vector == guard * adjacent == sp.zeros(len(gateways), 1)
    )
    return tuple(np.array(x, dtype=np.int64).ravel() for x in (vector, adjacent))


def clock_cells(board, route):
    sources = sorted({p // 12 for p in route})
    a, b = (np.asarray(board.external) // 12).T
    network = sparse.csr_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])),
        shape=(board.cells, board.cells),
    )
    ranges = (
        shortest_path(network, directed=False, unweighted=True, indices=sources)
        .min(axis=0)
        .astype(int)
    )
    clocks = {}
    for i, degree in enumerate(board.degree.reshape(-1, 12)):
        if i in sources:
            continue
        basis = modes(tuple(np.flatnonzero(degree == 6)))
        if basis:
            clocks[i] = basis
    return clocks, {i: int(ranges[i]) for i in clocks}


def verify_columns(h, index, vectors):
    # All four coin columns, and every output including other cells, are checked.
    ids = np.arange(48 * index, 48 * (index + 1))
    local = h[:, ids]
    vec, adj = (np.kron(v[:, None], np.eye(4, dtype=np.int64)) for v in vectors)
    # The spectral check below uses signed int64 arithmetic.
    assert (
        int(abs(local).sum(axis=1).max())
        * max(int(abs(vec).max()), int(abs(adj).max()))
        < 2**60
    )
    slot = np.flatnonzero(vectors[0])[0]
    center = int(local[ids[4 * slot], 4 * slot])
    for x, ax in ((vec, adj), (adj, 5 * vec)):
        expected = np.zeros((h.shape[0], 4), dtype=np.int64)
        expected[ids] = center * x - ax
        np.testing.assert_array_equal(local @ x, expected)


def ramsey_modes(vectors):
    v, av = vectors
    pair = [av + s * np.sqrt(5) * v for s in (1, -1)]
    lo, hi = (a / np.linalg.norm(a) for a in pair)
    assert abs(lo @ hi) < 1e-13
    return lo, hi


def experiment(h, clocks, ranges):
    first = min(clocks, key=lambda i: (ranges[i], i))
    last = max(clocks, key=lambda i: (ranges[i], -i))
    assert ranges[first] < ranges[last]
    rng = np.random.default_rng(162753)
    samples = []
    shots = 2**18
    for i in (first, last):
        low, high = ramsey_modes(clocks[i])
        initial = compile_wave(cell(), (low + high) / np.sqrt(2))
        analyzers = np.stack(
            [
                compile_wave(cell(), (low + phase * high) / np.sqrt(2))
                for phase in (1.0, -1j)
            ]
        )
        # Full native propagation, independent of a two-state fit. The source
        # block and coin column are arbitrary by the exact operator checks.
        local = np.arange(48 * i, 48 * (i + 1), 4)
        w = np.zeros(h.shape[0], dtype=complex)
        w[local] = initial
        for factor in (1.0, 1.25) if i == first else (1.0,):
            # Calibrate sensitivity by changing internal exchange timing on this
            # clock alone. Add (factor-1) times twice its internal exchange.
            perturb = sparse.kron(
                5 * sparse.eye(12) - cell().adjacency(), sparse.eye(4)
            ).tocoo()
            embedded = sparse.csr_matrix(
                (
                    perturb.data * (factor - 1),
                    (48 * i + perturb.row, 48 * i + perturb.col),
                ),
                shape=h.shape,
            )
            generator = -0.5j * (h + embedded)
            times = np.arange(5) * 0.35
            wave = expm_multiply(
                generator,
                w,
                start=0,
                stop=times[-1],
                num=len(times),
                traceA=generator.diagonal().sum(),
            )
            probability = np.abs(wave[:, local] @ analyzers.conj().T) ** 2
            phase = np.unwrap(
                np.arctan2(2 * probability[:, 1] - 1, 2 * probability[:, 0] - 1)
            )
            rate = float((times @ phase) / (times @ times))
            np.testing.assert_allclose(phase, np.sqrt(5) * times * factor, atol=2e-13)
            counts = rng.binomial(shots, np.clip(probability, 0, 1))
            counted_phase = np.unwrap(
                np.arctan2(2 * counts[:, 1] / shots - 1, 2 * counts[:, 0] / shots - 1)
            )
            samples.append(
                dict(
                    cell=i,
                    relay_hops=ranges[i],
                    internal_rate=factor,
                    shots_per_setting=shots,
                    times=times.tolist(),
                    counts=counts.tolist(),
                    frequency=rate,
                    sampled_frequency=float((times @ counted_phase) / (times @ times)),
                )
            )
    return samples


def compare_at_one_receiver(board, route, h, clocks, ranges):
    # Unitary splitters before and after free flight run on blank contacts.
    first = min(clocks, key=lambda i: (ranges[i], i))
    last = max(clocks, key=lambda i: (ranges[i], -i))
    marked = {tuple(sorted(pair)) for pair in (route[:2], route[2:])}
    apparatus = SimpleNamespace(car=board, active=set(board.seams) - marked)
    # Every setting returns to the same material pointer by inverse preparation.
    root = 12 * first + int(np.flatnonzero(clocks[first][0])[0])
    times = np.arange(6) * 0.2
    shape = h.shape
    exchange = sparse.kron(
        5 * sparse.eye(12) - cell().adjacency(), sparse.eye(4)
    ).tocoo()
    perturb = sparse.csr_matrix(
        (exchange.data / 4, (48 * first + exchange.row, 48 * first + exchange.col)),
        shape=shape,
    )
    modes_at = {i: ramsey_modes(clocks[i]) for i in (first, last)}
    rng = np.random.default_rng(162754)
    shots = 2**19
    experiments = []
    for change in (False, True):
        phases, counts = [], []
        counted_phases = []
        for band in range(2):
            probes = []
            for phase in (1.0, -1j):
                target = np.zeros(board.n, dtype=complex)
                for i, weight in ((first, 1.0), (last, phase)):
                    target[12 * i : 12 * (i + 1)] = (
                        weight * modes_at[i][band] / np.sqrt(2)
                    )
                assert all(board.degree[q] == 5 for q in np.flatnonzero(target))
                vector, _ = prepare(apparatus, target, root=root)
                probes.append(vector)
            w = np.zeros(shape[0], dtype=complex)
            w[::4] = probes[0]
            generator = -0.5j * (h + perturb if change else h)
            flight = expm_multiply(
                generator,
                w,
                start=0,
                stop=times[-1],
                num=len(times),
                traceA=generator.diagonal().sum(),
            )
            prob = np.abs(flight[:, ::4] @ np.stack(probes).conj().T) ** 2
            recorded = rng.binomial(shots, np.clip(prob, 0, 1))
            phases.append(np.unwrap(np.arctan2(2 * prob[:, 1] - 1, 2 * prob[:, 0] - 1)))
            counted_phases.append(
                np.unwrap(
                    np.arctan2(
                        2 * recorded[:, 1] / shots - 1, 2 * recorded[:, 0] / shots - 1
                    )
                )
            )
            counts.append(recorded.tolist())
        difference = phases[1] - phases[0]
        slope = float((times @ difference) / (times @ times))
        expected = -np.sqrt(5) / 4 if change else 0.0
        np.testing.assert_allclose(difference, expected * times, atol=3e-13)
        counted = counted_phases[1] - counted_phases[0]
        experiments.append(
            dict(
                near_exchange_factor=1.25 if change else 1.0,
                counts=counts,
                differential_angular_frequency=slope,
                sampled_differential_frequency=float(
                    (times @ counted) / (times @ times)
                ),
            )
        )
    # Hoeffding for all 48 separate probability settings, then the angle
    # subtended by a disk around an exact unit quadrature vector.
    settings = len(experiments) * 2 * 2 * len(times)
    probability_error = np.sqrt(np.log(2 * settings / 0.01) / (2 * shots))
    angle_error = np.arcsin(2 * np.sqrt(2) * probability_error)
    fit_error = float(2 * angle_error * times.sum() / (times @ times))
    assert all(
        abs(
            row["sampled_differential_frequency"]
            - row["differential_angular_frequency"]
        )
        <= fit_error
        for row in experiments
    )
    return dict(
        near_cell=first,
        far_cell=last,
        pointer_port=root,
        shots=shots,
        simultaneous_99_percent_frequency_half_width=fit_error,
        times=times.tolist(),
        experiments=experiments,
    )


def report():
    rows = []
    for level in range(3):
        board = federation(level)
        route = native.iso.pairing.formation.region(board)
        clocks, ranges = clock_cells(board, route)
        checks = 0
        for charge in (1, -1):
            for odd in (False, True):
                h = native.generator(board, route, charge, odd)
                for i, pair in clocks.items():
                    verify_columns(h, i, pair)
                    checks += 1
                if level == 1 and charge == 1 and not odd:
                    measured = experiment(h, clocks, ranges)
                    joint = compare_at_one_receiver(board, route, h, clocks, ranges)
        rows.append(
            dict(
                cells=board.cells,
                stored_heads=list(route),
                clock_count=len(clocks),
                by_relay_distance={
                    str(k): list(ranges.values()).count(k)
                    for k in sorted(set(ranges.values()))
                },
                verified_clock_sectors=checks,
            )
        )
    _, p, _, radius = density_center()
    assert radius <= Fraction(2, 10**7)
    # The trace of the four-dimensional odd block is 2-Lipschitz in HS norm.
    lo, hi = (float(Fraction(1, 2) - p + sign * 2 * radius) for sign in (-1, 1))
    return dict(
        boards=rows,
        equilibrium_board=80,
        ordinary_parked_energy_slope=[lo, hi],
        unpair_parked_energy_slope=0,
        common_angular_frequency="sqrt(5)",
        fractional_frequency_response_ordinary=0,
        fractional_frequency_response_unpair=0,
        material_ramsey=measured,
        one_receiver=joint,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "exterior.json").write_text(output)
    else:
        print(output, end="")
