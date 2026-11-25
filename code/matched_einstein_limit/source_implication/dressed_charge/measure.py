"""A compact exactly conserved dressing of the parked collar logarithm."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.sparse.linalg import expm_multiply

from ..full_algebra import native_precision, station
from . import interferometer as meter


def dark_cell(board, route):
    """Maximal subspace of this cell avoiding every gateway and marked contact."""
    cell = range(12)
    adjacency = sp.Matrix(12, 12, lambda i, j: int(j in board.neighbours[i]))
    excluded = sorted(
        set(p for p in cell if any(q not in cell for q in board.neighbours[p]))
        | (set(cell) & set(route))
    )
    pick = sp.eye(12)[excluded, :]
    invisible = sp.Matrix.vstack(*(pick * adjacency**j for j in cell))
    basis = sp.Matrix.hstack(*invisible.nullspace())
    assert basis.shape == (12, 2) and adjacency * basis == -basis
    projector = basis * (basis.T * basis).inv() * basis.T
    # With M = 3 projector, M_ii=1 at each of the six permitted launches.
    twice_m = 6 * projector
    assert all(a.q == 1 for a in twice_m)
    assert twice_m * twice_m == 6 * twice_m and twice_m.T == twice_m
    return (
        excluded,
        np.asarray(basis, dtype=np.int64),
        np.asarray(twice_m, dtype=np.int64),
    )


def initial_edge_numerator(h, twice_m, n):
    """Eight times the head block controlling d<J_4→5>/dt at a sharp launch from 5."""
    x = slice(20, 24)  # port 5
    y = slice(16, 20)  # port 4
    column = np.zeros((n, 1), dtype=np.int64)
    column[:12, 0] = twice_m[:, 5]
    loaded = np.kron(column, np.eye(4, dtype=np.int64))

    def skew(v):
        # With h=2H, this is 2(P4 H P5 - P5 H P4)v.
        out = np.zeros_like(v)
        out[y] = h[y, x] @ v[x]
        out[x] = -(h[x, y] @ v[y])
        return out

    return -(h @ skew(loaded) - skew(h @ loaded))[x]


def read_source():
    record = json.loads(Path(station.__file__).with_name("results.json").read_text())
    chi = np.array(record["x_real"]) + 1j * np.array(record["x_imag"])
    left = np.trace(chi.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)
    right = np.trace(chi.sum(axis=0).reshape(2, 2, 2, 2), axis1=0, axis2=2)
    product = np.stack([np.kron(block, right) for block in left])
    return chi, station.log_blocks(product) - station.log_blocks(chi)


def wave_at(h, port, time):
    start = np.zeros((h.shape[0], 4), dtype=complex)
    start[4 * port : 4 * port + 4] = np.eye(4)
    return expm_multiply(-0.5j * time * h, start).reshape(-1, 4, 4)


def read_flux(h, rho, k, m, wave, board):
    """Compute symmetrized charge densities, their time derivatives and edge currents."""
    n = len(wave)
    carried = np.zeros_like(wave)
    carried[:12] = np.einsum("pq,ab,qbi->pai", m, k, wave[:12])
    velocity = (-0.5j * (h @ wave.reshape(4 * n, 4))).reshape(n, 4, 4)
    carried_velocity = (-0.5j * (h @ carried.reshape(4 * n, 4))).reshape(n, 4, 4)

    def inner(a, b):
        return np.einsum("ij,...ai,...aj->...", rho.T, a.conj(), b, optimize=True).real

    density = inner(wave, carried)
    rate = inner(velocity, carried) + inner(wave, carried_velocity)
    currents = np.zeros((n, n))
    active = set(np.flatnonzero(np.diag(m)))
    edges = sorted({tuple(sorted((x, y))) for x in active for y in board.neighbours[x]})
    for x, y in edges:
        xy = h[4 * x : 4 * x + 4, 4 * y : 4 * y + 4].toarray() / 2
        # i (P_y H P_x - P_x H P_y) is the contribution into x from y.
        currents[y, x] = inner(wave[y], 1j * xy.T @ carried[x]) - inner(
            wave[x], 1j * xy @ carried[y]
        )
        currents[x, y] = -currents[y, x]
    bare = inner(wave, np.einsum("ab,pbi->pai", k, wave))
    return density, rate, currents, bare


def run():
    form = station.iso.pairing.formation
    board = form.federation(1)
    route = form.region(board)
    excluded, basis, twice_m = dark_cell(board, route)
    gram = basis.T @ basis
    active = np.flatnonzero(twice_m.diagonal())
    assert active.tolist() == [4, 5, 6, 7, 8, 11]
    np.testing.assert_array_equal(twice_m.diagonal()[active], 2)
    position_vectors = np.zeros((board.n, 2), dtype=np.int64)
    position_vectors[:12] = basis
    lifted = np.kron(position_vectors, np.eye(4, dtype=np.int64))
    spectral = {}
    hs = []
    for sign in (-1, 1):
        labels = []
        for odd in (False, True):
            h = native_precision.generator(board, route, sign, odd)
            idx = np.flatnonzero(lifted[:, 0])[0]
            eigen = int((h @ lifted[:, 0])[idx] // lifted[idx, 0])
            np.testing.assert_array_equal(h @ lifted, eigen * lifted)
            np.testing.assert_array_equal(
                initial_edge_numerator(h, twice_m, board.n),
                -3 * np.eye(4, dtype=np.int64),
            )
            labels.append(eigen)
            if sign == 1:
                hs.append(h)
        spectral[str(sign)] = labels
    chi, k = read_source()
    parked = float(station.pair(k, chi))
    spatial = twice_m[np.ix_(meter.PORTS, meter.PORTS)] / 2
    settings = meter.terminal_settings(meter.register(board), spatial)
    np.testing.assert_allclose(
        sum(a * b for a, b in settings), spatial, rtol=0, atol=3e-14
    )
    records = []
    tick = 1e-4

    def empty_readings():
        return [
            np.zeros(board.n),
            np.zeros(board.n),
            np.zeros((board.n, board.n)),
            np.zeros(board.n),
        ]

    for time in (0.5, 2.0):
        data = empty_readings()
        measured = 0.0
        side = np.zeros((2, board.n))
        for h, rho, charge in zip(hs, chi, k):
            for i, t in enumerate((time, time - tick, time + tick)):
                wave = wave_at(h, 5, t)
                row = read_flux(h, rho, charge, twice_m / 2, wave, board)
                if i == 0:
                    score = meter.spatial_score(wave, rho, charge)
                    measured += sum(a * np.sum(b * score).real for a, b in settings)
                    for out, contribution in zip(data, row):
                        out += contribution
                else:
                    side[i - 1] += row[0]
        density, rate, current, bare = data
        assert abs(density.sum() - parked) < 3e-12
        assert abs(measured - parked) < 3e-12
        assert np.max(abs(current.sum(axis=0) - rate)) < 1e-13
        diff = (side[1] - side[0]) / (2 * tick)
        err = np.max(abs(diff - rate))
        assert err < 1e-7
        local = current[np.ix_(active, active)]
        cut = current[np.ix_(np.setdiff1d(np.arange(board.n), active), active)].sum()
        pair = np.unravel_index(np.argmax(abs(local)), local.shape)
        records.append(
            dict(
                time=time,
                total=round(float(density.sum()), 14),
                terminal_read=float(measured),
                bare_head_change=float(bare.sum() - parked),
                six_port_population_charge=float(bare[active].sum()),
                six_port_coherence_charge=float(density.sum() - bare[active].sum()),
                density_on_six=density[active].tolist(),
                current_into_six=float(cut),
                maximal_internal_edge=[int(active[i]) for i in pair],
                signed_internal_current=float(local[pair]),
                finite_difference_error=float(err),
            )
        )
    return dict(
        excluded=excluded,
        active=active.tolist(),
        modes=basis.T.tolist(),
        gram=gram.tolist(),
        twice_restricted_dressing=twice_m[np.ix_(active, active)].tolist(),
        eigenvalues_of_integer_flight=spectral,
        ideal_projector_success="1/3",
        terminal_settings=len(settings),
        initial_current_rate_per_unit_head_charge="-3/8",
        parked_expectation=parked,
        records=records,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(out)
    else:
        print(out, end="")
