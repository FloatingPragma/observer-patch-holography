"""Physical cell-velocity measurements and shell occupations in the unchanged free evolution."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from .. import measure as front
from ...shared_radar.aperture.measure import source_and_changes
from ...material_clocks.measure import blank_clock

HERE = Path(__file__).parent
formation = front.optics.pairing.formation


def cell_distances(board):
    a, b = np.asarray(board.external).T // 12
    matrix = sparse.csr_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])),
        shape=(board.cells, board.cells),
    )
    assert np.max(matrix.data) == 1
    return shortest_path(matrix, directed=False, unweighted=True).astype(int)


def cross_block(board, h):
    co = h.tocoo()
    take = (co.row // 48 != co.col // 48) & (co.data != 0)
    return sparse.csr_matrix(
        (co.data[take], (co.row[take], co.col[take])), shape=h.shape
    )


def bridge_audit(board):
    ports = np.array(board.external)
    assert len(np.unique(ports)) == 2 * len(ports)
    a, b = ports.T
    raw = sparse.csr_matrix(
        (-np.ones(2 * len(a), dtype=int), (np.r_[a, b], np.r_[b, a])),
        shape=(board.n, board.n),
    )
    expected = sparse.kron(raw, np.eye(4, dtype=int), format="csr")
    route = formation.region(board)
    for q in (-1, 1):
        for odd in (False, True):
            h = front.native.generator(board, route, q, odd)
            assert (h - h.T).count_nonzero() == 0
            assert (cross_block(board, h) - expected).count_nonzero() == 0
    dist = cell_distances(board)
    # Every cell pair has an explicit optimal dual potential 2*d(origin,cell).
    lo, hi = ports.T // 12
    assert np.max(abs(dist[:, lo] - dist[:, hi])) == 1
    indicators = np.arange(board.cells)
    assert np.all(dist[indicators, indicators] == 0)
    return dict(
        cells=board.cells,
        bridge_count=len(a),
        port_count=board.n,
        cell_diameter=int(dist.max()),
        velocity_dual_diameter=int(2 * dist.max()),
    )


def probability_bound(distance, time):
    # Optimized squared weighted-propagator bound. Rate is per native clock unit.
    if time == 0:
        return float(distance == 0)
    v = distance / time
    if v <= 0.5:
        return 1.0
    return float(
        np.exp(-2 * time * (v * np.arccosh(2 * v) - 0.5 * np.sqrt(4 * v * v - 1)))
    )


def launch(board, ports, weights):
    x = np.zeros((4 * board.n, 4), dtype=complex)
    for p, w in zip(ports, weights):
        x[4 * p : 4 * p + 4] = w * np.eye(4)
    return x


def flight(h, initial, time):
    ih = (-0.5j * time) * h
    return expm_multiply(ih, initial, traceA=ih.diagonal().sum()).reshape(-1, 48, 4)


def data(board, chi, tangents, rng):
    route = formation.region(board)
    a, b = route[1:3]
    source_cell = a // 12
    dist = cell_distances(board)[source_cell]
    t_short = 1 / 2048
    shots = 2**44
    # Four signed preparations for each charge and one independent clock.
    tail = np.log(18 / 0.01)
    clock_shots = 2**52
    clock_p = abs(expm(-1j * t_short * blank_clock(1, 8))[1, 0]) ** 2
    clock_hits = int(rng.binomial(clock_shots, clock_p))
    clock_err = (np.sqrt(2 * clock_hits * tail) + 2 * tail) / clock_shots

    def clock_time(probability):
        return 2 * np.arcsin(np.sqrt(probability))

    dwell = clock_time(clock_hits / clock_shots)
    low, high = [clock_time(clock_hits / clock_shots + s * clock_err) for s in (-1, 1)]
    assert 0 < low < t_short < high < np.pi / 2
    clock = dict(
        counts=[clock_hits, clock_shots - clock_hits],
        dwell_from_clicks=dwell,
        simultaneous_dwell_interval=[low, high],
        assumed_branch="0 < t < pi/2",
    )
    out = {}
    for q in (-1, 1):
        local = {
            variant: {label: 0.0 for label in ["plus", "minus"]}
            for variant in ("active", "half_bridge_hop")
        }
        terminal = {t: [] for t in [2.0, 4.0, 8.0]}
        generators = []
        for parity in (False, True):
            h = front.native.generator(board, route, q, parity)
            generators.append(h)
            rho = chi[int(parity)]
            for variant, generator in [
                ("active", h),
                ("half_bridge_hop", h - 0.5 * cross_block(board, h)),
            ]:
                for sign, label in [(1, "plus"), (-1, "minus")]:
                    wave = flight(
                        generator,
                        launch(board, [a, b], [1 / np.sqrt(2), sign * 1j / np.sqrt(2)]),
                        t_short,
                    )[b // 12]
                    local[variant][label] += np.trace(wave.conj().T @ wave @ rho).real
            x = launch(board, [route[0]], [1])
            for t, rows in terminal.items():
                w = flight(h, x, t)
                rows.append(w.conj().swapaxes(-1, -2) @ w)
        hnorm = max(float(abs(h).sum(axis=1).max() / 2) for h in generators)
        assert hnorm <= 32
        observations = {}
        for variant, probs in local.items():
            counts = {k: int(rng.binomial(shots, p)) for k, p in probs.items()}
            difference = (counts["plus"] - counts["minus"]) / (2 * shots)
            measured = difference / dwell
            p_width = np.sqrt(tail / (2 * shots))
            # [H,[H,F]] has norm <= 2*||H||*||[H,F]|| <= ||H||.
            truncation = hnorm * high / 2
            interval = [
                (difference - p_width) / high - truncation,
                (difference + p_width) / low + truncation,
            ]
            expected = 0.5 if variant == "active" else 0.25
            assert interval[0] < expected < interval[1]
            observations[variant] = dict(
                oriented_click_probability=probs,
                counts={k: [n, shots - n] for k, n in counts.items()},
                measured_speed=measured,
                speed_interval_99pct=interval,
                forward_difference_truncation_bound=truncation,
            )
        controls = {}
        for t, rows in terminal.items():
            op = np.array(rows)
            p = np.einsum("bcij,bji->c", op, chi).real
            dp = {
                key: np.einsum("bcij,bji->c", op, tangents[name]).real
                for name, key in [("ordinary", "ordinary"), ("unpairing", "pairing")]
            }
            np.testing.assert_allclose(p.sum(), 1, atol=1e-13)
            shell = {}
            for r in [1, 2, 3, 4]:
                take = dist >= r
                val = p[take].sum()
                upper = probability_bound(r, t)
                assert 0 < val <= upper + 1e-14
                shell[str(r)] = dict(
                    population=val,
                    cone_upper=upper,
                    matched_log_slopes={
                        k: float(v[take].sum() / val) for k, v in dp.items()
                    },
                )
            controls[str(t)] = shell
        out[str(q)] = dict(velocity=observations, shells=controls)
    return dict(
        bridge=[a, b],
        short_native_dwell=t_short,
        dimer_clock=clock,
        sharp_launch=route[0],
        charges=out,
    )


def run():
    boards = [formation.federation(level) for level in (0, 1, 2)]
    chi, tangents, *_ = source_and_changes()
    return dict(
        hardware=[bridge_audit(b) for b in boards],
        tested_station=data(boards[1], chi, tangents, np.random.default_rng(194226)),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")
