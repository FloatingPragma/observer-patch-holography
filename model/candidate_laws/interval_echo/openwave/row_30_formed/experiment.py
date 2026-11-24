"""Local empty repairs, their band count, and a joint spatial and temporal test."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import (
    breadth_first_order,
    connected_components,
    shortest_path,
)

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
K = candidate.KAPPA
Q = 256
PERIOD = 64
BINS = np.arange(Q) * PERIOD / Q
SHOTS = 2**20
SPACE_SHOTS = 2**34
# Six boards, six spatial preparations, and all pair, gate and return readings.
COUNT_EPS = np.sqrt(np.log(2 * 6 * 6 * (Q + 5) / 0.01) / (2 * SHOTS))
BANKS = (178121, 178122)
TIMES = (0, 8, 64, 512)


def native():
    law = candidate.native_law
    blank = np.eye(law.DIM)[:, [law.SRC.index(x) for x in ((0, 0), (1, 0), (0, 1))]]
    vacuum = law.echo() @ blank[:, 0]
    written = np.zeros_like(blank)
    for j, pos in enumerate(np.argmax(blank, axis=0)):
        written[pos :: law.N_SRC, j] = vacuum[:: law.N_SRC]
    exchange, contrast = law.local_term_parts()
    block = np.zeros((3, 3))
    block[1:, 1:] = np.array([[1, -1], [-1, 1]]) / 2
    for h, sub in ((exchange, block), (contrast, np.diag([0, 1, 1]))):
        np.testing.assert_allclose(h @ blank, blank @ sub, atol=2e-14)
        np.testing.assert_allclose(h @ written, 0, atol=2e-14)
    np.testing.assert_allclose(law.echo() @ vacuum, blank[:, 0], atol=2e-14)
    half = expm(-1j * np.pi / 2 * (block[1:, 1:] + K * np.eye(2)))
    p = float(abs(half[0, 0]) ** 2)
    np.testing.assert_allclose(p, 0.5, atol=1e-13)
    # One quarter-turn of a relative degree phase, then the native half contact:
    # the target-port amplitude measures the antisymmetric spatial component.
    effects = half @ np.diag([1, 1j])
    np.testing.assert_allclose(
        np.outer(effects[0].conj(), effects[0]), block[1:, 1:], atol=2e-14
    )
    return 0.5


def formation(board, seed, level, p):
    rng = np.random.default_rng([seed, level])
    edges = np.array(board.seams)
    coins = rng.binomial(1, p, (len(edges), 2)).astype(np.uint8)
    repaired = coins.all(axis=1)
    live = edges[~repaired]
    writes = np.bincount(edges[repaired].ravel(), minlength=board.n)
    a, b = live.T
    g = coo_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(board.n, board.n)
    ).tocsr()
    _, labels = connected_components(g, directed=False)
    ids = np.flatnonzero(labels == np.argmax(np.bincount(labels)))
    sub = g[ids][:, ids]
    lap = diags(np.asarray(sub.sum(axis=1)).ravel()) - sub
    return (
        ids,
        writes,
        lap.toarray(),
        sub,
        {
            "coin_pairs_packed_little_bit_first_hex": np.packbits(
                coins.ravel(), bitorder="little"
            )
            .tobytes()
            .hex(),
            "seams": len(edges),
            "written": int(sum(repaired)),
            "largest_component_ports": len(ids),
            "largest_component_cells": len(np.unique(ids // 12)),
            "largest_component_root": int(ids[0]),
            "full_board_band_multiplicities": np.bincount(writes, minlength=7).tolist(),
        },
    )


def synthesize(target, graph, degree, root):
    order, parent = breadth_first_order(graph, root, directed=False)
    # At each splitter both outgoing intensities are available from subtree
    # sums. This also compiles exact zero amplitudes on a restricted domain.
    mass = abs(target) ** 2
    retained = np.zeros(len(target))
    for child in order[:0:-1]:
        a = parent[child]
        retained[child] = mass[a]
        mass[a] += mass[child]
    amp = np.zeros_like(target, dtype=complex)
    amp[root] = 1
    for child in order[1:]:
        a = parent[child]
        theta = 2 * np.arctan2(np.sqrt(mass[child]), np.sqrt(retained[child]))
        gain = np.exp(-1j * (K + 0.5) * theta) * amp[a]
        amp[a] = gain * np.cos(theta / 2)
        amp[child] = 1j * gain * np.sin(theta / 2)
    # A full contact is SWAP on binary blank inputs. A tree route to the
    # degree-five root, phase there and inverse route implement any site phase.
    assert degree[root] == 5
    amp *= np.exp(-1j * ((np.angle(amp) - np.angle(target)) % (2 * np.pi)))
    residual = float(np.linalg.norm(amp - target))
    if residual > 1e-11:
        raise ArithmeticError(f"Tree source misses spatial preparation by {residual}")
    return amp, residual


def clock(energies):
    # A Q-arm clock delays the source release, then recombines timing arms at
    # a Fourier port. Common-reference storage has zero energy on a degree-six
    # isolated port. Bins retain the same load, and measurement labels are kept.
    x = (energies[:, None] - BINS) / PERIOD
    kernel = np.exp(-1j * np.pi * (Q - 1) * x) * np.sinc(Q * x) / np.sinc(x)
    np.testing.assert_allclose((abs(kernel) ** 2).sum(axis=1), 1, rtol=0, atol=2e-12)
    return kernel


def band_window(r):
    low = K * (6 - r)
    window = (BINS >= low - 1) & (BINS <= low + 7)
    # Foreign bands are at least one unit from each kept bin. Bound their
    # Fourier sidelobes using the nearest possible frequency on the circle.
    others = [(K * (6 - j), K * (6 - j) + 6) for j in range(7) if j != r]
    bins = BINS[window]
    distances = np.full(len(bins), np.inf)
    for lo, hi in others:
        for shift in (-PERIOD, 0, PERIOD):
            d = np.maximum(np.maximum(lo + shift - bins, bins - hi - shift), 0.0)
            distances = np.minimum(distances, d)
    leakage = float(np.sum(1 / (Q * np.sin(np.pi * distances / PERIOD)) ** 2))
    return window, leakage


def counts(rng, n, p):
    p = np.asarray(p, float)
    if p.min() < -1e-10 or abs(p.sum() - 1) > 1e-8:
        raise ValueError("Unnormalised pointer distribution")
    p = np.maximum(p, 0).round(12)
    return rng.multinomial(n, p / p.sum())


def bond_meter(graph, wave):
    a, b = graph.nonzero()
    keep = a < b
    # Reverse the order of a pair if its addressed degree phase is applied
    # on the other end. The antisymmetric projection is orientation independent.
    difference = wave[a[keep]] - wave[b[keep]]
    return float(np.vdot(difference, difference).real / (2 * keep.sum()))


def line(rng, kernel, coeff, window, leak, lam, edge_count, odd):
    weights = abs(coeff) ** 2
    hist = counts(rng, SHOTS, weights @ abs(kernel) ** 2)
    idx = np.flatnonzero(window)
    selected = hist[idx]
    two = int(max(selected[:-1] + selected[1:]))
    # Frequency coverage includes all frequencies in the exact band by a guard
    # of one unit. A pure line places >0.80 probability in its nearest two bins.
    # Bonferroni over both banks, three boards, six modes and at most 256 pairs.
    eps = COUNT_EPS
    margin = 0.80 * (selected.sum() / SHOTS - eps - leak) - (two / SHOTS + eps)
    hits = int(rng.binomial(SPACE_SHOTS, odd))
    return {
        "spatial_lambda": float(lam),
        "spatial_odd_counts": hits,
        "lambda_from_bond_reads": float(2 * edge_count * hits / SPACE_SHOTS),
        "clock_counts": hist.tolist(),
        "selected_fraction": float(selected.sum() / SHOTS),
        "strongest_two_adjacent_share_of_tag": float(two / selected.sum())
        if selected.sum()
        else None,
        "single_line_rejection_margin": float(margin),
        "excludes_one_line_in_band_at_simultaneous_99": bool(margin > 0),
    }


def screen(level, seed, p):
    board = federation(level)
    ids, writes, lap, graph, facility = formation(board, seed, level, p)
    r = int(np.argmax(np.bincount(writes[ids], minlength=7)))
    src = int(np.flatnonzero(writes[ids] == r)[0])
    root = int(np.flatnonzero(board.degree[ids] == 5)[0])
    lo = K * (6 - r)
    h = lap / 2 + np.diag(K * (6 - writes[ids]))
    e, v = eigh(h)
    for j in range(7):
        bottom = K * (6 - j)
        m = (e >= bottom - 1e-9) & (e <= bottom + 6 + 1e-9)
        assert m.sum() == sum(writes[ids] == j)
    w, leak = band_window(r)
    ef = clock(e)
    rng = np.random.default_rng([seed, level, 1])
    which = (e >= lo - 1e-9) & (e <= lo + 6 + 1e-9)
    core_weight = float((v[src, which] ** 2).sum())
    lambdas, spatial = eigh(lap, subset_by_index=(0, 4))
    edge_count = graph.nnz // 2
    lines = []
    controls = []
    control_residuals = []
    for lam, phi in zip(lambdas, spatial.T, strict=True):
        prepared, residual = synthesize(phi, graph, board.degree[ids], root)
        coef = v.T @ prepared
        mean = abs(coef[which]) ** 2 @ e[which] / sum(abs(coef[which]) ** 2)
        std = float(
            np.sqrt(
                abs(coef[which]) ** 2
                @ ((e[which] - mean) ** 2)
                / sum(abs(coef[which]) ** 2)
            )
        )
        row = line(rng, ef, coef, w, leak, lam, edge_count, bond_meter(graph, prepared))
        row.update(
            source_residual=residual,
            within_exact_band_mean=float(mean),
            within_exact_band_width=std,
        )
        lines.append(row)
        # Removing contrast leaves exactly L/2 on this same occupied component.
        free = lap @ prepared / 2
        energy = float(np.vdot(prepared, free).real)
        control_residual = float(np.linalg.norm(free - energy * prepared))
        assert control_residual < 1e-10
        control_residuals.append(control_residual)
        controls.append(
            counts(rng, SHOTS, abs(clock(np.array([energy]))[0]) ** 2).tolist()
        )
    ext = np.array(board.external) // 12
    i, j = ext.T
    cg = coo_matrix(
        (np.ones(2 * len(i)), (np.r_[i, j], np.r_[j, i])),
        shape=(board.cells, board.cells),
    ).tocsr()
    distance = shortest_path(cg, indices=int(ids[src] // 12), unweighted=True)[
        ids // 12
    ].astype(int)
    # Coarse graining the Fourier ports means summing their position probabilities;
    # coherent addition would erase the retained energy record and give another read.
    filtered = (ef[:, w].T * v[src]).T
    passes = float(np.linalg.norm(filtered) ** 2)
    flights = []
    for t in TIMES:
        wave = v @ (np.exp(-1j * t * e)[:, None] * filtered)
        positions = np.sum(abs(wave) ** 2, axis=1)
        reception = np.bincount(
            distance, weights=positions, minlength=int(distance.max()) + 1
        )
        c = counts(rng, SHOTS, np.r_[reception, 1 - passes])
        flights.append(
            {
                "exposure": t,
                "radial_counts_then_tag_failure": c.tolist(),
                "two_or_more_cells_all_launches": float(c[2:-1].sum() / SHOTS),
                "two_or_more_cells_given_tag": float(c[2:-1].sum() / c[:-1].sum()),
            }
        )
    # Fresh assays, with separate seeds so no first-pass count is reused.
    follow = np.random.default_rng([seed, level, 2])
    spectral_rechecks = []

    def reread(coef):
        joint = abs((abs(coef) ** 2) @ ef) ** 2
        total = (abs(coef) ** 2) @ abs(ef) ** 2
        results = counts(follow, SHOTS, np.r_[joint, total - joint]).reshape(2, Q)
        accepted = results[:, w].sum()
        recovered = results[0, w].sum()
        # If the chosen band holds a single frequency, its projected timing
        # vector must survive a reversal of this same spatial preparation.
        p_lower = max(0, accepted / SHOTS - COUNT_EPS)
        signal = p_lower - leak - (1 - p_lower + leak) * np.sqrt(leak)
        q_lower = float(max(0, signal) ** 2)
        margin = float(q_lower - (recovered / SHOTS + COUNT_EPS))
        return {
            "returned_then_elsewhere_by_clock_bin": results.tolist(),
            "mode_recovery_given_tag": float(recovered / accepted),
            "joint_recovery_lower_if_single_line": q_lower,
            "joint_single_line_rejection_margin": margin,
            "joint_excludes_one_line_in_band_at_simultaneous_99": margin > 0,
        }

    for phi in spatial.T:
        wave, _ = synthesize(phi, graph, board.degree[ids], root)
        spectral_rechecks.append(reread(v.T @ wave))
    resonant = np.flatnonzero(writes[ids] == r)
    domain_graph = graph[resonant][:, resonant]
    _, labels = connected_components(domain_graph, directed=False)
    members = resonant[labels == np.argmax(np.bincount(labels))]
    phi = np.zeros(len(ids))
    phi[members] = 1 / np.sqrt(len(members))
    prepared, _ = synthesize(phi, graph, board.degree[ids], root)
    coef = v.T @ prepared
    domain = graph[members][:, members]
    rest = line(
        follow,
        ef,
        coef,
        w,
        leak,
        0,
        domain.nnz // 2,
        bond_meter(domain, prepared[members]),
    )
    rest.update(
        site_count=len(members),
        cells=len(np.unique(ids[members] // 12)),
        root=int(ids[members[0]]),
        spatial_reread=reread(coef),
    )
    return {
        "cells": board.cells,
        "seed": seed,
        "formation": facility,
        "selected_writes": r,
        "sharp_source": int(ids[src]),
        "transmitting_edges_on_component": edge_count,
        "spatial_rereads": spectral_rechecks,
        "one_write_domain_rest": rest,
        "estimated_tag_foreign_leakage_upper": leak,
        "sharp_source_weight_in_exact_band": core_weight,
        "sharp_source_tag_probability": passes,
        "study": lines,
        "exchange_only_same_modes_clock_counts": controls,
        "exchange_only_same_modes_eigen_residuals": control_residuals,
        "flight": flights,
    }


def run():
    p = native()
    # sinc is decreasing on [0,1]. For a line up to half a bin from its closest
    # center, this lower Riemann bound covers both closest timing outputs.
    ls = np.linspace(0, 0.5, 1001)
    assert min(np.sinc(ls[1:]) ** 2 + np.sinc(1 - ls[:-1]) ** 2) > 0.80
    return {
        "shots": SHOTS,
        "spatial_shots": SPACE_SHOTS,
        "clock_bins": Q,
        "frequency_period": PERIOD,
        "half_contact_positive_probability": p,
        "studies": [screen(level, seed, p) for seed in BANKS for level in range(3)],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(record, stream, allow_nan=False, separators=(",", ":"))
            stream.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            json.loads(json.dumps(record, allow_nan=False)),
            json.loads((HERE / "results.json").read_text()),
            "formed medium",
        )
        print(
            "Reproduced the local repair words, spatial and clock counts, and tagged flights."
        )
