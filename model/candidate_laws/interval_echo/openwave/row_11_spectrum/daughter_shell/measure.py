"""Energy-shell test of an autonomous charged triplet and its separated products."""

import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_02_dressed import (
    measure as triplet,
)
from model.candidate_laws.interval_echo.openwave.row_01_charge.measure import (
    sharp_contrast,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
K = int(triplet.K)
EXPOSURES = (Q(0), Q(1, 16), Q(1, 8))
RADIUS = 4
ORDER = 64
SHOTS = 2**22
TICKS = 10**12


def neutral_potential(dh, dx, adjacent, strength=K):
    return strength * (12 + 2 * (int(dh) - int(dx) + int(adjacent)))


def pair_high(dh, dx, adjacent):
    return dh == 6 and dx == 5 and adjacent


def daughter_ceilings(board):
    """Integer row bounds on every single-vacancy and neutral-pair energy."""
    common = Counter(pair for nb in board.neighbours for pair in combinations(nb, 2))
    assert all(nb == sorted(nb) for nb in board.neighbours)
    common_limit = max(common.values())
    assert common_limit == 2 and set(board.degree) == {5, 6}
    weights = np.where(board.degree == 6, 8, 1)
    vacancy = max(
        K * (2 * int(d) - 6)
        + Q(int(d), 2)
        + Q(int(sum(weights[board.neighbours[i]])), 2 * int(weights[i]))
        for i, d in enumerate(board.degree)
    )
    # Only a shared neighbour can make a new adjacent pair. A target
    # regular/gateway requirement allows one offspring per such neighbour.
    # A direct interchange of the existing two ends is a further possibility.
    caps = []
    for dh, dx, adj in product((5, 6), (5, 6), (0, 1)):
        degree = dh + dx - adj
        is_high = pair_high(dh, dx, adj)
        higher_neighbours = min(degree, common_limit + int(adj and not is_high))
        weight = 7 if is_high else 1
        weighted_row_sum = Q(
            7 * higher_neighbours + degree - higher_neighbours, 2 * weight
        )
        caps.append(neutral_potential(dh, dx, adj) + Q(degree, 2) + weighted_row_sum)
    return vacancy, max(caps), common_limit


def pair_moves(board, word):
    hole, excess = word
    for p in board.neighbours[hole]:
        yield (p, hole if p == excess else excess)
    for p in board.neighbours[excess]:
        if p != hole:
            yield (hole, p)


def separated_cell_check(board, source, dest):
    """Compare the actual compression with the independent daughter sum."""
    c, n = dest // 12, source[2] // 12
    charge = list(range(12 * c, 12 * (c + 1)))
    neutral = [
        p
        for p in range(12 * n, 12 * (n + 1))
        if not set(board.neighbours[p]).intersection(charge)
    ]
    assert dest in charge and {source[1], source[2]}.issubset(neutral)
    pairs = list(product(neutral, repeat=2))
    pairs = [z for z in pairs if z[0] != z[1]]
    pairs_index = {z: i for i, z in enumerate(pairs)}
    charge_index = {v: i for i, v in enumerate(charge)}
    h1 = np.zeros((len(charge), len(charge)), dtype=int)
    h2 = np.zeros((len(pairs), len(pairs)), dtype=int)
    for i, u in enumerate(charge):
        d = board.degree[u]
        h1[i, i] = 2 * K * (2 * d - 6) + d
        for v in board.neighbours[u]:
            if v in charge_index:
                h1[i, charge_index[v]] -= 1
    for i, (u, x) in enumerate(pairs):
        neighbours = tuple(pair_moves(board, (u, x)))
        h2[i, i] = 2 * neutral_potential(
            board.degree[u], board.degree[x], x in board.neighbours[u]
        ) + len(neighbours)
        for z in neighbours:
            if z in pairs_index:
                h2[i, pairs_index[z]] -= 1
    basis = [(*sorted((v, u)), x) for v in charge for u, x in pairs]
    actual = triplet.restriction(board, basis, K) * 2
    predicted = sparse.kronsum(sparse.csr_matrix(h2), h1, format="csr")
    assert (actual != predicted).nnz == 0
    return {
        "charge_owner": c,
        "neutral_owner": n,
        "charge_ports": charge,
        "neutral_ports": neutral,
        "joint_calibration_dimension": len(basis),
        "maximum_twice_hamiltonian_entry_error": 0,
        "charged_sharp_mean": str(
            Q(int(h1[charge_index[dest], charge_index[dest]]), 2)
        ),
        "neutral_sharp_mean": str(
            Q(int(h2[pairs_index[source[1:]], pairs_index[source[1:]]]), 2)
        ),
    }


def free_clicks(board, source, target, strength):
    """Exact rational Taylor coefficient on a finite ball and unrestricted-flight errors."""
    states = triplet.ball(board, (source,), RADIUS)
    h = 2 * triplet.restriction(board, states, strength)
    i, j = states.index(source), states.index(target)
    twice_diag = h.diagonal().astype(int)
    center = int((min(twice_diag) + max(twice_diag)) // 2)
    diagonal = (twice_diag - center).tolist()
    h.setdiag(0)
    h.eliminate_zeros()
    assert np.all(h.data == -1)
    neighbours = [
        h.indices[h.indptr[k] : h.indptr[k + 1]].tolist() for k in range(len(states))
    ]
    norm = max(abs(d) + len(ns) for d, ns in zip(diagonal, neighbours, strict=True))
    vector = [0] * len(states)
    vector[i] = 1
    coefficients = []
    for power in range(ORDER + 1):
        coefficients.append(vector[j])
        if power < ORDER:
            vector = [
                d * v - sum(vector[n] for n in ns)
                for v, d, ns in zip(vector, diagonal, neighbours, strict=True)
            ]
    records = []
    for t in EXPOSURES:
        terms = [Q(c, factorial(m)) * (-t / 2) ** m for m, c in enumerate(coefficients)]
        real = sum(terms[::4]) - sum(terms[2::4])
        imag = sum(terms[1::4]) - sum(terms[3::4])
        central_probability = real**2 + imag**2
        taylor = (norm * t / 2) ** (ORDER + 1) / factorial(ORDER + 1)
        # Any omitted path must exit at R+1 then reach the one-hop target.
        first = 2 * RADIUS + 1
        hop = 9 * t
        dyson = 2 * hop**first / factorial(first) / (1 - hop / (first + 1))
        error = (taylor + dyson) * (2 * (abs(real) + abs(imag)) + taylor + dyson)
        lower, upper = (
            max(0, central_probability - error),
            min(1, central_probability + error),
        )
        lo_tick, hi_tick = int(lower * TICKS), -(-upper * TICKS // 1)
        records.append(
            {
                "exposure": str(t),
                "joint_probability_ticks": [lo_tick, int(hi_tick)],
                "central_probability_ticks": round(central_probability * TICKS),
            }
        )
    return {
        "ball_dimension": len(states),
        "maximum_half_units_norm": norm,
        "records": records,
    }


def energy_read(board, source, target, dest, strength, rng):
    """Read the relative diagonal energy with a material beam splitter on the crossing."""
    q = next(
        p
        for p in board.neighbours[dest]
        if board.degree[p] == 5 and p not in source and p != dest
    )
    assert board.degree[source[0]] == board.degree[dest] == 6
    labels = triplet.LAW.SRC
    ix = [labels.index(w) for w in ((0, 1), (1, 0))]
    seam = triplet.material_reference.anchor.free_seam(K)
    splitter = expm(-0.5j * np.pi * seam)[np.ix_(ix, ix)]
    # Route the second arm across dest--q, expose q's native degree-five
    # private contact, and undo the swap. Every other sharp load returns.
    swap = expm(-1j * np.pi * seam)
    permutation = [labels.index((b, a)) for a, b in labels]
    np.testing.assert_allclose(swap, np.eye(9)[permutation], atol=1e-13)
    inputs = []
    for sign, angle in ((1, -np.pi / 2), (-1, np.pi / 2)):
        phase = np.exp(1j * K * ((angle % (2 * np.pi)) / K))
        prepared = np.diag([1, phase]) @ splitter[:, 0]
        expected = np.array([1, sign]) / np.sqrt(2)
        np.testing.assert_allclose(abs(np.vdot(expected, prepared)), 1, atol=1e-13)
        inputs.append(prepared)
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    # For the native pi/2 beam splitter, output d minus p reads Y.
    np.testing.assert_allclose(splitter.conj().T @ (-z) @ splitter, y, atol=1e-13)
    states = triplet.ball(board, (source, target), 3)
    ix = [states.index(source), states.index(target)]
    h = triplet.restriction(board, states, strength)
    energy_difference = h[ix[0], ix[0]] - h[ix[1], ix[1]]
    h.setdiag(h.diagonal() - h[ix[0], ix[0]])
    seed = np.zeros((len(states), 2), dtype=complex)
    seed[ix, :] = np.asarray(inputs).T
    records = []
    for t in (Q(1, 64), Q(1, 128)):
        timed = expm_multiply(
            -1j * float(t) * h, seed, traceA=-1j * float(t) * h.diagonal().sum()
        )
        rates = abs(splitter @ timed[ix, :]) ** 2
        shots = []
        exact = []
        for output in rates.T:
            ticks = np.rint(output * TICKS).astype(np.int64)
            ticks = np.r_[ticks, TICKS - sum(ticks)]
            shots.append(rng.multinomial(SHOTS, ticks / TICKS).tolist())
            exact.append((ticks[1] - ticks[0]) / TICKS)
        measured = [(row[1] - row[0]) / SHOTS for row in shots]
        hop = 9 * t
        error = 2 * hop**8 / factorial(8) / (1 - hop / 9)
        records.append(
            {
                "exposure": str(t),
                "plus_minus_counts_p_d_other": shots,
                "slope_from_counts": round(
                    (measured[0] - measured[1]) / (2 * float(t)), 8
                ),
                "expected_finite_time_slope": round(
                    (exact[0] - exact[1]) / (2 * float(t)), 8
                ),
                "whole_board_amplitude_entry_error_upper": float(error),
            }
        )
    return {
        "phase_route": [dest, q],
        "exact_initial_slope": float(energy_difference),
        "source_ball_dimension": len(states),
        "records": records,
    }


def detector_controls(board, source, target, dest):
    receiver = {p for p in range(board.n) if p // 12 == dest // 12}
    cut = [(a, b) for a, b in board.seams if (a in receiver) != (b in receiver)]
    assert len(cut) == 3 and tuple(sorted((source[0], dest))) in cut

    def signed_charge(word):
        return sum(
            sharp_contrast(2 if p == word[2] else 0 if p in word[:2] else 1, 0) - 1
            for p in receiver
        )

    initial = signed_charge(source)
    wanted = signed_charge(target)
    assert (initial, wanted) == (0, -1)
    positive = all(sharp_contrast(load, 0) == load for load in (0, 1, 2))
    assert positive
    return {
        "sealed_owner_edges": cut,
        "sealed_initial_charge": initial,
        "required_receiver_charge": wanted,
        "sealed_click_and_other": [0, SHOTS],
        "target_prepared_at_zero_click_and_other": [SHOTS, 0],
    }


def run():
    native = json.loads((triplet.HERE / "results.json").read_text())
    rows = []
    for level, stored in enumerate(native["upper_band"]):
        board = federation(level)
        band = triplet.high_band(board)
        for field in (
            "board_cells",
            "highest_band_dimension",
            "low_block_upper_bound",
            "high_block_enclosure",
            "sharp_energy",
        ):
            assert band[field] == stored[field]
        vacancy, neutral, shared = daughter_ceilings(board)
        independent = vacancy + neutral
        parent = Q(band["high_block_enclosure"][0])
        # One sharp top address has 15 different-load swaps. The rest of
        # the spectrum is bounded by the complementary principal submatrix.
        s = min(triplet.top_addresses(board))
        swaps = len(triplet.transitions(board, s))
        variance = Q(swaps, 4)
        mean = K * int(triplet.potential(board, s)) + Q(swaps, 2)
        gap = mean - Q(band["low_block_upper_bound"])
        acceptance = 1 - variance / (variance + gap**2)
        assert parent > independent and acceptance > 0
        rows.append(
            {
                "cells": board.cells,
                "max_common_neighbours": shared,
                "certified_parent_band_rank": band["highest_band_dimension"],
                "parent_floor": str(parent),
                "parent_sharp_mean": str(mean),
                "parent_sharp_variance": str(variance),
                "sharp_upper_band_weight_lower": str(acceptance),
                "vacancy_ceiling": str(vacancy),
                "neutral_pair_ceiling": str(neutral),
                "independent_daughter_sum_ceiling": str(independent),
                "minimum_energy_deficit": str(parent - independent),
            }
        )
    board = federation(0)
    source = tuple(native["high_release"]["launch"])
    assert min(triplet.top_addresses(board)) == source
    (dest,) = (q for q in board.neighbours[source[0]] if q // 12 != source[0] // 12)
    target = (*sorted((dest, source[1])), source[2])
    assert target in triplet.transitions(board, source)
    assert source[2] in board.neighbours[source[1]]
    calibration = separated_cell_check(board, source, dest)
    changed = free_clicks(board, source, target, 0)
    native_clicks = free_clicks(board, source, target, K)
    rng = np.random.default_rng(18439)
    for apparatus in (native_clicks, changed):
        for row in apparatus["records"]:
            count = int(rng.binomial(SHOTS, row["central_probability_ticks"] / TICKS))
            row["simulated_click_and_other"] = [count, SHOTS - count]
    native_meter = energy_read(board, source, target, dest, K, rng)
    changed_meter = energy_read(board, source, target, dest, 0, rng)
    return {
        "energy_reference": "each formation referred to its unit-filled host",
        "boards": rows,
        "source": list(source),
        "target": list(target),
        "calibration": calibration,
        "release": {
            "radius": RADIUS,
            "taylor_order": ORDER,
            "probability_tick_denominator": TICKS,
            "shots_per_setting": SHOTS,
            "seed": 18439,
            "native": native_clicks,
            "exchange_only_control": changed,
            "native_energy_read": native_meter,
            "exchange_only_energy_read": changed_meter,
            **detector_controls(board, source, target, dest),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    receipt = HERE / "results.json"
    result = json.loads(json.dumps(run()))
    if parser.parse_args().write:
        receipt.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(receipt.read_text()):
        raise ValueError("Daughter shell receipt differs")
    print(json.dumps(result["boards"], indent=2))


if __name__ == "__main__":
    main()
