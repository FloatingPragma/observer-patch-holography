"""Compare strict-vertex and owned-cell records using a retained source reference."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from math import log, sqrt
from pathlib import Path
from runpy import run_path

import numpy as np

from .....field_charge_records.paths import joint_support, local_terms
from ...word_reader import word_pairs
from ..direct_read import TIMES
from ..measure import (
    ORIGIN,
    PAIR,
    PAULI,
    STEP,
    TARGET,
    TICKS,
    TRIALS,
    WAITS,
    FanPath,
    archived_fixture,
    capture,
    density,
    evolution,
    interval,
    probability,
    source_columns,
)


HERE = Path(__file__).resolve().parent
CELL = 84
FAN = 2
HUB = 1011
EARLY_LAG = 8
LATE = 497
VPORT = np.array([10, -1 + 8j, -1 - 3j, -10j])
VHUB = np.array(
    [13, -4 - 11j, 24 + 16j, -13 + 13j, -18 - 30j, 3 - 12j, 10 + 20j, 28 - 1j]
)
FAMILY = 6
READ_ERROR = 1e-10


def physical_boundary(path, source):
    indices, addresses, heads, edges, _ = joint_support(
        path, source, path.sign, local_terms(path)
    )
    np.testing.assert_array_equal(indices, np.arange(len(path.basis)))
    assert all(x // 12 == CELL for x in path.cells[FAN])
    assert ORIGIN // 12 != CELL and ORIGIN not in path.cells[FAN]
    vertices = set(range(12 * CELL, 12 * (CELL + 1)))
    board_edges = {
        tuple(sorted((x, y))) for x in vertices for y in path.board.neighbours[x]
    }
    assert sum(set(e) <= vertices for e in board_edges) == 30
    assert sum(bool(set(e) - vertices) for e in board_edges) == 3
    incident_1010 = sorted(e for e in board_edges if TARGET in e)
    assert len(incident_1010) == 6
    local_indices = [i for i, e in enumerate(edges) if set(e) <= vertices]
    boundary_indices = [i for i, e in enumerate(edges) if len(set(e) & vertices) == 1]
    point_indices = [i for i, e in enumerate(edges) if TARGET in e]
    assert (
        len(local_indices) == 9
        and len(boundary_indices) == 2
        and len(point_indices) == 3
    )
    assert all(
        tuple(sorted(edges[i])) in board_edges
        for i in (*local_indices, *boundary_indices)
    )
    selected = addresses == TARGET
    owned = addresses // 12 == CELL

    def distinct(a):
        return np.unique(a, axis=0)

    point_rows = distinct(heads[selected][:, point_indices])
    full_cell_rows = distinct(
        np.column_stack((addresses[owned], heads[owned][:, local_indices]))
    )
    point_owned = distinct(heads[selected][:, local_indices])
    boundary = distinct(heads[owned][:, boundary_indices])
    assert point_rows.shape == (1, 3)
    assert point_owned.shape == (2, 9) and boundary.shape == (1, 2)
    assert full_cell_rows.shape == (20, 10)
    return {
        "any_sharp_1010_only_vertex_instrument_has_one_local_word": True,
        "sharp_1010_complete_cell_word_count": len(point_owned),
        "address_unresolved_complete_cell_word_count": len(full_cell_rows),
        "sharp_1010_vertex_active_links_and_heads": [
            [list(edges[i]), int(end)]
            for i, end in zip(point_indices, point_rows[0], strict=True)
        ],
        "all_vertex_incident_links_including_frozen": incident_1010,
        "active_cell_exterior_links_and_constant_heads": [
            [list(edges[i]), int(end)]
            for i, end in zip(boundary_indices, boundary[0], strict=True)
        ],
        "all_cell_incident_external_links_including_frozen": sorted(
            e for e in board_edges if set(e) - vertices
        ),
        "cell_owned_internal_head_registers": 30,
        "reconstructed_support_rows": len(indices),
    }


def local_basis(path, pol):
    # All cells except the material cell have the same words at each of
    # these ports. Sorting against the excess chart aligns both signs.
    shared = None
    blocks = []
    for i, port in enumerate(path.cells[FAN]):
        block = path.blocks[FAN, i]
        spectators = block.words[:FAN] + block.words[FAN + 1 :]
        if shared is None:
            shared = spectators
        assert shared == spectators
        words = block.words[FAN]
        target = [w ^ pol for w in sorted(w ^ pol for w in words)]
        order = [words.index(w) for w in target]
        block_faces = []
        for face, a, b, _ in word_pairs(path, words):
            block_faces.append((face, tuple(sorted((order.index(a), order.index(b))))))
        if port == TARGET:
            assert block_faces == [(0, (0, 1))] and len(words) == 2
        if port == HUB:
            assert len(words) == 4
            left, right = (path.field.masks[j] for j in (3, 0))
            assert left & right == 0
            for j, word in enumerate(target):
                assert word ^ left == target[j ^ 2]
                assert word ^ right == target[j ^ 1]
                for f in (0, 3):
                    mask = path.field.masks[f]
                    clockwise = path.field.clockwise[f]
                    assert word & mask in (clockwise, clockwise ^ mask)
        blocks.append((port, block, target, order, block_faces))
    return blocks


def partial_transpose(rho):
    d = len(rho) // 2
    return rho.reshape(d, 2, d, 2).transpose(0, 3, 2, 1).reshape(rho.shape)


def min_checked(matrix):
    """Enclose the numerical eigensolve within the separate read allowance."""
    values, frame = np.linalg.eigh(matrix)
    dim = len(matrix)
    assert dim <= 40 and np.max(abs(values)) <= 1
    residual = np.linalg.norm((frame * values) @ frame.conj().T - matrix, ord="fro")
    defect = np.linalg.norm(frame.conj().T @ frame - np.eye(dim), ord="fro")
    assert residual < 1e-12 and defect < 1e-12
    return float(values[0])


def short_record(rho, pad):
    trace = float(rho.trace().real)
    return {
        "probability_per_original": trace,
        "partial_transpose_min_interval_per_original": interval(
            min_checked(partial_transpose(rho)), pad
        ),
    }


def matrices_at_ports(wave, blocks):
    columns = []
    for port, block, words, order, pairs in blocks:
        raw = np.moveaxis(wave[block.rows].reshape(*block.shape, 2), FAN, -2)
        tensor = raw.reshape(-1, len(words), 2)[:, order, :]
        columns.append((port, words, pairs, tensor))
    return columns


def receiver(wave, blocks, stock, initial_error):
    pad = READ_ERROR + (2 * sqrt(stock) + initial_error) * initial_error
    matrices = matrices_at_ports(wave, blocks)
    all_rho = density(np.concatenate([mt[3] for mt in matrices], axis=1))
    records = []
    for port, words, pairs, tensor in matrices:
        rho = density(tensor)
        row = short_record(rho, pad)
        row["address"] = port
        row["ordered_reachable_words"] = words
        row["one_face_reads"] = [
            {
                "face": face,
                "words": [words[i], words[j]],
                **short_record(density(tensor[:, [i, j], :]), pad),
            }
            for face, (i, j) in pairs
        ]
        records.append(row)
    overall = short_record(all_rho, pad)
    overall["refusals_first_gate_then_outside_cell"] = [
        1 - stock,
        max(0, stock - overall["probability_per_original"]),
    ]
    return {
        "numerical_probability_and_matrix_error_per_original": pad,
        "whole_cell_without_sharp_address": overall,
        "per_address": records,
    }, pad


def operator_basis(n, choices):
    matrix = np.array([[1]], dtype=complex)
    for axis in choices:
        matrix = np.kron(matrix, PAULI[axis])
    assert len(matrix) == 2**n
    return matrix


def qubit_chsh(rho, pad, p):
    correlations = np.array(
        [[np.trace(rho @ np.kron(a, b)).real / p for b in PAULI[1:]] for a in PAULI[1:]]
    )
    left, singular, right = np.linalg.svd(correlations)
    assert np.linalg.norm((left * singular) @ right - correlations) < 1e-13
    assert np.linalg.norm(left.T @ left - np.eye(3)) < 1e-13
    assert np.linalg.norm(right @ right.T - np.eye(3)) < 1e-13
    # Every conditional entry costs at most 2 pad/(p-pad).
    # The norm of the two largest singular values is 1-Lipschitz in the
    # Frobenius norm (it is the maximum of rank-two projection norms).
    allowance = 12 * pad / (p - pad) + 1e-12
    return {
        "normalized_traceless_Pauli_correlations": correlations.tolist(),
        "optimal_unfiltered_traceless_CHSH": interval(
            float(2 * np.linalg.norm(singular[:2])), allowance
        ),
    }


def fixed_witness(rho, vector, stock, pad, rng):
    n = len(rho).bit_length() - 1
    assert len(rho) == 2**n and len(rho) in (4, 8)
    p = float(rho.trace().real)
    assert p > pad
    witness = (
        partial_transpose(np.outer(vector, vector.conj()))
        / np.vdot(vector, vector).real
    )
    axes = tuple(product(range(4), repeat=n))
    weights = np.array(
        [np.trace(witness @ operator_basis(n, alpha)).real / len(rho) for alpha in axes]
    )
    mean = float(np.trace(witness @ rho).real / p)
    diag = np.eye(2)
    dephased = np.where(np.tile(diag, (2 ** (n - 1), 2 ** (n - 1))), rho, 0)
    erased = float(np.trace(witness @ dephased).real / p)
    assert min_checked(partial_transpose(dephased)) >= -pad
    assert erased >= -pad / p
    # Each moment is charged to one analyzer bank. Identity slots use X
    # as the simultaneous setting and ignore that spin outcome.
    signs = np.array(list(product((1, -1), repeat=n)))
    running = float(weights[0])
    bounds = 0.0
    results = []
    for setting in product((1, 2, 3), repeat=n):
        probabilities = []
        for outcome in signs:
            operator = np.array([[1]], dtype=complex)
            for local, value in zip(setting, outcome, strict=True):
                operator = np.kron(operator, (PAULI[0] + value * PAULI[local]) / 2)
            probabilities.append(max(0.0, float(np.trace(operator @ rho).real)))
        np.testing.assert_allclose(sum(probabilities), p, atol=1e-14)
        histogram = rng.multinomial(
            TRIALS, [1 - stock, max(0, stock - p), *probabilities]
        )
        accepted = int(histogram[2:].sum())
        assert accepted > 0 and int(histogram.sum()) == TRIALS
        scores = np.zeros(len(signs))
        for alpha, weight in zip(axes[1:], weights[1:], strict=True):
            if setting == tuple(a or 1 for a in alpha):
                scores += weight * np.prod(signs[:, [bool(a) for a in alpha]], axis=1)
        running += float(histogram[2:] @ scores / accepted)
        bounds += float(np.ptp(scores) ** 2 / accepted)
        results.append(
            {
                "settings": "".join("IXYZ"[i] for i in setting),
                "counts": histogram.tolist(),
            }
        )
    model_allowance = (float(abs(weights).sum()) + abs(mean)) * pad / (p - pad)
    radius = sqrt(log(2 * FAMILY / 0.01) * bounds / 2) + model_allowance
    subsets = []
    if n == 3:
        for size in (2, 3):
            for sub in combinations(range(4), size):
                indices = [2 * i + j for i, j in product(sub, (0, 1))]
                subsets.append(
                    {
                        "word_indices": sub,
                        **short_record(rho[np.ix_(indices, indices)], pad),
                    }
                )
    return {
        "selected_probability_per_original": p,
        "fixed_score_conditional_enclosure": interval(mean, model_allowance),
        "selected_partial_transpose": short_record(rho, pad),
        "witness_after_reference_dephasing_conditional": interval(
            erased, model_allowance
        ),
        "measured_simultaneous_99_percent_interval": interval(running, radius),
        "axes_order": ["cell_factor_" + str(j) for j in range(n - 1)]
        + ["reference_98"],
        "counter_order": ["first_gate_refusal", "outside_selected_port"]
        + ["".join("+" if s == 1 else "-" for s in a) for a in signs],
        "all_word_subsets_at_hub": subsets,
        "unfiltered_qubit_CHSH": qubit_chsh(rho, pad, p) if n == 2 else None,
        "banks": results,
    }


def case(path, rng):
    source = source_columns(path)
    support = physical_boundary(path, source)
    mask = sum(path.field.bits.values()) if path.sign < 0 else 0
    blocks = local_basis(path, mask)
    direct = evolution(path.h, source, TIMES, sqrt(2), scipy_ticks=(LATE,))
    arms = []
    certified = []
    for first in (None, *TICKS):
        if first is None:
            passages = direct
            stock = 1.0
        else:
            wave, error = direct[first]
            _, gated, _ = capture(path, wave, tuple(w ^ mask for w in PAIR))
            stock = probability(gated)
            passages = evolution(
                path.h,
                gated,
                WAITS,
                sqrt(2) * (sqrt(stock) + error),
                error,
                scipy_ticks=(EARLY_LAG,) if first == TICKS[0] else (),
            )
        arm = {
            "first_gate": float(first * STEP) if first else None,
            "first_gate_weight": stock,
        }
        rows = []
        for s, (wave, error) in passages.items():
            report, pad = receiver(wave, blocks, stock, error)
            rows.append({"wait": float(s * STEP), **report})
            if (first is None and s == LATE) or (first == TICKS[0] and s == EARLY_LAG):
                targets = (TARGET, HUB) if first is None else (TARGET,)
                for port, words, _, tensor in matrices_at_ports(wave, blocks):
                    if port in targets:
                        test = fixed_witness(
                            density(tensor),
                            VPORT if port == TARGET else VHUB,
                            stock,
                            pad,
                            rng,
                        )
                        certified.append(
                            {
                                "first_gate": arm["first_gate"],
                                "wait": float(s * STEP),
                                "charge_address": port,
                                "ordered_cell_words": words,
                                **test,
                            }
                        )
        arm["appointments"] = rows
        arms.append(arm)
    return {
        "sign": path.sign,
        "exact_local_census": support,
        "complete_timings": arms,
        "counted_local_witnesses": certified,
    }


def run():
    with archived_fixture() as archive:
        paths = [FanPath(archive, 4, s, asymmetric=True) for s in (1, -1)]
        a, b = paths
        mask = sum(a.field.bits.values())
        complement = [b.ids[p, tuple(w ^ mask for w in words)] for p, words in a.basis]
        assert (a.h != b.h[complement][:, complement]).nnz == 0
        np.testing.assert_array_equal(a.source(), -b.source()[complement])
        rng = np.random.default_rng(840101098)
        cases = [case(path, rng) for path in paths]
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "site_cell_and_reference": [TARGET, CELL, 98],
        "both_signed_graphs_are_exact_head_complements": True,
        "number_of_original_starts_per_simulated_setting": TRIALS,
        "count_family_cohort_number": FAMILY,
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, allow_nan=False, indent=2) + "\n")
    else:
        compare = run_path(str(HERE.parents[3] / "replay.py"))["compare"]
        compare(result, json.loads((HERE / "records.json").read_text()), "fresh_port")
        print("Replayed strict-vertex, complete cell and six local count witnesses.")


if __name__ == "__main__":
    main()
