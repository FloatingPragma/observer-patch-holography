"""One signed defect, two written collars and the complete free Interval Echo generator."""

import argparse
import json
from fractions import Fraction
from itertools import product
from math import factorial
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.sparse import coo_matrix, diags, eye
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh, expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.spectral import (
    verified_inertia,
)
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator

from .receiver import receiver_control

# Owner 6 is the archived magnetic preparation cell. Choose its lexicographically
# first external seam and each gateway's first internal degree-five neighbor.
SOURCE_OWNER = 6
STRENGTH = int(KAPPA)
WAITS = (0, 1, 2, 4, 8, 16, 32)
ACCEPTANCE = 0.25


def geometry():
    board = federation(0)
    bridge = min(
        edge for edge in board.external if SOURCE_OWNER in {p // 12 for p in edge}
    )
    launch = next(p for p in bridge if p // 12 == SOURCE_OWNER)
    receive = next(p for p in bridge if p != launch)
    marked = tuple(
        min(
            edge
            for edge in board.seams
            if gateway in edge
            and edge[0] // 12 == edge[1] // 12
            and any(board.degree[p] == 5 for p in edge)
        )
        for gateway in (launch, receive)
    )
    return board, launch, receive, marked


def local_certificate():
    """Check the reduced integer columns against the exact 99-state echo."""

    def identity(size):
        return np.eye(size, dtype=np.int64)

    write = law._write_pattern()  # Scaled by sqrt(2).
    reversal, contrast = law._label_ops()
    written = np.kron(
        identity(law.N_LAB) - reversal + 2 * STRENGTH * contrast, identity(law.N_SRC)
    )
    # U acts as I-WW* on the written output rows. The remaining U entries
    # project to blank, and the native W intertwines the two records.
    p = write @ write.T
    assert np.array_equal(written @ p, p @ written)
    reduced = written @ (2 * identity(law.N_SRC * law.N_LAB) - p)  # Equals 4H.
    assert not np.any(write.T @ reduced)
    chart = np.kron(reversal, identity(law.N_SRC))
    assert np.array_equal(chart @ reduced, reduced @ chart)
    for sign, parity in product((-1, 1), repeat=2):
        pairs = ((1 + sign, 1), (1, 1 + sign), (1, 1))
        states = [(coin, positions) for coin in range(2) for positions in pairs]
        basis = np.zeros((90, len(states)), dtype=np.int64)
        for i, (coin, reads) in enumerate(states):
            for d in (-1, 1):
                row = 9 * law.LABELS.index((coin, d)) + law.SRC.index(reads)
                basis[row, i] = parity if d < 0 else 1
        # Of the four positions coincident with the marked seam, only the two
        # with the correct coin for that orientation share a written projector.
        a, b = (0, 4) if sign == -1 else (3, 1)
        twice_local = 2 * identity(6)
        twice_local[a, a] = twice_local[b, b] = 1
        twice_local[a, b] = twice_local[b, a] = -1
        gap = STRENGTH + (1 - parity) // 2
        assert np.array_equal(reduced @ basis, 2 * gap * (basis @ twice_local))
    # A native write of the appropriate coherent difference prepares the
    # same contrast parity with two equiprobable sharp local charge positions.
    for sign, parity in product((-1, 1), repeat=2):
        first, second = [law.SRC.index(t) for t in ((1 + sign, 1), (1, 1 + sign))]
        preparation = write[:, first] + parity * write[:, second]
        assert np.array_equal(chart @ preparation, parity * preparation)
        assert preparation @ preparation == 4
        for position in (first, second):
            selected = preparation[np.arange(90) % 9 == position]
            assert selected @ selected == 2


def native_generator(board, marked, sign, parity):
    """The integer generator 2H in the fixed parity block, with all seams open."""
    blank = defect_generator(board, sign, KAPPA, reference=0)
    occupancy_diagonal = np.rint(2 * blank.diagonal()).astype(np.int64)
    assert np.array_equal(occupancy_diagonal / 2, blank.diagonal())
    baseline = occupancy_diagonal.copy()
    for a, b in marked:
        baseline[[a, b]] -= 2 * STRENGTH + 1
    diagonal = np.tile(baseline, 1 << len(marked))
    edges = tuple(edge for edge in board.seams if edge not in marked)
    assert len(edges) + len(marked) == len(board.seams)
    rows, cols, values = [], [], []
    n = board.n

    def link(a, b, weight):
        rows.extend((a, b))
        cols.extend((b, a))
        values.extend((weight, weight))

    for word in range(1 << len(marked)):
        for a, b in edges:
            link(word * n + a, word * n + b, -1)
    for j, ((a, b), direction) in enumerate(zip(marked, parity, strict=True)):
        gap = STRENGTH + (1 - direction) // 2
        diagonal += 2 * gap
        bit = 1 << j
        for word in range(1 << len(marked)):
            if word & bit:
                continue
            first, second = (a, b) if sign == -1 else (b, a)
            left = word * n + first
            right = (word | bit) * n + second
            link(left, right, -gap)
            diagonal[[left, right]] -= gap
    pairs = coo_matrix((values, (rows, cols)), shape=(len(diagonal),) * 2).tocsr()
    assert connected_components(pairs, directed=False, return_labels=False) == 1
    return pairs + diags(diagonal, dtype=np.int64, format="csr")


def sharp(frame, n, pos, sbit, rbit):
    return frame[(sbit + 2 * rbit) * n + pos]


def source_slots(board, launch, receive, marked, sign):
    receive_coin = int((receive == marked[1][1]) == (sign == -1))
    return [(sbit + 2 * receive_coin) * board.n + launch for sbit in range(2)]


def joint_outcomes(frame, board, receiver_address):
    """Unnormalized Choi record, after tracing the sender's coin and charge port."""
    choi = np.zeros((4, 4), complex)
    for pos in receiver_address:
        for sender_coin in range(2):
            k = np.stack([sharp(frame, board.n, pos, sender_coin, c) for c in range(2)])
            flat = k.T.reshape(-1)
            choi += np.outer(flat, flat.conj()) * (ACCEPTANCE / 2)
    minimum = np.linalg.eigvalsh(
        choi.reshape(2, 2, 2, 2).transpose(2, 1, 0, 3).reshape(4, 4)
    )[0]
    operator = np.trace(choi.reshape(2, 2, 2, 2), axis1=1, axis2=3) * 2
    return {
        "probability_per_original_bank_start": round(float(np.trace(choi).real), 12),
        "accepted_input_operator_eigenvalues": np.linalg.eigvalsh(operator)
        .round(12)
        .tolist(),
        "reference_partial_transpose_minimum": round(float(minimum), 12),
        "preparation_refusal": 1 - ACCEPTANCE,
        "no_qualified_receiving_charge_given_prep_original_start": round(
            ACCEPTANCE - float(np.trace(choi).real), 12
        ),
    }


def flight(board, launch, receive, marked, sign, parity=(1, 1)):
    h = native_generator(board, marked, sign, parity) / 2
    # The receiver writes its collar first and the post-write gateway read
    # selects its native matching coin. A separate logical input is swapped
    # into the same coin coordinate on the sending collar.
    source = np.zeros((h.shape[0], 2), complex)
    source[source_slots(board, launch, receive, marked, sign), range(2)] = 1
    counts = {}
    moving = source.copy()
    diagonal_trace = float(h.diagonal().sum())
    before = 0
    for wait in WAITS:
        if wait > before:
            interval = -1j * (wait - before)
            moving = expm_multiply(
                interval * h, moving, traceA=interval * diagonal_trace
            )
        before = wait
        np.testing.assert_allclose(moving.conj().T @ moving, np.eye(2), atol=3e-10)
        counts[str(wait)] = {
            "receiver_gateway": joint_outcomes(moving, board, (receive,)),
            "other_marked_receiver_endpoint": joint_outcomes(
                moving, board, (next(p for p in marked[1] if p != receive),)
            ),
        }
    return counts


def reversal_control(board, launch, receive, marked, sign):
    """Coherent sender tags in two conserved external parity sectors."""
    sharp_slot = source_slots(board, launch, receive, marked, sign)[
        int((launch == marked[0][1]) == (sign == -1))
    ]
    branches = {}
    for parity in (1, -1):
        h = native_generator(board, marked, sign, (parity, 1)) / 2
        source = np.zeros((h.shape[0], 1), complex)
        source[sharp_slot, 0] = 1
        wait = -4j
        arrival = expm_multiply(wait * h, source, traceA=wait * h.diagonal().sum())
        rows = []
        for sender_coin in range(2):
            rows.append(
                [sharp(arrival, board.n, receive, sender_coin, c)[0] for c in range(2)]
            )
        response = np.asarray(rows).T
        local = ACCEPTANCE * (response @ response.conj().T)
        branches[str(parity)] = round(float(np.trace(local).real), 12)
    return {
        "gateway_arrival_per_original_start_with_each_sharp_parity": branches,
        "offdiagonal_reference_block_after_tracing_source_parity_exact": 0,
        "always_separable_by_conserved_orthogonal_exterior_labels": True,
    }


def integer_flight(h2, slots, wait):
    """Dyadic Taylor enclosure, returning integer real and imaginary amplitudes."""
    scale, slices, order = 1 << 53, 64, 16
    dimension = h2.shape[0]
    assert dimension < 2048  # sqrt(2 * dimension) < 64 per complex column.
    diagonal = h2.diagonal()
    radius = np.asarray(abs(h2).sum(axis=1)).ravel() - abs(diagonal)
    center = (int(min(diagonal - radius)) + int(max(diagonal + radius))) // 2
    shifted = h2 - center * eye(dimension, dtype=np.int64, format="csr")
    row_bound = int(max(np.asarray(abs(shifted).sum(axis=1)).ravel()))
    assert row_bound <= slices
    steps = wait * slices
    # Every term recurrence is contractive by at least 1/2 before rounding.
    # Each complex-vector rounding adds at most 32/scale to the l2 norm.
    per_step_rounding = Fraction(64 * order, scale)
    per_step_tail = 4 * Fraction(1, 2) ** (order + 1) / factorial(order + 1)
    error = steps * (per_step_tail + per_step_rounding)
    assert error < Fraction(1, 10**9)  # Also keeps every step input norm < 2.
    real = np.zeros((dimension, 2), np.int64)
    real[slots, range(2)] = scale
    imaginary = np.zeros_like(real)
    for _ in range(steps):
        term_real, term_imaginary = real.copy(), imaginary.copy()
        for index in range(1, order + 1):
            divisor = 2 * slices * index
            # A product's integer accumulator is bounded by row_bound times
            # the maximum term size. Partial sums obey the same absolute bound.
            assert np.max(abs(term_real)) < 2 * scale
            assert np.max(abs(term_imaginary)) < 2 * scale
            assert 2 * scale * row_bound + divisor < 2**62
            term_real, term_imaginary = (
                (shifted @ term_imaginary + divisor // 2) // divisor,
                (-(shifted @ term_real) + divisor // 2) // divisor,
            )
            real += term_real
            imaginary += term_imaginary
    return real, imaginary, scale, error


def certified_witness(h2, board, slots, receive):
    real, imaginary, scale, frame_error = integer_flight(h2, slots, wait=4)
    approximate = sp.zeros(4, 4)
    for sender_coin in range(2):
        re = [sharp(real, board.n, receive, sender_coin, c) for c in range(2)]
        im = [sharp(imaginary, board.n, receive, sender_coin, c) for c in range(2)]
        vector = sp.Matrix(
            [
                int(a) + sp.I * int(b)
                for a, b in zip(
                    np.ravel(np.transpose(re)), np.ravel(np.transpose(im)), strict=True
                )
            ]
        )
        approximate += vector * vector.H
    denominator = 8 * scale**2
    transpose = sp.zeros(4, 4)
    for i, j in product(range(4), repeat=2):
        transpose[i, j] = approximate[(j // 2) * 2 + i % 2, (i // 2) * 2 + j % 2]
    # Float selects a witness only. Its expectation is verified as a rational.
    _, modes = np.linalg.eigh(np.asarray((transpose / denominator).evalf(), complex))
    target = modes[:, 0] / (
        modes[np.argmax(abs(modes[:, 0])), 0] / max(abs(modes[:, 0]))
    )
    integers = np.rint(2**18 * np.column_stack((target.real, target.imag))).astype(int)
    vector = sp.Matrix([int(a) + sp.I * int(b) for a, b in integers])
    raw_expectation = sp.expand((vector.H * transpose * vector)[0])
    norm_squared = (vector.H * vector)[0].expand()
    assert raw_expectation.is_Integer and norm_squared.is_Integer
    exact_probe = Fraction(int(raw_expectation), denominator * int(norm_squared))
    # Contract the unnormalized Choi state before the local read. The
    # partial transpose of a rank-one witness has operator norm at most 1.
    choi_error = 3 * frame_error / 4
    upper = exact_probe + choi_error
    assert upper < 0
    return {
        "integer_witness_real_imaginary": integers.tolist(),
        "frame_error_bound": str(frame_error),
        "witness_expectation_on_approximate_frame": str(exact_probe),
        "true_witness_expectation_upper": float(upper),
        "witness_expectation_additive_error_bound": float(choi_error),
    }


def isolation_certificate(h2, sign, parity):
    cuts = (
        (Fraction(591, 16), Fraction(2367, 64), Fraction(37))
        if sign == -1
        else (Fraction(53), Fraction(1697, 32), Fraction(213, 4))
    )
    counts = [verified_inertia(h2 / 2, threshold)["count_below"] for threshold in cuts]
    assert counts == [0, int(parity == (1, 1)), int(parity == (1, 1))]
    return {"dyadic_cutoffs": list(map(str, cuts)), "level_counts_below": counts}


def ground_predictions(h2, board, controls, source, receive, positions):
    _, frame = eigsh(h2 / 2, k=1, which="SA", v0=np.ones(h2.shape[0]))
    ground = frame[:, 0]
    retained = {}
    for name, gate in zip(("x", "z"), controls, strict=True):
        # Three full marked phases give -1 on X's outside words. Z has
        # just the single uncompensated marked phase on the same words.
        output = -ground.astype(complex)
        for sender_coin in range(2):
            indices = [
                (sender_coin + 2 * c) * board.n + p for p in positions for c in range(2)
            ]
            output[indices] = gate @ ground[indices]
        retained[name] = round(float(abs(np.vdot(ground, output)) ** 2), 12)
        assert retained[name] < 1
    admission = ACCEPTANCE * float(sum(abs(ground[source]) ** 2)) / 2
    arrival = admission * float(sum(abs(ground[receive :: board.n]) ** 2))
    return {
        "spectral_admission_per_original_start_numeric": round(admission, 12),
        "ground_and_gateway_capture_per_original_start_numeric": round(arrival, 12),
        "energy_retention_given_projected_ground_numeric": retained,
        "reference_entanglement_after_rank_one_filter_exact": 0,
    }


def run():
    assert KAPPA == 8
    local_certificate()
    board, launch, receive, marked = geometry()
    results = {}
    for sign in (-1, 1):
        h0 = defect_generator(board, sign, KAPPA, reference=0)
        blank = eigsh(
            h0, k=4, which="SA", v0=np.ones(board.n), return_eigenvectors=False
        )
        levels = {}
        isolated = {}
        for parity in product((1, -1), repeat=2):
            h2 = native_generator(board, marked, sign, parity)
            isolated[str(parity)] = isolation_certificate(h2, sign, parity)
            h = h2 / 2
            low = eigsh(
                h, k=4, which="SA", v0=np.ones(h.shape[0]), return_eigenvectors=False
            )
            levels[str(parity)] = np.sort(low).round(10).tolist()
        h2 = native_generator(board, marked, sign, (1, 1))
        slots = source_slots(board, launch, receive, marked, sign)
        readout, ux, uz = receiver_control(board, marked[1], receive, sign, STRENGTH)
        flight_records = flight(board, launch, receive, marked, sign)
        readout["x_z_against_z_x_dark_probability_per_original_start_at_wait_four"] = (
            flight_records["4"]["receiver_gateway"][
                "probability_per_original_bank_start"
            ]
        )
        results[str(sign)] = {
            "lowest_four_blank_levels_numeric": np.sort(blank).round(10).tolist(),
            "lowest_four_written_levels_by_reversal_parity_numeric": levels,
            "verified_first_written_line_isolation": isolated,
            "same_even_parity_arrival_banks": flight_records,
            "wait_four_same_parity_entanglement_certificate": certified_witness(
                h2, board, slots, receive
            ),
            "wait_four_source_parity_tag_control": reversal_control(
                board, launch, receive, marked, sign
            ),
            "registered_protected_receiver_pulses": readout,
            "same_start_when_coarsely_selecting_the_ground_line": ground_predictions(
                h2, board, (ux, uz), slots, receive, readout["local_positions"]
            ),
        }
    return {
        "schema": "oph.material_spin_band.v1",
        "board_cells": board.cells,
        "matter_ports": board.n,
        "native_seams": len(board.seams),
        "sending_gateway": launch,
        "receiving_gateway": receive,
        "written_native_seams": [list(edge) for edge in marked],
        "contrast_magnitude": 1,
        "clock_waits": WAITS,
        "signed_banks": results,
    }


def main():
    if not __debug__:
        raise RuntimeError("Run without optimizations so the exact assertions execute.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    path = Path(__file__).with_name("results.json")
    result = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if parser.parse_args().record:
        path.write_text(result)
    elif path.read_text() != result:
        raise RuntimeError("Independent replay differs from the archived receipt.")
    print(
        "Full registered one-defect dynamics, native collar columns and arrival banks verified."
    )


if __name__ == "__main__":
    main()
