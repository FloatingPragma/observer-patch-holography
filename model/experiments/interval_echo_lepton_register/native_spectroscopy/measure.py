"""Complete Fourier-clock and held arrival records on registered Interval Echo."""

import argparse
import json
from math import log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation

from ..phase_feedback import circuit
from .neutral import blank_neutral_certificate, filled_memory_receiver
from .source import (
    blank_generator,
    energy_certificate,
    geometry,
    preparations,
    written_generator,
)

STEP = 1 / 16
COARSE = 64
FINE = 512
ARRIVAL_WAIT = 8
SHOTS = 2**26
SEED = 726111
CASES = 8


def clock_probabilities(moments):
    """Read the finite clock in its Fourier basis, including all side lobes."""
    n = len(moments)
    fejer = (1 - np.arange(n) / n)[:, None] * moments
    bins = 2 * np.fft.ifft(fejer, axis=0).real - 1 / n
    assert bins.min() >= -2e-12
    np.testing.assert_allclose(bins.sum(axis=0), 1, atol=2e-12)
    bins = np.maximum(bins, 0)
    return (bins / bins.sum(axis=0)).T


def draw(rng, probabilities):
    return [rng.multinomial(SHOTS, row).tolist() for row in probabilities]


def filtered_arrivals(matrix, history, chosen, head_count, rng):
    bins = sorted(set(chosen))
    times = np.arange(len(history))
    signs = np.exp(2j * pi * times[:, None] * bins / len(history))
    # Each selected bin is measured on every input, regardless of its source tag.
    filtered = np.einsum("tb,tij->ijb", signs, history) / len(history)
    filtered = filtered.reshape(matrix.shape[0], 3 * len(bins))
    held = expm_multiply((-1j * ARRIVAL_WAIT / 8) * matrix, filtered)
    outcomes = []
    for time, outputs in ((0, filtered), (ARRIVAL_WAIT, held)):
        local = np.sum(abs(outputs[: 12 * head_count]) ** 2, axis=0).reshape(3, -1)
        outside = np.sum(abs(outputs[12 * head_count :]) ** 2, axis=0).reshape(3, -1)
        chosen_outcomes = np.stack((outside, local), axis=2).reshape(3, -1)
        rejected = 1 - chosen_outcomes.sum(axis=1)
        detector = np.column_stack((chosen_outcomes, rejected))
        outcomes.append(
            dict(
                wait_after_selection=time,
                probabilities=detector.tolist(),
                outside_local_other_counts=draw(rng, detector),
            )
        )
    return bins, outcomes


def receiver_comparison(train, test):
    """Read the held counts using only the calibration bank to choose the decoder."""
    bin_to_tag = train[:, :-1].argmax(axis=0)
    table = np.column_stack(
        [test[:, :-1][:, bin_to_tag == tag].sum(axis=1) for tag in range(3)]
        + [test[:, -1]]
    )
    arrived = table[:, :3]
    right = np.diag(arrived)
    missed_tag = arrived.sum(axis=1) - right

    def count_interval(count):
        # Invert the Bernstein bound for the six rare-event probabilities.
        frequency = count / SHOTS
        r = log(12 * CASES / 0.01) / SHOTS
        slack = sqrt(r / 2)
        lo = max(0.0, sqrt(max(0.0, frequency - r / 6)) - slack) ** 2
        hi = min(1.0, (sqrt(frequency + 7 * r / 6) + slack) ** 2)
        return lo, hi

    right_lower = [count_interval(count)[0] for count in right]
    wrong_upper = [count_interval(count)[1] for count in missed_tag]
    lower_vs_any_constant = min(
        (sum(right_lower) - right_lower[tag] - wrong_upper[tag]) / 3 for tag in range(3)
    )
    # The best constant guess on the held arrivals always chooses the most
    # numerous true tag. This baseline gets to choose after the held bank.
    excess = (right.sum() - arrived.sum(axis=1).max()) / (3 * SHOTS)

    training_arrivals = train[:, :-1]
    normalized = training_arrivals / training_arrivals.sum(axis=1)[:, None]
    partition = normalized[1] > normalized[2]
    held_arrivals = test[:, :-1]
    totals = held_arrivals.sum(axis=1)
    opposite = held_arrivals[:, partition].sum(axis=1) / totals
    radius = sum(sqrt(log(4 * CASES / 0.01) / (2 * n)) for n in totals[1:])
    return dict(
        bin_to_receiver_tag=bin_to_tag.tolist(),
        guessed_tags_and_no_arrival=table.tolist(),
        correct_outside_fraction=float(right.sum() / arrived.sum()),
        outside_after_clock_fraction=float(arrived.sum() / (3 * SHOTS)),
        excess_over_constant_guess_all_starts=float(excess),
        excess_lower_99=float(lower_vs_any_constant),
        opposite_tag_partition_given_arrival=np.flatnonzero(partition).tolist(),
        opposite_tag_contrast_given_arrival=float(opposite[1] - opposite[2]),
        opposite_tag_radius_given_arrival_99=float(radius),
    )


def all_bin_receiver(held_history, fine, head_count, rng):
    """Train and test a receiving observer on every clock-bin and arrival outcome."""
    # The fixed clock filter commutes with the extra free dwell before the port read.
    delivered = np.fft.ifft(held_history[:, 12 * head_count :], axis=0)
    out = np.einsum("bij,bij->jb", delivered.conj(), delivered).real
    assert np.max(out - fine) < 2e-11
    probs = np.column_stack((out, 1 - out.sum(axis=1)))
    train = np.array(draw(rng, probs))
    test = np.array(draw(rng, probs))
    return dict(
        outside_probabilities=out.tolist(),
        outside_bin_training_counts=train.tolist(),
        outside_bin_held_counts=test.tolist(),
        **receiver_comparison(train, test),
    )


def measure(matrix, states, head_count, rng):
    step_generator = (-1j * STEP / 8) * matrix
    offset = int(ARRIVAL_WAIT / STEP)
    assert offset * STEP == ARRIVAL_WAIT
    full_history = expm_multiply(
        step_generator, states, start=0, stop=FINE + offset - 1, num=FINE + offset
    )
    history = full_history[:FINE]
    moments = np.einsum("ij,tij->tj", states.conj(), history)
    coarse = clock_probabilities(moments[:COARSE])
    fine = clock_probabilities(moments)
    train = np.array(draw(rng, fine))
    held = np.array(draw(rng, fine))
    # The tag comparison uses a partition fixed by the first bank.
    partition = train[1] > train[2]
    difference = (int(held[1, partition].sum()) - int(held[2, partition].sum())) / SHOTS
    radius = 2 * sqrt(log(4 * CASES / 0.01) / (2 * SHOTS))
    chosen = train.argmax(axis=1).tolist()
    bins, arrivals = filtered_arrivals(matrix, history, chosen, head_count, rng)
    target_weight = fine[:, bins]
    for row in arrivals:
        np.testing.assert_allclose(
            np.array(row["probabilities"])[:, :-1].reshape(3, len(bins), 2).sum(axis=2),
            target_weight,
            atol=2e-11,
        )
    receiver_tags = train[:, bins].argmax(axis=0)
    outside_counts = np.array(arrivals[-1]["outside_local_other_counts"])[:, :-1:2]
    correct = sum(outside_counts[tag, j] for j, tag in enumerate(receiver_tags))
    receiver = all_bin_receiver(full_history[offset:], fine, head_count, rng)
    np.testing.assert_allclose(
        np.array(receiver["outside_probabilities"])[:, bins],
        np.array(arrivals[-1]["probabilities"])[:, :-1:2],
        atol=2e-11,
    )
    return dict(
        coarse_counts=draw(rng, coarse),
        fine_training_counts=train.tolist(),
        fine_held_counts=held.tolist(),
        fine_probabilities=fine.tolist(),
        selected_bins=chosen,
        selected_bin_centres=(np.array(chosen) * (2 * pi / (FINE * STEP))).tolist(),
        held_read_bins=bins,
        calibrated_receiver_tag_per_bin=receiver_tags.tolist(),
        held_correct_tag_and_arrival=int(correct),
        held_selected_arrivals=int(outside_counts.sum()),
        opposite_tag_partition=np.flatnonzero(partition).tolist(),
        held_opposite_tag_contrast=difference,
        held_opposite_tag_radius_99=radius,
        arrivals=arrivals,
        all_bins_at_held_receiver=receiver,
    )


def run():
    face, tags, band_tags = geometry()
    boards = [federation(level) for level in range(3)]
    board = boards[0]
    blank = blank_generator(board, -1, private=True)
    np.testing.assert_array_equal(
        blank.toarray(), 8 * circuit.free_generator().toarray()
    )
    written = written_generator(board, face, -1, private=True)
    silent = written_generator(board, face, -1, private=False)
    conjugate = written_generator(board, face, 1, private=False)
    assert (silent - conjugate).nnz == 0

    rng = np.random.default_rng(SEED)
    cases = {
        "blank": measure(blank, preparations(board, tags), 1, rng),
        "written_real": measure(written, preparations(board, tags, 0), 8, rng),
        "written_positive_phase": measure(
            written, preparations(board, tags, pi / 2), 8, rng
        ),
        "written_negative_phase": measure(
            written, preparations(board, tags, -pi / 2), 8, rng
        ),
        "written_silent_reference": measure(
            silent, preparations(board, tags, pi / 2), 8, rng
        ),
        "face_band_real": measure(silent, preparations(board, band_tags, 0), 8, rng),
        "face_band_positive_phase": measure(
            silent, preparations(board, band_tags, pi / 2), 8, rng
        ),
        "face_band_negative_phase": measure(
            silent, preparations(board, band_tags, -pi / 2), 8, rng
        ),
    }
    assert len(cases) == CASES
    for name in ("blank", "written_real", "face_band_real"):
        a, b = cases[name]["fine_probabilities"][1:]
        np.testing.assert_allclose(a, b, atol=2e-12, rtol=0)
    for stem in ("written", "face_band"):
        np.testing.assert_allclose(
            cases[stem + "_positive_phase"]["fine_probabilities"],
            np.array(cases[stem + "_negative_phase"]["fine_probabilities"])[[0, 2, 1]],
            atol=2e-12,
            rtol=0,
        )

    certificates = []
    for host in boards:
        certificates.append(
            dict(
                cells=host.cells,
                blank=energy_certificate(blank_generator(host, -1, private=True)),
                chiral=[
                    dict(
                        sign=sign,
                        ordinary_private_reference=private,
                        support=energy_certificate(
                            written_generator(host, face, sign, private), face, 8
                        ),
                    )
                    for sign in (-1, 1)
                    for private in (False, True)
                ],
            )
        )
    return dict(
        native_cells=board.cells,
        marked_face=face,
        contrast=8,
        flight_physical_seams=len(board.seams),
        clock_step=STEP,
        clock_bins=[COARSE, FINE],
        shots=SHOTS,
        seed=SEED,
        cases=cases,
        support=certificates,
        blank_neutral_formation=blank_neutral_certificate(),
        filled_memory_receiver=filled_memory_receiver(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    path = parser.parse_args().output
    result = run()
    if path:
        with path.open("x") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result["support"], indent=2))
