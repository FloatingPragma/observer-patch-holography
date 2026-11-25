"""Prepare the rest line directly, impart momentum, and read energy and signed current."""

import argparse
import json

import numpy as np
from scipy.linalg import expm
from .measure import (
    HERE,
    LAW,
    binomial_interval,
    layout,
    prepare_seams,
    ring_hamiltonian,
    seam_data,
)
from .pulses import native_word

SHOTS = 2**30


def terminal(state, word, k, sign, rng):
    n = len(state)
    readings = np.zeros(3)
    pair = [LAW.SRC.index((1 + sign, 1)), LAW.SRC.index((1, 1 + sign))]
    local = seam_data(k)[0][np.ix_(pair, pair)]
    quarter = expm(-1j * np.pi / 2 * local)
    for a in range(n):
        b, c = (a + 1) % n, (a + 2) % n
        tested = native_word(state[[a, b, c]], word, k, sign, inverse=True)
        readings[0] += abs(tested[1]) ** 2 / n
        readings[1:] += abs(quarter @ state[[a, b]]) ** 2 / n
    # Take the occupation word on separate copies to measure the contrast,
    # and choose a uniformly random single addressed bond in the two other banks.
    location = rng.multinomial(SHOTS, abs(state) ** 2 / np.vdot(state, state).real)
    words = np.ones((n, n), int)
    words[np.diag_indices(n)] += sign
    contrast = (
        k
        * sum(
            np.sum((words[j] - np.roll(words[j], 1)) ** 2) * int(c)
            for j, c in enumerate(location)
        )
        / SHOTS
    )
    parity = int(rng.binomial(SHOTS, np.clip(readings[0], 0, 1)))
    left, right = readings[1:]
    cur = rng.multinomial(SHOTS, [left, right, 1 - left - right])
    # All eight native energy intervals use a union bound. The altered-law word
    # is a fixed-pulse counterfactual and its inverse measures a different basis.
    interval = contrast + n * binomial_interval(parity, SHOTS, 0.01 / 8)
    return dict(
        position_counts=location.tolist(),
        bond_minus_and_other=[parity, SHOTS - parity],
        current_left_right_other=cur.tolist(),
        coordinate_current=n * (int(cur[1]) - int(cur[0])) / (2 * SHOTS),
        load_contrast_energy=contrast,
        energy_from_native_bond=contrast + n * parity / SHOTS,
        simultaneous_energy_interval=interval.tolist(),
    )


def run():
    board, graph, ring, cuts, ref, donor = layout()
    n = len(ring)
    words = json.loads((HERE / "pulses.json").read_text())
    sharp = np.eye(n, dtype=complex)[:, 1]
    rng = np.random.default_rng(70733)
    cases = []
    for sign in (-1, 1):
        prep = prepare_seams(board, graph, ring, cuts, ref, donor, sign, True)
        gen, _, _ = ring_hamiltonian(board, ring, prep, 8, sign)
        first = native_word(sharp, words["0"], 8, sign)
        first = (
            expm(-1j * gen * 3) @ first
        )  # the selected charge waits freely at rest before the impulse
        records = []
        for m in (0, 1, 2):
            state = first
            if m:
                state = native_word(first, words["0"], 8, sign, inverse=True)
                state = native_word(state, words[str(m)], 8, sign)
            expected = np.exp(2j * np.pi * m * np.arange(n) / n) / np.sqrt(n)
            error = np.linalg.norm(state - expected * np.exp(-3j * 16))
            assert error < 3e-11
            test = native_word(
                np.eye(3, dtype=complex)[:, 1], words["bond_minus"], 8, sign
            )
            np.testing.assert_allclose(
                test, np.array([1, -1, 0]) / np.sqrt(2), atol=3e-12
            )
            record = terminal(state, words["bond_minus"], 8, sign, rng)
            record.update(
                mode=m,
                word_norm_error=error,
                free_eigen_residual=float(
                    np.linalg.norm(
                        gen @ state - (17 - np.cos(2 * np.pi * m / n)) * state
                    )
                ),
            )
            records.append(record)
        # Use the independently acquired energy intervals, with the same measured
        # zero-momentum level, to test one speed at both lattice momenta.
        E0_lo, E0_hi = records[0]["simultaneous_energy_interval"]
        for r in records[1:]:
            p2 = 2 * (1 - np.cos(2 * np.pi * r["mode"] / n))
            lo, hi = r["simultaneous_energy_interval"]
            r["mass_shell_speed_squared_interval"] = [
                (lo**2 - E0_hi**2) / p2,
                (hi**2 - E0_lo**2) / p2,
            ]
        changed = native_word(sharp, words["0"], 4, sign)
        # Complex-conjugating each real-Hamiltonian pulse time in the SAME
        # order launches the opposite momentum at exactly the same energy.
        opposite = native_word(first, words["0"], 8, sign, inverse=True)
        opposite = native_word(
            opposite, [(-t) % (2 * np.pi) for t in words["1"]], 8, sign
        )
        direction_ref = np.exp(-2j * np.pi * np.arange(n) / n - 3j * 16) / np.sqrt(n)
        np.testing.assert_allclose(opposite, direction_ref, atol=3e-11)
        negative_read = terminal(
            opposite, words["bond_minus"], 8, sign, np.random.default_rng(93032 + sign)
        )
        # Reuse the source word exactly and interrogate raw positions, independent of bond calibration.
        altered = rng.multinomial(
            SHOTS, abs(changed) ** 2 / np.vdot(changed, changed).real
        )
        cases.append(
            dict(
                sign=sign,
                rest_wait=3,
                records=records,
                opposite_first_momentum=negative_read,
                unchanged_pulses_k4_position=altered.tolist(),
                changed_law_uniform_overlap=float(abs(sum(changed) / np.sqrt(n)) ** 2),
            )
        )
    return dict(shots=SHOTS, ring=ring, cases=cases)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2) + "\n"
    if p.parse_args().write:
        (HERE / "rest_kick.json").write_text(data)
    else:
        print(data, end="")
