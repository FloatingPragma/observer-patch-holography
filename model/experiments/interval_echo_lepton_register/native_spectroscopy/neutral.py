"""Exact formation floor for detectable neutral departures in the blank host."""

from fractions import Fraction

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_11_spectrum.measure import (
    word_table,
)
from model.executable_model import cell

from .source import KAPPA


def blank_neutral_certificate():
    """Use the owner bound and sharpen the low inventories by exact enumeration."""
    _, _, words = word_table(12)
    loads = words.astype(np.int16) - 1
    counts_zero = np.count_nonzero(loads == -1, axis=1)
    counts_two = np.count_nonzero(loads == 1, axis=1)
    contrast = sum((loads[:, a] - loads[:, b]) ** 2 for a, b in cell().seams)
    # Dropping the three largest private costs minimizes over all gateway layouts.
    removed = 3 * np.minimum(counts_two, 3) - np.maximum(0, counts_zero - 9)
    cost = contrast + 3 * counts_two - counts_zero - removed
    vacancies, excesses = counts_zero, counts_two
    signed = excesses - vacancies
    defects = vacancies + excesses
    assert np.min(4 * cost - 6 * signed - 3 * defects) == 0
    local = np.full((13, 13), 10_000, dtype=int)
    np.minimum.at(local, (vacancies, excesses), cost)

    # Six pairs already cost at least 72 by the linear bound. Enumerate smaller
    # neutral inventories across any number of nonempty owners exactly.
    size = 6
    floor = np.full((size, size), 10_000, dtype=int)
    floor[0, 0] = 0
    for n in range(1, 2 * size - 1):
        for z in range(size):
            g = n - z
            if 0 <= g < size:
                floor[z, g] = min(
                    local[i, j] + floor[z - i, g - j]
                    for i in range(z + 1)
                    for j in range(g + 1)
                    if i + j
                )
    low_inventory = [KAPPA * int(floor[i, i]) for i in range(1, size)]
    high_inventory = 2 * size * 3 * KAPPA // 4
    neutral_floor = min(*low_inventory, high_inventory)
    # One sharp vacancy at a degree-five port bounds the lightest charge above.
    negative_trial = Fraction(8 * KAPPA + 5, 2)
    assert neutral_floor > 2 * negative_trial
    return dict(
        owner_words=len(words),
        small_neutral_formation_floors=low_inventory,
        floor_from_six_pairs=high_inventory,
        neutral_formation_floor=neutral_floor,
        lowest_negative_formation_trial_ceiling=str(negative_trial),
        neutral_to_lightest_negative_ratio_floor=str(neutral_floor / negative_trial),
    )


def filled_memory_receiver():
    """Check whether a source write can change a neighboring coin on filled material."""
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    local = 8 * (exchange + KAPPA * contrast)
    integer = np.rint(local).astype(np.int64)
    np.testing.assert_allclose(local, integer, atol=1e-12, rtol=0)
    loads = law.SRC.index((1, 1))
    filled = loads + len(law.SRC) * np.arange(11)
    other = np.setdiff1d(np.arange(99), filled)
    assert not integer[np.ix_(other, filled)].any()
    memory = integer[np.ix_(filled, filled)] / 8
    two_seams = np.kron(memory, np.eye(11)) + np.kron(np.eye(11), memory)
    coin_0 = 1 + law.LABELS.index((0, 2))
    coin_1 = 1 + law.LABELS.index((1, 2))
    read_coin_1 = [1 + i for i, (coin, _) in enumerate(law.LABELS) if coin == 1]
    first, second = np.eye(11)[:, [coin_0, coin_1]].T
    # The third input is the detector's independent local-response control.
    inputs = np.stack(
        (np.kron(first, first), np.kron(second, first), np.kron(first, second)), axis=1
    )
    rng = np.random.default_rng(726112)
    shots = 2**20
    banks = []
    for time in (0, 1 / 32, 1 / 8, 3 / 8, 7 / 8):
        released = (expm(-1j * time * two_seams) @ inputs).reshape(11, 11, 3)
        probability = np.sum(abs(released[:, read_coin_1]) ** 2, axis=(0, 1))
        np.testing.assert_allclose(probability[0], probability[1], atol=2e-14)
        counts = rng.binomial(shots, np.clip(probability, 0, 1))
        banks.append(
            dict(
                wait=time,
                probabilities=probability.tolist(),
                coin_one_counts=counts.tolist(),
            )
        )
    assert abs(banks[1]["probabilities"][0] - banks[1]["probabilities"][2]) > 0.1
    return dict(
        filled_columns_checked=len(filled),
        unfilled_matrix_elements=0,
        shots=shots,
        seed=726112,
        cases=["source coin 0", "source coin 1", "receiver switched locally"],
        banks=banks,
    )
