"""All-time energy-spread bound for the charged-level drop into either neutral orbital."""

from math import sqrt

import numpy as np
from .apparatus import frame_columns, local_frames, sector


def certificate():
    words, h = sector(8, 2, 1, 8)
    blocks, frames = local_frames(8)
    initial, si = frame_columns(words, [(0, 1), (3, 4)], frames)
    final, fi = frame_columns(words, [(0, 1), (5, 6)], frames)
    expected = np.kron(blocks[0], np.eye(2)) + np.kron(np.eye(2), blocks[1])
    for indices in (si, fi):
        np.testing.assert_array_equal(h[np.ix_(indices, indices)], expected)
        # Each pair configuration can hop out through two outer edges. A vacancy
        # on the second charged rail can also hop into the empty buffer.
        leaked = h[:, indices].copy()
        leaked[indices, :] = 0
        np.testing.assert_array_equal(
            leaked.T @ leaked, np.diag([0.5, 0.5, 0.75, 0.75])
        )
    cap = final[:, :2]
    gap = sqrt(293) / 2
    p = (1 - 17 / sqrt(293)) / 2
    variance = 0.5 + (1 - p) / 4
    escape_sq = 0.5 + p / 4
    s = initial[:, 3]
    np.testing.assert_allclose(np.linalg.norm(h @ s - (s @ h @ s) * s) ** 2, variance)
    hd = cap.T @ h @ cap
    residual = h @ cap - cap @ hd
    np.testing.assert_allclose(residual.T @ residual, escape_sq * np.eye(2), atol=1e-12)
    np.testing.assert_allclose(s @ h @ s - np.linalg.eigvalsh(hd)[-1], gap)
    upper = (sqrt(variance) + sqrt(escape_sq)) ** 2 / gap**2
    return dict(
        source_full_variance=round(variance, 12),
        daughter_leakage_squared=round(escape_sq, 12),
        local_gap=round(gap, 12),
        all_time_lower_charge_and_either_neutral_upper=round(upper, 12),
    )


def receiver_envelope():
    """All terminal neutral words, and the receiver-length-independent energy envelope."""
    wn, hn = sector(8, 1, 1, 8)
    neutral_ids = [j for j, w in enumerate(wn) if w[:5] == (1,) * 5]
    remote = hn[np.ix_(neutral_ids, neutral_ids)]
    # Weighted Gershgorin, using positive integer weights. Each orientation
    # of the high adjacent pair has the same weight.
    weights = np.array(
        [
            200 if x[-1] == 1 else 5 if x[-2] == 1 else 1
            for x in (wn[j] for j in neutral_ids)
        ]
    )
    # AC (weight 5) sees both the high pair (200) and the end pair (1).
    # The end pair (1) sees AC (5) and its other orientation.
    twice = np.rint(2 * remote).astype(int)
    row2 = np.diag(twice) * weights + (abs(twice) - np.diag(np.diag(twice))) @ weights
    from fractions import Fraction

    ceiling = max(
        Fraction(int(v), 2 * int(z)) for v, z in zip(row2, weights, strict=True)
    )
    w, h = sector(8, 2, 1, 8)
    _, frames = local_frames(8)
    source, _ = frame_columns(w, [(0, 1), (3, 4)], frames)
    cap = np.zeros((len(w), 6))
    for col, j in enumerate(neutral_ids):
        for pos in (0, 1):
            word = list(wn[j])
            word[pos] = 0
            cap[w.index(tuple(word)), col] = frames[0][pos, 0]
    charge = (51 - sqrt(293)) / 4
    hd = cap.T @ h @ cap
    np.testing.assert_allclose(hd, remote + charge * np.eye(6), atol=1e-10)
    leakage = h @ cap - cap @ hd
    slow_weight = (1 - 17 / sqrt(293)) / 2
    # Leaving the receiving owner's first port requires a non-unit load there.
    gram = np.diag(
        [(0.25 if wn[i][5] != 1 else 0) + slow_weight / 4 for i in neutral_ids]
    )
    np.testing.assert_allclose(leakage.T @ leakage, gram, atol=1e-12)
    spread_sq = 0.5 + (1 - slow_weight) / 4
    leak_sq = 0.25 + slow_weight / 4
    double_gap = sqrt(293) / 2

    def limit(maximum):
        return (sqrt(spread_sq) + sqrt(leak_sq)) ** 2 / (50 + double_gap - maximum) ** 2

    s = source[:, 3]
    return dict(
        two_times_neutral_matrix=twice.tolist(),
        integer_row_weights=weights.tolist(),
        neutral_upper_fraction=str(ceiling),
        neutral_energies=np.linalg.eigvalsh(remote).round(11).tolist(),
        raised_source_mean=round(float(s @ h @ s), 11),
        all_word_leakage_squared=round(leak_sq, 12),
        all_three_port_words_probability_upper=round(limit(float(ceiling)), 12),
        any_length_single_boundary_upper=round(limit(51), 12),
    )
