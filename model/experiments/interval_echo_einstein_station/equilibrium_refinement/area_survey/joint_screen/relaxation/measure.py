"""Certify an information transfer produced by one interval of native monitored evolution."""

import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from mpmath import iv
from scipy.sparse.csgraph import connected_components

from .. import acquire, certificate as cert, measure, wave

HERE = Path(__file__).resolve().parent
DISPLACEMENT = Fraction(1, 2**19)
LOTTERY = 2**50
SHOTS = 2**56
BATCH = 2**46


def instruments(s):
    im = measure.instrument(s)
    numerator, error = wave.kernel(s.h)
    old, old_error, _, joint, joint_error = cert.estimates(s, im, numerator, error)
    rays = (*measure.EDGES, (1, 0), (2, 0), (3, 0))
    starts = s.starts[:, [measure.clock.RAYS.index(ray) for ray in rays]].astype(int)
    unit = cert.M
    energy = np.rint(2 * s.h.diagonal()).astype(int).astype(object) * (unit // 512)
    rows = np.vstack(
        [
            energy,
            starts.T.astype(object) * unit,
            [old[ray] for ray in measure.clock.RAYS],
            256 * joint,
            (s.cut_count == 2).astype(object) * unit,
        ]
    )
    return rows, max(old_error, 256 * joint_error), numerator, error


def pushed_rows(s, rows, row_error, kernel, kernel_error):
    """Integer enclosures for the original instrument applied to (P-I)p."""
    unit, den = cert.M, wave.UNIT**2
    shifted = (rows @ kernel + den // 2) // den - rows
    # P symmetric stochastic, so ||(P-I)a||_infty <= 2||a||_infty.
    error = (
        2 * row_error
        + Fraction(max(abs(rows).flat), unit) * kernel_error
        + Fraction(1, 2 * unit)
    )
    return np.vstack([np.full(len(s.words), unit, dtype=object), shifted]), error


def certified_direction(matrix, error):
    unit, scale = cert.M, cert.B_UNIT
    b = cert.quantize(
        np.linalg.pinv(np.asarray(matrix, dtype=float) / unit, rcond=1e-13), scale
    )
    identity = np.eye(len(matrix), dtype=int).astype(object) * unit * scale
    residual = matrix @ b - identity
    defect = (
        Fraction(max(abs(residual).sum(axis=1)), unit * scale)
        + Fraction(sum(abs(b).flat), scale) * error
    )
    if defect >= 1:
        raise ArithmeticError(
            "The proposed right inverse does not certify full row rank"
        )
    last = b[:, -1]
    last_error = (
        Fraction(max(abs(residual[:, -1])), unit * scale)
        + Fraction(sum(abs(last)), scale) * error
    )
    correction = Fraction(max(abs(b).sum(axis=1)), scale) * last_error / (1 - defect)
    return last, defect, correction


def info_interval(s, tickets, kernel, kernel_error):
    """Bulk Shannon change for the deposited tickets, directly enclosing every marginal."""
    total, den = int(sum(tickets)), wave.UNIT**2
    next_num = kernel @ tickets
    groups = {}
    for idx, word in enumerate(s.words):
        key = tuple(v for v in word if v not in s.rim)
        groups.setdefault(key, []).append(idx)
    delta = iv.mpf(0)
    for group in groups.values():
        old = iv.mpf(int(sum(tickets[group]))) / total
        after = cert.exact_read(sum(next_num[group]), total * den, kernel_error)
        delta += old * iv.log(old) - after * iv.log(after)
    return acquire.ends(delta)


def record_counts(s, p, rng):
    """Fresh interior word and first-ray launch flags; all attempts enter their denominators."""
    labels = {}
    word_ids = [
        labels.setdefault(tuple(v for v in w if v not in s.rim), len(labels))
        for w in s.words
    ]
    probs = np.bincount(word_ids, weights=p)
    bulk = rng.multinomial(BATCH, probs / probs.sum(), size=SHOTS // BATCH).sum(axis=0)
    rate = float(p @ s.starts[:, measure.clock.RAYS.index((0, 1))])
    return dict(
        bulk=bulk.tolist(),
        herald=int(rng.binomial(BATCH, rate, size=SHOTS // BATCH).sum()),
    )


def run():
    iv.dps = 42
    s = measure.clock.survey(1)
    rows, row_error, kernel, kernel_error = instruments(s)
    matrix, error = pushed_rows(s, rows, row_error, kernel, kernel_error)
    last, defect, correction = certified_direction(matrix, error)
    n = len(last)
    max_direction = Fraction(max(abs(last)), cert.B_UNIT) + correction
    halfwidth = 1 / (n * max_direction)
    assert DISPLACEMENT < halfwidth
    # Pin the lab recipe independently of platform-specific least-squares roundoff.
    path = HERE / "lottery.json"
    if not path.exists():
        raise FileNotFoundError("The deposited 50-bit source lottery is required")
    tickets = np.asarray(json.loads(path.read_text()), dtype=object)
    assert len(tickets) == n and sum(tickets) == LOTTERY and min(tickets) > 0
    ticket_error = max(
        abs(
            Fraction(int(k), LOTTERY)
            - Fraction(1, n)
            - DISPLACEMENT * Fraction(int(v), cert.B_UNIT)
        )
        for k, v in zip(tickets, last, strict=True)
    )
    exact_source_error = ticket_error + DISPLACEMENT * correction
    # The measured difference on the ticket arm uses mean-zero weights. P preserves uniform.
    centered = n * tickets - LOTTERY
    scale = n * LOTTERY
    delta = rows @ ((kernel @ centered) - wave.UNIT**2 * centered)
    uncertainty = (
        2 * row_error + Fraction(max(abs(rows).flat), cert.M) * kernel_error
    ) * Fraction(sum(abs(centered)), scale)
    bounds = [
        Fraction(abs(v), cert.M * wave.UNIT**2 * scale) + uncertainty for v in delta
    ]
    rate_change = (
        s.starts[:, 0].astype(object)
        @ kernel
        @ (kernel @ centered - wave.UNIT**2 * centered)
    )
    # Two applications of the stochastic kernel; the approximate product can add an e^2 term.
    lag_err = (3 * kernel_error + kernel_error**2) * Fraction(sum(abs(centered)), scale)
    late = cert.exact_read(rate_change, wave.UNIT**4 * scale, lag_err)
    entropy = info_interval(s, tickets, kernel, kernel_error)
    p = np.array(tickets, dtype=float) / LOTTERY
    rng = np.random.default_rng(49235)
    banks = [record_counts(s, p, rng)]
    transition = np.asarray(kernel, dtype=float) / float(wave.UNIT**2)
    for _ in range(2):
        p = transition @ p
        banks.append(record_counts(s, p, rng))
    age_delta = (banks[2]["herald"] - banks[1]["herald"]) / SHOTS
    # Simultaneous binomial Hoeffding error over every bin on all three fresh cohorts.
    n_bins = 3 * (1 + len(banks[0]["bulk"]))
    eps = np.sqrt(np.log(2 * n_bins / 0.01) / (2 * SHOTS)) + 1e-12
    radius = 2 * eps
    entropies = []
    for bank in banks[:2]:
        rates = np.array(bank["bulk"]) / SHOTS
        assert min(rates) > eps and max(rates) + eps < 1 / np.e
        entropies.append(
            [-float(np.sum((rates + d) * np.log(rates + d))) for d in (-eps, eps)]
        )
    bulk_interval = [
        entropies[1][0] - entropies[0][1],
        entropies[1][1] - entropies[0][0],
    ]
    # Both stocks are normalized. Centering the entropy gradient costs half its
    # oscillation on each arm; marginalization and P contract the 1-norm.
    min_word = Fraction(min(tickets), LOTTERY) - exact_source_error
    assert min_word > Fraction(1, 1000)
    max_bulk = 6 * (Fraction(max(tickets), LOTTERY) + exact_source_error)
    ratio = max_bulk / min_word
    width = iv.log(iv.mpf(ratio.numerator) / ratio.denominator)
    entropy_correction = acquire.ends(width * cert.up(n * exact_source_error))[1]
    native2 = np.rint(2 * s.h).astype(int)
    edges = native2 < 0
    # Certified positive diagonal and every native exchange edge imply a single aperiodic class.
    assert connected_components(edges, directed=False)[0] == 1
    diagonal_lower = Fraction(min(kernel.diagonal()), wave.UNIT**2) - kernel_error
    exchange_lower = Fraction(min(kernel[edges]), wave.UNIT**2) - kernel_error
    assert diagonal_lower > 0 and exchange_lower > 0
    return dict(
        number_of_words=n,
        step=float(DISPLACEMENT),
        native_monitor_interval=0.5,
        native_exchange_probability_lower=-cert.up(-exchange_lower),
        stay_probability_lower=-cert.up(-diagonal_lower),
        rank=len(matrix),
        right_inverse_defect_upper=cert.up(defect),
        reachable_halfwidth_lower=-cert.up(-halfwidth),
        max_correction_per_ticket_probability=cert.up(exact_source_error),
        fixed_record_changes_on_deposited_stock=dict(
            mean_energy_upper=256 * cert.up(bounds[0]),
            herald_upper=3 * cert.up(max(bounds[1:10])),
            return_upper=cert.up(max(bounds[10:22])),
            pair_pixel_upper=cert.up(max(bounds[22:34])) / 256,
        ),
        interior_entropy_change_deposited_nats=entropy,
        interior_entropy_change_exact_stock_lower_nats=float(
            np.nextafter(entropy[0] - entropy_correction, -np.inf)
        ),
        first_variation_nats=float(4 * np.log(4) / 3 - np.log(6)),
        after_one_further_interval_first_ray_change=acquire.ends(late),
        estimated_late_herald_change=float(age_delta),
        simultaneous_99_counted_entropy_change_nats=bulk_interval,
        simultaneous_99_late_herald_interval=[
            float(age_delta - radius),
            float(age_delta + radius),
        ],
        shots_per_setting=SHOTS,
        ages_in_half_ticks=[0, 1, 2],
        counts=banks,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
