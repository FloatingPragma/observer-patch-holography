"""Paired reunion counts on completed transfer words using the candidate closing instrument."""
import hashlib
import itertools
import math
import random

import numpy as np

from .clock import closing_tick


# The first schedules reproduce the calibrated domain. The others expose the two age boundaries.
SCHEDULES = {
    "alternating_2_6": (2, 6) * 8,
    "alternating_6_2": (6, 2) * 8,
    "uniform_8": (8,) * 8,
    "mixed_8_15": (8, 15) * 4,
    "mixed_15_8": (15, 8) * 4,
    "overflow_16_then_2": (16, 2) * 4,
}
PREPARATIONS = 8192
SEED = 2026112127
FIELDS = ("home", "traveller_3bit", "traveller_4bit")


def suppression(q, width):
    """Probability mass at the first return, propagated independently on signed walk positions."""
    if q >= 2**width:
        return 0.
    if q == 1:
        return 1.
    choices = 2 ** (q - 1).bit_length()
    match = (1 - (1 - q / choices)**12) / q
    live = {-1: .5, 1: .5}
    suppressed = 0.
    for _ in range(11):
        following = {}
        for position, mass in live.items():
            for destination in (position - 1, position + 1):
                branch = mass * match * .5
                if destination == 0:
                    suppressed += branch
                else:
                    following[destination] = following.get(destination, 0.) + branch
        live = following
    return suppressed


def ratio_with_error(readings, home):
    ratio = float(readings.sum() / home.sum())
    # Independent rows, with covariance of numerator and denominator retained within each row.
    residual = readings - ratio * home
    error = np.std(residual, ddof=1) / (math.sqrt(len(home)) * home.mean())
    return {"ratio": ratio, "standard_error": float(error)}


def sample(seed=SEED):
    records, panels = {}, {}
    for index, (name, schedule) in enumerate(SCHEDULES.items()):
        source = np.random.default_rng(np.random.SeedSequence([seed, index, 0]))
        coins = random.Random(int(np.random.SeedSequence([seed, index, 1]).generate_state(1)[0]))

        def fair_stream():
            while True:
                yield coins.getrandbits(1)

        words = np.zeros((PREPARATIONS, len(FIELDS)), dtype="<u2")
        changed = 0
        for q in schedule:
            # Independent primary outcomes. The two widths see the same completed and missed services.
            home = source.integers(2, size=(PREPARATIONS, q), dtype=np.uint8)
            moving = source.integers(2, size=(PREPARATIONS, q), dtype=np.uint8)
            words[:, 0] += home.sum(axis=1, dtype=np.uint16)
            words[:, 1:] += moving[:, :-1].sum(axis=1, dtype=np.uint16)[:, None]
            for row, primary in enumerate(moving[:, -1]):
                # The same fair bit at each draw position, even if a branch consumes fewer bits.
                short_tape, long_tape = itertools.tee(fair_stream())
                narrow = closing_tick(q, int(primary), short_tape.__next__, age_bits=3)
                wide = closing_tick(q, int(primary), long_tape.__next__, age_bits=4)
                words[row, 1] += int(narrow)
                words[row, 2] += int(wide)
                changed += int(narrow != wide)

        records[name] = words
        h, narrow, wide = words.astype(np.int32).T
        total = sum(schedule)
        panels[name] = {
            "intervals": list(schedule), "completed_preparations": len(words),
            "paired_closing_decisions_changed": changed,
            "three_bits": ratio_with_error(narrow, h), "four_bits": ratio_with_error(wide, h),
            "paired_four_minus_three": ratio_with_error(wide - narrow, h),
            "finite_expectation": {str(w): 1 - sum(suppression(q, w) for q in schedule) / total
                                   for w in (3, 4)},
            "calibration_curve": sum(math.sqrt(q*q - 1) for q in schedule) / total,
            "words_sha256": hashlib.sha256(words.tobytes()).hexdigest(),
        }
    return {"seed": seed, "fields": list(FIELDS), "panels": panels}, records
