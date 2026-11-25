"""Integrity checks for archived receipts and independent high-exposure counts."""

import hashlib
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR

import numpy as np


def pinned_receipt(path, expected_sha256):
    payload = path.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == expected_sha256
    return json.loads(payload)


def verify_recorded_band(frozen, live, significant_digits=8):
    """Compare independently certified intervals at stated decimal precision."""
    lo, hi = json.loads(frozen, parse_float=Decimal)
    lower, upper = json.loads(live, parse_float=Decimal)
    assert lower <= upper and lo <= hi and max(lo, lower) <= min(hi, upper)

    def quantum(value):
        return Decimal(1).scaleb(value.adjusted() + 1 - significant_digits)

    assert lo.quantize(quantum(lo), rounding=ROUND_FLOOR) <= lower
    assert upper <= hi.quantize(quantum(hi), rounding=ROUND_CEILING)


def check_binomial(count, trials, probability, log_level):
    """Two-sided Bernstein check for one marginal of a predeclared family."""
    assert isinstance(count, int) and 0 <= count <= trials and 0 <= probability <= 1
    if probability in (0, 1):
        assert count == (trials if probability == 1 else 0)
        return
    expected = np.longdouble(trials) * probability
    variance = expected * (1 - probability)
    radius = np.sqrt(2 * variance * log_level) + 2 * log_level / 3
    assert abs(np.longdouble(count) - expected) <= radius
