"""A finite observer measures a regional entropy response across a physical cut."""

import math
from numbers import Integral, Real

# Counts at or above this size are multinomial draws over 2**38 to 2**48 shots.
# NumPy samples them with floating-point arithmetic, so their last units follow
# the platform's libm; smaller integers (cells, ports, classes, flags) are exact.
SAMPLED_COUNT = 2**32


def _mismatches(actual, stored, path, rel, absolute):
    if isinstance(stored, dict):
        if not isinstance(actual, dict) or actual.keys() != stored.keys():
            yield f"{path}: keys differ"
            return
        for key, value in stored.items():
            yield from _mismatches(actual[key], value, f"{path}.{key}", rel, absolute)
    elif isinstance(stored, list):
        if not isinstance(actual, (list, tuple)) or len(actual) != len(stored):
            yield f"{path}: length differs"
            return
        for index, (left, right) in enumerate(zip(actual, stored)):
            yield from _mismatches(left, right, f"{path}[{index}]", rel, absolute)
    elif isinstance(stored, float) and isinstance(actual, Real) and not isinstance(actual, (bool, Integral)):
        if not math.isclose(float(actual), stored, rel_tol=rel, abs_tol=absolute):
            yield f"{path}: {actual!r} != {stored!r}"
    elif (isinstance(stored, int) and not isinstance(stored, bool) and abs(stored) >= SAMPLED_COUNT
          and isinstance(actual, Integral) and not isinstance(actual, bool)):
        if abs(int(actual) - stored) > max(1, 1e-9 * abs(stored)):
            yield f"{path}: {actual!r} != {stored!r}"
    elif actual != stored:
        yield f"{path}: {actual!r} != {stored!r}"


def assert_same_replay(actual, stored, rel=1e-6, absolute=1e-12):
    """Require a fresh replay to reproduce its stored record.

    Structure, keys, labels, flags and ordinary counts agree exactly. Floating
    readouts agree to a relative 1e-6 (absolute 1e-12 near zero), and sampled
    counts of at least 2**32 to a relative 1e-9: both follow the platform's
    floating-point libraries in their last digits.
    """
    found = list(_mismatches(actual, stored, "$", rel, absolute))
    assert not found, f"{len(found)} replay mismatch(es): " + "; ".join(found[:12])
