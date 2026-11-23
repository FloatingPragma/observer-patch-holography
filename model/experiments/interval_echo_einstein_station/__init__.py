"""A finite observer measures a regional entropy response across a physical cut."""

import math
from numbers import Integral, Real

# Record fields holding multinomial draws (or sums of draws) from seeded
# generators. NumPy's samplers use floating-point arithmetic, so a draw can
# move by a unit between platforms while the seed, the law and the total stay
# fixed. Every other integer (cells, ports, classes, flags, shot budgets) is exact.
SAMPLED_FIELDS = frozenset({"reads", "readings", "counts", "selected_shots", "total_shots",
                            "flag_counts", "interior_counts"})


def _mismatches(actual, stored, path, sampled, rel, absolute):
    if isinstance(stored, dict):
        if not isinstance(actual, dict) or actual.keys() != stored.keys():
            yield f"{path}: keys differ"
            return
        for key, value in stored.items():
            yield from _mismatches(actual[key], value, f"{path}.{key}",
                                   sampled or key in SAMPLED_FIELDS, rel, absolute)
    elif isinstance(stored, list):
        if not isinstance(actual, (list, tuple)) or len(actual) != len(stored):
            yield f"{path}: length differs"
            return
        for index, (left, right) in enumerate(zip(actual, stored)):
            yield from _mismatches(left, right, f"{path}[{index}]", sampled, rel, absolute)
    elif isinstance(stored, float) and isinstance(actual, Real) and not isinstance(actual, (bool, Integral)):
        if not math.isclose(float(actual), stored, rel_tol=rel, abs_tol=absolute):
            yield f"{path}: {actual!r} != {stored!r}"
    elif (sampled and isinstance(stored, int) and not isinstance(stored, bool)
          and isinstance(actual, Integral) and not isinstance(actual, bool)):
        if abs(int(actual) - stored) > max(3, 1e-9 * abs(stored)):
            yield f"{path}: {actual!r} != {stored!r}"
    elif actual != stored:
        yield f"{path}: {actual!r} != {stored!r}"


def assert_same_replay(actual, stored, rel=1e-6, absolute=1e-12):
    """Require a fresh replay to reproduce its stored record.

    Structure, keys, labels, flags and counts agree exactly, except that a
    sampled draw (SAMPLED_FIELDS) may move by max(3, 1e-9 * count) units.
    Floating readouts agree to a relative 1e-6 (absolute 1e-12 near zero).
    Both allowances cover the last digits of the platform's floating-point
    libraries; every mismatch is listed.
    """
    found = list(_mismatches(actual, stored, "$", False, rel, absolute))
    assert not found, f"{len(found)} replay mismatch(es): " + "; ".join(found[:12])
