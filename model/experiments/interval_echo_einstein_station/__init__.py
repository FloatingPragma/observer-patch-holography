"""A finite observer measures a regional entropy response across a physical cut."""

import math
from numbers import Integral, Real


def assert_same_replay(actual, stored, path="$", rel=1e-9, absolute=1e-12):
    """Require a fresh replay to reproduce its stored record.

    Structure, keys, counts, labels and flags must agree exactly. Floating
    readouts may differ in their last digits between BLAS and libm builds, so
    they agree to a relative 1e-9 (absolute 1e-12 near zero).
    """
    if isinstance(stored, dict):
        assert isinstance(actual, dict) and actual.keys() == stored.keys(), f"{path}: keys differ"
        for key, value in stored.items():
            assert_same_replay(actual[key], value, f"{path}.{key}", rel, absolute)
    elif isinstance(stored, list):
        assert isinstance(actual, (list, tuple)) and len(actual) == len(stored), f"{path}: length differs"
        for index, (left, right) in enumerate(zip(actual, stored)):
            assert_same_replay(left, right, f"{path}[{index}]", rel, absolute)
    elif isinstance(stored, float) and isinstance(actual, Real) and not isinstance(actual, (bool, Integral)):
        assert math.isclose(float(actual), stored, rel_tol=rel, abs_tol=absolute), f"{path}: {actual!r} != {stored!r}"
    else:
        assert actual == stored, f"{path}: {actual!r} != {stored!r}"
