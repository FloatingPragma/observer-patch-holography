"""Require exact receipts, allowing only controlled roundoff in computed real fields."""

from math import isclose


def compare(actual, expected, path):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for key in expected:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        assert isclose(actual, expected, rel_tol=1e-10, abs_tol=5e-9), (
            path,
            actual,
            expected,
        )
    else:
        assert type(actual) is type(expected) and actual == expected, (
            path,
            actual,
            expected,
        )
