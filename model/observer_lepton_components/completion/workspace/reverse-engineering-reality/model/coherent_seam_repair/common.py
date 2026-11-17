"""Exact local source operators shared by the reversible and committed rules."""
import sympy as s


def require(ok, message):
    if not bool(ok):
        raise ArithmeticError(message)


def local_pump_kraus():
    """Complete pump on (00,10,01,11,20,02), including both double loads."""
    plus = s.Matrix([0, 1, 1, 0, 0, 0]) / s.sqrt(2)
    minus = s.Matrix([0, 1, -1, 0, 0, 0]) / s.sqrt(2)
    keep = s.diag(1, 0, 0, 1, 0, 0) + plus * plus.T
    jump = plus * minus.T
    split_a, split_b = s.zeros(6), s.zeros(6)
    split_a[3, 4] = split_b[3, 5] = 1
    return keep, jump, split_a, split_b


def local_repair_unitary():
    """Reversible 24-dimensional completion on source times four record values."""
    source, record = s.eye(6), s.eye(4)
    plus = (source[:, 1] + source[:, 2]) / s.sqrt(2)
    minus = (source[:, 1] - source[:, 2]) / s.sqrt(2)
    pairs = [
        (s.kronecker_product(minus, record[:, 0]), s.kronecker_product(plus, record[:, 1])),
        (s.kronecker_product(source[:, 4], record[:, 0]), s.kronecker_product(source[:, 3], record[:, 2])),
        (s.kronecker_product(source[:, 5], record[:, 0]), s.kronecker_product(source[:, 3], record[:, 3])),
    ]
    result = s.eye(24)
    for left, right in pairs:
        difference = left - right
        result -= difference * difference.H
    return result
