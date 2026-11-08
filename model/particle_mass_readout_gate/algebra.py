"""Exact Weyl-symbol identity for every retained kinetic valley."""
from itertools import product

import sympy as s


def casimir_checks():
    momentum = s.symbols('p1 p2 p3', real=True)
    speed = s.symbols('v', positive=True)
    pauli = (s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -s.I], [s.I, 0]]),
             s.diag(1, -1))
    output = []
    for valley in product((0, 1), repeat=3):
        symbol = -speed*sum(((-1)**b*p*sigma for b, p, sigma in
                            zip(valley, momentum, pauli, strict=True)), s.zeros(2))
        residual = (symbol*symbol-speed**2*sum(p*p for p in momentum)*s.eye(2)).applyfunc(s.expand)
        if residual != s.zeros(2):
            raise ArithmeticError('Weyl Casimir identity failed')
        output.append({'valley': list(valley), 'symbol_squared_equals_v2_p2': True})
    return output
