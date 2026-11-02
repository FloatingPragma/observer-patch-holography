"""Finite exact checks of the spin and rotation identities used in the proof."""
from functools import lru_cache

import sympy as sp


@lru_cache(maxsize=1)
def _checks():
    i = sp.I
    spins = (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -i], [i, 0]]) / 2,
        sp.diag(1, -1) / 2,
    )
    left = [sp.kronecker_product(s, sp.eye(2)) for s in spins]
    right = [sp.kronecker_product(sp.eye(2), s) for s in spins]
    bond = sum((left[a] * right[a] for a in range(3)), sp.zeros(4))
    x, y = sp.symbols('x y', real=True)
    comm = lambda a, b: a * b - b * a
    double = []
    for axis in range(3):
        a = x * left[axis] + y * right[axis]
        expected = -(x - y)**2 * sum(
            (left[b] * right[b] for b in range(3) if b != axis), sp.zeros(4)
        )
        double.append(sp.simplify(comm(a, comm(bond, a)) - expected) == sp.zeros(4))
    transverse = left[0] * right[0] + left[2] * right[2]
    spectrum = transverse.eigenvals()
    c, d = sp.symbols('c d', real=True)
    covariance = sp.diag(c, d, d)
    generators = (
        sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]]),
        sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]]),
    )
    rotation = sp.Matrix([[sp.cos(x), 0, sp.sin(x)], [0, 1, 0],
                          [-sp.sin(x), 0, sp.cos(x)]])
    one_axis = sp.simplify(sp.trace(covariance * rotation) - sp.trace(covariance)
                           - (-(c + d)) * (1 - sp.cos(x))) == 0
    transverse_quadratic = [sp.trace(covariance * g.T * g) == c + d
                            for g in generators[1:]]
    cross = sp.trace(covariance * (generators[1].T * generators[2]
                                  + generators[2].T * generators[1])) == 0
    canonical = comm(spins[2], -spins[1]) == i * spins[0]
    return {
        'all_three_double_commutators': all(double),
        'transverse_bond_eigenvalues': {str(k): int(v) for k, v in spectrum.items()},
        'transverse_bond_norm_one_half': max(abs(k) for k in spectrum) == sp.Rational(1, 2),
        'exact_single_axis_rotation_energy': bool(one_axis),
        'both_transverse_quadratic_coefficients': all(transverse_quadratic),
        'mixed_quadratic_rotation_term_zero': bool(cross),
        'canonical_sign': bool(canonical),
    }


def exact_checks():
    """Return detached results so callers cannot poison the cached computation."""
    import copy
    result = copy.deepcopy(_checks())
    if not all(v for k, v in result.items() if k != 'transverse_bond_eigenvalues'):
        raise ArithmeticError('an exact spin identity failed')
    return result
