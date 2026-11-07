"""Exact finite witnesses for the countable-effect extension theorem.

The all-domain theorem is proved in DENSE_EFFECTS.md and uses the existing
finite Busch--Gleason theorem. These finite calculations do not establish
an executable source's positivity, additivity or implementation equivalence.
"""
from __future__ import annotations

from collections.abc import Callable

import sympy as sp


def rational(value) -> sp.Rational:
    """Accept exact real rationals; do not rationalize measured floats."""
    if isinstance(value, bool):
        raise ValueError("an exact rational, not a Boolean, is required")
    result = sp.sympify(value)
    if result.is_Rational is not True:
        raise ValueError("an exact real rational is required")
    return result


def hermitian(value) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(value)
    if matrix.rows < 1 or matrix.rows != matrix.cols:
        raise ValueError("a nonempty square matrix is required")
    for entry in matrix:
        if entry in (sp.true, sp.false):
            raise ValueError("Boolean matrix entries are not admitted")
        real, imaginary = sp.expand_complex(entry).as_real_imag()
        rational(real)
        rational(imaginary)
    if matrix != matrix.H:
        raise ValueError("a Hermitian matrix over Q(i) is required")
    return matrix


def positive_semidefinite(value) -> bool:
    """Exact Schur elimination; zero diagonal pivots keep their full row."""
    matrix = sp.Matrix(hermitian(value))
    while matrix.rows:
        diagonal = [matrix[i, i] for i in range(matrix.rows)]
        if any(entry < 0 for entry in diagonal):
            return False
        positive = next((i for i, entry in enumerate(diagonal) if entry > 0), None)
        if positive is None:
            return matrix == sp.zeros(matrix.rows)
        indices = [positive] + [i for i in range(matrix.rows) if i != positive]
        matrix = matrix.extract(indices, indices)
        column = matrix[1:, :1]
        matrix = sp.simplify(matrix[1:, 1:] - column * column.H / matrix[0, 0])
    return True


def effect(value) -> sp.ImmutableMatrix:
    matrix = hermitian(value)
    if not positive_semidefinite(matrix) or not positive_semidefinite(sp.eye(matrix.rows) - matrix):
        raise ValueError("the matrix must satisfy 0 <= E <= I")
    return matrix


def order_sandwich(value, radius) -> bool:
    """Certify -radius I <= value <= radius I exactly."""
    matrix = hermitian(value)
    radius = rational(radius)
    if radius < 0:
        raise ValueError("a nonnegative radius is required")
    identity = sp.eye(matrix.rows)
    return positive_semidefinite(radius * identity - matrix) and positive_semidefinite(radius * identity + matrix)


def lipschitz_order_witness(first, second, delta):
    """The two rational order comparisons used to prove 1-Lipschitzness."""
    first, second = effect(first), effect(second)
    delta = rational(delta)
    if first.shape != second.shape or not order_sandwich(first - second, delta):
        raise ValueError("the two effects need the stated common order bound")
    identity = sp.eye(first.rows)
    scale = 1 + delta
    forward = (effect(first / scale), effect((second + delta * identity) / scale))
    reverse = (effect(second / scale), effect((first + delta * identity) / scale))
    if not all(positive_semidefinite(right - left) for left, right in (forward, reverse)):
        raise ArithmeticError("the rational order implication failed")
    return forward, reverse


def regularize_hermitian(approximation, delta) -> sp.ImmutableMatrix:
    """Turn -delta I <= B <= (1+delta)I into a rational effect."""
    approximation = hermitian(approximation)
    delta = rational(delta)
    if delta < 0:
        raise ValueError("a nonnegative approximation radius is required")
    identity = sp.eye(approximation.rows)
    if not positive_semidefinite(approximation + delta * identity) or not positive_semidefinite((1 + delta) * identity - approximation):
        raise ValueError("the approximation is outside the stated relaxed effect bounds")
    return effect((approximation + delta * identity) / (1 + 2 * delta))


def regularize_coexistent(first, second, delta):
    """Regularize two approximants jointly, keeping their sum an effect."""
    first, second = hermitian(first), hermitian(second)
    delta = rational(delta)
    if first.shape != second.shape or delta < 0:
        raise ValueError("matching matrices and nonnegative radius are required")
    identity = sp.eye(first.rows)
    if not all(positive_semidefinite(matrix + delta * identity) for matrix in (first, second)):
        raise ValueError("an approximant violates the lower order bound")
    if not positive_semidefinite((1 + 2 * delta) * identity - first - second):
        raise ValueError("the approximated sum violates the upper order bound")
    result = tuple(effect((matrix + delta * identity) / (1 + 4 * delta)) for matrix in (first, second))
    effect(result[0] + result[1])
    return result


def tomographic_candidate(valuation: Callable, dimension: int) -> sp.ImmutableMatrix:
    """Reconstruct a candidate; finite tomography does not certify additivity."""
    if type(dimension) is not int or dimension < 1:
        raise ValueError("a positive integer dimension is required")
    units = []
    for i in range(dimension):
        unit = sp.zeros(dimension)
        unit[i, i] = 1
        units.append(sp.ImmutableMatrix(unit))
    diagonal = [rational(valuation(unit)) for unit in units]
    candidate = sp.diag(*diagonal)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            real_probe = sp.Matrix((units[i] + units[j]) / 2)
            real_probe[i, j] = real_probe[j, i] = sp.Rational(1, 2)
            imaginary_probe = sp.Matrix((units[i] + units[j]) / 2)
            imaginary_probe[i, j], imaginary_probe[j, i] = sp.I / 2, -sp.I / 2
            centre = (diagonal[i] + diagonal[j]) / 2
            real = rational(valuation(effect(real_probe))) - centre
            imaginary = rational(valuation(effect(imaginary_probe))) - centre
            candidate[i, j], candidate[j, i] = real + sp.I * imaginary, real - sp.I * imaginary
    return hermitian(candidate)


def _entries(matrix):
    return [[str(sp.simplify(entry)) for entry in row] for row in matrix.tolist()]


def build():
    """Stable exact fixture data; no source-premise admission is inferred."""
    density = sp.ImmutableMatrix([[sp.Rational(2, 3), (1 + sp.I) / 6], [(1 - sp.I) / 6, sp.Rational(1, 3)]])
    assert positive_semidefinite(density) and sp.trace(density) == 1
    valuation = lambda matrix: sp.simplify(sp.trace(density * matrix))
    recovered = tomographic_candidate(valuation, 2)
    assert recovered == density
    first = effect(sp.diag(0, 1))
    delta = sp.Rational(1, 100)
    approximation = first + delta * sp.Matrix([[-sp.Rational(1, 2), sp.Rational(1, 3)], [sp.Rational(1, 3), sp.Rational(1, 2)]])
    assert order_sandwich(approximation - first, delta)
    regularized = regularize_hermitian(approximation, delta)
    approximation_bound = 2 * delta / (1 + 2 * delta)
    assert order_sandwich(regularized - first, approximation_bound)
    comparisons = lipschitz_order_witness(first, regularized, approximation_bound)
    joint = regularize_coexistent(approximation, sp.eye(2) - approximation, delta)
    # A positive normalized complement-respecting assignment can agree on all
    # rank-one qubit projectors and still fail a coexistent rational sum.
    cubic = lambda t: 3 * t - 6 * t**2 + 4 * t**3
    impostor = lambda matrix: cubic(sp.trace(matrix) / 2)
    false_candidate = tomographic_candidate(impostor, 2)
    witness = sp.eye(2) / 8
    discrepancy = sp.simplify(impostor(2 * witness) - 2 * impostor(witness))
    assert false_candidate == sp.eye(2) / 2 and discrepancy == -sp.Rational(9, 64)
    return {
        'schema': 'oph.dense_rational_effect_witnesses.v1',
        'theorem_scope': 'Countable rational-effect positivity and compatible additivity imply a unique Born extension; proof is mathematical, not inferred from these finite fixtures.',
        'existing_finite_theorem': 'Lean/EventAlgebra/FiniteBuschGleason.lean:finite_busch_gleason',
        'density_reconstruction': _entries(recovered),
        'regularization': {'delta': str(delta), 'effect': _entries(regularized), 'operator_error_bound': str(approximation_bound)},
        'lipschitz_comparisons': [[_entries(left), _entries(right)] for left, right in comparisons],
        'coexistent_regularization': [_entries(matrix) for matrix in joint],
        'finite_tomography_counterexample': {
            'valuation': 'f(Tr(E)/2), f(t)=3t-6t^2+4t^3',
            'candidate': _entries(false_candidate),
            'coexistent_pair': 'E=F=I/8',
            'value_of_sum': str(impostor(2 * witness)),
            'sum_of_values': str(2 * impostor(witness)),
            'additivity_defect': str(discrepancy)},
        'source_dense_effect_menu_realized': False,
        'source_equal_effect_record_equivalence_established': False,
        'canonical_physical_premises_discharged': [],
    }
