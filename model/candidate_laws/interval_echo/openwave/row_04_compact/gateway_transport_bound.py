"""Uniform home-detection bound for the registered Poisson swap control.

With rate 1/2 on each seam, the continuous generator is -L/2 = 3(T-I),
where T = I-L/6 is stochastic. For a home indicator h, max(T**k h)
decreases with k. After 24 discrete steps the integer walk count bounds
every home simultaneously. At exposure 12 the step count is Poisson(36).
The early-count probability at later exposures can only decrease. A convex
mixture of starting ports satisfies the same bound, including a release
whose total mass is less than one because preparation can exhaust.
"""

from fractions import Fraction
from math import ceil, factorial

import numpy as np
from scipy import sparse

from model.executable_model import federation


STEPS = 24
EXPOSURE = 12


def exchange_step_counts(board):
    """Return 6I-L with nonnegative integer entries and row and column sums six."""
    lap = board.laplacian()
    integer_lap = lap.astype(np.int64)
    if (lap != integer_lap).nnz or board.n != 12 * board.cells:
        raise ArithmeticError("Expected an integer port Laplacian and complete home cells")
    step = 6 * sparse.eye(board.n, format="csr", dtype=np.int64) - integer_lap
    if ((step != step.T).nnz or np.any(step.data < 0)
            or np.any(np.asarray(step.sum(axis=1)) != 6)):
        raise ArithmeticError("Uniformized exchange must be symmetric and stochastic")
    return step


def _walk_ceiling(board):
    """Maximum home probability after STEPS, counted over all sharp starts and homes."""
    step = exchange_step_counts(board)
    denominator = 6 ** STEPS
    if denominator > np.iinfo(np.int64).max:
        raise ArithmeticError("Walk counts exceed the integer accumulator")
    # Home indicators partition the ports. Each nonnegative row has sum 6**k
    # at step k, which bounds every intermediate sum in the multiplication.
    counts = np.repeat(np.eye(board.cells, dtype=np.int64), 12, axis=0)
    for _ in range(STEPS):
        counts = step @ counts
    if np.any(counts.sum(axis=1) != denominator):
        raise ArithmeticError("Walk count lost a port or a home")
    return Fraction(int(counts.max()), denominator)


def swap_home_ceiling():
    """Return the discrete ceiling, early-count ceiling and full-flight ceiling."""
    walk = _walk_ceiling(federation(0))
    mean = 3 * EXPOSURE
    # Truncating exp(mean) below bounds P(N < STEPS) above, using integers only.
    weights = [Fraction(mean**k, factorial(k)) for k in range(49)]
    early = sum(weights[:STEPS]) / sum(weights)
    return walk, early, walk + (1 - walk) * early


if __name__ == "__main__":
    walk, early, ceiling = swap_home_ceiling()
    if ceiling >= Fraction(4, 5):
        raise ArithmeticError("Poisson swap law exceeds the registered home ceiling")
    print(f"20 cells, 240 starts, all exposures >= {EXPOSURE}")
    print(f"Discrete walk ceiling {walk}, early-count ceiling {early}")
    print(f"Joint home probability <= {ceil(10**12 * ceiling)}/10^12")
