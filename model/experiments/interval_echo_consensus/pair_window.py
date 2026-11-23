"""Bound the four-appointment M10.3 mean over all initial two-load density matrices."""
from fractions import Fraction
from math import ceil

import numpy as np
from scipy import sparse

from .alias_window import setup
from .bounded_clock_flight import DetectorEvolution


TIMES = (0, 5, 10, 15)
SCALE = 1 << 20


def gram_ceiling(returns, column_error, trial):
    """Certify an upper bound for the mean effect by an integer Gram factorization."""
    duration_count, width, columns = returns.shape
    n = duration_count+1
    if (n != len(TIMES) or tuple(TIMES) != tuple(i*TIMES[1] for i in range(n))
            or TIMES[1] <= 0 or width != columns or not 1 <= width <= 30
            or not np.isfinite(returns).all() or np.max(abs(returns)) > 1.00001
            or not 0 <= column_error <= 1e-5 or not 0 < trial < 1):
        raise ValueError("Expected finite compressed unitary amplitudes and a column norm budget")
    # F stacks P U(t_j)/sqrt(n); the nonzero eigenvalues of F*F and FF* coincide.
    # The block (i,j) uses U(t_i-t_j). Centering H conjugates the Gram
    # by a diagonal unitary on its four time labels.
    dim = n*width
    gram = np.eye(dim, dtype=complex)
    for i in range(n):
        for j in range(i):
            a, b = slice(i*width, (i+1)*width), slice(j*width, (j+1)*width)
            delay = TIMES[i]-TIMES[j]
            gram[a,b] = returns[TIMES[1:].index(delay)]
            gram[b,a] = gram[a,b].conj().T
    real = np.rint(SCALE*gram.real).astype(np.int64)
    imag = np.rint(SCALE*gram.imag).astype(np.int64)
    block = np.block([[real, -imag], [imag, real]])
    # Dividing by n*SCALE gives the rounded real mean Gram. A floating factor
    # merely proposes a PSD witness. Its integer residual supplies the proof.
    scale = n*SCALE
    proposed = int(trial * scale)
    gap = proposed*np.eye(2*dim, dtype=np.int64)-block
    try:
        lower = np.rint(SCALE*np.linalg.cholesky(gap.astype(float)/scale))
    except np.linalg.LinAlgError as error:
        raise ArithmeticError("Trial does not enclose the mean Gram") from error
    if (not np.isfinite(lower).all() or np.max(abs(lower)) > 2*SCALE
            or (2*dim)*(2*SCALE)**2 + 4*SCALE*SCALE > np.iinfo(np.int64).max):
        raise ArithmeticError("Uncertified integer factor range")
    lower = lower.astype(np.int64)
    # n=4 divides SCALE. Check the absolute row sum with Python integers.
    residual = gap*(SCALE//n) - lower @ lower.T
    defect = max(sum(abs(int(x)) for x in row) for row in residual)
    # Per-column error and complex dyadic rounding each bound every entry.
    # The n-normalized off-diagonal block rows have total operator error
    # <= (n-1)*width/n times the maximum entry error.
    entry_error = Fraction(ceil(column_error*SCALE)+1, SCALE)
    return Fraction(proposed, scale) + Fraction(defect, SCALE*SCALE) + (n-1)*width*entry_error/n


def window_ceiling(strength):
    """Largest possible mean home-and-adjacent probability on the adopted board."""
    if strength not in (0, Fraction(1, 2)):
        raise ValueError("Expected exchange alone or the positive clock alias 1/2")
    exchange, adjacent, home, _ = setup()
    h = exchange - float(2*strength)*sparse.diags(adjacent, format="csr")
    flight = DetectorEvolution(h, TIMES[1:])
    returns = np.stack([flight.amplitudes(j, home) for j in home], axis=-1)
    trial = Fraction(933, 1000) if not strength else Fraction(9375, 10000)
    return gram_ceiling(returns, flight.amplitude_error, trial)


if __name__ == "__main__":
    for strength in (0, Fraction(1, 2)):
        bound = window_ceiling(strength)
        if bound >= Fraction(94, 100):
            raise ArithmeticError("The four-read mean does not exclude the continuous-window claim")
        print(f"kappa={strength}, four-appointment mean <= {ceil(bound*10**9)}/10^9")
