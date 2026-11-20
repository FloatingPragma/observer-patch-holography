"""The local contrast-memory law on one seam and its exact identities.

One seam joins two ports `a` and `b`. Every port holds a load in {0, 1, 2}, so the seam sees nine load states
`|n_a, n_b>`. The seam carries a private memory of eleven states: one blank and ten labels `(c, d)` with coin
`c` in {0, 1} and contrast `d` in {-2, ..., 2}. The local space has 9 x 11 = 99 states, ordered as the nine
load states with blank memory first, then label index times nine plus load index.

- Native event on loads `(a, b)`: `s = a + b`, `p = s mod 2`, output for coin `c`
  `y_c = (floor(s/2) + c p, floor(s/2) + (1 - c) p)`.
- Contrast write, an isometry: `W|a, b> = 2^(-1/2) sum_c |c, a - b>|y_c>`.
- Echo: `U = [[0, W*], [W, I - W W*]]`, a Hermitian unitary on all 99 states.
- Record operators: `R|c, d> = |c, -d>`, `D|c, d> = d|c, d>`.
- Memory generator `G_kappa = 0_blank (+) [(I - R)/2 + kappa D^2] (x) I` and local term `h_kappa = U G_kappa U*`.

On the blank memory the local term reduces to the seam term on loads,
`h_kappa J = J [(I - SWAP)/2 + kappa (n_a - n_b)^2]`, for every real `kappa`. `exact_checks` verifies this and
the other identities with exact SymPy arithmetic.
"""
from __future__ import annotations

import numpy as np

SRC = tuple((a, b) for a in range(3) for b in range(3))
LABELS = tuple((c, d) for c in (0, 1) for d in range(-2, 3))
N_SRC, N_LAB = len(SRC), len(LABELS)
DIM = N_SRC * (1 + N_LAB)


def native_event(a: int, b: int, c: int) -> tuple[int, int]:
    """Output loads of the native seam event for coin `c`."""
    s = a + b
    p = s % 2
    return (s // 2 + c * p, s // 2 + (1 - c) * p)


def _write_pattern():
    """Integer matrix `sqrt(2) W` (90 x 9): record index first, loads second."""
    Wp = np.zeros((N_LAB * N_SRC, N_SRC), dtype=np.int64)
    for col, (a, b) in enumerate(SRC):
        for c in (0, 1):
            y = native_event(a, b, c)
            Wp[LABELS.index((c, a - b)) * N_SRC + SRC.index(y), col] += 1
    return Wp


def _label_ops():
    R = np.zeros((N_LAB, N_LAB), dtype=np.int64)
    for i, (c, d) in enumerate(LABELS):
        R[LABELS.index((c, -d)), i] = 1
    D2 = np.diag([d * d for c, d in LABELS]).astype(np.int64)
    return R, D2


def _load_ops():
    T = np.zeros((N_SRC, N_SRC), dtype=np.int64)
    for i, (a, b) in enumerate(SRC):
        T[SRC.index((b, a)), i] = 1
    C = np.diag([(a - b) ** 2 for a, b in SRC]).astype(np.int64)
    return T, C


def write() -> np.ndarray:
    """The contrast write `W` as a float matrix (90 x 9)."""
    return _write_pattern() / np.sqrt(2.0)


def echo() -> np.ndarray:
    """The echo unitary `U` on the 99 local states."""
    W = write()
    U = np.zeros((DIM, DIM))
    U[:N_SRC, N_SRC:] = W.T
    U[N_SRC:, :N_SRC] = W
    U[N_SRC:, N_SRC:] = np.eye(N_LAB * N_SRC) - W @ W.T
    return U


def generator_parts():
    """`(G_exchange, G_contrast)` with `G_kappa = G_exchange + kappa G_contrast` on the 99 local states."""
    R, D2 = _label_ops()
    G0 = np.zeros((DIM, DIM))
    G1 = np.zeros((DIM, DIM))
    G0[N_SRC:, N_SRC:] = np.kron((np.eye(N_LAB) - R) / 2, np.eye(N_SRC))
    G1[N_SRC:, N_SRC:] = np.kron(D2, np.eye(N_SRC))
    return G0, G1


def local_term_parts():
    """`(h_exchange, h_contrast)` with `h_kappa = h_exchange + kappa h_contrast`."""
    U = echo()
    G0, G1 = generator_parts()
    return U @ G0 @ U.T, U @ G1 @ U.T


def seam_term(kappa: float) -> np.ndarray:
    """The reduced seam term `(I - SWAP)/2 + kappa (n_a - n_b)^2` on the nine load states."""
    T, C = _load_ops()
    return (np.eye(N_SRC) - T) / 2 + kappa * C


def exact_checks() -> dict:
    """Exact SymPy checks of the local law, for symbolic `kappa`.

    Matrices are built over the field Q(sqrt 2) with the symbol `kappa`; every identity is checked entry by
    entry after expansion, so it holds as a polynomial identity in `kappa`, hence for every real `kappa`.
    """
    import sympy as sp

    k = sp.Symbol('kappa', real=True)
    s2 = sp.sqrt(2)
    Wp = _write_pattern()
    R_, D2_ = _label_ops()
    T_, C_ = _load_ops()
    W = sp.Matrix(Wp.tolist()) * s2 / 2
    R = sp.Matrix(R_.tolist())
    D2 = sp.Matrix(D2_.tolist())
    T = sp.Matrix(T_.tolist())
    C = sp.Matrix(C_.tolist())
    I9, I10, I90 = sp.eye(N_SRC), sp.eye(N_LAB), sp.eye(N_LAB * N_SRC)
    U = sp.BlockMatrix([[sp.zeros(N_SRC), W.T], [W, I90 - W * W.T]]).as_explicit()
    G = sp.diag(sp.zeros(N_SRC), sp.kronecker_product((I10 - R) / 2 + k * D2, I9))
    J = sp.zeros(DIM, N_SRC)
    J[:N_SRC, :] = I9
    reduced = (I9 - T) / 2 + k * C

    def zero(M):
        return all(sp.expand(x) == 0 for x in M)

    hJ = U * (G * (U * J))                      # h_kappa J with U = U*
    number = sp.diag(*[a + b for a, b in SRC])
    squares = sp.diag(*[a * a + b * b for a, b in SRC])
    checks = {
        'native_event_conserves_load': all(sum(native_event(a, b, c)) == a + b for a, b in SRC for c in (0, 1)),
        'W_isometry': zero(W.T * W - I9),
        'U_hermitian': zero(U - U.T),
        'U_unitary_on_all_99_states': zero(U * U - sp.eye(DIM)),
        'R_intertwines_swap': zero(sp.kronecker_product(R, I9) * W - W * T),
        'D2_intertwines_contrast_square': zero(sp.kronecker_product(D2, I9) * W - W * C),
        'R_D2_commute': zero(R * D2 - D2 * R),
        'h_kappa_J_equals_J_seam_term': zero(hJ - J * reduced),
        'blank_code_reducing': zero(hJ[N_SRC:, :]),
        'seam_term_hermitian': zero(reduced - reduced.T),
        'seam_term_conserves_load_number': zero(reduced * number - number * reduced),
        'seam_term_conserves_sum_of_squared_loads': zero(reduced * squares - squares * reduced),
    }
    # binary loads: the seam term is the ferromagnetic XXZ interaction
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Z = sp.Matrix([[1, 0], [0, -1]])
    kp = sp.kronecker_product
    xxz = ((1 + 2 * k) * (sp.eye(4) - kp(Z, Z)) - kp(X, X) - kp(Y, Y)) / 4
    binary = [SRC.index((a, b)) for a in (0, 1) for b in (0, 1)]
    checks['binary_seam_term_is_ferromagnetic_XXZ'] = zero(reduced.extract(binary, binary) - xxz)
    # degree term: D^2 alone on a load/reference pair (n_v, n_ref)
    G_ref = sp.diag(sp.zeros(N_SRC), sp.kronecker_product(D2, I9))
    hJ_ref = U * (G_ref * (U * J))
    checks['reference_pair_D2_reduces_to_contrast_square'] = zero(hJ_ref - J * C)
    zero_ref = [SRC.index((n, 0)) for n in range(3)]
    checks['reference_pair_gives_n_squared_at_reference_zero'] = all(
        C[i, i] == SRC[i][0] ** 2 for i in zero_ref)
    checks['reference_stays_at_zero'] = zero(C.extract([i for i in range(N_SRC) if i not in zero_ref], zero_ref))
    exchange_on_reference = (I9 - T) / 2
    checks['exchange_on_reference_pair_would_move_load'] = not zero(
        exchange_on_reference.extract([i for i in range(N_SRC) if i not in zero_ref], zero_ref))
    if not all(checks.values()):
        raise ValueError({key: v for key, v in checks.items() if not v})
    return {'local_dimension': DIM, 'load_states': N_SRC, 'memory_states': 1 + N_LAB,
            'checks': checks,
            'element_02_20_of_exchange': str(((I9 - T) / 2)[SRC.index((0, 2)), SRC.index((2, 0))])}


def numeric_checks(tol: float = 1e-12) -> dict:
    """The same identities in floating point (used by the tests for fast mutation checks)."""
    U = echo()
    h0, h1 = local_term_parts()
    J = np.zeros((DIM, N_SRC))
    J[:N_SRC] = np.eye(N_SRC)
    out = {}
    for kappa in (0.0, 0.5, 8.0, -3.0):
        h = h0 + kappa * h1
        out[kappa] = float(np.abs(h @ J - J @ seam_term(kappa)).max())
    return {'U_unitary': float(np.abs(U @ U - np.eye(DIM)).max()) < tol,
            'reduction_residual': max(out.values())}
