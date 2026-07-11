#!/usr/bin/env python3
"""Machine receipts for the null-net standardness packet (GitHub #524).

Implements finite-stage witnesses for the compact paper's subsection
`subsec:null-net-standardness`:

* stagewise standardness (Theorem `thm:null-net-standardness`, finite part):
  a faithful state on a finite-dimensional algebra has a cyclic and separating
  GNS vector, separation passes to subalgebras, and cyclicity of a subalgebra
  of the commutant implies separation (the mechanism used for half-lines);
* Markov modular locality (Theorem `thm:markov-modular-locality`): for exact
  quantum Markov states across a cut, the reduced modular Hamiltonian splits
  into one-sided terms plus a bounded interface block, giving the
  endpoint-Lipschitz control, with the interface bound independent of the
  interval length;
* the counterexample boundary
  (Proposition `prop:gibbs-modular-nonlocality-witness`): an explicit
  four-qubit nearest-neighbour Gibbs chain whose reduced modular Hamiltonian
  contains a nonzero endpoint-coupling term while the collar CMI is strictly
  positive, showing that finite-range Gibbs locality does NOT imply reduced
  modular locality without the Markov clause;
* the Lie-algebra consistency receipts of the four-translation assembly
  (Theorem `thm:four-translation-assembly`): in the defining affine
  representation of the Poincare group, per-null-direction translation
  generators commute, transform with the null-vector weight under Lorentz
  conjugation, satisfy the linearity relation over dependent null quadruples,
  and have future-cone dual positivity.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.linalg import expm, logm

TOL = 1e-9

I2 = np.eye(2)
SX = np.array([[0.0, 1.0], [1.0, 0.0]])
SY = np.array([[0.0, -1.0j], [1.0j, 0.0]])
SZ = np.array([[1.0, 0.0], [0.0, -1.0]])
PAULI = {"I": I2, "X": SX, "Y": SY, "Z": SZ}


def kron(*ops: np.ndarray) -> np.ndarray:
    out = np.array([[1.0 + 0.0j]])
    for op in ops:
        out = np.kron(out, op)
    return out


def site_op(op: np.ndarray, site: int, n: int) -> np.ndarray:
    return kron(*[op if i == site else I2 for i in range(n)])


def partial_trace(rho: np.ndarray, dims: list[int], keep: list[int]) -> np.ndarray:
    n = len(dims)
    keep = sorted(keep)
    rho = rho.reshape(dims + dims)
    traced = [i for i in range(n) if i not in keep]
    for cnt, i in enumerate(sorted(traced, reverse=True)):
        rho = np.trace(rho, axis1=i, axis2=i + n - cnt)
        n_now = n - cnt  # bookkeeping only
        del n_now
    d = int(np.prod([dims[i] for i in keep])) if keep else 1
    return rho.reshape(d, d)


def von_neumann_entropy(rho: np.ndarray) -> float:
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-14]
    return float(-np.sum(evals * np.log(evals)))


# ---------------------------------------------------------------------------
# stagewise standardness (finite part of thm:null-net-standardness)
# ---------------------------------------------------------------------------

def gns_vector(rho: np.ndarray) -> np.ndarray:
    """Omega = vec(sqrt(rho)) in H (x) H; the algebra acts as a (x) 1."""
    evals, vecs = np.linalg.eigh(rho)
    evals = np.clip(evals, 0.0, None)
    sqrt_rho = vecs @ np.diag(np.sqrt(evals)) @ vecs.conj().T
    return sqrt_rho.reshape(-1)


def algebra_basis(d: int) -> list[np.ndarray]:
    return [np.eye(d)[:, [i]] @ np.eye(d)[[j], :] for i in range(d) for j in range(d)]


def left_action(a: np.ndarray, d: int) -> np.ndarray:
    """a (x) 1 acting on vec(m) = (a m) as a matrix on C^{d^2}."""
    return np.kron(a, np.eye(d))


def right_action(a: np.ndarray, d: int) -> np.ndarray:
    """1 (x) a^T acting on vec(m) = (m a): commutant of the left action."""
    return np.kron(np.eye(d), a.T)


def cyclicity_dimension(ops: list[np.ndarray], omega: np.ndarray) -> int:
    mat = np.stack([op @ omega for op in ops], axis=1)
    return int(np.linalg.matrix_rank(mat, tol=1e-10))


def separating_modulus(ops_basis: list[np.ndarray], omega: np.ndarray) -> float:
    """min over unit-HS-norm algebra elements a of |a Omega|: > 0 iff separating.

    Assumes `ops_basis` is HS-orthogonal with equal norms (matrix units,
    possibly tensored with an identity), so a = sum_i c_i b_i / beta has unit
    HS norm exactly when |c| = 1 and the modulus is sigma_min(M) / beta for
    M = [b_1 Omega | b_2 Omega | ...]."""
    beta = float(np.linalg.norm(ops_basis[0], "fro"))
    mat = np.stack([op @ omega for op in ops_basis], axis=1)
    sv = np.linalg.svd(mat, compute_uv=False)
    smallest = sv[-1] if len(sv) == len(ops_basis) else 0.0
    if mat.shape[0] < len(ops_basis):
        smallest = 0.0  # more basis elements than dimensions: kernel exists
    return float(smallest) / beta


def full_algebra_standard(rho: np.ndarray) -> tuple[bool, bool]:
    """(cyclic, separating) for the full matrix algebra in its GNS space."""
    d = rho.shape[0]
    omega = gns_vector(rho)
    ops = [left_action(a, d) for a in algebra_basis(d)]
    cyclic = cyclicity_dimension(ops, omega) == d * d
    separating = separating_modulus(ops, omega) > 1e-10
    return cyclic, separating


def subalgebra_separating(rho: np.ndarray, dims: list[int], region: list[int]) -> bool:
    """Separation passes to subalgebras: the half-line algebra on `region`."""
    d = int(np.prod(dims))
    omega = gns_vector(rho)
    sub_ops = []
    d_reg = int(np.prod([dims[i] for i in region]))
    for a_small in algebra_basis(d_reg):
        # embed: operators supported on `region` (contiguous prefix assumed)
        rest = d // d_reg
        a_full = np.kron(a_small, np.eye(rest))
        sub_ops.append(left_action(a_full, d))
    return separating_modulus(sub_ops, omega) > 1e-10


def commutant_cyclicity_implies_separating(d: int, seed: int = 7) -> tuple[bool, bool]:
    """The separation mechanism of thm:null-net-standardness: on H_L (x) H_R
    with a full-Schmidt-rank vector, B = 1 (x) M_d lies in A' = (M_d (x) 1)',
    B Omega is dense, and A is separated by Omega."""
    rng = np.random.default_rng(seed)
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    # full Schmidt rank vector
    omega = (m @ m.conj().T + np.eye(d)).reshape(-1)
    omega = omega / np.linalg.norm(omega)
    ops_b = [right_action(a, d) for a in algebra_basis(d)]
    ops_a = [left_action(a, d) for a in algebra_basis(d)]
    b_cyclic = cyclicity_dimension(ops_b, omega) == d * d
    a_separated = separating_modulus(ops_a, omega) > 1e-10
    return b_cyclic, a_separated


# ---------------------------------------------------------------------------
# Markov modular locality (thm:markov-modular-locality)
# ---------------------------------------------------------------------------

def markov_chain_state(dims_a: int, dims_bl: int, dims_br: int, dims_c: int,
                       weights: list[float], seed: int = 11) -> np.ndarray:
    """Exact quantum Markov state rho = sum_k p_k rho^k_{A bL} (x) rho^k_{bR C}
    on A (x) (bL (x) bR) (x) C (single-block collar factorization per k, with
    the direct sum realized as classical mixing on a block label carried by
    the collar; here one block suffices for the locality receipt, several
    weights give the general mixed case on aligned supports)."""
    rng = np.random.default_rng(seed)

    def random_state(d: int) -> np.ndarray:
        m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        rho = m @ m.conj().T + 0.1 * np.eye(d)
        return rho / np.trace(rho)

    total = None
    for w in weights:
        left = random_state(dims_a * dims_bl)
        right = random_state(dims_br * dims_c)
        term = w * np.kron(left, right)
        total = term if total is None else total + term
    return total / np.trace(total)


def modular_hamiltonian(rho: np.ndarray) -> np.ndarray:
    return -logm(rho + 1e-300 * np.eye(rho.shape[0]))


def one_sided_split_defect(k_matrix: np.ndarray, d_left: int, d_right: int) -> float:
    """Distance of K to h_L (x) 1 + 1 (x) h_R (the Markov split of the theorem):
    zero iff K contains no cross terms between the two factors."""
    d = d_left * d_right
    k = k_matrix.reshape(d_left, d_right, d_left, d_right)
    h_l = np.trace(k, axis1=1, axis2=3) / d_right
    h_r = np.trace(k, axis1=0, axis2=2) / d_left
    trace_part = np.trace(k_matrix) / d
    split = (
        np.kron(h_l - np.trace(h_l) / d_left * np.eye(d_left), np.eye(d_right))
        + np.kron(np.eye(d_left), h_r - np.trace(h_r) / d_right * np.eye(d_right))
        + trace_part * np.eye(d)
    )
    return float(np.linalg.norm(k_matrix - split, 2))


def markov_modular_split_receipt(seed: int = 11, beta: float = 1.0,
                                 g: float = 1.0) -> dict[str, float]:
    """Receipts for Theorem thm:markov-modular-locality.

    * `markov_split_defect`: for an exact Markov (collar-factorized) state,
      the modular Hamiltonian splits exactly across the cut (defect ~ 0);
    * `gibbs_split_defect`: for the nearest-neighbour TFI Gibbs chain at the
      same sizes, the SAME functional is far from zero, so the split is a
      consequence of the Markov clause, not of finite-range locality;
    * `endpoint_lipschitz_growth`: extending the interval on the Markov
      branch changes K by exactly the added local block (interface-bounded,
      independent of the interval length)."""
    out = {}
    dims = (2, 2, 2, 2)  # A, bL, bR, C
    rho = markov_chain_state(*dims, weights=[1.0], seed=seed)
    d_a, d_bl, d_br, d_c = dims
    k_full = modular_hamiltonian(rho)
    out["markov_split_defect"] = one_sided_split_defect(k_full, d_a * d_bl, d_br * d_c)
    h = tfi_chain_hamiltonian(4, g)
    rho_gibbs = expm(-beta * h)
    rho_gibbs = rho_gibbs / np.trace(rho_gibbs)
    out["gibbs_split_defect"] = one_sided_split_defect(
        modular_hamiltonian(rho_gibbs), d_a * d_bl, d_br * d_c
    )
    # endpoint-Lipschitz on the Markov branch: interface bound stable in |A|
    rho_i = partial_trace(rho, list(dims), keep=[0, 1])
    k_i = modular_hamiltonian(rho_i)
    lips = []
    for extra in (1, 2):
        d_a2 = 2 ** extra
        rng = np.random.default_rng(seed + extra)
        m = rng.normal(size=(d_a2, d_a2))
        pad = m @ m.T + 0.1 * np.eye(d_a2)
        pad = pad / np.trace(pad)
        rho_ext = np.kron(pad, rho_i)
        k_ext = modular_hamiltonian(rho_ext)
        diff = k_ext - (
            np.kron(modular_hamiltonian(pad), np.eye(d_a * d_bl))
            + np.kron(np.eye(d_a2), k_i)
        )
        lips.append(float(np.linalg.norm(diff, 2)))
    out["endpoint_lipschitz_growth"] = max(lips)
    return out


# ---------------------------------------------------------------------------
# counterexample boundary (prop:gibbs-modular-nonlocality-witness)
# ---------------------------------------------------------------------------

def tfi_chain_hamiltonian(n: int, g: float = 1.0) -> np.ndarray:
    h = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(n - 1):
        h += site_op(SZ, i, n) @ site_op(SZ, i + 1, n)
    for i in range(n):
        h += g * site_op(SX, i, n)
    return h


def pauli_coefficients(k_matrix: np.ndarray, n: int) -> dict[str, complex]:
    coeffs = {}
    dim = 2 ** n
    for labels in itertools.product("IXYZ", repeat=n):
        op = kron(*[PAULI[c] for c in labels])
        coeffs["".join(labels)] = np.trace(op @ k_matrix) / dim
    return coeffs


def endpoint_coupling_norm(k_matrix: np.ndarray, n: int) -> float:
    """Sum of |coefficients| of Pauli words nontrivial on BOTH endpoint sites
    of the interval (sites 0 and n-1): the irreducibly nonlocal content."""
    coeffs = pauli_coefficients(k_matrix, n)
    total = 0.0
    for word, c in coeffs.items():
        if word[0] != "I" and word[-1] != "I":
            total += abs(c)
    return float(total)


def collar_cmi(rho: np.ndarray, dims: list[int], a: list[int], b: list[int],
               c: list[int]) -> float:
    """I(A:C|B) = S(AB) + S(BC) - S(B) - S(ABC)."""
    s_ab = von_neumann_entropy(partial_trace(rho, dims, a + b))
    s_bc = von_neumann_entropy(partial_trace(rho, dims, b + c))
    s_b = von_neumann_entropy(partial_trace(rho, dims, b))
    s_abc = von_neumann_entropy(partial_trace(rho, dims, a + b + c))
    return float(s_ab + s_bc - s_b - s_abc)


def gibbs_nonlocality_witness(beta: float = 1.0, g: float = 1.0) -> dict[str, float]:
    """The finite witness of prop:gibbs-modular-nonlocality-witness.

    Four-qubit nearest-neighbour transverse-field Ising Gibbs state; interval
    I = sites {0,1,2}. Reports the endpoint-coupling norm of K_I (nonzero:
    modular nonlocality), the collar CMI I(01 : 3 | 2) (nonzero: outside the
    Markov branch), and, for contrast, the same functionals for an exactly
    Markov product-comparison state (both zero)."""
    n = 4
    h = tfi_chain_hamiltonian(n, g)
    rho = expm(-beta * h)
    rho = rho / np.trace(rho)
    dims = [2] * n
    rho_i = partial_trace(rho, dims, keep=[0, 1, 2])
    k_i = modular_hamiltonian(rho_i)
    witness = {
        "endpoint_coupling_norm": endpoint_coupling_norm(k_i, 3),
        "collar_cmi": collar_cmi(rho, dims, a=[0, 1], b=[2], c=[3]),
    }
    # exactly Markov comparison: product of two-site Gibbs blocks
    h2 = tfi_chain_hamiltonian(2, g)
    blk = expm(-beta * h2)
    blk = blk / np.trace(blk)
    rho_markov = np.kron(blk, blk)
    rho_mi = partial_trace(rho_markov, dims, keep=[0, 1, 2])
    k_mi = modular_hamiltonian(rho_mi)
    witness["markov_endpoint_coupling_norm"] = endpoint_coupling_norm(k_mi, 3)
    witness["markov_collar_cmi"] = collar_cmi(rho_markov, dims, a=[0, 1], b=[2], c=[3])
    return witness


# ---------------------------------------------------------------------------
# four-translation assembly receipts (thm:four-translation-assembly)
# ---------------------------------------------------------------------------

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def null_vector(direction: np.ndarray) -> np.ndarray:
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    return np.concatenate(([1.0], d))


def translation_generator(q: np.ndarray) -> np.ndarray:
    """P(q) in the 5x5 affine representation of the Poincare group."""
    gen = np.zeros((5, 5))
    gen[:4, 4] = q
    return gen


def boost_generator(direction: np.ndarray) -> np.ndarray:
    """Lorentz boost generator along a spatial unit direction, 5x5 affine."""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    gen = np.zeros((5, 5))
    gen[0, 1:4] = d
    gen[1:4, 0] = d
    return gen


def lorentz_matrix(boost_dir: np.ndarray, rapidity: float) -> np.ndarray:
    return expm(rapidity * boost_generator(boost_dir))


def assembly_receipts(seed: int = 3) -> dict[str, float]:
    """Verify the algebraic clauses of Theorem thm:four-translation-assembly
    in the defining affine representation:

    * commutativity of null translations across directions;
    * Lorentz covariance with the null-vector weight
      Lambda P(q) Lambda^{-1} = P(Lambda q);
    * linearity over dependent null families (the MI-branch relation);
    * dual-cone positivity: eta(q, p) >= 0 for all future null q iff p is in
      the closed future cone (sampled check on both sides of the boundary).
    """
    rng = np.random.default_rng(seed)
    out = {}
    dirs = [rng.normal(size=3) for _ in range(6)]
    qs = [null_vector(d) for d in dirs]
    ps = [translation_generator(q) for q in qs]
    out["commutators"] = max(
        float(np.linalg.norm(p1 @ p2 - p2 @ p1)) for p1, p2 in itertools.combinations(ps, 2)
    )
    lam = lorentz_matrix(rng.normal(size=3), 0.7)
    lam4 = lam[:4, :4]
    out["covariance"] = max(
        float(np.linalg.norm(lam @ translation_generator(q) @ np.linalg.inv(lam)
                             - translation_generator(lam4 @ q)))
        for q in qs
    )
    # linearity: express q5 in the basis q1..q4 and compare generators
    basis = np.stack(qs[:4], axis=1)
    coeff = np.linalg.solve(basis, qs[4])
    assembled = sum(c * p for c, p in zip(coeff, ps[:4]))
    out["linearity"] = float(np.linalg.norm(assembled - ps[4]))
    # dual-cone positivity
    inside = np.array([2.0, 0.3, -0.4, 0.5])       # future timelike
    outside = np.array([0.5, 1.0, 0.0, 0.0])       # spacelike
    sample_qs = [null_vector(rng.normal(size=3)) for _ in range(200)]
    pos_inside = min(-float(q @ ETA @ inside) for q in sample_qs)
    pos_outside = min(-float(q @ ETA @ outside) for q in sample_qs)
    out["dual_cone_inside_min"] = pos_inside      # must be > 0
    out["dual_cone_outside_min"] = pos_outside    # must be < 0
    return out
