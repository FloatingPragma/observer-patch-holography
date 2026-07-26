#!/usr/bin/env python3
"""Numerical verification of the Markov-split alignment (MSA) characterizations.

Implements the objects of The spacetime and Einstein paper's Proposition `prop:msachar`
(operational characterizations of Markov-split alignment) and Corollary
`cor:msareduction` (axiom-side reduction), introduced for paper-audit issue 001
(GitHub #543).

Collar model: H = oplus_alpha  H_A (x) H_{bL^alpha} (x) H_{bR^alpha} (x) H_D.
A state is represented blockwise as a list of (weight, block density matrix,
(dA, dL, dR, dD)) triples, which enforces [rho, P_alpha] = 0 by construction;
`embed_blocks` produces the full direct-sum matrix when needed.

Checked characterizations, for faithful blockwise states:

1. EC-aligned normal form  rho = oplus p_a rho_{A bL}^a (x) rho_{bR D}^a
2. modular splitting       log rho in M_L + M_R (blockwise K = h_L (x) 1 + 1 (x) h_R)
3. Takesaki criterion      [log rho, x (x) 1] in B(H_{A bL}) (x) 1  for all x
4. entropic criterion      I(A bL : bR D) = 0 on every block

together with the implication (any of 1-4)  =>  I(A:D|B) = 0, and the
Bell-pair counterexample showing that I(A:D|B) = 0 alone implies none of 1-4.
"""

from __future__ import annotations

import numpy as np

LN2 = float(np.log(2.0))


# ---------------------------------------------------------------------------
# linear-algebra helpers
# ---------------------------------------------------------------------------

def dagger(m: np.ndarray) -> np.ndarray:
    return m.conj().T


def partial_trace(rho: np.ndarray, dims: list[int], keep: list[int]) -> np.ndarray:
    """Trace out all subsystems not in `keep` (indices into `dims`)."""
    n = len(dims)
    keep = sorted(keep)
    rho = rho.reshape(dims + dims)
    traced = [i for i in range(n) if i not in keep]
    for count, t in enumerate(sorted(traced, reverse=True)):
        cur_n = n - count
        rho = np.trace(rho, axis1=t, axis2=t + cur_n)
    d = int(np.prod([dims[k] for k in keep])) if keep else 1
    return rho.reshape(d, d)


def von_neumann_entropy(rho: np.ndarray) -> float:
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-14]
    return float(-np.sum(evals * np.log(evals)))


def mutual_information(rho: np.ndarray, dims: list[int],
                       part_x: list[int], part_y: list[int]) -> float:
    s_x = von_neumann_entropy(partial_trace(rho, dims, part_x))
    s_y = von_neumann_entropy(partial_trace(rho, dims, part_y))
    s_xy = von_neumann_entropy(partial_trace(rho, dims, sorted(part_x + part_y)))
    return s_x + s_y - s_xy


def conditional_mutual_information(rho: np.ndarray, dims: list[int],
                                   part_a: list[int], part_b: list[int],
                                   part_d: list[int]) -> float:
    """I(A:D|B) = S(AB) + S(BD) - S(B) - S(ABD)."""
    s_ab = von_neumann_entropy(partial_trace(rho, dims, sorted(part_a + part_b)))
    s_bd = von_neumann_entropy(partial_trace(rho, dims, sorted(part_b + part_d)))
    s_b = von_neumann_entropy(partial_trace(rho, dims, part_b))
    s_abd = von_neumann_entropy(
        partial_trace(rho, dims, sorted(part_a + part_b + part_d)))
    return s_ab + s_bd - s_b - s_abd


def logm_psd(rho: np.ndarray) -> np.ndarray:
    """Matrix logarithm of a faithful density matrix."""
    evals, vecs = np.linalg.eigh(rho)
    if np.min(evals) <= 0.0:
        raise ValueError("state is not faithful; log rho undefined")
    return (vecs * np.log(evals)) @ dagger(vecs)


# ---------------------------------------------------------------------------
# blockwise state representation
# ---------------------------------------------------------------------------

Block = tuple[float, np.ndarray, tuple[int, int, int, int]]


def block_dims(block: Block) -> list[int]:
    return list(block[2])


def embed_blocks(blocks: list[Block]) -> tuple[np.ndarray, list[np.ndarray]]:
    """Direct-sum embedding; returns (rho_full, center projectors P_alpha)."""
    sizes = [b[1].shape[0] for b in blocks]
    total = sum(sizes)
    rho = np.zeros((total, total), dtype=complex)
    projectors = []
    offset = 0
    for (p, rho_a, _), size in zip(blocks, sizes):
        rho[offset:offset + size, offset:offset + size] = p * rho_a
        proj = np.zeros((total, total))
        proj[offset:offset + size, offset:offset + size] = np.eye(size)
        projectors.append(proj)
        offset += size
    return rho, projectors


def collar_cmi(blocks: list[Block]) -> float:
    """I(A:D|B) of the blockwise state; entropies add over blocks."""
    total = 0.0
    for p, rho_a, dims in blocks:
        if p <= 0.0:
            continue
        d = list(dims)
        total += p * conditional_mutual_information(rho_a, d, [0], [1, 2], [3])
    # Classical block label contributes equally to S(AB), S(BD), S(B), S(ABD)
    # and cancels in the CMI combination, so the blockwise sum is exact.
    return total


# ---------------------------------------------------------------------------
# characterization 4: entropic criterion (no faithfulness needed)
# ---------------------------------------------------------------------------

def entropic_alignment_defect(blocks: list[Block]) -> float:
    """max_alpha I(A bL^alpha : bR^alpha D); zero iff EC-aligned (item 4)."""
    worst = 0.0
    for p, rho_a, dims in blocks:
        if p <= 0.0:
            continue
        d = list(dims)
        worst = max(worst, mutual_information(rho_a, d, [0, 1], [2, 3]))
    return worst


def is_ec_aligned(blocks: list[Block], tol: float = 1e-9) -> bool:
    return entropic_alignment_defect(blocks) < tol


# ---------------------------------------------------------------------------
# characterization 2: modular splitting  log rho in M_L + M_R
# ---------------------------------------------------------------------------

def one_sided_projection(k: np.ndarray, d_left: int, d_right: int) -> np.ndarray:
    """Hilbert-Schmidt projection of K onto B(H_L) (x) 1 + 1 (x) B(H_R)."""
    k4 = k.reshape(d_left, d_right, d_left, d_right)
    k_l = np.trace(k4, axis1=1, axis2=3) / d_right       # Tr_R K / d_R
    k_r = np.trace(k4, axis1=0, axis2=2) / d_left        # Tr_L K / d_L
    trace_part = np.trace(k) / (d_left * d_right)
    return (np.kron(k_l, np.eye(d_right)) + np.kron(np.eye(d_left), k_r)
            - trace_part * np.eye(d_left * d_right))


def modular_splitting_defect(blocks: list[Block]) -> float:
    """max_alpha || log rho^alpha - Proj_{M_L + M_R}(log rho^alpha) ||_inf.

    Zero iff log rho in M_L + M_R blockwise (item 2). Requires faithful blocks.
    """
    worst = 0.0
    for p, rho_a, dims in blocks:
        if p <= 0.0:
            continue
        d_a, d_l, d_r, d_d = dims
        k = logm_psd(rho_a)
        proj = one_sided_projection(k, d_a * d_l, d_r * d_d)
        worst = max(worst, float(np.linalg.norm(k - proj, ord=2)))
    return worst


# ---------------------------------------------------------------------------
# characterization 3: Takesaki commutator criterion
# ---------------------------------------------------------------------------

def takesaki_defect(blocks: list[Block], rng: np.random.Generator | None = None,
                    n_probes: int = 8) -> float:
    """max over blocks/probes of the off-subalgebra part of [log rho, x (x) 1].

    Zero iff Ad rho^{it} preserves B(H_{A bL}) (x) 1 for all t, i.e. iff a
    rho-preserving conditional expectation onto M_L exists (item 3, via the
    derivative form of Takesaki's criterion used in the paper proof).
    """
    if rng is None:
        rng = np.random.default_rng(7)
    worst = 0.0
    for p, rho_a, dims in blocks:
        if p <= 0.0:
            continue
        d_a, d_l, d_r, d_d = dims
        dl_tot, dr_tot = d_a * d_l, d_r * d_d
        k = logm_psd(rho_a)
        for _ in range(n_probes):
            x = rng.normal(size=(dl_tot, dl_tot)) + 1j * rng.normal(size=(dl_tot, dl_tot))
            x = (x + dagger(x)) / 2.0
            x_full = np.kron(x, np.eye(dr_tot))
            comm = k @ x_full - x_full @ k
            # component of comm orthogonal to B(H_L) (x) 1
            comm4 = comm.reshape(dl_tot, dr_tot, dl_tot, dr_tot)
            comm_l = np.trace(comm4, axis1=1, axis2=3) / dr_tot
            resid = comm - np.kron(comm_l, np.eye(dr_tot))
            worst = max(worst, float(np.linalg.norm(resid, ord=2)))
    return worst


# ---------------------------------------------------------------------------
# state constructors
# ---------------------------------------------------------------------------

def random_density(dim: int, rng: np.random.Generator, floor: float = 0.05) -> np.ndarray:
    """Random faithful density matrix with spectral floor."""
    g = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    rho = g @ dagger(g)
    rho = rho / np.trace(rho).real
    rho = (1.0 - floor) * rho + floor * np.eye(dim) / dim
    return (rho + dagger(rho)) / 2.0


def random_aligned_blocks(rng: np.random.Generator,
                          spec: list[tuple[int, int, int, int]] | None = None
                          ) -> list[Block]:
    """Random EC-aligned exact Markov state: blockwise product across the cut."""
    if spec is None:
        spec = [(2, 2, 2, 2), (2, 2, 2, 2)]
    weights = rng.dirichlet(np.ones(len(spec)))
    blocks: list[Block] = []
    for w, dims in zip(weights, spec):
        d_a, d_l, d_r, d_d = dims
        left = random_density(d_a * d_l, rng)
        right = random_density(d_r * d_d, rng)
        blocks.append((float(w), np.kron(left, right), dims))
    return blocks


def random_generic_blocks(rng: np.random.Generator,
                          spec: list[tuple[int, int, int, int]] | None = None
                          ) -> list[Block]:
    """Random faithful blockwise state with generic cross-cut correlation."""
    if spec is None:
        spec = [(2, 2, 2, 2)]
    weights = rng.dirichlet(np.ones(len(spec)))
    return [(float(w), random_density(int(np.prod(d)), rng), d)
            for w, d in zip(weights, spec)]


def bell_state(dim: int = 2) -> np.ndarray:
    vec = np.eye(dim).reshape(-1) / np.sqrt(dim)
    return np.outer(vec, vec.conj())


def bell_counterexample() -> list[Block]:
    """rho = Phi(A, bR) (x) Phi(bL, D): exactly Markov, maximally misaligned.

    Single center block, qubits ordered (A, bL, bR, D).
    """
    phi = bell_state(2)
    rho = np.kron(phi, phi)  # order (A, bR, bL, D)
    rho = rho.reshape([2] * 8)
    perm = [0, 2, 1, 3]      # -> (A, bL, bR, D)
    rho = rho.transpose(perm + [p + 4 for p in perm]).reshape(16, 16)
    return [(1.0, rho, (2, 2, 2, 2))]


def gibbs_blocks(hamiltonians: list[tuple[np.ndarray, tuple[int, int, int, int]]],
                 central_energies: list[float], beta: float = 1.0) -> list[Block]:
    """Gibbs state of H = oplus_alpha (H^alpha + e_alpha P_alpha) at inverse temp beta."""
    unnorm = []
    for (h_a, dims), e_a in zip(hamiltonians, central_energies):
        evals, vecs = np.linalg.eigh(h_a)
        block = (vecs * np.exp(-beta * (evals + e_a))) @ dagger(vecs)
        unnorm.append((block, dims))
    z_total = sum(np.trace(b).real for b, _ in unnorm)
    return [(float(np.trace(b).real / z_total), b / np.trace(b).real, dims)
            for b, dims in unnorm]


def one_sided_hamiltonian(rng: np.random.Generator,
                          dims: tuple[int, int, int, int]) -> np.ndarray:
    """H = h_L (x) 1 + 1 (x) h_R: every cross-cut coupling routed through the center."""
    d_a, d_l, d_r, d_d = dims
    h_l = rng.normal(size=(d_a * d_l, d_a * d_l))
    h_l = (h_l + h_l.T) / 2.0
    h_r = rng.normal(size=(d_r * d_d, d_r * d_d))
    h_r = (h_r + h_r.T) / 2.0
    return np.kron(h_l, np.eye(d_r * d_d)) + np.kron(np.eye(d_a * d_l), h_r)


def cross_cut_coupling(rng: np.random.Generator,
                       dims: tuple[int, int, int, int]) -> np.ndarray:
    """A traceless non-central bL--bR interaction term crossing the cut."""
    d_a, d_l, d_r, d_d = dims
    g_l = rng.normal(size=(d_l, d_l))
    g_l = (g_l + g_l.T) / 2.0
    g_l -= np.trace(g_l) / d_l * np.eye(d_l)
    g_r = rng.normal(size=(d_r, d_r))
    g_r = (g_r + g_r.T) / 2.0
    g_r -= np.trace(g_r) / d_r * np.eye(d_r)
    return np.kron(np.kron(np.eye(d_a), g_l), np.kron(g_r, np.eye(d_d)))
