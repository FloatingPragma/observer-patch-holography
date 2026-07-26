#!/usr/bin/env python3
"""Tests for the MSA characterization package (paper-audit issue 001, GitHub #543).

Each test verifies one clause of The spacetime and Einstein paper's Proposition `prop:msachar`
or Corollary `cor:msareduction`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from msa_characterizations import (  # noqa: E402
    LN2,
    bell_counterexample,
    collar_cmi,
    cross_cut_coupling,
    embed_blocks,
    entropic_alignment_defect,
    gibbs_blocks,
    is_ec_aligned,
    modular_splitting_defect,
    mutual_information,
    one_sided_hamiltonian,
    partial_trace,
    random_aligned_blocks,
    random_generic_blocks,
    takesaki_defect,
)

TOL = 1e-9
RNG = np.random.default_rng(20260711)


def test_bell_counterexample_is_markov_but_not_aligned():
    """The audit's Bell-pair state: I(A:D|B)=0 yet every alignment criterion fails."""
    blocks = bell_counterexample()
    assert collar_cmi(blocks) < TOL
    # item 4 fails: one-sided mutual information is large
    defect = entropic_alignment_defect(blocks)
    assert defect > 2.0 * LN2 - 1e-6
    assert not is_ec_aligned(blocks)
    # the specific witness used in the papers: I(A:B_R) = 2 ln 2
    _, rho, dims = blocks[0]
    i_a_br = mutual_information(rho, list(dims), [0], [2])
    assert abs(i_a_br - 2.0 * LN2) < 1e-9


def test_aligned_states_satisfy_all_characterizations_and_are_markov():
    """Items 1 => 2, 3, 4 and => I(A:D|B)=0 on random EC-aligned states."""
    for seed in range(5):
        rng = np.random.default_rng(seed)
        blocks = random_aligned_blocks(rng)
        assert is_ec_aligned(blocks), "item 4 must hold for aligned states"
        assert modular_splitting_defect(blocks) < 1e-8, "item 2 must hold"
        assert takesaki_defect(blocks, rng) < 1e-8, "item 3 must hold"
        assert collar_cmi(blocks) < 1e-8, "alignment must imply exact Markovianity"


def test_generic_states_fail_all_characterizations_consistently():
    """Items 2, 3, 4 fail together on generic (non-aligned) faithful states."""
    for seed in range(5):
        rng = np.random.default_rng(100 + seed)
        blocks = random_generic_blocks(rng)
        e_defect = entropic_alignment_defect(blocks)
        m_defect = modular_splitting_defect(blocks)
        t_defect = takesaki_defect(blocks, rng)
        assert e_defect > 1e-3
        assert m_defect > 1e-3
        assert t_defect > 1e-3


def test_characterizations_agree_on_mixed_ensembles():
    """Items 2 and 4 vanish together (equivalence), including multi-block states."""
    rng = np.random.default_rng(42)
    for _ in range(4):
        aligned = random_aligned_blocks(
            rng, spec=[(2, 2, 2, 2), (2, 2, 2, 2), (2, 2, 2, 2)])
        assert entropic_alignment_defect(aligned) < 1e-8
        assert modular_splitting_defect(aligned) < 1e-8
        generic = random_generic_blocks(rng, spec=[(2, 2, 2, 2), (2, 2, 2, 2)])
        assert entropic_alignment_defect(generic) > 1e-4
        assert modular_splitting_defect(generic) > 1e-4


def test_block_diagonal_embedding_is_consistent():
    """The direct-sum embedding commutes with the center projectors."""
    rng = np.random.default_rng(3)
    blocks = random_aligned_blocks(rng)
    rho, projectors = embed_blocks(blocks)
    assert abs(np.trace(rho).real - 1.0) < 1e-12
    for proj in projectors:
        assert np.linalg.norm(rho @ proj - proj @ rho) < 1e-12


def test_gibbs_with_central_interface_is_aligned():
    """Corollary cor:msareduction, positive direction.

    H_eff = H_L + H_R + H_Z with the cross-cut energy central (a function of
    the block label alpha) produces an EC-aligned exact Markov Gibbs state.
    """
    rng = np.random.default_rng(9)
    dims = (2, 2, 2, 2)
    hams = [(one_sided_hamiltonian(rng, dims), dims) for _ in range(2)]
    central_energies = [0.0, 0.7]  # H_Z = f(alpha): central interface energy
    blocks = gibbs_blocks(hams, central_energies)
    assert is_ec_aligned(blocks, tol=1e-8)
    assert modular_splitting_defect(blocks) < 1e-8
    assert takesaki_defect(blocks, rng) < 1e-8
    assert collar_cmi(blocks) < 1e-8


def test_gibbs_with_noncentral_cross_cut_coupling_breaks_alignment():
    """Corollary cor:msareduction, negative direction.

    Adding a non-central bL--bR coupling to H_eff destroys every alignment
    criterion (and, generically, exact Markovianity).
    """
    rng = np.random.default_rng(11)
    dims = (2, 2, 2, 2)
    base = one_sided_hamiltonian(rng, dims)
    coupled = base + 1.5 * cross_cut_coupling(rng, dims)
    blocks = gibbs_blocks([(coupled, dims)], [0.0])
    assert entropic_alignment_defect(blocks) > 1e-3
    assert modular_splitting_defect(blocks) > 1e-3
    assert takesaki_defect(blocks, rng) > 1e-3
    assert collar_cmi(blocks) > 1e-6


def test_no_modulus_controls_alignment_distance():
    """Remark rem:hjpwalignment: CMI -> 0 does not force alignment defect -> 0.

    Interpolating toward the Bell counterexample keeps CMI small while the
    alignment defect stays bounded away from zero, so no vanishing modulus
    delta(eps) for the EC-aligned class exists.
    """
    bell = bell_counterexample()[0][1]
    dims = (2, 2, 2, 2)
    ident = np.eye(16) / 16.0
    for eps in (1e-3, 1e-4, 1e-5):
        rho = (1.0 - eps) * bell + eps * ident
        blocks = [(1.0, rho, dims)]
        cmi = collar_cmi(blocks)
        defect = entropic_alignment_defect(blocks)
        assert cmi < 0.05, "state stays near-Markov"
        assert defect > 1.0, "alignment defect stays O(2 ln 2)"


def test_bell_smoothing_is_excluded_from_central_interface_gibbs_class():
    """Theorem thm:msaderivation, exclusion clause.

    Every Gibbs state of a central-interface package is EC-aligned (and
    faithful). Faithful smoothings of the Bell counterexample keep a
    cross-cut component of log rho outside M_L + M_R and keep
    I(A:B_R) > 0, so neither the counterexample nor its faithful
    neighborhood arises on the declared central-interface branch.
    """
    bell = bell_counterexample()[0][1]
    dims = (2, 2, 2, 2)
    ident = np.eye(16) / 16.0
    rng = np.random.default_rng(17)
    for eps in (1e-2, 1e-3):
        rho = (1.0 - eps) * bell + eps * ident
        blocks = [(1.0, rho, dims)]
        # faithful, so the modular criteria are defined -- and they fail:
        assert modular_splitting_defect(blocks) > 1e-2, \
            "log rho must keep a cross-cut term outside M_L + M_R"
        assert takesaki_defect(blocks, rng) > 1e-2, \
            "no rho-preserving conditional expectation onto M_L exists"
        # the entropic witness of the exclusion:
        i_a_br = mutual_information(rho, list(dims), [0], [2])
        assert i_a_br > LN2, "I(A:B_R) stays bounded away from zero"
        assert not is_ec_aligned(blocks)


def _pre_quotient_model(rng: np.random.Generator):
    """Z_2 boundary-action collar model for Theorem thm:msaderivation.

    H~_BL = C^2(x)C^2 with u_L = diag(1,-1)(x)I_2 (irrep label (x) multiplicity),
    H~_BR carries the contragredient (= same, for Z_2) action. The invariant
    subspace of H~_BL (x) H~_BR is the two-block collar model with
    (b_L, b_R) = (2, 2) per block. A and D are qubits.

    Returns (U, P_iso_plus, P_iso_minus, isometries) where U is the lifted
    diagonal boundary action on H_A (x) H~_BL (x) H~_BR (x) H_D (dim 64) and
    isometries[s] maps block s of the quotient model (dim 16) into the
    invariant subspace.
    """
    d = 64  # 2 * 4 * 4 * 2
    u = np.diag([1.0, 1.0, -1.0, -1.0])  # diag(1,-1) (x) I_2 on C^4
    U = np.kron(np.kron(np.eye(2), np.kron(u, u)), np.eye(2))

    # isotypic projectors of u on C^4 (flux sectors of the left action)
    p_plus = np.diag([1.0, 1.0, 0.0, 0.0])
    p_minus = np.diag([0.0, 0.0, 1.0, 1.0])

    # block isometries: |A> (x) |s,m>_L (x) |s,m'>_R (x) |D>  ->  block s
    isometries = {}
    for s, offset in ((0, 0), (1, 2)):  # s=0: +, s=1: -
        v = np.zeros((d, 16))
        col = 0
        for a in range(2):
            for m in range(2):
                for mp in range(2):
                    for dd in range(2):
                        row = ((a * 4 + (offset + m)) * 4 + (offset + mp)) * 2 + dd
                        v[row, col] = 1.0
                        col += 1
        isometries[s] = v
    return U, p_plus, p_minus, isometries


def _group_average(h: np.ndarray, U: np.ndarray) -> np.ndarray:
    return 0.5 * (h + U @ h @ U.conj().T)


def test_descent_central_interface_pre_quotient_gibbs_is_aligned():
    """Theorem thm:msaderivation via Lemma lem:onesideddescent, end to end.

    Pre-quotient H~ = H~_AL + H~_RD + H~_Sigma with H~_AL, H~_RD one-sided
    K-invariant and H~_Sigma a flux (isotypic-projector) function. The Gibbs
    state of the descended Hamiltonian on the invariant subspace is EC-aligned
    exact Markov on the preselected blocks.
    """
    rng = np.random.default_rng(23)
    U, p_plus, p_minus, isometries = _pre_quotient_model(rng)

    def rand_herm(n):
        x = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        return 0.5 * (x + x.conj().T)

    # one-sided K-invariant terms, lifted to the 64-dim pre-quotient space
    u4 = np.diag([1.0, 1.0, -1.0, -1.0])
    h_al = _group_average(rand_herm(8), np.kron(np.eye(2), u4))     # on A(x)BL
    h_rd = _group_average(rand_herm(8), np.kron(u4, np.eye(2)))     # on BR(x)D
    H_AL = np.kron(h_al, np.eye(8))
    H_RD = np.kron(np.eye(8), h_rd)
    # central interface energy: a function of the flux through Sigma
    H_SIG = np.kron(np.eye(2), np.kron(0.9 * p_plus - 0.4 * p_minus, np.eye(8)))
    H_pre = H_AL + H_RD + H_SIG

    # descend: restrict the Gibbs state of H_pre to the invariant subspace
    blocks = []
    weights = []
    for s in (0, 1):
        v = isometries[s]
        h_s = v.conj().T @ H_pre @ v
        rho_s = np.asarray(np.linalg.eigh(h_s)[1], dtype=complex)
        evals = np.linalg.eigh(h_s)[0]
        w = rho_s @ np.diag(np.exp(-evals)) @ rho_s.conj().T
        weights.append(np.trace(w).real)
        blocks.append(w)
    z = sum(weights)
    blocks = [(weights[s] / z, blocks[s] / weights[s], (2, 2, 2, 2))
              for s in (0, 1)]

    assert is_ec_aligned(blocks, tol=1e-8)
    assert modular_splitting_defect(blocks) < 1e-8
    assert collar_cmi(blocks) < 1e-8


def test_descent_invariant_but_noncentral_interface_breaks_alignment():
    """Sharp negative direction of the central-interface collar clause.

    A K-invariant cross-cut term that couples the multiplicity factors
    (rather than acting through the flux sectors) destroys alignment even
    though it respects the boundary symmetry. Invariance alone is not the
    clause; centrality is.
    """
    rng = np.random.default_rng(29)
    U, p_plus, p_minus, isometries = _pre_quotient_model(rng)

    def rand_herm(n):
        x = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        return 0.5 * (x + x.conj().T)

    u4 = np.diag([1.0, 1.0, -1.0, -1.0])
    h_al = _group_average(rand_herm(8), np.kron(np.eye(2), u4))
    h_rd = _group_average(rand_herm(8), np.kron(u4, np.eye(2)))
    H_AL = np.kron(h_al, np.eye(8))
    H_RD = np.kron(np.eye(8), h_rd)
    # K-invariant but NON-central cross-cut coupling on BL (x) BR
    cross = _group_average(rand_herm(16), np.kron(u4, u4))
    H_CROSS = np.kron(np.eye(2), np.kron(cross, np.eye(2)))
    H_pre = H_AL + H_RD + 1.5 * H_CROSS

    blocks = []
    weights = []
    for s in (0, 1):
        v = isometries[s]
        h_s = v.conj().T @ H_pre @ v
        vecs = np.asarray(np.linalg.eigh(h_s)[1], dtype=complex)
        evals = np.linalg.eigh(h_s)[0]
        w = vecs @ np.diag(np.exp(-evals)) @ vecs.conj().T
        weights.append(np.trace(w).real)
        blocks.append(w)
    z = sum(weights)
    blocks = [(weights[s] / z, blocks[s] / weights[s], (2, 2, 2, 2))
              for s in (0, 1)]

    assert entropic_alignment_defect(blocks) > 1e-3
    assert modular_splitting_defect(blocks) > 1e-3


# ---------------------------------------------------------------------------
# Independence of the central-interface collar clause (issue #544,
# Proposition prop:clauseindependence): a finite regulator package that
# satisfies the checked repair/consensus axiom-side conditions while its
# retained Axiom-3 family carries a K-invariant NON-central cross-cut density.
# ---------------------------------------------------------------------------


def _rand_herm(rng: np.random.Generator, n: int) -> np.ndarray:
    x = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return 0.5 * (x + x.conj().T)


def _noncentral_package(rng: np.random.Generator, g_cross: float = 1.5):
    """Retained Axiom-3 density family with an invariant noncentral coupling.

    Returns (U, isometries, densities) where densities = (H_AL, H_RD, H_X)
    are the three retained gauge-invariant local densities on the 64-dim
    pre-quotient space: two one-sided K-invariant terms and one group-averaged
    (hence K-invariant) but NON-central bL--bR cross-cut density, scaled by
    g_cross. This is exactly the coupling class of
    test_descent_invariant_but_noncentral_interface_breaks_alignment.
    """
    U, _, _, isometries = _pre_quotient_model(rng)
    u4 = np.diag([1.0, 1.0, -1.0, -1.0])
    h_al = _group_average(_rand_herm(rng, 8), np.kron(np.eye(2), u4))
    h_rd = _group_average(_rand_herm(rng, 8), np.kron(u4, np.eye(2)))
    H_AL = np.kron(h_al, np.eye(8))
    H_RD = np.kron(np.eye(8), h_rd)
    cross = _group_average(_rand_herm(rng, 16), np.kron(u4, u4))
    H_X = g_cross * np.kron(np.eye(2), np.kron(cross, np.eye(2)))
    return U, isometries, (H_AL, H_RD, H_X)


def _gibbs(h: np.ndarray) -> np.ndarray:
    evals, vecs = np.linalg.eigh(h)
    w = (vecs * np.exp(-(evals - evals.min()))) @ vecs.conj().T
    return w / np.trace(w).real


def _descend_gibbs(h_pre: np.ndarray, isometries) -> list:
    """Descend the Gibbs state of a pre-quotient Hamiltonian to collar blocks."""
    blocks, weights = [], []
    for s in (0, 1):
        v = isometries[s]
        h_s = v.conj().T @ h_pre @ v
        evals, vecs = np.linalg.eigh(h_s)
        w = (vecs * np.exp(-evals)) @ vecs.conj().T
        weights.append(np.trace(w).real)
        blocks.append(w)
    z = sum(weights)
    return [(weights[s] / z, blocks[s] / weights[s], (2, 2, 2, 2))
            for s in (0, 1)]


def _relative_entropy(rho: np.ndarray, sigma: np.ndarray) -> float:
    from msa_characterizations import logm_psd as _logm
    return float(np.trace(rho @ (_logm(rho) - _logm(sigma))).real)


def test_noncentral_model_satisfies_axioms():
    """Independence witness, quantum layer (Proposition prop:clauseindependence).

    The retained family {H_AL, H_RD, H_X} with H_X the invariant-but-noncentral
    cross-cut density passes every axiom-side condition checked here:

    (a) Axiom-3 admissibility: each density is self-adjoint and gauge
        (K_Sigma-)invariant, and the family plus the identity is linearly
        independent, so the MaxEnt multipliers are unique
        (lem:closure-residual, item i);
    (b) Axiom 2 (overlap consistency): restrictions of the realized Gibbs
        state to the two overlapping patches {A,BL,BR} and {BL,BR,D} agree
        on the shared half-collar pair;
    (c) refinement closure (Axiom-3 clause): the declared ancilla-extension
        refinement channel (extra regulator cells disjoint from the collar,
        coarse-graining = partial trace) has closure defect exactly zero at
        every tested multiplier vector, with identity induced multiplier map,
        so the noncentral coupling survives to the finer stage unchanged;

    while the central-interface collar clause FAILS: the descended MaxEnt
    state is not EC-aligned. Repair/consensus and closure therefore do not
    force the clause; it is an independent axiom-level input.
    """
    rng = np.random.default_rng(31)
    U, isometries, densities = _noncentral_package(rng)
    H_AL, H_RD, H_X = densities

    # (a) admissibility: invariance and linear independence
    for h in densities:
        assert np.linalg.norm(h - h.conj().T) < 1e-12
        assert np.linalg.norm(h @ U - U @ h) < 1e-12, \
            "every retained density must be gauge-invariant"
    ops = [H_AL, H_RD, H_X, np.eye(64)]
    gram = np.array([[np.trace(a.conj().T @ b).real for b in ops] for a in ops])
    assert np.linalg.matrix_rank(gram, tol=1e-8) == 4, \
        "densities plus identity must be linearly independent"

    dims = [2, 4, 4, 2]
    for lam in ([1.0, 1.0, 1.0], [0.7, 1.3, 0.9], [1.1, 0.6, 1.4]):
        h_eff = lam[0] * H_AL + lam[1] * H_RD + lam[2] * H_X
        omega = _gibbs(h_eff)

        # (b) overlap consistency of the realized state
        rho_p1 = partial_trace(omega, dims, [0, 1, 2])       # patch {A,BL,BR}
        rho_p2 = partial_trace(omega, dims, [1, 2, 3])       # patch {BL,BR,D}
        shared_1 = partial_trace(rho_p1, [2, 4, 4], [1, 2])
        shared_2 = partial_trace(rho_p2, [4, 4, 2], [0, 1])
        assert np.linalg.norm(shared_1 - shared_2) < 1e-12

        # (c) refinement closure: fine stage adds two ancilla cells away from
        # the collar; the channel is the partial trace over them.
        h_fine = np.kron(h_eff, np.eye(4))
        omega_fine = _gibbs(h_fine)
        coarse_grained = partial_trace(omega_fine, dims + [4], [0, 1, 2, 3])
        assert _relative_entropy(coarse_grained, omega) < 1e-9, \
            "closure defect must vanish (def:closure-defect)"
        for h in densities:  # identity induced multiplier map: moments match
            m_fine = np.trace(omega_fine @ np.kron(h, np.eye(4))).real
            m_coarse = np.trace(omega @ h).real
            assert abs(m_fine - m_coarse) < 1e-10

        # the clause fails at every stage: the MaxEnt state is not aligned
        blocks = _descend_gibbs(h_eff, isometries)
        assert entropic_alignment_defect(blocks) > 1e-3
        assert modular_splitting_defect(blocks) > 1e-3
        assert not is_ec_aligned(blocks)


def _rooted_extension(boundary: int, edge_labels: tuple) -> tuple:
    """E(b) of thm:functional-selected-fiber on the Z_2 seam chain."""
    vals = [boundary]
    for c in edge_labels:
        vals.append(vals[-1] ^ c)
    return tuple(vals)


def _enabled_repairs(state: tuple, boundary: int, edge_labels: tuple):
    """Touched-overlap transactions enabled at `state`.

    Transaction T_i proposes s_i := s_{i-1} XOR c_i (the recovery-derived
    seam proposal), reads the two seams meeting register i
    (validation-complete read set), writes only register i >= 1 (boundary
    preserved), and commits only under strict lexicographic descent of
    mu = (Phi, N_unresolved), the declared well-founded measure of
    def:transactional-quotient-repair.
    """
    full = (boundary,) + state
    ext = _rooted_extension(boundary, edge_labels)

    def mu(f):
        phi = sum(1 for i, c in enumerate(edge_labels) if f[i] ^ c != f[i + 1])
        unresolved = sum(1 for a, b in zip(f, ext) if a != b)
        return (phi, unresolved)

    moves = []
    for i in range(1, len(full)):
        proposal = full[i - 1] ^ edge_labels[i - 1]
        if proposal == full[i]:
            continue
        new = full[:i] + (proposal,) + full[i + 1:]
        if mu(new) < mu(full):  # touched-overlap acceptance: strict descent
            moves.append((i, new[1:]))
    return moves


def test_noncentral_model_repair_confluence():
    """Independence witness, classical repair layer.

    The flux-sector patch net carrying the noncentral package is a layered
    functional boundary carrier (thm:layered-carrier-HB-Hfib). Exhaustive
    enumeration over every initial state and every schedule of accepted
    touched-overlap transactions verifies: strict descent on accepted moves,
    termination, repair completeness (terminal states are exactly the
    consistent fiber), and schedule-independent confluence to the rooted
    extension E(b) (thm:confluence, thm:boundary-conditioned-uniqueness).
    The repair layer reads only interface sector data; the settled net then
    realizes the noncentral quantum package with no enabled transaction, so
    overlap-consistent repair is inert on it and cannot exclude the coupling.
    """
    from itertools import product

    for boundary, edge_labels in ((0, (0, 1, 1)), (1, (1, 0, 1)), (0, (0, 0, 0))):
        ext = _rooted_extension(boundary, edge_labels)
        for init in product((0, 1), repeat=len(edge_labels)):
            terminals = set()
            stack = [tuple(init)]
            seen = set()
            while stack:
                state = stack.pop()
                if state in seen:
                    continue
                seen.add(state)
                moves = _enabled_repairs(state, boundary, edge_labels)
                if not moves:
                    terminals.add(state)
                    continue
                for _, new_state in moves:
                    stack.append(new_state)
            # confluence + completeness: unique terminal = consistent fiber
            assert terminals == {ext[1:]}, \
                f"schedules from {init} must settle on E(b)"

    # the settled net realizes the noncentral package; repair has no move
    # left, yet the realized MaxEnt state breaks alignment.
    rng = np.random.default_rng(37)
    _, isometries, (H_AL, H_RD, H_X) = _noncentral_package(rng)
    settled = _rooted_extension(0, (0, 1, 1))
    assert not _enabled_repairs(settled[1:], 0, (0, 1, 1))
    blocks = _descend_gibbs(H_AL + H_RD + H_X, isometries)
    assert entropic_alignment_defect(blocks) > 1e-3
    assert not is_ec_aligned(blocks)


def test_noncentral_coupling_passes_small_cmi_budget():
    """The small-CMI clause of Axiom 4 does not exclude the coupling class.

    At coupling strength g the collar CMI of the descended Gibbs state
    scales quadratically while the modular-splitting defect scales linearly.
    Every finite CMI budget therefore admits invariant noncentral couplings
    whose distance from the EC-aligned class is parametrically larger than
    the budget, so collar recoverability cannot force the clause.
    """
    rng = np.random.default_rng(41)
    _, isometries, (H_AL, H_RD, H_X) = _noncentral_package(rng, g_cross=1.0)

    strengths = (0.4, 0.2, 0.1)
    cmi_vals, msd_vals = [], []
    for g in strengths:
        blocks = _descend_gibbs(H_AL + H_RD + g * H_X, isometries)
        cmi_vals.append(collar_cmi(blocks))
        msd_vals.append(modular_splitting_defect(blocks))

    # quadratic CMI decay versus linear splitting-defect decay
    assert cmi_vals[1] / cmi_vals[0] < 0.35      # expect ~ (1/2)^2
    assert cmi_vals[2] / cmi_vals[1] < 0.35
    assert msd_vals[1] / msd_vals[0] > 0.4       # expect ~ 1/2
    assert msd_vals[2] / msd_vals[1] > 0.4
    # the misalignment outruns the CMI budget as g decreases
    assert msd_vals[2] / cmi_vals[2] > msd_vals[0] / cmi_vals[0]
    assert msd_vals[2] > 10.0 * cmi_vals[2]


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-v"]))
