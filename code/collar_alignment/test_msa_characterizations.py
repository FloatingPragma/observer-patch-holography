#!/usr/bin/env python3
"""Tests for the MSA characterization package (paper-audit issue 001, GitHub #543).

Each test verifies one clause of the compact paper's Proposition `prop:msachar`
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


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-v"]))
