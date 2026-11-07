"""Complete-source and independent instrument checks for probe uncertainty."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import uncertainty as u


@pytest.fixture(scope="module")
def probes():
    return u.load_probes()


def test_exact_certificate_covers_all_six_complete_probes():
    result = u.certificate()
    assert len(result["probes"]) == 6
    for row in result["probes"]:
        assert (row["neutral_dimension"], row["charged_dimension"]) == (438, 1596)
        assert row["rank_each_signed_outcome"] == 384
        assert row["rank_blocked_outcome"] == 1266
        assert row["neutral_blocked_dimension"] == 54
        assert row["neutral_tracial_probabilities_minus_zero_plus"] == ["32/73", "9/73", "32/73"]
    assert not result["scope"]["physical_position_momentum_identified"]
    assert not result["scope"]["new_native_compiler_for_projective_instrument"]


def test_actual_neutral_states_saturate_with_every_blocked_branch(probes):
    rng = np.random.default_rng(2063)
    for probe in probes:
        neutral = rng.normal(size=(438, 3)) + 1j*rng.normal(size=(438, 3))
        neutral /= np.linalg.norm(neutral)
        state = np.vstack((neutral, np.zeros((1596, 3))))
        b = probe.creation()
        w = np.linalg.norm(b @ neutral)**2
        assert 0 < w < 1
        matrices = {}
        for axis in ("X", "Y"):
            outcomes = probe.outcomes(state, axis)
            probabilities = {k: np.linalg.norm(v)**2 for k, v in outcomes.items()}
            np.testing.assert_allclose([probabilities[k] for k in (-1, 0, 1)],
                                       [w/2, 1-w, w/2], rtol=0, atol=2e-14)
            effects = probe.projectors(axis)
            matrices[axis] = effects[1] - effects[-1]
            assert abs(np.vdot(state, matrices[axis] @ state)) < 2e-14
            assert abs(np.vdot(matrices[axis] @ state, matrices[axis] @ state).real-w) < 2e-14
        x, y = matrices["X"], matrices["Y"]
        np.testing.assert_allclose(np.vdot(state, x @ (y @ state) - y @ (x @ state)),
                                   2j*w, rtol=0, atol=2e-14)
        # The failure to retain zero would lose real probability on this input.
        assert probabilities[-1] + probabilities[1] < 1


def test_complete_instrument_preserves_arbitrary_reference_marginal(probes):
    rng = np.random.default_rng(991)
    for probe in probes:
        state = rng.normal(size=(probe.dimension, 4)) + 1j*rng.normal(size=(probe.dimension, 4))
        state /= np.linalg.norm(state)
        for axis in ("X", "Y"):
            branches = probe.outcomes(state, axis)
            marginal = sum((v.conj().T @ v for v in branches.values()), np.zeros((4, 4), complex))
            np.testing.assert_allclose(marginal, state.conj().T @ state, rtol=0, atol=2e-14)


def test_spectral_projectors_against_independent_diagonalization():
    probe = u.Probe("control", 3, 4, ((2, 0, -1), (0, 2, 1)))
    b = probe.creation().toarray()
    for axis, phase in (("X", 1), ("Y", -1j)):
        full = np.block([[np.zeros((3, 3)), phase*b.conj().T],
                         [np.conjugate(phase)*b, np.zeros((4, 4))]])
        eigenvalues, vectors = np.linalg.eigh(full)
        for label, effect in probe.projectors(axis).items():
            eigenspace = vectors[:, np.isclose(eigenvalues, label, atol=1e-12)]
            np.testing.assert_allclose(effect.toarray(), eigenspace @ eigenspace.conj().T,
                                       rtol=0, atol=2e-14)
    blocked = np.eye(7)[:, 1]
    for axis in ("X", "Y"):
        out = probe.outcomes(blocked, axis)
        np.testing.assert_array_equal(out[0], blocked)
        assert np.linalg.norm(out[-1]) == np.linalg.norm(out[1]) == 0


@pytest.mark.parametrize("entries", [((0, 0, 1), (0, 1, 1)), ((0, 0, 1), (1, 0, 1)),
                                    ((0, 0, 2),), ((4, 0, 1),)])
def test_invalid_partial_permutations_rejected(entries):
    with pytest.raises(ValueError):
        u.Probe("bad", 3, 4, entries)


def test_nonunit_exact_normalization_is_rejected():
    with pytest.raises(ValueError):
        u._probe_from_rows("bad", [(0, 0, "1")], [F(2)], [F(1)])


def test_instrument_input_domains(probes):
    probe = probes[0]
    with pytest.raises(ValueError):
        probe.projectors("Z")
    for state in (np.ones(3), np.full(probe.dimension, np.nan), np.zeros((probe.dimension, 2, 2))):
        with pytest.raises(ValueError):
            probe.outcomes(state)
