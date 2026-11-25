"""Independent checks of native block structure and the energy-reference boundary."""

import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

from model.candidate_laws.interval_echo.openwave.row_03_clock.neutral_reference import measure


def test_exact_local_uniform_block_and_zero_generators():
    root = sp.sqrt(2) / 2
    echo = sp.Matrix([[0, root, root], [root, sp.Rational(1, 2), -sp.Rational(1, 2)],
                      [root, -sp.Rational(1, 2), sp.Rational(1, 2)]])
    assert echo * echo == sp.eye(3)
    u, h, _, indices, _, _ = measure.primitives()
    np.testing.assert_allclose(u[np.ix_(indices, indices)], np.array(echo, dtype=float), atol=1e-15)
    assert np.max(np.abs(h[:, indices])) < 1e-15


def test_complete_adaptive_native_records_ignore_neutral_source_coherence():
    u, h, _, _, filled, pair = measure.primitives()
    phi = np.exp(-0.371j) * pair
    pure = measure.transcript(u, h, (filled + phi) / np.sqrt(2))
    mixture = (measure.transcript(u, h, filled) + measure.transcript(u, h, phi)) / 2
    assert pure.shape == (99,)
    assert sum(pure) == pytest.approx(1, abs=2e-13)
    assert np.max(np.abs(pure - mixture)) < 2e-13


def test_block_mixing_analyzer_is_a_detectable_change():
    _, _, _, _, filled, pair = measure.primitives()
    plus = (filled + pair) / np.sqrt(2)
    assert abs(np.vdot(plus, plus)) ** 2 == pytest.approx(1)
    mixed_plus = (abs(np.vdot(plus, filled))**2 + abs(np.vdot(plus, pair))**2) / 2
    assert mixed_plus == pytest.approx(0.5)
    projector = np.outer(plus, plus.conj())
    filled_projector = np.outer(filled, filled)
    assert np.linalg.norm(projector @ filled_projector - filled_projector @ projector) > 0.7


def test_uniform_shift_and_controlled_drive_are_distinct():
    _, h, charge, _, _, _ = measure.primitives()
    result = measure.branch_boundary(h, charge)
    assert result['native_plus'] == pytest.approx((2 + np.sqrt(2)) / 4, abs=2e-13)
    assert result['shifted_every_branch_plus'] == pytest.approx(result['native_plus'], abs=2e-13)
    assert result['conditional_shift_only_plus'] == pytest.approx((2 - np.sqrt(2)) / 4, abs=2e-13)


def test_balancing_reference_and_adaptive_branch_phases():
    # A sample can exchange charge with an explicit reference. In this fixed
    # total-charge subspace its two basis vectors have sample charges -1,0
    # and reference charges 0,-1. The total charge is -I, not sample charge.
    root = sp.sqrt(2) / 2
    qsample, qreference = sp.diag(-1, 0), sp.diag(0, -1)
    qtotal = qsample + qreference
    assert qtotal == -sp.eye(2)
    splitter = sp.Matrix([[root, sp.I * root], [sp.I * root, root]])
    h = sp.diag(35, 0)
    # Arbitrary classical read outcomes admit different waits on each branch.
    t, lam = sp.symbols('t lam', real=True)
    phases = [sp.exp(-sp.I * 35 * t), 1]
    native = splitter * sp.diag(*phases) * splitter
    shifted = splitter * sp.diag(*[sp.exp(sp.I * lam * t) * p for p in phases]) * splitter
    assert sp.simplify(shifted - sp.exp(sp.I * lam * t) * native) == sp.zeros(2)
    # Shifting only the sample instead changes its relative transfer energy.
    assert qsample != qtotal
    assert h - 16 * qsample == sp.diag(51, 0)


def test_two_convention_minimax_bound_and_neutral_sum():
    assert 35 + 67 == 51 + 51 == 102
    for estimate in np.linspace(20, 70, 101):
        assert max(abs(estimate - 35), abs(estimate - 51)) >= 8
    assert max(abs(43 - 35), abs(43 - 51)) == 8


def test_negative_and_unnormalized_probabilities_fail_instead_of_being_hidden():
    rng = np.random.default_rng(12)
    with pytest.raises(ArithmeticError, match="negative"):
        measure.sample_counts(rng, [-1e-14, 1 + 1e-14])
    with pytest.raises(ArithmeticError, match="normalization"):
        measure.sample_counts(rng, [0.2, 0.2])
    counts, bound = measure.sample_counts(rng, [0.5, 0.5 + 1e-14])
    assert sum(counts) == measure.SHOTS
    assert 0 < bound < 2e-12


def test_public_replay_command_uses_json_compatible_scalars():
    result = subprocess.run([sys.executable, '-B', '-m', measure.__package__ + '.measure'],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert 'Replayed native neutral-record' in result.stdout


def test_same_source_pair_reader_detects_charge_and_swap_preserving_changed_law():
    _, h, charge, _, filled, pair = measure.primitives()
    result = measure.changed_creation_control(h, charge, filled, pair)
    assert result['native_pair_census_probability'] == 0
    assert result['changed_law_pair_census_probability'] == pytest.approx(0.5, abs=2e-13)
