"""Independent circuit, cutoff-complement and complete-record checks."""
from fractions import Fraction as F
import numpy as np
import sympy as sp

from . import model


def test_exact_finite_car_and_normalized_rotation():
    result = model.verify_full_car()
    assert result["occupation_words"] == 8192
    assert result["checked_mixed_CAR_columns"] == 1384448
    assert result["checked_creator_anticommutator_columns"] == 1384448


def test_kinematic_boundary_rotation_requires_explicit_physical_extension():
    # Two matter modes, R in0..2, one longitudinal bit. All24 states retained.
    states = [(word, r, ell) for word in range(4) for r in range(3) for ell in range(2)]
    index = {state: j for j, state in enumerate(states)}
    size = len(states)
    physical = [r == word.bit_count() and ell == bool(word & 2) for word, r, ell in states]
    P = sp.diag(*map(int, physical))
    creators = []
    for mode in range(2):
        D = sp.zeros(size)
        for col, (word, r, ell) in enumerate(states):
            action = model.car_step(word, mode, True, 2)
            if action is None or r == 2 or (mode == 1 and ell):
                continue
            word2, sign = action
            target = (word2, r+1, ell+int(mode == 1))
            D[index[target], col] = sign
        creators.append(D)
    D = sp.Rational(3, 5)*creators[0] + sp.Rational(4, 5)*creators[1]
    A = D-D.T
    naive = sp.sqrt(3)/2*sp.eye(size)+A/2
    assert naive.T*naive != sp.eye(size)
    U = P*naive*P + sp.eye(size)-P
    assert sp.simplify(U.T*U) == sp.eye(size)
    # Gauss at source1 is ell-n1; sum Gauss is R-N.
    G1 = sp.diag(*[ell-int(bool(word & 2)) for word, r, ell in states])
    Gsum = sp.diag(*[r-word.bit_count() for word, r, ell in states])
    assert U*G1 == G1*U and U*Gsum == Gsum*U
    vacuum = sp.zeros(size, 1); vacuum[index[(0, 0, 0)]] = 1
    output = U*vacuum
    assert sum(output[j]**2 for j, (w, _, _) in enumerate(states) if w.bit_count() == 0) == sp.Rational(3, 4)
    assert sum(output[j]**2 for j, (w, _, _) in enumerate(states) if w.bit_count() == 1) == sp.Rational(1, 4)


def test_fourier_sign_from_literal_controlled_circuit_and_all_outcomes():
    size = 16
    labels = np.arange(size)
    Fplus = np.exp(2j*np.pi*np.outer(labels, labels)/size)/np.sqrt(size)
    for energy in (0, 3, 7, 11, 1.8738097549661659):
        controlled_state = np.exp(-2j*np.pi*labels*energy/16)/np.sqrt(size)
        output = Fplus @ controlled_state
        direct = model.phase_kernel(energy, size)
        np.testing.assert_allclose(output, direct, atol=2e-15)
        assert abs(np.vdot(output, output)-1) < 2e-15
    correct = Fplus @ (np.exp(-2j*np.pi*labels*3/16)/np.sqrt(size))
    reversed_sign = Fplus.conj() @ (np.exp(-2j*np.pi*labels*3/16)/np.sqrt(size))
    assert np.argmax(abs(correct)) == 3
    assert np.argmax(abs(reversed_sign)) == 13


def test_uniform_grid_tail_includes_phase_alignment_and_wrapping():
    for size in (16, 64, 257):
        grid = np.arange(size)/size
        for theta in (0, .003, .499, .501, .999, np.sqrt(2)/4):
            probabilities = abs(model.phase_kernel(16*theta, size))**2
            assert abs(probabilities.sum()-1) < 3e-14
            distances = abs(grid-theta)
            distances = np.minimum(distances, 1-distances)
            for distance in (F(1, 100), F(1, 20), F(1, 5), F(2, 5)):
                actual = probabilities[distances >= float(distance)].sum()
                assert actual <= float(model.tail_bound(size, distance)) + 1e-14


def test_full_phase_pointer_completeness_diagnostic():
    # All2^20 pointer outcomes, on a specified source-energy eigenphase.
    energy = float(model.full.CENTER)
    probabilities = abs(model.phase_kernel(energy))**2
    assert probabilities.shape == (2**20,)
    assert abs(probabilities.sum()-1) < 1e-12
    lo, hi = model.acceptance_indices(model.full.CENTER)
    actual_rejection = probabilities[:lo].sum()+probabilities[hi+1:].sum()
    assert actual_rejection < float(model.tail_bound(model.M, (model.HALF_WINDOW-model.full.RADIUS)/16))


def test_rational_clock_global_identity_and_error():
    identities = model.rational_clock_identities()
    assert "100800" in identities["single_step_norm_error"]
    assert model.clock_error() < F(7409, 10**15)
    assert 2*model.SUBSTEPS*(model.M-1) == 8589926400
    # Actual spectral values at both signs; no small-step assumption is used
    # for the theorem. Large x controls unitarity and the exact denominator.
    for x in (-20., -1., 0., .3, 3., 20.):
        a, b = 1-x*x/10, x/2-x**3/120
        R = (a-1j*b)/(a+1j*b)
        assert abs(abs(R)-1) < 4e-16
        assert abs(R-np.exp(-1j*x)) <= abs(x)**7/100800 + 1e-15


def test_source_replay_and_complete_conditional_unconditional_bounds():
    saved = model.source()
    result = model.record_bounds(saved)
    lo, hi = map(F, result["conditional_detector_acceptance_interval"])
    ulo, uhi = map(F, result["unconditional_accepted_record_interval"])
    rlo, rhi = map(F, result["unconditional_rejected_record_interval"])
    assert lo > F("0.99967332932") and hi <= 1
    assert ulo <= lo/4 and uhi >= hi/4
    assert ulo+rlo <= F(1, 4) <= uhi+rhi
    assert F(result["unconditional_excited_false_acceptance_upper"]) < F("1.400761e-9")
    assert F(result["accepted_source_ground_fidelity_lower"]) > F("0.99999999439512")
    assert result["accepted_pointer_outcomes"]+result["rejected_pointer_outcomes"] == model.M
    assert result["exact_ground_projector_implemented"] is False


def test_public_receipt_preserves_frozen_science_and_portable_closure():
    from . import certificate as cert
    result = cert.certificate()
    assert result["independent_research_comparison"]["all_scientific_fields_equal"]
    assert result["schema"] == "oph.affine_gauss_charged_clock.v1"
    assert all(not name.startswith("/") and ".." not in name.split("/") for name in result["input_sha256"])
    cert.validate_pins(result["input_sha256"])
    assert "model/affine_gauss_charged_clock/research_reference.json" in result["input_sha256"]
    assert result["finite_clock"]["large_source_QPE_matrix_numerically_executed"] is False
