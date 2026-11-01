"""Exact full-field eta lift, actual complete repair and interacting response."""
import math

import numpy as np
import pytest
from scipy import sparse
import sympy as sy

from . import cube, repair, spectrum


def test_complete_3d_action_eta_injection_and_both_intertwiners():
    result = cube.certificate()
    assert result['protected_eta_image_dimension'] == math.comb(16, 4) == 1820
    assert result['eta3_parent_dimension'] == 119
    assert result['spin_frame_alignment_used_for_protection'] is False
    assert cube.eta_image().shape == (8008, 1820)
    assert sum(bool(group) for group in cube.GEOMETRY.groups) == 3


def test_actual_four_dimer_support_all_seams_and_all_outcomes():
    result = repair.certificate()
    assert result['complete_outcomes_including_zero_and_outside'] == 130
    assert result['admitted_native_seams'] == [[0, 1], [0, 5], [0, 7], [1, 5], [1, 7]]
    assert result['repair_support_dimension'] == 24
    assert result['complete_full_Fock_domain'] == 65536


def test_complete_repair_preserves_arbitrary_reference_cross_terms_in_image():
    rng = np.random.default_rng(3701)
    # Include complete parent columns with singly occupied sites; not only eta condensates.
    parent = np.zeros((1820, 3), complex)
    parent[[0, 3, 280, 801, 1701]] = rng.normal(size=(5, 3))+1j*rng.normal(size=(5, 3))
    image = cube.eta_image()@parent
    image /= np.linalg.norm(image)
    state = np.zeros((65536, 3), complex); state[cube.weights(6)] = image
    for law in ('C', 'E'):
        for seam in repair.SEAMS:
            out = repair.apply(state, law, seam)
            assert sum(np.linalg.norm(value)**2 for label, value in out[:-1]) < 1e-30
            np.testing.assert_allclose(out[-1][1], state, rtol=0, atol=1e-15)
            actual = sum(value.conj().T@value for _, value in out)
            np.testing.assert_allclose(actual, state.conj().T@state, rtol=0, atol=1e-14)


def test_unprotected_native_support_is_not_misreported_as_the_dark_space():
    state = repair.embedding()[:, 0].toarray()
    outcomes = repair.apply(state, 'E', repair.SEAMS[0])
    assert sum(np.linalg.norm(value)**2 for label, value in outcomes[:-1]) == pytest.approx(1)
    assert np.linalg.norm(outcomes[-1][1]) < 1e-14
    with pytest.raises(ValueError): repair.matrices('E', (2, 3))


def test_interacting_119_state_spectrum_and_positive_eight_line_probe():
    result = spectrum.polynomial_certificate()
    assert result['eta3_distinct_energies'] == 11
    assert result['genuinely_contact_shifted_dimensions'] == 20
    report = spectrum.local_probe_certificate()
    lines = report['eight_interacting_energy_lines']
    assert len(lines) == 8
    assert all(row['normalized_probe_weight'] > 0 for row in lines)
    assert sum(row['normalized_probe_weight'] for row in lines) == pytest.approx(1)
    assert report['normalized_probe_H_variance'] == '3/112'
    assert report['charged_vacuum_superposition_required'] is False
    assert report['exact_full_exterior_reduced_state_unchanged'] is True
    assert report['phase_probe_energy_variance'] == '45/1792'


def test_exact_probe_resolvent_agrees_off_spectrum_with_direct_linear_solve():
    numerator, denominator, determinant = spectrum.probe_resolvent()
    _, _, gram, h, _, q = spectrum.singlet_data()
    assert determinant == 1488395072912184839188800532906554156908544
    for z in (sy.Integer(20), sy.Integer(30)):
        direct = (sy.Matrix(q).T*sy.diag(*map(int, gram))*(z*sy.eye(36)-sy.Matrix(h)).inv()*sy.Matrix(q))[0]/56
        assert direct == numerator.subs(spectrum.Z, z)/denominator.subs(spectrum.Z, z)


def test_local_phase_has_the_claimed_fixed_number_elastic_and_inelastic_parts():
    from .source import ideal_reference
    reference = ideal_reference()
    n = np.array([(mask & 3).bit_count() for mask in range(65536)])
    probe = (4*n-3)*reference/math.sqrt(15)
    kicked = np.exp(1j*math.pi*n/2)*reference
    np.testing.assert_allclose(kicked, (reference-math.sqrt(15)*probe)/4, atol=2e-16, rtol=0)
    assert np.vdot(reference, kicked) == pytest.approx(.25)
    assert np.linalg.norm(probe) == pytest.approx(1)
    assert sum(abs(reference[j])**2 for j in range(65536) if j.bit_count() != 6) == 0


def test_changed_contact_action_does_not_satisfy_the_retained_119_polynomial():
    _, _, gram, h, _, _ = spectrum.singlet_data()
    changed = h.copy()
    for j, x in enumerate(np.diag(h)):
        if x: changed[j, j] += 16
    assert sy.Matrix(changed).charpoly(spectrum.Z).as_expr() != sy.Matrix(h).charpoly(spectrum.Z).as_expr()
