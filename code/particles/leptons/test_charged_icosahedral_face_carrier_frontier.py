"""Tests for the source-only icosahedral face-corner carrier frontier."""

from __future__ import annotations

import numpy as np

import derive_charged_icosahedral_face_carrier_frontier as lane


def test_icosahedral_incidence_is_12_30_20_with_euler_two():
    incidence = lane.icosahedral_incidence()
    assert (incidence["vertices"], incidence["edges"], incidence["faces"]) == (12, 30, 20)
    assert incidence["euler_characteristic"] == 2
    assert set(incidence["vertex_degrees"]) == {5}
    assert set(incidence["face_sizes"]) == {3}
    assert set(incidence["edge_face_counts"]) == {2}
    assert incidence["face_corner_flag_count"] == 60


def test_face_stabilizer_has_regular_c3_corner_action():
    shift = lane.cyclic_shift()
    assert np.allclose(np.linalg.matrix_power(shift, 3), np.eye(3))
    assert not np.allclose(shift, np.eye(3))
    assert 60 // lane.icosahedral_incidence()["faces"] == 3


def test_hermitian_circulant_commutes_with_face_shift_and_has_real_spectrum():
    shift = lane.cyclic_shift().astype(complex)
    matrix = lane.hermitian_circulant(1.0, 0.2, 0.1)
    assert np.allclose(matrix, np.conjugate(matrix.T))
    assert np.allclose(matrix @ shift, shift @ matrix)
    spectrum = lane.circulant_spectrum(1.0, 0.2, 0.1)
    assert len(spectrum) == 3
    assert min(spectrum) > 0.0
    assert np.allclose(
        spectrum,
        lane.analytic_circulant_spectrum(1.0, 0.2, 0.1),
    )
    assert np.allclose(
        spectrum,
        lane.circulant_spectrum(1.0, 0.2, -0.1),
    )


def test_c3_symmetry_does_not_select_amplitude_or_phase():
    left = lane.circulant_spectrum(1.0, 0.20, 0.10)
    right = lane.circulant_spectrum(1.0, 0.27, 0.23)
    assert left != right
    assert min(left) > 0.0
    assert min(right) > 0.0


def test_artifact_promotes_geometry_but_not_charged_masses():
    screen = {
        "status": "closed_on_declared_triangulated_screen_branch",
        "checks": {"icosahedral_orbit_has_twelve_vertices": True},
        "orbit_stabilizer": {
            "orbit_size": 12,
            "group_order": 60,
            "group": "I ~= A5",
            "fivefold_stabilizer_order": 5,
        },
    }
    artifact = lane.build_artifact(
        screen,
        screen_sha256="frozen-screen-hash",
    )
    assert artifact["checks_pass"] is True
    assert artifact["charged_reference_data_consumed"] is False
    assert artifact["public_charged_mass_promotion_allowed"] is False
    assert artifact["screen_source_sha256"] == "frozen-screen-hash"
    assert artifact["status"] == (
        "CLOSED_GEOMETRIC_C3_FACE_CARRIER_PHYSICAL_ATTACHMENT_AND_VALUE_LAWS_OPEN"
    )


def test_malformed_screen_certificate_is_rejected():
    with np.testing.assert_raises(ValueError):
        lane.build_artifact({"orbit_stabilizer": {"orbit_size": 12}})
