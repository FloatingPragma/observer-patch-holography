"""Reject physical-spin transfers of scoped static-flux identities."""
import json
import numpy as np
from oph_openwave.flux_controls import run, static_flux_laplacian, transport_rotation


def test_finite_controls_and_scope():
    result = run()
    assert all(result['checks'].values())
    assert all(result['mutations'].values())
    assert result['full_criterion_met'] is False
    assert json.loads(json.dumps(result))['full_criterion_met'] is False
    assert result['measurements']['compact_residual'] < 1e-12


def test_rephasing_preserves_symmetry_and_changes_cube():
    h = static_flux_laplacian({(0, 1): 3, (0, 5): 3, (6, 8): 3,
                               (7, 8): 3, (1, 7): 3})
    r = transport_rotation()
    assert np.max(abs(r @ h - h @ r)) < 1e-12
    assert np.max(abs((-r) @ h - h @ (-r))) < 1e-12
    np.testing.assert_array_equal(np.linalg.matrix_power(r, 3), -np.eye(12))
    np.testing.assert_array_equal(np.linalg.matrix_power(-r, 3), np.eye(12))


def test_doublet_needs_flux_placement():
    adjacent = np.linalg.eigvalsh(static_flux_laplacian({(0, 1): 3}))
    assert adjacent[1] - adjacent[0] > 2
    bare = np.linalg.eigvalsh(static_flux_laplacian({}))
    assert abs(bare[0]) < 1e-12
