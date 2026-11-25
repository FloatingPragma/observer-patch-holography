"""Independent physical geometry, exact isometry and complete-history checks."""

import json
import subprocess
import sys

import numpy as np
import sympy as sp

from model.executable_model.carrier import cell, icosahedron_from_coordinates, isomorphism
from model.experiments.interval_echo_angular.charged_orientation.projected_turns import measure as m


def test_actual_half_turn_geometry_and_perpendicular_axes():
    R, S, _, _, _, _ = m.setup()
    points, edges = icosahedron_from_coordinates()
    positions = points[isomorphism(cell().seams, edges)]
    axes = []
    for op in (R, S):
        permutation = np.argmax(op, axis=0)
        physical = np.linalg.lstsq(positions, positions[permutation], rcond=None)[0]
        np.testing.assert_allclose(physical.T @ physical, np.eye(3), atol=1e-13)
        assert abs(np.linalg.det(physical) - 1) < 1e-13
        assert abs(np.trace(physical) + 1) < 1e-13
        _, vectors = np.linalg.eigh(physical)
        axes.append(vectors[:, -1])
    assert abs(np.dot(*axes)) < 1e-13


def test_exact_symbolic_isometry_saturates_threshold_and_positive_literal_full_turn():
    a = sp.sqrt(2 + sp.sqrt(2)) / (2 * sp.sqrt(2))
    b = sp.sqrt(2 - sp.sqrt(2)) / (2 * sp.sqrt(2))
    n = sp.Matrix([[a, b], [b, a], [-a, b], [-b, a]])
    r, s = sp.diag(1, 1, -1, -1), sp.diag(1, -1, 1, -1)
    A, B = n.T * r * n, n.T * s * n
    assert sp.simplify(n.T * n) == sp.eye(2)
    assert sp.simplify(A * B + B * A) == sp.zeros(2)
    assert sp.simplify(A * A) == sp.eye(2) / 2
    assert sp.simplify(B * B) == sp.eye(2) / 2
    assert sp.simplify(A**4) == sp.eye(2) / 4


def test_commutator_leakage_identity_on_arbitrary_coherent_codes():
    R, S, _, _, _, _ = m.setup()
    rng = np.random.default_rng(711)
    for _ in range(12):
        V, _ = np.linalg.qr(rng.normal(size=(12, 3)) + 1j * rng.normal(size=(12, 3)))
        P = V @ V.conj().T
        Q = np.eye(12) - P
        A, B = V.conj().T @ R @ V, V.conj().T @ S @ V
        actual = A @ B - B @ A
        identity = -V.conj().T @ R @ Q @ S @ V + V.conj().T @ S @ Q @ R @ V
        np.testing.assert_allclose(actual, identity, atol=1e-13)
        assert np.linalg.norm(actual, 2) <= 2 * np.linalg.norm(Q @ R @ V, 2) * np.linalg.norm(Q @ S @ V, 2) + 1e-13


def test_untouched_reference_and_both_rejection_stages():
    R, S, V, P, _, _ = m.setup()
    bell = V.reshape(-1) / np.sqrt(2)
    result = m.complete([np.kron(R, np.eye(2)), np.kron(S, np.eye(2))],
                        [np.kron(S, np.eye(2)), np.kron(R, np.eye(2))],
                        np.kron(P, np.eye(2)), bell, True, np.random.default_rng(2))
    np.testing.assert_allclose(result['probabilities'][:2], [0, 0.25], atol=2e-13)
    assert len(result['counts']) == 8
    assert sum(result['counts']) == m.SHOTS
    histories = {tuple(row['code_acceptance']) for row in result['labels']}
    assert histories == {(True, True), (True, False), (False, True), (False, False)}


def test_deposited_full_turn_differs_from_selected_order_sign():
    data = json.loads((m.HERE / 'results.json').read_text())
    for source in data['records']:
        for panel in source['panels'].values():
            assert sum(panel['counts']) == m.SHOTS
            assert abs(sum(panel['probabilities']) - 1) < 2e-13
        p = source['panels']
        np.testing.assert_allclose(p['filtered_order']['probabilities'][:2], [0, 0.25], atol=2e-13)
        np.testing.assert_allclose(p['filtered_full_turn']['probabilities'][:2], [0.25, 0], atol=2e-13)
        np.testing.assert_allclose(p['unfiltered_order']['probabilities'], [1, 0], atol=2e-13)
    even, changed = data['changed_law_controls']
    assert max(even['calibration']['spatial_pi_calibration_errors']) < 1e-12
    assert max(changed['calibration']['spatial_pi_calibration_errors']) > 0.1
    for control in data['changed_law_controls']:
        for source in control['records']:
            for key, panel in source.items():
                if isinstance(panel, dict):
                    assert sum(panel['counts']) == m.SHOTS


def test_optimized_receipt_checker_fails_on_a_changed_count():
    code = f"from {m.__package__} import measure; measure.compare({{'count':1}},{{'count':2}})"
    result = subprocess.run([sys.executable, '-O', '-B', '-c', code], capture_output=True, text=True)
    assert result.returncode != 0 and 'ArithmeticError' in result.stderr


def test_actual_optimized_replay_preserves_the_archive():
    import hashlib
    path = m.HERE / 'results.json'
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    result = subprocess.run([sys.executable, '-O', '-B', '-W', 'error', '-m',
                             m.__package__ + '.measure'], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert 'Replayed complete projected native half-turn' in result.stdout
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
