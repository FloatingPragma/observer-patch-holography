"""Independent graph sum rule, full reference circuit and raw-record inference."""

from copy import deepcopy
import json
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import measure as m, readback


def receipt():
    return json.loads((m.HERE / 'results.json').read_text())


def test_scalar_sum_rule_with_arbitrary_inhomogeneity_and_all_lines():
    # A connected nonregular five-site graph, unrelated to the native solver.
    edges = [(0, 1), (1, 2), (1, 3), (3, 4)]
    incidence = np.zeros((4, 5))
    for i, (a, b) in enumerate(edges):
        incidence[i, a], incidence[i, b] = 1, -1
    lap = incidence.T @ incidence
    psi = np.array([1., 2., 3., 2., 4.])
    psi /= np.linalg.norm(psi)
    A = .37
    potential = -A * (lap @ psi) / psi
    h = A * lap + np.diag(potential)
    energies, basis = np.linalg.eigh(h)
    assert abs(energies[0]) < 1e-14
    for f in (np.array([1., -1., .3, .2, -.1]), np.array([.3, .7, -.4, 1., 0.])):
        spectral = np.sum((energies - energies[0]) * abs(basis.T @ (f * psi))**2)
        rhs = A * sum(psi[a] * psi[b] * (f[a] - f[b])**2 for a, b in edges)
        np.testing.assert_allclose(spectral, rhs, rtol=2e-13)
        assert np.count_nonzero(abs(basis.T @ (f * psi)) > 1e-5) > 3


def test_fringe_is_full_reference_interference_with_retained_loss():
    h = np.array([[0., -.1, .03], [-.1, .2, .06], [.03, .06, .4]])
    energy, vectors = np.linalg.eigh(h)
    f = np.diag([.8, -.3, .5])
    k = vectors.T @ f @ vectors[:, :1]
    g = np.ones(1)
    tau = 2.7
    plus, minus, operator = m.fringe(energy, k, g, tau)
    # Independent direct three-site controlled-arm construction.
    source = vectors[:, 0]
    u = expm(-1j * (h - energy[0] * np.eye(3)) * tau)
    a, b = u @ f @ source, f @ u @ source
    reference = np.r_[a, b] / np.sqrt(2)
    # Y output |0>=(arm0-i arm1)/sqrt(2), |1>=(arm0+i arm1)/sqrt(2).
    analyzer = np.kron(np.array([[1, -1j], [1, 1j]]) / np.sqrt(2), np.eye(3))
    final = analyzer @ reference
    np.testing.assert_allclose([plus, minus], [np.linalg.norm(final[:3])**2, np.linalg.norm(final[3:])**2], atol=1e-14)
    assert 0 < 1 - plus - minus < 1
    np.testing.assert_allclose((plus - minus) / tau, operator[0, 0].real, atol=1e-14)
    moment = np.sum((energy - energy[0]) * abs(k[:, 0])**2)
    sinc = np.sinc(np.ptp(energy) * tau / np.pi)
    assert (plus - minus) / tau <= moment <= (plus - minus) / tau / sinc


def test_raw_counts_reproduce_all_inference_and_changed_controls():
    data = receipt()
    readback.clock(data)
    for r in data['clock_standard']['records']:
        assert sum(r['counts']) == m.SHOTS
        np.testing.assert_allclose(r['raw_probabilities'][1], np.sin(r['commanded_appointment'] / 2)**2, atol=2e-15)
        assert sum(r['compensated_input']) == 3
    readback.close(readback.infer(data, 'native_sources'), data['native_inference'])
    readback.close(readback.infer(data, 'scalar_controls'), data['scalar_inference'])
    for source in data['native_sources'] + data['scalar_controls']:
        assert len(source['population']['counts']) == 21
        assert sum(source['population']['counts']) == m.SHOTS
        for law in source['response'].values():
            for panel in law:
                for r in panel['records']:
                    assert len(r['counts']) == 5 and sum(r['counts']) == m.SHOTS
                    assert min(r['raw_probabilities']) >= 0
        null = source['which_arm_control']['first_moment_interval']
        assert null[0] == 0 and null[1] < 1e-8
    native = data['native_inference']['held']
    assert sum(r['outcome'] == 'reject' and r['shell'] == 3 and r['center'] == 0
               and r['law'] == 'native' for r in native) == 3
    assert all(r['outcome'] != 'reject' for r in native if r['law'] == 'native' and r['shell'] == 1)
    assert all(r['outcome'] == 'reject' for r in native if r['law'] == 'half_rate')
    assert all(r['outcome'] != 'reject' for r in data['scalar_inference']['held'] if r['law'] == 'native')


def test_falsifiable_held_signal_and_optimized_receipt_rejection():
    data = receipt()
    changed = deepcopy(data)
    # Alter a held fringe by a resolvable amount at fixed total, without changing metadata.
    r = changed['native_sources'][0]['response']['native'][4]['records'][0]
    transfer = sum(r['counts'][:4]) // 100
    r['counts'][0] += transfer
    r['counts'][1] -= transfer
    with pytest.raises(ArithmeticError):
        readback.close(readback.infer(changed, 'native_sources'), data['native_inference'])
    code = f'from {m.__package__} import measure; measure.compare({{"a":0}},{{"a":1}})'
    result = subprocess.run([sys.executable, '-O', '-B', '-c', code], capture_output=True, text=True)
    assert result.returncode != 0 and 'ArithmeticError' in result.stderr
    with pytest.raises(ArithmeticError):
        m.sample([-1e-14, 1 + 1e-14], np.random.default_rng(0))
