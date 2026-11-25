"""Check the physical coupling, energy accounting and complete detector records."""

from itertools import product
import json

import numpy as np
from scipy.linalg import expm

from . import measure as m


def test_local_generator_and_full_source_blocks_from_port_loads():
    for k, cut in m.parent.ARMS:
        full = m.walls.operators(k, cut, 3).reshape(16, 64, 16, 64)
        _, _, ordered = m.parent.reference()
        boundary = m.parent.reference()[0]
        for phase, index in enumerate(ordered):
            hs = []
            for b in boundary:
                b = int(b if phase == 0 else 2 - b)
                words = list(product((0, 1), repeat=3))
                independent = np.zeros((8, 8))
                cost = 4 * k + .5
                for j, bits in enumerate(words):
                    loads = [x for bit in bits for x in (2 * bit, 2 - 2 * bit)]
                    independent[j, j] = cost * (3 + (not cut) * (b != loads[0])
                        + (loads[1] != loads[2]) + (loads[3] != loads[4]))
                    for q in range(3):
                        changed = tuple(bit ^ (r == q) for r, bit in enumerate(bits))
                        independent[words.index(changed), j] = -.5
                h = m.walls.wing_ref(k, b, cut, 3)[0]
                np.testing.assert_allclose(h, independent, atol=1e-13)
                hs.append(h)
            sum_wings = np.kron(hs[0], np.eye(8)) + np.kron(np.eye(8), hs[1])
            residual = full[index, :, index, :] - sum_wings
            np.testing.assert_allclose(residual, residual[0, 0] * np.eye(64), atol=1e-12)


def test_detector_is_static_local_and_conserves_coupled_energy():
    h = m.ex.wing(8, 3, 0)
    levels, vectors = np.linalg.eigh(h)
    omega = levels[2] - levels[0]
    meter = m.receiver(h, omega, 32)
    x = np.array([[0., 1.], [1., 0.]])
    direct = (np.kron(h - levels[0] * np.eye(8), np.eye(2))
              + np.kron(np.eye(8), np.diag([0., omega]))
              + m.COUPLING * np.kron(np.kron(np.eye(4), x), x))
    np.testing.assert_allclose(meter['generator'], direct, atol=1e-12)
    np.testing.assert_allclose(meter['unitary'], expm(-32j * direct), atol=4e-12)
    state = np.outer(vectors[:, 2], vectors[:, 2])
    ledger = m.ledger(state, [meter, meter])
    assert ledger[1] > 1.9 * omega
    assert abs(ledger[0] - ledger[1] - ledger[2]) < 1e-10
    # At zero receiver coupling the same hard source transfers no energy.
    initial = np.kron(vectors[:, 2], [1., 0.])
    decoupled = direct - meter['operators'][2]
    final = expm(-32j * decoupled) @ initial
    assert np.sum(abs(final.reshape(8, 2)[:, 1])**2) < 1e-25


def test_complete_cohorts_controls_and_energy_budget():
    data = json.loads((m.HERE / 'results.json').read_text())
    for arm in data['arms']:
        for panel in arm['readings']:
            sham, walls = panel['cohorts']
            assert abs(walls['initial_total_energy'] - sham['initial_total_energy']
                       - (8 * arm['contrast'] + 1)) < 1e-9
            for cohort in panel['cohorts']:
                assert abs(cohort['before_ring_read_energy'] - cohort['initial_total_energy']) < 2e-7
                for record in cohort['records']:
                    assert len(record['counts']) == 9
                    assert sum(record['counts']) == data['shots']
                    assert min(record['probabilities']) >= 0
                    assert abs(sum(record['probabilities']) - 1) < 3e-10
                    e = record['energy_per_original_launch']
                    assert abs(e['wing_loss'] - e['detector_gain'] - e['interaction_final']) < 2e-9
                    assert abs(record['switch_off_work_per_original_launch'] + e['interaction_final']) < 1e-12
            if arm['attachments_cut']:
                for s, w in zip(sham['records'], walls['records'], strict=True):
                    ps = np.array(s['probabilities'][:-1]).reshape(2, 4).sum(axis=0) / sham['source_order_probability']
                    pw = np.array(w['probabilities'][:-1]).reshape(2, 4).sum(axis=0) / walls['source_order_probability']
                    np.testing.assert_allclose(ps, pw, atol=1e-10)
    native = data['arms'][0]['readings'][1]['cohorts']
    sham, wall = [c['records'][2] for c in native]
    both = lambda r: r['probabilities'][3] + r['probabilities'][7]
    assert both(wall) > .0005 and both(sham) < 1e-7
    assert data['native_late_double_click_excess'] > data['simultaneous_99_difference_halfwidth']
    assert wall['detector_gain_per_ordered_source'] > 32
    changed = data['arms'][2]['readings'][1]['cohorts'][1]['records'][2]
    assert both(changed) < 1e-9


def test_independent_detector_count_energy_and_input_calibration():
    data = json.loads((m.HERE / 'results.json').read_text())
    for arm in data['arms']:
        for panel in arm['readings']:
            for cohort in panel['cohorts']:
                for record in cohort['records']:
                    p = np.array(record['probabilities'][:-1]).reshape(2, 4).sum(axis=0)
                    detector_energy = data['target_frequency'] * (p[1] + p[2] + 2 * p[3])
                    assert abs(detector_energy - record['energy_per_original_launch']['detector_gain']) < 1e-10
    for calibration in data['calibration'][:2]:
        response = calibration['settings'][2]['click_probabilities_by_level']
        assert response[0] < 2e-6 and response[2] > .97
