import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.law import seam_term
from .guide import Flight, physical_rings, relative_hamiltonian, signed_pulse_audit

HERE = Path(__file__).resolve().parent


def test_unreduced_physical_seams_and_orientation_lines():
    size = 8
    words = [(x, y) for x in range(size) for y in range(size) if x != y]
    h = np.zeros((len(words), len(words)))
    local = seam_term(8)
    for col, (x, y) in enumerate(words):
        loads = [1] * size
        loads[x], loads[y] = 0, 2
        for i in range(size):
            j = (i + 1) % size
            n, m = loads[i], loads[j]
            h[col, col] += local[3 * n + m, 3 * n + m]
            if n != m:
                moved = loads.copy()
                moved[i], moved[j] = m, n
                row = words.index((moved.index(0), moved.index(2)))
                h[row, col] += local[3 * m + n, 3 * n + m]
    for k in range(size):
        basis = np.zeros((len(words), size - 1), complex)
        for row, (x, y) in enumerate(words):
            basis[row, (y - x) % size - 1] = np.exp(2j * np.pi * k * x / size) / np.sqrt(size)
        np.testing.assert_allclose(h @ basis, basis @ relative_hamiltonian(size, k, 8), atol=1e-12)
    for sign, energy in [(1, 49), (-1, 50)]:
        rest = Flight(size, size // 2, 8, sign)
        np.testing.assert_allclose(rest.h @ rest.initial, energy * rest.initial, atol=1e-12)


def test_native_impulse_on_every_ring_and_signed_controls():
    _, rings, _ = physical_rings()
    data = json.loads((HERE / 'results.json').read_text())
    for size, ring in rings.items():
        assert len(ring['ports']) == size
        checks = signed_pulse_audit(size)
        assert max(row['max_column_entry_error'] for row in checks.values()) < 1e-11
    for panel in data['impulses']:
        delta = {'pair': (0, 2), 'vacancy': (-1, 1), 'excess': (1, 1)}[panel['stock']]
        for mode, shift in zip(('inactive', 'signed', 'unsigned'), (0, *delta)):
            counts = panel['terminal_character_counts'][mode]
            assert sum(counts) == counts[(panel['incoming'] + shift) % panel['size']] == 2**20


def test_held_reader_uses_calibration_counts_and_position_histograms():
    data = json.loads((HERE / 'results.json').read_text())
    native_errors, changed_errors = [], []
    for row in data['positions']:
        clock = data['calibration'][str(row['orientation'])]['displacement_quadratures']
        z = complex(*[(p - m) / (p + m) for m, p in clock])
        inverse_mass = -np.angle(z) / (256 * (1 - np.cos(2 * np.pi / 20)))
        q = 2 * np.pi * row['offset'] / row['size']
        for line in row['steps']:
            counts = np.asarray(line['vacancy_position_counts'])
            assert sum(counts) == 2**20
            signal = np.sum(np.exp(-1j * q * np.arange(len(counts))) * counts)
            if line['time'] == 0:
                baseline = signal
            error = abs(signal / baseline - np.exp(-1j * line['time'] * (1 - np.cos(q)) * inverse_mass))
            np.testing.assert_allclose(error, line['mass_discrepancy'], atol=2e-14)
            (native_errors if row['contrast'] == 8 else changed_errors).append(error)
    assert max(native_errors) < .07
    assert max(changed_errors) > .15


def test_unrestricted_release_counts_and_independent_receivers():
    release = json.loads((HERE / 'results.json').read_text())['full_release']
    filtered = json.loads((HERE / 'free_results.json').read_text())
    for panel in release['sources']:
        final = panel['exposures'][-1]
        assert sum(final['joint_counts']) == sum(final['ring_character_and_outside_counts']) == 2**20
        assert final['ring_character_probability'] < .99
    for panel in filtered['mounted_sources']:
        assert panel.get('unresolved') or panel['admission'] > 0
        if 'unresolved' in panel:
            continue
        for record in panel['signed_remote']['records'].values():
            assert sum(record['counts']) == 2**58
        for row in panel['response']['native']:
            for record in row['records']:
                assert sum(record['counts']) == 2**58


def test_free_clock_aliases_current_and_refusals_from_raw_counts():
    data = json.loads((HERE / 'rest_results.json').read_text())
    assert data['absolute_energy_upper_bound'] / 128 < np.pi
    assert max(data['all_time_doublet_bridge_current_bounds']) < 1e-5
    for label, i in [('right', 0), ('left', 1), ('host', 2)]:
        assert data['bridge_calibration_counts'][label][i] == 2**58
    for panel in data['sources']:
        coarse = None
        for clock in panel['energies']:
            t = clock['tau']
            for counts in clock['quadratures']:
                assert sum(counts) == 2**58
            z = complex(*[(b - a) / (a + b) for a, b, _ in clock['quadratures']])
            angle = -np.angle(z)
            if coarse is None:
                coarse = angle / t
            angle += 2 * np.pi * round((coarse * t - angle) / (2 * np.pi))
            np.testing.assert_allclose(angle / t, clock['measured_cost'])
            assert 132.11 < angle / t < 132.12
        for currents in panel['currents']:
            assert len(currents['bridge_counts']) == 30
            assert all(sum(counts) == 2**58 for counts in currents['bridge_counts'])
            assert currents['largest_abs_measured_current'] < 1e-5
