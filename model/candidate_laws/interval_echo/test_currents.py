import json
from pathlib import Path

import numpy as np
import pytest

from model.executable_model import cell
from ..contrast_memory_spherical_v1.sectors import BinarySector
from .currents import SOURCE_PORT, STRENGTH, _geometry, pulse, run


def test_current_menu_and_control_replay():
    result = run()
    assert result == json.loads(Path(__file__).with_name('currents.json').read_text())
    for row in (result['one_load'], result['two_load']):
        assert row['guarded']['rank'] == row['guarded']['closure'] == row['guarded']['visible_in_records'] == 12
        assert row['with_bare_port_zero']['rank'] == 13
        assert row['with_bare_port_zero']['closure'] == 18
        assert row['with_bare_port_zero']['visible_in_records'] == 17
        assert row['one_edge_rotations'] == 12 and row['port_rotation_checks'] == 720
        assert row['brackets_checked'] == 66
        assert [p['dimension'] for p in row['source_poles'] if p['retained']] == [3, 3]


def test_public_pulse_and_complete_detector():
    from scipy.linalg import expm

    car = cell()
    sector = BinarySector(car, 2)
    antipodes, _ = _geometry(car)
    source = np.zeros(len(sector.states))
    for sign, p in ((1, SOURCE_PORT), (-1, antipodes[SOURCE_PORT])):
        for q in car.neighbours[p]:
            source[sector.index[tuple(sorted((p, q)))]] = sign / np.sqrt(2 * len(car.neighbours[p]))
    h = sector.hamiltonian(8)
    forward = pulse(0, STRENGTH)
    backward = pulse(0, -STRENGTH)
    assert np.allclose(forward @ backward, np.eye(len(sector.states)), atol=1e-12, rtol=0)
    probabilities = abs(expm(-16j * h) @ forward @ source)**2
    seam = [sector.index[tuple(e)] for e in car.seams]
    reading = np.r_[probabilities[seam], probabilities[sector.agreed == 0].sum()]
    report = json.loads(Path(__file__).with_name('currents.json').read_text())
    rows = report['two_load']['detector']['rows']
    assert np.allclose(reading, rows[0]['positive'], atol=1e-11, rtol=0)
    assert 0.11 < rows[0]['total_variation'] < 0.12
    for row in rows:
        assert abs(sum(row['positive']) - 1) < 1e-9 and abs(sum(row['negative']) - 1) < 1e-9
    assert pulse(0, 0, loads=1).shape == (12, 12)
    with pytest.raises(ValueError):
        pulse(12, STRENGTH)
