"""Cross-check the detector probabilities against the original tour and cut identities."""

import json
import math
from pathlib import Path

import pytest

from model.candidate_laws.interval_echo.openwave.quark_winding import measure, phase_control


HERE = Path(__file__).parent


@pytest.fixture(scope='module')
def parity():
    archived = json.loads((HERE / 'results.json').read_text())
    assert json.loads(json.dumps(measure.run())) == archived
    return archived


@pytest.fixture(scope='module')
def third():
    archived = json.loads((HERE / 'third_turn.json').read_text())
    assert json.loads(json.dumps(phase_control.run())) == archived
    return archived


def test_full_return_records(parity):
    parent = json.loads((HERE.parent / 'electron_field/flux_string/results.json').read_text())
    primary, held = parity['settings']
    assert primary['faces'] == [0, parent['opposite_face']]
    assert held['faces'] == [1, 8]
    previous = next(t for t in parent['scan']['trials'] if t['ends'] == primary['faces'])
    size = parity['field_words']
    for setting in (primary, held):
        assert setting['cut_pair_multiplicities'] == [{'cut': [0, 0], 'words': size}]
        assert {f for f, e in enumerate(setting['ring_energy_by_face']) if e > 0} == set(setting['faces'])
        for row in setting['panels']:
            counts = row['path_counts']
            assert sum(counts.values()) == size
            inside = len(set(setting['faces']).intersection(row['disk']))
            sign = (-1) ** inside
            expected = (sign * counts['winding'] + counts['idle']) / size
            assert row['string']['inclusive_bright'] == pytest.approx((1 + expected) / 2)
            assert row['fully_dephased_bright_both_inputs'] == (1 + counts['idle'] / size) / 2
            for prep in ('bare', 'string'):
                bins = row[prep]['by_cuts_and_final_port']
                assert sum(b['bright'] + b['dark'] for b in bins) == 4 * size
                assert all(b['cut'] == [0, 0] for b in bins)
                voltage = row['pre_tour_material_voltage_by_region'][prep]
                assert all(axis == {'bright_x': 1., 'bright_y': 0.5}
                           for region in voltage for axis in region.values())
            if setting is primary and row['contour'] == 'enclosing':
                record = previous['tour'][str(row['sign'])]
                assert row['string']['inclusive_bright'] == pytest.approx(record['inclusive_bright_port'])
            if row['contour'] == 'deformed':
                assert counts['winding'] == (0 if setting is primary else 18)


@pytest.mark.parametrize('end', ['near', 'far'])
def test_charge_controls(parity, end):
    for setting in parity['settings']:
        for control in setting['calibrations_by_region'][end]:
            sign, fraction = control['sign'], control['tips_inside'] / 3
            written, unwritten = (control['reads'][k] for k in ('written', 'unwritten'))
            assert written['sharp']['joint'] == [{'material': sign, 'cut': sign, 'probability': 1}]
            assert unwritten['sharp']['joint'] == [{'material': sign, 'cut': 0, 'probability': 1}]
            probabilities = {(b['material'], b['cut']): b['probability'] for b in written['spread']['joint']}
            assert probabilities == pytest.approx({(0, 0): 1 - fraction, (sign, sign): fraction})
            assert unwritten['spread']['mean_cut'] == 0
            voltage = written['spread']['material_voltage']
            assert voltage['pi']['bright_x'] == pytest.approx(1 - fraction)
            assert voltage['pi_over_two']['bright_y'] == pytest.approx((1 + sign * fraction) / 2)
            assert voltage != written['spread']['hypothetical_sharp_charge_at_same_mean']


def test_oriented_phase_from_receiver_ports(parity, third):
    assert third['singleton_cut_values'] == [[0]] * 12
    for binary, oriented in zip(parity['settings'], third['settings'], strict=True):
        assert binary['faces'] == oriented['faces']
        assert oriented['energy_by_face'] == pytest.approx([0.75 * e for e in binary['ring_energy_by_face']])
        for base, panel in zip(binary['panels'], oriented['panels'], strict=True):
            assert (base['contour'], base['sign']) == (panel['contour'], panel['sign'])
            for bins in panel['outcomes'].values():
                for axis in ('x', 'y'):
                    assert sum(b['bright_' + axis] + b['dark_' + axis] for b in bins) == pytest.approx(1)
            enclosed = panel['contour'] != 'exterior'
            winding = base['path_counts']['winding'] if enclosed else 0
            assert panel['conjugate_total_variation_y'] == pytest.approx(
                math.sqrt(3) / 2 * winding / parity['field_words'], abs=3e-12)
            if panel['contour'] == 'enclosing':
                # Reconstruct the turn from the measured ports, without using the input string.
                angles = []
                for prep in ('positive', 'negative'):
                    forward = next(b for b in panel['outcomes'][prep] if b['port'] == panel['cycle'][0]
                                   and b['heads'] == 'forward')
                    angles.append(math.atan2(forward['bright_y'] - forward['dark_y'],
                                             forward['bright_x'] - forward['dark_x']))
                assert sorted(angles) == pytest.approx([-2 * math.pi / 3, 2 * math.pi / 3])
