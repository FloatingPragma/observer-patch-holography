"""Independent CAR/reflection oracles and complete source-detector controls."""
from fractions import Fraction as F
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import expm

from ..full_fock_field.field import FieldGeometry
from ..observer.source_scalar import compiler
from ..particle_response import dimer
from . import detector as d


def annihilator(modes, mode):
    result = np.zeros((1 << modes, 1 << modes), complex)
    for mask in range(1 << modes):
        if mask & (1 << mode):
            result[mask ^ (1 << mode), mask] = (-1)**sum(
                (mask >> earlier) & 1 for earlier in range(mode))
    return result


def independent_reflection(index, phase):
    a, r = annihilator(5, index), annihilator(5, 4)
    return np.eye(32)-a.conj().T@a-r.conj().T@r+(-1)**phase*(a.conj().T@r+r.conj().T@a)


@pytest.mark.parametrize('index', range(4))
@pytest.mark.parametrize('phase', (0, 1))
def test_reflection_matches_literal_occupation_signs_and_all_sectors(index, phase):
    actual = d.reflection(index, phase)
    assert np.array_equal(actual, independent_reflection(index, phase))
    assert np.array_equal(actual@actual, np.eye(32))
    for mask in range(32):
        nj, nr = (mask >> index) & 1, (mask >> 4) & 1
        if nj == nr:
            assert actual[mask, mask] == (-1 if nj else 1)
        else:
            other = mask ^ (1 << index) ^ (1 << 4)
            between = sum((mask >> j) & 1 for j in range(index+1, 4))
            assert actual[other, mask] == (-1)**(between+phase)
    with pytest.raises(ValueError):
        actual[0, 0] = 0


def test_full_operator_and_actual_source_word_admission():
    report = d.operator_certificate()
    assert report['full_controlled_input_columns_checked'] == 512
    assert report['complete_field_operator_identities_checked'] == 48
    assert report['arbitrary_untouched_reference_extension']
    assert not report['primitive_pulse_CAR_number_conservation_claimed']


@pytest.mark.parametrize('axis', range(3))
@pytest.mark.parametrize('reservoir', (0, 1))
def test_nonstationary_entangled_number_two_input_obeys_operator_green_identity(axis, reservoir):
    sides = [1, 1, 1]; sides[axis] = 2
    spectrum = dimer.setup(FieldGeometry(tuple(sides)))
    unitary = expm(-.37j*spectrum.generator)
    evolved = np.kron(np.eye(2), unitary)
    # Three arbitrary untouched-reference columns, with coherent N=2 input.
    columns = np.zeros((16, 3), complex)
    for mask in range(16):
        if mask.bit_count() == 2:
            columns[mask] = [(mask+1)*(j+1)+1j*(mask-j*j) for j in range(3)]
    columns /= np.linalg.norm(columns)
    joint = np.zeros((32, 3), complex)
    joint[16*reservoir:16*(reservoir+1)] = columns
    for index in range(4):
        a = annihilator(4, index)
        at = unitary.conj().T@a@unitary
        correlations = []
        for phase in (0, 1):
            out = evolved.conj().T@independent_reflection(index, phase)@evolved@independent_reflection(index, 0)@joint
            correlations.append(complex(np.vdot(joint, out)))
        expected = (at@a.conj().T) if reservoir else at.conj().T@a
        assert abs((correlations[0]-correlations[1])/2-np.vdot(columns, expected@columns)) < 2e-13


@pytest.fixture(scope='module')
def complete_readout():
    return d.readout()


def test_same_candidate_green_functions_from_all_complete_source_settings(complete_readout):
    row = complete_readout
    assert row['independent_preparations_required'] == 8
    assert len(row['complete_records']) == 8
    assert {(x['reservoir_occupation'], x['reflection_phase_in_pi'], x['pointer_quadrature'])
            for x in row['complete_records']} == {(r, p, q) for r in (0, 1) for p in (0, 1) for q in (0, 1)}
    for name in ('forward', 'reverse'):
        assert abs(complex(*row['observed_'+name])-complex(*row['complete_Lehmann_'+name])) < 2e-6
    assert row['maximum_complete_vector_error'] < 1e-5
    assert row['maximum_work_bits'] == 6
    assert not row['controlled_odd_field_operators_used']
    assert row['retained_purifier_bits'] == 0
    assert not row['conditioned_on_success_or_particle_outcome']
    totals = {key: 0 for key in d.COUNT_KEYS}
    for record in row['complete_records']:
        assert len(record['pointer_probabilities']) == 2
        assert abs(sum(record['pointer_probabilities'])+record['preparation_abort_probability']
                   +record['continuation_abort_probability']-1) < 2e-9
        assert record['first_export_event_count'] == record['resources']['native_pulses']
        assert record['preparation_export_event_count'] == record['preparation_resources']['native_pulses']
        assert record['retained_continuation_complement_count'] > record['first_export_event_count']
        for key in totals:
            totals[key] += record['resources'][key]+record['preparation_resources'][key]
    assert totals == row['total_resources_including_preparations']


def test_zero_time_requires_subtraction_not_mislabeling_the_fswaP_correlation():
    spectrum = dimer.setup()
    direct = d.correlator(spectrum, duration=F(0), steps=1)
    opposite = d.correlator(spectrum, phase=1, duration=F(0), steps=1)
    assert abs(direct['measured_mean']-1) < 2e-10
    assert abs(opposite['measured_mean']) < 2e-10
    assert abs((direct['measured_mean']-opposite['measured_mean'])/2-.5) < 2e-10
    # Keeping the uncycled reflection alone would wrongly report G(0)=1.
    assert abs(direct['measured_mean']-.5) > .49


def test_pointer_dephasing_erases_signal_without_erasing_records():
    spectrum = dimer.setup()
    ground = spectrum.vectors[:, 0]
    initial = np.zeros(32, complex); initial[16:] = ground
    acted = independent_reflection(0, 0)@initial
    # Retain both dephasing outcomes. Their incoherent mixture has no path
    # off-diagonal, whereas summing vectors would invent interference.
    a = np.concatenate((initial, np.zeros(32)))/np.sqrt(2)
    b = np.concatenate((np.zeros(32), acted))/np.sqrt(2)
    x = np.kron(np.array([[0, 1], [1, 0]]), np.eye(32))
    assert np.vdot(a, x@a)+np.vdot(b, x@b) == 0
    assert abs(np.vdot(a, a)+np.vdot(b, b)-1) < 2e-13


@pytest.mark.parametrize('damage', ('drop_scalar', 'omit_jw_string'))
def test_complete_source_gate_audit_rejects_phase_or_routing_corruption(monkeypatch, damage):
    original = d.controlled_reflection_word
    def changed(index=0, phase=0):
        word = original(3 if damage == 'omit_jw_string' else index, phase)
        if damage == 'drop_scalar':
            word = [event for event in word if event['kind'] != 'I']
        return word
    monkeypatch.setattr(d, 'controlled_reflection_word', changed)
    with pytest.raises(ArithmeticError, match='controlled reflection'):
        d.operator_certificate()


@pytest.mark.parametrize('damage', ('drop_event', 'drop_complement'))
def test_zero_export_complements_cannot_disappear(monkeypatch, damage):
    execute = d.execute_with_exports
    def changed(word, bits, initial):
        state, failures = execute(word, bits, initial)
        failures = list(failures)
        if damage == 'drop_event':
            failures.pop()
        else:
            failures[0] = failures[0][:-1]
        return state, tuple(failures)
    monkeypatch.setattr(d, 'execute_with_exports', changed)
    with pytest.raises(ArithmeticError, match='complete'):
        d.correlator(dimer.setup(), duration=F(0), steps=1)


def test_actual_export_leakage_is_rejected(monkeypatch):
    native = compiler.native_pulse
    def leaking(kind, angle):
        full = native(kind, angle).copy()
        rotation = np.eye(6, dtype=complex)
        angle = .001
        rotation[0, 0] = rotation[5, 5] = np.cos(angle)
        rotation[0, 5] = rotation[5, 0] = 1j*np.sin(angle)
        return rotation@full
    monkeypatch.setattr(compiler, 'native_pulse', leaking)
    with pytest.raises(ArithmeticError):
        d.correlator(dimer.setup(), duration=F(0), steps=1)


def test_preparation_and_counters_are_not_replaceable(monkeypatch):
    spectrum = dimer.setup()
    with pytest.raises(TypeError, match='prepared'):
        d.correlator(spectrum, prepared=np.eye(16)[0])
    prepare = d.source.prepare
    def wrong(*args, **kwargs):
        original = prepare(*args, **kwargs)
        attributes = {key: getattr(original, key) for key in
                      ('state', 'target', 'failures', 'resources', 'initial', 'word')}
        attributes['resources'] = dict(attributes['resources'])
        attributes['resources']['native_pulses'] = 0
        return SimpleNamespace(**attributes)
    monkeypatch.setattr(d.source, 'prepare', wrong)
    with pytest.raises(ArithmeticError, match='resource account'):
        d.correlator(spectrum, duration=F(0), steps=1)


@pytest.mark.parametrize('bad', (True, 0., -1, 4))
def test_strict_mode_domain_before_any_execution(bad):
    with pytest.raises(ValueError): d.reflection(bad)
    with pytest.raises(ValueError): d.controlled_reflection_word(bad)


@pytest.mark.parametrize('key,value', (('reservoir', True), ('phase', .0),
    ('quadrature', 2), ('duration', .25), ('duration', F(2)), ('steps', False), ('steps', 0)))
def test_detector_domain_guards(key, value):
    with pytest.raises(ValueError):
        d.correlator(dimer.setup(), **{key: value})
