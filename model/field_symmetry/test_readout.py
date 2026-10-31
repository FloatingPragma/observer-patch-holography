"""Independent Hamiltonian, tensor-output and hostile source controls."""
from fractions import Fraction as F
import math
import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.full_fock_field.field import FullFock
from model.field_symmetry import readout as r


def test_preparation_macro_preserves_parity_on_every_input():
    word = r.preparation(r.geometry())
    matrix = r.compiler.execute(word, 8)
    parity = np.array([(-1)**i.bit_count() for i in range(256)])
    assert np.linalg.norm(parity[:, None]*matrix-matrix*parity[None, :]) < 2e-11


@pytest.mark.parametrize('plane', [0, 1, 2])
@pytest.mark.parametrize('centered', [False, True])
def test_literal_complete_source_against_full_hamiltonian(plane, centered):
    g = r.geometry(plane)
    out = r.execute(plane, centered=centered)
    space = FullFock(8)
    h = space.hamiltonian(g)
    if centered:
        h -= sum(space.number(i) for i in range(8))/2
    # Construct preparation/analyzer independently as tensor matrices.
    encoder = np.ones((1, 1), complex)
    for point in reversed(g.sites):
        sign = (-1)**sum(point)
        local = np.array([[1, 0, 0, -sign], [0, 1, sign, 0],
                          [0, -sign, 1, 0], [sign, 0, 0, 1]], complex)/math.sqrt(2)
        encoder = np.kron(encoder, local)
    blank = np.zeros(256); blank[0] = 1
    expected = encoder.conj().T@expm_multiply(-1j*h, encoder@blank)
    np.testing.assert_allclose(out['output'], expected, atol=2e-12)
    assert out['report']['first_failure_event_slots'] == out['report']['resources']['native_pulses']
    assert len(out['output']) == 256
    expected_bit = 0 if centered else math.sin(.5)**2
    np.testing.assert_allclose(out['report']['read_one_probabilities'], [expected_bit]*4, atol=2e-12)
    assert out['report']['state_specific_ideal_split_error'] == '0'
    assert F(out['report']['arbitrary_input_ideal_split_bound']) > F(1, 2)


@pytest.mark.parametrize('duration,steps', [(F(0), 1), (F(-1), 1), (F(1), 2), (F(1, 3), 1)])
def test_zero_negative_time_and_split_refinement(duration, steps):
    a = r.execute(duration=duration, steps=steps)
    assert a['report']['complete_output_error'] < 5e-12
    assert r.independent_evolution(duration=duration) < 2e-14


def test_free_control_and_exact_signal_enclosure():
    free = r.execute(interaction=0.)
    assert max(free['report']['read_one_probabilities']) < 1e-20
    assert F(529, 2304) < math.sin(.5)**2 < F(1, 4)


@pytest.mark.parametrize('change', ['drop_event', 'drop_complement', 'nan', 'wrong_phase', 'wrong_address'])
def test_complete_output_gate_refuses_corruption(monkeypatch, change):
    original = r.execute_with_exports
    def corrupt(*args):
        output, failures = original(*args)
        if change == 'drop_event': failures = failures[:-1]
        elif change == 'drop_complement': failures = (failures[0][:-1],)+failures[1:]
        elif change == 'nan':
            row = failures[0][0].copy(); row[0] = float('nan')
            failures = ((row,)+failures[0][1:],)+failures[1:]
        elif change == 'wrong_phase':
            output = output.copy(); output[3] *= -1
        else: output = np.roll(output, 2)
        return output, failures
    monkeypatch.setattr(r, 'execute_with_exports', corrupt)
    with pytest.raises(ArithmeticError): r.execute()


@pytest.mark.parametrize('change', ['nan', 'wrong_unitary', 'leak'])
def test_primitive_gate_refuses_wrong_complete_action(monkeypatch, change):
    original = r.compiler.native_pulse
    def corrupt(kind, angle):
        u = original(kind, angle)
        if change == 'nan': u[0, 0] = float('nan')
        elif change == 'wrong_unitary': u[0] *= -1
        else:
            u[[0, 5]] = u[[5, 0]]
        return u
    monkeypatch.setattr(r.compiler, 'native_pulse', corrupt)
    with pytest.raises(ArithmeticError): r.execute()


@pytest.mark.parametrize('kwargs', [{'plane': True}, {'centered': 1}, {'duration': True},
                                   {'duration': 1.}, {'steps': True}, {'steps': 0},
                                   {'interaction': float('nan')}])
def test_strict_public_domain(kwargs):
    with pytest.raises(ValueError): r.execute(**kwargs)


def test_occupation_read_alone_does_not_read_eta_phase():
    g = r.geometry()
    np.testing.assert_allclose(abs(r.tower_state(g, .7))**2, abs(r.tower_state(g))**2, atol=1e-16)
    # A local wrong staggered sign is not a globally protected eta state.
    state = r.tower_state(g).copy()
    state[(np.arange(256) & 3) == 3] *= -1
    space = FullFock(8)
    kinetic = space.hamiltonian(r.geometry(interaction=0.))
    assert np.linalg.norm(kinetic@state) > .01


def test_complete_source_retains_small_nonzero_export_failures(monkeypatch):
    original = r.compiler.native_pulse
    angle = 3e-12
    rotation = np.eye(6, dtype=complex)
    rotation[0, 0] = rotation[5, 5] = math.cos(angle)
    rotation[0, 5], rotation[5, 0] = -math.sin(angle), math.sin(angle)
    monkeypatch.setattr(r.compiler, 'native_pulse', lambda kind, value: rotation@original(kind, value))
    out = r.execute()
    assert out['report']['first_export_failure_probability'] > 1e-22
    assert sum(len(event) for event in out['failures']) == 8160
    assert out['report']['complete_output_error'] < float(F(out['report']['numerical_dilation_budget']))
