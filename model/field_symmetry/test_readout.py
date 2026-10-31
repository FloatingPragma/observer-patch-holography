"""Independent Hamiltonian, tensor-output and hostile source controls."""
from fractions import Fraction as F
import math
import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.full_fock_field.field import FullFock
from model.field_symmetry import readout as r


def test_preparation_macro_preserves_parity_on_every_input():
    word = r.preparation(r.geometry())
    matrix = r.compiler.execute(word, 8)
    parity = np.array([(-1)**i.bit_count() for i in range(256)])
    assert np.linalg.norm(parity[:, None]*matrix-matrix*parity[None, :]) < 2e-11


@pytest.mark.parametrize('centered', [False, True])
def test_complete_strang_target_against_independent_factor_exponentials(centered):
    g = r.geometry(plane=1, interaction=.3)
    duration, steps = F(-2, 3), 2
    space = FullFock(8)
    a = [space.annihilator(mode) for mode in range(8)]
    hops = []
    for group in g.groups:
        for p, q, z in group:
            directed = z*a[p].conj().T@a[q]
            generator = (directed+directed.conj().T).toarray()
            hops.append(expm(-1j*float(duration)/(2*steps)*generator))
    contact = space.hamiltonian(r.field.FieldGeometry(g.sides, speed=0., interaction=g.interaction))
    middle = expm(-1j*float(duration)/steps*contact.toarray())
    target = np.eye(256, dtype=complex)
    for _ in range(steps):
        for factor in hops+[middle]+list(reversed(hops)):
            target = factor@target
    if centered:
        number = sum(space.number(mode) for mode in range(8))
        target = expm(.5j*g.interaction*float(duration)*number.toarray())@target
    actual = r._ideal_strang(g, duration, steps, centered)
    assert np.linalg.norm(actual-target, 2) < 2e-14


@pytest.mark.parametrize('change', ['empty', 'hopping_free', 'reversed', 'dark_phase'])
def test_full_word_binding_refuses_wrong_action_hidden_by_eta_state(monkeypatch, change):
    original = r.field.source_word

    def corrupt(g, duration, steps):
        if change == 'empty':
            return ()
        if change == 'hopping_free':
            return original(r.field.FieldGeometry(g.sides, speed=0., interaction=g.interaction), duration, steps)
        word = original(g, duration, steps)
        if change == 'reversed':
            return tuple(reversed(word))
        # Eta states have no singly occupied sites: this wrong extra action
        # is invisible on the executed tower, even below the coarse H bound.
        phase = r.compiler.compile_diagonal(np.array([0., 1e-4, 0., 0.]))
        return tuple(word)+tuple(phase)

    monkeypatch.setattr(r.field, 'source_word', corrupt)
    with pytest.raises(ArithmeticError, match='declared Strang action'):
        r.execute(steps=4 if change == 'hopping_free' else 1)


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
    assert out['report']['full_field_word_input_columns'] == 256
    bound = F(out['report']['full_field_word_success_budget'])
    assert out['report']['full_field_word_success_error'] <= float(bound)
    assert F(out['report']['full_field_word_complement_budget']) == bound
    assert F(out['report']['full_field_word_dilation_budget']) == 2*bound


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


@pytest.mark.parametrize('numerator,denominator', [
    (True, 1), (1, True), (1, -1), (1, 0), (0, 0),
    (2, 2), (0, 2), (1.5, 1), (1, 1.5),
])
def test_exact_duration_rejects_malformed_fraction_storage(numerator, denominator):
    duration = F(1)
    object.__setattr__(duration, '_numerator', numerator)
    object.__setattr__(duration, '_denominator', denominator)
    with pytest.raises(ValueError):
        r.execute(duration=duration)
    with pytest.raises(ValueError):
        r.independent_evolution(duration=duration)


@pytest.mark.parametrize('duration', [F(0), F(-2, 3), F(7, 14)])
def test_canonical_exact_duration_remains_admitted(duration):
    result = r.execute(plane=1, duration=duration)
    expected = math.sin(float(duration)/2)**2
    np.testing.assert_allclose(result['report']['read_one_probabilities'], [expected]*4, atol=2e-12)
    assert F(result['report']['duration']) == duration


def test_success_array_subclass_cannot_hide_wrong_terminal_state(monkeypatch):
    original = r.execute_with_exports

    class BlindState(np.ndarray):
        def __sub__(self, other):
            return np.zeros(self.shape, complex)

        def __array_function__(self, function, types, args, kwargs):
            if function is np.vdot:
                return 1.+0j
            return super().__array_function__(function, types, args, kwargs)

    def corrupt(*args):
        _, failures = original(*args)
        blank = np.zeros(256, complex)
        blank[0] = 1
        return blank.view(BlindState), failures

    monkeypatch.setattr(r, 'execute_with_exports', corrupt)
    with pytest.raises(ArithmeticError):
        r.execute()


@pytest.mark.parametrize('dtype,amplitude', [(np.int64, 1 << 32), (np.float16, 1e-4)])
def test_export_probability_cannot_disappear_in_unsupported_arithmetic(monkeypatch, dtype, amplitude):
    original = r.execute_with_exports

    def corrupt(*args):
        output, failures = original(*args)
        branch = np.zeros(failures[0][0].shape, dtype=dtype)
        branch[0] = amplitude
        return output, ((branch,)+failures[0][1:],)+failures[1:]

    monkeypatch.setattr(r, 'execute_with_exports', corrupt)
    with pytest.raises(ArithmeticError):
        r.execute()


def test_small_floating_export_is_preserved_in_trusted_arithmetic(monkeypatch):
    original = r.execute_with_exports
    amplitude = np.float32(2e-12)

    def perturbed(*args):
        output, failures = original(*args)
        branch = np.zeros(failures[0][0].shape, dtype=np.float32)
        branch[0] = amplitude
        return output, ((branch,)+failures[0][1:],)+failures[1:]

    monkeypatch.setattr(r, 'execute_with_exports', perturbed)
    result = r.execute()
    branch = result['failures'][0][0]
    assert type(branch) is np.ndarray and branch.dtype == np.dtype(complex)
    assert branch[0] == float(amplitude)
    assert result['report']['first_export_failure_probability'] >= float(amplitude)**2
    assert result['report']['complete_output_error'] < float(F(result['report']['numerical_dilation_budget']))


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
