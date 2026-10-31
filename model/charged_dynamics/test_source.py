"""Independent full-input, referenced, complete-failure source controls."""
from dataclasses import replace
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import source, chain
from ..field_thermal.thermal import execute_with_exports, failure_probability
from ..full_fock_field import field
from ..observer.source_scalar import compiler


@pytest.mark.parametrize('axis', range(3))
@pytest.mark.parametrize('kind', ('hole', 'particle'))
def test_preparation_preserves_unknown_complex_spin_reference(kind, axis):
    spin = np.array([[1+2j, -3j, .5], [2j, 1-1j, -2]], complex)/math.sqrt(24.25)
    actual = source.prepare(kind, axis, spin)
    target = chain.bare_columns(kind, axis)[:, :2]@spin
    assert np.linalg.norm(actual.success-target) < 1e-10
    assert abs(np.linalg.norm(actual.success)**2+actual.failure_probability-np.linalg.norm(spin)**2) < 1e-10
    assert actual.failure_probability >= 0
    assert not actual.initial.flags.writeable
    # No postselection or normalization is hidden in this linear map.
    scaled = source.prepare(kind, axis, 2*spin)
    assert np.linalg.norm(scaled.success-2*actual.success) < 1e-12
    assert scaled.failure_probability == pytest.approx(4*actual.failure_probability, rel=1e-12, abs=0)


@pytest.mark.parametrize('axis', range(3))
def test_encoder_complete_columns_charge_and_phase(axis):
    target, word = source._encoder_data(axis)
    actual = source.encoder(axis)
    assert np.linalg.norm(actual.success-target, 2) < 1e-10
    assert np.linalg.norm(actual.success-compiler.execute(word, 4), 2) < 1e-12
    numbers = np.array([j.bit_count() for j in range(16)])
    assert np.linalg.norm((numbers[:, None]-numbers[None, :])*target, 2) < 2e-12
    inverse = source.encoder(axis, inverse=True)
    assert np.linalg.norm(inverse.success@actual.success-np.eye(16), 2) < 1e-10


@pytest.mark.parametrize('axis', range(3))
def test_one_ulp_ground_perturbations_leave_the_complete_control_word_fixed(axis, monkeypatch):
    original = chain.ground_column
    target, baseline = source._encoder_data(axis)
    for address in (3, 5, 10, 12):
        for direction in (-np.inf, np.inf):
            def changed(which_axis, address=address, direction=direction):
                value = original(which_axis).copy()
                value[address] = complex(np.nextafter(value[address].real, direction), value[address].imag)
                return value
            monkeypatch.setattr(chain, 'ground_column', changed)
            nearby, word = source._encoder_data(axis)
            assert word == baseline
            assert np.linalg.norm(nearby-target) < 2e-13
    monkeypatch.setattr(chain, 'ground_column', original)


def test_control_canonicalization_cannot_hide_a_substantive_wrong_phase(monkeypatch):
    target, word = source._encoder_data(0)
    primitive = compiler.native_pulse
    monkeypatch.setattr(compiler, 'native_pulse', lambda kind, angle: np.exp(2e-7j)*primitive(kind, angle))
    with pytest.raises(ArithmeticError, match='complete local source instrument'):
        source.compile_macro(word, 4, target, source.ENCODER_GATE)


@pytest.mark.parametrize('axis', range(3))
def test_independent_fsum_scalar_arithmetic_preserves_control_schedule(axis, monkeypatch):
    target, baseline = source._encoder_data(axis)
    old_dot, old_norm = np.dot, np.linalg.norm
    def independent_dot(a, b, *args, **kwargs):
        if np.ndim(a) == np.ndim(b) == 1 and np.isrealobj(a) and np.isrealobj(b):
            return math.fsum(float(x)*float(y) for x, y in zip(a, b, strict=True))
        return old_dot(a, b, *args, **kwargs)
    def independent_norm(a, *args, **kwargs):
        if np.ndim(a) == 1 and not args and not kwargs:
            return math.sqrt(math.fsum(abs(complex(x))**2 for x in a))
        return old_norm(a, *args, **kwargs)
    monkeypatch.setattr(np, 'dot', independent_dot)
    monkeypatch.setattr(np.linalg, 'norm', independent_norm)
    nearby, alternative = source._encoder_data(axis)
    assert alternative == baseline
    assert old_norm(nearby-target) < 2e-13


@pytest.mark.parametrize('axis', range(3))
def test_every_field_input_against_independent_existing_executor(axis):
    duration, steps = F(1, 8), 2
    actual = source.evolve(np.eye(256), axis, duration, steps)
    geometry = chain.geometry(axis)
    other = field.evolve(geometry, np.eye(256), float(duration), steps)
    assert np.linalg.norm(actual.success-other, 2) < 2e-12
    # The literal pulse word also binds all scalar-phase addresses, not just
    # a success matrix modulo a global phase or an ideal gate target.
    assert source.literal_word(actual) == field.source_word(geometry, float(duration), steps)
    assert actual.resources['native_pulses'] == field.resources(geometry, float(duration), steps)['native_pulses']
    exact = expm(-1j*float(duration)*field.FullFock(8).hamiltonian(geometry).toarray())
    bound = source.split_bound(geometry.sides, duration, steps, refined=True)
    assert np.linalg.norm(actual.success-exact, 2) <= float(bound)


def test_complete_failure_archive_matches_literal_nonadjacent_remaps(monkeypatch):
    angle = 1e-8
    full = np.eye(6, dtype=complex)
    full[np.ix_((0, 2), (0, 2))] = [[math.cos(angle), 1j*math.sin(angle)],
                                     [1j*math.sin(angle), math.cos(angle)]]
    monkeypatch.setattr(compiler, 'native_pulse', lambda kind, theta: full.copy())
    word = [compiler.pulse('X', 0, [0]), compiler.pulse('X', 0, [1])]
    macro = source.compile_macro(word, 2, np.eye(4), source.ENCODER_GATE)
    # Interleaved and reversed local placement makes untouched-output ordering
    # observable. Include genuinely complex entanglement, not only basis data.
    initial = np.arange(24).reshape(8, 3)+1j*np.arange(24, 48).reshape(8, 3)
    initial = initial/np.linalg.norm(initial)
    schedule = (source.Placement(macro, (2, 0)), source.Placement(macro, (1, 2)))
    archive = source.execute_schedule(schedule, 3, initial)
    expected, failures = execute_with_exports(source.literal_word(archive), 3, initial)
    assert np.linalg.norm(archive.success-expected) < 2e-14
    assert archive.failure_probability > 1e-17
    assert archive.failure_probability == pytest.approx(failure_probability(failures), rel=2e-14)
    complete_reference = archive.success.conj().T@archive.success
    for macro_index in range(2):
        for event in range(2):
            actual = archive.failure(macro_index, event)
            for a, b in zip(actual, failures[2*macro_index+event], strict=True):
                assert np.linalg.norm(a-b) < 1e-22
                complete_reference += a.conj().T@a
    assert np.linalg.norm(complete_reference-initial.conj().T@initial) < 2e-14
    # Later source mutations and caller-buffer changes cannot alter the owned
    # archive's quantum outputs.
    baseline = archive.failure(1, 1)
    initial[:] = 0; full[:] = 0
    monkeypatch.setattr(compiler, 'native_pulse', lambda *_: (_ for _ in ()).throw(RuntimeError('not owned')))
    for a, b in zip(archive.failure(1, 1), baseline, strict=True):
        assert np.array_equal(a, b)


@pytest.mark.parametrize('basis', ('occupation', 'candidate'))
def test_full_read_law_keeps_all_outputs_and_reference(basis):
    indices = np.arange(512).reshape(256, 2)
    state = np.sin(indices+1)+1j*np.cos(3*indices+2)
    state /= np.linalg.norm(state)
    decoder, pointer, outputs = source.read(state, 2, basis=basis)
    ideal = state if decoder is None else np.kron(source._encoder_data(2)[0].conj().T,
                                                source._encoder_data(2)[0].conj().T)@state
    assert len(outputs) == 256
    for record, output in enumerate(outputs):
        expected = np.zeros_like(state); expected[record] = ideal[record]
        assert np.linalg.norm(output-expected) < 1e-10
    assert abs(sum(np.linalg.norm(a)**2 for a in outputs)+pointer.failure_probability
               +(0 if decoder is None else decoder.failure_probability)-1) < 1e-10
    # Two refinements agree only on block number bins, not on each raw flag.
    original_bins = [sum(np.linalg.norm(state[j])**2 for j in range(256)
                    if (j & 15).bit_count() == a and (j >> 4).bit_count() == b)
                    for a in range(5) for b in range(5)]
    decoded_bins = [sum(np.linalg.norm(outputs[j])**2 for j in range(256)
                   if (j & 15).bit_count() == a and (j >> 4).bit_count() == b)
                   for a in range(5) for b in range(5)]
    assert np.linalg.norm(np.array(original_bins)-decoded_bins) < 1e-10


def test_canonical_long_history_positive_transfer_and_complete_counts():
    result = source.run(basis='candidate')
    ideal = expm(-4j*chain.setup().K)@chain.bare_columns('hole')[:, 0]
    exact_probability = np.linalg.norm(chain.bare_columns('hole')[:, 2:].conj().T@ideal)**2
    assert exact_probability == pytest.approx(.005851526872527, abs=2e-14)
    assert result['candidate_location'][1] == pytest.approx(exact_probability, abs=1e-7)
    assert result['resources']['native_pulses'] == 307625
    assert result['resources']['proper_code_service_events'] == 615250
    assert result['ideal_split_operator_bound'] == F(841, 1572864)
    assert abs(sum(result['probabilities'])+result['failure_probability']-1) < 1e-7
    assert result['candidate_location'][2] > .03
    # Removing the middle bond is an independent no-transfer control.
    decoupled = expm(-4j*chain.setup(link_strength=F(0)).K)@chain.bare_columns('hole')[:, 0]
    assert np.linalg.norm(chain.bare_columns('hole')[:, 2:].conj().T@decoupled) < 2e-13


@pytest.mark.parametrize('spin', ([0, 0], [2, 0], [1e308, 1e308], [[1, 0], [0, 1]]))
def test_probability_experiment_rejects_nonunit_supplied_input(spin):
    with pytest.raises(ValueError, match='normalized supplied spin/reference'):
        source.run(spin=spin, duration=F(0), steps=1)


def test_probability_experiment_accepts_complex_entangled_reference_without_normalizing():
    spin = np.array([[1+2j, -3j, .5], [2j, 1-1j, -2]], complex)/math.sqrt(24.25)
    result = source.run(spin=spin, duration=F(0), steps=1, basis='candidate')
    target = chain.bare_columns('hole')[:, :2]@spin
    assert np.linalg.norm(result['preparation'].success-target) < 1e-10
    assert result['outcomes'][49].shape == (256, 3)
    assert abs(math.fsum(result['probabilities'])+result['failure_probability']-1) < 1e-10
    assert result['candidate_location'][0] == pytest.approx(1, abs=1e-10)


def test_composition_uses_norm_growth_not_linear_residual_sum():
    archive = source.execute_schedule((source.Placement(source._local('x'), (0,)),)*10, 1, [1, 0])
    s = 10*source.SMALL_GATE
    assert source.composition_bound([archive]) == s/(1-s)
    assert source.composition_bound([archive]) >= (1+source.SMALL_GATE)**10-1 > s


@pytest.mark.parametrize('damage', ('schedule', 'initial', 'success', 'failure', 'resources', 'bits', 'targets'))
@pytest.mark.parametrize('in_place', (False, True))
def test_owned_archive_rejects_forged_or_warm_mutated_evidence(damage, in_place, monkeypatch):
    # A valid nonidentity instrument with a genuinely nonzero abort output:
    # the independent complement rotation makes omission measurable.
    angle = 1e-8
    full = np.eye(6, dtype=complex)
    full[:2, :2] = [[0, 1], [1, 0]]
    turn = np.eye(6, dtype=complex)
    turn[np.ix_((0, 2), (0, 2))] = [[math.cos(angle), 1j*math.sin(angle)],
                                    [1j*math.sin(angle), math.cos(angle)]]
    monkeypatch.setattr(compiler, 'native_pulse', lambda *_: turn@full)
    macro = source.compile_macro([compiler.pulse('X', 0, [0])], 1,
                                 np.array([[0, 1], [1, 0]]), source.ENCODER_GATE)
    initial = np.array([[1, 2j], [-3j, 4]], complex)/math.sqrt(30)
    archive = source.execute_schedule((source.Placement(macro, (0,)),), 1, initial)
    baseline = archive.failure(0, 0)
    assert archive.failure_probability > 0
    assert sum(np.linalg.norm(row)**2 for row in baseline) == pytest.approx(
        archive.failure_probability, rel=2e-14)
    assert source.composition_bound([archive]) > 0
    assert len(source.literal_word(archive)) == 1
    complete = archive.success.conj().T@archive.success
    for row in baseline: complete += row.conj().T@row
    assert np.linalg.norm(complete-initial.conj().T@initial) < 2e-14
    assert not archive.initial.flags.writeable and not archive.success.flags.writeable

    field, changed = {
        'schedule': ('schedule', ()),
        'initial': ('initial', np.zeros_like(archive.initial)),
        'success': ('success', np.zeros_like(archive.success)),
        'failure': ('failure_probability', -1.),
        'resources': ('resources', {**archive.resources, 'native_pulses': 0}),
        'bits': ('bits', True),
        'targets': ('schedule', (source.Placement(macro, (False,)),)),
    }[damage]
    if in_place:
        # Frozen dataclasses stop ordinary assignment, but verification also
        # rejects deserialization/object-level tampering after warm use.
        object.__setattr__(archive, field, changed)
        forged = archive
    else:
        forged = replace(archive, **{field: changed})
    for consumer in (lambda: source.composition_bound([forged]),
                     lambda: source.literal_word(forged),
                     lambda: forged.failure(0, 0)):
        with pytest.raises(ValueError, match='unchanged, owned complete archive'):
            consumer()


@pytest.mark.parametrize('change', ({'failure_gram': -np.eye(2)}, {'success': np.zeros((2, 2))},
                                  {'gate': F(0)}, {'resources': {}}, {'bits': True}))
def test_forged_public_macro_is_rejected(change):
    genuine = source._local('x')
    forged = replace(genuine, **change)
    with pytest.raises(ValueError, match='owned|distinct'):
        source.execute_schedule((source.Placement(forged, (0,)),), 1, [1, 0])


@pytest.mark.parametrize('axis', (True, 0., np.int64(0), -1, 3))
def test_warm_axis_guards(axis):
    source.encoder(0)
    with pytest.raises(ValueError): source.encoder(axis)
    with pytest.raises(ValueError): source.prepare(axis=axis)
    with pytest.raises(ValueError): source.evolve(np.eye(256)[:, 0], axis)


@pytest.mark.parametrize('value', (True, 1., float('nan'), '1', F(-1), F(5)))
def test_exact_duration_guards(value):
    with pytest.raises(ValueError): source.evolution_schedule(duration=value)


def test_malformed_fraction_and_nonfinite_inputs_fail_closed():
    bad = F(1, 2); bad._denominator = -2
    with pytest.raises(ValueError): source.evolution_schedule(duration=bad)
    with pytest.raises(ValueError): source.compile_macro([], 1, np.eye(2), gate=bad)
    with pytest.raises(ValueError): source.compile_macro([], 1, np.eye(2), gate=True)
    for value in (np.array([1, np.nan]), np.ma.array([1, 0]), np.empty((2, 0))):
        with pytest.raises(ValueError): source.prepare(spin=value)
    with pytest.raises(ValueError): source.read(np.zeros(256), basis='postselect')
    with pytest.raises(ValueError): source.read(np.zeros((256, 17)))


def test_macro_gate_is_owned_and_bound_as_original_integer_components():
    supplied = F(1, 10**10)
    macro = source.compile_macro([], 1, np.eye(2), supplied)
    archive = source.execute_schedule((source.Placement(macro, (0,)),), 1, [1, 0])
    supplied._numerator = -1
    assert macro.gate == F(1, 10**10)
    macro.gate._numerator = -1
    with pytest.raises(ValueError, match='unchanged'):
        source.execute_schedule((source.Placement(macro, (0,)),), 1, [1, 0])
    with pytest.raises(ValueError, match='unchanged'):
        source.composition_bound([archive])


def test_wrong_phase_or_leaking_source_primitive_fails(monkeypatch):
    for full in (np.exp(.001j)*np.eye(6), np.full((6, 6), np.nan), np.zeros((6, 6))):
        monkeypatch.setattr(compiler, 'native_pulse', lambda *_, full=full: full)
        with pytest.raises((ValueError, ArithmeticError)):
            source.compile_macro([compiler.pulse('X', 0, [0])], 1, np.eye(2))


def test_certificate_bound_is_positive_and_scopes_are_explicit():
    value = source.certificate()
    assert len(value['rows']) == 4
    for row in value['rows']:
        assert row['raw_occupation_outcomes'] == 256
        assert row['all_number_sectors_and_export_failures_retained']
        assert row['initialization_resources']['supplied_field_blank_bits'] == 7
        assert row['initialization_resources']['fresh_pointer_blank_bits'] == 8
        assert row['initialization_resources']['maximum_live_work_bits'] == 16
        if row['read_basis'] == 'candidate':
            lower = F(row['conditional_positive_candidate_transfer_probability_lower'])
            assert .003 < float(lower) < row['candidate_location_probabilities'][1]
    assert not value['initial_reference_is_coupled_ground']
    assert not value['native_passive_protection_or_topological_defect_claimed']


def test_fresh_certificate_rejects_changed_native_pulse_after_warm_execution(monkeypatch):
    source.run(duration=F(1, 8), steps=1, basis='candidate')
    native = compiler.native_pulse
    monkeypatch.setattr(compiler, 'native_pulse', lambda kind, angle: np.exp(.001j)*native(kind, angle))
    with pytest.raises(ArithmeticError, match='complete local source instrument'):
        source.certificate()


def test_fresh_certificate_reloads_actual_control_files_after_warm_execution(tmp_path, monkeypatch):
    source.run(duration=F(1, 8), steps=1, basis='candidate')
    for relative in (compiler.RESPONSE, compiler.CONTROL_RECEIPT):
        destination = tmp_path/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((compiler.ROOT/relative).read_bytes())
    target = tmp_path/compiler.CONTROL_RECEIPT
    mutated = json.loads(target.read_text(encoding='utf-8'))
    mutated['evidence']['controls']['rows'][0]['coefficients'][0] = ['123', '0']
    target.write_text(json.dumps(mutated), encoding='utf-8')
    monkeypatch.setattr(compiler, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match='exact source control mismatch'):
        source.certificate()
