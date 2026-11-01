"""Complete source histories and adversarial controls outside the cat code."""
from dataclasses import replace
from fractions import Fraction
import math
from types import MappingProxyType
import warnings

import numpy as np
import pytest

from . import occupation_service as s
from ..exchange_sector import many_body
from ..observer.source_scalar import compiler
from ..reference_observation.observation import frozen


@pytest.fixture(scope='module')
def histories():
    return {(law, sign, family, control): s.experiment(law, sign, family=family,
             control=control, return_run=True)
            for law in ('C', 'E') for sign in (1, -1)
            for family, control in (('extreme', 'none'), ('fixed_number', 'none'))}


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('sign', [1, -1])
def test_complete_word_on_every_native_input_including_signing_spectators(law, sign):
    result = s.certify_repair(law, sign)
    assert result['all_native_input_columns'] == 4096
    assert result['compact_active_input_columns'] == 64
    assert result['complete_raw_slots_per_coin'] == (4 if law == 'C' else 8)
    assert Fraction(result['full_input_complete_amplitude_bound']) > Fraction(result['ideal_word_amplitude_ceiling'])
    # A moving one-particle state and an intervening spectator produce the
    # independent CAR contrast, even though both extreme code columns hold.
    target = s.target_columns(law, sign, (1, 3))
    for column, mask in enumerate((1, 3)):
        for coin in (0, 1):
            raw, moved, phase, _ = many_body.column(mask, law, s.SEAM, coin, sign)
            assert target[moved+(coin << 12)+(raw << 13), column] == phase/math.sqrt(2)


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('sign', [1, -1])
def test_extreme_code_unknown_complex_reference_survives_all_complete_outputs(histories, law, sign):
    report, run, target = histories[law, sign, 'extreme', 'none']
    assert report['encoded_occupation_numbers'] == [0, 12]
    assert report['cat_cross_number_phase_apparatus_supplied']
    assert report['program_class_v2_execution_claimed'] is False
    assert report['encoded_preparation_preserves_declared_total'] is False
    # Arbitrary reference columns, with complex and unequal Schmidt weights.
    reference = np.array([[1, 2j], [3+1j, -2]], complex)/math.sqrt(19)
    assert np.linalg.norm((run.success-target)@reference) < float(Fraction(report['complete_output_amplitude_bound']))
    assert np.linalg.norm(run.success.conj().T@run.success+run.failure_gram-np.eye(2)) < 2e-11
    for column in range(2):
        support = np.flatnonzero(target[:, column])
        assert set(support) == {column, column+(1 << 12)}
    # No input-dependent tape distribution, and no hidden normalization.
    assert np.allclose(np.sum(abs(run.success)**2, axis=0), 1, rtol=0, atol=2e-11)


def _source_density(target, column):
    keep = np.flatnonzero(np.any(abs(target.reshape(-1, 4096, 2)) > 0, axis=(0, 2)))
    pieces = target.reshape(-1, 4096, 2)[:, keep, column]
    return pieces.T@pieces.conj()


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('sign', [1, -1])
def test_fixed_six_particle_control_loses_phase_into_complete_native_records(histories, law, sign):
    report, run, target = histories[law, sign, 'fixed_number', 'none']
    assert report['supplied_initial_field_occupations'] == [6, 6]
    assert report['encoded_occupation_numbers'] == [6, 6]
    assert report['encoded_preparation_preserves_declared_total']
    assert report['primitive_encoded_number_conservation_claimed'] is False
    assert report['program_class_v2_execution_claimed'] is False
    assert report['cat_cross_number_phase_apparatus_supplied'] is False
    assert report['postselection'] is False
    assert all((i & s.FULL).bit_count() == 6 for i in np.flatnonzero(np.any(target != 0, axis=1)))
    assert np.array_equal(abs(target[:, 0])**2, abs(target[:, 1])**2)
    # Complete output is an isometry; committing/ignoring its native record
    # removes the distinguishing cat phase from the source output.
    assert np.linalg.norm(target.conj().T@target-np.eye(2)) < 1e-14
    assert np.linalg.norm(_source_density(target, 0)-_source_density(target, 1)) < 1e-14
    addresses = np.flatnonzero(np.any(target != 0, axis=1))
    assert {int(i >> 13) for i in addresses} == ({0, 1, 2} if law == 'C' else {4, 5, 6})
    assert run.failure_weights.shape == (report['resources']['native_pulses'], 4)


@pytest.mark.parametrize('control', ['local_z', 'which_path'])
def test_one_port_phase_and_retained_pointer_controls(control):
    report, run, target = s.experiment('E', -1, control=control, return_run=True)
    assert report['work_bits'] == (17 if control == 'which_path' else 16)
    if control == 'local_z':
        assert {int(i & s.FULL) for i in np.flatnonzero(target[:, 0])} == {1}
        assert {int(i & s.FULL) for i in np.flatnonzero(target[:, 1])} == {0}
    else:
        assert np.linalg.norm(_source_density(target, 0)-np.eye(2)/2) < 1e-14
        assert np.linalg.norm(_source_density(target, 1)-np.eye(2)/2) < 1e-14
        assert {int(i >> 16) for i in np.flatnonzero(np.any(target != 0, axis=1))} == {0, 1}
    assert run.failure_probability >= 0


def test_preparation_macros_preserve_number_on_all_their_input_columns():
    report = s.preparation_certificate()
    assert report['fixed_number_preparation_macro_columns'] == 12
    assert report['ideal_macros_preserve_total_number_on_all_inputs']


def test_identity_repair_word_cannot_pass_via_extreme_code_invariance(monkeypatch):
    monkeypatch.setattr(s, 'repair_word', lambda law, sign: [compiler.pulse('Y', -math.pi/4, [12])])
    with pytest.raises(ArithmeticError, match='general-input'):
        s.certify_repair('C', -1)


def test_wrong_spectator_phase_cannot_pass_general_input_guard(monkeypatch):
    original = s.repair_word
    monkeypatch.setattr(s, 'repair_word', lambda law, sign: original(law, 1))
    with pytest.raises(ArithmeticError, match='general-input'):
        s.certify_repair('E', -1)


def test_bad_but_unitary_primitive_is_rejected(monkeypatch):
    original = compiler.native_pulse
    def wrong(kind, angle):
        if kind == 'H':
            return np.eye(6, dtype=complex)
        return original(kind, angle)
    monkeypatch.setattr(compiler, 'native_pulse', wrong)
    with pytest.raises(ArithmeticError, match='primitive complete-output'):
        s.pulse_budget(compiler.cnot(0, 1))


def test_small_nonzero_complement_is_retained_and_corrupt_archive_rejected(monkeypatch):
    original = compiler.native_pulse
    def leaky(kind, angle):
        full = original(kind, angle)
        turn = np.eye(6, dtype=complex)
        theta = 1e-12
        turn[np.ix_([0, 2], [0, 2])] = [[math.cos(theta), -math.sin(theta)],
                                       [math.sin(theta), math.cos(theta)]]
        return turn@full
    monkeypatch.setattr(compiler, 'native_pulse', leaky)
    word = [compiler.pulse('Y', .2, [0])]
    initial = np.eye(2, dtype=complex)
    target = s._ideal_pulse(word[0])[:2]
    run = s.execute_complete(word, 1, initial)
    s.admit(run, target, word=word, bits=1, initial=initial)
    assert run.failure_probability > 1e-26
    assert np.any(run.failure_weights > 0)
    failed = run.failure(0)
    assert len(failed) == 4
    assert sum(np.vdot(x, x).real for x in failed) == pytest.approx(run.failure_probability, rel=1e-14, abs=0)
    corrupted = run.failure_weights.copy()
    corrupted[0, 0] = 0
    with pytest.raises(ArithmeticError, match='failure-prefix'):
        s.admit(replace(run, failure_weights=frozen(corrupted)), target, word=word, bits=1, initial=initial)
    altered = run.success.copy(); altered[0, 0] += 1e-12
    with pytest.raises(ArithmeticError, match='failure-prefix'):
        s.admit(replace(run, success=frozen(altered)), target, word=word, bits=1, initial=initial)
    with pytest.raises(ArithmeticError, match='initial'):
        s.admit(run, target, word=word, bits=1, initial=initial*1j)


def _mutable_snapshot(value, kind):
    if kind == 'buffer':
        result = np.frombuffer(bytearray(value.tobytes()), dtype=value.dtype).reshape(value.shape)
    else:
        result = value.copy()
    if kind != 'writable':
        result.flags.writeable = False
    return result


@pytest.mark.parametrize('field', ['initial', 'success', 'failure_weights', 'failure_gram', 'primitive'])
@pytest.mark.parametrize('kind', ['writable', 'readonly_owner', 'buffer'])
def test_archive_rejects_mutable_backing_even_when_values_and_readonly_flags_match(field, kind):
    word = [compiler.pulse('Y', .2, [0])]
    initial = np.eye(2, dtype=complex)
    run = s.execute_complete(word, 1, initial)
    target = s._ideal_pulse(word[0])[:2]
    with pytest.raises(ValueError, match='immutable bytes-backed'):
        if field == 'primitive':
            entries = tuple((key, _mutable_snapshot(value, kind)) for key, value in run.primitives.items())
            forged = replace(run, primitives=s._FrozenPrimitives(entries))
        else:
            forged = replace(run, **{field: _mutable_snapshot(getattr(run, field), kind)})
        s.admit(forged, target, word=word, bits=1, initial=initial)


def test_archive_rejects_mapping_proxy_with_caller_retained_mutable_dictionary():
    word = [compiler.pulse('Y', .2, [0])]
    initial = np.eye(2, dtype=complex)
    run = s.execute_complete(word, 1, initial)
    target = s._ideal_pulse(word[0])[:2]
    catalogue = dict(run.primitives)
    forged = replace(run, primitives=MappingProxyType(catalogue))
    with pytest.raises(ValueError, match='literal complete source archive'):
        s.admit(forged, target, word=word, bits=1, initial=initial)
    # Accepted archives have neither mutable array owners nor dictionary aliases.
    s.admit(run, target, word=word, bits=1, initial=initial)
    for array in (run.initial, run.success, run.failure_weights, run.failure_gram, *run.primitives.values()):
        with pytest.raises(ValueError):
            array.flags.writeable = True
    with pytest.raises((AttributeError, TypeError)):
        run.primitives.entries = ()


def test_archive_envelope_has_no_mutable_instance_dictionary():
    word = [compiler.pulse('Y', .2, [0])]
    initial = np.eye(2, dtype=complex)
    run = s.execute_complete(word, 1, initial)
    target = s._ideal_pulse(word[0])[:2]
    s.admit(run, target, word=word, bits=1, initial=initial)
    assert not hasattr(run, '__dict__')
    with pytest.raises((AttributeError, TypeError)):
        run.initial = frozen(2*initial)
    # The earlier shared frozen dataclass still has a mutable envelope and
    # cannot be passed off as this independently bounded archive.
    legacy = s.recovery.CompleteRun(run.word, run.bits, run.initial, run.primitives,
                                   run.success, run.failure_weights, run.failure_gram)
    assert hasattr(legacy, '__dict__')
    with pytest.raises(ValueError, match='literal complete source archive'):
        s.admit(legacy, target, word=word, bits=1, initial=initial)


@pytest.mark.parametrize('field', ['initial', 'success', 'failure_weights', 'failure_gram',
                                  'primitive', 'primitive_entry'])
def test_public_array_metadata_and_base_mutations_do_not_change_archive(field):
    word = [compiler.pulse('Y', .2, [0])]
    initial = np.eye(2, dtype=complex)
    run = s.execute_complete(word, 1, initial)
    target = s._ideal_pulse(word[0])[:2]
    def exposed():
        if field == 'primitive':
            return run.primitives['Y', .2]
        if field == 'primitive_entry':
            return run.primitives.entries[0][1]
        return getattr(run, field)
    original = exposed().copy()
    output = exposed()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        output.shape = (output.size,)
        output.dtype = np.uint8
        if type(output.base) is np.ndarray:
            output.base.shape = (1, output.base.size)
            output.base.dtype = np.uint8
    assert exposed().shape == original.shape
    assert exposed().dtype == original.dtype
    assert np.array_equal(exposed(), original)
    s.admit(run, target, word=word, bits=1, initial=initial)


def test_unrepresentable_input_gram_is_refused_before_execution():
    with pytest.raises(ValueError, match='reference Gram'):
        s.execute_complete([], 1, 1e200*np.eye(2, dtype=complex))
    # Finite, unnormalized linear input is not rejected merely for its scale.
    initial = 1e100*np.eye(2, dtype=complex)
    run = s.execute_complete([], 1, initial)
    error, bound = s.admit(run, initial, word=[], bits=1, initial=initial)
    assert error == 0 and bound == 0
    assert np.array_equal(run.success, initial)


@pytest.mark.parametrize('dtype', [np.float32, np.float64, np.longdouble, np.clongdouble])
def test_exactly_representable_input_and_target_values_are_preserved(dtype):
    initial = np.array([[1., .5], [0., -.25]], dtype=dtype)
    original = initial.copy()
    run = s.execute_complete([], 1, initial)
    assert s.admit(run, initial, word=[], bits=1, initial=initial) == (0., Fraction(0))
    assert run.initial.dtype == np.complex128
    assert np.array_equal(run.initial.real, original)
    assert np.array_equal(initial, original)


@pytest.mark.skipif(np.finfo(np.longdouble).nmant <= np.finfo(np.float64).nmant,
                    reason='platform longdouble has no extra mantissa precision')
@pytest.mark.parametrize('dtype', [np.longdouble, np.clongdouble])
@pytest.mark.parametrize('field', ['initial', 'target'])
def test_wider_precision_values_cannot_be_silently_rounded(dtype, field):
    initial = np.eye(2, dtype=complex)
    wider = np.eye(2, dtype=dtype)
    wider[0, 0] = np.nextafter(np.longdouble(1), np.longdouble(2))
    assert wider[0, 0] != 1
    with pytest.raises(ValueError, match=f'{field} requires a value-preserving complex128'):
        if field == 'initial':
            s.execute_complete([], 1, wider)
        else:
            run = s.execute_complete([], 1, initial)
            s.admit(run, wider, word=[], bits=1, initial=initial)


@pytest.mark.skipif(np.finfo(np.longdouble).minexp >= -1074,
                    reason='platform longdouble has no exponent range below binary64 subnormals')
@pytest.mark.parametrize('field', ['initial', 'target'])
def test_wider_nonzero_values_cannot_be_silently_erased(field):
    initial = np.eye(2, dtype=complex)
    wider = np.eye(2, dtype=np.longdouble)
    wider[0, 1] = np.longdouble('1e-400')
    assert wider[0, 1] > 0
    with pytest.raises(ValueError, match=f'{field} requires a value-preserving complex128'):
        if field == 'initial':
            s.execute_complete([], 1, wider)
        else:
            run = s.execute_complete([], 1, initial)
            s.admit(run, wider, word=[], bits=1, initial=initial)


@pytest.mark.parametrize('bits,columns', [(0, 1), (18, 1), (True, 1), (1, 0), (1, 3)])
def test_executor_strict_bounded_domain(bits, columns):
    initial = np.zeros((2, columns), complex)
    with pytest.raises(ValueError):
        s.execute_complete([], bits, initial)


def test_nonfinite_input_and_oversized_word_rejected():
    initial = np.eye(2, dtype=complex); initial[0, 0] = np.nan
    with pytest.raises(ValueError):
        s.execute_complete([], 1, initial)
    with pytest.raises(ValueError, match='pulse domain'):
        s.execute_complete([compiler.pulse('I', 0., [0])]*(s.MAX_PULSES+1), 1, np.eye(2, dtype=complex))


@pytest.mark.parametrize('kwargs', [dict(law='bad', character=1), dict(law='C', character=True),
                                  dict(law='C', character=-1, control='bad'),
                                  dict(law='C', character=-1, family='fixed_number', control='which_path')])
def test_experiment_rejects_ambiguous_or_unsupported_requests(kwargs):
    with pytest.raises(ValueError):
        s.experiment(**kwargs)
