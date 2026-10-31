"""Independent full-carrier and hostile controls for supplied candidate transport."""
from copy import deepcopy
import math

import numpy as np
import pytest

from . import algebra, source
from ..particle_response import dimer
from ..observer.source_scalar import compiler
from ..full_fock_field.field import FieldGeometry


def literal_permutation(modes=8, fermionic=True):
    result = np.zeros((1 << modes, 1 << modes), complex)
    for occupation in range(1 << modes):
        moved = [(mode+4) % 8 for mode in range(modes) if occupation >> mode & 1]
        inversions = sum(a > b for i, a in enumerate(moved) for b in moved[i+1:])
        result[sum(1 << p for p in moved), occupation] = (-1)**inversions if fermionic else 1
    return result


@pytest.mark.parametrize("first,second", [(1, 2), (3, 5), (3, 12), (7, 14)])
def test_equal_number_selective_rotations_on_every_Fock_column(first, second):
    angle = .317
    actual = compiler.execute(source.number_rotation_word(first, second, angle), 4)
    expected = np.eye(16, dtype=complex)
    expected[np.ix_([first, second], [first, second])] = [
        [math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]]
    np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=0)
    number = np.diag([i.bit_count() for i in range(16)])
    np.testing.assert_allclose(actual@number, number@actual, atol=2e-11, rtol=0)


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("kind", ["ground", "hole", "particle"])
def test_preparation_has_actual_fixed_charge_and_predecessor_spectral_support(axis, kind):
    prepared = source.prepare(kind, axis)
    state = np.asarray(prepared.state).reshape(16, -1)
    initial = np.asarray(prepared.initial).reshape(16, -1)
    number = source.SEEDS[kind][0].bit_count()
    values = np.array([i.bit_count() for i in range(16)])
    assert np.count_nonzero(initial[values != number]) == 0
    np.testing.assert_allclose(state[values != number], 0, atol=2e-11, rtol=0)
    sides = tuple(2 if j == axis else 1 for j in range(3))
    spectrum = dimer.setup(FieldGeometry(sides))
    projector = spectrum.reference if kind == "ground" else dimer.candidate_projectors(spectrum)[
        "charge_minus_one_hole" if kind == "hole" else "charge_plus_one_particle"]
    # The predecessor's independently assembled spectral projector checks the
    # phase-fixed analytic columns used to request the encoder.
    np.testing.assert_allclose(projector@state, state, atol=2e-11, rtol=0)
    source._gram(prepared.state, prepared.failures)
    assert len(prepared.failures) == prepared.resources["native_pulses"]
    assert source.failure_probability(prepared.failures) < 1e-20
    assert prepared.error < 2e-11


def test_full_encoder_other_sectors_are_identity_and_no_caller_state_substitution():
    initial, target, unitary, word = source._encoder("ground", 0)
    actual = compiler.execute(word, 4)
    for number in (0, 1, 3, 4):
        addresses = [i for i in range(16) if i.bit_count() == number]
        np.testing.assert_allclose(actual[:, addresses], np.eye(16)[:, addresses], atol=2e-11, rtol=0)
    np.testing.assert_allclose(actual@initial, target, atol=2e-11, rtol=0)
    with pytest.raises(TypeError):
        source.prepare("ground", state=target)
    prepared = source.prepare("ground")
    with pytest.warns(DeprecationWarning, match="Setting the shape"):
        prepared.state.shape = (4, 4)
    prepared.word[0]["angle"] = 1.234
    later = source.prepare("ground")
    assert later.state.shape == (16,) and later.error < 2e-11


def test_actual_controlled_route_on_all512_columns_and_allnumbersectors():
    expected = np.eye(512, dtype=complex)
    expected[256:, 256:] = literal_permutation()
    actual = compiler.execute(source.block_swap_word(controlled=True), 9)
    np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=0)
    wrong = compiler.execute(source.block_swap_word(controlled=True, ordinary=True), 9)
    expected_wrong = np.eye(512, dtype=complex)
    expected_wrong[256:, 256:] = literal_permutation(fermionic=False)
    np.testing.assert_allclose(wrong, expected_wrong, atol=2e-11, rtol=0)
    assert np.linalg.norm(actual-wrong, 2) > 1.99


@pytest.mark.parametrize("kind", ["hole", "particle"])
def test_transport_unknown_complex_spin_entangled_reference_and_exterior(kind):
    result = source.transport(kind)
    coefficients = np.array([[1, 2j, -1j], [2, -1, 3j]], complex)
    coefficients /= np.linalg.norm(coefficients)
    actual, target = result["state"]@coefficients, result["target"]@coefficients
    np.testing.assert_allclose(actual, target, atol=2e-11, rtol=0)
    np.testing.assert_allclose(result["state"].conj().T@result["state"], np.eye(2), atol=2e-11, rtol=0)
    # After transport the low four modes are precisely the same ground exterior,
    # even when the charged internal state and an arbitrary reference are mixed.
    tensor = actual.reshape(16, 16, 3)
    exterior = np.einsum("bar,bcr->ac", tensor, tensor.conj())
    ground = algebra.ground_column()
    np.testing.assert_allclose(exterior, np.outer(ground, ground.conj()), atol=2e-11, rtol=0)
    assert result["resources"]["native_pulses"] == 6226
    assert result["initialized_field_number"] == (3 if kind == "hole" else 5)
    assert result["first_failure_kraus_slots"] > 0


@pytest.mark.parametrize("kind,ordinary,dephase,expected", [
    ("hole", False, False, [0, 1]), ("particle", False, False, [0, 1]),
    ("ground", False, False, [1, 0]), ("hole", True, False, [1, 0]),
    ("hole", False, True, [.5, .5]),
])
def test_exchange_and_complete_countercontrols(kind, ordinary, dephase, expected):
    result = source.exchange(kind, ordinary=ordinary, dephase=dephase)
    values = [sum(r["probability"] for r in result["outcomes"] if r["pointer"] == bit)
              for bit in (0, 1)]
    np.testing.assert_allclose(values, expected, atol=2e-11, rtol=0)
    assert abs(result["complete_probability"]-1) < 2e-11
    assert len(result["outcomes"]) == (4 if dephase else 2)
    if dephase:
        np.testing.assert_allclose([r["probability"] for r in result["outcomes"]], [.25]*4,
                                   atol=2e-11, rtol=0)
    assert result["first_failure_kraus_slots"] > 0
    assert result["neutral_pointer_reads"] == (2 if dephase else 1)


def test_injected_valid_unitary_export_leakage_is_retained_and_rejected(monkeypatch):
    original = compiler.native_pulse
    def leaked(kind, angle):
        full = np.array(original(kind, angle), copy=True)
        turn = np.eye(6, dtype=complex)
        turn[np.ix_([0, 5], [0, 5])] = [[math.cos(.001), -math.sin(.001)],
                                       [math.sin(.001), math.cos(.001)]]
        return turn@full
    monkeypatch.setattr(compiler, "native_pulse", leaked)
    with pytest.raises(ArithmeticError, match="Stinespring"):
        source.prepare("ground")


def test_global_source_phase_is_part_of_the_preparation_admission(monkeypatch):
    original = compiler.native_pulse
    def shifted(kind, angle):
        value = original(kind, angle)
        return np.exp(.0001j)*value if kind == "I" else value
    monkeypatch.setattr(compiler, "native_pulse", shifted)
    with pytest.raises(ArithmeticError, match="Stinespring"):
        source.prepare("hole")


@pytest.mark.parametrize("kwargs", [dict(kind=True), dict(kind="unknown"), dict(axis=True),
                                   dict(axis=-1), dict(axis=3), dict(axis=1.0)])
def test_preparation_refuses_coercions_and_unknown_domains(kwargs):
    with pytest.raises(ValueError):
        source.prepare(**kwargs)


@pytest.mark.parametrize("args", [(1, 3, .2), (True, 2, .2), (1, 2, float("nan")),
                                 (1, 2, True), (2, 1, .2), (1, 32, .2)])
def test_number_rotations_refuse_wrong_support_and_bad_angles(args):
    with pytest.raises(ValueError):
        source.number_rotation_word(*args)


def test_certificate_binds_emitted_route_instead_of_only_known_witness(monkeypatch):
    original = source.block_swap_word
    def corrupted(*, controlled=False, ordinary=False):
        word = deepcopy(original(controlled=controlled, ordinary=ordinary))
        # Extra phase on an exterior full-carrier direction must fail the
        # global route check even if a selected prepared input missed it.
        word += compiler.pauli_rotation("Z"+"I"*(8 if controlled else 7), .01)
        return word
    monkeypatch.setattr(source, "block_swap_word", corrupted)
    with pytest.raises(ArithmeticError, match="full-Fock route"):
        source._swap_certificate()


@pytest.mark.parametrize("kind", ["hole", "particle"])
def test_complete_presence_read_preserves_spin_and_rejects_every_other_sector(kind):
    # Independently assembled predecessor projector, including the upper doublet
    # and every other occupation sector in its orthogonal complement.
    spectrum = dimer.setup(FieldGeometry((2, 1, 1)))
    projector = dimer.candidate_projectors(spectrum)[
        "charge_minus_one_hole" if kind == "hole" else "charge_plus_one_particle"]
    complement = np.eye(16)-projector
    actual = compiler.execute(source.presence_readout_word(kind), 5)
    expected = np.block([[complement, projector], [projector, complement]])
    np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=0)
    ground = algebra.ground_column()
    assert np.linalg.norm(actual[16:, :16]@ground) < 2e-11
    upper = algebra.candidate_columns(kind, band="upper")
    assert np.linalg.norm(actual[16:, :16]@upper) < 2e-11
    number = source.SEEDS[kind][0].bit_count()
    other = [j for j in range(16) if j.bit_count() != number]
    assert np.linalg.norm(actual[16:, other]) < 2e-11
    assert source._presence_certificate(kind)["full_Fock_and_pointer_columns"] == 32


@pytest.mark.parametrize("kind", ["hole", "particle"])
def test_joined_source_transport_and_presence_read_keeps_unknown_reference(kind):
    result = source.transport_readout(kind)
    coefficients = np.array([[1, 2j, -1j], [2, -1, 3j]], complex)
    coefficients /= np.linalg.norm(coefficients)
    actual = result["state"]@coefficients
    expected = result["target"]@coefficients
    np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=0)
    assert len(result["outcomes"]) == 2
    assert np.linalg.norm(result["outcomes"][0]["state"]) < 2e-11
    source._gram(result["state"], result["failures"])
    assert result["neutral_pointer_reads"] == 1
    assert result["resources"]["native_pulses"] > 2*1642+6226
    assert len(result["failures"]) == result["resources"]["native_pulses"]


@pytest.mark.parametrize("api", [source.prepare, source.transport, source.exchange,
                                source.transport_readout])
@pytest.mark.parametrize("corruption", ["event", "complement", "shape", "nan", "masked"])
def test_public_services_reject_missing_zero_outcomes_and_malformed_exports(monkeypatch, api, corruption):
    original = source.execute_with_exports
    def damaged(word, bits, initial):
        state, failures = original(word, bits, initial)
        if corruption == "event":
            return state, failures[:-1]
        events = list(failures)
        first = list(events[0])
        if corruption == "complement":
            first = first[:-1]
        elif corruption == "shape":
            first[0] = np.zeros((1,), complex)
        elif corruption == "nan":
            first[0] = np.full(np.asarray(first[0]).shape, np.nan)
        else:
            first[0] = np.ma.array(first[0], mask=False)
        events[0] = tuple(first)
        return state, tuple(events)
    monkeypatch.setattr(source, "execute_with_exports", damaged)
    with pytest.raises(ArithmeticError, match="export|failure"):
        api()


def test_stable_resource_record_keeps_integer_costs_and_admits_diagnostics():
    resources = compiler.resources(source.block_swap_word())
    counts = source._counts(resources)
    assert all(type(value) in (bool, int) for value in counts.values())
    assert counts["native_pulses"] == 416
    assert counts["integrated_processor_norm_upper_integer"] >= resources["integrated_processor_norm_upper"]
    resources["maximum_numerical_valid_code_leakage"] = float("nan")
    with pytest.raises(ArithmeticError, match="diagnostics"):
        source._counts(resources)
