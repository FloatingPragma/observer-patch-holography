"""Independent instrument and hostile controls for charged transport/readout."""
import numpy as np
import pytest

from . import source
from ..full_fock_field.field import FieldGeometry
from ..observer.source_scalar import compiler
from ..particle_response import dimer


def predecessor_projector(kind, axis=0):
    sides = tuple(2 if j == axis else 1 for j in range(3))
    spectrum = dimer.setup(FieldGeometry(sides))
    name = "charge_minus_one_hole" if kind == "hole" else "charge_plus_one_particle"
    return dimer.candidate_projectors(spectrum)[name]


@pytest.mark.parametrize("kind,axis", [("hole", 1), ("particle", 2)])
def test_full_source_presence_instrument_preserves_both_reference_branches(kind, axis):
    # The predecessor spectral calculation does not use the new encoder or its
    # analytic candidate frame. Nondefault axes exercise complex spin frames.
    projector = predecessor_projector(kind, axis)
    other = np.eye(16)-projector
    word = source.presence_readout_word(kind, axis)
    actual, failures = source.execute_with_exports(word, 5, np.eye(32, dtype=complex))
    expected = np.block([[other, projector], [projector, other]])
    np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=0)

    # A three-dimensional untouched reference is entangled with every field
    # occupation sector. Both readout branches have nonzero weight.
    rng = np.random.default_rng(721+axis)
    joint = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    joint /= np.linalg.norm(joint)
    outputs = (actual[:16, :16]@joint, actual[16:, :16]@joint)
    for output, effect in zip(outputs, (other, projector), strict=True):
        np.testing.assert_allclose(output, effect@joint, atol=2e-11, rtol=0)
        np.testing.assert_allclose(effect@output, output, atol=2e-11, rtol=0)
        assert np.linalg.norm(output) > .1

    # Complete-reference trace preservation includes each failed export column.
    gram = sum(output.conj().T@output for output in outputs)
    assert len(failures) == len(word)
    for event, branches in zip(word, failures, strict=True):
        dimension = 1 << len(event["targets"])
        assert len(branches) == 6-dimension
        for branch in branches:
            assert branch.shape == (32//dimension, 32)
            failed = branch[:, :16]@joint
            gram += failed.conj().T@failed
    np.testing.assert_allclose(gram, joint.conj().T@joint, atol=2e-11, rtol=0)

    # Completed readout commutes with total field number on both pointer arms.
    number = np.diag([j.bit_count() for j in range(16)])
    full_number = np.kron(np.eye(2), number)
    np.testing.assert_allclose(actual@full_number, full_number@actual,
                               atol=2e-11, rtol=0)


def test_complement_only_phase_passes_selected_transport_but_fails_full_certificate(monkeypatch):
    original = source.presence_readout_word

    def changed(kind="hole", axis=0):
        # Vacuum belongs to the rejected complement for either charged doublet.
        # Its phase cannot be detected by a transported candidate witness.
        phases = np.zeros(16)
        phases[0] = .17
        return original(kind, axis)+compiler.compile_diagonal(phases)

    monkeypatch.setattr(source, "presence_readout_word", changed)
    selected = source.transport_readout("hole")
    assert selected["error"] < 2e-11
    with pytest.raises(ArithmeticError, match="Stinespring"):
        source._presence_certificate("hole")


def test_correct_presence_probabilities_cannot_hide_internal_reference_rotation(monkeypatch):
    original = source.presence_readout_word
    # For axis0 this is a twisted internal Sz rotation. It commutes with the
    # lower-doublet projector and number, but changes its unknown internal state.
    spin_word = []
    for label, angle in (("ZIII", .17), ("IZII", -.17),
                         ("IIZI", -.17), ("IIIZ", .17)):
        spin_word += compiler.pauli_rotation(label, angle)
    spin = compiler.execute(spin_word, 4)
    projector = predecessor_projector("hole")
    np.testing.assert_allclose(spin@projector, projector@spin, atol=2e-11, rtol=0)
    assert np.linalg.norm((spin-np.eye(16))@projector) > .1

    def changed(kind="hole", axis=0):
        return original(kind, axis)+spin_word

    actual = compiler.execute(changed(), 5)
    # Both presence effects are still exactly correct despite the disturbance.
    for block, effect in ((actual[:16, :16], np.eye(16)-projector),
                          (actual[16:, :16], projector)):
        np.testing.assert_allclose(block.conj().T@block, effect, atol=2e-11, rtol=0)
    monkeypatch.setattr(source, "presence_readout_word", changed)
    with pytest.raises(ArithmeticError, match="Stinespring"):
        source.transport_readout("hole")


@pytest.mark.parametrize("kind", ["hole", "particle"])
def test_canonical_join_counts_preparation_transport_and_complete_readout(kind):
    result = source.transport_readout(kind)
    resources = result["resources"]
    # 1642 charged +4168 ground +416 route +3620 complete presence read.
    assert resources["native_pulses"] == 9846
    assert resources["proper_code_load_events"] == 9846
    assert resources["proper_code_export_events"] == 9846
    assert resources["proper_code_service_events"] == 19692
    assert resources["maximum_work_qubits"] == 9
    assert len(result["failures"]) == 9846
    assert sum(map(len, result["failures"])) == 34336
    assert result["first_failure_kraus_slots"] == 34336
    assert [row["pointer"] for row in result["outcomes"]] == [0, 1]
    assert result["neutral_pointer_bits_initialized"] == 1
    assert result["neutral_pointer_reads"] == 1
    assert result["supplied_occupation_register_bits"] == 8
    assert result["initialized_field_number"] == (3 if kind == "hole" else 5)
