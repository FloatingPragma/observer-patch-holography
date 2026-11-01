from fractions import Fraction as F
import copy

import numpy as np
import pytest

from model.interacting_vacuum import ground
from model.cube_eta_protection import cube
from model.field_symmetry.symmetry import onsite_frame
from model.protected_particle_readout import execute
from . import preparation as p, circuits as c


def test_source_derived_signed_weights_are_exact():
    weights, signs = p.signed_weights()
    assert sum(weights) == 1296
    assert sum(bool(w) for w in weights) == 62
    assert set(weights) == {0, 9, 27, 36, 81}
    assert all(not w or mask.bit_count() == 4 for mask, w in enumerate(weights))
    for mask, (a, b) in ground.reference().items():
        j = sum(((mask >> (2*x)) & 1) << x for x in range(8))
        assert weights[j] == a*a+3*b*b and signs[j]*(a+b) > 0


def test_signed_tree_certificate_and_negative_sign_control():
    word = p.tree_word(); weights, signs = p.signed_weights()
    checked = p.signed_tree_certificate(word, weights, signs)
    assert checked['conditional_nodes'] == 255
    assert checked['signed_tree_error']+checked['fixed_angle_error'] < F(1, 10**11)
    wrong = p.signed_tree_certificate(word, weights, (1,)*256)
    assert wrong['signed_tree_error'] > 1


def test_tape_mutation_and_remaining_event_are_detected():
    weights, signs = p.signed_weights(); word = p.tree_word()
    changed = copy.deepcopy(word); changed[0]['angle'] += .01
    result = p.signed_tree_certificate(changed, weights, signs)
    assert result['signed_tree_error']+result['fixed_angle_error'] > F(1, 1000)
    with pytest.raises(ArithmeticError): p.signed_tree_certificate(word+word[:1], weights, signs)
    with pytest.raises(ValueError): p.signed_tree_certificate(word, tuple(float(x) for x in weights), signs)
    with pytest.raises(ValueError): p.signed_tree_certificate(word, weights, (0,)*256)


def test_inverse_frame_full_four_columns_including_scalar_phases():
    assert p.exact_frame_certificate()['exact_local_Fock_columns']==32
    for site in cube.GEOMETRY.sites:
        word = p.emit(p.inverse_frame_plan(site))
        actual = execute.compile_block(execute.key(word), 2)
        expected = np.eye(4, dtype=complex)
        expected[1:3, 1:3] = onsite_frame(site).conj().T
        assert np.linalg.norm(actual.success-expected) < 2e-13
        assert actual.complete_error < 2e-13


def test_complete_native_ground_prep_no_normalization():
    tree, frame = p.runs(); target = p.ideal_ground()
    bound = F(p.certificate()['complete_preparation_state_error'])
    complete_error = np.sqrt(np.linalg.norm(frame.success[:, 0]-target)**2+
                             tree.failure_probability+frame.failure_probability)
    assert complete_error < float(bound)+1e-14
    assert abs(np.linalg.norm(frame.success)**2+tree.failure_probability+frame.failure_probability-1) < 1e-11
    assert p.certificate()['all_abort_outputs_retained']


def test_spectator_embedding_preserves_all_rows_and_reference_columns():
    rng = np.random.default_rng(724)
    vector = rng.normal(size=(256, 3))+1j*rng.normal(size=(256, 3))
    output = p.embed_blank_modes(vector)
    assert np.linalg.norm(output.conj().T@output-vector.conj().T@vector) < 1e-11
    assert np.count_nonzero(output) == np.count_nonzero(vector)
    with pytest.raises(ValueError): p.embed_blank_modes(np.zeros(128))


def test_car_majorana_and_missing_jw_negative_control():
    state = p.ideal_ground(); result = c.majorana(state)
    assert np.linalg.norm(c.majorana(result)-state) == 0
    assert np.linalg.norm(result[cube.weights(5)]) > 1-1e-14
    wrong = state[np.arange(65536)^2]
    assert abs(np.vdot(result, wrong)) < 1e-14
    assert abs(np.linalg.norm(result-wrong)-np.sqrt(2)) < 1e-14
    assert c.exact_probe_certificate()['native_qubit_parity_coherence_service_required']


def test_origin_creation_preserves_the_initial_exterior_and_number_flag():
    result=c.localized_defect_certificate()
    assert result['creation_is_unitary_and_origin_supported']
    assert result['initial_exterior_CAR_restriction_unchanged']
    assert (result['input_number'],result['created_number'])==(4,5)
    assert result['mean_added_energy']=='1/2'
    assert not result['persistent_spatial_localization_proved']


def test_controlled_majorana_all_eight_local_columns():
    word = p.emit(c.majorana_plan())
    actual = execute.compile_block(execute.key(word), 3)
    ideal = np.zeros((8, 8), complex)
    for mask in range(8):
        target = mask^2 if mask & 4 else mask
        ideal[target, mask] = (-1)**(mask & 1) if mask & 4 else 1
    assert np.linalg.norm(actual.success-ideal) < 1e-13
    assert actual.complete_error < 1e-13


def test_complete_pointer_prep_has_ground_and_charged_arms():
    ground_vector = p.ideal_ground()
    initial = np.zeros((131072, 1), complex); initial[:65536, 0] = ground_vector
    result = execute.complete(c.pointer_blocks(), 17, initial)
    target = np.concatenate((ground_vector, c.majorana(ground_vector)))[:, None]/np.sqrt(2)
    assert np.linalg.norm(result.success-target) < 1e-13
    assert np.linalg.norm(result.success.conj().T@result.success+result.failure_gram-np.eye(1)) < 1e-12


@pytest.mark.parametrize('phase', [0, .5, F(1, 3), '1/2', True])
def test_invalid_analyzer_phase_rejected(phase):
    with pytest.raises(ValueError): c.pointer_word(phase)
