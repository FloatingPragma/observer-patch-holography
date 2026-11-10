"""Exact invariants and nontrivial complete detector controls of the joint clock."""
from fractions import Fraction as F
import pytest
import sympy as s
from .model import JointObserver, SparseState, combine
from .certificate import interaction_certificate, structure_certificate


def test_joint_controls_and_trace():
    result = structure_certificate()
    assert result['dimension'] == 72000
    assert result['unitary_trace'] == '-47850'
    assert result['all_30_link_signs_used_in_wave_repair']


def test_executed_reciprocal_interaction_and_coherence():
    result = interaction_certificate()
    assert result['wave_to_defect']['deterministic_final_face_indices'] == [17, 0]
    assert result['defect_to_wave']['yes_and_complement_probabilities'] == [['1/120', '119/120'], ['0', '1']]
    assert result['coherence']['coherent_inverse_controls_return_and_complement'] == ['1', '0']
    assert result['coherence']['after_complete_face_readout_return_and_complement'] == ['1/2', '1/2']


def test_multilink_q_matches_literal_central_sign_on_each_address():
    observer = JointObserver()
    # Exhaust all60 orientedaddresses and10 actualfields; both code charges.
    for word in observer.link_words:
        for address, (a, b) in enumerate(observer.addresses):
            for internal in (0, 3):
                initial = observer.basis_state(word, a, internal, 0, address)
                output = observer.repair(initial)
                sign = -1 if internal >= 3 and word >> (address//2) & 1 else 1
                assert output.as_dict() == {
                    (word, a, internal, 0, address): F(1, 2),
                    (word, b, internal, 0, address): F(sign, 2),
                    (word, a, internal, 1, address): F(1, 2),
                    (word, b, internal, 1, address): F(sign, 2),
                }
                assert observer.repair(output) == initial


def test_complex_mixed_charge_state_preserves_norm_inverse_and_charges():
    observer = JointObserver()
    state = SparseState({(1, 0, 0, 0, 0): F(1), (3, 5, 3, 1, 4): s.I,
                         (6553635, 11, 4, 0, 59): F(2)}, 6)
    output = observer.step(state)
    assert state.norm_squared() == output.norm_squared() == 1
    assert observer.inverse_step(output) == state
    for internal in range(5):
        before = SparseState({k: v for k, v in state.amplitudes if k[2] == internal}, 6)
        after = SparseState({k: v for k, v in output.amplitudes if k[2] == internal}, 6)
        assert before.norm_squared() == after.norm_squared()
    assert all(observer.law.region_charge(key[0]) == 1 for key, _ in output.amplitudes)
    assert sum(branch.norm_squared() for _, branch in observer.region_face_measurement(output)) == 1


def test_sparse_input_validation():
    observer = JointObserver()
    with pytest.raises(ValueError):
        SparseState({(1, 0, 0, 0, 0): 0.5})
    with pytest.raises(ValueError):
        SparseState({}, 0)
    with pytest.raises(ValueError):
        observer.step(SparseState({(0, 0, 0, 0, 0): 1}))
    with pytest.raises(ValueError):
        observer.embedded_ports(1, (0, 0))
    with pytest.raises(ValueError):
        observer.defect_control(observer.embedded_port(1, 0), True)
    with pytest.raises(ValueError):
        combine(SparseState({}, 1), SparseState({}, 2))


def test_constrained_complete_return_apparatus_and_recording_contrast():
    from .certificate import invariant_coherence_certificate
    from .gauss import covariant_preparation_return, gamma_amplitudes
    result = invariant_coherence_certificate()
    assert result['coherent_return_and_complement'] == ['1', '0']
    assert result['after_complete_face_readout_return_and_complement'] == ['1/2', '1/2']
    assert result['return_projector_commutes_with_Gamma']
    observer = JointObserver()
    state = SparseState({(1, 5, 3, 0, 0): s.I, (6553635, 6, 3, 0, 11): F(2),
                         (3, 0, 0, 1, 8): F(1)}, 6)
    branches = covariant_preparation_return(observer, state, (5, 6))
    assert sum(branch.norm_squared() for _, branch in branches) == 1
    for index, (_, branch) in enumerate(branches):
        again = covariant_preparation_return(observer, branch, (5, 6))
        assert again[index][1] == branch
        assert again[1-index][1].norm_squared() == 0
    transformed = SparseState(gamma_amplitudes(state.as_dict()), 6)
    gamma_branches = covariant_preparation_return(observer, transformed, (5, 6))
    for (_, branch), (_, other) in zip(branches, gamma_branches):
        assert SparseState(gamma_amplitudes(branch.as_dict()), 6) == other
