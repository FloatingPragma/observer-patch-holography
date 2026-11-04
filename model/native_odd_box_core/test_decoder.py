"""Native full-CAR parity, complete instrument and resource checks."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from ..native_spin_observer import circuits, execute
from . import decoder


def test_crossing_Givens_keeps_spectator_parity_and_double_occupations():
    modes = [decoder._annihilator(4, mode) for mode in range(4)]
    generator = modes[0].T@modes[3]-modes[3].T@modes[0]
    target = expm(np.pi*generator/4)/4
    template = decoder.givens(4, 0, 3, F(1, 4))
    run = execute.complete(circuits.emit(template), 4, np.eye(16)/4)
    assert np.linalg.norm(run.success-target) < 1e-12
    assert run.failure_probability < 1e-24
    # The same one-particle rotation without the intervening CAR parity has
    # opposite action on a sector with one occupied spectator.
    assert abs(target[3, 10]+target[1, 8]) < 1e-12
    assert abs(target[3, 10]-target[1, 8]) > .3
    assert abs(target[9, 9]-.25) < 1e-12


def test_chain_four_outcome_projectors_form_a_complete_QND_instrument():
    initial, target = decoder.chain_target()
    projectors = [8*target[j*64:(j+1)*64] for j in range(4)]
    assert np.linalg.norm(sum(projectors)-np.eye(64)) < 1e-12
    for j, p in enumerate(projectors):
        assert np.linalg.norm(p-p.conj().T) < 1e-12
        assert np.linalg.norm(p@p-p) < 1e-12
        for k, q in enumerate(projectors):
            if j != k:
                assert np.linalg.norm(p@q) < 1e-12
    assert np.linalg.norm(target.conj().T@target-initial.conj().T@initial) < 1e-12


def test_complete_native_chain_readout_with_full_quantum_reference():
    result = decoder.native_execution()
    assert result['quantum_reference_columns'] == 64
    assert result['native_first_failure_slots'] == 980
    assert result['complete_output_error_diagnostic'] < 2e-11
    assert result['complete_first_failure_probability'] < 1e-24
    assert not result['success_postselected_or_renormalized']
    assert not result['reference_columns_dephased']


def test_exact_cube_word_and_refinable_commands_have_complete_resource_bounds():
    value = decoder.certificate()
    assert value['cube_Givens_rotations_per_decode'] == 108
    assert value['cube_decode_native_pulses'] == 12312
    assert value['cube_complete_four_outcome_readout_native_pulses'] == 24630
    assert value['refined_dyadic_fraction_bits'] == 51
    assert F(value['refined_complete_unhalved_trace_error_upper']) < F(1, 10**9)
    for command in value['distinct_refined_commands'].values():
        assert len(command['integrated_port_coefficients']) == 12
        for coefficient in command['integrated_port_coefficients']:
            q = F(coefficient)
            assert q.denominator & (q.denominator-1) == 0
        assert F(command['generator_and_complete_output_error']) <= F(29, 2**51)
    assert not value['cube_full_Fock_native_readout_executed']
    assert not value['interleaved_native_repair_protection_proved']


@pytest.mark.parametrize('arguments', [(4, 2, 1, F(1, 4)), (4, 0, 4, F(1, 4)),
                                      (True, 0, 1, F(1, 4)), (4, 0, 1, .25)])
def test_invalid_native_CAR_targets_or_inexact_angles_rejected(arguments):
    with pytest.raises(ValueError):
        decoder.givens(*arguments)
