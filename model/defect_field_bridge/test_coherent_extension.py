"""Independent domain, coherence, and complete source checks."""
import math

import numpy as np
import pytest

from . import coherent_extension as coherent, full_domain


@pytest.mark.parametrize('law', ['C', 'E'])
def test_all_binary_columns_equal_old_law_and_nonbinary_columns_hold_coherently(law):
    for mask in range(256):
        binary = (mask & 3).bit_count() <= 1 and ((mask >> 6) & 3).bit_count() <= 1
        for coin in (0, 1):
            actual = coherent.branch_action(mask, law, coin, sites=4, flavors=2, seam=(0, 3))
            if binary:
                assert actual == full_domain.branch_action(mask, law, coin, sites=4, flavors=2, seam=(0, 3))
            else:
                assert actual == (0, mask, 1, 0)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_all_sectors_and_unital_effect_sums(law):
    rows = coherent.finite_branches(sites=4, flavors=2, law=law, seam=(0, 3))
    inputs, outputs = np.zeros(256, int), np.zeros(256, int)
    for row in rows.values():
        assert len({x[1] for x in row}) == len(row)
        for origin, target, sign in row:
            inputs[origin] += 1
            outputs[target] += 1
            assert origin.bit_count() == target.bit_count()
            assert sign in (-1, 1)
    assert np.array_equal(inputs, np.full(256, 2))
    assert np.array_equal(outputs, np.full(256, 2))
    assert not rows[0, 8] and not rows[1, 8]


@pytest.mark.parametrize('law', ['C', 'E'])
def test_coherent_control_preserves_vacuum_pair_superposition_flagged_control_dephases(law):
    psi = np.zeros(16, complex)
    psi[0], psi[3] = 1/math.sqrt(2), 1j/math.sqrt(2)
    rho = np.outer(psi, psi.conj())
    actual = coherent.channel(rho, sites=2, flavors=2, law=law, seam=(0, 1))
    flagged = full_domain.channel(rho, sites=2, flavors=2, law=law, seam=(0, 1))
    assert np.array_equal(actual, rho)
    assert flagged[0, 3] == 0
    assert abs(np.trace(flagged@flagged).real-.5) < 1e-14


@pytest.mark.parametrize('law', ['C', 'E'])
def test_identity_on_even_endpoint_subspace_with_arbitrary_reference(law):
    rng = np.random.default_rng(114001)
    support = coherent.even_endpoint_identity(2, 2, law, (0, 1))
    value = np.zeros((16, 3), complex)
    value[list(support)] = rng.normal(size=(len(support), 3))+1j*rng.normal(size=(len(support), 3))
    value /= np.linalg.norm(value)
    actual, _, _ = coherent.source_isometry(2, 2, law)
    evolved = actual@value
    # Complete output ordering field+16coin+32raw. No outcome is removed.
    expected = np.zeros((512, 3), complex)
    expected[:16] = expected[16:32] = value/math.sqrt(2)
    assert np.linalg.norm(evolved-expected) < coherent.SOURCE_BOUND


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('sites,flavors', [(2, 2), (3, 1)])
def test_source_dilation_retains_every_input_output_and_failure_record(sites, flavors, law):
    actual, target, resources = coherent.source_isometry(sites, flavors, law)
    assert actual.shape == (512, 1 << (sites*flavors))
    assert np.linalg.norm(actual-target) < coherent.SOURCE_BOUND
    assert np.linalg.norm(actual.conj().T@actual-np.eye(actual.shape[1])) < coherent.SOURCE_BOUND
    assert resources['proper_code_complement_outcome_slots'] == (
        4*resources['single_code_pulses']+2*resources['pair_code_pulses'])
    assert resources['proper_code_complement_outcome_slots'] > 0


def test_design_does_not_remove_backaction_for_a_moving_binary_state():
    psi = np.zeros(16, complex)
    psi[1] = psi[3] = 1/math.sqrt(2)
    rho = np.outer(psi, psi.conj())
    result = coherent.channel(rho, sites=2, flavors=2, law='C', seam=(0, 1))
    assert np.linalg.norm(result-rho) > .1
    assert abs(result[1, 3]-.25) < 1e-14


def test_source_inverse_drive_echo_returns_one_for_both_laws():
    for row in coherent.source_echo().values():
        assert row['all_record_field_amplitudes'].shape == (32, 16)
        assert abs(row['return_probability']-1) < 4*coherent.SOURCE_BOUND
