"""Independent source controls, complete colored outcomes and phase checks."""
from copy import deepcopy
import math

import numpy as np
import pytest

from . import dynamics as d


def test_controlled_swap_complete_truth_table_and_encoded_exchange():
    certificate = d.encoder_certificate()
    assert certificate['controlled_swap_full_input_columns'] == 8
    matrices, _, _ = d.prepare_internal(0, 1, 2)
    np.testing.assert_allclose(matrices[1], -matrices[0].T, atol=2e-9, rtol=0)
    assert abs(sum(np.linalg.norm(x)**2 for x in matrices)-1) < 2e-9
    # An orientation flip alone acts as -SWAP after decoding E_col.
    changed = matrices[::-1]
    decoded = (changed[0]-changed[1].T)/math.sqrt(2)
    np.testing.assert_allclose(decoded, -math.sqrt(2)*matrices[0].T, atol=2e-9, rtol=0)


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_full_source_phase_gate_retains_vacuum_and_padding(current):
    energies = np.array([.2, .7, 1.1, 2.])
    charges = np.array([-1., -1., 1., 1.])
    U, word = d.phase_gate(energies, charges, .5, .125, current, vacuum_energy=-3.)
    phase = -.5*(1+current/4)*(energies-1.5)+.125*current*charges
    expected = np.ones(128, complex)
    expected[:4] = np.exp(1j*phase)
    np.testing.assert_allclose(U, np.diag(expected), atol=2e-9, rtol=0)
    assert any(event['kind'] == 'I' for event in word)
    # Tensoring the two one-particle gates includes exactly one E_vac.
    target = np.exp(-.5j*(1+current/4)*(energies[0]+energies[2]-3.))
    assert abs(U[0, 0]*U[2, 2]-target) < 2e-9


def test_factor_contraction_matches_full_dense_state_without_truncation():
    rng = np.random.default_rng(1093)
    vectors = rng.normal(size=(3, 4))+1j*rng.normal(size=(3, 4))
    matrices = rng.normal(size=(4, 128, 128))+1j*rng.normal(size=(4, 128, 128))
    dense = sum(np.einsum('s,ab->sab', vectors[:, i], matrices[i]) for i in range(4))
    np.testing.assert_allclose(d.pointer_probabilities(vectors, matrices),
                               np.sum(abs(dense)**2, axis=0), atol=1e-10, rtol=0)
    assert abs(d.tensor_norm_squared(vectors, matrices)-np.linalg.norm(dense)**2) < 1e-8


def test_complete_charge_probe_distinguishes_all_seven_charges_and_preserves_them():
    charges = np.array([-1., -2/3, -1/3, 0., 1/3, 2/3, 1.])
    actual, _, certificate = d.charge_probe(charges)
    assert certificate['source_full_input_columns'] == 256
    assert len(certificate['charge_values']) == 7
    for index, charge in enumerate(charges):
        output = actual[:, index]
        support = output.copy()
        support[[index, 128+index]] = 0
        assert np.linalg.norm(support) < 2e-9
        assert abs(abs(output[index])**2-(1-math.sin(charge))/2) < 2e-9
        assert abs(np.linalg.norm(output)-1) < 2e-9
    values = [row['probability'] for row in certificate['pointer_zero_probability_by_charge']]
    assert all(a > b for a, b in zip(values, values[1:]))


@pytest.fixture(scope='module')
def executed():
    return d.execute([0., .2, .5, .9], [-1., -1., 1., 1.],
                     ramsey_indices=(0, 1), partner=2, vacuum_energy=-2.)


def test_every_native_and_internal_pointer_outcome_retained(executed):
    report, arrays = executed
    bound = report['ideal_contracted_source_factors']
    assert bound['decoder_factors']+bound['internal_pointer_copies'] == 23
    delta = math.expm1(23*math.log1p(d.ADMISSION))
    assert bound['normalized_input_operator_error_below'] == delta
    assert bound['normalized_input_probability_TV_error_below'] == delta+delta**2/2
    assert delta < 4.601e-6
    for law, count in [('C', 8), ('E', 16)]:
        assert len(report['laws'][law]) == count
        total = 0.
        for row in report['laws'][law]:
            key = row['array_prefix']
            P = arrays[key+'_all_pointer_probabilities']
            V = arrays[key+'_source_factors']
            M = arrays[key+'_internal_factors']
            assert P.shape == (2, 2, 128, 128)
            assert V.shape == (2, 512, 4) and M.shape == (2, 4, 128, 128)
            assert np.min(P) >= 0
            for orientation in (0, 1):
                for charge in (0, 1):
                    np.testing.assert_allclose(P[orientation, charge],
                                               d.pointer_probabilities(V[orientation], M[charge]),
                                               atol=1e-12, rtol=0)
            assert abs(P.sum()-row['ideal_native_weight']) < 2e-8
            assert P[1].sum() < 2e-16  # retained, not projected away
            total += P.sum()
        assert abs(total-1) < 2e-8


def test_measured_ramsey_gap_is_native_record_conditioned(executed):
    report, arrays = executed
    values = {}
    for row in report['laws']['C']:
        if not row['ideal_native_weight']:
            continue
        P = arrays[row['array_prefix']+'_all_pointer_probabilities'][0].sum(axis=0)
        selected = P[0].sum() if row['internal_Ramsey_slot'] == 0 else P[:, 0].sum()
        expected = .5*math.cos(.2*row['effective_duration']/2)**2
        assert abs(selected-expected) < 2e-8
        charge_zero = arrays[row['array_prefix']+'_all_pointer_probabilities'][:, 0].sum()
        assert abs(charge_zero-.5*(1+math.sin(1.))/2) < 2e-8
        values[row['current']] = selected
    assert set(values) == {0, 1} and abs(values[0]-values[1]) > 1e-3
    assert report['record_infers_energy_difference_not_absolute_rest_energy']
    assert not report['physical_species_identification_complete']


def test_complete_output_phase_mutation_is_detectable_even_with_same_probabilities(executed):
    from .dynamics_packet import validate_arrays, compare_arrays
    report, arrays = executed
    validate_arrays(report, arrays)
    altered = {k: v.copy() for k, v in arrays.items()}
    key = next(row['array_prefix'] for row in report['laws']['C'] if row['ideal_native_weight'])
    # Change a quantum output phase while retaining its saved probabilities.
    altered[key+'_internal_factors'][0, 0] *= 1j
    with pytest.raises(ValueError):
        compare_arrays(altered, arrays)
    missing = dict(arrays)
    missing.pop(key+'_all_pointer_probabilities')
    with pytest.raises(ValueError):
        validate_arrays(report, missing)


@pytest.mark.parametrize('energies,charges', [([0., -1.], [0., 0.]),
    ([0., float('nan')], [0., 0.]), ([0., 1.], [0.]), ([True, False], [0., 0.])])
def test_invalid_spectral_inputs_are_rejected(energies, charges):
    with pytest.raises(ValueError):
        d.phase_word(energies, charges, 1, 0, 0)
