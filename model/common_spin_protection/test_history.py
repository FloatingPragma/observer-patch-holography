"""Arbitrary-reference factor test, with every outside and native outcome."""
import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from . import history, spectrum, spin


def test_arbitrary_correlated_multiplicity_spin_reference_survives_joined_instrument():
    chart = spectrum.chart(4, 2)
    rng = np.random.default_rng(1141)
    tensor = rng.normal(size=(105, 3, 2))+1j*rng.normal(size=(105, 3, 2))
    tensor /= np.linalg.norm(tensor)
    state = sum(chart.complete_columns(j, original=True)@tensor[:, j, :] for j in range(3))
    expected = np.einsum('mar,mbs->arbs', tensor, tensor.conj()).reshape(6, 6)
    generator = spin.operators().K.matrix(spin.operators().scale)
    evolved = expm_multiply(-.37j*generator, state)
    observed = np.zeros_like(expected)
    masses = []
    for label, operator in spin.repair_outcomes('hole', 'E', (1, 2)):
        output = operator@evolved
        observed += history.spin_reference_marginal(output[:, :, None], chart)
        masses.append(float(np.vdot(output, output).real))
    np.testing.assert_allclose(observed, expected, atol=1e-14, rtol=0)
    assert abs(sum(masses)-1) < 1e-14
    # Discarding or success-normalizing the outside event would be material.
    assert masses[-1] > .9
    assert sum(masses[:-1]) < .1


def test_complete_two_repair_history_retains_small_carrier_leakage_and_reference():
    receipt = history.report()
    assert len(receipt['first_complete_probabilities']) == 17
    assert len(receipt['complete_joint_outcomes']) == 153
    assert receipt['first_outside_small_carrier_probability'] > .005
    assert receipt['final_spin_reference_Frobenius_error'] < 1e-13
    assert receipt['initial_spin_reference_purity'] == pytest.approx(1, abs=1e-13)
    assert receipt['spin_z_record_control_reference_trace_distance'] == pytest.approx(2/3, abs=1e-13)
    assert receipt['postselection_used'] is False
    assert receipt['physical_pulse_source_executed'] is False


@pytest.mark.parametrize('bad', [np.zeros((4096, 2)), np.zeros((10, 2, 1)),
                               np.zeros((4096, 0, 1)), np.ones((4096, 2, 1))*np.nan])
def test_malformed_complete_history_columns_reject(bad):
    with pytest.raises(ValueError):
        history.spin_reference_marginal(bad)
