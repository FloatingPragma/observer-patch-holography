"""Independent spectral determinant and complete original-field Ramsey checks."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from ..common_spin_protection import spin
from . import exact, experiment as e


def test_all_fifteen_lines_and_thirteen_energies_are_visible_without_thresholds():
    value = exact.spectral_certificate()
    rows, groups = value['all15_mode_pairs'], value['all13_distinct_energy_groups']
    assert value['exact_mean_H_energy'] == '1'
    assert value['exact_H_energy_variance'] == '1/256'
    assert len(rows) == 15 and len(groups) == 13
    assert all(row['weight'] > 0 for row in rows)
    assert sum(row['weight'] for row in rows) == pytest.approx(1, abs=2e-15)
    central = next(row for row in groups if row['group'] == 'central')
    assert central['mode_pairs'] == [[1, 6], [2, 5], [3, 4]]
    assert central['H_energy'] == pytest.approx(1)
    assert all(row['H_energy']-row['K_energy'] == pytest.approx(2) for row in rows)


@pytest.mark.parametrize('time', [0., .5, 4., 17.])
def test_full_fock_spectral_and_independent_two_by_two_slater_determinant(time):
    state = exact.target()
    h = spin.operators().H.matrix(spin.operators().scale)
    direct = np.vdot(state, expm_multiply(-1j*time*h, state))
    adjacent = np.diag(np.ones(5), 1)+np.diag(np.ones(5), -1)
    single = expm(-1j*time*adjacent/16)
    determinant = np.exp(-1j*time)*np.linalg.det(single[:2, :2])
    assert direct == pytest.approx(exact.amplitude(time), abs=3e-15)
    assert determinant == pytest.approx(direct, abs=3e-15)


@pytest.fixture(scope='module')
def history():
    return e.report()


def test_complete_native_history_keeps_two_pointer_and_all_other_outputs(history):
    assert len(history['all40_native_repair_records']) == 40
    assert len(history['common_prefix_source_records']) == 3
    assert history['success_state_renormalized'] is False
    for setting in history['pointer_settings']:
        assert len(setting['pointer_records']) == len(setting['pointer_probabilities']) == 6
        assert setting['complete_probability'] == pytest.approx(1, abs=1e-11)
        assert setting['complete_output_reference_error_diagnostic'] < 3e-7
        assert setting['pointer_probabilities'][0] == pytest.approx(setting['analytic_pointer_zero_probability'], abs=3e-7)
    assert history['pointer_settings'][0]['pointer_probabilities'][0] > .93
    assert .73 < history['pointer_settings'][1]['pointer_probabilities'][0] < .75


def test_source_specification_and_strang_errors_have_separate_exact_budgets(history):
    certificate = e.source_certificate()
    assert certificate['Hamiltonian_duration'] == F(1, 2)
    assert certificate['each_field_segment']['pulses'] == 3360
    assert float(certificate['preparation']['complete_source_operator_error']) < 2e-13
    for bound, setting in zip(certificate['complete_output_error_bounds'], history['pointer_settings'], strict=True):
        assert type(bound) is F
        assert setting['complete_output_reference_error_diagnostic'] < float(bound) < .0073
    # Omitting the real field produces a completely different quadrature.
    assert abs(history['pointer_settings'][1]['pointer_probabilities'][0]-.5) > .2
