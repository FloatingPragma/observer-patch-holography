from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from model.native_spin_observer import algebra
from model.protected_defect_catalog import charge
from . import energy, clock


def test_exact_record_blindness_includes_different_energies():
    row = energy.repair_blindness()
    assert row['complete_protected_basis_branches_checked'] == 12288
    assert row['two_protected_product_input_energies'] == ['0', '-6']
    assert row['both_inputs_have_identical_complete_repair_history_laws'] is True


def test_povm_span_obstruction_and_positive_score_identity():
    # Exact separating density matrices for an energy-blind complete POVM.
    first, second = np.diag([0.75, 0.25]), np.diag([0.25, 0.75])
    effects = [np.eye(2)/4]*4
    assert all(np.trace(first@f) == np.trace(second@f) for f in effects)
    assert np.trace(first@np.diag([1, -1])) != np.trace(second@np.diag([1, -1]))
    constant, coefficients = energy.physical_terms()
    reconstructed = np.zeros((256, 256), complex)
    settings = list(coefficients.items())+[('I'*8, F(0))]*20
    for label, coefficient in settings:
        p = algebra.tensor(label)
        for sign in (1, -1):
            effect = (np.eye(256)+sign*p)/128
            reconstructed += float(constant+64*coefficient*sign)*effect
    op = charge.operators()
    diagonal = np.array([(-1)**sum((mask >> j) & 1 for j, sign in enumerate(op['signs']) if sign == -1)
                         for mask in range(256)])
    physical = diagonal[:, None]*op['H']*diagonal[None, :]/op['scale']
    assert np.array_equal(reconstructed, physical)


def test_local_work_against_full_independent_ground_energy():
    op = charge.operators()
    h = op['H']/op['scale']
    eigenvalues, vectors = np.linalg.eigh(h)
    ground, e0 = vectors[:, 0], eigenvalues[0]
    value = energy.work()
    total = 0.0
    for row in value['channels']:
        a = op['twice_creation'][row['charge']]/2
        image = a@ground
        probability = float(np.vdot(image, image).real)
        moment = float(np.vdot(image, (h-e0*np.eye(256))@image).real)
        commutator = np.vdot(ground, a.T@(h@a-a@h)@ground).real
        assert abs(moment-commutator) < 2e-14
        assert float(row['probability_lower']) <= probability <= float(row['probability_upper'])
        assert float(row['unnormalized_excitation_lower']) <= moment <= float(row['unnormalized_excitation_upper'])
        assert float(row['mean_excitation_lower']) <= moment/probability <= float(row['mean_excitation_upper'])
        total += moment
    z = op['local_Z']
    assert abs(total-np.vdot(ground, (z@h@z-h)@ground).real) < 2e-14
    assert value['unitary_Z_regional_entropy_change_exact'] == 0
    assert value['unitary_Z_modular_relative_entropy_lower'] > F(17, 10)


def test_source_signed_port_action_uses_full_precision_and_orientation():
    rows = energy.primitive_action()['primitive_records']
    assert all(row['pulse_count'] == 1 and len(row['twelve_literal_controls']) == 12 for row in rows)
    for row in rows:
        lo, hi = map(F, row['signed_action_interval'])
        target = -F(row['angle']) if row['kind'] in ('Z', 'I') else F(0)
        error = F(row['full_generator_error'])
        assert lo-error <= target <= hi+error
        assert hi-lo < F(1, 10**40)
    assert [F(x) for x in rows[0]['signed_action_Q_sqrt5']] == [-F(x) for x in rows[1]['signed_action_Q_sqrt5']]
    assert rows[2]['signed_action_interval'] == ['0', '0']


def test_maintenance_action_and_zero_power_on_actual_created_state():
    op = charge.operators()
    h = op['H']/op['scale']
    e, v = np.linalg.eigh(h)
    psi = op['local_Z']@v[:, 0]
    energy0 = np.vdot(psi, h@psi).real
    for time in (0.1, 0.7, 2.0):
        evolved = expm(-1j*time*h)@psi
        assert abs(np.vdot(evolved, h@evolved).real-energy0) < 1e-13
    report = clock.maintenance(F(7, 3))
    by_charge = {row['charge']: row for row in energy.work()['channels']}
    for row in report['created_channel_action']:
        assert row['autonomous_system_power_exact'] == 0
        assert row['bare_protected_repair_energy_input_exact'] == 0
        assert row['conditional_relative_action_lower'] == F(7, 3)*by_charge[row['charge']]['mean_excitation_lower']


@pytest.mark.parametrize('time', [-1, F(-1), 0.2, True])
def test_maintenance_requires_declared_exact_nonnegative_clock(time):
    with pytest.raises(ValueError): clock.maintenance(time)
