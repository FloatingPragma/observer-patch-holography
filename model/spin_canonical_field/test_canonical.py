from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import algebra, bounds, packet


def unitary(matrix, time=1):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors*np.exp(1j*time*values))@vectors.conj().T


def matrices():
    spin = [np.asarray(a, complex) for a in algebra.spin()]
    eye = np.eye(2)
    h = sum(np.kron(a, a) for a in spin)-np.eye(4)/4
    return spin, eye, h+np.eye(4)


def test_exact_all_axis_currents_and_canonical_orientation():
    result = algebra.exact_checks()
    assert result['all_three_axis_double_commutators']
    assert result['second_transverse_momentum_commutator_is_nonzero']


@pytest.mark.parametrize('a,b', [(F(1, 50), F(3, 100)),
                                (F(-2, 5), F(1, 3)), (F(3), F(-2))])
def test_finite_weyl_remainder_keeps_orientation_and_both_terms(a, b):
    spin, _, _ = matrices()
    u, v = unitary(spin[1], float(a)), unitary(spin[2], float(b))
    actual = u@v@u.conj().T@v.conj().T
    target = unitary(spin[0], -float(a*b))
    assert np.linalg.norm(actual-target, 2) <= float(bounds.site_weyl_error(a, b))+2e-14
    if abs(a)+abs(b) < F(1, 10):
        wrong_sign = unitary(spin[0], float(a*b))
        assert np.linalg.norm(actual-wrong_sign, 2) > float(bounds.site_weyl_error(a, b))


@pytest.mark.parametrize('axis', [0, 1, 2])
def test_energy_kick_controls_arbitrary_vectors_not_only_ground(axis):
    spin, eye, k = matrices()
    left, right = F(2, 5), F(-3, 7)
    a = float(left)*np.kron(spin[axis], eye)+float(right)*np.kron(eye, spin[axis])
    u = unitary(a)
    delta2 = bounds.energy_kick(1, (left-right)**2)
    vectors = [np.array([1, 0, 0, 0], complex),
               np.array([1, 2j, -3, 4j], complex),
               np.array([0, 1, -1, 0], complex)]
    for vector in vectors:
        vector /= np.linalg.norm(vector)
        energy = float(np.vdot(vector, k@vector).real)
        moved = u@vector
        actual = float(np.vdot(moved, k@moved).real)
        assert actual <= (np.sqrt(max(0, energy))+np.sqrt(float(delta2)))**2+2e-14


def test_transport_comparison_has_a_nontrivial_finite_control():
    spin, eye, k = matrices()
    a = .07*np.kron(spin[0], eye)-.03*np.kron(eye, spin[0])
    v = unitary(a)
    kp = v.conj().T@k@v
    first = k@a-a@k
    nested = a@first-first@a
    delta2 = np.linalg.norm(nested, 2)/2
    psi = np.array([1, 2j, -3, 4j], complex)
    psi /= np.linalg.norm(psi)
    energy = np.vdot(psi, k@psi).real
    time = .3
    actual = np.linalg.norm((unitary(k, time)-unitary(kp, time))@psi)**2
    bound = 4*time*np.sqrt(delta2*energy)+10*time*delta2
    assert 0 < actual <= bound < 1


def test_positive_order_alone_does_not_force_a_scalar_center():
    # Exact mixture of the two product Neel states: finite averages are
    # already +/-1/2 on the two sectors, at every volume.
    plus, minus = F(3, 4), F(1, 4)
    mean = (plus-minus)/2
    variance = (plus+minus)/4-mean**2
    assert mean == F(1, 4) and variance == F(3, 16)
    assert mean < F(1, 2)  # This state does not attain the maximal order.


def test_transport_examples_have_actual_decay_without_claiming_order_rate():
    rows = [bounds.scaling_example(n) for n in (2, 4, 8, 16, 64)]
    values = [F(row['transport_distance_squared_upper']) for row in rows]
    assert all(a > b for a, b in zip(values, values[1:]))
    assert values[-1] < F(1, 1000)
    assert all(row['order_self_averaging_rate_claimed'] is False for row in rows)


@pytest.mark.parametrize('value', [F(0), F(1, 3), F(2), F(7, 10**80)])
def test_outward_root(value):
    result = bounds.sqrt_upper(value)
    assert result**2 >= value
    assert result == 0 or (result-F(1, 2**80))**2 <= value


@pytest.mark.parametrize('bad', [True, 1.0, -1, '1'])
def test_nonrational_or_negative_bounds_rejected(bad):
    with pytest.raises(ValueError):
        bounds.transport_error_squared(bad, 1, 1)


def test_second_moment_and_energy_constants_are_distinct():
    result = bounds.constants()
    assert result['uniform_second_moment_coefficient'] == '1/4'
    assert result['staggered_second_moment_coefficient'] == '3'
    assert result['source_L1_cone_speed_upper'] == '72'
    assert bounds.energy_after(F(4), F(1)) == 9


def test_receipt_rejects_species_shell_and_domain_promotions():
    value = packet.build()
    packet.validate(value)
    for key in ('sharp_mass_shell_proved', 'physical_species_identified',
                'surviving_connected_interaction_proved',
                'common_Wightman_field_domain_proved',
                'two_independent_canonical_pairs_proved'):
        altered = deepcopy(value)
        altered['theorem'][key] = True
        with pytest.raises(ValueError):
            packet.validate(altered)


def test_boolean_alias_cannot_launder_an_unproved_scope():
    value = packet.build()
    value['theorem']['sharp_mass_shell_proved'] = 0
    with pytest.raises(ValueError):
        packet.validate(value)


def test_source_change_during_replay_is_rejected(monkeypatch):
    original = packet.bindings()
    calls = 0
    def changed():
        nonlocal calls
        calls += 1
        result = dict(original)
        if calls > 1:
            result[next(iter(result))] = '0'*64
        return result
    monkeypatch.setattr(packet, 'bindings', changed)
    with pytest.raises(RuntimeError, match='source changed'):
        packet.build()
