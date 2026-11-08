"""Independent exact two-spin state, clocks and complete projector controls."""
import sympy as s


def controls():
    eye = s.eye(2)
    x = s.Matrix([[0, 1], [1, 0]])
    y = s.Matrix([[0, -s.I], [s.I, 0]])
    z = s.diag(1, -1)
    h0 = sum((s.kronecker_product(a, a) for a in (x, y, z)), s.zeros(4))/4-s.eye(4)/4
    selection = (s.kronecker_product(x, eye)-s.kronecker_product(eye, x))/2
    b = s.Rational(3, 8)
    hb = h0-b*selection
    vector = s.Matrix([-1, 3, -3, 1])/s.sqrt(20)
    rho = vector*vector.H
    singlet = s.Matrix([0, 1, -1, 0])/s.sqrt(2)
    ground = singlet*singlet.H
    excited = s.eye(4)-ground
    assert (vector.H*vector)[0] == 1
    assert hb*vector == -s.Rational(9, 8)*vector
    assert hb.eigenvals() == {-s.Rational(9, 8): 1, s.Rational(1, 8): 1, s.Integer(0): 2}
    assert h0 == -ground
    assert ground+excited == s.eye(4)
    weights = [s.trace(projector*rho) for projector in (ground, excited)]
    assert weights == [s.Rational(9, 10), s.Rational(1, 10)]
    assert h0*rho-rho*h0 != s.zeros(4)
    assert hb*rho-rho*hb == s.zeros(4)
    # At t=pi, exp(-itH0)=I-2 P_ground. No floating exponential is used.
    zero_evolution = s.eye(4)-2*ground
    rho_after = zero_evolution*rho*zero_evolution.H
    return_effects = [rho, s.eye(4)-rho]
    before = [s.trace(e*rho) for e in return_effects]
    after = [s.trace(e*rho_after) for e in return_effects]
    assert before == [1, 0] and after == [s.Rational(16, 25), s.Rational(9, 25)]
    delta = rho_after-rho
    assert delta.eigenvals() == {s.Rational(3, 5): 1, -s.Rational(3, 5): 1, s.Integer(0): 2}
    variance = s.trace(h0*h0*rho)-s.trace(h0*rho)**2
    assert variance == s.Rational(9, 100)
    # Torque is checked on the actual source spin and selected field.
    a = s.kronecker_product(y, eye)/2
    torque = -s.I*b*s.kronecker_product(z, eye)/2
    assert (hb*a-a*hb)-(h0*a-a*h0) == torque
    return {'field': '3/8', 'sites': 2, 'bonds': 1,
            'HB_ground_energy': '-9/8', 'HB_first_gap': '9/8',
            'HJ_ground_energy': '-1', 'HJ_first_gap': '1',
            'HJ_excitation_energies': ['0', '1'],
            'complete_HJ_energy_probabilities': [str(w) for w in weights],
            'HJ_excitation_mean': '1/10', 'HJ_energy_variance': str(variance),
            'HB_ground_is_HJ_stationary': False,
            'return_time': 'pi', 'HB_return_probabilities': [str(w) for w in before],
            'HJ_return_probabilities': [str(w) for w in after],
            'complete_return_record_TV': '9/25',
            'full_state_halftrace_at_pi': '3/5',
            'selection_torque_exact': True,
            'source_pulse_tape_executed_here': False}
