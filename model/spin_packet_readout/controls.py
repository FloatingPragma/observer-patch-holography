"""Exact all-output LCU interference and retained-reference countercontrols."""
import sympy as s


def exact_control():
    eye = s.eye(2)
    x = s.Matrix([[0, 1], [1, 0]])
    z = s.diag(1, -1)
    unitaries = (eye, x, z, x*z)
    phases = (1, s.I, -1, -s.I)
    kraus = []
    for outcome in range(4):
        kraus.append(sum(((-1)**((outcome & term).bit_count())*phases[term]*unitaries[term]
                          for term in range(4)), s.zeros(2))/4)
    complete = sum((k.adjoint()*k for k in kraus), s.zeros(2))
    if complete != eye:
        raise ArithmeticError('a rejected output was lost')
    target = (eye+s.I*x-z-s.I*x*z)/4
    if kraus[0] != target:
        raise ArithmeticError('blank-ancilla interference has the wrong phase')
    probabilities = [s.simplify(s.trace(k.adjoint()*k)/2) for k in kraus]
    retained_reference = sum(((k.adjoint()*k).T/2 for k in kraus), s.zeros(2))
    if retained_reference != eye/2 or sum(probabilities) != 1:
        raise ArithmeticError('the complete instrument changed an untouched reference')
    probe_effects = ((eye+z)/2, (eye-z)/2)
    refined = sum((k.adjoint()*effect*k for k in kraus for effect in probe_effects), s.zeros(2))
    if refined != eye:
        raise ArithmeticError('composition lost an instrument output')
    evolution = s.diag(1, s.I)
    conjugated = evolution*x*evolution.adjoint()
    shifted = (-evolution)*x*(-evolution).adjoint()
    if shifted != conjugated or evolution*x == conjugated:
        raise ArithmeticError('scalar-clock cancellation negative control failed')
    wrong = (eye-s.I*x-z+s.I*x*z)/4
    if wrong == target:
        raise ArithmeticError('the filter phase-sign countercontrol disappeared')
    bell_plus = s.Matrix([1, 0, 0, 1])/s.sqrt(2)
    bell_minus = s.Matrix([1, 0, 0, -1])/s.sqrt(2)
    bell_difference = bell_plus*bell_plus.adjoint()-bell_minus*bell_minus.adjoint()
    marginal_difference = s.Matrix(2, 2, lambda i, j:
        sum(bell_difference[2*i+r, 2*j+r] for r in range(2)))
    joint_distance = sum(abs(e)*multiplicity
                         for e, multiplicity in bell_difference.eigenvals().items())
    if marginal_difference != s.zeros(2) or joint_distance != 2:
        raise ArithmeticError('marginal distance cannot bound unknown reference correlations')
    return {'all_four_Kraus_effects_sum_to_identity': True,
            'reference_marginal_exactly_preserved': True,
            'subsequent_two_bin_probe_all_eight_effects_complete': True,
            'maximally_entangled_input_output_probabilities': [str(p) for p in probabilities],
            'scalar_clock_shift_cancels_in_Heisenberg_SELECT': True,
            'one_sided_evolution_countercontrol_rejected': True,
            'opposite_filter_phase_countercontrol_rejected': True,
            'same_marginal_joint_reference_trace_distance': str(joint_distance),
            'control_kind': 'exact logical LCU matrix control; not native pulse execution'}
