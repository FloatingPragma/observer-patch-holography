"""Complete finite Pauli-score instrument for an actual bond coefficient."""
from fractions import Fraction as F
from math import ceil


def plan(coefficients, excluded_axis, score_accuracy=F(1,100), failure_probability=F(1,100)):
    if excluded_axis not in (1,2,3):
        raise ValueError('one physical spin axis is required')
    coefficients = tuple(F(x) for x in coefficients)
    total = sum(map(abs, coefficients), F())
    accuracy, failure = F(score_accuracy), F(failure_probability)
    if accuracy <= 0 or not 0 < failure < 1:
        raise ValueError('valid accuracy budget required')
    if total == 0:
        return {'complete_fine_records':[{'record':'identity','score':'0'}],
                'score_absolute_bound':'0',
                'independent_repetitions_for_Chebyshev_bound':0,
                'score_accuracy':str(accuracy),'failure_probability_upper':str(failure),
                'observable':'zero','score_is_spectral_measurement_of_entire_sum':True,
                'all_setting_and_local_Pauli_records_retained':True,
                'native_coherent_control_and_local_Pauli_pointer_services_required':False,
                'large_native_tape_executed':False}
    rows = []
    for edge, value in enumerate(coefficients):
        if not value:
            continue
        for axis in (1,2,3):
            if axis == excluded_axis:
                continue
            probability = abs(value)/(2*total)
            for left in (-1,1):
                for right in (-1,1):
                    rows.append({'edge_index':edge, 'Pauli_axis':axis,
                        'setting_probability':str(probability),
                        'left_record':left, 'right_record':right,
                        'score':str(total*F(1 if value > 0 else -1,2)*left*right)})
    scale = total/2
    shots = ceil(scale**2/(failure*accuracy**2))
    return {'complete_fine_records':rows, 'score_absolute_bound':str(scale),
            'independent_repetitions_for_Chebyshev_bound':shots,
            'score_accuracy':str(accuracy), 'failure_probability_upper':str(failure),
            'observable':'sum_b z_b sum_(a != excluded_axis) S_x^a S_y^a',
            'score_is_spectral_measurement_of_entire_sum':False,
            'all_setting_and_local_Pauli_records_retained':True,
            'native_coherent_control_and_local_Pauli_pointer_services_required':True,
            'large_native_tape_executed':False}
