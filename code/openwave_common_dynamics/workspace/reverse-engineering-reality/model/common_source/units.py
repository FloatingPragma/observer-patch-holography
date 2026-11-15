"""Conditional unit conversion, separate from any physical model calibration.

c fixes L/T and hbar fixes E*T. Even together they leave one dimensional
scale free until an independent clock, energy or length calibration is given.
No propagation speed is inferred from the present observer construction.
"""
from fractions import Fraction as F
import math


C_SI = F(299792458)
H_SI = F('6.62607015e-34')
ROWS = {'c': (1, -1, 0), 'hbar': (0, 1, 1), 'time_calibration': (0, 1, 0)}
SOURCES = {
    'SI': 'https://www.bipm.org/en/measurement-units/',
    'metre_definition': 'https://www.bipm.org/en/si-base-units/metre',
    'defining_constants': 'https://www.bipm.org/en/measurement-units/si-defining-constants',
}


def constraints(names):
    """Exact rank and logarithmic scale freedoms for the selected constraints."""
    if (not isinstance(names, (tuple, list)) or any(name not in ROWS for name in names)
            or len(set(names)) != len(names)):
        raise ValueError('distinct declared dimensional constraints required')
    matrix = [[F(value) for value in ROWS[name]] for name in names]
    pivots = []
    for col in range(3):
        pivot = next((i for i in range(len(pivots), len(matrix)) if matrix[i][col]), None)
        if pivot is None:
            continue
        row = len(pivots)
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        scale = matrix[row][col]
        matrix[row] = [value/scale for value in matrix[row]]
        for i in range(len(matrix)):
            if i != row:
                factor = matrix[i][col]
                matrix[i] = [a-factor*b for a, b in zip(matrix[i], matrix[row], strict=True)]
        pivots.append(col)
    basis = []
    for free in (col for col in range(3) if col not in pivots):
        vector = [F(0)]*3
        vector[free] = F(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -matrix[row][free]
        if next(value for value in vector if value) < 0:
            vector = [-value for value in vector]
        basis.append(vector)
    return {'variables': ['log length scale', 'log time scale', 'log energy scale'],
            'constraints': list(names), 'exponent_rows': [list(ROWS[name]) for name in names],
            'rank_exact': len(pivots), 'independent_scale_freedoms': 3-len(pivots),
            'nullspace_exponent_basis_exact': [[str(value) for value in row] for row in basis],
            'nullspace_meaning': 'each vector v permits multiplying scales by (lambda^vL,lambda^vT,lambda^vE)'}


def conditional_si_example(time_seconds):
    """Externally chosen seconds per model-time unit, imposing c_model=hbar_model=1."""
    if type(time_seconds) not in (int, F) or time_seconds <= 0:
        raise ValueError('a positive exact rational external time calibration required')
    time = F(time_seconds)
    length = C_SI*time
    energy_coefficient = H_SI/(2*time)
    return {'external_time_scale_seconds_exact': str(time),
            'external_time_scale_seconds_approx': float(time),
            'conditional_c_model': '1', 'conditional_hbar_model': '1',
            'length_scale_metres_exact': str(length),
            'length_scale_metres_approx': float(length),
            'energy_scale_joules_coefficient_of_inverse_pi_exact': str(energy_coefficient),
            'energy_scale_joules_approx': float(energy_coefficient)/math.pi,
            'c_reconstructed_metres_per_second_exact': str(length/time),
            'hbar_reconstructed_joule_seconds_coefficient_of_inverse_pi_exact': str(energy_coefficient*time),
            'external_time_scale_selected_by_model': False}


def report():
    examples = [conditional_si_example(F(1, 10**power)) for power in (15, 12)]
    return {'schema': 'oph.common_source.conditional_units.v1',
            'sources': SOURCES,
            'SI_defining_constants': {'c_metres_per_second_exact': str(C_SI),
                                      'h_joule_seconds_exact': str(H_SI),
                                      'hbar_definition': 'h/(2*pi)',
                                      'hbar_joule_seconds_coefficient_of_inverse_pi_exact': str(H_SI/2)},
            'dimensional_constraints': [constraints(['c']), constraints(['c', 'hbar']),
                                        constraints(['c', 'hbar', 'time_calibration'])],
            'conditional_examples': examples,
            'example_scale_ratios_second_over_first': {'time': '1000', 'length': '1000', 'energy': '1/1000'},
            'scope': {'emitted_model_propagation_speed_derived': False,
                      'spatial_readout_linked_to_modular_clock': False,
                      'external_calibration_selected': False,
                      'dimensionless_physical_predictions_established_by_unit_choice': False,
                      'meaning': 'unit constraints and conditional conversions only; SI constants are externally supplied'},
            'all_pass': True}
