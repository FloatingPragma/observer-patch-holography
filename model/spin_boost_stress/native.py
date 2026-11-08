"""Measured weighted-energy flow on the selected source ground."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.native_joint_spin_spectroscopy.circuits import controlled_pauli
from model.native_ordered_ground import experiment, source
from model.native_spectral_detector.circuits import hadamard_template
from model.native_spin_observer import algebra as spin, circuits
from model.certified_reference import certificate as cert
from model.spin_modular_geometry import native as engine, regional
from . import algebra


def terms(axis):
    op = algebra.operators(axis)
    coordinates = op['twice_coordinates']
    edges, odd = source.cube()
    coefficients = defaultdict(F)
    for x, y in edges:
        weight = F(coordinates[x]+coordinates[y], 16)
        coefficients['I'*8] -= weight
        for axis_name in 'XYZ':
            label = ['I']*8
            label[x] = label[y] = axis_name
            coefficients[''.join(label)] += weight
    for x in range(8):
        label = ['I']*8
        label[x] = 'X'
        coefficients[''.join(label)] += F((1 if x in odd else -1)*coordinates[x], 16)
    coefficients = {k: v for k, v in sorted(coefficients.items()) if v}
    matrix = sum((float(value)*spin.tensor(label) for label, value in coefficients.items()),
                 np.zeros((256, 256), complex))
    if not np.array_equal(matrix, op['B16']/16):
        raise ArithmeticError('native weighted-energy generator differs from the source density')
    return coefficients


@lru_cache(maxsize=1)
def run():
    axis, rapidity, subdivisions = 0, F(1, 1024), 4
    coefficients = terms(axis)
    source_row = algebra.source_controls()[axis]
    data = regional.state()
    full, _, _, _, preparation = experiment.plan()
    old_probe = circuits.probe('X', 8, 0, 16)
    prefix = full[:-len(old_probe)]
    end, _, variables = cert.consume(full, len(prefix), old_probe)
    if end != len(full) or variables:
        raise ArithmeticError('selected source preparation boundary changed')
    ordered = list(coefficients.items())+list(reversed(list(coefficients.items())))
    piece = []
    for label, value in ordered:
        physical = ''.join(axis_name+('X' if axis_name in 'XY' else 'I') for axis_name in label)+'I'
        piece += controlled_pauli(physical, -rapidity*value/(2*subdivisions), 16)
    template = hadamard_template(16)+piece*subdivisions+hadamard_template(16)
    psi = np.asarray(data['vector'], complex)/math.sqrt(data['vector_norm_squared'])
    evolved = psi.copy()
    for _ in range(subdivisions):
        for label, value in ordered:
            angle = -float(rapidity*value/(2*subdivisions))
            evolved = math.cos(angle)*evolved+1j*math.sin(angle)*(spin.tensor(label)@evolved)
    indices = np.array([sum(1 << (2*j+(0 if word >> j & 1 else 1)) for j in range(8))
                        for word in range(256)])
    target = np.zeros((1 << 17, 1), complex)
    target[indices, 0] = (psi+evolved)/2
    target[indices+(1 << 16), 0] = (psi-evolved)/2
    blank = np.zeros_like(target)
    blank[0, 0] = 1
    actual, report = engine.complete(template, 17, blank, target,
                                     prefix_word=prefix,
                                     prefix_error=F(preparation['complete_ground_isometry_error']))
    norm_sum = sum(abs(value) for value in coefficients.values())
    splitting = norm_sum**3*rapidity**3/(3*subdivisions**2)
    total_error = (F(report['complete_isometry_error_bound'])+data['ground_vector_error']+splitting)
    ideal_lower = (rapidity**2*source_row['ground_boost_second_moment_lower']/4
                   -rapidity**4*source_row['boost_norm_upper']**4/48)
    lower = ideal_lower-total_error
    probability = float(np.linalg.norm(actual.success[1 << 16:])**2)
    if lower <= F(5, 10**9) or probability < float(lower)-2e-11:
        raise ArithmeticError('source weighted-energy flow lost its measured ground response')
    report.update({'axis': axis, 'rapidity_parameter': str(rapidity),
                   'subdivisions': subdivisions, 'source_generator_terms': {k: str(v) for k, v in coefficients.items()},
                   'symmetric_product_error_upper': str(splitting),
                   'complete_true_source_output_TV_error_upper': str(total_error),
                   'ideal_minus_pointer_probability_lower': str(ideal_lower),
                   'native_unconditioned_minus_pointer_probability_lower': str(lower),
                   'native_unconditioned_minus_pointer_probability': probability,
                   'native_plus_pointer_probability': float(np.linalg.norm(actual.success[:1 << 16])**2),
                   'reserved_six_level_pointer_probabilities': [0.0]*4,
                   'input_is_blank': True, 'complete_work_rows': 1 << 17,
                   'physical_spacetime_boost_identified': False,
                   'selected_ground_is_invariant_under_this_weighted_energy_flow': False})
    return report
