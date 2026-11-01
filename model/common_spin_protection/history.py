"""Complete interacting evolution/repair histories on one spin factor.

This evaluates the declared ideal full-Fock operators.  It is not source
pulse execution or a proof that physical encoders and global decoders exist.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..charged_defect import branches, dictionary
from . import spin, spectrum


ADMISSION = 2e-10


def spin_reference_marginal(state, chart=None):
    """Trace multiplicity, keeping every reference column and spin weight.

    Input state has shape (4096, reference, earlier_record).  Its norm need
    not be one: complete outcomes are combined without success conditioning.
    """
    if np.ma.isMaskedArray(state):
        raise ValueError('masked state columns are not a complete input')
    value = np.asarray(state)
    if (value.ndim != 3 or value.shape[0] != 4096 or min(value.shape[1:]) < 1
            or not np.isfinite(value).all() or value.dtype.kind not in 'fci'):
        raise ValueError('finite full-Fock/reference/record columns required')
    chart = spectrum.chart(4, 2) if chart is None else chart
    flattened = value.reshape(4096, -1)
    components = np.stack([chart.complete_columns(j, original=True).conj().T@flattened
                           for j in range(chart.twice_spin+1)], axis=1)
    components = components.reshape(chart.dimension, chart.twice_spin+1, *value.shape[1:])
    # m = multiplicity, a/b = protected spin, r/s = reference, h = old record.
    density = np.einsum('marh,mbsh->arbs', components, components.conj()).reshape(
        (chart.twice_spin+1)*value.shape[1], -1)
    projected_mass = float(np.vdot(components, components).real)
    total_mass = float(np.vdot(value, value).real)
    if abs(projected_mass-total_mass) > ADMISSION:
        raise ArithmeticError('history left its declared complete spin/charge block')
    return density


def report():
    chart = spectrum.chart(4, 2)
    triplet = branches.COUPLING[:, 1:]/np.sqrt([1., 2., 1.])
    initial = dictionary.embedded_columns(3, 'hole')[:, :4]@triplet/math.sqrt(3)
    expected = spin_reference_marginal(initial[:, :, None], chart)
    generator = spin.operators().K.matrix(spin.operators().scale)
    first = expm_multiply(-1j*generator, initial)
    # Retain all E slots and its full-Fock outside outcome, including zeros.
    first_operators = spin.repair_outcomes('hole', 'E', (0, 2))
    first_outputs = np.stack([operator@first for _, operator in first_operators], axis=2)
    first_probabilities = np.sum(abs(first_outputs)**2, axis=(0, 1))
    continued = expm_multiply(-.25j*generator, first_outputs.reshape(4096, -1)).reshape(first_outputs.shape)
    final_marginal = np.zeros_like(expected)
    probabilities = []
    second_operators = spin.repair_outcomes('hole', 'C', (0, 1))
    for second_label, operator in second_operators:
        output = (operator@continued.reshape(4096, -1)).reshape(continued.shape)
        final_marginal += spin_reference_marginal(output, chart)
        masses = np.sum(abs(output)**2, axis=(0, 1))
        for (first_label, _), mass in zip(first_operators, masses, strict=True):
            probabilities.append({'first': list(first_label), 'second': list(second_label),
                                  'probability': float(mass)})
    error = float(np.linalg.norm(final_marginal-expected))
    total = sum(row['probability'] for row in probabilities)
    if error > ADMISSION or abs(total-1) > ADMISSION:
        raise ArithmeticError('complete joined history failed its protected reference marginal')
    if first_probabilities[-1] < 1e-5:
        raise ArithmeticError('the old small-carrier leakage control disappeared')
    # Merely retaining S^2 is insufficient.  A complete S_z measurement
    # preserves this spin-one sector but dephases its reference entanglement.
    dephased = np.zeros_like(expected)
    for weight in range(3):
        sl = slice(3*weight, 3*weight+3)
        dephased[sl, sl] = expected[sl, sl]
    control_distance = float(np.sum(abs(np.linalg.eigvalsh(expected-dephased)))/2)
    if abs(control_distance-2/3) > ADMISSION:
        raise ArithmeticError('spin-sector-only control failed to destroy the encoded qutrit')
    return {'initial': 'old localized hole-pair triplet maximally entangled with a three-state reference',
            'number': 4, 'total_spin': 1,
            'old_triplet_carrier_dimension': 9,
            'expanded_multiplicity_dimension': chart.dimension,
            'expanded_complete_spin_sector_dimension': chart.dimension*3,
            'field_durations_in_supplied_Hamiltonian_units': [1., .25],
            'first_law': 'E', 'first_seam': [0, 2],
            'second_law': 'C', 'second_seam': [0, 1],
            'first_complete_probabilities': first_probabilities.tolist(),
            'first_outside_small_carrier_probability': float(first_probabilities[-1]),
            'complete_joint_outcomes': probabilities,
            'total_probability': total,
            'final_spin_reference_Frobenius_error': error,
            'initial_spin_reference_purity': float(np.trace(expected@expected).real),
            'spin_z_record_control_reference_trace_distance': control_distance,
            'exact_spin_z_control_distance': '2/3',
            'numerical_admission': ADMISSION,
            'postselection_used': False, 'physical_pulse_source_executed': False}
