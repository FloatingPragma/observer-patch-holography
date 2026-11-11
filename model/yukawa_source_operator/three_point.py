"""Complete finite source traces and normalized connected three-point records.

The source contains binary occupation observations and a faithful reference.
The declared insertions below are records of the same protected-core charge.
Their static connected amplitudes have no proved CAR/Higgs interpretation.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

import sympy as s

from ..common_source import slow_geometry
from ..native_modular_reference import local_preparation
from ..observer import native


def clean(value):
    return s.ImmutableMatrix(value).applyfunc(s.simplify)


def entries(value):
    return [[str(x) for x in row] for row in value.tolist()]


@lru_cache(maxsize=1)
def source_data():
    states = tuple(native.states(binary=True))
    probabilities = tuple(s.Rational(x) for x in local_preparation.reference())
    basis = s.Matrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a, b in row]
                      for row in slow_geometry.raw_frame()])
    if (len(states) != 66 or len(probabilities) != 66 or sum(probabilities) != 1
            or any(p <= 0 for p in probabilities)):
        raise ArithmeticError('the complete faithful 66-state source reference is required')
    gram = clean(basis.T*basis)
    if gram != (10+2*s.sqrt(5))*s.eye(3):
        raise ArithmeticError('source family frame Gram changed')
    values = tuple(s.ImmutableMatrix(basis.T*s.Matrix(x)) for x in states)
    core = tuple(s.Integer(sum(x[i] for i in local_preparation.CORE)) for x in states)
    return states, probabilities, basis, values, core


def expectation(values):
    probabilities = source_data()[1]
    if len(values) != len(probabilities):
        raise ValueError('one value for every source state required')
    initial = s.zeros(*values[0].shape) if isinstance(values[0], s.MatrixBase) else s.Integer(0)
    return sum((p*v for p, v in zip(probabilities, values, strict=True)), initial)


@lru_cache(maxsize=2)
def source_moments(core_power=1):
    """Exact first, second and connected third traces over every source state."""
    if type(core_power) is not int or core_power not in (1, 2):
        raise ValueError('the stated core-charge or squared-core-charge insertion is required')
    _, probabilities, _, values, core = source_data()
    h = tuple(x**core_power for x in core)
    mean, hmean = clean(expectation(values)), s.simplify(expectation(h))
    centered = tuple(clean(x-mean) for x in values)
    covariance = clean(expectation(tuple(x*x.T for x in centered)))
    variance = s.simplify(expectation(tuple((x-hmean)**2 for x in h)))
    mixed = clean(expectation(tuple(x*(y-hmean) for x, y in zip(centered, h, strict=True))))
    third = clean(expectation(tuple(x*x.T*(y-hmean) for x, y in zip(centered, h, strict=True))))
    if (not covariance.is_diagonal() or any(x.is_positive is not True for x in covariance.diagonal())
            or variance.is_positive is not True):
        raise ArithmeticError('positive diagonal source two-point form required')
    normalized = clean(s.Matrix(3, 3, lambda i, j:
        third[i, j]/s.sqrt(covariance[i, i]*covariance[j, j]*variance)))
    # Total expectation through actual core-charge records is an independent
    # aggregation of all66 labels; every charge stratum is retained.
    strata = defaultdict(list)
    for index, charge in enumerate(core):
        strata[charge].append(index)
    conditional = []
    reconstructed = s.zeros(3)
    for charge, indices in sorted(strata.items()):
        weight = sum(probabilities[i] for i in indices)
        conditional_second = clean(sum((probabilities[i]*centered[i]*centered[i].T
                                       for i in indices), s.zeros(3))/weight)
        reconstructed += weight*(charge**core_power-hmean)*conditional_second
        conditional.append({'core_charge': str(charge), 'state_count': len(indices),
                            'probability': str(weight),
                            'conditional_centered_family_second_moment': entries(conditional_second)})
    if clean(reconstructed-third) != s.zeros(3):
        raise ArithmeticError('complete conditional-expectation reconstruction failed')
    return {'mean': mean, 'H_mean': hmean, 'covariance': covariance,
            'H_variance': variance, 'H_family_covariance': mixed,
            'connected_third': third, 'normalized_third': normalized,
            'conditional': conditional}


def certificate():
    rows = []
    for power in (1, 2):
        data = source_moments(power)
        rows.append({'insertion': f'N_core**{power}',
                     'insertion_mean_exact': str(data['H_mean']),
                     'insertion_variance_exact': str(data['H_variance']),
                     'family_mean_exact': entries(data['mean']),
                     'family_covariance_exact': entries(data['covariance']),
                     'insertion_family_covariance_exact': entries(data['H_family_covariance']),
                     'connected_third_exact': entries(data['connected_third']),
                     'normalized_connected_third_exact': entries(data['normalized_third']),
                     'normalized_connected_eigenvalues_diagnostic':
                         [float(x) for x in data['normalized_third'].diagonal()],
                     'complete_core_conditional_records': data['conditional']})
    if source_moments(1)['normalized_third'] == source_moments(2)['normalized_third']:
        raise ArithmeticError('the distinct source-insertion witness collapsed')
    return {
        'source_states': 66, 'all_states_retained': True,
        'source_reference_normalization_exact': '1',
        'protected_core': list(local_preparation.CORE),
        'family_observable': 'q=B.T*n in the actual native rank-three coordinate frame',
        'three_point_definition': 'T_ab=E[(q_a-E[q_a])*(q_b-E[q_b])*(H-E[H])]',
        'normalized_definition': 'A_ab=T_ab/sqrt(C_aa*C_bb*Var(H)); C is exactly diagonal',
        'response_identity': 'd/dh Cov_(rho_h)(q_a,q_b)|h=0=T_ab, rho_h=rho*exp(hH)/E_rho[exp(hH)]',
        'positive_rescaling_of_each_observable_cancels_in_normalized_amplitude': True,
        'declared_source_insertions': rows,
        'same_reference_and_same_modular_flow_at_zero_insertion': True,
        'different_insertions_give_different_normalized_response_exact': True,
        'source_intervention_is_selected_by_zero_field_reference_alone': False,
        'finite_native_preparation_errors_propagated': False,
        'static_two_point_covariance_is_proved_physical_kinetic_residue': False,
        'source_scalar_insertions_are_proved_Higgs_doublet': False,
        'occupation_observables_are_proved_charged_CAR_fields': False,
        'physical_Yukawa_coefficients_derived': False,
    }
