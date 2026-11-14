"""Complete finite-environment history budgets using the existing theorems."""
from fractions import Fraction as F
import math

from . import continuation_bounds as continuation
from . import open_environment as noise, protocol


def budget(program, environment, *, prepared_noiseless_error=None):
    """Bounds for distinct, explicitly named complete output comparisons.

    Every phase commutes with the number moment and cutoff projection, so
    inserting the matched phase channel adds no occupation-bound increment.
    Temporal and Gaussian-quadrature errors telescope through arbitrary
    complete later instruments, including reference-entangled inputs.
    """
    program = protocol.validate_program(program)
    environment = noise.validated(program, environment)
    inherited = continuation.program_bound(program)
    channels = noise.schedule_budget(program, environment)
    temporal = F(inherited['ideal_temporal_instrument_error_upper'])
    quadrature = F(channels['ideal_finite_to_gaussian_half_diamond_upper'])
    cutoff = F(inherited['occupation_certificate']['complete_instrument_trace_distance_upper'])
    effect = F(channels['ideal_finite_to_identity_half_diamond_upper'])
    finite = min(F(1), temporal+quadrature)
    infinite = min(F(1), cutoff+temporal+quadrature)
    comparison = None
    if prepared_noiseless_error is not None:
        if type(prepared_noiseless_error) not in (int, F) or not 0 <= prepared_noiseless_error <= 1:
            raise ValueError('an exact separately proved prepared noiseless error is required')
        comparison = {
            'prepared_noiseless_ideal_product_to_infinite_error_upper': str(prepared_noiseless_error),
            'ideal_finite_noise_to_noiseless_finite_error_upper': str(effect),
            'ideal_noisy_finite_to_NOISELESS_infinite_error_upper': str(
                min(F(1), F(prepared_noiseless_error)+effect)),
            'target_contains_phase_noise': False,
            'prepared_input_contract': 'independent blank gauge and declared scalar vacuum',
            'supplied_prepared_certificate_recomputed_by_this_function': False,
        }
    required = (*continuation.required_service_keys(program.operations),
                *(f'environment_{index}:site_{event.site}@{event.after}'
                  for index, event in enumerate(environment.events)))
    return {
        'schema': 'oph.common_source.open_history_budget.v1',
        'canonical_program': program.specification(), 'environment': environment.specification(),
        'environment_budget': channels,
        'ideal_temporal_instrument_error_upper': str(temporal),
        'ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper': str(finite),
        'finite_comparison_includes_arbitrary_work_reference_inputs': True,
        'matched_phase_noise_occupation_certificate': inherited['occupation_certificate'],
        'ideal_finite_source_to_infinite_GAUSSIAN_history_error_upper': str(infinite),
        'occupation_comparison_requires_declared_preparation': True,
        'number_moment_growth_and_cutoff_leakage_from_noise_exact': '0',
        'prepared_noiseless_reference_comparison': comparison,
        'physical_service_coverage': continuation.service_error_budget(finite, required_stage_keys=required),
        'floating_execution_error_certified': False,
        'simultaneous_Hamiltonian_Lindblad_flow_claimed': False,
        'scope': 'complete quantum/reference/record output; independent phase channels at declared boundaries',
    }


def gaussian_refinement(dimension, variances):
    """Analytic quadrature orders for the Gaussian occupation limit.

    M!/(2M)! <= (M+1)^(-M). With M >= 2*v*(d-1)^2 and M >= d,
    each local half-diamond error is <= d*2^(-M-1). The finite executor
    currently supports only d=2,4 and M<=64; this does not expand that domain.
    """
    if type(dimension) is not int or dimension < 2:
        raise ValueError('integer occupation dimension at least two required')
    if type(variances) is not tuple or any(type(v) not in (int, F) or v < 0 for v in variances):
        raise ValueError('tuple of exact nonnegative phase variances required')
    orders, errors = [], []
    # Round outward on a dimension-dependent grid. The rounding allowance
    # tends to zero and avoids constructing/serializing enormous powers of 2.
    scale = 10**12*(dimension+2)**6
    for variance in variances:
        variance = F(variance)
        threshold = 2*variance*(dimension-1)**2
        order = max(dimension, math.ceil(threshold))
        orders.append(order)
        if variance == 0:
            error = F(0)
        elif order+1 > dimension.bit_length()+scale.bit_length():
            error = F(1, scale)
        else:
            exact = F(dimension, 2**(order+1))
            error = F(math.ceil(exact*scale), scale)
        errors.append(error)
    return {'dimension': dimension, 'orders': orders,
            'ideal_exact_parameter_quadrature_error_upper': str(min(F(1), sum(errors, F(0)))),
            'outward_rounding_grid': str(F(1, scale)),
            'exact_unrounded_formula_per_nonzero_variance': 'd*2^(-M-1)',
            'source_execution_available': dimension in (2, 4) and all(m <= 64 for m in orders),
            'parameter_rounding_control_required': True,
            'scope': 'analytic ideal family, not an execution beyond the implemented cutoffs'}
