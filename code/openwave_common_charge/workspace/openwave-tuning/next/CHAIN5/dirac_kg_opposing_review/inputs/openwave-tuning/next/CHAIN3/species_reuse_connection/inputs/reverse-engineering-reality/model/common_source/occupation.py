"""Rational occupation bounds for the actual interacting common history.

The scalar preparation is the existing projected free coupled vacuum; the
subsequent state need not be Gaussian. A weighted number-moment inequality
allows arbitrary matter/reference inputs and the retained intermediate
matter measurement/feedback. The final discontinuous scalar-sign readout
is deliberately a separate comparison: this bound ends immediately before
that readout and cannot be passed through two different cutoff sign maps.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

from ..observer.source_scalar.bounds import root_upper
from ..observer.source_scalar.cutoff import exponential_upper


KAPPA = F(9, 32)
MAXIMUM_COUPLING = F(1, 4)
MAXIMUM_TIME = F(3, 2)
INITIAL_WEIGHT = F(8)
YOUNG_PARAMETER = F(1, 2)
CHARACTERISTIC_RATE = 2 * KAPPA + MAXIMUM_COUPLING * YOUNG_PARAMETER
MOMENT_RATE = KAPPA + MAXIMUM_COUPLING / YOUNG_PARAMETER


def _time(value):
    if type(value) is not F or not 0 <= value <= MAXIMUM_TIME:
        raise ValueError("an exact Fraction model time between zero and 3/2 required")
    return value


def _dimension(value):
    if type(value) is not int or value < 2 or value % 2:
        raise ValueError("an even integer onsite occupation dimension at least two required")
    return value


def _root_lower(value, scale):
    return F(math.isqrt(value.numerator * scale * scale // value.denominator), scale)


def preparation_tail(d):
    """Omitted norm squared of the exact preparation's total-N projection.

    The actual vector keeps N=0,2,...,d-2 and is then normalized. Only this
    initial state is Gaussian; no Gaussian propagation hypothesis is used.
    Precision grows with d, so directed rounding creates no positive floor.
    """
    d = _dimension(d)
    term, retained_series = F(1), F(0)
    for pair in range(d // 2):
        retained_series += term
        term *= F(2 * pair + 1, 81 * (2 * pair + 2))
    return 1 - F(4, 9) * _root_lower(F(5), 10 ** (2 * d + 30)) * retained_series


def number_moment_bound(time=MAXIMUM_TIME):
    """Uniform <q**N> bound through the interacting retained first read.

    For M(s,t)=<exp(sN)>, weighted Cauchy and Young inequalities give
    M_t <= a*sinh(s)*M_s + b*sinh(s)*M, with a=11/16, b=25/32.
    Along tanh(s/2)=exp(-a*t)*tanh(log(8)/2), the bound is
    sqrt(80/17)*(8/q)**(25/22). We choose a rational q below this
    characteristic and bound its fractional power by concavity. These
    constants cover both g=1/4 and the disconnected g=0 control.
    """
    time = _time(time)
    growth = exponential_upper(CHARACTERISTIC_RATE * time, digits=40)
    initial_tanh = (INITIAL_WEIGHT - 1) / (INITIAL_WEIGHT + 1)
    final_tanh_lower = initial_tanh / growth
    characteristic_weight_lower = (1 + final_tanh_lower) / (1 - final_tanh_lower)
    # A modest denominator keeps the exported q**d certificate manageable.
    q = F((64 * characteristic_weight_lower).numerator
          // (64 * characteristic_weight_lower).denominator, 64)
    if not 1 < q <= INITIAL_WEIGHT:
        raise ArithmeticError("the fixed time window must retain an exponential number weight")
    # q <= characteristic q follows from this exact rational inequality.
    if (q - 1) / (q + 1) * growth > initial_tanh:
        raise ArithmeticError("the rational weight does not lie below its characteristic")
    initial_moment = root_upper(F(80, 17), scale=10**40)
    ratio = INITIAL_WEIGHT / q
    # 25/22 = 1+3/22. Use tangents to x**(3/22) at 1 and 4.
    # The second tangent has an irrational factor 4**(3/22), bounded
    # by 121/100 because (121/100)**22 >= 4**3. Every operation here
    # is rational; taking the smaller of two upper bounds is valid.
    fourth_power_upper = F(121, 100)
    if fourth_power_upper**22 < 4**3:
        raise ArithmeticError("the rational tangent coefficient must be an upper bound")
    fractional_power_upper = min(
        1 + F(3, 22) * (ratio - 1),
        fourth_power_upper * (1 + F(3, 88) * (ratio - 4)),
    )
    power_upper = ratio * fractional_power_upper
    moment = initial_moment * power_upper
    return {
        "total_model_time": str(time),
        "number_convention": "N=N0+N1 in the original unit-frequency onsite Fock basis",
        "maximum_absolute_coupling": str(MAXIMUM_COUPLING),
        "initial_exponential_weight": str(INITIAL_WEIGHT),
        "exponential_weight": str(q),
        "characteristic_rate": str(CHARACTERISTIC_RATE),
        "moment_rate": str(MOMENT_RATE),
        "characteristic_exponential_upper": str(growth),
        "initial_number_moment_upper": str(initial_moment),
        "characteristic_power_upper": str(power_upper),
        "number_moment_upper": str(moment),
        "interacting_state_assumed_gaussian": False,
        "complete_intermediate_matter_records_retained": True,
        "arbitrary_matter_and_untouched_reference": True,
    }


def occupation_cutoff_bound(d, time=MAXIMUM_TIME):
    """Complete ideal history immediately before final scalar-sign readout.

    Compare the infinite scalar vacuum to the normalized total-N<=d-2
    preparation already used by field.vacuum_vector, and compare H with
    P_d H P_d for P_d={n0,n1<d}. Matter preparation, the first complete
    matter measurement, and its controlled matter feedback commute with
    this projection and preserve the number bound. Total evolution time
    is at most 3/2, split into the declared intervals 1 and 1/2.

    The bound is uniform over the initial matter plus any untouched
    reference. It is NOT a diamond bound over arbitrary scalar inputs.
    Finite splitting, numerical execution and physical errors are separate.
    """
    d, time = _dimension(d), _time(time)
    moment_certificate = number_moment_bound(time)
    q = F(moment_certificate["exponential_weight"])
    moment = F(moment_certificate["number_moment_upper"])
    tail = moment / q**d
    scale = 10 ** (d + 30)
    tail_norm = root_upper(tail, scale=scale)
    initial_tail = preparation_tail(d)
    initial_error = root_upper(2 * initial_tail, scale=scale)
    # P H (1-P) receives two onsite lowering-pair terms, three bond
    # monomials and the one force-lowering term. Matter terms give zero.
    cross_norm = (KAPPA * root_upper(F(d * (d + 1)), scale=scale) / 2
                  + 3 * KAPPA * d / 2
                  + MAXIMUM_COUPLING * root_upper(F(d), scale=scale))
    duhamel = time * cross_norm * tail_norm
    projected_error = initial_error + duhamel
    vector_error = root_upper(projected_error**2 + tail, scale=scale)
    trace_error = min(F(1), vector_error)
    return {
        "dimension": d,
        "total_model_time": str(time),
        "number_moment_certificate": moment_certificate,
        "preparation_omitted_probability_upper": str(initial_tail),
        "preparation_vector_error_upper": str(initial_error),
        "rectangular_outside_probability_upper": str(min(F(1), tail)),
        "unclipped_rectangular_outside_probability_upper": str(tail),
        "projected_hamiltonian_cross_norm_upper": str(cross_norm),
        "projected_duhamel_error_upper": str(duhamel),
        "projected_history_vector_error_upper": str(projected_error),
        "complete_history_vector_error_upper": str(vector_error),
        "pre_final_readout_trace_distance_upper": str(trace_error),
        "common_postprocessing_probability_error_upper": str(trace_error),
        "first_matter_record_and_feedback_included": True,
        "all_unnormalized_first_record_blocks_included": True,
        "arbitrary_matter_and_untouched_reference": True,
        "arbitrary_scalar_input_diamond_bound": False,
        "final_scalar_sign_readout_comparison_included": False,
        "finite_splitting_error_included": False,
        "numerical_execution_or_physical_error_included": False,
        "rate": "O(d*q**(-d/2)) for fixed time, with q>1; initial preparation converges faster",
        "scope": "ideal interacting history before the final scalar sign readout, retaining the first matter record and quantum state; fixed scalar vacuum, arbitrary matter/reference input",
    }


def certificate_table(dimensions=(32, 64, 128)):
    """Concrete rational comparisons; these dimensions are not source runs."""
    if type(dimensions) not in (tuple, list) or not dimensions:
        raise ValueError("a nonempty list or tuple of even onsite dimensions required")
    return {str(_dimension(d)): occupation_cutoff_bound(d) for d in dimensions}
