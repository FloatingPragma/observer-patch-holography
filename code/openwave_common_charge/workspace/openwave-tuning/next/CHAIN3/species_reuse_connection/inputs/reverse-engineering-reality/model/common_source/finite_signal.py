"""Exact finite detector-response enclosure for the common interacting model.

The d=2 scalar projection makes the complete 16-dimensional Hamiltonian
rational. Taylor actions and retained branch probabilities are computed in
Gaussian rationals. An operator-norm integral remainder encloses the ideal
continuous-model-time instrument; no numerical execution supplies a bound.
The optional inherited splitting budget also encloses its ideal source word.
Spatial/occupation limits, numerical controls and physical services are not
certified here.
"""
from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
import math

from . import gauge_matter as gauge


Vector = tuple[gauge.Gaussian, ...]
DIMENSION = 16
HAMILTONIAN_NORM_UPPER = Q(109, 16)
DEFAULT_DEGREE = 64
INTERVAL_SCALE = 10**15


def _degree(value):
    if type(value) is not int or not 0 <= value <= 256:
        raise ValueError("Taylor degree must be an integer from 0 through 256")
    return value


def _time(value):
    if type(value) not in (int, Q) or value < 0:
        raise ValueError("model time must be a nonnegative exact rational")
    return Q(value)


def _vector(value):
    if (not isinstance(value, tuple) or len(value) != DIMENSION
            or any(not isinstance(cell, tuple) or len(cell) != 2
                   or any(type(part) is not Q for part in cell) for cell in value)):
        raise ValueError("a 16-component immutable Gaussian-rational vector is required")
    return value


def _sum(a, b):
    return a[0]+b[0], a[1]+b[1]


def _scale_real(a, scale):
    return a[0]*scale, a[1]*scale


def initial_vector() -> Vector:
    """The normalized d=2 projected scalar preparation and blank gauge code."""
    return (gauge.ONE,)+(gauge.ZERO,)*(DIMENSION-1)


def apply_hamiltonian(vector: Vector) -> Vector:
    """Independent exact H on b0,b1,n0,n1, in little-endian order.

    H = X0+X1+E² +(73/64)(n0+n1+1)
        -(9/64)X2 X3 +(1/4)X2 Z0,
    where E=-b0+b1. The scalar squares use the actual projected infinite
    oscillator squares, including their finite-boundary correction.
    """
    vector = _vector(vector)
    output = []
    for index in range(DIMENSION):
        b0, b1, n0, n1 = ((index >> bit) & 1 for bit in range(4))
        diagonal = Q((b1-b0)**2)+Q(73, 64)*(n0+n1+1)
        result = _scale_real(vector[index], diagonal)
        result = _sum(result, vector[index ^ 1])
        result = _sum(result, vector[index ^ 2])
        result = _sum(result, _scale_real(vector[index ^ 12], -Q(9, 64)))
        result = _sum(result, _scale_real(vector[index ^ 4], Q(1-2*b0, 4)))
        output.append(result)
    return tuple(output)


def hamiltonian_matrix() -> gauge.Matrix:
    columns = [apply_hamiltonian(tuple(gauge.ONE if i == j else gauge.ZERO
                                       for i in range(DIMENSION)))
               for j in range(DIMENSION)]
    return tuple(tuple(column[row] for column in columns) for row in range(DIMENSION))


def taylor_action(vector: Vector, time: Q, degree=DEFAULT_DEGREE):
    """Return polynomial action and an operator-norm remainder upper bound.

    The remainder is epsilon times the norm of the supplied vector. It is
    not an absolute vector bound for an arbitrary unnormalized input.
    """
    vector, time, degree = _vector(vector), _time(time), _degree(degree)
    total, term = vector, vector
    for power in range(1, degree+1):
        applied = apply_hamiltonian(term)
        factor = time/power
        term = tuple((imag*factor, -real*factor) for real, imag in applied)
        total = tuple(_sum(a, b) for a, b in zip(total, term, strict=True))
    epsilon = (HAMILTONIAN_NORM_UPPER*time)**(degree+1)/math.factorial(degree+1)
    return total, epsilon


def first_projection(vector: Vector, outcome: int) -> Vector:
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError("first outcome must be 0 or 1")
    return tuple(cell if (index & 1) == outcome else gauge.ZERO
                 for index, cell in enumerate(_vector(vector)))


def controlled_feedback(vector: Vector, outcome: int, enabled: bool) -> Vector:
    vector = _vector(vector)
    if type(outcome) is not int or outcome not in (0, 1) or type(enabled) is not bool:
        raise ValueError("a binary committed record and Boolean feedback switch required")
    if not enabled or outcome == 0:
        return vector
    # (3I-4iX0)/5, the same source operation as the joint-history runner.
    return tuple((Q(3, 5)*real+Q(4, 5)*vector[index ^ 1][1],
                  Q(3, 5)*imag-Q(4, 5)*vector[index ^ 1][0])
                 for index, (real, imag) in enumerate(vector))


def scalar_sign_projection(vector: Vector, outcome: int) -> Vector:
    """B1=(I+X2)/2 reads positive Q0; B0=(I-X2)/2 reads negative Q0."""
    vector = _vector(vector)
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError("scalar sign outcome must be 0 or 1")
    sign = 2*outcome-1
    return tuple(_scale_real(_sum(cell, _scale_real(vector[index ^ 4], sign)), Q(1, 2))
                 for index, cell in enumerate(vector))


def norm_squared(vector: Vector) -> Q:
    return sum((real*real+imag*imag for real, imag in _vector(vector)), Q(0))


def branch_vectors(feedback_enabled: bool, degree=DEFAULT_DEGREE):
    """All four UNNORMALIZED approximate branches and their shared error.

    The approximation does not renormalize, postselect, combine coherent
    record values, or use a supplied detector probability as an input.
    """
    if type(feedback_enabled) is not bool:
        raise ValueError("feedback switch must be Boolean")
    degree = _degree(degree)
    return _branch_vectors(feedback_enabled, degree)


@lru_cache(maxsize=12)
def _branch_vectors(feedback_enabled: bool, degree: int):
    first, epsilon1 = taylor_action(initial_vector(), Q(1), degree)
    branches = []
    epsilon2 = Q(0)
    for a in (0, 1):
        selected = first_projection(first, a)
        controlled = controlled_feedback(selected, a, feedback_enabled)
        second, epsilon2 = taylor_action(controlled, Q(1, 2), degree)
        for b in (0, 1):
            branches.append(((a, b), scalar_sign_projection(second, b)))
    eta = epsilon2+(1+epsilon2)*epsilon1
    return tuple(branches), epsilon1, epsilon2, eta


def _outward(center: Q, radius: Q, limits):
    """Outward rational rounding only, to keep retained fractions readable."""
    low, high = max(Q(limits[0]), center-radius), min(Q(limits[1]), center+radius)
    lower = (low.numerator*INTERVAL_SCALE)//low.denominator
    upper = -((-high.numerator*INTERVAL_SCALE)//high.denominator)
    return [str(Q(lower, INTERVAL_SCALE)), str(Q(upper, INTERVAL_SCALE))]


def _source_word_errors(steps):
    # The inherited finite proof uses exact rational root enclosures.
    # It does not call any source unitary or floating target evolution.
    from .scalar_coupling import time_error_bound
    if type(steps) is not int or steps < 1:
        raise ValueError("a positive integer source step count is required")
    segments = [time_error_bound(2, steps, time, Q(1, 4))
                for time in (Q(1), Q(1, 2))]
    probability_error = sum((Q(row["ideal_channel_half_diamond_error_upper"])
                             for row in segments), Q(0))
    return segments, min(Q(1), probability_error)


def certificate(degree=DEFAULT_DEGREE, steps=512):
    """JSON-compatible rigorous ideal-probability and response intervals.

    Each model law is compared with and without the same record-controlled
    operation. The contrast is p(positive scalar | feedback enabled) minus
    p(positive scalar | feedback disabled), without conditioning on either
    first outcome. Both first records and both scalar outcomes are retained.
    """
    degree = _degree(degree)
    source_segments, source_probability_error = _source_word_errors(steps)
    configurations = []
    centers = []
    probability_error = Q(0)
    epsilon1 = epsilon2 = eta = Q(0)
    for enabled in (True, False):
        branches, epsilon1, epsilon2, eta = branch_vectors(enabled, degree)
        branch_error = 2*eta+eta*eta
        probability_error = 2*branch_error
        probabilities = [(record, norm_squared(vector)) for record, vector in branches]
        center = sum((value for record, value in probabilities if record[1] == 1), Q(0))
        centers.append(center)
        configurations.append({
            "feedback_enabled": enabled,
            "retained_branches": [
                {"record": list(record),
                 "continuous_probability_interval": _outward(value, branch_error, (0, 1))}
                for record, value in probabilities],
            "positive_scalar_continuous_probability_interval": _outward(center, probability_error, (0, 1)),
            "positive_scalar_ideal_source_word_probability_interval": _outward(
                center, probability_error+source_probability_error, (0, 1)),
        })
    contrast = centers[0]-centers[1]
    continuous_interval = _outward(contrast, 2*probability_error, (-1, 1))
    source_interval = _outward(contrast, 2*(probability_error+source_probability_error), (-1, 1))
    nonzero = lambda interval: Q(interval[1]) < 0 or Q(interval[0]) > 0
    return {
        "method": "exact Gaussian-rational Taylor actions with unitary integral remainder",
        "onsite_scalar_dimension": 2, "system_dimension": DIMENSION,
        "logical_bits": ["species0_right", "species1_right", "scalar0_occupation", "scalar1_occupation"],
        "Hamiltonian_exact": "X0+X1+(-b0+b1)^2+(73/64)*(n0+n1+1)-(9/64)*X2*X3+(1/4)*X2*Z0",
        "norm_term_bounds": ["2", "1", "219/64", "9/64", "1/4"],
        "Hamiltonian_norm_upper": str(HAMILTONIAN_NORM_UPPER),
        "initial_state": "|0000>, the d=2 normalized projected scalar preparation and gauge blank",
        "model_time_segments": ["1", "1/2"],
        "first_projector": "species0 right occupation (bit0)",
        "feedback_if_first_record_one": "(3I-4iX0)/5",
        "second_projectors": {"0": "(I-X2)/2", "1": "(I+X2)/2"},
        "Taylor_degree": degree,
        "Taylor_operator_error_upper": [str(epsilon1), str(epsilon2)],
        "each_branch_vector_error_upper": str(eta),
        "positive_scalar_probability_error_upper": str(probability_error),
        "configurations": configurations,
        "contrast_definition": "unconditional scalar-positive probability with feedback minus without feedback",
        "continuous_contrast_interval": continuous_interval,
        "continuous_nonzero_response_certified": nonzero(continuous_interval),
        "source_steps_per_segment": steps,
        "source_segment_certificates": source_segments,
        "each_configuration_source_probability_error_upper": str(source_probability_error),
        "ideal_source_word_contrast_interval": source_interval,
        "ideal_source_word_nonzero_response_certified": nonzero(source_interval),
        "scope": {
            "same_declared_finite_interacting_instrument": True,
            "both_first_records_and_scalar_outcomes_retained": True,
            "exact_rational_probability_enclosures": True,
            "floating_execution_error_included": False,
            "physical_service_error_included": False,
            "occupation_or_spatial_continuum_limit": False,
            "comparison_with_nature": False,
        },
    }
