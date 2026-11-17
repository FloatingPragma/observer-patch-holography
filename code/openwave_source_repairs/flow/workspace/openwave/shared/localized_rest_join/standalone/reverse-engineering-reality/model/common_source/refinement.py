"""Occupation refinement of the declared interacting two-site history.

The fourth-number-moment/Duhamel certificate is quantitative up to the last
field readout. The actual finite sign(Q_d) instrument also converges strongly,
as proved in OCCUPATION_REFINEMENT.md, but no finite-d rate for that sharp
readout is supplied by this module. The separate sharp_output module gives
state-specific complete-output finite-cutoff enclosures. None is a spatial continuum or
physical-service certificate. The source execution remains at d=2 and d=4.
"""
from fractions import Fraction as F
import json
import math

from ..observer.source_scalar.bounds import root_upper
from ..observer.source_scalar.cutoff import exponential_upper


KAPPA = F(9, 32)
COUPLING = F(1, 4)
MOMENT_RATE = 8 * KAPPA + 4 * abs(COUPLING)


def _dimension(d):
    if type(d) is not int or d < 2 or d % 2:
        raise ValueError("an even onsite occupation dimension at least two is required")
    return d


def _time(time):
    if type(time) not in (int, F) or not 0 <= time <= F(3, 2):
        raise ValueError("exact rational total model time between zero and 3/2 required")
    return F(time)


def initial_fourth_moment():
    """Exact <(N+2)^4> in the actual free coupled-oscillator vacuum.

    N=2R, R has generating function sqrt(80/(81-z)), and its factorial
    moments are (1/2)_k / 80^k. Conditioning on N <= d-2 only decreases
    the expectation of this increasing weight.
    """
    factorial = [F(1)]
    for k in range(1, 5):
        factorial.append(factorial[-1] * F(2 * k - 1, 160))
    return 16 * sum(coefficient * moment for coefficient, moment in
                    zip((1, 15, 25, 10, 1), factorial, strict=True))


def preparation_tail_upper(d):
    """Rational omitted mass for the actual total-occupation preparation.

    The first omitted pair is R=d/2. Consecutive vacuum pair probabilities
    have ratio <1/81; dropping their normalization sqrt(80/81) increases
    this geometric tail bound. This only concerns the initial preparation.
    """
    r = _dimension(d) // 2
    return F(81, 80) * F(math.comb(2 * r, r), 4**r * 81**r)


def occupation_evolution_bound(d, time=F(3, 2)):
    """Certified finite/infinite error immediately BEFORE sign(Q0) readout.

    Both histories use the same declared interacting Hamiltonian, the same
    gauge-bit measurement and its conditional gauge rotation. The comparison
    retains every unnormalized branch and its quantum output; it never
    conditions on a rare outcome. A native-current gauge preparation kick
    commutes with N, so the bound is uniform over its three declared values.

    The finite evolution here is exp(-it Pi_d H Pi_d). Finite product-word
    error, floating contraction, final sign(Q_d) readout and physical-service
    error must not be silently included in this number.
    """
    d, t = _dimension(d), _time(time)
    # Polynomially increasing precision avoids a positive rounding floor as
    # d increases without creating exponentially long returned root bounds.
    scale = 10**30 * (d + 2)**6
    moment = initial_fourth_moment()
    amplitude = root_upper(moment, scale=scale)
    growth = exponential_upper(MOMENT_RATE * t / 2)
    final_amplitude = amplitude * growth
    integral = amplitude * (growth - 1) / (MOMENT_RATE / 2)
    omitted = preparation_tail_upper(d)
    preparation = root_upper(2 * omitted, scale=scale)
    inverse_root = root_upper(F(1, d + 2), scale=scale)
    leakage_coefficient = (KAPPA + abs(COUPLING) * inverse_root) / (d + 2)
    residual = integral * leakage_coefficient
    projected = preparation + residual
    outside = final_amplitude / (d + 2)**2
    vector = root_upper(projected**2 + outside**2, scale=scale)
    return {
        "dimension": d,
        "total_model_time": str(t),
        "coupling_exact": str(COUPLING),
        "initial_fourth_number_moment_exact": str(moment),
        "fourth_number_moment_growth_rate": str(MOMENT_RATE),
        "final_fourth_number_moment_upper": str(final_amplitude**2),
        "preparation_omitted_probability_upper": str(omitted),
        "preparation_vector_error_upper": str(preparation),
        "projected_duhamel_residual_integral_upper": str(residual),
        "projected_history_vector_error_upper": str(projected),
        "final_outside_vector_norm_upper": str(outside),
        "pre_readout_stinespring_vector_error_upper": str(vector),
        "pre_readout_trace_distance_upper": str(min(F(1), vector)),
        "includes_intermediate_gauge_record_and_feedback": True,
        "uniform_over_declared_native_currents": True,
        "includes_final_sharp_field_readout": False,
        "includes_source_splitting_or_floating_execution": False,
        "includes_spatial_or_physical_service_error": False,
        "scope": "ideal interacting fixed two-site history before its final field readout",
    }


def sharp_readout_status():
    """Theorem scope, not a replacement for a numerical error certificate."""
    return {
        "complete_history_occupation_limit_proved": True,
        "proof": "OCCUPATION_REFINEMENT.md",
        "limiting_instrument": "positive/negative spectral projectors of infinite Q0",
        "finite_instrument": "positive/negative spectral projectors of compressed Q0",
        "pre_readout_quantitative_error_bound_supplied": True,
        "complete_history_finite_error_bound_supplied": False,
        "sharp_readout_rate": "not supplied",
        "compiled_dimensions": [2, 4],
        "limit_extends_to_arbitrary_even_occupation_dimensions": True,
        "spatial_continuum_or_empirical_result": False,
        "post_sharp_readout_fourth_moment_preservation_claimed": False,
    }


def report(dimensions=(4, 16, 64, 256)):
    return {"status": sharp_readout_status(),
            "pre_readout_bounds": [occupation_evolution_bound(d) for d in dimensions]}


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
