"""Blind signed phase inversion from complete native two-read records.

For the two source analyzers, the score in phase_readout has ideal means
A*cos(s*Delta) and A*sin(s*Delta), A=(9/2)*(b1-b0).  The earlier native
reference records certify A>0 and fix the ENACTED learned gap
Delta=log(bhat1/bhat0).  A is not fitted, nor is s passed to this inverse.

After any past records, a complete-trial half-trace allowance T gives mean
bias <=4T and second moment <=9q/20+8T for X in {-2,0,2}.  Thus
|X-E[X|past]|<=4.  For each fixed-length quadrature batch, martingale
Bernstein bounds both tails by 2 exp(-N*r^2/(2*(V+4r/3))).
The rational radius satisfies N*r^2>=2*k*(V+4r/3), so the two batches cost
at most 4*2^-k.  This conditional statement requires neither independent
input copies nor independence from training.  Source and sampler error
are conditional mean bias, not an error multiplied by the number of trials.

A declared parameter window must lie inside the right phase half-plane.
The confidence rectangle must also have strictly positive cosine.  Exact
rational atan enclosures then give an oriented phase interval; dividing
by the learned gap gives the implemented parameter.  Dividing separately
by the reference gap box is a comparison in reference units, not a claim
that the enacted generator was the unknown reference generator.

The proof is conditional on the declared source/sampler admission and
fixed design.  Counts do not authenticate execution.  This is active
source-word calibration, not an autonomous clock or a physical time unit.
"""
from __future__ import annotations

from fractions import Fraction as F

from . import common_caps as common, phase_readout as readout
from . import reference_tomography as tomography
from ..common_source.halfline import pi_interval

SCHEMA = "oph.native_modular_reference.phase_inference.v1"
QUADRATURES = ("cos", "sin")
ROUND_BITS = 64
ARCTAN_TERMS = 64
SCORE_CENTERED_BOUND = F(4)


def _exact(value, *, nonnegative=False):
    if type(value) is int:
        result = F(value)
    elif type(value) is F:
        numerator, denominator = value.numerator, value.denominator
        if (type(numerator) is not int or type(denominator) is not int
                or denominator <= 0):
            raise ValueError("an original exact rational with valid integer slots is required")
        result = F(numerator, denominator)
    else:
        raise ValueError("an original exact rational is required")
    if nonnegative and result < 0:
        raise ValueError("an original nonnegative exact rational is required")
    return result


def _outward(interval):
    lo, hi = interval
    return -tomography._up(-lo, ROUND_BITS), tomography._up(hi, ROUND_BITS)


def _small_atan(value):
    """Signed alternating-series enclosure, valid for |value|<=1/2."""
    if value < 0:
        lo, hi = _small_atan(-value)
        return -hi, -lo
    if value > F(1, 2):
        raise ValueError("the reduced atan argument exceeds one half")
    square, power, total = value*value, value, F()
    for k in range(ARCTAN_TERMS):
        total += (-1)**k*power/(2*k+1)
        power *= square
    next_term = power/(2*ARCTAN_TERMS+1)
    # 64 terms end with a negative term; the next term is positive.
    return total, total+next_term


def atan_interval(value):
    """Outward dyadic enclosure from exact reductions and Taylor tails."""
    value = _exact(value)
    if value < 0:
        lo, hi = atan_interval(-value)
        return -hi, -lo
    p_lo, p_hi = pi_interval()
    if value < F(1, 2):
        return _outward(_small_atan(value))
    if value <= 2:
        lo, hi = _small_atan((value-1)/(value+1))
        return _outward((p_lo/4+lo, p_hi/4+hi))
    lo, hi = _small_atan(1/value)
    return _outward((p_lo/2-hi, p_hi/2-lo))


def _divide(numerator, denominator):
    if denominator[0] <= 0 or denominator[0] > denominator[1]:
        raise ValueError("a positive denominator interval is required")
    corners = tuple(x/y for x in numerator for y in denominator)
    return _outward((min(corners), max(corners)))


def complete_alphabet(law="C", index=0):
    """All source and preparation branches, including zero observed counts."""
    return readout.record_alphabet(law, index, prepared=True)+(('numerical_abort',),)


def _counts(counts, law, index):
    if (type(counts) is not dict or set(counts) != set(QUADRATURES)
            or any(type(key) is not str for key in counts)):
        raise ValueError("both complete cos and sin count dictionaries are required")
    alphabet = complete_alphabet(law, index)
    result = {}
    for quadrature in QUADRATURES:
        column = counts[quadrature]
        if (type(column) is not dict or set(column) != set(alphabet)
                or any(type(key) is not tuple
                       or any(type(part) not in (int, str) for part in key)
                       for key in column)):
            raise ValueError("every original native record key, including zero branches, is required")
        if any(type(value) is not int or value < 0 for value in column.values()):
            raise ValueError("complete counts must be nonnegative original integers")
        n = sum(column.values())
        tomography._integer(n, 1, 1 << 40, "quadrature trial count")
        scores = {key: (0 if key == ('numerical_abort',) else readout.score(key))
                  for key in alphabet}
        total = sum(column[key]*scores[key] for key in alphabet)
        squared = sum(column[key]*scores[key]**2 for key in alphabet)
        result[quadrature] = {
            "trials": n, "score_sum": total, "score_square_sum": squared,
            "mean": F(total, n), "empirical_second_moment": F(squared, n),
            "complete_counts": tuple((key, column[key]) for key in alphabet)}
    return result


def score_budget(systematic_halftrace_upper, confidence_bits=12):
    """Source-uniform exact constants; observed variance cannot lower them."""
    tomography._integer(confidence_bits, 4, 64, "phase confidence bits")
    total = _exact(systematic_halftrace_upper, nonnegative=True)
    source = readout.complete_trial_budget()
    if total < source["complete_trial_halftrace_upper"] or total > 1:
        raise ValueError("the exact systematic allowance must cover the complete source budget and be at most one")
    return {
        "systematic_halftrace_upper": total,
        "source_complete_trial_halftrace_upper": source["complete_trial_halftrace_upper"],
        "score_mean_bias_upper": 4*total,
        "score_conditional_second_moment_upper": source["ideal_score_second_moment_upper"]+8*total,
        "score_centered_absolute_bound": F(SCORE_CENTERED_BOUND.numerator, SCORE_CENTERED_BOUND.denominator),
        "confidence_bits": confidence_bits,
        "phase_confidence_failure_upper": F(4, 2**confidence_bits),
        "fixed_lengths_and_analyzers_required": True,
        "uniform_conditional_source_and_sampler_admission_required": True}


def _radius(n, budget):
    variance, k = budget["score_conditional_second_moment_upper"], budget["confidence_bits"]
    radius = tomography._sqrt_upper(2*variance*k/n)+2*SCORE_CENTERED_BOUND*k/(3*n)
    if n*radius**2 < 2*k*(variance+SCORE_CENTERED_BOUND*radius/3):
        raise ArithmeticError("the exact martingale Bernstein exponent is insufficient")
    return tomography._up(radius, ROUND_BITS)


def trial_budget(score_radius, *, systematic_halftrace_upper, confidence_bits=12):
    """A sufficient finite trial count per quadrature for a total mean radius.

    The two Bernstein terms are bounded separately by half the requested
    statistical radius.  This is an a priori sufficient budget, not a fit
    to successful or unusually quiet records.
    """
    budget = score_budget(systematic_halftrace_upper, confidence_bits)
    requested = _exact(score_radius, nonnegative=True)
    remaining = requested-budget["score_mean_bias_upper"]
    if remaining <= 0:
        raise ValueError("the requested radius must exceed systematic mean bias")
    variance, k = budget["score_conditional_second_moment_upper"], budget["confidence_bits"]
    # Leave room for the upward square-root and final dyadic roundings.
    margin = F(1, 1 << 46)
    if remaining <= margin:
        raise ValueError("the requested statistical radius is below the exact rounding margin")
    target = remaining-margin
    requirements = (8*variance*k/target**2, 4*SCORE_CENTERED_BOUND*k/(3*target))
    n = max(1, *(-(-x.numerator//x.denominator) for x in requirements))
    tomography._integer(n, 1, 1 << 40, "sufficient quadrature trial count")
    attained = _radius(n, budget)+budget["score_mean_bias_upper"]
    if attained > requested:
        raise ArithmeticError("the sufficient finite budget did not attain its requested radius")
    return dict(budget, trials_per_quadrature=n, complete_trials=2*n,
                requested_score_radius=requested, attained_score_radius_upper=attained)


def _window(value):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError("a declared two-endpoint rational parameter window is required")
    lo, hi = (_exact(x) for x in value)
    if lo >= hi:
        raise ValueError("the declared parameter window must have positive width")
    return lo, hi


def _squared_modulus(box):
    lows, highs = [], []
    for lo, hi in box:
        lows.append(F() if lo <= 0 <= hi else min(lo*lo, hi*hi))
        highs.append(max(lo*lo, hi*hi))
    return sum(lows, F()), sum(highs, F())


def infer_counts(counts, reference_counts, *, law="C", index=0,
                 systematic_halftrace_upper, reference_sampling_bias=F(),
                 confidence_bits=12, reference_confidence_bits=10,
                 parameter_window=(F(-1), F(1))):
    """Return a blind oriented interval, or an explicit unresolved reason.

    The learned center ratio must be the ratio used by the admitted source
    word.  Replaying only this function does not certify that hardware or
    the numerical sampler implemented that word; orchestration checks it.
    No hidden parameter, fitted amplitude, normalized successful frequency
    or caller-supplied score enters this inverse.
    """
    observations = _counts(counts, law, index)
    window = _window(parameter_window)
    budget = score_budget(systematic_halftrace_upper, confidence_bits)
    reference = tomography.infer_counts(reference_counts,
        sampling_bias=reference_sampling_bias, confidence_bits=reference_confidence_bits)
    for observed in observations.values():
        statistical = _radius(observed["trials"], budget)
        radius = statistical+budget["score_mean_bias_upper"]
        observed.update(statistical_radius_upper=statistical,
                        total_mean_radius_upper=radius,
                        mean_interval=_outward((observed["mean"]-radius, observed["mean"]+radius)))
    result = {
        "schema": SCHEMA, "status": "unresolved", "law": law, "cap_index": index,
        "declared_parameter_window": window, "observed_scores": observations,
        "reference_inference": reference, "budget": budget,
        "confidence_failure_upper": reference["confidence_failure_upper"]+budget["phase_confidence_failure_upper"],
        "independent_input_copies_required": False, "unknown_amplitude_fitted": False,
        "successful_trials_conditioned_on": False,
        "source_or_hardware_execution_authenticated_by_counts": False,
        "scope": "active learned-source phase calibration on a declared window; no autonomous clock or physical time unit",
        "reference_unit_scope": "same measured phase divided by the true-reference gap box; the implemented generator uses the learned center"}
    if reference["status"] != "faithful":
        result["reason"] = "reference records do not certify the positive three-sector family"
        return result
    centers = reference["central_subnormalized_weights"]
    boxes = reference["subnormalized_weight_intervals"]
    amplitude = (F(9, 2)*(boxes[1][0]-boxes[0][1]),
                 F(9, 2)*(boxes[1][1]-boxes[0][0]))
    result["positive_amplitude_interval"] = amplitude
    if amplitude[0] <= 0:
        result["reason"] = "reference records do not certify the positive quadrature orientation"
        return result
    gap = _outward(common.logarithm_interval(centers[1]/centers[0]))
    reference_gap = _outward((common.logarithm_interval(boxes[1][0]/boxes[0][1])[0],
                               common.logarithm_interval(boxes[1][1]/boxes[0][0])[1]))
    result.update(implemented_gap_interval=gap, reference_gap_interval=reference_gap)
    if gap[0] <= 0 or reference_gap[0] <= 0:
        result["reason"] = "a positive implemented and reference modular gap is not certified"
        return result
    pi = _outward(pi_interval())
    result["score_mean_alias_period_interval"] = _divide((2*pi[0], 2*pi[1]), gap)
    result["alias_scope"] = "the selected score means repeat after 2*pi/gap; equality of complete instruments is not claimed"
    phase_window = (min(x*y for x in window for y in gap), max(x*y for x in window for y in gap))
    result["declared_phase_window_interval"] = _outward(phase_window)
    if phase_window[0] <= -pi[0]/2 or phase_window[1] >= pi[0]/2:
        result["reason"] = "the declared window is not certified inside one oriented right half-plane"
        return result
    cosine, sine = (observations[q]["mean_interval"] for q in QUADRATURES)
    modulus = _squared_modulus((cosine, sine))
    result["quadrature_squared_modulus_interval"] = modulus
    if modulus[0] > amplitude[1]**2 or modulus[1] < amplitude[0]**2:
        result["reason"] = "the quadrature box is incompatible with the trained amplitude family"
        return result
    if cosine[0] <= 0:
        result["reason"] = "the complete quadrature records do not certify a positive cosine"
        return result
    ratio = _divide(sine, cosine)
    phase = _outward((atan_interval(ratio[0])[0], atan_interval(ratio[1])[1]))
    raw_parameter = _divide(phase, gap)
    clipped = (max(raw_parameter[0], window[0]), min(raw_parameter[1], window[1]))
    result.update(phase_interval=phase, sine_cosine_ratio_interval=ratio,
                  unclipped_parameter_interval=raw_parameter,
                  reference_unit_parameter_interval=_divide(phase, reference_gap))
    if clipped[0] > clipped[1]:
        result["reason"] = "the phase confidence interval misses the declared parameter window"
        return result
    result.update(status="identified", parameter_interval=clipped,
                  parameter_interval_width=clipped[1]-clipped[0],
                  orientation=("positive" if phase[0] > 0 else "negative" if phase[1] < 0 else "unresolved"),
                  confidence_statement="simultaneous reference and phase coverage under the declared conditional source and fixed-design contracts")
    return result

