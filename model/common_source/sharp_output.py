"""Certified finite-code comparisons for complete sharp-sign output states.

All arithmetic is rational apart from alpha=1/sqrt(2*pi), enclosed by the
existing exact half-line certificate. Polynomial cancellations are performed
before interval evaluation. Old records, both signs and unnormalized branch
masses are retained. These are mathematical occupation/output comparisons,
not floating source-execution or physical-service error certificates.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

from . import halfline as h
from . import infinite_history as history
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT_SCALE = 10**60
DIMENSIONS = (2, 4, 8, 16, 64, 256)
_ZERO = (F(0), F(0))
SOURCE_PATHS = (
    "model/common_source/sharp_output.py",
    "model/common_source/sharp_pointer.py",
    "model/common_source/halfline.py",
    "model/common_source/infinite_history.py",
    "model/common_source/finite_signal.py",
    "model/common_source/experiment.py",
    "model/common_source/protocol.py",
    "model/common_source/scalar_coupling.py",
    "model/common_source/gauge_matter.py",
    "model/common_source/occupation.py",
    "model/observer/source_scalar/compiler.py",
    "model/observer/source_scalar/experiment.py",
    "model/observer/source_scalar/field.py",
    "model/observer/source_scalar/bounds.py",
    "model/observer/source_scalar/cutoff.py",
    "code/source_selection_model/response.json",
    "code/m1_source_realization/receipt.json",
    "model/requirements.txt",
)


def _dimension(value):
    if type(value) is not int or not 1 <= value <= h.MAX_OCCUPATION + 1:
        raise ValueError("a finite code dimension from 1 through 513 required")
    return value


def _outcome(value):
    if type(value) is not int or value not in (0, 1):
        raise ValueError("an integer sign outcome zero or one required")
    return value


def _normalization(value):
    if type(value) is not F or value <= 0:
        raise ValueError("a positive exact rational initial normalization required")
    return value


def _polynomial(value):
    if (type(value) not in (tuple, list) or not value
            or any(type(coefficient) is not F for coefficient in value)):
        raise ValueError("nonempty exact rational polynomial coefficients required")
    return tuple(value)


def _add(left, right):
    return tuple((left[k] if k < len(left) else F(0))
                 + (right[k] if k < len(right) else F(0))
                 for k in range(max(len(left), len(right))))


def _scale(value, factor):
    return tuple(factor * coefficient for coefficient in value)


def _multiply(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return tuple(result)


def _absolute_squared(value):
    real, imag = tuple(part[0] for part in value), tuple(part[1] for part in value)
    return _add(_multiply(real, real), _multiply(imag, imag))


def polynomial_interval(coefficients):
    """Enclose sum_k coefficients[k]*alpha**k, alpha=1/sqrt(2*pi).

    Exact coefficients may have either sign. In particular, an outside
    mass or a rank-one distance radicand is subtracted symbolically before
    this function is called; no floating subtraction certifies a tail.
    """
    coefficients = _polynomial(coefficients)
    low, high = h.inverse_sqrt_two_pi_interval()
    lower = upper = F(0)
    for power, coefficient in enumerate(coefficients):
        a, b = low**power, high**power
        if coefficient < 0:
            a, b = b, a
        lower += coefficient*a
        upper += coefficient*b
    return lower, upper


def _nonnegative_interval(coefficients):
    lower, upper = polynomial_interval(coefficients)
    if upper < 0:
        raise ArithmeticError("proved nonnegative polynomial has a negative enclosure")
    return max(F(0), lower), upper


def root_interval(value):
    """Outward rational square root, retaining zero and arbitrarily small mass."""
    if (type(value) not in (tuple, list) or len(value) != 2
            or any(type(part) is not F for part in value)
            or not 0 <= value[0] <= value[1]):
        raise ValueError("a nonnegative ordered exact rational interval required")
    lower = F(math.isqrt(value[0].numerator * ROOT_SCALE**2 // value[0].denominator),
              ROOT_SCALE)
    upper = F(math.isqrt(value[1].numerator * ROOT_SCALE**2 // value[1].denominator),
              ROOT_SCALE)
    if upper*upper < value[1]:
        upper += F(1, ROOT_SCALE)
    return lower, upper


def interval_strings(value, places=15):
    """Short outward-rounded receipt values; internal calculations stay exact."""
    if (type(value) not in (tuple, list) or len(value) != 2
            or any(type(part) is not F for part in value) or value[0] > value[1]
            or type(places) is not int or not 0 <= places <= 60):
        raise ValueError("an ordered rational interval and decimal places 0..60 required")
    scale = 10**places
    low = F(value[0].numerator*scale // value[0].denominator, scale)
    high = F(-((-value[1].numerator*scale) // value[1].denominator), scale)
    return [str(low), str(high)]


def _groups(vector):
    groups = {}
    for (n0, n1, gauge), coefficient in h._vector(vector).items():
        groups.setdefault((n1, gauge), {})[n0] = coefficient
    return groups


@lru_cache(maxsize=65536)
def _projection_row(coefficients, m):
    """Immutable exact <m_B|P_positive|u> coefficients, reused across cutoffs."""
    amplitude = [[F(0), F(0)], [F(0), F(0)]]
    for n, z in coefficients:
        r, s = h.halfline_coefficient(m, n)
        for power, factor in enumerate((r, s)):
            for part in (0, 1):
                amplitude[power][part] += factor*z[part]
    return tuple(tuple(part) for part in amplitude)


def halfline_overlap_coefficients(left, right, outcome=1):
    """Return Gaussian coefficients (z0,z1) of <left|P_b|right>=z0+alpha*z1.

    Both vectors use the unnormalized Bargmann basis. Spectator occupations
    and gauge labels are orthogonal; factorial weights are retained.
    """
    outcome = _outcome(outcome)
    first, second = _groups(left), _groups(right)
    result = [[F(0), F(0)], [F(0), F(0)]]
    for (n1, gauge), rows in first.items():
        coefficients = tuple(sorted(rows.items()))
        for m, (c, d) in second.get((n1, gauge), {}).items():
            for power, (a, b) in enumerate(_projection_row(coefficients, m)):
                factor = math.factorial(n1)*(1 if power == 0 else 2*outcome-1)
                result[power][0] += factor*(a*c+b*d)
                result[power][1] += factor*(a*d-b*c)
    return tuple(tuple(part) for part in result)


def projected_masses(vector, dimension, outcome, normalization=F(1)):
    """Exact coefficient triples for the total, inside and outside sign mass.

    The code projection keeps n0,n1 < dimension and every gauge value.
    P_b is the infinite coordinate projector, never sign of compressed Q.
    Division uses the supplied initial norm only; branch masses are untouched.
    """
    dimension, outcome = _dimension(dimension), _outcome(outcome)
    normalization = _normalization(normalization)
    groups = _groups(vector)
    constant, cross = h.positive_weight_coefficients(vector, normalization)
    total = (constant, (2*outcome-1)*cross, F(0))
    inside = (F(0), F(0), F(0))
    for (n1, _), coefficients in groups.items():
        if n1 >= dimension:
            continue
        entries = tuple(sorted(coefficients.items()))
        for m in range(dimension):
            positive = _projection_row(entries, m)
            amplitude = (positive[0], tuple((2*outcome-1)*part for part in positive[1]))
            # <m_B|P_b|u>/m! is the output Bargmann coefficient.
            # Multiplying its squared norm by m! n1! leaves n1!/m!.
            inside = _add(inside, _scale(_absolute_squared(amplitude),
                         F(math.factorial(n1), math.factorial(m))/normalization))
    outside = _add(total, _scale(inside, F(-1)))
    return {"total": total, "inside": inside, "outside": outside}


def support_witness(vector, dimension, outcome, normalization=F(1)):
    """An unconditional lower bound against every code-supported output block.

    If a,t are the inside/outside squared norms of v=P_b u, the norm-one
    witness coupling those two normalized components gives
    (1/2)|| |v><v|-sigma ||_1 >= sqrt(a*t), for every sigma in the code.
    Sigma may have any branch mass and need not be pure.
    """
    masses = projected_masses(vector, dimension, outcome, normalization)
    enclosed = {name: _nonnegative_interval(value) for name, value in masses.items()}
    inside, outside = enclosed["inside"], enclosed["outside"]
    witness = root_interval((inside[0]*outside[0], inside[1]*outside[1]))
    return {"masses": enclosed, "half_trace_distance_witness_interval": witness}


def pure_output_distance_interval(infinite_vector, finite_vector, outcome,
                                  normalization=F(1), finite_normalization=F(1)):
    """Exact READOUT distance between P_b u and an already post-read vector w.

    Return the half trace norm of |P_b u><P_b u|-|w><w|, with outward
    arithmetic. Relative phase, unequal raw mass, and the infinite output
    tail are all included. The finite vector need not be P_b-invariant.
    """
    outcome = _outcome(outcome)
    normalization = _normalization(normalization)
    finite_normalization = _normalization(finite_normalization)
    r, s = h.positive_weight_coefficients(infinite_vector, normalization)
    p = (r, (2*outcome-1)*s)
    q = h.norm_squared(finite_vector)/finite_normalization
    overlap = halfline_overlap_coefficients(infinite_vector, finite_vector, outcome)
    summed = _add(p, (q,))
    radicand = _add(_multiply(summed, summed),
                    _scale(_absolute_squared(overlap), -F(4)/(normalization*finite_normalization)))
    low, high = root_interval(_nonnegative_interval(radicand))
    return low/2, high/2


def retained_history(certificate):
    """Validate complete retained polynomial records and their imported budgets.

    The fixed public run regenerates this certificate from its source. This
    parser also binds the supplied coefficient hashes, exact masses, record
    grid and recomputed error summary; it does not accept missing branches.
    """
    if (type(certificate) is not dict
            or certificate.get("schema") != "oph.common_source.infinite_history.v1"):
        raise ValueError("the complete infinite-history polynomial certificate required")
    normalization = F(certificate["initial_polynomial_norm_squared"])
    _normalization(normalization)
    configurations = certificate.get("configurations")
    if type(configurations) is not list or len(configurations) != 2:
        raise ValueError("both feedback configurations required")
    preparation = F(certificate["initial_trace_distance_upper"])
    first_errors = {key: F(value) for key, value in certificate["first_segment_vector_errors"].items()}
    error_keys = {"Taylor", "discarded_support", "coefficient_rounding"}
    if (set(first_errors) != error_keys or preparation < 0
            or any(value < 0 for value in first_errors.values())):
        raise ValueError("nonnegative complete inherited vector errors required")
    result = {}
    for configuration in configurations:
        enabled = configuration.get("feedback_enabled")
        rows = configuration.get("retained_first_records")
        if (type(enabled) is not bool or enabled in result
                or type(rows) is not list or len(rows) != 2):
            raise ValueError("two unique feedback settings with both first records required")
        branches, second_errors = {}, []
        for row in rows:
            first = row.get("first_record")
            data = row.get("pre_readout_polynomial")
            if (type(first) is not int or first not in (0, 1) or first in branches
                    or type(data) is not list):
                raise ValueError("unique first records and retained polynomial coefficients required")
            digest = hashlib.sha256(json.dumps(data, separators=(",", ":")).encode()).hexdigest()
            if digest != row.get("pre_readout_polynomial_sha256"):
                raise ValueError("retained polynomial hash mismatch")
            vector = {}
            for item in data:
                if (type(item) is not list or len(item) != 5
                        or any(type(value) is not int for value in item[:3])
                        or any(type(value) is not str for value in item[3:])):
                    raise ValueError("exact retained Bargmann coefficient rows required")
                key = tuple(item[:3])
                if key in vector:
                    raise ValueError("duplicate retained polynomial coordinate")
                vector[key] = tuple(map(F, item[3:]))
            h._vector(vector)
            if (len(data) != row.get("pre_readout_nonzero_coefficients")
                    or h.norm_squared(vector)/normalization != F(row["approximate_pre_readout_norm_squared"])):
                raise ValueError("retained polynomial mass or coefficient count mismatch")
            errors = {key: F(value) for key, value in row["second_segment_vector_errors"].items()}
            if set(errors) != error_keys or any(value < 0 for value in errors.values()):
                raise ValueError("all nonnegative second-segment errors required")
            second_errors.append(sum(errors.values(), F(0)))
            branches[first] = vector
        error = sum(first_errors.values(), F(0)) + history._root_upper(
            sum((value*value for value in second_errors), F(0)))
        complete_error = preparation + error + error*error/2
        if (error != F(configuration["normalized_preparation_to_polynomial_vector_error_upper"])
                or complete_error != F(configuration["complete_output_half_trace_norm_error_upper"])):
            raise ValueError("inherited complete-output error summary mismatch")
        result[enabled] = {"branches": branches, "complete_output_error": complete_error}
    return normalization, result


def history_support_bounds(certificate, dimensions=DIMENSIONS):
    """Finite-code lower bounds for all four records of both actual histories."""
    if (type(dimensions) not in (tuple, list) or not dimensions
            or len(set(dimensions)) != len(dimensions)):
        raise ValueError("nonempty distinct code dimensions required")
    for dimension in dimensions:
        _dimension(dimension)
    normalization, configurations = retained_history(certificate)
    output = []
    for enabled in (True, False):
        configuration = configurations[enabled]
        rows = []
        for dimension in dimensions:
            witnesses, lower = [], F(0)
            for first in (0, 1):
                for outcome in (0, 1):
                    result = support_witness(configuration["branches"][first], dimension,
                                             outcome, normalization)
                    local = result["half_trace_distance_witness_interval"][0]
                    lower += local
                    witnesses.append({"record": [first, outcome],
                        "polynomial_sign_mass_interval": interval_strings(result["masses"]["total"]),
                        "polynomial_inside_mass_interval": interval_strings(result["masses"]["inside"]),
                        "polynomial_outside_mass_interval": interval_strings(result["masses"]["outside"]),
                        "polynomial_block_half_trace_distance_lower": str(interval_strings((local, local))[0])})
            physical_lower = max(F(0), lower-configuration["complete_output_error"])
            if physical_lower > 1:
                raise ArithmeticError("normalized-history distance lower bound exceeds one")
            rows.append({"onsite_code_dimension": dimension,
                         "retained_records": witnesses,
                         "polynomial_complete_output_lower": interval_strings((lower, lower))[0],
                         "infinite_history_complete_output_lower": interval_strings((physical_lower, physical_lower))[0]})
        output.append({"feedback_enabled": enabled,
                       "inherited_complete_output_error_upper": str(configuration["complete_output_error"]),
                       "dimensions": rows})
    return output


def vacuum_certificate(dimensions=DIMENSIONS):
    """Independent vacuum tail identity, retaining both equally likely outcomes."""
    vacuum = {(0, 0, 0): (F(1), F(0))}
    rows = []
    for dimension in dimensions:
        _dimension(dimension)
        if dimension % 2:
            raise ValueError("the vacuum closed sum uses an even code dimension")
        coefficient = sum((F(math.comb(2*k, k), 4**k*(2*k+1))
                           for k in range(dimension//2)), F(0))
        tail_coefficients = (F(1, 4), F(0), -coefficient)
        tail = _nonnegative_interval(tail_coefficients)
        positive = support_witness(vacuum, dimension, 1)
        negative = support_witness(vacuum, dimension, 0)
        if projected_masses(vacuum, dimension, 1)["outside"] != tail_coefficients:
            raise ArithmeticError("general projection disagrees with independent vacuum tail")
        lower = (positive["half_trace_distance_witness_interval"][0]
                 + negative["half_trace_distance_witness_interval"][0])
        # The actual finite sign(Q_d) vacuum effects have equal probabilities.
        # This countercontrol concerns complete outputs, not probability error.
        rows.append({"onsite_code_dimension": dimension,
                     "positive_and_negative_probabilities": ["1/2", "1/2"],
                     "each_sign_outside_mass_interval": interval_strings(tail),
                     "complete_output_half_trace_distance_lower": interval_strings((lower, lower))[0]})
    return rows


def finite_two_level_comparison(certificate, degree=64, steps=512):
    """Bound actual complete d=2 continuous and ideal source-word discrepancies.

    Every approximate finite branch is compared directly with the infinite
    half-line polynomial branch, including different initial preparations.
    Inherited approximation errors widen the resulting complete-output
    interval; probabilities are never substituted for quantum output states.
    """
    from . import finite_signal
    normalization, configurations = retained_history(certificate)
    if certificate["parameters"]["coupling"] != "1/4":
        raise ValueError("the finite d=2 comparison uses the declared coupling 1/4")
    segments, split_error = finite_signal._source_word_errors(steps)
    output = []
    for enabled in (True, False):
        branches, _, _, eta = finite_signal.branch_vectors(enabled, degree)
        # Gauge-record writing, controlled feedback and the complete sign
        # pointer are isometries on the coherently flagged history. The two
        # Taylor errors therefore give the same GLOBAL vector bound eta,
        # before classical commitment. Its half-density error is at most
        # eta+eta²/2; committing every record contracts that distance.
        finite_error = eta+eta*eta/2
        total = (F(0), F(0))
        rows = []
        for (first, outcome), vector in branches:
            finite = {((index >> 2) & 1, (index >> 3) & 1, index & 3): amplitude
                      for index, amplitude in enumerate(vector)}
            interval = pure_output_distance_interval(configurations[enabled]["branches"][first],
                                                    finite, outcome, normalization)
            total = tuple(a+b for a, b in zip(total, interval, strict=True))
            rows.append({"record": [first, outcome],
                         "polynomial_block_half_trace_distance_interval": interval_strings(interval)})
        error = finite_error+configurations[enabled]["complete_output_error"]
        continuous = (max(F(0), total[0]-error), min(F(1), total[1]+error))
        ideal_word = (max(F(0), continuous[0]-split_error), min(F(1), continuous[1]+split_error))
        output.append({"feedback_enabled": enabled, "retained_records": rows,
                       "finite_polynomial_complete_output_error_upper": interval_strings((finite_error, finite_error), 60)[1],
                       "infinite_polynomial_complete_output_error_upper": str(configurations[enabled]["complete_output_error"]),
                       "polynomial_complete_output_distance_interval": interval_strings(total),
                       "continuous_complete_output_distance_interval": interval_strings(continuous),
                       "ideal_source_word_complete_output_distance_interval": interval_strings(ideal_word)})
    return {"onsite_code_dimension": 2, "finite_Taylor_degree": degree,
            "source_steps_per_segment": steps, "finite_source_segment_certificates": segments,
            "finite_source_complete_output_error_upper": str(split_error),
            "configurations": output,
            "comparison": "same declared current-zero histories, both first records and both signs, including different finite/infinite preparations",
            "numerical_execution_or_physical_service_error_included": False}


def run():
    """Regenerate the actual infinite histories before every fixed comparison."""
    from . import sharp_pointer
    certificate = history.certificate(retain_polynomials=True)
    root = Path(__file__).resolve().parents[2]
    return {"schema": "oph.common_source.sharp_output.v1",
            "source_sha256": bound_pins({path: hashlib.sha256((root/path).read_bytes()).hexdigest()
                              for path in SOURCE_PATHS}),
            "infinite_polynomial_certificate_sha256": hashlib.sha256(
                history.canonical(certificate).encode()).hexdigest(),
            "infinite_history_parameters": certificate["parameters"],
            "vacuum_countercontrol": vacuum_certificate(),
            "actual_history_code_support_bounds": history_support_bounds(certificate),
            "actual_two_level_complete_output_comparison": finite_two_level_comparison(certificate),
            "larger_code_complete_output_comparisons": [
                sharp_pointer.complete_cutoff_comparison(certificate, d) for d in (32, 64, 128)],
            "scope": {"native_preparation_current": 0,
                      "all_first_and_sign_records_retained": True,
                      "branch_renormalization": False,
                      "lower_bound_applies_to_every_output_in_the_named_occupation_code": True,
                      "optimality_of_finite_sign_instrument_claimed": False,
                      "all_cutoff_asymptotic_rate_claimed": False,
                      "source_execution_above_existing_dimensions_claimed": False,
                      "physical_services_or_spatial_refinement_certified": False}}


def main():
    import argparse
    from tools.strict_json import load
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create an exact receipt without overwriting")
    group.add_argument("--check", type=Path, help="recompute and compare every exact receipt field")
    args = parser.parse_args()
    expected = load(args.check) if args.check else None
    result = run()
    if args.check:
        history.compare_receipt(expected, result)
    if args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(history.canonical(result))
    print(history.canonical({"schema": result["schema"],
                             "complete_two_level_comparisons": [row["ideal_source_word_complete_output_distance_interval"]
                                 for row in result["actual_two_level_complete_output_comparison"]["configurations"]]}))


if __name__ == "__main__":
    main()
