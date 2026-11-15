"""Finite native preparation of a faithful Palm modular reference.

The existing law-uniform binary warm-up gives D(rho_nR,tau tensor rho_R)
<= delta, tau=I/66. Apply the COMPLETE random-seam movement instrument,
keeping success A and failure C. Its equilibrium success is B=q*sigma,
q=5/33, sigma=diag(degree/600). Contractivity gives
||A-B||_1 + ||C-D||_1 <= 2 delta, hence |p-q|<=delta and
D(A/p,B/q)<=delta/q. This remains valid with an untouched reference.
The success probability is never hidden in a renormalized whole experiment.

For the source marginal, its traceless perturbation obeys operator norm
<=delta. Put t=66delta. Then (1-t)tau<=rho_n<=(1+t)tau. Positivity of
the actual herald gives (1-t)sigma<=A/q<=(1+t)sigma. If t<1, A/p is
faithful on M66. Operator monotonicity of log gives
||log(A/p)-log(sigma)+log(p/q) I|| <= -log(1-t) <= t/(1-t).
The scalar term cancels in modular conjugation. Duhamel's unitary bound
therefore gives half-diamond distance <= min(1, |s| t/(1-t)) between
their modular state actions, also with arbitrary probe reference systems.

Operator monotonicity here follows from
log X=int_0^infinity[(1+u)^-1 I-(X+u I)^-1]du and reversal of positive
matrix order by inversion. This proof does not require the two states to
commute. It instantiates finitePrivateStarAutomorphism_inner and the finite
Gibbs/KMS identity with H=-log(A/p), beta=1. The Gibbs Heisenberg parameter
is -s for sigma^(is) A sigma^(-is); no new dynamics is selected.

This is the full random-seam Palm herald, NOT the fixed-seam signed probe
herald, whose successful work support is only twenty dimensional. The
classical success/failure flag is retained separately; faithfulness refers
to the successful work marginal, never to every flag or external-reference
block. Warm-up counts native attempts; s is dimensionless modular parameter.
Programmed modular services remain admitted, and are not passive repair or
physical seconds. No full-78-state or continuum faithfulness is asserted.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from ..event_circulation import construction as events
from ..native_relational_readout import construction as relational
from ..observer import native
from ..predictions.repair_response import warmup_prior
from ..record_laws import burn_in

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "preparation_results.json"
Q = F(5, 33)
DIMENSION = 66


@lru_cache(maxsize=2)
def marked(law):
    return relational.source(law)


def moving_columns(law):
    columns = defaultdict(lambda: defaultdict(F))
    for pairs in marked(law).moving_maps:
        for i, out_i in pairs:
            for j, out_j in pairs:
                columns[i, j][out_i, out_j] += F(1, 60)
    return {pair: dict(row) for pair, row in columns.items()}


def anchor():
    """Bind the CP herald and its full faithful output to actual C/E events."""
    degree, sigma, _, _ = events.embedded()
    columns = moving_columns("C")
    if columns != moving_columns("E"):
        raise ArithmeticError("binary moving-current C/E maps differ")
    counts = marked("C").no_current_counts
    prepared = [F(0)] * DIMENSION
    for i in range(DIMENSION):
        for (a, b), value in columns.get((i, i), {}).items():
            if a != b:
                raise ArithmeticError("diagonal native herald created coherence")
            prepared[a] += value / DIMENSION
        if sum(columns.get((i, i), {}).values()) != F(degree[i], 60):
            raise ArithmeticError("native movement effect lost its degree count")
        for j in range(DIMENSION):
            total_trace = (F(counts[i][j], 60) if i == j else F(0))
            total_trace += sum(value for (a, b), value in columns.get((i, j), {}).items() if a == b)
            if total_trace != int(i == j):
                raise ArithmeticError("complete native movement/failure instrument lost trace")
    if tuple(prepared) != tuple(Q * p for p in sigma) or min(sigma) != F(1, 75):
        raise ArithmeticError("global movement herald does not prepare the full faithful Palm reference")
    return degree, sigma


def budget(attempts, modular_window=F(16), coverage=F(1)):
    """Exact conditional-reference and active-channel bounds from native work."""
    if type(modular_window) is not F or modular_window < 0:
        raise ValueError("nonnegative exact modular-parameter window required")
    delta = burn_in.binary_bound(attempts, coverage)
    t = DIMENSION * delta
    # The actual source effect has eigenvalues 8/60 and 10/60, even for
    # coherent inputs. These bounds also keep success strictly positive.
    p_lower = max(F(2, 15), Q - delta)
    p_upper = min(F(1, 6), Q + delta)
    faithful = t < 1
    return {
        "warmup_native_attempts": attempts,
        "total_native_attempts_per_trial": attempts + 1,
        "modular_parameter_window": modular_window,
        "exogenous_seam_coverage": F(coverage),
        "warmup_half_trace_error_upper": delta,
        "herald_success_interval": (p_lower, p_upper),
        "herald_failure_interval": (1-p_upper, 1-p_lower),
        "conditional_reference_half_trace_error_upper": min(F(1), delta / Q),
        "relative_operator_sandwich_radius": t,
        "successful_work_reference_faithful": faithful,
        "minimum_work_eigenvalue_lower": (1-t)*Q/(75*p_upper) if faithful else F(0),
        "generator_error_modulo_scalar_upper": t/(1-t) if faithful else None,
        "modular_channel_half_diamond_error_upper": min(F(1), modular_window*t/(1-t)) if faithful else None,
        "expected_attempts_with_full_warmup_repeated_each_trial_upper": F(attempts+1)/p_lower,
    }


def first_sufficient(error, modular_window=F(16), coverage=F(1)):
    if type(error) is not F or not 0 < error < 1:
        raise ValueError("exact modular-channel error in (0,1) required")
    for attempts in range(10001):
        row = budget(attempts, modular_window, coverage)
        bound = row["modular_channel_half_diamond_error_upper"]
        if bound is not None and bound < error:
            return attempts
    raise ValueError("no certificate within 10000 native attempts")


def diagonal_trial(initial, attempts):
    """Actual finite warm-up and specified-attempt herald, retaining failure."""
    prior = warmup_prior(initial, attempts)
    degree, sigma = anchor()
    columns = moving_columns("C")
    success = [F(0)] * DIMENSION
    failure = []
    for i, weight in enumerate(prior):
        for (a, b), value in columns.get((i, i), {}).items():
            if a != b:
                raise ArithmeticError("population trial created an off-diagonal output")
            success[a] += weight * value
        failure.append(weight * (1-F(degree[i], 60)))
    p = sum(success)
    normalized = tuple(value/p for value in success)
    delta = sum(abs(value-F(1, DIMENSION)) for value in prior)/2
    flagged = (sum(abs(value-Q*target) for value, target in zip(success, sigma, strict=True))
               + sum(abs(value-F(60-d, 60*DIMENSION))
                     for value, d in zip(failure, degree, strict=True)))/2
    conditional = sum(abs(value-target) for value, target in zip(normalized, sigma, strict=True))/2
    if sum(success)+sum(failure) != 1 or flagged > delta or conditional > delta/Q:
        raise ArithmeticError("finite flagged preparation violates its trace-distance bound")
    contract = budget(attempts)
    if delta > contract["warmup_half_trace_error_upper"]:
        raise ArithmeticError("actual native warm-up exceeds inherited certificate")
    if contract["successful_work_reference_faithful"]:
        t = contract["relative_operator_sandwich_radius"]
        if (min(normalized) < contract["minimum_work_eigenvalue_lower"]
                or any(not (1-t)*target <= value/Q <= (1+t)*target
                       for value, target in zip(success, sigma, strict=True))):
            raise ArithmeticError("actual herald violates the inherited operator sandwich")
    return {"prior": prior, "success": tuple(success), "failure": tuple(failure),
            "success_probability": p, "conditional": normalized,
            "warmup_error": delta, "flagged_error": flagged,
            "conditional_error": conditional}


def _upper(value, digits=12):
    scale = 10**digits
    return F(-(-(value*scale).numerator//(value*scale).denominator), scale)


def _lower(value, digits=12):
    return -_upper(-value, digits)


def _compact(row):
    # Exact rational budgets can have thousands of digits. Every displayed
    # bound is rounded outward; the replay recomputes the exact value first.
    lower_keys = {"minimum_work_eigenvalue_lower"}
    output = {}
    for key, value in row.items():
        if isinstance(value, F):
            output[key] = str((_lower if key in lower_keys else _upper)(value))
        elif isinstance(value, tuple):
            output[key] = [str(_lower(value[0])), str(_upper(value[1]))]
        else:
            output[key] = value
    return output


def run():
    degree, _ = anchor()
    n = first_sufficient(F(1, 1000))
    if n != 369:
        raise ArithmeticError("bounded-window preparation threshold changed")
    initial = tuple(int(port in (0, 2)) for port in range(12))
    trials = []
    for attempts in (0, 8, 221, n):
        row = diagonal_trial(initial, attempts)
        trials.append({"warmup_native_attempts": attempts,
                       "success_probability_interval": [str(_lower(row["success_probability"])),
                                                         str(_upper(row["success_probability"]))],
                       "source_half_trace_error_upper": str(_upper(row["warmup_error"])),
                       "complete_flagged_half_trace_error_upper": str(_upper(row["flagged_error"])),
                       "conditional_work_half_trace_error_upper": str(_upper(row["conditional_error"])),
                       "conditional_work_minimum_eigenvalue_lower": str(_lower(min(row["conditional"]))),
                       "complete_exact_trial_sha256": hashlib.sha256(canonical(row).encode()).hexdigest()})
    # A fixed seam movement herald cannot replace the full Palm instrument.
    states = relational.quantum.STATES
    fixed = {native.native(state, (0, 1), coin)[0]
             for state in states for coin in (0, 1)
             if native.native(state, (0, 1), coin)[1][1]}
    if len(fixed) != 20:
        raise ArithmeticError("fixed-seam probe-herald support changed")
    return {
        "schema": "oph.native_modular_reference.finite_preparation.v1",
        "global_C_E_moving_maps_identical": True,
        "complete_flag_instrument_matrix_units_checked": DIMENSION**2,
        "degree_census": {str(d): degree.count(d) for d in sorted(set(degree))},
        "ideal_success_probability": str(Q), "ideal_work_eigenvalue_floor": "1/75",
        "binary_quantum_warmup_uniform_over_untouched_references": True,
        "budget": _compact(budget(n)), "exact_diagonal_trials": trials,
        "fixed_seam_probe_herald_work_rank": len(fixed),
        "fixed_seam_probe_herald_is_full_faithful_reference_preparation": False,
        "finite_theorems": ["OPH.Dynamics.finitePrivateStarAutomorphism_inner",
                            "OPH.Dynamics.hamiltonianPropagator_mem_unitary",
                            "OPH.QFT.gibbsState_kms"],
        "scope": {"actual_successful_work_state_prepared_by_finite_native_events": True,
                  "full_success_failure_law_retained": True,
                  "faithfulness_claimed_for_external_reference_or_flag_algebra": False,
                  "warmup_tape_and_controller_unavailable_to_continuation": True,
                  "record_law_or_scheduler_selected_by_this_result": False,
                  "native_passive_repair_identified_with_modular_flow": False,
                  "active_modular_control_services_remain_admitted": True,
                  "physical_seconds_or_temperature_identified": False},
        "all_pass": True,
    }


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2,
                      default=lambda x: str(x) if isinstance(x, F) else x) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    encoded = canonical(run())
    if args.check is not None:
        if args.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit("finite native reference preparation receipt drift")
    else:
        with OUTPUT.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    print("finite native faithful-reference preparation: all checks passed")


if __name__ == "__main__":
    main()
