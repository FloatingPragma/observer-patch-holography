"""Reference identification by complete native records on one cross seam.

After the existing local preparation and one unread core attempt, read the
native signed current on (1,7).  Write b_q=t*a_q for the three charge-sector
weights of the ideal successful reference, retaining its success mass t.
The three success probabilities obey u=(9*b1+b2)/2, v=(9*b0+b1)/2,
t=45*b0+20*b1+b2.  Consequently

 b0=(22*v+2*u-t)/54, b1=(t-10*v-2*u)/6,
 b2=5*u+15*v-3*t/2.

Neither t nor the reference weights are inputs to the inverse.  Its ratios
determine the modular action; the positive normalized weights determine the
two-moment MaxEnt family in common_variation.  This is reduced-family
identification, not tomography of an unconstrained coherent state.

For each score, |X-E[X|past]|<=L and E[X^2|past]<=V.  The existing complete,
arbitrary-reference preparer bounds every conditional record distribution,
even on a reused carrier with an initially correlated auxiliary.  Bernstein's
martingale bound gives 2 exp(-N*r^2/(2*(V+L*r/3))).  We choose a rational r
with N*r^2>=2*k*(V+L*r/3), so three scores fail with probability <=6*2^-k.
Systematic source error is added once to the mean, not N times.  Individual
records need not be independent, identically distributed, or monotone.

The retained execution samples the actual finite population-map controller,
core repair and current instrument, keeping every terminal flag and state.
Its binary CDF is converted to an explicit dyadic Markov kernel.  An exact
column audit certifies THIS sampler's bias and variance against the ideal
target; floating population maps are separately numerical diagnostics.  A
common-part decomposition of that same kernel accelerates the continuing
sample without substituting an ideal multinomial or an independent reset.
No microscopic repair tape or physical service execution is claimed.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import common_caps as common, local_preparation as local
from ..native_equilibration import common_variation as variation
from ..observer import native

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("reference_tomography_results.json")
SEAM = (1, 7)
LABELS = ("success_-1", "success_0", "success_1", "herald_failure",
          "warmup_timeout", "completion_timeout", "numerical_abort")
DIMENSIONS = (45, 20, 1)
GRID = 1 << 53
SAMPLER_GATE = F(1, 10**10)
SCORES = ((F(7, 18), -F(1, 54), F(1, 54), F(), F(), F(), F()),
          (-F(3, 2), F(1, 6), -F(1, 6), F(), F(), F(), F()),
          (F(27, 2), -F(3, 2), F(7, 2), F(), F(), F(), F()))
WIDTHS = (F(11, 27), F(5, 3), F(15))
SECTOR_SECOND_MOMENT_MAX = (F(5, 324), F(1, 4), F(81, 4))


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an original integer in [{low},{high}]")


def _rational(value, *, nonnegative=True):
    if type(value) not in (int, F) or (nonnegative and value < 0):
        raise ValueError("an exact nonnegative rational is required")
    return F(value)


def _up(value, bits=48):
    value = F(value)
    scale = 1 << bits
    return F(-(-value.numerator*scale//value.denominator), scale)


def _sqrt_upper(value):
    scale = 1 << 48
    floor = math.isqrt(value.numerator*scale*scale//value.denominator)
    return F(floor+1, scale)


@lru_cache(maxsize=1)
def _population_data():
    """Twice each actual native population map, with all current branches."""
    basis = native.states(binary=True)
    lookup = {x: i for i, x in enumerate(basis)}
    core = np.zeros((66, 66), dtype=np.int64)
    reads = {j: np.zeros((66, 66), dtype=np.int64) for j in (-1, 0, 1)}
    for i, x in enumerate(basis):
        for coin in (0, 1):
            y, _ = native.native(x, local.CORE, coin)
            core[lookup[y], i] += 1
            y, (_, j) = native.native(x, SEAM, coin)
            reads[j][lookup[y], i] += 1
    if not np.array_equal(core.sum(axis=0), np.full(66, 2)) or not np.array_equal(
            sum(reads.values()).sum(axis=0), np.full(66, 2)):
        raise ArithmeticError("the complete native instrument lost a column")
    return core.tobytes(), tuple(reads[j].tobytes() for j in (-1, 0, 1))


def population_counts():
    """Fresh arrays: callers cannot mutate the cached native proof object."""
    core, reads = _population_data()
    unpack = lambda data: np.frombuffer(data, dtype=np.int64).reshape(66, 66).copy()
    return unpack(core), {j: unpack(data) for j, data in zip((-1, 0, 1), reads, strict=True)}


@lru_cache(maxsize=1)
def _record_contract():
    """Source budget scalars without evaluating any target reference weights.

    These are the record-bound and resource parts of local.budget at its
    declared default settings.  Its additional log/faithfulness sandwich
    needs the known reference; record inference needs none of that data.
    Equality to those existing budget fields is separately tested.
    """
    target, deadline = F(1, 2**24), 256
    selection = local.selectors(32)
    _, completion_fraction, gamma = local._schedule()
    warm_fraction = F(len(local.incident.incident_edges(local.CORE)), 30)
    warm_blocks = complete_blocks = 0
    while F(33)/warm_fraction/2**warm_blocks > target:
        warm_blocks += 1
    while F(33)/completion_fraction/2**complete_blocks > target:
        complete_blocks += 1
    warm_ticks, complete_ticks = 16*warm_blocks, math.ceil(1/gamma)*complete_blocks
    warm_slots, complete_slots = (warm_ticks+1)*deadline, complete_ticks*deadline
    warm_loss = local._dyadic_upper((warm_ticks+1)*(deadline*F(selection["warmup"]["timeout_probability"])
                                                  +(1-warm_fraction)**deadline))
    complete_loss = local._dyadic_upper(complete_ticks*(deadline*selection["completion_timeout_probability"]
                                                        +(1-completion_fraction)**deadline))
    q = F(5, 33)
    bound = local.incident.local_bound(warm_ticks, local.CORE)+warm_loss
    bound += q*(F(33)/completion_fraction/2**complete_blocks+complete_loss)
    total = warm_slots*(1+selection["warmup"]["maximum_auxiliary_native_attempts"])
    total += complete_slots*(1+selection["completion_auxiliary_calls_per_slot_upper"])
    return bound, q, total


def exact_identification():
    """Check the actual effects on each charge sector, not a fitted table."""
    _, reads = population_counts()
    charges = variation.charges()
    effects = {j: tuple(F(int(x), 2) for x in a.sum(axis=0)) for j, a in reads.items()}
    responses, second = [], []
    for q, g in enumerate(DIMENSIONS):
        law = tuple(sum((effects[j][i] for i, qi in enumerate(charges) if qi == q), F())/g
                    for j in (-1, 0, 1))+(F(),)*4
        estimates = tuple(sum((s*p for s, p in zip(score, law, strict=True)), F())
                          for score in SCORES)
        if estimates != tuple(F(int(k == q), g) for k in range(3)):
            raise ArithmeticError("native current effects do not invert the charge family")
        moments = tuple(sum((s*s*p for s, p in zip(score, law, strict=True)), F())
                        for score in SCORES)
        responses.append(law[:3]); second.append(moments)
    if tuple(max(row[i] for row in second) for i in range(3)) != SECTOR_SECOND_MOMENT_MAX:
        raise ArithmeticError("source-uniform second moment certificate changed")
    # Equal cross-seam occupations can conceal different protected charges.
    basis = native.states(binary=True)
    hidden = [tuple(int(p in occupied) for p in range(12)) for occupied in ((0, 1), (1, 2))]
    hidden_laws = [tuple(effects[j][basis.index(x)] for j in (-1, 0, 1)) for x in hidden]
    if hidden_laws[0] != hidden_laws[1] or sum(hidden[0][p] for p in local.CORE) == sum(hidden[1][p] for p in local.CORE):
        raise ArithmeticError("uncompleted-state nonidentifiability countercontrol changed")
    return {"seam": SEAM, "charge_sector_current_probabilities": responses,
            "score_second_moment_maxima": SECTOR_SECOND_MOMENT_MAX,
            "unknown_success_mass_cancels_in_modular_ratios": True,
            "source_laws": ["C", "E"], "full_state_tomography": False,
            "uncompleted_countercontrol": {"occupied_ports": ((0, 1), (1, 2)),
                                           "same_current_law": hidden_laws[0], "different_core_charges": (2, 1)}}


def infer_counts(counts, *, sampling_bias=F(), confidence_bits=10):
    """Target-blind exact confidence box from all complete-trial counts.

    sampling_bias=0 is the physical finite-instrument statement.  A supplied
    positive value must be established by a separate sampler audit; the
    retained producer computes it from every exact dyadic transition column.
    The source herald probability enters only a variance bound, not inversion.
    """
    if type(counts) is not dict or set(counts) != set(LABELS):
        raise ValueError("all seven complete terminal count categories are required")
    if any(type(v) is not int or v < 0 for v in counts.values()):
        raise ValueError("terminal counts must be nonnegative original integers")
    n = sum(counts.values())
    _integer(n, 1, 1 << 40, "trial count")
    _integer(confidence_bits, 4, 64, "confidence bits")
    sampling_bias = _rational(sampling_bias)
    if sampling_bias > 1:
        raise ValueError("sampling bias must be at most one")
    source_error, q, native_cost = _record_contract()
    error = source_error+sampling_bias
    estimates = tuple(sum((score[i]*counts[label] for i, label in enumerate(LABELS)), F())/n
                      for score in SCORES)
    radii, variances = [], []
    for score, width, vmax in zip(SCORES, WIDTHS, SECTOR_SECOND_MOMENT_MAX, strict=True):
        variance = q*vmax+error*max(s*s for s in score)
        radius = _sqrt_upper(2*variance*confidence_bits/n)+2*width*confidence_bits/(3*n)
        if n*radius*radius < 2*confidence_bits*(variance+width*radius/3):
            raise ArithmeticError("rational Bernstein exponent is insufficient")
        radii.append(_up(radius+width*error)); variances.append(variance)
    intervals = tuple((value-radius, value+radius) for value, radius in zip(estimates, radii, strict=True))
    total = sum((g*w for g, w in zip(DIMENSIONS, estimates, strict=True)), F())
    if total != F(sum(counts[x] for x in LABELS[:3]), n):
        raise ArithmeticError("inversion does not retain observed success mass")
    result = {"status": "unresolved", "trials": n,
              "central_subnormalized_weights": estimates,
              "subnormalized_weight_intervals": intervals,
              "observed_success_mass": total, "systematic_halftrace_upper": error,
              "score_conditional_second_moment_upper": tuple(variances),
              "confidence_failure_upper": F(6, 2**confidence_bits),
              "confidence_bits": confidence_bits,
              "source_native_attempts_per_trial_upper": native_cost+2,
              "source_native_attempts_all_trials_upper": n*(native_cost+2),
              "independent_input_copies_required": False,
              "source_or_hardware_execution_authenticated_by_counts": False}
    if total <= 0 or any(lo <= 0 for lo, _ in intervals):
        result["reason"] = "complete records do not yet certify three positive weights"
        return result
    normalized = tuple(w/total for w in estimates)
    mean, pair = 20*normalized[1]+2*normalized[2], normalized[2]
    if variation.family(mean, pair) != tuple(normalized[q] for q in variation.charges()):
        raise ArithmeticError("the inferred positive reference does not inhabit the MaxEnt family")
    log_error = F()
    for center, (lo, hi) in zip(estimates, intervals, strict=True):
        log_error = max(log_error, -common.logarithm_interval(lo/center)[0],
                        common.logarithm_interval(hi/center)[1])
    result.update({"status": "faithful", "central_weights": normalized,
                   "moments": (mean, pair),
                   "adjacent_modular_ratio_intervals": tuple(
                       (intervals[i+1][0]/intervals[i][1], intervals[i+1][1]/intervals[i][0])
                       for i in range(2)),
                   "logarithm_operator_error_upper": _up(log_error),
                   "modular_halfdiamond_upper_at_one": min(F(1), _up(log_error)),
                   "normalization_scalar_ignored_only_in_modular_channel": True})
    return result


def _raw_kernel(law):
    maps = local.population_maps(law=law)
    core, reads = population_counts()
    success = (core/2)@maps["success"]
    arrays = [(reads[j]/2)@success for j in (-1, 0, 1)]
    arrays += [maps[flag] for flag in LABELS[3:6]]
    kernel = np.vstack(arrays)
    if kernel.shape != (396, 66) or not np.isfinite(kernel).all() or np.any(kernel < 0):
        raise ArithmeticError("original complete population kernel is invalid")
    residual = float(np.max(np.abs(kernel.sum(axis=0)-1)))
    if residual > 5e-10:
        raise ArithmeticError("finite population-kernel mass diagnostic failed")
    return kernel, residual


def quantized_kernel(law="C"):
    """Explicit 53-bit CDF sampler, including underfull-CDF numerical abort.

    Saturating a CDF at the grid endpoint is part of this declared numerical
    sampler, not correction of the raw source array.  Its exact law is audited
    below against an independently constructed ideal complete instrument.
    """
    raw, residual = _raw_kernel(law)
    cdf = np.cumsum(raw, axis=0)
    result = np.zeros((462, 66), dtype=np.uint64)
    for x in range(66):
        earlier = 0
        for y in range(396):
            scaled = F(float(cdf[y, x]))*GRID
            threshold = min(GRID, -(-scaled.numerator//scaled.denominator))
            if threshold < earlier:
                raise ArithmeticError("the original binary CDF is not monotone")
            result[y, x] = threshold-earlier
            earlier = threshold
        result[396+x, x] = GRID-earlier
    if not np.all(result.sum(axis=0) == GRID):
        raise ArithmeticError("dyadic sampler does not retain unit mass")
    return result, {"binary_population_mass_residual_diagnostic": residual,
                    "grid": GRID, "numerical_abort_is_separate_terminal": True}


def audit_kernel(kernel):
    """Exact bias/variance of all retained sampler columns, not a float norm."""
    _validate_kernel(kernel)
    _, reads = population_counts()
    q = local.budget()["uniform_ideal_herald_probability"]
    sigma = local.reference()  # Validation only; infer_counts never reads this.
    raw_palm = local.incident.reference(local.CORE)[0]
    ideal = []
    for j in (-1, 0, 1):
        ideal.extend(sum((F(int(c), 2)*q*w for c, w in zip(row, sigma, strict=True)), F())
                     for row in reads[j])
    ideal.extend(F(1, 66)-q*w for w in raw_palm)
    ideal.extend([F()]*(3*66))
    if len(ideal) != 462 or sum(ideal) != 1 or min(ideal) < 0:
        raise ArithmeticError("independent exact target instrument is incomplete")
    weights = tuple(next(w for w, charge in zip(sigma, variation.charges(), strict=True) if charge == i)
                    for i in range(3))
    worst, means, seconds = F(), [F()]*3, [F()]*3
    for x in range(66):
        column = tuple(F(int(c), GRID) for c in kernel[:, x])
        worst = max(worst, sum((abs(p-r) for p, r in zip(column, ideal, strict=True)), F())/2)
        terminal = tuple(sum(column[66*k:66*(k+1)], F()) for k in range(7))
        for i, score in enumerate(SCORES):
            m = sum((s*p for s, p in zip(score, terminal, strict=True)), F())
            v = sum((s*s*p for s, p in zip(score, terminal, strict=True)), F())
            means[i] = max(means[i], abs(m-q*weights[i])); seconds[i] = max(seconds[i], v)
    if any(m > width*worst for m, width in zip(means, WIDTHS, strict=True)):
        raise ArithmeticError("exact sampler score bias exceeds complete variation")
    return {"exact_dyadic_columns_checked": 66, "complete_total_variation_upper": worst,
            "score_mean_bias_upper": tuple(means), "score_second_moment_upper": tuple(seconds),
            "entire_multitrial_path_distance_claimed": False,
            "physical_finite_instrument_bound_is_separate": True}


def _validate_kernel(kernel):
    if (type(kernel) is not np.ndarray or kernel.shape != (462, 66) or kernel.dtype != np.uint64
            or np.any(kernel > GRID) or any(sum(map(int, kernel[:, x])) != GRID for x in range(66))):
        raise ValueError("a complete original uint64 dyadic terminal/state kernel is required")


def sample_kernel(kernel, trials, *, seed=20261024, initial=0, chunk_size=1 << 16):
    """Exact dyadic common/residual sampling of one continuing carrier.

    Every common draw includes both its public record and its terminal state.
    Residual draws use that immediately preceding terminal state.  Batching
    only exploits the identical subkernel obtained by a columnwise minimum.
    """
    _validate_kernel(kernel)
    _integer(trials, 1, 1 << 30, "trial count")
    _integer(seed, 0, (1 << 63)-1, "seed")
    _integer(initial, 0, 65, "initial state")
    _integer(chunk_size, 1, 1 << 20, "chunk size")
    common_part = kernel.min(axis=1)
    common_mass = int(common_part.sum())
    common_cdf = np.cumsum(common_part)
    residual_cdf = np.cumsum(kernel-common_part[:, None], axis=0)
    rng = np.random.default_rng(seed)
    counts = np.zeros(462, dtype=np.int64)
    digest, prefix, state, exceptional = hashlib.sha256(), [], initial, 0
    for start in range(0, trials, chunk_size):
        draws = rng.integers(0, GRID, size=min(chunk_size, trials-start), dtype=np.uint64)
        selected = np.searchsorted(common_cdf, draws, side="right")
        for position in np.flatnonzero(draws >= common_mass):
            prior = state if position == 0 else int(selected[position-1]) % 66
            selected[position] = np.searchsorted(residual_cdf[:, prior],
                                                   int(draws[position])-common_mass, side="right")
            exceptional += 1
        if np.any(selected >= 462):
            raise ArithmeticError("dyadic sampling omitted terminal mass")
        state = int(selected[-1]) % 66
        counts += np.bincount(selected, minlength=462)
        digest.update(selected.astype("<u2").tobytes())
        if len(prefix) < 32:
            prefix.extend((LABELS[int(i)//66], int(i)%66) for i in selected[:32-len(prefix)])
    terminal_counts = {label: int(counts[66*i:66*(i+1)].sum()) for i, label in enumerate(LABELS)}
    return {"counts": terminal_counts, "terminal_state_counts": counts.reshape(7, 66).tolist(),
            "record_and_terminal_state_sha256": digest.hexdigest(), "prefix": prefix,
            "trials": trials, "seed": seed, "initial_state": initial, "final_state": state,
            "common_subkernel_probability": F(common_mass, GRID), "residual_draws": exceptional,
            "same_carrier_reused": True, "microscopic_attempt_tapes_executed": False,
            "seeded_PRNG_demonstration_proves_independent_uniform_randomness": False}


def execute_trials(trials=1 << 24, *, seed=20261024, law="C"):
    kernel, numerical = quantized_kernel(law)
    audit = audit_kernel(kernel)
    sample = sample_kernel(kernel, trials, seed=seed)
    certificate = _portable_sampler_certificate()
    if (audit["complete_total_variation_upper"] > certificate["complete_total_variation_upper"]
            or any(actual > bound for field in ("score_mean_bias_upper", "score_second_moment_upper")
                   for actual, bound in zip(audit[field], certificate[field], strict=True))
            or sample["common_subkernel_probability"] < 1-SAMPLER_GATE):
        raise ArithmeticError("actual sampler fails its portable source-comparison gates")
    # Export fixed certified caps, not platform-dependent exact diagnostics.
    # The latter must pass this run's audit before any portable receipt exists.
    del sample["common_subkernel_probability"]
    sample["common_subkernel_probability_lower"] = 1-SAMPLER_GATE
    inferred = infer_counts(sample["counts"], sampling_bias=SAMPLER_GATE)
    return {"law": law, "numerical_kernel": numerical, "sampler_certificate": certificate,
            "execution": sample, "inference": inferred}


def _portable_sampler_certificate():
    """Fixed rational claims independently checked against every exact column."""
    q = _record_contract()[1]
    return {"exact_dyadic_columns_checked": 66,
            "complete_total_variation_upper": SAMPLER_GATE,
            "score_mean_bias_upper": tuple(width*SAMPLER_GATE for width in WIDTHS),
            "score_second_moment_upper": tuple(
                q*vmax+SAMPLER_GATE*max(s*s for s in score)
                for score, vmax in zip(SCORES, SECTOR_SECOND_MOMENT_MAX, strict=True)),
            "entire_multitrial_path_distance_claimed": False,
            "physical_finite_instrument_bound_is_separate": True}


def run(trials=1 << 24, *, seed=20261024):
    execution = execute_trials(trials, seed=seed)
    inferred = execution["inference"]
    if inferred["status"] == "faithful":
        q = local.budget()["uniform_ideal_herald_probability"]
        weights = tuple(next(w for w, c in zip(local.reference(), variation.charges(), strict=True) if c == i)
                        for i in range(3))
        execution["validation_only_known_target_in_confidence_box"] = all(
            lo <= q*w <= hi for w, (lo, hi) in zip(weights, inferred["subnormalized_weight_intervals"], strict=True))
    paths = tuple(sorted(set((*local.SOURCE_PATHS, "model/native_modular_reference/reference_tomography.py",
                             "model/native_equilibration/common_variation.py"))))
    return {"schema": "oph.native_modular_reference.reference_tomography.v1",
            "identification": exact_identification(), "continuing_execution": execution,
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            "scope": {"fixed_known_reference_weights_used_by_estimator": False,
                      "full_state_tomography": False, "independent_copies_required": False,
                      "finite_physical_service_accuracy_established": False,
                      "strict_monotone_individual_clock_readout": False,
                      "confidence_theorem_domain": "actual native records or ideal uniform dyadic-kernel draws",
                      "fixed_seed_PRNG_run": "reproducible numerical illustration, not an independence certificate"}}


def _json(value):
    return json.loads(json.dumps(value, default=str))


def compare_receipt(supplied, actual):
    """Check exact records and source custody; revalidate binary diagnostics.

    All sampler proof fields are fixed portable rational caps, certified by
    this run's exact all-column audit. Only the raw floating mass residual may
    vary between platforms, within its declared diagnostic gate.
    """
    cleaned = []
    for original in (supplied, actual):
        value = deepcopy(original)
        try:
            row = value["continuing_execution"]
            diagnostic = row["numerical_kernel"]["binary_population_mass_residual_diagnostic"]
            if type(diagnostic) is not float or not math.isfinite(diagnostic) or not 0 <= diagnostic <= 5e-10:
                return False
            row["numerical_kernel"]["binary_population_mass_residual_diagnostic"] = "checked"
            if json.dumps(row["sampler_certificate"], sort_keys=True, allow_nan=False) != json.dumps(
                    _json(_portable_sampler_certificate()), sort_keys=True, allow_nan=False):
                return False
            if row["execution"]["common_subkernel_probability_lower"] != str(1-SAMPLER_GATE):
                return False
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            return False
        cleaned.append(value)
    try:
        return json.dumps(cleaned[0], sort_keys=True, allow_nan=False, separators=(",", ":")) == json.dumps(
            cleaned[1], sort_keys=True, allow_nan=False, separators=(",", ":"))
    except (TypeError, ValueError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    payload = _json(run())
    if args.check:
        if not compare_receipt(strict_json.load(args.check), payload):
            raise SystemExit("native reference tomography receipt drift")
    else:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print("Native complete-record reference identification verified.")


if __name__ == "__main__":
    main()
