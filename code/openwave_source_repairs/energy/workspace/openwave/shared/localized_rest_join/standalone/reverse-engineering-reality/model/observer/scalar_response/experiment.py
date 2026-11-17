"""Execute a conditional scalar instrument between two native repair events.

The scalar action is supplied apparatus. Exact rational bounds belong to the
mathematical source; floating spectral execution is a separately checked
diagnostic, not a new interval or physical-realization certificate.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import random

from model.observer import native

ROOT = Path(__file__).resolve().parents[3]
PARENT = "code/source_scalar_packet/source_common_scalar_receipt.json"
PINS = {
    PARENT: "c050889ea4b8e46b4635ed6ab7b850f2b8af24f3cfcac51e4a2b349d865ba7de",
    "model/observer/native.py": "dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d",
}
GRIDS = (13, 34, 89, 233)
TIMES = (Q(19, 20), Q(39, 40), Q(1))
TRIALS = 200_000
SEED = 20261012


def sources():
    for name, expected in PINS.items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source pin mismatch: {name}")
    parent = json.loads((ROOT / PARENT).read_text())
    for name, pin in parent["source_pins"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != pin["sha256"]:
            raise ValueError(f"Scalar parent source pin mismatch: {name}")
    return parent


def detector_probability(response, current, outcome):
    from .field import _current, _outcome, _real
    response = _real(response, "scalar response")
    current, outcome = _current(current), _outcome(outcome)
    if abs(response) > 0.5:
        raise ValueError("Invalid scalar response or detector symbol")
    return 0.5 + outcome * current * response


def decoder_seam(detector_bit):
    """Only the committed field-detector bit is accessible to the actuator."""
    from .field import _outcome
    detector_bit = _outcome(detector_bit)
    return native.F if detector_bit == 1 else native.E


@dataclass(frozen=True)
class Record:
    first: tuple[int, int]
    detector: int
    second: tuple[int, int] | None = None


def feedback(loads, record, coin):
    if record.second is not None:
        return loads, record  # full tape halts
    result, second = native.native(loads, decoder_seam(record.detector), coin)
    return result, Record(record.first, record.detector, second)


def affine_histories():
    """All 528 paths, with exact weight a+b*r for the field response r."""
    rows = []
    for initial in native.states(binary=True):
        for c1 in (0, 1):
            first, r1 = native.native(initial, native.E, c1)
            for bit in (-1, 1):
                for c2 in (0, 1):
                    final, record = feedback(first, Record(r1, bit), c2)
                    rows.append((Q(1, 528), Q(bit*r1[1], 264),
                                 initial, first, final, record))
    return rows


def exact_statistics():
    rows = affine_histories()
    functions = {
        "total_probability": lambda row: 1,
        "second_current": lambda row: row[-1].second[1],
        "final_port0": lambda row: row[-2][0],
        "detector_current_correlation": lambda row: row[-1].detector * row[-1].first[1],
        "detector_plus": lambda row: int(row[-1].detector == 1),
    }
    coefficients = {
        name: [sum((row[k]*fn(row) for row in rows), Q()) for k in (0, 1)]
        for name, fn in functions.items()
    }
    expected = {"total_probability": [Q(1), Q(0)],
                "second_current": [Q(0), Q(5, 132)],
                "final_port0": [Q(1, 6), Q(-5, 132)],
                "detector_current_correlation": [Q(0), Q(10, 33)],
                "detector_plus": [Q(1, 2), Q(0)]}
    if coefficients != expected:
        raise ValueError("Native/scalar coefficient mismatch")
    tv_coefficient = sum((abs(row[1]) for row in rows), Q()) / 2
    if tv_coefficient != Q(5, 33):
        raise ValueError("History variation coefficient mismatch")
    return {
        "path_count": len(rows),
        "reachable_record_words": len({row[-1] for row in rows}),
        "affine_coefficients_in_response": {k: list(map(str, v)) for k, v in coefficients.items()},
        "history_total_variation_per_response_difference": str(tv_coefficient),
    }


def transferred_bounds(parent):
    """Exact postprocessing of the pinned *mathematical* source certificate."""
    rows = []
    for level in parent["levels"]:
        epsilon = Q(level["graph_vs_compact_continuum_error_upper"][1])
        response = list(map(Q, level["graph_compact_induced_response"]))
        rows.append({
            "q": level["q"],
            "field_probability_error_upper": str(epsilon),
            "history_total_variation_upper": str(Q(5, 33)*epsilon),
            "native_current_error_upper": str(Q(5, 132)*epsilon),
            "native_current_interval": [str(Q(5, 132)*x) for x in response],
            "native_current_magnitude_lower": str(Q(5, 132)*Q(level["resolved_graph_signal_lower"])),
        })
    return rows


def mean_summary(total, squared, n, expected):
    mean = total / n
    variance = max(0.0, (squared - total*total/n) / (n-1)) if n > 1 else 0.0
    return {"sum": total, "sum_squares": squared, "mean": mean,
            "standard_error": math.sqrt(variance/n), "forecast": expected}


def sample(response, trials=TRIALS, seed=SEED):
    """Fresh field preparations, one Born outcome and two native events each.

    Field evolution and its exact Gaussian instrument are evaluated before
    sampling. The postmeasurement field is terminal; it is never silently
    reset and reused as an undisturbed field for a later measurement.
    """
    if type(trials) is not int or trials < 2 or type(seed) is not int:
        raise ValueError("At least two trials and an integer seed are required")
    detector_probability(response, 1, 1)
    rng = random.Random(seed)
    configurations = native.states(binary=True)
    totals, squares, words = Counter(), Counter(), Counter()
    traces = []
    forecast = 5*response/132  # fixed before outcome sampling
    for trial in range(trials):
        initial = rng.choice(configurations)
        first, r1 = native.native(initial, native.E, rng.randrange(2))
        u = rng.random()
        bit = 1 if u < detector_probability(response, r1[1], 1) else -1
        reverse_bit = 1 if u < detector_probability(-response, r1[1], 1) else -1
        # A donor detector has exactly the fair marginal on the declared prior.
        donor_bit = 1 if rng.random() < 0.5 else -1
        c2 = rng.randrange(2)
        final, record = feedback(first, Record(r1, bit), c2)
        _, reverse = feedback(first, Record(r1, reverse_bit), c2)
        _, donor = feedback(first, Record(r1, donor_bit), c2)
        _, disconnected = native.native(first, native.E, c2)
        if record.first != r1 or record.detector != bit or feedback(final, record, c2) != (final, record):
            raise ValueError("Committed record or halt failure")
        values = {"current": record.second[1], "reversed": reverse.second[1],
                  "independent_donor": donor.second[1], "disconnected": disconnected[1],
                  "paired_current_minus_donor": record.second[1]-donor.second[1],
                  "final_port0": final[0], "detector_plus": int(bit == 1),
                  "detector_current_correlation": bit*r1[1]}
        for name, value in values.items():
            totals[name] += value
            squares[name] += value*value
        words[(r1, bit, record.second)] += 1
        if trial < 8:
            traces.append({"first_record": list(r1), "detector_bit": bit,
                           "selected_seam": list(decoder_seam(bit)),
                           "second_record": list(record.second), "halted": True,
                           "field_branch": {"preparation_current": r1[1], "kraus_outcome": bit}})
    expected = {"current": forecast, "reversed": -forecast,
                "independent_donor": 0.0, "disconnected": 0.0,
                "paired_current_minus_donor": forecast,
                "final_port0": 1/6-forecast, "detector_plus": 0.5,
                "detector_current_correlation": 10*response/33}
    return {"trials": trials, "seed": seed,
            "statistics": {k: mean_summary(totals[k], squares[k], trials, v) for k, v in expected.items()},
            "record_counts": [{"first": list(a), "detector": b, "second": list(c), "count": n}
                              for (a, b, c), n in sorted(words.items())],
            "first_eight_trials": traces}


def run(trials=TRIALS, seed=SEED):
    from .field import ScalarField
    parent = sources()
    exact = exact_statistics()
    bounds = transferred_bounds(parent)
    executions = []
    checks = {"all_528_native_paths_normalized_affine_in_field_response": True,
              "exact_history_error_transfer": exact["history_total_variation_per_response_difference"] == "5/33"}
    for q, level in zip(GRIDS, parent["levels"]):
        field = ScalarField(q)
        responses = []
        for time in TIMES:
            result = field.response(float(time))
            # The field implementation computes mode evolution, not this receipt.
            p = result["plus_probability"]
            r = p-0.5
            lo, hi = map(lambda x: float(Q(x)), level["graph_compact_induced_response"])
            checks[f"q{q}_t{time}_inside_parent_response_interval"] = lo-1e-10 <= r <= hi+1e-10
            responses.append({"model_time_fraction": str(time), **result, "response": r,
                              "predicted_native_current": 5*r/132})
        primary = responses[1]
        branches = [field.branch_descriptor(float(TIMES[1]), y, current=j)
                    for j in (-1, 0, 1) for y in (-1, 1)]
        for branch in branches:
            c = [complex(*x) for x in branch["coefficients"]]
            norm = sum(c[a].conjugate()*complex(*branch["gram"][a][b])*c[b]
                       for a in (0, 1) for b in (0, 1))
            checks[f"q{q}_j{branch['current']}_y{branch['outcome']}_conditional_field_norm"] = abs(norm-1) < 2e-12
        sampled = sample(primary["response"], trials, seed)
        for name, item in sampled["statistics"].items():
            checks[f"q{q}_{name}_sample_within_six_standard_errors"] = (
                abs(item["mean"]-item["forecast"]) <= 6*item["standard_error"]+1/trials)
        executions.append({"q": q, "responses": responses,
                           "terminal_field_branches": branches, "sample": sampled})
    result = {
        "schema": "oph.observer.scalar_response.v1", "source_pins": PINS,
        "scope": {
            "construction": "native current -> supplied reversible scalar apparatus -> binary record -> native feedback",
            "time": "declared scalar model coordinate c*t/L, not native count or laboratory time",
            "field": "full Gaussian mode evolution and one centered Weyl instrument; no Fock cutoff",
            "postmeasurement_field": "two-coherent-component branch K_y|alpha>; terminal after read, not reused",
            "numerics": "floating diagnostics; no certified roundoff bound or executed native scalar pulse word",
            "bounds": "exact inherited ideal-model continuum bounds transferred to complete classical histories",
            "supplied": "golden geometry, scalar action, vacuum, coherent preparation, Born readout, ideal controls and model clock",
            "physical_prediction": "conditional apparatus prediction, not empirical validation or an OPH-specific new constant",
        },
        "exact_native_statistics": exact, "refinement_bounds": bounds,
        "executions": executions, "checks": checks, "all_pass": all(checks.values()),
    }
    if not result["all_pass"]:
        raise ValueError("Scalar observer check failed: " + ", ".join(k for k, v in checks.items() if not v))
    return result


def compare_receipt(expected, actual, path="root"):
    """Exact structure/counts/bounds, tolerance only for floating diagnostics."""
    if type(expected) is not type(actual):
        raise ValueError(f"Receipt type mismatch at {path}")
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError(f"Receipt keys mismatch at {path}")
        for key in expected:
            compare_receipt(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f"Receipt length mismatch at {path}")
        for i, (a, b) in enumerate(zip(expected, actual)):
            compare_receipt(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        if not math.isfinite(expected) or not math.isfinite(actual) or abs(expected-actual) > 2e-10:
            raise ValueError(f"Floating replay mismatch at {path}")
    elif expected != actual:
        raise ValueError(f"Receipt mismatch at {path}")


def canonical(result):
    return json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n"
