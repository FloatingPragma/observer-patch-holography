"""Frozen, source-bound finite experiments; --record never overwrites evidence."""
import argparse
from collections import Counter
from dataclasses import replace
from fractions import Fraction
import hashlib
import gzip
from itertools import product
import json
from math import comb, log, pi, sqrt, ulp
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.linalg import expm

from tools.pin_scope import bound, bound_pins

from .clock import (Construction, clock_matrix,
                    clock_probabilities, evolve, measured_return,
                    prepared_span, readout, require, step)
from .native import (ADDRESSES, BACKGROUNDS, decode_history, loads,
                     omitted_record_overlap,
                     reused_record_counterexample)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOL = 2e-10
BINDINGS = tuple(bound("model/native_memory_clock/" + name for name in (
    "__init__.py", "native.py", "clock.py", "experiment.py", "protocol.json",
    "provenance.json")))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_event(source, address, background, law):
    """Direct incidence/flux oracle, not the reversible gate implementation."""
    q = list(background)
    for edge, (a, b) in enumerate(ADDRESSES[::2]):
        flux = ((source >> (2 * edge + 1)) & 1) - ((source >> (2 * edge)) & 1)
        q[a] += flux; q[b] -= flux
    a, b = ADDRESSES[address]
    bit = (source >> address) & 1
    sign = -1 if law == "uphill" else 1
    changes = law != "hold" and ((not bit and sign * (q[a] - q[b]) > 0)
                                or (bit and sign * (q[a] - q[b]) < 0))
    return source ^ (1 << address) if changes else source


def native_records(source, n, background, law):
    """All address histories, exact integer counts, including no-operations."""
    out = Counter()
    for addresses in product(range(12), repeat=n):
        current, old = source, []
        for k in addresses:
            old.append((current >> k) & 1)
            current = independent_event(current, k, background, law)
        out[(n, current, addresses, tuple(old))] += 1
    return {k: Fraction(v, 12 ** n) for k, v in out.items()}


def compare(expected, actual, path="result"):
    """Full typed scientific comparison, with no ignored scientific fields."""
    require(type(expected) is type(actual), path + ": type")
    if isinstance(expected, dict):
        require(expected.keys() == actual.keys(), path + ": keys")
        for key in expected:
            compare(expected[key], actual[key], path + "." + key)
    elif isinstance(expected, list):
        require(len(expected) == len(actual), path + ": length")
        field = path.rsplit(".", 1)[-1]
        if field in ("outcomes", "complete_classical_records", "endpoint_source_law"):
            check_probability_vectors([row[-1] for row in expected], [row[-1] for row in actual], path)
        elif field in ("x_record_law", "uniform_clock_same_read_time_law", "interior_measured_return_clock_law"):
            check_probability_vectors(expected, actual, path)
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, path + f"[{i}]")
    elif type(expected) is float:
        require(np.isfinite(expected) and np.isfinite(actual)
                and abs(expected - actual) <= TOL, path + ": numerical value")
    else:
        require(expected == actual, path + ": value")


def check_probability_vectors(left, right, path):
    for values in (left, right):
        require(all(type(x) is float and np.isfinite(x) and x >= 0 for x in values),
                path + ": invalid probability")
        require(abs(sum(values) - 1) <= TOL, path + ": incomplete outcome law")
    require(sum(abs(x - y) for x, y in zip(left, right)) <= TOL,
            path + ": aggregate probability discrepancy")


def law_error(actual, expected):
    return sum(abs(actual.get(k, 0) - expected.get(k, 0)) for k in actual.keys() | expected.keys())


def binary_coupling_bound(c):
    """Exact rational upper bound on the actual binary coupling error.

    |x-sqrt(r)|=|x*x-r|/(x+sqrt(r)) <= |x*x-r|/x. Row sums bound
    the full dressed Hamiltonian norm. This does not certify SciPy roundoff.
    """
    errors = []
    for j in range(c.horizon):
        x = Fraction.from_float(c.coupling(j))
        r = (Fraction(1, 4) if c.profile == "uniform" else
             Fraction((j + 1) * (c.horizon - j), c.horizon ** 2))
        errors.append(abs(x*x-r)/x)
    norm_bound = max((errors[j-1] if j else 0) + (errors[j] if j < c.horizon else 0)
                     for j in range(c.horizon + 1))
    return {"operator_norm_error_upper_rational": str(norm_bound),
            "state_error_bound_at_binary_time_tau": "abs(tau) * operator_norm_error_upper_rational",
            "complete_record_l1_error_bound": "2 * state_error_bound",
            "numerical_exponential_residual_gate": TOL,
            "numerical_roundoff_is_not_certified_by_coupling_bound": True}


def pack_law(probabilities):
    # Every computed label is retained, including numerically zero labels.
    return [[n, source, list(addresses), list(old), probability]
            for (n, source, addresses, old), probability in sorted(probabilities.items())]


def sampled_indices(law, count, seed):
    """Independent preparations; no claimed sequentially measured trajectory."""
    rng = np.random.default_rng(seed)
    probabilities = np.array([item[-1] for item in law])
    require(np.all(probabilities >= 0) and abs(probabilities.sum() - 1) < TOL,
            "incomplete sampled outcome law")
    # No clipping or renormalization. Last-bit floating CDF error is checked.
    cdf = np.cumsum(probabilities)
    draws = rng.random(count)
    indices = np.searchsorted(cdf, draws, side="right")
    require(np.all(indices < len(law)), "sampler encountered missing probability")
    return indices.tolist()


def sampling_budget(table, indices, seed):
    probabilities = np.array([row[-1] for row in table])
    k, count = len(table), len(indices)
    # Audit the actual binary CDF against exact rational cumulative sums.
    cdf = np.cumsum(probabilities)
    exact_sum, delta = Fraction(0), Fraction(0)
    for probability, cumulative in zip(probabilities, cdf):
        exact_sum += Fraction.from_float(float(probability))
        delta = max(delta, abs(Fraction.from_float(float(cumulative)) - exact_sum))
    eta = abs(exact_sum - 1)
    exact_bias = 2 * delta + eta + Fraction(1, 2 ** 52)
    bias = float(np.nextafter(float(exact_bias), np.inf))
    defect = float(abs(cdf[-1] - 1))
    missing = max(Fraction(0), 1 - Fraction.from_float(float(cdf[-1]))) + Fraction(1, 2 ** 53)
    radius = sqrt(log(2 * k * 12 / 0.01) / (2 * count)) if count else None
    error = (float(np.max(abs(np.bincount(indices, minlength=k) / count - probabilities)))
             if count else None)
    bound = radius + bias if count else None
    require(not count or error <= bound, "registered simultaneous sampling gate")
    return {"performed": bool(count), "count": count, "seed": seed,
            "hoeffding_radius": radius, "cdf_grid_bias_bound": bias,
            "cdf_endpoint_defect": defect,
            "missing_probability_upper_bound": float(np.nextafter(float(missing), np.inf)),
            "simultaneous_per_outcome_error_bound": bound,
            "largest_sample_frequency_error": error}


def validate_sampling_record(record, path="record"):
    """Recheck a retained budget against its own table, before replay comparison.

    Tiny conservative bounds do not receive the general absolute replay
    tolerance. Eight binary64 ULPs permit platform log/sqrt differences in
    the Hoeffding radius and final addition; this is a comparison allowance,
    not a certified transcendental-function error estimate.
    """
    table = record["outcomes"]
    indices = record["independent_shot_outcome_indices"]
    supplied = record["sampling"]
    require(type(table) is list and len(table) > 0, path + ": empty sampling table")
    probabilities = [row[-1] for row in table]
    check_probability_vectors(probabilities, probabilities, path + ".outcomes")
    require(type(indices) is list and all(type(index) is int and 0 <= index < len(table)
                                         for index in indices), path + ": invalid shot indices")
    require(type(supplied) is dict and type(supplied.get("seed")) is int
            and supplied["seed"] >= 0, path + ": invalid sampling seed")
    derived = sampling_budget(table, indices, supplied["seed"])
    require(supplied.keys() == derived.keys(), path + ": sampling fields")
    require(sampled_indices(table, len(indices), supplied["seed"]) == indices,
            path + ": recorded shots disagree with seed and table")
    conservative = {"cdf_grid_bias_bound", "missing_probability_upper_bound"}
    platform_arithmetic = {"hoeffding_radius", "simultaneous_per_outcome_error_bound"}
    for field, expected in derived.items():
        actual = supplied[field]
        require(type(actual) is type(expected), path + ".sampling." + field + ": type")
        if type(expected) is float:
            require(np.isfinite(actual) and actual >= 0,
                    path + ".sampling." + field + ": negative or nonfinite budget")
            if field in conservative:
                require(actual >= expected,
                        path + ".sampling." + field + ": understated conservative bound")
            else:
                allowance = 8 * max(ulp(actual), ulp(expected)) if field in platform_arithmetic else 0.0
                require(abs(actual - expected) <= allowance,
                        path + ".sampling." + field + ": inconsistent retained budget")
        else:
            require(actual == expected, path + ".sampling." + field + ": value")
    if indices:
        total = supplied["hoeffding_radius"] + supplied["cdf_grid_bias_bound"]
        allowance = 8 * max(ulp(total), ulp(supplied["simultaneous_per_outcome_error_bound"]))
        require(supplied["simultaneous_per_outcome_error_bound"] + allowance >= total,
                path + ": total sampling bound omits retained bias")
        require(supplied["largest_sample_frequency_error"]
                <= supplied["simultaneous_per_outcome_error_bound"],
                path + ": retained sampling gate fails")


def combined_sampling_failure_bound(records):
    """Upward binary enclosure of the exact sum of retained upper bounds."""
    exact = Fraction(1, 100) + sum(
        row["sampling"]["count"] * Fraction.from_float(row["sampling"]["missing_probability_upper_bound"])
        for row in records)
    return float(np.nextafter(float(exact), np.inf))


def validate_sampling(science):
    """Validate complete retained sampling metadata without ignoring fields."""
    records = []
    for case_index, case in enumerate(science["cases"]):
        for record_index, record in enumerate(case["records"]):
            validate_sampling_record(record, f"cases[{case_index}].records[{record_index}]")
            records.append(record)
    summary = science["sampling"]
    performed = any(row["sampling"]["count"] for row in records)
    expected = {"performed": performed, "ideal_iid_assumption": True,
                "fixed_seed_prng_does_not_prove_iid": True,
                "combined_frequency_and_overflow_failure_bound":
                combined_sampling_failure_bound(records) if performed else None}
    require(type(summary) is dict and summary.keys() == expected.keys(), "sampling summary fields")
    for key, value in expected.items():
        actual = summary[key]
        require(type(actual) is type(value), "sampling summary type: " + key)
        if type(value) is float:
            require(np.isfinite(actual) and actual >= value,
                    "sampling summary understates combined failure bound")
        else:
            require(actual == value, "sampling summary value: " + key)


def finite_run(c, shots):
    states, h, initial = prepared_span(c, {1: 1.0})
    t = c.horizon
    records = []
    errors = []
    energies = []
    for i, fraction in enumerate((0.25, 0.5, 1.0)):
        tau = pi * t * fraction
        psi = evolve(h, initial, tau)
        law = readout(states, psi)
        if i == 0:
            for (_, endpoint, addresses, old) in law:
                recovered, _ = decode_history(endpoint, addresses, old, c.background, c.law)
                require(recovered == 1, "terminal records fail native history decoding")
        clock = clock_probabilities(t, tau)
        expected = {}
        for n in range(t + 1):
            for key, probability in native_records(1, n, c.background, c.law).items():
                expected[key] = clock[n] * float(probability)
        error = float(law_error(law, expected))
        require(error < TOL, "complete native joint instrument disagrees with independent oracle")
        energies.append(float(np.vdot(psi, h @ psi).real))
        errors.append(error)
        table = pack_law(law)
        indices = sampled_indices(table, shots, 6100 + i)
        records.append({"time_over_pi_T": fraction, "outcomes": table,
                        "independent_shot_outcome_indices": indices,
                        "sampling": sampling_budget(table, indices, 6100 + i)})
    require(max(abs(v) for v in energies) < TOL, "free energy not conserved")
    end = {k[1]: 0.0 for k in native_records(1, t, c.background, c.law)}
    for key, value in native_records(1, t, c.background, c.law).items():
        end[key[1]] += float(value)
    result = {"horizon": t, "background": list(c.background), "law": c.law,
              "literal_reducing_dimension": len(states), "resources": c.resources(),
              "decoded_distinct_histories": len(law),
              "binary_operator_budget": binary_coupling_bound(c),
              "complete_joint_l1_errors": errors, "free_mean_energy": energies,
              "records": records, "endpoint_source_law": [[k, v] for k, v in sorted(end.items())]}
    if c.law == "native":
        measured = measured_return(c, states, h, initial, pi * t / 4, pi * t)
        check_probability_vectors(measured.tolist(), measured.tolist(), "measured return")
        predicted = comb(2 * t, t) / 4 ** t
        require(abs(measured[0] - predicted) < TOL, "clock dephasing return control")
        require(measured[0] < 0.9, "coherence-insensitive return detector")
        result["interior_measured_return_clock_law"] = measured.tolist()
        result["predicted_measured_return"] = predicted
    return result


def coherent_control():
    """Same coherent decoder for native and phase-twisted completions.

    At endpoint T=1, undo the native gate and read X on input bit 0. The
    inverse circuit and X reference are declared terminal apparatus.
    The complete classical readouts of both laws are also compared.
    """
    c = Construction(1)
    states, h, initial = prepared_span(c, {0: 1 / sqrt(2), 1: 1 / sqrt(2)})
    alternatives = []
    baseline = None
    for phase in (False, True):
        cp = replace(c, phase=phase)
        other, hp, start = prepared_span(cp, {0: 1 / sqrt(2), 1: 1 / sqrt(2)})
        require(other == states, "phase control changed the apparatus basis")
        psi = evolve(hp, start, pi / 2)
        law = readout(states, psi)
        if baseline is None:
            baseline = law
        require(law_error(law, baseline) < TOL, "phase control changed classical output law")
        # Retain unfinished clock branch, and every address/record outcome.
        decoded = {}
        for b, amp in zip(states, psi):
            dest = step(c, b, -1)[0] if b.clock == 1 else b
            # A terminal pointer separately retains the original clock label.
            key = (b.clock, dest.buffers, dest.addresses, dest.records)
            decoded[key] = decoded.get(key, 0j) + amp
        plus = minus = 0.0
        groups = {(n, buffers[0] & ~1, buffers[1:], addresses, old)
                  for n, buffers, addresses, old in decoded}
        for n, base, rest, addresses, old in groups:
            a = decoded.get((n, (base,) + rest, addresses, old), 0j)
            b = decoded.get((n, (base | 1,) + rest, addresses, old), 0j)
            plus += abs((a + b) / sqrt(2)) ** 2
            minus += abs((a - b) / sqrt(2)) ** 2
        require(abs(plus + minus - 1) < TOL, "incomplete coherent detector")
        alternatives.append({"phase_twist": phase, "x_record_law": [float(plus), float(minus)],
                             "classical_l1_from_native": float(law_error(law, baseline)),
                             "complete_classical_records": pack_law(law)})
    require(abs(alternatives[0]["x_record_law"][1]) < TOL, "native coherence lost")
    require(abs(alternatives[1]["x_record_law"][1] - 1 / 12) < TOL, "coherent alternative invisible")
    return alternatives


def compute(extended=False):
    horizons = (3, 4) if extended else (2,)
    cases = [finite_run(Construction(t, background, law), 0 if extended else 256)
             for t in horizons for background in BACKGROUNDS
             for law in ("native", "uphill")]
    separations = []
    for native, uphill in zip(cases[::2], cases[1::2]):
        a, b = dict(native["endpoint_source_law"]), dict(uphill["endpoint_source_law"])
        tv = float(law_error(a, b) / 2)
        require(tv > 0.1, "nonzero competing law fails to separate source records")
        separations.append(tv)
    spectra = []
    for t in (1, 2, 3, 4, 8, 16):
        c = Construction(t)
        values = np.linalg.eigvalsh(clock_matrix(c))
        expected = np.linspace(-1, 1, t + 1)
        require(np.max(abs(values - expected)) < TOL, "complete clock spectrum")
        spectra.append({"horizon": t, "all_distinct_eigenvalues": values.tolist(),
                        "multiplicity_each": c.resources()["full_data_dimension"]})
    uniform = Construction(3, profile="uniform")
    uv = expm(-1j * 3 * pi / 2 * clock_matrix(uniform))[:, 0]
    require(abs(uv[-1]) ** 2 < 0.99, "clock alternative not discriminated")
    no_memory = omitted_record_overlap(0, 1)
    require(no_memory == Fraction(1, 12), "old-bit omission control")
    witness = reused_record_counterexample()
    require(witness["reused"] != witness["fresh"], "dirty reuse control")
    before = sum(q * q for q in loads(1))
    after = sum(sum(q * q for q in loads(independent_event(1, k, BACKGROUNDS[0], "native")))
                for k in range(12)) / 12
    require(after < before, "strict marginal recurrence obstruction not witnessed")
    all_records = [record for case in cases for record in case["records"]]
    return {"schema": "native-memory-clock.v1", "extended": extended,
            "sampling": {"performed": not extended, "ideal_iid_assumption": True,
                         "fixed_seed_prng_does_not_prove_iid": True,
                         "combined_frequency_and_overflow_failure_bound":
                         combined_sampling_failure_bound(all_records) if not extended else None},
            "cases": cases, "competing_source_total_variation": separations,
            "complete_spectrum": spectra, "coherent_completion_controls": coherent_control(),
            "uniform_clock_same_read_time_law": (abs(uv) ** 2).tolist(),
            "omitted_record_gram": str(no_memory), "dirty_reuse": witness,
            "strict_disagreement_drop": [before, after],
            "physical_rows_validated": [], "empirical_mass_inputs": []}


def frozen():
    freeze = json.loads((HERE / "freeze.json").read_text(encoding="utf8"))
    require(freeze["schema"] == "native-memory-clock-freeze.v1", "freeze schema")
    files = bound_pins(freeze.get("files"))
    require(type(files) is dict and set(files) == set(BINDINGS),
            "freeze must bind the entire declared source and inputs")
    for name, sha in files.items():
        path = ROOT / name
        require(path.is_file() and digest(path) == sha, "frozen source changed: " + name)
    return digest(HERE / "freeze.json")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--extended", action="store_true")
    args = parser.parse_args()
    path = HERE / ("extended_results.json.gz" if args.extended else "results.json")
    require(not args.record or not path.exists(), "refusing to overwrite retained evidence")
    pin = frozen()
    start = time.perf_counter()
    result = json.loads(json.dumps(compute(args.extended), allow_nan=False))
    validate_sampling(result)
    if args.record:
        payload = {"freeze_sha256": pin, "science": result,
                   "execution": {"python": platform.python_version(), "numpy": np.__version__,
                                 "scipy": scipy.__version__, "seconds": time.perf_counter() - start}}
        raw = (json.dumps(payload, separators=(",", ":"), allow_nan=False) + "\n").encode("utf8")
        with path.open("xb") as stream:
            stream.write(gzip.compress(raw, mtime=0) if args.extended else raw)
    else:
        raw = gzip.decompress(path.read_bytes()) if args.extended else path.read_bytes()
        saved = json.loads(raw)
        require(saved["freeze_sha256"] == pin, "receipt freeze mismatch")
        validate_sampling(saved["science"])
        compare(saved["science"], result)
    print(json.dumps({"verified": True, "extended": args.extended,
                      "seconds": time.perf_counter() - start, "source_freeze": pin,
                      "physical_rows_validated": []}))


if __name__ == "__main__":
    main()
