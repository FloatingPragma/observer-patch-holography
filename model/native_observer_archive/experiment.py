"""Source-bound exhaustive observer laws, with read-only scientific replay."""
import argparse
from collections import Counter
from dataclasses import replace
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from .oracle import expected_records
from .packet import (Rail, evolve_clock, factored_vector, fraction_record,
                     interval_certificate, literal_evolve, packet, prepared_span,
                     preparation_error_bound, require, spectrum)
from .program import Program, prefixes, preparation, public_record, resources

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOL = 2e-10
BINDINGS = tuple("model/native_observer_archive/" + name for name in (
    "__init__.py", "program.py", "oracle.py", "packet.py", "experiment.py",
    "protocol.json", "provenance.json", "PROOF.md", "README.md",
    "test_program.py", "test_channel.py", "test_packet.py", "test_experiment.py",
)) + tuple("model/native_memory_clock/" + name for name in (
    "__init__.py", "native.py", "PROOF.md", "provenance.json", "freeze.json")) + (
    "model/__init__.py",
) + tuple("model/native_observer_archive/receipts/v1/" + name for name in (
    "freeze.json", "results.json.gz", "extended_results.json.gz"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def config_for(variant, background):
    config = Program(tuple(background))
    if variant == "uphill":
        return replace(config, law="uphill")
    if variant == "no_communication":
        return replace(config, communication=False)
    if variant == "no_feedback":
        return replace(config, feedback=False)
    if variant == "no_publish":
        return replace(config, write_a=False)
    require(variant == "native", "undeclared program variant")
    return config


def tv(left, right):
    return sum(abs(Fraction(left[k] - right[k], 12)) for k in left.keys() | right.keys()) / 2


def complete_case(config, source):
    """Every branch at every prefix. Each address has exact weight 1/12."""
    records = []
    for address in range(12):
        actual = tuple(public_record(data) for data in prefixes(config, preparation(source, address)))
        require(actual == expected_records(config, source, address), "complete native oracle disagreement")
        records.append(actual)
    return records


def public_law(records, stage_weights):
    """Marginalize actual registers, retaining invalid and disagreeing pointers."""
    law = np.zeros((4, 4))
    for branch in records:
        for stage, weight in enumerate(stage_weights):
            a_valid, a_value, b_valid, b_value = branch[stage][:4]
            law[2 * a_valid + a_value, 2 * b_valid + b_value] += weight / 12
    return law.reshape(-1).tolist()


def clock_record(rail, tau, direction):
    wave = evolve_clock(rail, tau, direction)
    raw = abs(wave) ** 2
    norm = float(raw.sum())
    require(abs(norm - 1) <= TOL, "packet norm residual")
    # Store a positive normalized numerical law and its unnormalized norm.
    probabilities = raw / norm
    stages = [float(probabilities[[rail.stage(n) == j for n in range(rail.sites)]].sum())
              for j in range(5)]
    energy = float(np.vdot(wave[:-1], wave[1:]).real)
    require(abs(energy) <= TOL, "declared packet energy residual")
    return {"time": tau, "direction": direction,
            "probabilities": probabilities.tolist(), "raw_norm": norm,
            "stage_probabilities": stages, "mean_hamiltonian": energy}


def literal_check(spec, source_amplitudes, config):
    rail = Rail(spec["sites"], spec["offset"], spec["order"])
    states, h, initial = prepared_span(rail, config, source_amplitudes)
    final = literal_evolve(h, initial, spec["time"])
    expected = factored_vector(rail, config, source_amplitudes, states, spec["time"])
    discrepancy = float(np.linalg.norm(final - expected))
    require(discrepancy <= TOL, "literal full-column propagation differs from gauge prediction")
    initial_energy = float(np.vdot(initial, h @ initial).real)
    final_energy = float(np.vdot(final, h @ final).real)
    pointer_energy = 0.0
    for a in range(4):
        for b in range(4):
            branch = final * np.array([data.pointer_a == a and data.pointer_b == b for _, data in states])
            pointer_energy += float(np.vdot(branch, h @ branch).real)
    require(abs(final_energy - initial_energy) <= TOL and abs(pointer_energy - final_energy) <= TOL,
            "declared preparation's energy/meter residual")
    return {"rail": {"sites": rail.sites, "offset": rail.offset, "order": rail.order},
            "time": spec["time"], "source_amplitudes": [[s, [float(complex(a).real), float(complex(a).imag)]]
                                                        for s, a in source_amplitudes.items()],
            "basis_dimension": len(states), "hamiltonian_nonzeros": h.nnz,
            "vector_l2_error": discrepancy, "norm_squared": float(np.vdot(final, final).real),
            "initial_energy": initial_energy, "final_energy": final_energy,
            "pointer_measurement_energy": pointer_energy}


def compute(extended=False):
    protocol = json.loads((HERE / "protocol.json").read_text(encoding="utf8"))
    rail = Rail(**protocol["rail"])
    certificate = interval_certificate(rail, *protocol["window"])
    require(certificate["failure_below_one_in_thousand"], "uniform interval gate")
    preparation_error = preparation_error_bound(rail)
    ideal_bound = Fraction(int(certificate["failure_upper"]["numerator"]),
                           int(certificate["failure_upper"]["denominator"]))
    binary_bound = ideal_bound + 2 * preparation_error
    require(binary_bound < Fraction(1, 1000), "normalized binary packet preparation budget")
    certificate["normalized_binary_preparation_extra_failure_upper"] = fraction_record(2 * preparation_error)
    certificate["normalized_binary_packet_failure_upper"] = fraction_record(binary_bound)
    clocks = [clock_record(rail, tau, direction) for direction in protocol["packet_directions"]
              for tau in protocol["read_times"]]
    lookup = {(row["direction"], row["time"]): row for row in clocks}
    require(lookup[(-1, 256)]["stage_probabilities"][4] <= .01, "reversed packet must challenge readiness")
    require(lookup[(1, 5120)]["stage_probabilities"][4] <= .5, "finite memory loss control")
    cases, completed_sources = [], {}
    for background in protocol["backgrounds"]:
        for source in protocol["sources"]:
            for variant in protocol["variants"]:
                config = config_for(variant, background)
                records = complete_case(config, source)
                laws = []
                for clock in clocks:
                    probabilities = public_law(records, clock["stage_probabilities"])
                    ready = float(sum(probabilities[4 * a + b] for a in (2, 3) for b in (2, 3)))
                    disagreement = float(probabilities[4 * 2 + 3] + probabilities[4 * 3 + 2])
                    if variant in ("no_communication", "no_publish"):
                        require(ready == 0, "invalid blank pointers counted as agreement")
                    else:
                        require(abs(ready - clock["stage_probabilities"][4]) <= TOL
                                and disagreement == 0, "actual autonomous observer record mismatch")
                    laws.append({"time": clock["time"], "direction": clock["direction"],
                                 "probabilities": probabilities, "both_valid": ready,
                                 "valid_disagreement": disagreement})
                key = (tuple(background), source, variant)
                completed_sources[key] = Counter(branch[4][5][4] for branch in records)
                cases.append({"background": background, "source": source, "variant": variant,
                              "address_weight": {"numerator": "1", "denominator": "12"},
                              "complete_prefix_records": records, "public_laws": laws})
    native = completed_sources[((0, 0, 0, 0), 85, "native")]
    law_tv = tv(native, completed_sources[((0, 0, 0, 0), 85, "uphill")])
    feedback_tv = tv(native, completed_sources[((0, 0, 0, 0), 85, "no_feedback")])
    require(law_tv >= Fraction(1, 2) and feedback_tv >= Fraction(1, 4), "nonvacuous same-apparatus controls")
    spec = protocol["extended_literal_check"] if extended else protocol["literal_check"]
    source_amplitudes = {spec["source"]: 1.0} if extended else {1: 1 / np.sqrt(2), 85: 1j / np.sqrt(2)}
    literal = literal_check(spec, source_amplitudes, Program())
    return {"construction": protocol["construction"], "extended": extended,
            "protocol_sha256": digest(HERE / "protocol.json"), "interval_certificate": certificate,
            "binary_packet_vector_error_upper": fraction_record(preparation_error),
            "resources": {**resources(), "clock_sites": rail.sites, "clock_sector": "exactly one excitation",
                          "clock_edges": rail.sites - 1, "edge_coupling": "1/2",
                          "preparation": "supplied compact coherent packet, source, address and blank memories",
                          "full_carrier_dimension": str(rail.sites * int(resources()["full_data_dimension"]))},
            "spectrum": {"values": spectrum(rail).tolist(),
                         "each_multiplicity": resources()["full_data_dimension"],
                         "definition": "cos(k*pi/(N+1)), k=1,...,N; all modes retained"},
            "clock_laws": clocks, "cases": cases,
            "controls": {"completed_source_tv_native_uphill": fraction_record(law_tv),
                         "completed_source_tv_native_feedback_off": fraction_record(feedback_tv)},
            "literal_operator_check": literal,
            "physical_rows_validated": [], "empirical_calibration_inputs": []}


def compare(expected, actual, path="science"):
    require(type(expected) is type(actual), path + ": type")
    if isinstance(expected, dict):
        require(expected.keys() == actual.keys(), path + ": keys")
        for key in expected:
            compare(expected[key], actual[key], path + "." + key)
    elif isinstance(expected, list):
        require(len(expected) == len(actual), path + ": length")
        if path.endswith("probabilities"):
            for values in (expected, actual):
                require(all(type(x) is float and np.isfinite(x) and x >= 0 for x in values), path + ": positivity")
                require(abs(sum(values) - 1) <= TOL, path + ": normalization")
            require(sum(abs(x - y) for x, y in zip(expected, actual)) <= TOL, path + ": aggregate error")
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, path + f"[{i}]")
    elif type(expected) is float:
        field = path.rsplit(".", 1)[-1]
        if field in ("both_valid", "valid_disagreement"):
            require(all(0 <= x <= 1 + 8 * np.spacing(1.0) for x in (expected, actual)),
                    path + ": probability domain")
        elif field == "vector_l2_error":
            require(expected >= 0 and actual >= 0, path + ": nonnegative error")
        require(np.isfinite(expected) and np.isfinite(actual) and abs(expected - actual) <= TOL,
                path + ": numerical value")
    else:
        require(expected == actual, path + ": value")


def validate_science(science):
    """Check each result's frozen gates before numerical replay comparison.

    A tolerance between two executions does not admit a stored result that
    fails its own gate. Derived pointer summaries also describe the stored
    probability law itself, with a separate small summation allowance.
    """
    protocol = json.loads((HERE / "protocol.json").read_text(encoding="utf8"))
    rail = Rail(**protocol["rail"])
    arithmetic = 8 * np.spacing(1.0)

    def number(value, label):
        require(type(value) is float and np.isfinite(value), label + ": finite float required")
        return value

    def probabilities(values, size, label):
        require(type(values) is list and len(values) == size, label + ": complete outcome menu")
        require(all(number(value, label) >= 0 for value in values), label + ": positivity")
        require(abs(sum(values) - 1) <= TOL, label + ": normalization")

    def rational(value, label):
        require(type(value) is dict and value.keys() == {"numerator", "denominator"}
                and all(type(item) is str for item in value.values()), label + ": rational fields")
        out = Fraction(int(value["numerator"]), int(value["denominator"]))
        require(str(out.numerator) == value["numerator"] and str(out.denominator) == value["denominator"],
                label + ": canonical rational required")
        return out

    require(type(science) is dict and type(science["extended"]) is bool, "scientific result type")
    require(science["construction"] == protocol["construction"]
            and science["protocol_sha256"] == digest(HERE / "protocol.json"), "scientific protocol binding")
    clocks = {}
    for row in science["clock_laws"]:
        require(type(row["direction"]) is int and type(row["time"]) is int, "clock coordinates")
        key = (row["direction"], row["time"])
        require(key not in clocks, "duplicate clock law")
        clocks[key] = row
        probabilities(row["probabilities"], rail.sites, "clock probabilities")
        probabilities(row["stage_probabilities"], 5, "clock stage probabilities")
        require(abs(number(row["raw_norm"], "raw_norm") - 1) <= TOL, "stored packet norm residual")
        require(abs(number(row["mean_hamiltonian"], "mean_hamiltonian")) <= TOL,
                "stored packet energy residual")
        for stage in range(5):
            derived = float(np.array(row["probabilities"])[[rail.stage(n) == stage for n in range(rail.sites)]].sum())
            require(abs(derived - row["stage_probabilities"][stage]) <= arithmetic,
                    "stored clock stages disagree with their probability law")
    require(set(clocks) == {(direction, tau) for direction in protocol["packet_directions"]
                           for tau in protocol["read_times"]}, "complete clock settings")
    require(clocks[-1, 256]["stage_probabilities"][4] <= .01, "stored reversed packet readiness gate")
    require(clocks[1, 5120]["stage_probabilities"][4] <= .5, "stored finite memory loss gate")

    seen = set()
    for case in science["cases"]:
        key = (tuple(case["background"]), case["source"], case["variant"])
        require(key not in seen, "duplicate source case")
        seen.add(key)
        require(len(case["complete_prefix_records"]) == 12
                and all(len(branch) == 5 for branch in case["complete_prefix_records"]),
                "complete address and prefix menu")
        require(rational(case["address_weight"], "address weight") == Fraction(1, 12), "uniform address weight")
        records = case["complete_prefix_records"]
        support = set()
        for branch in records:
            for record in branch:
                require(all(type(bit) is int and bit in (0, 1) for bit in record[:4]),
                        "stored pointer record bit domain")
                support.add(4 * (2 * record[0] + record[1]) + 2 * record[2] + record[3])
        settings = set()
        for row in case["public_laws"]:
            setting = (row["direction"], row["time"])
            require(setting in clocks and setting not in settings, "public clock setting")
            settings.add(setting)
            law = row["probabilities"]
            probabilities(law, 16, "public probabilities")
            ready = sum(law[4 * a + b] for a in (2, 3) for b in (2, 3))
            disagreement = law[11] + law[14]
            stored_ready = number(row["both_valid"], "both_valid")
            stored_disagreement = number(row["valid_disagreement"], "valid_disagreement")
            require(0 <= stored_ready <= 1 + arithmetic and abs(stored_ready - ready) <= arithmetic,
                    "stored readiness disagrees with its pointer law")
            require(disagreement == 0 and stored_disagreement == 0, "stored valid pointer disagreement")
            if case["variant"] in ("no_communication", "no_publish"):
                require(ready == 0 and stored_ready == 0, "stored invalid blank pointers counted as agreement")
                require(sum(law[4 * a + b] for a in range(4) for b in (2, 3)) == 0,
                        "stored disabled receiver has a valid pointer")
                if case["variant"] == "no_publish":
                    require(sum(law[4 * a + b] for a in (2, 3) for b in range(4)) == 0,
                            "stored disabled publisher has a valid pointer")
            else:
                require(abs(ready - clocks[setting]["stage_probabilities"][4]) <= TOL,
                        "stored readiness differs from actual program completion")
            require(all(value == 0 for index, value in enumerate(law) if index not in support),
                    "stored public law has mass outside its complete prefix records")
            derived_law = public_law(records, clocks[setting]["stage_probabilities"])
            require(sum(abs(left - right) for left, right in zip(law, derived_law)) <= arithmetic,
                    "stored public law disagrees with its clock and record factors")
        require(settings == set(clocks), "complete public clock settings")
    require(seen == {(tuple(background), source, variant) for background in protocol["backgrounds"]
                     for source in protocol["sources"] for variant in protocol["variants"]}, "complete source cases")

    literal = science["literal_operator_check"]
    require(0 <= number(literal["vector_l2_error"], "vector_l2_error") <= TOL,
            "stored literal propagation residual gate")
    require(abs(number(literal["norm_squared"], "norm_squared") - 1) <= TOL,
            "stored literal norm residual gate")
    initial = number(literal["initial_energy"], "initial_energy")
    final = number(literal["final_energy"], "final_energy")
    measured = number(literal["pointer_measurement_energy"], "pointer_measurement_energy")
    require(abs(initial) <= TOL and abs(final - initial) <= TOL and abs(measured - final) <= TOL,
            "stored literal energy/meter residual gate")

    certificate = science["interval_certificate"]
    require(certificate["window"] == protocol["window"], "stored analytic interval")
    cantelli = rational(certificate["cantelli_failure"], "Cantelli bound")
    boundary = rational(certificate["boundary_upper"], "boundary bound")
    ideal = rational(certificate["failure_upper"], "ideal failure bound")
    preparation_error = rational(science["binary_packet_vector_error_upper"], "preparation error")
    extra = rational(certificate["normalized_binary_preparation_extra_failure_upper"], "preparation allowance")
    binary = rational(certificate["normalized_binary_packet_failure_upper"], "binary failure bound")
    require(0 <= cantelli and 0 < boundary and ideal == cantelli + boundary
            and 0 <= preparation_error and extra == 2 * preparation_error
            and binary == ideal + extra and binary < Fraction(1, 1000)
            and certificate["failure_below_one_in_thousand"] is True, "stored analytic interval budget gate")
    controls = science["controls"]
    require(Fraction(1, 2) <= rational(controls["completed_source_tv_native_uphill"], "law control") <= 1
            and Fraction(1, 4) <= rational(controls["completed_source_tv_native_feedback_off"], "feedback control") <= 1,
            "stored same-apparatus control gates")
    require(science["physical_rows_validated"] == [] and science["empirical_calibration_inputs"] == [],
            "undeclared physical award or calibration input")


def frozen():
    freeze = json.loads((HERE / "freeze.json").read_text(encoding="utf8"))
    require(freeze.keys() == {"schema", "files"}
            and freeze["schema"] == "native-observer-archive-freeze.v2", "freeze schema")
    require(type(freeze["files"]) is dict and set(freeze["files"]) == set(BINDINGS), "complete source binding")
    for name, sha in freeze["files"].items():
        require(type(sha) is str and digest(ROOT / name) == sha, "frozen source changed: " + name)
    return digest(HERE / "freeze.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--extended", action="store_true")
    args = parser.parse_args()
    target = HERE / ("extended_results.json.gz" if args.extended else "results.json.gz")
    require(not args.record or not target.exists(), "refusing to overwrite retained evidence")
    pin = frozen()
    started = time.perf_counter()
    result = json.loads(json.dumps(compute(args.extended), allow_nan=False))
    validate_science(result)
    if args.record:
        payload = {"freeze_sha256": pin, "science": result,
                   "execution": {"python": platform.python_version(), "numpy": np.__version__,
                                 "scipy": scipy.__version__, "seconds": time.perf_counter() - started}}
        raw = (json.dumps(payload, separators=(",", ":"), allow_nan=False) + "\n").encode("utf8")
        with target.open("xb") as stream:
            stream.write(gzip.compress(raw, mtime=0))
    else:
        saved = json.loads(gzip.decompress(target.read_bytes()))
        require(saved.keys() == {"freeze_sha256", "science", "execution"}
                and saved["freeze_sha256"] == pin, "receipt schema or source pin")
        validate_science(saved["science"])
        compare(saved["science"], result)
    print(json.dumps({"verified": True, "extended": args.extended,
                      "seconds": time.perf_counter() - started, "source_freeze": pin,
                      "physical_rows_validated": []}))


if __name__ == "__main__":
    main()
