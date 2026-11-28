"""Record each donor-release clock disposition before attempting a flight."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import log, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.executable_model import federation
from tools.strict_json import load, loads

from . import transfer_clock as clock
from .donor_release import (
    BENCHMARK,
    NUMERICAL_ALLOWANCE,
    SHOTS,
    commissioning_failure,
    port_records,
    witness_floor,
)
from .experiment import BOARD_SHA, CARRIER_DIR, HERE, TIMES, compare
from .graph import DONOR, History, donor_selection
from .local import check_local, require
from .phase_reader import effects
from .source import install, triangles

STRENGTHS = tuple(
    sign * magnitude for magnitude in (0.5, 1, 2, 4, 8, 16, 64, 128) for sign in (1, -1)
) + (0,)


def clock_outcomes(receipt):
    """Return one disposition per scheduled bank, including every refusal."""
    if (
        not isinstance(receipt, dict)
        or receipt.get("shots_per_appointment") != clock.SHOTS
        or receipt.get("appointments") != clock.T.tolist()
        or not isinstance(receipt.get("runs"), list)
    ):
        raise ValueError("Clock receipt has an invalid schedule or trial count")
    runs = receipt["runs"]
    if len(runs) > len(STRENGTHS):
        raise ValueError("Clock receipt contains unexpected extra settings")
    for index, strength in enumerate(STRENGTHS):
        recorded = runs[index] if index < len(runs) else {}
        row = {"strength_in_simulator": strength, "clock_calibration": None}
        if not isinstance(recorded, dict):
            yield row | {"disposition": "invalid_clock_bank"}
            continue
        if "strength_in_simulator" in recorded and (
            type(recorded["strength_in_simulator"]) not in (int, float)
            or recorded["strength_in_simulator"] != strength
        ):
            yield row | {"disposition": "invalid_clock_bank", "clock_record": recorded}
            continue
        if "counts" not in recorded:
            yield row | {"disposition": "missing_clock_bank"}
            continue
        counts = recorded["counts"]
        row["clock_counts"] = counts
        if not isinstance(counts, list):
            yield row | {"disposition": "invalid_clock_bank"}
            continue
        try:
            reading = clock.infer(counts)
        except ValueError:
            yield row | {"disposition": "invalid_clock_bank"}
            continue
        if reading is None:
            yield row | {"disposition": "inconclusive_clock"}
            continue
        row["clock_calibration"] = reading
        _, lower, upper, _ = reading
        # Simulator consistency check. The clock reader above only takes counts.
        if not lower <= abs(strength) <= upper:
            yield row | {"disposition": "clock_rejected"}
            continue
        yield row | {"disposition": "calibrated"}


def run(receipt):
    decisions = list(clock_outcomes(receipt))
    check_local()
    board = federation(0)
    fixture = CARRIER_DIR / "carriers" / "FEDERATION_LEVEL0.npz"
    require(sha256(fixture.read_bytes()).hexdigest() == BOARD_SHA[0], "Board differs")
    heads, operations = install(board, triangles(board)[0])
    tree = History(board, heads)
    choices = donor_selection(tree)["choices"]
    readers = {
        choice["analyzer"]: effects(tree, choice["paths"], -1) for choice in choices
    }
    analyzer = load(DONOR / "analyzers.json")
    addresses = np.array([tree.sites == port for port in analyzer["ports"]], dtype=int)
    angles = np.array([word["angles"] for word in analyzer["pulses"]]) @ addresses
    operations_count = sum(name in ("write", "head_hop") for name, _, _ in operations)
    require(operations_count == len(board.seams) + 1, "Commissioning count differs")
    failure = commissioning_failure(operations_count, 96)
    energies, basis = eigh(tree.twice_scaled)
    rates, directions = eigh(tree.lap)
    ports = np.unique(tree.read_ports).tolist()
    coverage = sqrt(log(2 * len(choices) * len(STRENGTHS) / 0.01) / (2 * SHOTS))

    for index, row in enumerate(decisions):
        if row["disposition"] != "calibrated":
            continue
        strength = row["strength_in_simulator"]
        states = clock.run_word(
            tree, angles, analyzer["step"], strength / row["clock_calibration"][0]
        )
        evolutions = {
            time: (basis * np.exp(-0.5j * strength * time * energies)) @ basis.T
            for time in TIMES
        }
        stochastic = {
            time: (directions * np.exp(-0.5 * abs(strength) * time * rates))
            @ directions.T
            for time in TIMES
        }
        row["launches"] = []
        row["disposition"] = "passed"
        for choice in choices:
            column = choice["analyzer"]
            prepared = states[:, column]
            effect = readers[column]
            phase_law = np.einsum("i,kij,j->k", prepared.conj(), effect, prepared).real
            require(min(phase_law) >= -NUMERICAL_ALLOWANCE, "Invalid witness law")
            outcomes = np.r_[(1 - failure) * phase_law, failure]
            require(
                abs(outcomes.sum() - 1) <= NUMERICAL_ALLOWANCE, "Incomplete witness"
            )
            rounded = outcomes.round(12)
            require(min(rounded) >= 0, "Witness sampling grid contains a negative bin")
            counts = np.random.default_rng([441304, index, column]).multinomial(
                SHOTS, rounded / rounded.sum()
            )
            computed = witness_floor(float(outcomes[-2]))
            sampled = witness_floor(float(counts[-2] / SHOTS), coverage)
            violations = []
            if min(computed, sampled) < BENCHMARK:
                violations.append("launch_witness")
            records = port_records(tree, prepared, evolutions, stochastic, failure)
            tips = [ports.index(path[-1]) for path in choice["paths"]]
            for record in records:
                home = {
                    arm: sum(probabilities[tip] for tip in tips)
                    for arm, probabilities in record["port_probabilities"].items()
                }
                if home["coherent"] < computed - NUMERICAL_ALLOWANCE:
                    violations.append("protected_weight")
                if (
                    record["exchange_exposure"] == 37
                    and max(home["dephased"], home["rate_matched_swaps"]) >= BENCHMARK
                ):
                    violations.append("changed_law_control")
                record["two_port_probabilities"] = home
            row["launches"].append(
                {
                    "column": column,
                    "launch_witness_probabilities": rounded.tolist(),
                    "simulated_witness_counts": counts.tolist(),
                    "computed_all_time_floor": computed,
                    "sampling_all_time_floor": sampled,
                    "terminal": records,
                    "violations": violations,
                }
            )
            if violations:
                row["disposition"] = "scientific_rejection"
    return {
        "schema": "oph.interval_echo.row_04.donor_release_checked.v1",
        "physical_port_order": ports,
        "simultaneous_99_half_width": coverage,
        "fair_commissioning_operations": operations_count,
        "commissioning_exhaustion": failure,
        "launches": decisions,
    }


def verify_historical(receipt):
    """Check the resolved results and zero-clock refusal against the sealed replay."""
    archived = load(HERE / "donor_release.json")
    compare(receipt["physical_port_order"], archived["physical_port_order"])
    compare(
        receipt["simultaneous_99_half_width"], archived["simultaneous_99_half_width"]
    )
    compare(
        receipt["commissioning_exhaustion"],
        archived["exhaustion_by_cap"]["96"],
        "root/refused/96",
    )
    for actual, expected in zip(receipt["launches"], archived["launches"], strict=True):
        measured = {name: actual[name] for name in expected}
        if "launches" in expected:
            require(
                actual["disposition"] == "passed", "Resolved historical launch failed"
            )
            measured["launches"] = [
                {name: tested[name] for name in original}
                for tested, original in zip(
                    actual["launches"], expected["launches"], strict=True
                )
            ]
        compare(measured, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clock-receipt", type=Path)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    if args.clock_receipt and not args.record:
        parser.error("--clock-receipt requires a separate --record path")
    source = args.clock_receipt or HERE / "clock.json"
    snapshot = source.read_bytes()
    data = run(loads(snapshot.decode("utf-8")))
    data["clock_receipt_sha256"] = sha256(snapshot).hexdigest()
    content = json.dumps(data, indent=2, allow_nan=False) + "\n"
    if not args.clock_receipt:
        verify_historical(json.loads(content))
    if args.record:
        with args.record.open("x", encoding="utf-8") as file:
            file.write(content)
    dispositions = Counter(row["disposition"] for row in data["launches"])
    print(", ".join(f"{count} {label}" for label, count in dispositions.items()))


if __name__ == "__main__":
    main()
