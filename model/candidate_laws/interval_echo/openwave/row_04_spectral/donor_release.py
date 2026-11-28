"""Connect the original spectroscopy preparation to the free two-port read."""

import argparse
import json
from hashlib import sha256
from math import expm1, log, log1p, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.executable_model import federation
from tools.strict_json import load

from . import transfer_clock as clock
from .experiment import BOARD_SHA, CARRIER_DIR, HERE, TIMES, compare
from .graph import DONOR, History, donor_selection
from .local import check_local, require
from .phase_reader import effects
from .source import install, triangles


SHOTS = 2**24
NUMERICAL_ALLOWANCE = 2e-11
BENCHMARK = 0.98


def witness_floor(click_probability, half_width=0.0):
    """Transfer a lower bound on the original-start click probability."""
    return max(0.0, 8 * (click_probability - half_width - NUMERICAL_ALLOWANCE))


def commissioning_failure(operations, cap):
    return -expm1(operations * log1p(-(2.0**-cap)))


def port_records(tree, prepared, evolutions, stochastic, failure):
    arms = []
    read = np.searchsorted(np.unique(tree.read_ports), tree.read_ports)
    for time in TIMES:
        unitary = evolutions[time]
        population = abs(prepared) ** 2
        laws = {
            "coherent": abs(unitary @ prepared) ** 2,
            "dephased": abs(unitary) ** 2 @ population,
            "rate_matched_swaps": stochastic[time] @ population,
        }
        ports = {}
        for name, probability in laws.items():
            require(np.min(probability) >= -NUMERICAL_ALLOWANCE, "Negative probability")
            require(
                abs(probability.sum() - 1) <= NUMERICAL_ALLOWANCE,
                "Port law loses probability before commissioning losses",
            )
            bins = (1 - failure) * np.bincount(read, weights=probability)
            ports[name] = bins.round(12).tolist()
        arms.append({"exchange_exposure": time, "port_probabilities": ports})
    return arms


def run():
    check_local()
    board = federation(0)
    fixture = CARRIER_DIR / "carriers" / "FEDERATION_LEVEL0.npz"
    require(sha256(fixture.read_bytes()).hexdigest() == BOARD_SHA[0], "Board differs")
    heads, operations = install(board, triangles(board)[0])
    tree = History(board, heads)
    choices = donor_selection(tree)["choices"]
    readers = {item["analyzer"]: effects(tree, item["paths"], -1) for item in choices}
    analyzer = load(DONOR / "analyzers.json")
    addresses = np.array([tree.sites == port for port in analyzer["ports"]], dtype=int)
    angles = np.array([word["angles"] for word in analyzer["pulses"]]) @ addresses
    fair_operations = sum(op in ("write", "head_hop") for op, _, _ in operations)
    require(fair_operations == len(board.seams) + 1, "Commissioning count differs")
    budgets = {
        str(cap): commissioning_failure(fair_operations, cap)
        for cap in (1, 8, 16, 32, 96)
    }
    failure = budgets["96"]
    energies, basis = eigh(tree.twice_scaled)
    rates, directions = eigh(tree.lap)
    ports = np.unique(tree.read_ports).tolist()
    clock_receipt = load(HERE / "clock.json")
    settings = clock_receipt["runs"]
    coverage = sqrt(log(2 * len(choices) * len(settings) / 0.01) / (2 * SHOTS))
    result = []
    for index, setting in enumerate(settings):
        strength = setting["strength_in_simulator"]
        calibration = clock.infer(setting["counts"])
        row = {"strength_in_simulator": strength, "clock_calibration": calibration}
        result.append(row)
        if calibration is None:
            require(strength == 0, "Unresolved nonzero clock")
            row["disposition"] = "inconclusive_clock"
            continue
        value, lower, upper, _ = calibration
        require(lower <= abs(strength) <= upper, "Contact calibration misses strength")
        states = clock.run_word(tree, angles, analyzer["step"], strength / value)
        evolutions = {
            t: (basis * np.exp(-0.5j * strength * t * energies)) @ basis.T
            for t in TIMES
        }
        stochastic = {
            t: (directions * np.exp(-0.5 * abs(strength) * t * rates)) @ directions.T
            for t in TIMES
        }
        launches = []
        for choice in choices:
            column = choice["analyzer"]
            prepared = states[:, column]
            effect = readers[column]
            phase_law = np.einsum("i,kij,j->k", prepared.conj(), effect, prepared).real
            require(min(phase_law) >= -NUMERICAL_ALLOWANCE, "Invalid witness law")
            outcome = np.r_[(1 - failure) * phase_law, failure]
            require(abs(outcome.sum() - 1) <= NUMERICAL_ALLOWANCE, "Incomplete witness")
            rounded = outcome.round(12)
            require(min(rounded) >= 0, "Witness sampling grid contains a negative bin")
            count = np.random.default_rng([441304, index, column]).multinomial(
                SHOTS, rounded / rounded.sum()
            )
            require(count.sum() == SHOTS, "Incomplete original-start count bank")
            floor = witness_floor(float(count[-2] / SHOTS), coverage)
            computed = witness_floor(float(outcome[-2]))
            require(min(floor, computed) >= BENCHMARK, "Launch witness below benchmark")
            receipts = port_records(tree, prepared, evolutions, stochastic, failure)
            tips = [ports.index(path[-1]) for path in choice["paths"]]
            for receipt in receipts:
                home = {
                    arm: sum(values[tip] for tip in tips)
                    for arm, values in receipt["port_probabilities"].items()
                }
                require(
                    home["coherent"] >= computed - NUMERICAL_ALLOWANCE,
                    "Native flight violates its protected weight",
                )
                if receipt["exchange_exposure"] == 37:
                    require(
                        max(home["dephased"], home["rate_matched_swaps"]) < BENCHMARK,
                        "A changed-law control does not discriminate persistence",
                    )
                receipt["two_port_probabilities"] = home
            launches.append(
                {
                    "column": column,
                    "launch_witness_probabilities": rounded.tolist(),
                    "simulated_witness_counts": count.tolist(),
                    "computed_all_time_floor": computed,
                    "sampling_all_time_floor": floor,
                    "terminal": receipts,
                }
            )
        row["launches"] = launches
    return {
        "protocol_sha256": sha256((HERE / "DONOR_RELEASE.md").read_bytes()).hexdigest(),
        "clock_receipt_sha256": sha256((HERE / "clock.json").read_bytes()).hexdigest(),
        "carrier_sha256": BOARD_SHA[0],
        "physical_port_order": ports,
        "witness_order": [
            "wrong_prefix",
            "reverse_refusal_1",
            "reverse_refusal_2",
            "reverse_refusal_3",
            "reverse_refusal_4",
            "reverse_refusal_5",
            "parent_vacant",
            "signed_ray_click",
            "commissioning_exhaustion",
        ],
        "fair_operations": fair_operations,
        "exhaustion_by_cap": budgets,
        "shots_per_witness": SHOTS,
        "simultaneous_99_half_width": coverage,
        "terminal_exposures": TIMES,
        "launches": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.record:
        with args.record.open("x", encoding="utf-8", newline="\n") as output:
            json.dump(result, output, separators=(",", ":"), allow_nan=False)
            output.write("\n")
    else:
        archived = load(HERE / "donor_release.json")
        compare(result, archived)
        compare(
            result["exhaustion_by_cap"], archived["exhaustion_by_cap"], "root/refused"
        )
    resolved = [item for row in result["launches"] for item in row.get("launches", ())]
    print(
        f"{len(resolved)} launches, minimum computed free-persistence bound "
        f"{min(item['computed_all_time_floor'] for item in resolved):.10f}, "
        "99% sampling lower bound "
        f"{min(item['sampling_all_time_floor'] for item in resolved):.10f}"
    )


if __name__ == "__main__":
    main()
