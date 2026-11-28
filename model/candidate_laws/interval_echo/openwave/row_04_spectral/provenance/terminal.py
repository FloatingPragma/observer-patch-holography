"""Re-read the donor's signed material history after independent free releases."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from tools.strict_json import load
from ..donor_release import BENCHMARK, commissioning_failure
from ..donor_release_checked import clock_outcomes
from ..experiment import HERE as ARCHIVE, TIMES
from ..graph import DONOR, History, donor_selection
from ..local import check_local
from ..phase_reader import effects
from ..source import install, triangles
from ..transfer_clock import run_word


HERE = Path(__file__).resolve().parent
# Summed roundoff allowances for numerical source and flight comparisons.
# The all-exposure result uses the exact generator and detector identities.
PAD = 2e-9
ARMS = (
    "native_coherent",
    "dephased_at_release",
    "classical_swaps_same_graph",
    "one_tip_phase_flipped",
)
RECORDS = (
    "wrong_prefix",
    "reverse_query_refusal_1",
    "reverse_query_refusal_2",
    "reverse_query_refusal_3",
    "reverse_query_refusal_4",
    "reverse_query_refusal_5",
    "final_parent_vacant",
    "signed_ray_hit",
    "commissioning_exhaustion",
)


def launch_fixture():
    board = federation(0)
    heads, operations = install(board, triangles(board)[0])
    fair = sum(job in ("write", "head_hop") for job, _, _ in operations)
    assert fair == len(board.seams) + 1 == 631
    return History(board, heads), commissioning_failure(fair, 96)


def complete_bank(effect, rho, failure):
    law = np.array([np.trace(e @ rho).real for e in effect])
    assert min(law) > -PAD and abs(sum(law) - 1) < PAD
    return np.r_[(1 - failure) * law, failure]


def read_after_wait(tree, prep, choice, u, stochastic, failure):
    bank = effects(tree, choice["paths"], -1)
    reference = np.trace(bank[-1])
    assert reference == 0.125
    left, right = (tree.index[tuple(path)] for path in choice["paths"])
    d = np.zeros(len(tree.paths), dtype=int)
    d[left], d[right] = 1, -1
    np.testing.assert_array_equal(tree.twice_scaled @ d, d)
    np.testing.assert_array_equal(bank[-1] * 16, np.outer(d, d))
    sharp = abs(prep) ** 2
    opposite = prep.copy()
    # On this one-load component a native local phase on a physical tip
    # address multiplies all configurations at that port by minus one.
    opposite[tree.sites == tree.sites[right]] *= -1
    coherent = u @ prep
    flipped = u @ opposite
    states = (
        np.outer(coherent, coherent.conj()),
        (u * sharp) @ u.conj().T,
        np.diag(stochastic @ sharp),
        np.outer(flipped, flipped.conj()),
    )
    bins = [complete_bank(bank, rho, failure) for rho in states]
    return {label: vector.tolist() for label, vector in zip(ARMS, bins, strict=True)}


def run():
    check_local()
    tree, failure = launch_fixture()
    choices = donor_selection(tree)["choices"]
    words = load(DONOR / "analyzers.json")
    site_select = np.array([tree.sites == p for p in words["ports"]], dtype=int)
    angle_bank = np.array([entry["angles"] for entry in words["pulses"]]) @ site_select
    report = []
    for record in clock_outcomes(load(ARCHIVE / "clock.json")):
        if record["disposition"] != "calibrated":
            assert record["disposition"] == "inconclusive_clock"
            report.append(
                {k: record[k] for k in ("strength_in_simulator", "disposition")}
            )
            continue
        strength = record["strength_in_simulator"]
        calibrated = record["clock_calibration"][0]
        source = run_word(tree, angle_bank, words["step"], strength / calibrated)
        # Recompute small full-flight matrices by exponential, independently
        # of the eigensystem propagation in the archived release driver.
        native = {t: expm(-0.5j * t * strength * tree.twice_scaled) for t in TIMES}
        swaps = {t: expm(-0.5 * t * abs(strength) * tree.lap) for t in TIMES}
        out = []
        for choice in choices:
            column = choice["analyzer"]
            marked = tuple(tree.index[tuple(path)] for path in choice["paths"])
            prep = source[:, column]
            target_weight = (
                (1 - failure) * abs(prep[marked[0]] - prep[marked[1]]) ** 2 / 2
            )
            floor = float(target_weight) - PAD
            assert floor > BENCHMARK
            appointments = []
            for t in TIMES:
                read = read_after_wait(tree, prep, choice, native[t], swaps[t], failure)
                sign = {law: value[-2] * 8 for law, value in read.items()}
                assert abs(sign["native_coherent"] - target_weight) < PAD
                assert sign["dephased_at_release"] <= 0.5 + PAD
                assert sign["classical_swaps_same_graph"] <= 0.5 + PAD
                assert sign["one_tip_phase_flipped"] <= 1 - floor + PAD
                appointments.append({"exposure": t, "complete_probability_banks": read})
            out.append(
                {
                    "column": column,
                    "tips": [int(tree.sites[i]) for i in marked],
                    "all_exposure_signed_ray_floor_per_original": floor,
                    "terminal_appointments": appointments,
                }
            )
        report.append(
            {
                "strength_in_simulator": strength,
                "disposition": "verified_positive_signed_return",
                "calibrated_strength_magnitude": calibrated,
                "source_commissioning_refusal": failure,
                "launches": out,
            }
        )
    return {
        "records": RECORDS,
        "probability_banks_are_predictions_not_draws": True,
        "proof_uncertainty_allowance": PAD,
        "benchmark": BENCHMARK,
        "launches": report,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.record:
        with args.record.open("x") as output:
            output.write(json.dumps(data, separators=(",", ":")) + "\n")
    else:
        from ..experiment import compare

        compare(data, json.loads((HERE / "terminal.json").read_text()), "terminal")
    resolved = [s for e in data["launches"] for s in e.get("launches", [])]
    smallest = min(x["all_exposure_signed_ray_floor_per_original"] for x in resolved)
    print(f"Signed return floor {smallest:.12f} on {len(resolved)} resolved launches")


if __name__ == "__main__":
    main()
