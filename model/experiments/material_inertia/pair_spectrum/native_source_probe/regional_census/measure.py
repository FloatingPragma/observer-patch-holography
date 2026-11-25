"""Recompute whole-band charge-event bounds and complete adaptive census records."""

import argparse
from fractions import Fraction
import json
from math import ceil, floor, log
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import measure as probe
from ..count_readback.measure import interval
from ..measure import check, flight, previous
from . import geometry, meter

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "results.json"
SEED = 7191931


def operator_bound(basis, charged, projector_error):
    """Maximum of the charge-event probability on an arbitrary cluster state."""
    sub = basis[charged]
    gram = sub.T @ sub
    values, vectors = np.linalg.eigh(gram)
    residual = np.linalg.norm(gram @ vectors - vectors * values, 2)
    ortho = np.linalg.norm(vectors.T @ vectors - np.eye(len(values)), 2)
    # Standard dot-product error bound for the entire 24x24 Gram matrix.
    eps = np.finfo(float).eps
    gamma = len(sub) * eps / (1 - len(sub) * eps)
    products = gamma * np.linalg.norm(sub, "fro") ** 2
    # The producer checks the orthonormal basis defect below 1e-11.
    allowance = 2 * projector_error + products + residual / (1 - ortho) + ortho + 2e-11
    maximum = float(values[-1])
    return dict(
        maximum=maximum,
        numerical_allowance=float(allowance),
        upper=float(min(1, maximum + allowance)),
        lower=float(max(0, maximum - allowance)),
    )


def acquisition(weights, bins, admission, rng, charged_mask, fallback_mask):
    law = np.bincount(bins, weights=weights, minlength=geometry.BINS) * admission
    law[geometry.FILTER_MISS] = 1 - admission
    result = probe.sample(law, rng)
    counts = np.array(result["counts"])
    success = int(counts.sum() - counts[geometry.FILTER_MISS])
    charged = int(counts[charged_mask].sum())
    bounds = interval(charged, success, log(2 * 540 / 0.01))
    return dict(
        accepted=success,
        charged=charged,
        fallback=int(counts[fallback_mask].sum()),
        simultaneous_interval=[
            floor(bounds[0] * 1e10) / 1e10,
            ceil(bounds[1] * 1e10) / 1e10,
        ],
        counts=[[int(i), int(counts[i])] for i in np.flatnonzero(counts)],
        sampling_l1_complement_bound=result["sampling_l1_complement_bound"],
    )


def run():
    check(probe.SHOTS == 2**26, "source count denominator changed")
    local = meter.local_blocks()
    calibration = meter.two_target_calibration()
    decoder = calibration["decoder"]
    board = flight.federation(0)
    profiles = geometry.survey(board)
    x, y, _, kinetic, h, cost = flight.defect.generators(board)
    original = json.loads((flight.HERE / "results.json").read_text())
    reduction = flight.PairReduction(h, cost)
    energies, basis, residual, _ = reduction.lowest(h, 24)
    check(
        reduction.inertia(flight.CUT)["upper_lines_below"] == 24,
        "native cluster count differs",
    )
    projector_error = max(
        original["cluster_projector_error"], residual / (flight.CUT - energies[-1])
    )
    check(
        np.linalg.norm(basis.T @ basis - np.eye(24), 2) < 1e-11, "invalid cluster basis"
    )
    any_separation = operator_bound(basis, x // 12 != y // 12, projector_error)
    labels = []
    for profile in profiles:
        bins, q = geometry.labels(x, y, profile, decoder)
        labels.append(bins)
        profile["charge_on_entire_cluster"] = operator_bound(
            basis, q != 0, projector_error
        )
    masks = geometry.event_masks(decoder)
    rng = np.random.default_rng(SEED)
    launches = []
    for launch in original["flights"]:
        if launch["owner"] not in (0, 5, 10) or launch["admission"] < 0.001:
            continue
        address = flight.defect.ordered_index(*launch["ports"], board.n)
        overlap = basis[address]
        admission = float(overlap @ overlap)
        state = basis @ (
            np.exp(-1j * (energies - energies[0]) * probe.TIME)
            * overlap
            / np.sqrt(admission)
        )
        weight, enorm = previous.normalized_weights(state)
        matched, match_info = previous.hard_core_control(weight, x, y, board.n)
        changed, cnorm = previous.normalized_weights(
            expm_multiply(-4j * kinetic, state)
        )
        panels = {}
        for label, law in [
            ("native", weight),
            ("matched_marginals", matched),
            ("exchange_only", changed),
        ]:
            panels[label] = [
                acquisition(law, b, admission, rng, *masks) for b in labels
            ]
        launches.append(
            dict(
                owner=launch["owner"],
                ports=launch["ports"],
                admission=admission,
                native_per_event_error=float(
                    2 * projector_error / np.sqrt(admission)
                    + (probe.TIME + 4) * residual
                    + max(enorm, cnorm)
                    + 1e-9
                ),
                marginal_residual=match_info["maximum_marginal_residual"],
                panels=panels,
            )
        )
    check(
        len(launches) == 3 and len(profiles) == 60,
        "missing an enrolled source or region",
    )
    native_upper = max(
        row["simultaneous_interval"][1]
        for cohort in launches
        for row in cohort["panels"]["native"]
    )
    control_lower = min(
        cohort["panels"][mode][i]["simultaneous_interval"][0]
        for cohort in launches
        for mode in ("matched_marginals", "exchange_only")
        for i, p in enumerate(profiles)
        if p["center"] == cohort["owner"] and p["radius"] == 0
    )
    uniform = max(p["charge_on_entire_cluster"]["upper"] for p in profiles)
    return dict(
        candidate="interval_echo",
        shots=probe.SHOTS,
        seed=SEED,
        exposure=probe.TIME,
        full_material_input_words=board.n * (board.n - 1),
        spectral_projector_error=float(projector_error),
        total_output_bins=geometry.BINS,
        sparse_encoding="[index,count]; unlisted indices have zero count",
        pointer_decoder=decoder,
        local_pulse_checks=local,
        calibration=calibration,
        omitted_y=meter.two_target_calibration(omit_y=True),
        half_exposure=meter.two_target_calibration(Fraction(1, 2)),
        universal_distinct_owner_bound=any_separation,
        profiles=profiles,
        launches=launches,
        tests=dict(
            maximum_native_upper=native_upper,
            minimum_launch_cell_control_lower=control_lower,
            all_time_all_vectors_upper=uniform,
            below_native_002=native_upper < 0.002,
            controls_above_20=control_lower > 0.20,
            whole_cluster_below_01=uniform < 0.01,
        ),
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        probe.compare(result, json.loads(RECEIPT.read_text()))
        print(
            "Verified native regional charge on the entire filtered moving-pair ensemble."
        )


if __name__ == "__main__":
    main()
