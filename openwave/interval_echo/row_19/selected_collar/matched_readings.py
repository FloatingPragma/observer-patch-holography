"""Test four-tick focusing at matched range, collar-charge, work and current means."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import measure, observer_records

RECEIPT = Path(__file__).with_name("matched_readings.json")
SEED = 190819


def prepare():
    ref, *_ = measure.jacobi.beam.aperture.source_and_changes()
    basis = measure.joint_basis()
    left, right = measure.marginals(ref)
    product = np.array(
        [
            np.kron(left_block, right_block) / np.trace(left_block)
            for left_block, right_block in zip(left, right, strict=True)
        ]
    )
    meters = {
        "radar": measure.station.radar_menus()[6],
        "modular": [
            measure.station.log_blocks(product) - measure.station.log_blocks(ref)
        ],
    }
    formation = measure.jacobi.beam.front.optics.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    rings = measure.jacobi.beam.screens(board, route)
    for charge in (-1, 1):
        for name, observable in measure.stress_effects(
            board, route, rings, charge
        ).items():
            meters[f"{charge}_{name}"] = [observable]

    rows = np.array(
        [
            measure.station.pair(basis, effect)
            for effects in meters.values()
            for effect in effects
        ]
    )
    normalized = rows / np.linalg.norm(rows, axis=1)[:, None]
    _, singular, right_vectors = np.linalg.svd(normalized, full_matrices=True)
    rank = int(np.count_nonzero(singular > 1e-10 * singular[0]))
    kernel = right_vectors[rank:]
    projector = kernel.T @ kernel
    lengths = np.linalg.norm(projector, axis=0)
    coordinate = int(np.argmax(lengths))
    vector = projector[:, coordinate] / lengths[coordinate]
    direction = np.einsum("i,ibjk->bjk", vector, basis)

    for power in range(7, 28):
        step = 2.0**-power
        ends = [ref + sign * step * direction for sign in (-1, 1)]
        if min(np.linalg.eigvalsh(x).min() for x in ends) <= 0:
            continue
        recipes = [measure.station.preparation(x) for x in ends]
        if min(min(p["min_kept"], p["min_dephased"]) for p in recipes) > 0:
            break
    else:
        raise ArithmeticError("No positive two-sided preparation")

    levels = {
        key: [measure.station.pair(np.asarray(effects), x).tolist() for x in ends]
        for key, effects in meters.items()
    }
    gaps = {k: float(np.max(np.abs(np.subtract(*v)))) for k, v in levels.items()}
    marginal_gap = float(
        np.max(
            np.abs(np.subtract(measure.marginals(ends[0]), measure.marginals(ends[1])))
        )
    )
    if max(*gaps.values(), marginal_gap, *(p["max_error"] for p in recipes)) > 1e-11:
        raise ArithmeticError("Unmatched source or invalid preparation")
    forward, reverse = [
        measure.focusing(board, route, ring, -1, ref, {"joint": direction})["joint"]
        for ring in (rings, rings[::-1])
    ]
    summary = dict(
        coordinate_zero_based=coordinate,
        rank=rank,
        singular_values=singular.tolist(),
        selected_coefficients=vector.tolist(),
        half_step=step,
        max_normalized_direction_residual=float(np.max(np.abs(normalized @ vector))),
        max_marginal_difference=marginal_gap,
        endpoint_meter_means=levels,
        max_endpoint_differences=gaps,
        realizations=recipes,
        differential_odd_focusing=float(forward - reverse),
    )
    return direction, summary


def run(precision):
    if precision:
        raw = RECEIPT.read_bytes()
        source = json.loads(raw)["source"]
        direction = np.einsum(
            "i,ibjk->bjk", source["selected_coefficients"], measure.joint_basis()
        )
        record = observer_records.collect(
            direction, source["half_step"], 4, SEED + 1, precision=True
        )
        return dict(
            primary_sha256=sha256(raw).hexdigest(),
            counts=record,
            decoded=observer_records.decode(record),
        )
    direction, summary = prepare()
    record = observer_records.collect(direction, summary["half_step"], 4, SEED)
    decoded = observer_records.decode(record)
    return dict(source=summary, counts=record, decoded=decoded)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="create the acquisition receipt"
    )
    parser.add_argument(
        "--counts-only", action="store_true", help="decode the existing counts"
    )
    parser.add_argument(
        "--precision", action="store_true", help="repeat on the frozen primary source"
    )
    args = parser.parse_args()
    if args.write and args.counts_only:
        parser.error("--write requires a fresh acquisition")
    target = RECEIPT.with_name("matched_precision.json") if args.precision else RECEIPT
    if args.counts_only:
        data = json.loads(target.read_text())
        if (
            args.precision
            and data["primary_sha256"] != sha256(RECEIPT.read_bytes()).hexdigest()
        ):
            raise ValueError("Precision bank refers to a different primary source")
        if observer_records.decode(data["counts"]) != data["decoded"]:
            raise ValueError("Decoded intervals differ from the receipt")
    else:
        data = run(args.precision)
        if args.write:
            with target.open("x") as handle:
                json.dump(data, handle, allow_nan=False, indent=2)
                handle.write("\n")
        else:
            saved = json.loads(target.read_text())
            if args.precision:
                if data["primary_sha256"] != saved["primary_sha256"]:
                    raise ValueError(
                        "Primary source differs from the saved acquisition"
                    )
            else:
                if data["source"] != saved["source"]:
                    raise ValueError(
                        "Reconstructed input meters or source recipe differ"
                    )
            if data["counts"] != saved["counts"] or data["decoded"] != saved["decoded"]:
                raise ValueError("Simulated records differ from the saved acquisition")
    print(
        json.dumps(
            {
                k: data["decoded"][k]
                for k in (
                    "predicted_odd_focusing",
                    "measured_odd_focusing",
                    "measured_minus_predicted",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
