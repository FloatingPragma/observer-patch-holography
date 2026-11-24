"""Read a labelled reflection using only the sending station's clock."""

import argparse
import json

import numpy as np

from .capacity import path_records
from .experiment import HERE, layout


def echoes(level, closures, origin, targets, emission, transmitting):
    outward = path_records(level, closures, origin, emission, targets, transmitting)
    returns = []
    for target, hit in zip(targets, outward, strict=True):
        if hit is None:
            returns.append(None)
            continue
        (back,) = path_records(
            level, closures, target, hit["arrival"], [origin], transmitting
        )
        returns.append(None if back is None else back["arrival"])
    return returns


def run():
    cases = json.loads((HERE / "homogeneity.json").read_text())
    records = []
    with (
        np.load(HERE / "signals.npz") as signals,
        np.load(HERE / "records.npz") as survey,
    ):
        for case in cases:
            if case["indices"] != [0, 1]:
                continue
            level, seed, dilute = case["level"], case["seed"], case["dilute"]
            closures = (
                signals[f"{seed}_1_commits"] if dilute else survey[f"epoch_{seed}"]
            )
            witness = case["witness"]
            origin, targets = witness["origin"], witness["destinations"]
            rounds = (0, 4, 48, 52)
            native = [
                echoes(level, closures, origin, targets, t, False) for t in rounds
            ]
            control = [
                echoes(level, closures, origin, targets, t, True) for t in rounds
            ]
            assert all(
                row == [t + 2 * witness["initial_range"]] * 2
                for t, row in zip(rounds, control, strict=True)
            )
            records.append(
                dict(
                    level=level,
                    seed=seed,
                    dilute=dilute,
                    origin=origin,
                    targets=targets,
                    emissions=list(rounds),
                    returns=native,
                    transmitting_control=control,
                    timeout_ticks=2 * layout(level)[0],
                )
            )
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "radar.json"
    records = run()
    if parser.parse_args().write:
        path.write_text(json.dumps(records, indent=1) + "\n")
    else:
        assert records == json.loads(path.read_text())
        print(f"Verified source-clock echo records on {len(records)} repair histories.")


if __name__ == "__main__":
    main()
