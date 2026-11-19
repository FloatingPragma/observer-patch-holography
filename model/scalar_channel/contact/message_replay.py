"""Send a receiver's private sign through repeated use of the same time-bin switch."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..carrier import require
from . import material_radar, radar_checks, remote_message, tagged_radar
from .__main__ import HERE, check_sources, digest

SOURCES = {"model/scalar_channel/contact/" + name for name in (
    "radar_checks_sources.json", "message_replay.py", "remote_message.py")}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--record", type=Path, help="acquire private-message records at a new path")
    mode.add_argument("--read-archive", action="store_true", help="grade original timed histograms directly")
    cli.add_argument("--details", action="store_true", help="print every calibration, decision and resource bound")
    args = cli.parse_args()
    spec = material_radar.protocol()
    check_sources(HERE / "radar_checks_sources.json", radar_checks.SOURCES)
    check_sources(HERE / "message_sources.json", SOURCES)
    _, tag = tagged_radar.run(spec)
    path = HERE / "message_counts.npz"
    if args.record is None:
        require(digest(path) == path.with_suffix(".sha256").read_text().strip(), "private-message archive differs")
    saved = None
    if args.read_archive:
        with np.load(path, allow_pickle=False) as archive:
            saved = archive["remote_message_time_counts"]
    weight_path = HERE / "message_weights.npz"
    require(digest(weight_path) == "568416179b0117fc4c709322cb2d312a20ab9bfad02075a959041a05e42795ee",
            "message sampling distribution differs from the original acquisition")
    with np.load(weight_path, allow_pickle=False) as weights:
        require(weights.files == ["weights"], "invalid message sampling distribution archive")
        report, arrays = remote_message.run(spec, tag["tagged_time_probabilities"],
                                            saved_records=saved, frozen_weights=weights["weights"])
    counts = {key: arrays[key] for key in ("remote_message_time_counts", "remote_message_private_bits")}
    if args.record is not None:
        with args.record.open("xb") as f:
            np.savez_compressed(f, **counts)
        material_radar.verify_counts(args.record, counts)
    else:
        material_radar.verify_counts(path, counts)
    if not args.details:
        report = {key: report[key] for key in ("total_correct_bits", "total_incorrect_bits",
                  "total_abstentions", "weakest_source_dyad", "trials_per_calibration_or_bit")}
    print(json.dumps(report, allow_nan=False, indent=2))


if __name__ == "__main__":
    main()
