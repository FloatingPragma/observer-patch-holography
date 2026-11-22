"""Replay patch clocks and reunion cohorts, with additional finite patch measurements in the full run."""
import argparse
import json
from pathlib import Path
import tempfile
import sys

import numpy as np

from model.executable_model import carrier
from ..contrast_memory_spherical_v1.verify import compare
from . import candidate, clock, cohorts

HERE = Path(__file__).resolve().parent


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def generate(directory, full, seed):
    levels = (0, 1, 2) if full else (0,)
    material = {}
    for level in levels:
        sphere = carrier.federation(level)
        boards = {}
        for width in (3, 4):
            # Keep the common domain, both capacity boundaries and their local resets in one physical trip.
            panels = (2, 6, 8, 15, 16, 2)
            boards[str(width)] = clock.material_replay(sphere, panels, age_bits=width, seed=seed + level)
        material[str(level)] = {"seam_sha256": carrier.seam_digest(sphere), "layouts": boards}
    cohort, words = cohorts.sample(seed)
    rates = {str(width): [clock.interval_response(q, age_bits=width) for q in range(1, 18)]
             for width in (3, 4)}
    for width, rows in rates.items():
        for row in rows:
            direct = cohorts.suppression(row["interval"], int(width))
            if abs(direct - row["closing_candidate_suppression"]) > 1e-13:
                raise ArithmeticError("first-return enumeration disagrees with finite walk")

    # Seed 31 exhibits a positive primary and a first return at q=8; it is a selected mechanism witness.
    witness = {}
    for width in (3, 4):
        trip = clock.material_replay(carrier.federation(0), (8, 1), age_bits=width, seed=31)
        witness[str(width)] = {"intervals": trip["intervals"], "work_bits": trip["work_bits"],
                               "moving_cells": trip["moving_cells"], "closing": trip["local_events"][7],
                               "reunion": trip["reunion"]}
    result = {"schema": "oph.interval_echo.v1", "candidate": candidate.NAME, "seed": seed,
              "material": material, "selected_eighth_service_seed_31": witness,
              "finite_clock": rates, "cohorts": cohort}
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / "RESULTS.json", result)
    np.savez_compressed(directory / "reunions.npz", **words)
    if full:
        from . import consensus, currents, observations, projection
        observations.run(seed=seed + 4, output=directory)
        write_json(directory / "consensus.json", consensus.run())
        write_json(directory / "currents.json", currents.run())
        write_json(directory / "projection.json", projection.run())


def verify(generated, reference, full):
    failures = []
    for path in sorted(generated.iterdir()):
        target = reference / path.name
        if not target.is_file():
            failures.append(f"missing {path.name}")
        elif path.suffix == ".json":
            compare(json.loads(target.read_text()), json.loads(path.read_text()), "", failures, full)
        elif path.suffix == ".npz":
            with np.load(target, allow_pickle=False) as saved, np.load(path, allow_pickle=False) as replay:
                if set(saved.files) != set(replay.files):
                    failures.append(f"{path.name}: archive keys")
                    continue
                for key in replay.files:
                    if saved[key].dtype != replay[key].dtype or not np.array_equal(saved[key], replay[key]):
                        failures.append(f"{path.name}: {key}")
    return failures


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ['--clock-trip']:
        return clock.trip_cli(argv[1:], default_bits=5)
    if argv[:1] == ['--repair-standard']:
        from model.experiments.interval_echo_consensus.modular_standard import main as repair_standard
        return repair_standard(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__,
                                     epilog="Use --clock-trip [--level 0|1|2] [--age-bits 3|4|5] OUT BACK ... for a material journey. Use --repair-standard to reconstruct and test the seam generator from local offers.")
    parser.add_argument("--full", action="store_true",
                        help="also replay 80/320-cell trips, bulk, repair forks, port currents and the two-patch information projection")
    parser.add_argument("--seed", type=int, default=cohorts.SEED)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--output", type=Path, help="write observations to a directory")
    mode.add_argument("--verify", type=Path, help="compare with a reference directory (default: candidate package)")
    args = parser.parse_args(argv)
    if args.seed < 0:
        parser.error("the seed must be nonnegative")
    if args.output:
        if args.output.resolve() == HERE and (not args.full or args.seed != cohorts.SEED):
            parser.error("a partial or alternative-seed run needs its own output directory")
        args.output.mkdir(parents=True, exist_ok=True)
        # Complete the measurements before replacing any published data. Publish the run summary last.
        with tempfile.TemporaryDirectory(dir=args.output) as scratch:
            stage = Path(scratch)
            generate(stage, args.full, args.seed)
            for artifact in sorted(stage.iterdir(), key=lambda f: f.name == "RESULTS.json"):
                artifact.replace(args.output / artifact.name)
        print(json.dumps({"candidate": candidate.NAME, "written": str(args.output), "full": args.full}))
        return 0
    with tempfile.TemporaryDirectory() as scratch:
        directory = Path(scratch)
        generate(directory, args.full, args.seed)
        failures = verify(directory, args.verify or HERE, args.full)
    print(json.dumps({"candidate": candidate.NAME, "verified": not failures, "full": args.full,
                      "failures": failures}))
    return int(bool(failures))
