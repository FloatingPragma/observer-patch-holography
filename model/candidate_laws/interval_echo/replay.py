"""Replay patch clocks and reunion cohorts, with additional finite patch measurements in the full run."""
import argparse
import json
from pathlib import Path
import shutil
import tempfile
import sys
from zipfile import BadZipFile

import numpy as np

from model.executable_model import carrier
from tools.strict_json import load
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


def artifacts(full):
    names = ("reunions.npz",)
    if full:
        names += ("observations.json", "patch_records.npz", "consensus.json", "currents.json", "projection.json")
    return (*names, "RESULTS.json")


def artifact_errors(generated, full):
    expected = set(artifacts(full))
    present = {path.name for path in generated.iterdir()} if generated.is_dir() else set()
    failures = [f"missing {name} in replay" for name in sorted(expected) if not (generated / name).is_file()]
    failures.extend(f"unexpected {name} in replay" for name in sorted(present - expected))
    return failures


def check_output(directory, full):
    if not full and any((directory / name).exists() for name in set(artifacts(True)) - set(artifacts(False))):
        raise ValueError("a quick replay needs an output directory without full-only measurements")
    for name in artifacts(full):
        target = directory / name
        if (target.exists() or target.is_symlink()) and not target.is_file():
            raise ValueError(f"output measurement is not a file: {target}")


def publish(stage, directory, full):
    # The backup directory also excludes simultaneous publishers. Backups
    # remain available if the filesystem prevents a rollback.
    backup = directory / ".interval-echo-publication"
    backup.mkdir()
    replaced = []
    try:
        check_output(directory, full)
        for name in artifacts(full):
            if (directory / name).is_file():
                shutil.copy2(directory / name, backup / name, follow_symlinks=False)
        for name in artifacts(full):
            (stage / name).replace(directory / name)
            replaced.append(name)
    except BaseException:
        for name in reversed(replaced):
            prior = backup / name
            if prior.exists() or prior.is_symlink():
                prior.replace(directory / name)
            else:
                (directory / name).unlink()
        shutil.rmtree(backup)
        raise
    shutil.rmtree(backup)


def verify(generated, reference, full):
    failures = artifact_errors(generated, full)
    for name in artifacts(full):
        path, target = generated / name, reference / name
        if not path.is_file():
            continue
        if not target.is_file():
            failures.append(f"missing {name} in reference")
            continue
        try:
            if path.suffix == ".json":
                expected, actual = load(target), load(path)
                if not full and name == "RESULTS.json" and isinstance(expected, dict):
                    # Only the 80- and 320-cell material trips are absent in the quick measurement.
                    if isinstance(expected.get("material"), dict):
                        expected["material"] = {k: v for k, v in expected["material"].items() if k not in ("1", "2")}
                compare(expected, actual, name, failures, True)
            else:
                with np.load(target, allow_pickle=False) as saved, np.load(path, allow_pickle=False) as replay:
                    keys = set(saved.files)
                    if keys != set(replay.files) or len(keys) != len(saved.files) or len(keys) != len(replay.files):
                        failures.append(f"{name}: archive keys")
                        continue
                    for key in replay.files:
                        if saved[key].dtype != replay[key].dtype or not np.array_equal(saved[key], replay[key]):
                            failures.append(f"{name}: {key}")
        except (ValueError, OSError, EOFError, TypeError, BadZipFile) as exc:
            failures.append(f"{name}: {exc}")
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
        try:
            check_output(args.output, args.full)
        except ValueError as exc:
            parser.error(str(exc))
        args.output.mkdir(parents=True, exist_ok=True)
        # Complete the measurements before replacing any published data. Publish the run summary last.
        with tempfile.TemporaryDirectory(dir=args.output) as scratch:
            stage = Path(scratch)
            generate(stage, args.full, args.seed)
            failures = artifact_errors(stage, args.full)
            if failures:
                parser.error(", ".join(failures))
            try:
                publish(stage, args.output, args.full)
            except (OSError, ValueError) as exc:
                parser.error(str(exc))
        print(json.dumps({"candidate": candidate.NAME, "written": str(args.output), "full": args.full}))
        return 0
    with tempfile.TemporaryDirectory() as scratch:
        directory = Path(scratch)
        generate(directory, args.full, args.seed)
        failures = verify(directory, args.verify or HERE, args.full)
    print(json.dumps({"candidate": candidate.NAME, "verified": not failures, "full": args.full,
                      "failures": failures}))
    return int(bool(failures))
