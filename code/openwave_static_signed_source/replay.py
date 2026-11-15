#!/usr/bin/env python3
"""Replay the signed-source evidence without changing its archived bytes.

The endpoint rule is j=(d-(2c-1)(abs(d) mod 2))/2. The compact
transfer is T=I-L/1260; its selected energy is -1260 log(1-6/1260).
The only execution-copy amendment removes an obsolete bytecode-cache pin.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"
CACHE = "openwave/shared/paired_source_moments/__pycache__/run.cpython-313.pyc"
PREREG = "openwave/shared/static_candidate_join/PREREG.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    """Verify closure, every published byte, and the one cache omission."""
    manifest = json.loads((HERE / "MANIFEST.json").read_text())["files"]
    actual = {
        p.relative_to(EVIDENCE).as_posix() for p in EVIDENCE.rglob("*") if p.is_file()
    }
    if actual != set(manifest):
        raise ValueError("Evidence file set differs from MANIFEST.json")
    for name, digest in manifest.items():
        path = EVIDENCE / name
        if path.is_symlink() or not path.resolve().is_relative_to(EVIDENCE):
            raise ValueError("Unsafe evidence path: " + name)
        if sha(path) != digest:
            raise ValueError("Evidence hash mismatch: " + name)
    amendment = json.loads((HERE / "PORTABILITY.json").read_text())
    if (
        amendment["omitted_cache"] != CACHE
        or sha(EVIDENCE / PREREG) != amendment["prereg_sha256"]
    ):
        raise ValueError(
            "Portability amendment does not match the original preregistration"
        )
    original = dict(
        line.split("  ", 1)[::-1]
        for line in (EVIDENCE / "MANIFEST.sha256").read_text().splitlines()
    )
    if original.pop(CACHE) != amendment["omitted_cache_sha256"]:
        raise ValueError("Cache omission is not the declared historical pin")
    if original != {k: v for k, v in manifest.items() if k != "MANIFEST.sha256"}:
        raise ValueError("Scientific export differs from the historical manifest")
    prereg = json.loads((EVIDENCE / PREREG).read_text())
    if prereg["files"].pop(CACHE) != amendment["omitted_cache_sha256"]:
        raise ValueError("Unexpected preregistration cache pin")
    for name, digest in prereg["files"].items():
        if manifest.get(name) != digest:
            raise ValueError("Preregistered scientific input is not bound: " + name)
    # Restore provider paths only from already verified, identical export bytes.
    provider = EVIDENCE / "openwave/shared/localized_rest_join"
    pins = json.loads((provider / "PREREG.txt").read_text())["files"]
    pins.update(json.loads((provider / "SOURCE_PINS.json").read_text()))
    by_hash = {digest: name for name, digest in manifest.items()}
    aliases = {}
    for name, digest in pins.items():
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("Unsafe provider path")
        if name in manifest:
            if manifest[name] != digest:
                raise ValueError("Provider input hash differs: " + name)
        else:
            aliases[name] = by_hash[digest]
    declared = json.loads((HERE / "PROVIDER_PATHS.json").read_text())
    if declared["unresolved"] or aliases != declared["aliases"]:
        raise ValueError("Provider aliases differ from the frozen input pins")
    return manifest, prereg, aliases


def compare(actual, expected, skip=()):
    """Compare every retained field; timings and separately hashed artifacts differ."""
    if not isinstance(actual, dict) or set(actual) - set(skip) != set(expected) - set(
        skip
    ):
        raise ValueError("Receipt fields differ")
    for key, value in expected.items():
        if key in skip:
            continue
        got = actual[key]
        if isinstance(value, dict):
            compare(got, value)
        elif isinstance(value, float):
            if (
                not isinstance(got, (int, float))
                or not math.isfinite(got)
                or abs(got - value) > 1e-10
            ):
                raise ValueError("Numeric mismatch: " + key)
        elif got != value:
            raise ValueError("Exact mismatch: " + key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    manifest, prereg, aliases = verify()
    if args.verify_only:
        print(
            json.dumps(
                {
                    "verified_files": len(manifest),
                    "omitted_cache_pins": 1,
                    "identical_provider_aliases": len(aliases),
                },
                sort_keys=True,
            )
        )
        return
    jobs = []
    with tempfile.TemporaryDirectory(prefix="oph-static-") as temporary:
        root = Path(temporary) / "evidence"
        shutil.copytree(EVIDENCE, root)
        for destination, source in aliases.items():
            target = root / destination
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / source, target)
        # The sealed original stays byte-identical; only this scratch copy changes.
        (root / PREREG).write_text(json.dumps(prereg, indent=2, sort_keys=True) + "\n")
        env = os.environ.copy()
        env.update(
            PYTHONDONTWRITEBYTECODE="1",
            OPENBLAS_NUM_THREADS="1",
            OMP_NUM_THREADS="1",
            PYTHONHASHSEED="0",
        )
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONOPTIMIZE", None)
        for name in (
            "compact_static_adjudication",
            "localized_rest_join",
            "static_candidate_join",
        ):
            command = [sys.executable, "-B", "openwave/shared/" + name + "/headline.py"]
            started = time.perf_counter()
            result = subprocess.run(
                command, cwd=root, env=env, text=True, capture_output=True, timeout=180
            )
            if result.returncode:
                # Tracebacks stay local; never include machine paths in the JSON receipt.
                sys.stderr.write(result.stderr)
                raise RuntimeError("Replay failed: " + name)
            data = json.loads(result.stdout)
            if data.get("all_pass") is not True:
                raise ValueError("Replay did not pass: " + name)
            if name == "compact_static_adjudication":
                expected = json.loads(
                    (root / ("openwave/shared/" + name + "/RECEIPT.json")).read_text()
                )
                compare(data, expected, skip=("elapsed_seconds", "artifact_sha256"))
            if name == "localized_rest_join":
                receipt = "openwave/shared/localized_rest_join/READER_RECEIPT.json"
                expected = json.loads((EVIDENCE / receipt).read_text())
                compare(data, expected, skip=("elapsed_seconds",))
                # This historical reader rewrites its timing receipt. Restore
                # the frozen bytes before the next script checks its input pin.
                shutil.copyfile(EVIDENCE / receipt, root / receipt)
            jobs.append(
                {
                    "name": name,
                    "command": "python3 -B " + command[-1],
                    "exit": result.returncode,
                    "wall_seconds": time.perf_counter() - started,
                    "result": data,
                }
            )
    verify()
    print(
        json.dumps(
            {
                "schema": "oph.static.replay.v1",
                "candidate": "A-static-signed-source-v1",
                "all_pass": True,
                "verified_files": len(manifest),
                "omitted_cache_pins": 1,
                "identical_provider_aliases": len(aliases),
                "scientific_inputs_changed": False,
                "jobs": jobs,
                "official_rows_solved": [],
                "external_review": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
