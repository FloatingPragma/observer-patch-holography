#!/usr/bin/env python3
"""Verify or replay the portable finite component review without editing it."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PACKET = HERE / "completion"
MANIFEST_SHA256 = "4b5d22aae764b3784f73abc239c152f1296ff8d145cc43bcc8f28a2b6ceb60a2"
REVIEW = "workspace/openwave/row_11_lepton_masses/reviews/2026-11-15_codex_completion_review"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_packet(packet=PACKET):
    manifest = packet / "MANIFEST.sha256"
    if digest(manifest) != MANIFEST_SHA256:
        raise ValueError("Portable completion manifest mismatch")
    expected = {}
    for line in manifest.read_text().splitlines():
        sha, name = line.split(None, 1)
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or name in expected:
            raise ValueError("Unsafe or duplicate portable path: " + name)
        expected[name] = sha
    actual = set()
    for path in packet.rglob("*"):
        if path.is_symlink():
            raise ValueError("Symlink in portable completion")
        if path.is_file() and path != manifest:
            name = path.relative_to(packet).as_posix()
            actual.add(name)
            if expected.get(name) != digest(path):
                raise ValueError("Changed or unexpected portable file: " + name)
    if actual != set(expected):
        raise ValueError("Missing portable completion files")
    return {"manifest_sha256": MANIFEST_SHA256, "files_verified": len(expected), "all_pass": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Fresh output directory for the four recorded review jobs")
    args = parser.parse_args()
    result = verify_packet()
    if args.output is not None:
        output = args.output.resolve()
        if output.exists() or output.is_relative_to(PACKET.resolve()):
            raise ValueError("Replay requires a fresh output directory outside the sealed closure")
        env = {key: value for key, value in os.environ.items() if not key.startswith("PYTHON")}
        command = [sys.executable, "-I", "-B", str(PACKET / REVIEW / "run_review.py"), "--output", str(output)]
        process = subprocess.run(command, cwd=PACKET / "workspace", env=env, text=True, capture_output=True)
        output.mkdir(parents=True, exist_ok=True)
        (output / "wrapper.stdout").write_text(process.stdout)
        (output / "wrapper.stderr").write_text(process.stderr)
        result["source_after"] = verify_packet()
        result["child_exit_code"] = process.returncode
        if process.returncode != 0:
            raise ValueError("Component replay failed; complete child output retained")
        summary = json.loads((output / "SUMMARY.json").read_text())
        jobs = summary.get("runs", [])
        if summary.get("all_pass") is not True or len(jobs) != 4 or any(
            job.get("exit_code") != 0 or job.get("all_pass") is not True for job in jobs
        ):
            raise ValueError("Component replay returned a rejected or incomplete verdict")
        result["replay"] = summary
        (output / "PORTABLE_REPLAY.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
