#!/usr/bin/env python3
"""Capture the regression suite and its actual process exit on either platform."""
import argparse
from pathlib import Path
import subprocess
import sys
import replay as r

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--run", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True, help="New log directory")
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=False)
completed = subprocess.run([sys.executable, "-B", "-X", "utf8", str(r.HERE / "test_replay.py"),
                            "--run", str(args.run.resolve()), "-v"], capture_output=True)
for name, data in [("stdout", completed.stdout), ("stderr", completed.stderr)]:
    (args.output / (name + ".txt")).write_bytes(data)
r.save(args.output / "TESTS.json", dict(exit_code=completed.returncode, finished_utc=r.utc(),
    registration_sha256=r.verify_registration(), test_source_sha256=r.sha(r.HERE / "test_replay.py"),
    stdout_sha256=r.sha_bytes(completed.stdout), stderr_sha256=r.sha_bytes(completed.stderr)))
print(completed.stderr.decode("utf-8"))
raise SystemExit(completed.returncode)
