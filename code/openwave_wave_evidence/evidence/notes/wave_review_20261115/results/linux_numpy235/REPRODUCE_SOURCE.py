#!/usr/bin/env python3
"""Review sealed wave evidence without modifying it. Run from any directory."""
import argparse
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUNDLE = "shared/unified_wave_submission_20261115"
COMMIT = "60519f4"
FREE = "evidence/free/evidence/openwave/row_30_klein_gordon/work/2026-11-14_codex_common_cubic"
ORBITAL = "evidence/orbital/evidence/openwave/row_31_orbital_levels/work/2026-11-14_codex_production_field"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode())


def checked_path(base, relative):
    path = (base / relative).resolve()
    if not path.is_relative_to(base.resolve()):
        raise ValueError(f"Path outside evidence: {relative}")
    return path


def verify_manifest(base, manifest):
    count = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split(maxsplit=1)
        path = checked_path(base, relative.lstrip("*"))
        if sha(path) != expected:
            raise ValueError(f"Manifest mismatch: {path}")
        count += 1
    return count


def extract(destination):
    commit = subprocess.check_output(
        ["git", "rev-parse", COMMIT + "^{commit}"], cwd=REPO, text=True).strip()
    data = subprocess.check_output(
        ["git", "-c", "core.autocrlf=false", "-c", "core.eol=lf",
         "archive", "--format=tar", commit, BUNDLE], cwd=REPO)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:") as archive:
        for member in archive:
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f"Non-regular archive entry: {member.name}")
            path = checked_path(destination, member.name)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(archive.extractfile(member).read())
    # Git archive can apply checkout conversions: verify every extracted blob.
    tree = subprocess.check_output(["git", "ls-tree", "-r", "-z", commit, "--", BUNDLE], cwd=REPO)
    for entry in tree.split(b"\0"):
        if not entry:
            continue
        meta, name = entry.split(b"\t", 1)
        _, kind, object_id = meta.split()
        if kind != b"blob":
            raise ValueError(f"Unexpected Git object: {entry!r}")
        content = checked_path(destination, name.decode("utf-8")).read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()
        if actual != object_id.decode():
            raise ValueError(f"Archive is not Git-blob-identical: {name!r}")
    return commit, hashlib.sha256(data).hexdigest(), destination / BUNDLE


def compare_json(expected, actual):
    """Diagnostic only: 1e-10 is NOT a replacement for the sealed byte gate."""
    result = dict(float_values=0, changed_float_values=0, max_absolute_difference=0.,
                  max_difference_path=None, above_1e_minus_10=0,
                  exact_mismatch_count=0, exact_mismatch_examples=[])

    def mismatch(path, a, b):
        result["exact_mismatch_count"] += 1
        if len(result["exact_mismatch_examples"]) < 12:
            result["exact_mismatch_examples"].append(dict(path=path, expected=str(a)[:160],
                                                           actual=str(b)[:160]))

    def walk(a, b, path):
        if type(a) is not type(b):
            mismatch(path, type(a).__name__, type(b).__name__)
        elif isinstance(a, dict):
            if a.keys() != b.keys():
                mismatch(path + "/keys", sorted(a), sorted(b))
            for key in a.keys() & b.keys():
                walk(a[key], b[key], path + "/" + key)
        elif isinstance(a, list):
            if len(a) != len(b):
                mismatch(path + "/length", len(a), len(b))
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, path + "/" + str(i))
        elif isinstance(a, float):
            result["float_values"] += 1
            if not math.isfinite(a) or not math.isfinite(b):
                mismatch(path, a, b)
                return
            error = abs(a - b)
            result["changed_float_values"] += int(error != 0)
            result["above_1e_minus_10"] += int(error > 1e-10)
            if error > result["max_absolute_difference"]:
                result["max_absolute_difference"] = error
                result["max_difference_path"] = path
        elif a != b:
            mismatch(path, a, b)

    walk(expected, actual, "")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default="windows")
    parser.add_argument("--wsl-python", help="Absolute Linux Python path; launches via Ubuntu WSL")
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--scratch", type=Path, default=Path(
        os.environ.get("LOCALAPPDATA", tempfile.gettempdir())) / "Temp/opencode")
    args = parser.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="oph_wave_review_", dir=args.scratch))
    output = args.output / args.label
    output.mkdir(parents=True, exist_ok=True)
    (output / "REPRODUCE_SOURCE.py").write_bytes(Path(__file__).read_bytes())
    commit, archive_hash, bundle = extract(work / "committed")
    report = dict(review_started_utc=datetime.now(timezone.utc).isoformat(),
                  subject_commit=commit, bundle=BUNDLE, git_archive_sha256=archive_hash,
                  scratch=str(work), reviewer_script_sha256=sha(Path(__file__)),
                  original_manifest_sha256=sha(bundle / "MANIFEST.sha256"))

    def linux(path):
        return subprocess.check_output(
            ["wsl", "-d", "Ubuntu", "--", "wslpath", "-a", path.as_posix()], text=True).strip()

    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1", PYTHONIOENCODING="utf-8",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")

    def run(name, script=None, flags=(), code=None, optimize=False):
        local_env = env.copy()
        if optimize:
            local_env["PYTHONOPTIMIZE"] = "1"
        if args.wsl_python:
            command = ["wsl", "-d", "Ubuntu", "--", "env", "-u", "PYTHONOPTIMIZE",
                       "PYTHONDONTWRITEBYTECODE=1", "PYTHONUTF8=1", "PYTHONIOENCODING=utf-8",
                       "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1",
                       "TMPDIR=" + linux(work)]
            if optimize:
                command.append("PYTHONOPTIMIZE=1")
            command += [args.wsl_python, "-B", "-X", "utf8"]
            command += ["-c", code] if code is not None else [linux(script), *flags]
        else:
            local_env["TMP"] = local_env["TEMP"] = str(work)
            command = [sys.executable, "-B", "-X", "utf8"]
            command += ["-c", code] if code is not None else [str(script), *flags]
        started = time.monotonic()
        completed = subprocess.run(command, cwd=work, env=local_env, capture_output=True,
                                   timeout=240)
        (output / (name + ".stdout.txt")).write_bytes(completed.stdout)
        (output / (name + ".stderr.txt")).write_bytes(completed.stderr)
        entry = dict(exit_code=completed.returncode, elapsed_seconds=time.monotonic() - started,
                     command=command, optimization_requested=optimize,
                     stdout_sha256=hashlib.sha256(completed.stdout).hexdigest(),
                     stderr_sha256=hashlib.sha256(completed.stderr).hexdigest())
        try:
            entry["json"] = json.loads(completed.stdout)
        except (ValueError, UnicodeError):
            pass
        return entry

    report["environment"] = run("environment", code=(
        "import sys,platform,json,numpy,zlib; "
        "print(json.dumps(dict(python=sys.version,numpy=numpy.__version__,"
        "platform=platform.platform(),optimize=sys.flags.optimize,"
        "zlib=zlib.ZLIB_RUNTIME_VERSION)))"))
    report["outer_manifest_entries"] = verify_manifest(bundle, bundle / "MANIFEST.sha256")
    pins = json.loads((bundle / "SOURCE_PINS.json").read_bytes())
    mismatches = []
    for relative, expected in pins.items():
        exact = checked_path(bundle, relative)
        if sha(exact) != expected:
            raise ValueError(f"Committed pin mismatch: {relative}")
        local = REPO / BUNDLE / relative
        if not local.exists() or sha(local) != expected:
            mismatches.append(dict(path=relative, crlf_only=local.exists() and
                                   local.read_bytes().replace(b"\r\n", b"\n") == exact.read_bytes()))
    report["working_tree_pin_mismatches"] = mismatches
    report["copied_pins_verified"] = len(pins)
    origins = json.loads((bundle / "ORIGINS.json").read_bytes())
    report["ancestor_manifests"] = {}
    for name, origin in origins.items():
        base = bundle / "evidence" / name
        if sha(base / "MANIFEST.sha256") != origin["manifest_sha256"]:
            raise ValueError(f"Origin mismatch: {name}")
        report["ancestor_manifests"][name] = verify_manifest(base, base / "MANIFEST.sha256")
    report["custody"] = run("custody", bundle / "replay.py", ("--custody-only",))
    report["original_replay"] = run("original_replay", bundle / "replay.py")
    save(output / "REPRODUCTION.json", report)

    # Leaf execution diagnoses a failed wrapper without weakening its byte gate.
    report["leaf_diagnostics"] = {}
    for name, relative, script, artifacts in [
        ("free", FREE, "probe.py", ("RECEIPT.json", "RECORDS.json.gz")),
        ("orbital", ORBITAL, "probe.py", ("RECEIPT.json", "RECORDS.json.gz")),
        ("unified", "evidence/unified", "headline.py", ("RECEIPT.json", "CONTROL_RECORDS.json.gz")),
        ("boost", "evidence/boost", "headline.py", ("RECEIPT.json",)),
    ]:
        source = bundle / relative
        target = work / "leaf_runs" / name
        shutil.copytree(source, target)
        entry = run("leaf_" + name, target / script)
        entry["artifacts"] = {}
        for artifact in artifacts:
            a, b = source / artifact, target / artifact
            if not b.exists():
                entry["artifacts"][artifact] = dict(missing=True)
                continue
            ab, bb = a.read_bytes(), b.read_bytes()
            comparison = dict(byte_equal=ab == bb, expected_sha256=sha(a), actual_sha256=sha(b))
            if artifact.endswith(".gz"):
                ab, bb = gzip.decompress(ab), gzip.decompress(bb)
                comparison["uncompressed_byte_equal"] = ab == bb
            comparison["json_diagnostic"] = compare_json(json.loads(ab), json.loads(bb))
            entry["artifacts"][artifact] = comparison
        report["leaf_diagnostics"][name] = entry
        save(output / "REPRODUCTION.json", report)

    report["optimized_original_replay"] = run("optimized_original_replay", bundle / "replay.py",
                                               optimize=True)
    # Change only a top-level origin digest in a disposable copy.
    mutated = work / "mutated_origin"
    shutil.copytree(bundle, mutated)
    changed = json.loads((mutated / "ORIGINS.json").read_bytes())
    changed["free"]["manifest_sha256"] = "0" * 64
    save(mutated / "ORIGINS.json", changed)
    report["origin_mutation_normal"] = run("origin_mutation_normal", mutated / "replay.py",
                                            ("--custody-only",))
    report["origin_mutation_optimized"] = run("origin_mutation_optimized", mutated / "replay.py",
                                               ("--custody-only",), optimize=True)
    # Verify that both the immutable scratch subject and workspace stayed intact.
    report["outer_manifest_entries_after"] = verify_manifest(bundle, bundle / "MANIFEST.sha256")
    report["subject_unchanged"] = all(sha(bundle / p) == h for p, h in pins.items())
    report["review_finished_utc"] = datetime.now(timezone.utc).isoformat()
    save(output / "REPRODUCTION.json", report)
    print(json.dumps(dict(report=str(output / "REPRODUCTION.json"), scratch=str(work),
                          original_exit=report["original_replay"]["exit_code"],
                          optimized_exit=report["optimized_original_replay"]["exit_code"],
                          leaf_exits={k:v["exit_code"] for k,v in report["leaf_diagnostics"].items()}),
                     indent=2))


if __name__ == "__main__":
    main()
