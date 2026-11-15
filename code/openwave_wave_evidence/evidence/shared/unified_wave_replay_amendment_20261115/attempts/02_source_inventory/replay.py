#!/usr/bin/env python3
"""Portable, fail-closed replay of the sealed wave evidence; see AMENDMENT.md."""
import argparse
import copy
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FREE = "evidence/free/evidence/openwave/row_30_klein_gordon/work/2026-11-14_codex_common_cubic"
ORBITAL = "evidence/orbital/evidence/openwave/row_31_orbital_levels/work/2026-11-14_codex_production_field"
LEAVES = {
    "free": (FREE, "probe.py", "RECORDS.json.gz", "raw_uncompressed_sha256", "actual_native_calls", 169364),
    "orbital": (ORBITAL, "probe.py", "RECORDS.json.gz", "raw_uncompressed_sha256", "native_calls", 347533),
    "unified": ("evidence/unified", "compare.py", "CONTROL_RECORDS.json.gz", "control_raw_uncompressed_sha256", "total_native_calls", 492837),
    "boost": ("evidence/boost", "check_boost.py", None, None, None, 0),
}
PATCHES = (
    ("u=proj[:,np.argmax(np.diag(proj))].copy()", "u=proj[:,0].copy()"),
    ("sv,su=np.linalg.eigh(sl);u=internal_vectors[band_index]",
     "sv,su=np.linalg.eigh(sl);u=internal_vectors[band_index]\n"
     "            if su[0,0]<0:su[:,0]*=-1"),
)


class ReplayError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise ReplayError(message)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return sha_bytes(path.read_bytes())


def utc():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8"))


def strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def finite(value):
        result = float(value)
        require(math.isfinite(result), f"Nonfinite JSON number: {value}")
        return result

    def invalid(value):
        raise ReplayError(f"Invalid JSON constant: {value}")

    return json.loads(data, object_pairs_hook=pairs, parse_float=finite, parse_constant=invalid)


def checked_path(root, relative):
    require(isinstance(relative, str) and relative and "\\" not in relative,
            f"Invalid relative path: {relative!r}")
    require(not Path(relative).is_absolute() and ":" not in relative,
            f"Absolute path: {relative}")
    path = (root / relative).resolve()
    require(path.is_relative_to(root.resolve()) and path != root.resolve(),
            f"Path escapes package: {relative}")
    return path


def verify_manifest(base, manifest):
    seen = set()
    for line in manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split(maxsplit=1)
        relative = relative.lstrip("*")
        require(relative not in seen, f"Duplicate manifest entry: {relative}")
        seen.add(relative)
        require(sha(checked_path(base, relative)) == expected, f"Manifest mismatch: {relative}")
    require(bool(seen), f"Empty manifest: {manifest}")
    return len(seen)


def extract_git(destination, policy):
    """Read raw blobs, bypassing checkout filters and line-ending settings."""
    commit, bundle = policy["source_commit"], policy["bundle"]
    tree = subprocess.check_output(["git", "ls-tree", "-r", "-z", commit, "--", bundle], cwd=REPO)
    entries = []
    for entry in tree.split(b"\0"):
        if not entry:
            continue
        meta, name = entry.split(b"\t", 1)
        mode, kind, oid = meta.split()
        require(mode == b"100644" or mode == b"100755", f"Unsupported Git mode: {mode!r}")
        require(kind == b"blob", "Non-blob subject")
        entries.append((name.decode("utf-8"), oid))
    require(bool(entries), "No subject blobs at pinned commit")
    raw = subprocess.check_output(["git", "cat-file", "--batch"], cwd=REPO,
                                  input=b"".join(oid + b"\n" for _, oid in entries))
    stream = io.BytesIO(raw)
    for name, oid in entries:
        actual_oid, kind, length = stream.readline().split()
        require(actual_oid == oid and kind == b"blob", f"Git batch mismatch: {name}")
        content = stream.read(int(length))
        require(stream.read(1) == b"\n", "Invalid Git batch framing")
        digest = hashlib.sha1(b"blob " + length + b"\0" + content).hexdigest().encode()
        require(digest == oid, f"Git blob mismatch: {name}")
        path = checked_path(destination, name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    require(stream.read() == b"", "Unexpected Git batch data")
    return destination / bundle


def custody(bundle, policy):
    require(sha(bundle / "MANIFEST.sha256") == policy["outer_manifest_sha256"],
            "Original outer manifest digest mismatch")
    outer = verify_manifest(bundle, bundle / "MANIFEST.sha256")
    pins = strict_json((bundle / "SOURCE_PINS.json").read_bytes())
    for relative, digest in pins.items():
        require(sha(checked_path(bundle, relative)) == digest, f"Copied pin mismatch: {relative}")
    origins = strict_json((bundle / "ORIGINS.json").read_bytes())
    ancestors = {}
    for name, origin in origins.items():
        base = checked_path(bundle, "evidence/" + name)
        manifest = base / "MANIFEST.sha256"
        require(sha(manifest) == origin["manifest_sha256"], f"Origin digest mismatch: {name}")
        ancestors[name] = verify_manifest(base, manifest)
    require(outer == 420 and len(pins) == 409 and len(ancestors) == 6, "Unexpected custody census")
    return dict(outer_manifest_entries=outer, copied_pins=len(pins), ancestor_manifests=ancestors)


def compare_json(expected, actual, policy):
    stats = dict(exact_scalars=0, float_scalars=0, changed_floats=0,
                 max_absolute_error=0.0, max_error_path=None, max_budget_fraction=0.0)

    def walk(a, b, path):
        require(type(a) is type(b), f"JSON type mismatch at {path}")
        if isinstance(a, dict):
            require(a.keys() == b.keys(), f"JSON keys mismatch at {path}")
            for key in sorted(a):
                walk(a[key], b[key], path + "/" + key)
        elif isinstance(a, list):
            require(len(a) == len(b), f"JSON length mismatch at {path}")
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, path + "/" + str(i))
        elif isinstance(a, float):
            require(math.isfinite(a) and math.isfinite(b), f"Nonfinite value at {path}")
            error = abs(a - b)
            budget = policy["float_absolute_tolerance"] + policy["float_relative_tolerance"] * max(abs(a), abs(b))
            require(error <= budget, f"Float mismatch at {path}: {a!r} versus {b!r}, budget {budget}")
            stats["float_scalars"] += 1
            stats["changed_floats"] += int(error != 0)
            stats["max_budget_fraction"] = max(stats["max_budget_fraction"], error / budget)
            if error > stats["max_absolute_error"]:
                stats["max_absolute_error"], stats["max_error_path"] = error, path
        else:
            require(a == b, f"Discrete value mismatch at {path}: {a!r} versus {b!r}")
            stats["exact_scalars"] += 1

    walk(expected, actual, "")
    return stats


def amend_orbital_reference(receipt, records):
    """Independent polynomial projector reference, using sealed records only."""
    import numpy as np

    result = copy.deepcopy(receipt)
    lap = np.zeros((12, 12))
    internal = [(a, b) for a, b, axis in records["outer_graph"]["edges"] if axis == -1 and a < 12 and b < 12]
    require(len(internal) == 30, "Unexpected internal graph")
    for a, b in internal:
        lap[a, a] += 1
        lap[b, b] += 1
        lap[a, b] -= 1
        lap[b, a] -= 1
    values = [0.0, 5 - math.sqrt(5), 6.0, 5 + math.sqrt(5)]
    changes = []
    for i, population in enumerate(result["population"]):
        for j, band in enumerate(population["bands"]):
            lam = values[j + 1]
            projector = np.eye(12)
            for other in values:
                if other != lam:
                    projector = projector @ (lap - other * np.eye(12)) / (lam - other)
            internal_mode = projector[:, 0] / math.sqrt(projector[0, 0])
            old_profile = band["ground"]["profile"]
            profile = np.array(old_profile)
            require(np.all(profile > 0) or np.all(profile < 0), "Ground profile is not strictly one-signed")
            profile = abs(profile)
            spatial = np.repeat(profile, 3) / math.sqrt(3)
            mode = np.kron(spatial, internal_mode)
            observed = []
            for channel in records["runs"][i]["channels"]:
                counts = np.asarray(channel["occupation_counts"])
                require(counts.size == 972 and int(sum(counts)) == 32768, "Unexpected orbital occupation")
                y = math.sqrt(972 / (4 * 6318 * 32768)) * (counts - 32768 / 972)
                observed.append(float(y @ mode))
            changes.append(dict(population=i, band=j, old_profile=old_profile,
                                amended_profile=profile.tolist(),
                                old_readouts=band["actual_ground_mode_record_readouts"],
                                amended_readouts=observed, reference_port=0))
            band["ground"]["profile"] = profile.tolist()
            band["actual_ground_mode_record_readouts"] = observed
    return result, changes


def patch_orbital(path):
    original = path.read_bytes()
    text = original.decode("utf-8")
    for before, after in PATCHES:
        require(text.count(before) == 1, f"Orbital patch anchor mismatch: {before}")
        text = text.replace(before, after)
    path.write_bytes(text.encode("utf-8"))
    return dict(original_sha256=sha_bytes(original), amended_sha256=sha(path),
                reference_port=0, ground_sign="positive first spatial coordinate")


def compare_artifacts(name, expected_dir, actual_dir, policy, amend_reference=True):
    _, _, archive, hash_key, count_key, count = LEAVES[name]
    expected = strict_json((expected_dir / "RECEIPT.json").read_bytes())
    actual = strict_json((actual_dir / "RECEIPT.json").read_bytes())
    status_key = "all_registered_analysis_checks" if name == "boost" else "all_registered_checks"
    require(actual.get(status_key) == "pass", f"Scientific checks failed: {name}")
    if count_key:
        require(type(actual.get(count_key)) is int and actual[count_key] == count, f"Native-call census failed: {name}")
    result = dict(artifacts={})
    filenames = ["RECEIPT.json"] + ([archive] if archive else [])
    for filename in filenames:
        a, b = expected_dir / filename, actual_dir / filename
        result["artifacts"][filename] = dict(original_sha256=sha(a), actual_sha256=sha(b), byte_equal=a.read_bytes() == b.read_bytes())
    if archive:
        a = gzip.decompress((expected_dir / archive).read_bytes())
        b = gzip.decompress((actual_dir / archive).read_bytes())
        require(expected.get(hash_key) == sha_bytes(a), f"Original embedded raw hash failed: {name}")
        require(actual.get(hash_key) == sha_bytes(b), f"Generated embedded raw hash failed: {name}")
        ar, br = strict_json(a), strict_json(b)
        result["records"] = compare_json(ar, br, policy)
        result["artifacts"][archive].update(uncompressed_byte_equal=a == b,
            original_raw_sha256=sha_bytes(a), actual_raw_sha256=sha_bytes(b), embedded_hashes_verified=True)
        if name == "orbital" and amend_reference:
            expected, result["amended_diagnostics"] = amend_orbital_reference(expected, ar)
        # Both bindings have just been independently verified; all other strings stay exact.
        expected[hash_key] = actual[hash_key] = "verified against corresponding decompressed artifact"
    result["receipt"] = compare_json(expected, actual, policy)
    result["native_calls"] = count
    return result


CHILD_BOOTSTRAP = """import json, platform, runpy, sys
if sys.flags.optimize != 0:
    raise RuntimeError('Scientific child requires assertions enabled')
import numpy, zlib
print(json.dumps(dict(python=sys.version, numpy=numpy.__version__, platform=platform.platform(),
                     optimize=sys.flags.optimize, zlib=zlib.ZLIB_RUNTIME_VERSION)), file=sys.stderr, flush=True)
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).parent))
runpy.run_path(sys.argv[1], run_name='__main__')
"""


def run_child(script, output):
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    command = [sys.executable, "-E", "-B", "-X", "utf8", "-c", CHILD_BOOTSTRAP, str(script)]
    started = time.monotonic()
    completed = subprocess.run(command, cwd=script.parent, env=env, capture_output=True, timeout=240)
    output.mkdir(parents=True, exist_ok=True)
    (output / "stdout.txt").write_bytes(completed.stdout)
    (output / "stderr.txt").write_bytes(completed.stderr)
    result = dict(exit_code=completed.returncode, elapsed_seconds=time.monotonic() - started,
                  stdout_sha256=sha_bytes(completed.stdout), stderr_sha256=sha_bytes(completed.stderr),
                  environment=None)
    try:
        result["environment"] = strict_json(completed.stderr.splitlines()[0])
    except (ValueError, ReplayError, IndexError):
        pass
    save(output / "PROCESS.json", result)
    require(completed.returncode == 0, f"Scientific child failed: {script.name}, exit {completed.returncode}; see {output}")
    require(isinstance(result["environment"], dict) and result["environment"].get("optimize") == 0,
            "Scientific child did not confirm assertions enabled")
    return result


def verify_registration():
    registration = strict_json((HERE / "REGISTRATION.json").read_bytes())
    for relative, digest in registration["inputs"].items():
        require(sha(checked_path(HERE, relative)) == digest, f"Amendment input changed since freeze: {relative}")
    return sha(HERE / "REGISTRATION.json")


def main():
    require(sys.flags.optimize == 0, "Optimized Python is unsupported; unset PYTHONOPTIMIZE and remove -O/-OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, help="Use an exact, byte-preserving bundle instead of extracting pinned Git blobs")
    parser.add_argument("--output", type=Path, help="New directory for all artifacts and the replay receipt")
    parser.add_argument("--custody-only", action="store_true")
    args = parser.parse_args()
    require(args.custody_only or args.output is not None, "Full replay requires --output NEW_DIRECTORY")
    registration_sha = verify_registration()
    policy = strict_json((HERE / "POLICY.json").read_bytes())
    if args.output:
        args.output = args.output.resolve()
        args.output.mkdir(parents=True, exist_ok=False)
    report = dict(status="running", started_utc=utc(), policy_sha256=sha(HERE / "POLICY.json"),
                  registration_sha256=registration_sha, source_commit=policy["source_commit"],
                  official="untested", historical_byte_gate="not passed by this amendment",
                  parent_environment=dict(python=sys.version, platform=platform.platform(), optimize=sys.flags.optimize))
    try:
        with tempfile.TemporaryDirectory(prefix="oph_wave_amendment_") as directory:
            work = Path(directory)
            bundle = args.bundle.resolve() if args.bundle else extract_git(work / "subject", policy)
            report["source_mode"] = "provided exact bundle" if args.bundle else "raw pinned Git blobs"
            report["custody"] = custody(bundle, policy)
            if not args.custody_only:
                report["leaves"] = {}
                for name, (relative, script, archive, _, _, _) in LEAVES.items():
                    source, target = bundle / relative, work / "runs" / name
                    shutil.copytree(source, target)
                    if name == "orbital":
                        report["orbital_patch"] = patch_orbital(target / script)
                    # Outputs must be newly generated, never inherited from the sealed copy.
                    for filename in ["RECEIPT.json"] + ([archive] if archive else []):
                        (target / filename).unlink()
                    destination = args.output / name
                    process = run_child(target / script, destination)
                    for filename in ["RECEIPT.json"] + ([archive] if archive else []):
                        shutil.copyfile(target / filename, destination / filename)
                    result = compare_artifacts(name, source, destination, policy)
                    result["process"] = process
                    report["leaves"][name] = result
                    save(args.output / "REPLAY.json", report)
                report["total_native_calls_in_overlapping_lineage"] = sum(x["native_calls"] for x in report["leaves"].values())
                require(report["total_native_calls_in_overlapping_lineage"] == 1009734, "Lineage census failed")
                report["independent_repetitions"] = False
            require(custody(bundle, policy) == report["custody"], "Source custody changed during replay")
            report["source_unchanged"] = True
        require(verify_registration() == registration_sha, "Amendment inputs changed during replay")
        report["status"] = "pass"
    except Exception as error:
        report.update(status="fail", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        report["finished_utc"] = utc()
        if args.output:
            save(args.output / "REPLAY.json", report)
    print(json.dumps(dict(status=report["status"], custody=report["custody"],
                          output=str(args.output) if args.output else None,
                          native_calls=report.get("total_native_calls_in_overlapping_lineage", 0)), indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Replay failed: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
