"""Run OW-30 from the author's byte-preserved modules shipped inside OpenWave."""

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
SCRIPTS = Path("openwave/xperiments/m10_oph/research/scripts/m10_5")
MATERIAL = Path("model/candidate_laws/interval_echo/openwave/row_30_collective")
TASK = "openwave/xperiments/m10_oph/research/tasks/m10_5_task_details.md"


def git(repo, *args):
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True
    ).stdout


def mounted_source(platform, source, revision):
    """Check every physical input and every imported project module against Git blobs."""
    if HERE != platform / SCRIPTS:
        raise ValueError("Run the copy in the OpenWave experiment directory")
    layout = json.loads((HERE / "layout.json").read_text())
    for target, original in layout.items():
        path = platform / target
        if path.read_bytes() != git(source, "show", f"{revision}:{original}"):
            raise ValueError("Platform bytes differ from source release: " + target)
    vendor = HERE / "vendor"
    expected = {p for p in layout if p.startswith(f"{SCRIPTS}/vendor/")}
    runtime_suffixes = {".py", ".pyc", ".so", ".pyd", ".dylib", ".npz", ".json"}
    actual = {
        p.relative_to(platform).as_posix() for p in vendor.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
        and (p.suffix in runtime_suffixes or p.relative_to(platform).as_posix() in expected)
    }
    if actual != expected:
        raise ValueError("Unexpected or missing platform input: " + str(sorted(actual ^ expected)))
    return vendor


def imports_used(vendor):
    origins = {}
    for name, module in sys.modules.items():
        if name == "model" or name.startswith("model."):
            file = getattr(module, "__file__", None)
            if file:
                path = Path(file).resolve()
                origins[name] = path.relative_to(vendor).as_posix()
            elif any(not Path(p).resolve().is_relative_to(vendor) for p in module.__path__):
                raise ValueError("External namespace package: " + name)
    return origins


def send_to_detector(protocol, platform, source, revision, vendor, args):
    # Git custody is checked on the official projects. No author-tree module is
    # added to the Python search path during that check.
    protocol.HERE = source / MATERIAL / "registration"
    protocol.ROOT = source
    anchor = protocol.official_anchor(platform)
    if not anchor:
        raise SystemExit("Prospective task absent from official OpenWave main; no settings derived.")
    if protocol.git(platform, "show", anchor + ":" + protocol.SOURCE_POINTER).decode().strip() != revision:
        raise ValueError("Registration refers to a different source release")
    # The owner runs the port from an OpenWave checkout descended from the
    # registered task. The code and resulting records may land together later.
    git(platform, "merge-base", "--is-ancestor", anchor, "HEAD")
    for rel in (protocol.TASK, protocol.SOURCE_POINTER):
        if (platform / rel).read_bytes() != git(platform, "show", anchor + ":" + rel):
            raise ValueError("Platform task or source pointer differs from registration: " + rel)
    # run_bank obtains its generator, inputs, meters and profile from this root.
    protocol.ROOT = vendor
    protocol.HERE = vendor / MATERIAL / "registration"
    if args.command == "acquire":
        return protocol.run_bank(anchor, args.release, revision)
    with tempfile.TemporaryDirectory() as staging:
        trial = Path(staging) / "run"
        passed = protocol.run_bank(anchor, trial, revision)
        protocol.compare_release(args.release, trial)
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("pilot", "acquire", "replay"))
    parser.add_argument("--source", required=True, type=Path, help="Git checkout for independent byte verification")
    parser.add_argument("--release", required=True, type=Path, help="New directory for a pilot or acquisition; saved acquisition for a replay")
    args = parser.parse_args()
    if not sys.flags.isolated or not sys.dont_write_bytecode:
        parser.error("use python -I -B to isolate imports and preserve the tested payload")
    platform = Path(git(HERE, "rev-parse", "--show-toplevel").decode().strip()).resolve()
    source = args.source.resolve()
    # The same pointer names both the pilot snapshot and the formal one. Only
    # the formal run additionally requires custody on the two public main refs.
    pointer = (platform / TASK).with_suffix(".source").read_bytes()
    revision = pointer.decode().strip()
    if len(revision) != 40 or any(ch not in "0123456789abcdef" for ch in revision) or pointer != (revision + "\n").encode():
        raise ValueError("Malformed RER Git revision in the platform source pointer")
    vendor = mounted_source(platform, source, revision)
    # -B suppresses writes, but Python can still read an existing bytecode
    # cache. A fresh empty prefix makes those ignored checkout caches inert.
    cache = tempfile.TemporaryDirectory()
    sys.pycache_prefix = cache.name
    sys.path.insert(0, str(vendor))
    base = "model.candidate_laws.interval_echo.openwave.row_30_collective"
    if args.command == "pilot":
        # The published local 927016 holdout is the only seed this command uses.
        holdout = importlib.import_module(base + ".holdout")
        assert holdout.SEED == 927016
        args.release.mkdir(parents=True, exist_ok=False)
        receipt = holdout.run()
        (args.release / "pilot.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        expected = json.loads((vendor / MATERIAL / "holdout.json").read_text())
        holdout.verify(receipt, expected)
        passed = all(b["score"]["operational_pass"] and b["score"]["shape_pass"] for b in receipt["boards"])
    else:
        protocol = importlib.import_module(base + ".registration.verify")
        passed = send_to_detector(protocol, platform, source, revision, vendor, args)
    origins = imports_used(vendor)
    if not any(name.endswith(".contrast_memory_spherical_v1.law") for name in origins):
        raise ValueError("The native seam generator was not loaded")
    custody = {
        "command": args.command, "source_commit": revision,
        "platform_head": git(platform, "rev-parse", "HEAD").decode().strip(),
        "execution_root": SCRIPTS.as_posix() + "/vendor",
        "task_sha256": hashlib.sha256((platform / TASK).read_bytes()).hexdigest(),
        "imported_model_files": origins, "passed": passed,
    }
    if args.command == "pilot":
        (args.release / "execution.json").write_text(json.dumps(custody, indent=2) + "\n")
    print(json.dumps(custody, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
