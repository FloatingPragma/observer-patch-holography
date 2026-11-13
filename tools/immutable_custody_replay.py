"""Run the particle-interface replay with an audited dependency-path cache.

Usage (from a dedicated, immutable checkout)::

    python -m tools.immutable_custody_replay --audit /tmp/write.json -- \
        --write /tmp/particle-interfaces.json
    python -m tools.immutable_custody_replay --audit /tmp/check.json -- \
        --check /tmp/particle-interfaces.json

Only four reviewed, deterministic source-inventory functions are adapted.
Their first calls perform the original discovery and validation; subsequent
calls reuse path names and rehash every live file into a fresh dictionary.
Original discovery runs again before acceptance. The scientific report and
its type-strict comparison execute unchanged. No computed result is cached.

Equivalence requires an immutable source snapshot. The guard compares names
and bytes of tracked/unignored files and all Python/Markdown discovery
candidates, including ignored candidates. Only the declared writer output is
excluded; no cached inventory may depend on it. Fixed retained-record checks
inside native_program.source_bindings execute at first and final discovery.
Use an isolated checkout: this is not an operating-system write lock.

The audit identifies an adapted execution, never an unmodified CLI run.
An independent --check process is required after a successful --write.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import importlib
import json
import os
import platform
import runpy
import subprocess
import sys
import tempfile
import time
from pathlib import Path


OWNER = "model.particle_identification.completion_certificate"
OUTPUT = "model/particle_identification/completion_results.json"
# Updating one of these modules requires reviewing its discovery contract.
TARGET_SOURCE_SHA256 = {
    "model.native_program.receipt":
        "9f507a9f2e14ee628fa3fa420e5b3ba393ccd8d9b5ca475314da0d12e3bc3cfd",
    "model.native_quantum_gauss.packet":
        "3a4b7e8b7e172d366fbe80e4db2bb6862fe00f478edcdef69a9f1ac9dd88d248",
    "model.native_sector_readout.packet":
        "914a3077bc6525e18740c247d9e5d50f63a85af8499b0133edfb68f44c22e4d1",
    "model.native_fair_local_preparation.packet":
        "880103c4330faf0813b82aa44c2391fab9a418ef2fb5183504effa2f137c71ac",
}
DISCOVERY_ROOTS = (
    "model", "tools", "evidence/observer_dynamics_20260925/oph-physics-sim",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def catalogue(root: Path, excluded: frozenset[str]) -> tuple[str, ...]:
    raw = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=root,
    )
    names = set(raw.decode().split("\0")) - {""}
    for directory in DISCOVERY_ROOTS:
        for suffix in ("*.py", "*.md"):
            names.update(str(p.relative_to(root))
                         for p in (root / directory).rglob(suffix) if p.is_file())
    return tuple(sorted(n for n in names - excluded if (root / n).is_file()))


class FrozenSnapshot:
    def __init__(self, root: Path, excluded: frozenset[str] = frozenset()):
        self.root = root.resolve()
        self.excluded = excluded
        self.names = catalogue(self.root, excluded)
        self.ensure_contained(self.names)
        self.bindings = {name: sha(self.root / name) for name in self.names}

    def ensure_contained(self, names):
        for name in names:
            if not (self.root / name).resolve().is_relative_to(self.root):
                raise ValueError("source escapes the immutable repository: " + name)

    def check(self):
        if catalogue(self.root, self.excluded) != self.names:
            raise ValueError("immutable source filename inventory changed")
        self.ensure_contained(self.names)
        for name, expected in self.bindings.items():
            if sha(self.root / name) != expected:
                raise ValueError("immutable source bytes changed: " + name)


class PathCache:
    def __init__(self, module, snapshot: FrozenSnapshot, forbidden: frozenset[str]):
        self.module = module
        self.snapshot = snapshot
        self.forbidden = forbidden
        self.original = module.source_bindings
        self.keys = None
        self.calls = 0
        self.original_discoveries = 0
        self.final_discovery_verified = False
        functools.update_wrapper(self, self.original)

    def fresh(self):
        if Path(self.module.ROOT).resolve() != self.snapshot.root:
            raise ValueError("metadata root changed")
        self.snapshot.ensure_contained(self.keys)
        return {name: sha(self.snapshot.root / name) for name in self.keys}

    def discover(self):
        value = self.original()
        self.original_discoveries += 1
        if (type(value) is not dict or not value
                or any(type(k) is not str or type(v) is not str
                       for k, v in value.items())):
            raise ValueError("metadata discovery is not a nonempty string hash map")
        keys = tuple(sorted(value))
        if (set(keys) & self.forbidden
                or not set(keys) <= self.snapshot.bindings.keys()):
            raise ValueError("discovery reads a mutable or uncaptured input")
        if self.keys is not None and keys != self.keys:
            raise ValueError("original discovered path keys changed")
        self.keys = keys
        if value != self.fresh():
            raise ValueError("original discovery differs from fresh source hashes")
        return value

    def __call__(self):
        self.calls += 1
        return self.discover() if self.keys is None else self.fresh()

    def verify(self):
        if self.keys is not None:
            self.discover()
            self.final_discovery_verified = True

    def evidence(self):
        return {
            "module": self.module.__name__,
            "calls": self.calls,
            "original_discoveries": self.original_discoveries,
            "path_count": len(self.keys or ()),
            "path_keys_sha256": hashlib.sha256(json.dumps(self.keys).encode()).hexdigest(),
            "final_original_discovery_verified": self.final_discovery_verified,
        }


def same_file(left: Path, right: Path) -> bool:
    return left == right or (left.exists() and right.exists()
                             and os.path.samefile(left, right))


def check_expected(path: Path, expected: str):
    if sha(path) != expected:
        raise ValueError("expected receipt changed during strict comparison")


def install(snapshot: FrozenSnapshot, forbidden: frozenset[str], slots: list):
    for name, expected in TARGET_SOURCE_SHA256.items():
        source = snapshot.root / (name.replace(".", "/") + ".py")
        if sha(source) != expected:
            raise ValueError("source-inventory contract requires review: " + name)
        module = importlib.import_module(name)
        if Path(module.__file__).resolve() != source.resolve():
            raise ValueError("metadata module loaded outside selected checkout: " + name)
        slot = PathCache(module, snapshot, forbidden)
        slots.append(slot)
        module.source_bindings = slot


def authenticate_owner(root: Path):
    spec = importlib.util.find_spec(OWNER)
    expected = root / (OWNER.replace(".", "/") + ".py")
    if spec is None or spec.origin is None or Path(spec.origin).resolve() != expected:
        raise ValueError("scientific owner loaded outside selected checkout")


def authenticate_loaded_modules(snapshot: FrozenSnapshot):
    for name, module in tuple(sys.modules.items()):
        if name != "model" and name != "tools" and not name.startswith(("model.", "tools.")):
            continue
        if name == "tools.immutable_custody_replay":
            continue  # The adapter is separately bound by its own source hash.
        source = getattr(module, "__file__", None)
        if source is None:
            continue  # Namespace packages have no executing initializer.
        path = Path(source).resolve()
        if not path.is_relative_to(snapshot.root):
            raise ValueError("repository module loaded outside selected checkout: " + name)
        relative = str(path.relative_to(snapshot.root))
        if snapshot.bindings.get(relative) != sha(path):
            raise ValueError("loaded repository module is not an immutable input: " + name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("original_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    original_args = args.original_args
    if original_args[:1] == ["--"]:
        original_args = original_args[1:]
    if len(original_args) != 2 or original_args[0] not in ("--write", "--check"):
        parser.error("pass -- followed by --write or --check and an explicit output path")
    output = Path(original_args[1]).resolve()
    root, audit, adapter = args.root.resolve(), args.audit.resolve(), Path(__file__).resolve()
    writing = original_args[0] == "--write"
    if writing and output.is_relative_to(root) and output != root / OUTPUT:
        parser.error("only the declared receipt may be written inside the repository")
    if (audit.is_relative_to(root) or same_file(audit, output)
            or same_file(audit, adapter) or same_file(output, adapter)):
        parser.error("audit must be external; audit, output and adapter must be distinct")
    original_args[1] = str(output)
    excluded = frozenset({OUTPUT}) if writing and output == root / OUTPUT else frozenset()
    os.chdir(root)
    sys.path.insert(0, str(root))
    started = time.monotonic()
    snapshot = FrozenSnapshot(root, excluded)
    slots = []
    report = {
        "schema": "oph.adapted_metadata_path_replay.v1",
        "mode": "audited_path_inventory_cache",
        "root": str(root), "owner": OWNER, "original_arguments": original_args,
        "adapter_sha256": sha(adapter), "python": sys.version,
        "platform": platform.platform(), "input_sha256": snapshot.bindings,
        "scientific_result_cache": False, "comparison_function_changed": False,
        "unmodified_cli_run": False,
    }
    staged = None
    try:
        if not writing:
            report["expected_output_sha256_before"] = sha(output)
        if writing and any(same_file(output, root / n) for n in snapshot.names):
            raise ValueError("writer output aliases an immutable source")
        if any(same_file(audit, root / n) for n in snapshot.names):
            raise ValueError("audit output aliases an immutable source")
        authenticate_owner(root)
        authenticate_loaded_modules(snapshot)
        install(snapshot, frozenset({OUTPUT}), slots)
        executed_args = list(original_args)
        if writing and output == root / OUTPUT:
            # A late audit failure must preserve the canonical receipt.
            handle, temporary = tempfile.mkstemp(prefix="oph-custody-replay-", suffix=".json")
            os.close(handle)
            staged = Path(temporary)
            executed_args[1] = str(staged)
        report["executed_original_arguments"] = executed_args
        sys.argv = [OWNER, *executed_args]
        runpy.run_module(OWNER, run_name="__main__")
        if not writing:
            check_expected(output, report["expected_output_sha256_before"])
        for slot in slots:
            slot.verify()
        authenticate_loaded_modules(snapshot)
        snapshot.check()
        if sha(adapter) != report["adapter_sha256"]:
            raise ValueError("adapter source changed during replay")
        if staged is not None:
            # A same-directory staging file makes the final replacement atomic
            # even when the external temporary directory is on another device.
            descriptor, local_name = tempfile.mkstemp(prefix=".oph-custody-", dir=output.parent)
            local = Path(local_name)
            try:
                with os.fdopen(descriptor, "wb") as stream:
                    stream.write(staged.read_bytes())
                local.chmod((output.stat().st_mode & 0o777) if output.exists() else 0o644)
                os.replace(local, output)
            finally:
                local.unlink(missing_ok=True)
        report.update(
            verdict=("PASS_ADAPTED_ORIGINAL_SCIENCE_WRITER" if writing else
                     "PASS_ADAPTED_ORIGINAL_SCIENCE_AND_STRICT_COMPARISON"),
            output_sha256=sha(output), output=str(output), original_report_executed=True,
            original_strict_comparison_executed=not writing,
            immutable_source_bytes_and_filenames_preserved=True,
        )
    except BaseException as error:
        report.update(verdict="FAIL", error=repr(error))
        raise
    finally:
        for slot in slots:
            slot.module.source_bindings = slot.original
        if staged is not None:
            staged.unlink(missing_ok=True)
        report["metadata_caches"] = [slot.evidence() for slot in slots]
        report["elapsed_seconds"] = time.monotonic() - started
        # Failed alias validation must not overwrite an immutable input.
        if not any(same_file(audit, root / n) for n in snapshot.names):
            audit.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({key: report.get(key) for key in
                          ("verdict", "mode", "elapsed_seconds", "output_sha256")}), flush=True)


if __name__ == "__main__":
    main()
