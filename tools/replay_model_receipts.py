#!/usr/bin/env python3
"""Plan or replay live model receipts after source-only maintenance.

Default is read-only and plans all registered live receipts. --execute supports
only the explicitly reviewed, source-pinned path-output writers below; it refuses
the complete operation before any writes if an affected writer is unsupported.
Registered checker commands alone do not establish a safe writer interface.
Supported producers and cold checkers run in dependency order; no receipt digest
is patched.
Frozen historical captures are excluded because they are not in REPLAYS.
Non-provenance differences stop execution for review, retaining original bytes
and command logs. Source/schema changes must be reviewed before this command;
this maintenance operation does not promote scientific claims.
Every invocation requires a new or empty output directory, so earlier audit
logs and original-byte backups cannot be replaced by a later invocation.
Run from the repository with its pinned Python environment.
"""
from __future__ import annotations
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading

REPOSITORY_FIELDS = {"sources", "source_sha256", "source_bindings", "input_sha256",
                     "source_binding_sha256", "source_control_bindings", "source_pins"}
PACKAGE_LEVELS = {"implementation_sha256": 0, "source_interface_hashes": 1}
ENDPOINT_PINS = "root.retrospective_endpoint_audit.source.pins"

# Reviewed interfaces: both --write PATH and --check PATH honor the supplied
# path; writers do not replace canonical outputs. A parser/source change needs
# renewed review. Other REPLAYS entries may have pathless or different writers.
SUPPORTED_WRITERS = {
    ('model.laboratory_clock', None): ('model/laboratory_clock/__main__.py',
        '4031dbfa8ab56b51e5de87f6371c9399ba864eb035967a5235637b906a9f37e7'),
    ('model.native_clock.host_clock', None): ('model/native_clock/host_clock.py',
        '0c803a72b0a0eede334fb77ee37c25d0e48195f6e66c2055b844163d8d62aa5f'),
    ('model.native_clock.clock_bridge', None): ('model/native_clock/clock_bridge.py',
        '466f9b99a3d105ccaa0c6981144b3087bafebcb9a5916af049c1ccbdede11be5'),
    ('model.interacting_field', None): ('model/interacting_field/__main__.py',
        '0e388f9f0deda353c3990af4157a8401038f5d7a1328c3d086630400112ef714'),
}


def writer_problem(root, row):
    module, command, _ = row
    contract = SUPPORTED_WRITERS.get((module, command))
    if contract is None:
        return "no reviewed path-output writer contract"
    path, expected = contract
    source = root / path
    if not source.is_file() or digest(source) != expected:
        return "reviewed writer source changed: " + path
    return None


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_maps(value, location="root"):
    if isinstance(value, dict):
        for key, child in value.items():
            here = f"{location}.{key}"
            recognized = key in REPOSITORY_FIELDS or key in PACKAGE_LEVELS or here == ENDPOINT_PINS
            if recognized:
                yield key, child, here
            yield from source_maps(child, here)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from source_maps(child, f"{location}[{index}]")


def resolve_map(root, receipt, key, child, here):
    local_only = key != "sources"
    if local_only and (not isinstance(child, dict) or not child):
        raise ValueError(f"Empty or malformed local source map: {receipt} at {here}")
    if not isinstance(child, dict):
        return
    base = (root / receipt).parents[PACKAGE_LEVELS[key]] if key in PACKAGE_LEVELS else root
    for relative, expected in child.items():
        if not isinstance(relative, str) or "://" in relative:
            if local_only:
                raise ValueError(f"Invalid local source path: {receipt} at {here}")
            continue
        path = Path(relative)
        if path.is_absolute():
            if local_only:
                raise ValueError(f"Absolute local source path: {receipt} at {here}[{relative}]")
            continue
        try:
            target = (base / path).resolve()
            local_file = target.is_relative_to(root) and target.is_file()
        except (OSError, ValueError):
            local_file = False
        if not local_file:
            if local_only:
                raise ValueError(f"Missing or escaping local source path: {receipt} at {here}[{relative}]")
            continue
        yield relative, expected, target


def load_receipt(root, receipt, strict_json):
    raw = (root / receipt).read_bytes()
    if receipt.endswith(".gz"):
        raw = gzip.decompress(raw)
    return strict_json.loads(raw.decode("utf-8"))


def owned_outputs(root, receipt, value):
    """Identify this producer's receipt and declared/conventional archives.

    Only top-level output metadata is traversed. An imported receipt nested
    inside a proof does not make its archives outputs of the importing job.
    """
    own = (root / receipt).resolve()
    paths = {own}
    stem = own.name.removesuffix(".gz").removesuffix(".json")
    paired = own.with_name(stem + ".npz")
    if paired.is_file():
        paths.add(paired)

    def collect(node):
        if isinstance(node, dict):
            for key in ("path", "filename"):
                if key not in node:
                    continue
                name = node[key]
                if not isinstance(name, str) or not name.endswith((".npz", ".json.gz")):
                    raise ValueError(f"Invalid archive path for {receipt}: {name!r}")
                relative = Path(name)
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError(f"Escaping archive path for {receipt}: {name!r}")
                path = ((root if relative.parts[0] == "model" else own.parent) / relative).resolve()
                if not path.is_relative_to(root) or not path.is_file():
                    raise ValueError(f"Missing or escaping archive for {receipt}: {name!r}")
                paths.add(path)
            for key, child in node.items():
                if key not in ("path", "filename"):
                    collect(child)
        elif isinstance(node, list):
            for child in node:
                collect(child)

    for key in ("artifact", "complete_record_archive", "array_archive"):
        if key in value:
            collect(value[key])
    return tuple(sorted(paths))


def output_owners(root, values):
    owners, outputs = {}, {}
    for receipt, value in values.items():
        outputs[receipt] = owned_outputs(root, receipt, value)
        for path in outputs[receipt]:
            if path in owners:
                raise ValueError(f"Multiple producers own {path}: {owners[path]}, {receipt}")
            owners[path] = receipt
    return owners, outputs


def bindings(root, receipt, value):
    for key, child, here in source_maps(value):
        for relative, expected, target in resolve_map(root, receipt, key, child, here):
            yield here, relative, expected, target


def science_payload(root, receipt, value):
    # Preserve all paths, map structure, external metadata and non-source data.
    # Only values which the authoritative source-map checker treats as local
    # provenance are masked. New/removed keys or any other change remain visible.
    value = json.loads(json.dumps(value, allow_nan=False))
    for key, child, here in source_maps(value):
        for relative, expected, target in list(resolve_map(root, receipt, key, child, here)):
            child[relative] = "<LOCAL_SOURCE_SHA256>"
    return value


def differences(left, right, path="root", limit=50):
    result = []
    if type(left) is not type(right):
        return [{"path": path, "old": repr(left)[:240], "new": repr(right)[:240]}]
    if isinstance(left, dict):
        for key in sorted(left.keys() | right.keys()):
            if key not in left or key not in right:
                result.append({"path": f"{path}.{key}", "old": repr(left.get(key))[:240],
                               "new": repr(right.get(key))[:240]})
            else:
                result.extend(differences(left[key], right[key], f"{path}.{key}", limit))
            if len(result) >= limit:
                break
    elif isinstance(left, list):
        if len(left) != len(right):
            result.append({"path": path, "old_length": len(left), "new_length": len(right)})
        for i, (old, new) in enumerate(zip(left, right)):
            result.extend(differences(old, new, f"{path}[{i}]", limit))
            if len(result) >= limit:
                break
    elif left != right:
        result.append({"path": path, "old": repr(left)[:240], "new": repr(right)[:240]})
    return result[:limit]


def topological(graph):
    pending = {key: set(value) for key, value in graph.items()}
    order = []
    while pending:
        ready = sorted(key for key, deps in pending.items() if not deps)
        if not ready:
            raise ValueError("Live receipt dependency cycle: " + repr(pending))
        order.extend(ready)
        for key in ready:
            pending.pop(key)
        for deps in pending.values():
            deps.difference_update(ready)
    return order


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--jobs", type=int, default=1,
                        help="maximum simultaneous independent receipts (1 through 32)")
    args = parser.parse_args()
    if not 1 <= args.jobs <= 32:
        parser.error("--jobs must be between 1 and 32")
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("--output must be new or empty; preserve previous audit logs and backups")
    root = args.root.resolve()
    sys.path.insert(0, str(root))
    from model import verify
    from tools import strict_json
    if verify.ROOT.resolve() != root:
        raise ValueError("Imported model does not belong to requested root")
    replay = {row[2]: row for row in verify.REPLAYS}
    if len(replay) != len(verify.REPLAYS):
        raise ValueError("Duplicate registered receipt")
    values = {name: load_receipt(root, name, strict_json) for name in replay}
    by_path, outputs = output_owners(root, values)
    graph = {}
    initial = {}
    digests = {}
    for name, value in values.items():
        graph[name] = set()
        for here, relative, expected, target in bindings(root, name, value):
            if target in by_path and by_path[target] != name:
                graph[name].add(by_path[target])
            if target not in digests:
                digests[target] = digest(target)
            if expected != digests[target]:
                initial.setdefault(name, []).append({"map": here, "path": relative,
                    "target": target.relative_to(root).as_posix(), "expected": expected,
                    "actual": digests[target]})
    order = topological(graph)
    affected = set(initial)
    for name in order:
        if graph[name] & affected:
            affected.add(name)
    report = {"root": str(root), "authoritative_verifier_sha256": digest(root / "model/verify.py"),
              "execute": args.execute, "jobs": args.jobs, "registered_count": len(replay),
              "initial_stale": initial,
              "unsupported_affected_writers": {name: problem for name in order
                  if name in affected and (problem := writer_problem(root, replay[name]))},
              "affected_topological_order": [name for name in order if name in affected],
              "producer_outputs": {name: [p.relative_to(root).as_posix() for p in paths]
                                   for name, paths in outputs.items()},
              "dependency_graph": {name: sorted(deps) for name, deps in graph.items()},
              "results": []}
    args.output.mkdir(parents=True, exist_ok=True)
    report_path = args.output / "report.json"
    report_lock = threading.RLock()
    def save():
        with report_lock:
            report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    save()
    print(f"Registered={len(replay)}, directly stale={len(initial)}, affected={len(affected)}", flush=True)
    if not args.execute:
        print(f"Read-only plan: {report_path}", flush=True)
        return
    if report["unsupported_affected_writers"]:
        raise SystemExit("No producers started: unsupported or changed writer interfaces; "
                         "inspect unsupported_affected_writers in " + str(report_path))
    environment = dict(os.environ)
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[var] = "1"
    def execute_one(name):
        module, command, _ = replay[name]
        prefix = [sys.executable, "-W", "error", "-m", module]
        if command is not None:
            prefix.append(command)
        before = load_receipt(root, name, strict_json)
        backed_up = {}
        for path in outputs[name]:
            relative = path.relative_to(root)
            backup = args.output / "before" / relative
            backup.parent.mkdir(parents=True, exist_ok=True)
            backup.write_bytes(path.read_bytes())
            backed_up[relative.as_posix()] = digest(path)
        result = {"receipt": name, "before_sha256": digest(root / name),
                  "backed_up_outputs": backed_up}
        with report_lock:
            report["results"].append(result)
            save()
        print(f"Replaying {name}", flush=True)
        candidate_root = args.output.resolve() / "candidate"
        candidate = candidate_root / name
        candidate.parent.mkdir(parents=True, exist_ok=True)
        for operation in ("write", "check"):
            # A number of evidence producers correctly refuse to overwrite an
            # existing path. Give every producer a fresh candidate path, and
            # do not replace canonical outputs until its cold checker passes.
            argv = prefix + ["--" + operation, str(candidate)]
            log = args.output / (module + ("_" + command.replace("-", "") if command else "") + "." + operation + ".log")
            with report_lock:
                result[operation] = {"argv": argv, "log": str(log)}
                save()
            with log.open("w") as stream:
                completed = subprocess.run(argv, cwd=root, env=environment, stdout=stream, stderr=subprocess.STDOUT)
            with report_lock:
                result[operation]["returncode"] = completed.returncode
                save()
            if completed.returncode:
                raise RuntimeError(f"{name} {operation} failed; inspect {log}")
            if operation == "write":
                after = load_receipt(candidate_root, name, strict_json)
                diff = differences(science_payload(root, name, before), science_payload(root, name, after))
                with report_lock:
                    result["after_sha256"] = digest(candidate)
                    result["non_provenance_differences"] = diff
                    save()
                if diff:
                    raise RuntimeError(f"{name} changed non-provenance payload; backups={args.output / 'before'}")
        candidate_outputs = owned_outputs(candidate_root, name, after)
        relative_outputs = {p.relative_to(candidate_root).as_posix(): p for p in candidate_outputs}
        if relative_outputs.keys() != backed_up.keys():
            raise RuntimeError(f"{name} changed its output catalogue; review candidate outputs")
        if any(digest(root / relative) != expected for relative, expected in backed_up.items()):
            raise RuntimeError(f"{name} canonical output changed during replay")
        # Commit archives first and the receipt last. Dependent DAG nodes are
        # released only after this entire operation returns successfully.
        for relative in sorted(relative_outputs, key=lambda p: p == name):
            destination = root / relative
            with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".replay-",
                                             delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(relative_outputs[relative].read_bytes())
            try:
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
        return name

    pending = set(affected)
    complete = set(replay) - affected
    running = {}
    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        while pending or running:
            if not failures:
                ready = [name for name in order if name in pending and graph[name] <= complete]
                for name in ready[:args.jobs - len(running)]:
                    pending.remove(name)
                    running[pool.submit(execute_one, name)] = name
            if not running:
                if pending and not failures:
                    failures.append({"receipt": None, "error": "No ready DAG nodes remain"})
                break
            finished, _ = wait(running, return_when=FIRST_COMPLETED)
            for future in finished:
                name = running.pop(future)
                try:
                    future.result()
                except Exception as error:
                    failures.append({"receipt": name, "error": str(error)})
                else:
                    complete.add(name)
            with report_lock:
                report["failures"] = list(failures)
                report["completed"] = sorted(complete & affected)
                report["not_scheduled"] = sorted(pending)
                report["running"] = sorted(running.values())
                save()
        # Exiting the executor waits for every already-running producer. After
        # the first failure no further node, dependent or independent, starts.
    if failures:
        raise SystemExit("Stopped after awaiting running nodes: " + repr(failures))
    # Fresh authoritative validation, not the local copy of discovery logic.
    verify.verify_receipt_source_bindings()
    report["authoritative_final_source_check"] = "PASS"
    save()
    print(f"All {len(affected)} affected receipts replayed and checked: {report_path}", flush=True)


if __name__ == "__main__":
    main()
