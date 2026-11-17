#!/usr/bin/env python3
"""Plan or replay live model receipts after source-only maintenance.

Default is read-only and plans all registered live receipts (model/verify.py
REPLAYS). With --discover it also plans the unregistered receipts whose owners
tools/audit_receipt_sources.py finds from literal OUTPUT declarations.
--execute runs the producer and the cold checker of every affected receipt in
dependency order; no receipt digest is patched. Each producer's write and
check interface is read from its own --help: a path-taking --write or
--output writes a fresh candidate, which replaces the canonical outputs only
after its cold checker passes; a flag-only --write writes in place, and the
backed-up outputs are restored when the payload changes or the checker fails.

Provenance is exactly the custody pins: the recognized source maps below and
any other {path: sha256} map or {path, sha256} record list, whose prose and
test entries (tools/pin_scope.py) may disappear, plus any SHA-256 value that
changes to the current digest of a repository file other than the
producer's own outputs. Every other difference stops the operation for
review, retaining original bytes and command logs. Frozen historical captures
are excluded because they are not in REPLAYS.
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
import re
import subprocess
import sys
import tempfile
import threading
import time
from subprocess import run as _git_run

try:
    from tools.pin_scope import exempt
except ImportError:  # executed as tools/replay_model_receipts.py or imported from tools/
    from pin_scope import exempt

REPOSITORY_FIELDS = {"sources", "source_sha256", "source_bindings", "input_sha256",
                     "source_binding_sha256", "source_control_bindings", "source_pins"}
PACKAGE_LEVELS = {"implementation_sha256": 0, "source_interface_hashes": 1}
ENDPOINT_PINS = "root.retrospective_endpoint_audit.source.pins"
HEX = re.compile(r"[0-9a-f]{64}\Z")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Registered rows whose checker is a router: the receipt's producer.
def writer_row(row):
    module, command, receipt = row
    if module == "model.replay_native_fixed":
        from model.replay_native_fixed import OWNERS
        return OWNERS[command], None
    if module == "model.common_source.replay":
        return f"model.common_source.{command}.run", None
    return module, command


def producer_help(root, module, command):
    argv = [sys.executable, "-m", module] + ([command] if command is not None else []) + ["--help"]
    completed = subprocess.run(argv, cwd=root, capture_output=True, text=True)
    if completed.returncode:
        return None, "--help failed: " + (completed.stdout + completed.stderr)[-400:]
    return completed.stdout, None


def writer_interface(root, row):
    """Read a producer's write/check interface from its own --help.

    write: ("path", flag) writes a fresh candidate at a given path,
    ("inplace", "--write") rewrites the canonical outputs, and ("bare", None)
    writes the canonical outputs when run without arguments (they are moved
    aside first, since such producers refuse to overwrite).
    """
    module, command = writer_row(row)
    text, problem = producer_help(root, module, command)
    if problem:
        return None, problem
    write_flag = re.search(r"--write\b(?! \[?[A-Z_])", text) is not None
    for flag in ("--write", "--output", "--out"):
        if re.search(re.escape(flag) + r" \[?[A-Z_]+", text):
            # A flag-only --write that directs its output with --output.
            write = ("path", "--write --output" if flag == "--output" and write_flag else flag)
            break
    else:
        write = ("inplace", "--write") if write_flag else ("bare", None)
    if (module, command) != row[:2]:
        text, problem = producer_help(root, row[0], row[1])
        if problem:
            return None, problem
    if re.search(r"--check \[?[A-Z_]+", text):
        check = "path"
    elif re.search(r"--check\b", text):
        check = "flag"
    elif re.search(r"--verify\b", text):
        check = "verify"
    else:
        # No separate checker: the exact payload comparison is the check.
        check = None
    return {"write": write, "check": check}, None


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
        if exempt(relative):
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


def output_owners(root, values, optional=()):
    """Map every output file to its producer; an optional producer that would
    share an output with another one is dropped and returned."""
    owners, outputs, dropped = {}, {}, []
    for receipt in sorted(values, key=lambda name: name in optional):
        try:
            paths = owned_outputs(root, receipt, values[receipt])
        except ValueError:
            if receipt not in optional:
                raise
            dropped.append(receipt)
            continue
        clash = [path for path in paths if path in owners]
        if clash and receipt in optional:
            dropped.append(receipt)
            continue
        if clash:
            raise ValueError(f"Multiple producers own {clash[0]}: {owners[clash[0]]}, {receipt}")
        outputs[receipt] = paths
        for path in paths:
            owners[path] = receipt
    return owners, outputs, dropped


def bindings(root, receipt, value):
    for key, child, here in source_maps(value):
        for relative, expected, target in resolve_map(root, receipt, key, child, here):
            yield here, relative, expected, target


def path_like(key):
    return isinstance(key, str) and "://" not in key and (
        "/" in key or re.search(r"\.[A-Za-z0-9]{1,6}\Z", key) is not None)


def pin_map(node):
    """A {path: sha256} map, whatever its field name."""
    return (isinstance(node, dict) and bool(node) and any(path_like(k) for k in node)
            and all(isinstance(v, str) and HEX.match(v) for v in node.values()))


def pin_records(node):
    """A [{path, sha256, ...}] record list."""
    return (isinstance(node, list) and bool(node)
            and all(isinstance(r, dict) and isinstance(r.get("path"), str)
                    and isinstance(r.get("sha256"), str) and HEX.match(r["sha256"]) for r in node))


def hex_values(value, found=None):
    found = set() if found is None else found
    if isinstance(value, dict):
        for child in value.values():
            hex_values(child, found)
    elif isinstance(value, list):
        for child in value:
            hex_values(child, found)
    elif isinstance(value, str) and HEX.match(value):
        found.add(value)
    return found


def science_payload(root, receipt, value):
    """The receipt without its custody pins.

    Removed: every local entry of a recognized source map (external metadata
    in mixed ``sources`` maps stays), every {path: sha256} map under any field
    name, and the digests and prose/test rows of {path, sha256} record lists.
    Which files a map pins is custody; every other value remains visible.
    """
    value = json.loads(json.dumps(value, allow_nan=False))

    def strip(node, location="root"):
        if isinstance(node, dict):
            for key, child in list(node.items()):
                here = f"{location}.{key}"
                recognized = key in REPOSITORY_FIELDS or key in PACKAGE_LEVELS or here == ENDPOINT_PINS
                if recognized and isinstance(child, dict):
                    base = (root / receipt).parents[PACKAGE_LEVELS[key]] if key in PACKAGE_LEVELS else root
                    for relative in list(child):
                        local = False
                        if isinstance(relative, str) and "://" not in relative and not Path(relative).is_absolute():
                            try:
                                target = (base / relative).resolve()
                                local = target.is_relative_to(root) and target.is_file()
                            except (OSError, ValueError):
                                local = False
                        hexed = isinstance(child[relative], str) and HEX.match(child[relative])
                        if local or exempt(relative) or (key != "sources" and hexed):
                            del child[relative]
                    node[key] = "<SOURCE_MAP>" if key != "sources" else child
                    if key == "sources":
                        strip(child, here)
                elif pin_map(child):
                    node[key] = "<PIN_MAP>"
                elif pin_records(child):
                    node[key] = [dict(row, sha256="<SHA256>") for row in child if not exempt(row["path"])]
                else:
                    strip(child, here)
        elif isinstance(node, list):
            for index, child in enumerate(node):
                if pin_map(child):
                    node[index] = "<PIN_MAP>"
                else:
                    strip(child, f"{location}[{index}]")
    if pin_map(value):
        return "<PIN_MAP>"
    strip(value)
    return value


def differences(left, right, path="root", limit=50, file_digest=None, float_tolerance=None):
    """Differences of two payloads. A SHA-256 value that changes to a value for
    which ``file_digest`` holds (the current digest of a repository file) is a
    custody pin and no difference. With ``float_tolerance`` (relative, with an
    absolute floor of 1e-12), floats that agree to it are no difference; this
    classifies differences for a report and never admits a changed payload."""
    result = []
    if (float_tolerance is not None and isinstance(left, float) and isinstance(right, float)
            and abs(left - right) <= max(1e-12, float_tolerance * max(abs(left), abs(right)))):
        return result
    if type(left) is not type(right):
        return [{"path": path, "old": repr(left)[:240], "new": repr(right)[:240]}]
    if isinstance(left, dict):
        for key in sorted(left.keys() | right.keys()):
            if key not in left or key not in right:
                result.append({"path": f"{path}.{key}", "old": repr(left.get(key))[:240],
                               "new": repr(right.get(key))[:240]})
            else:
                result.extend(differences(left[key], right[key], f"{path}.{key}", limit, file_digest, float_tolerance))
            if len(result) >= limit:
                break
    elif isinstance(left, list):
        if len(left) != len(right):
            result.append({"path": path, "old_length": len(left), "new_length": len(right)})
        for i, (old, new) in enumerate(zip(left, right)):
            result.extend(differences(old, new, f"{path}[{i}]", limit, file_digest, float_tolerance))
            if len(result) >= limit:
                break
    elif left != right:
        pin = (file_digest is not None and isinstance(left, str) and isinstance(right, str)
               and HEX.match(left) and HEX.match(right) and file_digest(right))
        if not pin:
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


def touched_since(root, start_ns, skip=()):
    """Files under model/ modified at or after ``start_ns``."""
    found = set()
    for directory, names, files in os.walk(root / "model"):
        names[:] = [n for n in names if n != "__pycache__"]
        for file_name in files:
            path = Path(directory) / file_name
            try:
                if path.stat().st_mtime_ns >= start_ns and path.resolve() not in skip:
                    found.add(path.resolve())
            except OSError:
                pass
    return found


def restore_from_git(root, paths):
    """Return side-effect files to their committed bytes; remove new ones."""
    for path in sorted(paths):
        relative = path.relative_to(root).as_posix()
        tracked = _git_run(["git", "cat-file", "-e", "HEAD:" + relative], cwd=root,
                           capture_output=True).returncode == 0
        if tracked:
            _git_run(["git", "checkout", "HEAD", "--", relative], cwd=root, capture_output=True)
        else:
            path.unlink(missing_ok=True)


def repository_digests(root):
    """Current SHA-256 of every repository file (tracked files when in Git)."""
    try:
        listing = subprocess.run(["git", "ls-files", "-z"], cwd=root, capture_output=True, check=True).stdout
        names = [name for name in listing.decode().split("\0") if name]
    except (OSError, subprocess.CalledProcessError):
        names = [p.relative_to(root).as_posix() for p in root.rglob("*")
                 if p.is_file() and ".git" not in p.relative_to(root).parts]
    digests = {}
    for name in names:
        path = root / name
        if path.is_file() and not path.is_symlink():
            digests[path.resolve()] = digest(path)
    return digests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--discover", action="store_true",
                        help="also plan unregistered receipts with a literal OUTPUT owner")
    parser.add_argument("--keep-going", action="store_true",
                        help="after a failure, still run receipts that do not depend on it")
    parser.add_argument("--stop-after", type=float, metavar="MINUTES",
                        help="start no producer after this many minutes; running ones finish")
    parser.add_argument("--float-tolerance", type=float, metavar="REL",
                        help="report only: also record which payload differences remain "
                             "beyond this relative float tolerance; any payload difference "
                             "still stops the receipt and keeps its original bytes")
    parser.add_argument("--jobs", type=int, default=1,
                        help="maximum simultaneous independent receipts (1 through 32)")
    args = parser.parse_args()
    if not 1 <= args.jobs <= 32:
        parser.error("--jobs must be between 1 and 32")
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("--output must be new or empty; preserve previous audit logs and backups")
    started = time.monotonic()
    root = args.root.resolve()
    sys.path.insert(0, str(root))
    from model import verify
    from tools import strict_json
    if verify.ROOT.resolve() != root:
        raise ValueError("Imported model does not belong to requested root")
    replay = {row[2]: row for row in verify.REPLAYS}
    if len(replay) != len(verify.REPLAYS):
        raise ValueError("Duplicate registered receipt")
    discovered = set()
    if args.discover:
        from tools.audit_receipt_sources import discover_owners
        owners, _ = discover_owners(root)
        for receipt, owner in owners.items():
            if receipt not in replay:
                replay[receipt] = (owner.module, owner.command, receipt)
                discovered.add(receipt)
        # A package's results.json next to its __main__.py is that package's
        # conventional output (python -m model.<package> --write/--check).
        for main_file in sorted((root / "model").rglob("__main__.py")):
            relative = main_file.parent.relative_to(root)
            if {"historical", "predecessor"} & set(relative.parts):
                continue
            receipt = (relative / "results.json").as_posix()
            if receipt not in replay and (root / receipt).is_file():
                replay[receipt] = (".".join(relative.parts), None, receipt)
                discovered.add(receipt)
        # Any other receipt with a source map whose package holds exactly one
        # module that names it and has a write or check option.
        for path in sorted((root / "model").rglob("*.json")):
            relative = path.relative_to(root)
            receipt = relative.as_posix()
            if (receipt in replay or {"historical", "predecessor"} & set(relative.parts)
                    or path.name in ("manifest.json", "MANIFEST.json")
                    or "campaign" in path.name or "contract" in path.name):
                continue
            try:
                value = strict_json.loads(path.read_text())
            except (OSError, ValueError, UnicodeDecodeError):
                continue
            if not any(isinstance(child, dict) and child for _, child, _ in source_maps(value)):
                continue
            owners_found = []
            for module_file in sorted(path.parent.glob("*.py")):
                if module_file.name.startswith("test_") or module_file.name == "__init__.py":
                    continue
                text = module_file.read_text(errors="ignore")
                if path.name in text and re.search(r"--(write|output|check)\b", text):
                    owners_found.append(module_file)
            if len(owners_found) == 1:
                module_path = owners_found[0].relative_to(root).with_suffix("")
                parts = module_path.parts[:-1] if module_path.name == "__main__" else module_path.parts
                replay[receipt] = (".".join(parts), None, receipt)
                discovered.add(receipt)
    values = {}
    for name in list(replay):
        try:
            values[name] = load_receipt(root, name, strict_json)
        except (OSError, ValueError, EOFError):
            if name not in discovered:
                raise
            replay.pop(name)
            discovered.discard(name)
    by_path, outputs, dropped = output_owners(root, values, discovered)
    for name in dropped:
        replay.pop(name)
        values.pop(name)
        discovered.discard(name)
    receipt_digests = {name: digest(root / name) for name in replay}
    by_digest = {}
    for name, value in receipt_digests.items():
        by_digest.setdefault(value, set()).add(name)
    graph = {}
    initial = {}
    digests = {}
    unreadable = {}
    for name in sorted(discovered):
        # A discovered receipt whose maps the authoritative checker would
        # reject keeps its owning verifier and is not planned here.
        try:
            list(bindings(root, name, values[name]))
        except ValueError as error:
            unreadable[name] = str(error)
            replay.pop(name)
            values.pop(name)
            outputs.pop(name)
            discovered.discard(name)
    by_path = {path: name for path, name in by_path.items() if name in replay}
    by_digest = {found: names & set(replay) for found, names in by_digest.items()}
    for name, value in values.items():
        graph[name] = set()
        try:
            found = list(bindings(root, name, value))
        except ValueError as error:
            # A pinned file that no longer exists makes the receipt stale.
            initial.setdefault(name, []).append({"error": str(error)})
            found = []
            for key, child, here in source_maps(value):
                if not isinstance(child, dict):
                    continue
                base = (root / name).parents[PACKAGE_LEVELS[key]] if key in PACKAGE_LEVELS else root
                for relative, expected in child.items():
                    if isinstance(relative, str) and not exempt(relative) and (base / relative).is_file():
                        found.append((here, relative, expected, (base / relative).resolve()))
        for here, relative, expected, target in found:
            if target in by_path and by_path[target] != name:
                graph[name].add(by_path[target])
            if target not in digests:
                digests[target] = digest(target)
            if expected != digests[target]:
                initial.setdefault(name, []).append({"map": here, "path": relative,
                    "target": target.relative_to(root).as_posix(), "expected": expected,
                    "actual": digests[target]})
        # A receipt that records another receipt's digest outside a source map
        # depends on it as well.
        for found in hex_values(value):
            graph[name].update(by_digest.get(found, set()) - {name})
    order = topological(graph)
    affected = set(initial)
    for name in order:
        if graph[name] & affected:
            affected.add(name)
    report = {"root": str(root), "authoritative_verifier_sha256": digest(root / "model/verify.py"),
              "execute": args.execute, "jobs": args.jobs, "registered_count": len(verify.REPLAYS),
              "discovered_count": len(discovered), "discovered": sorted(discovered),
              "discovered_with_unreadable_source_maps": unreadable,
              "initial_stale": initial,
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
    print(f"Registered={len(verify.REPLAYS)}, discovered={len(discovered)}, "
          f"directly stale={len(initial)}, affected={len(affected)}", flush=True)
    if not args.execute:
        print(f"Read-only plan: {report_path}", flush=True)
        return
    environment = dict(os.environ)
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[var] = "1"
    current = repository_digests(root)
    holders = {}
    for path, found in current.items():
        holders.setdefault(found, set()).add(path)
    current_lock = threading.Lock()
    def file_digest_for(own):
        def check(value):
            with current_lock:
                return bool(holders.get(value, set()) - own)
        return check

    def execute_one(name):
        module, command, _ = replay[name]
        interface, problem = writer_interface(root, replay[name])
        if problem:
            raise RuntimeError(f"{name}: {problem}")
        prefix = [sys.executable, "-W", "error", "-m", module]
        if command is not None:
            prefix.append(command)
        writer_module, writer_command = writer_row(replay[name])
        writer = [sys.executable, "-W", "error", "-m", writer_module]
        if writer_command is not None:
            writer.append(writer_command)
        before = load_receipt(root, name, strict_json)
        backed_up = {}
        for path in outputs[name]:
            relative = path.relative_to(root)
            backup = args.output / "before" / relative
            backup.parent.mkdir(parents=True, exist_ok=True)
            backup.write_bytes(path.read_bytes())
            backed_up[relative.as_posix()] = digest(path)
        result = {"receipt": name, "before_sha256": digest(root / name),
                  "backed_up_outputs": backed_up, "interface": interface}
        with report_lock:
            report["results"].append(result)
            save()
        print(f"Replaying {name}", flush=True)
        log_stem = module + ("_" + command.replace("-", "") if command else "")
        mode, flag = interface["write"]
        candidate_root = args.output.resolve() / "candidate" if mode == "path" else root
        candidate = candidate_root / name
        candidate.parent.mkdir(parents=True, exist_ok=True)

        def restore():
            for relative in backed_up:
                (root / relative).write_bytes((args.output / "before" / relative).read_bytes())

        def run(operation, argv):
            log = args.output / (log_stem + "." + operation + ".log")
            with report_lock:
                result[operation] = {"argv": argv, "log": str(log)}
                save()
            printed = args.output / (log_stem + "." + operation + ".stdout")
            if operation == "write" and mode == "bare":
                # A bare producer may print its receipt instead of writing it.
                with log.open("w") as stream, printed.open("w") as out:
                    completed = subprocess.run(argv, cwd=root, env=environment, stdout=out, stderr=stream)
                receipt_path = root / name
                if not completed.returncode and not receipt_path.exists():
                    text = printed.read_text()
                    try:
                        strict_json.loads(text)
                    except ValueError:
                        pass
                    else:
                        receipt_path.write_text(text)
                stream_text = printed.read_text()
                with log.open("a") as stream:
                    stream.write(stream_text)
            else:
                with log.open("w") as stream:
                    completed = subprocess.run(argv, cwd=root, env=environment, stdout=stream, stderr=subprocess.STDOUT)
            with report_lock:
                result[operation]["returncode"] = completed.returncode
                save()
            if completed.returncode:
                if mode != "path":
                    restore()
                raise RuntimeError(f"{name} {operation} failed; inspect {log}")

        # A number of evidence producers correctly refuse to overwrite an
        # existing path. A path-taking producer writes a fresh candidate, and
        # canonical outputs are replaced only after its cold checker passes.
        if mode == "bare":
            for relative in backed_up:
                (root / relative).unlink()
        started_ns = time.time_ns()
        declared = {path.resolve() for path in by_path}
        failure = None
        try:
            run("write", writer + (flag.split() + [str(candidate)] if mode == "path" else [flag] if flag else []))
        except RuntimeError as error:
            failure = error
        # Files the producer wrote besides its declared outputs.
        side_effects = touched_since(root, started_ns, skip=declared)
        with report_lock:
            result["side_effect_files"] = sorted(p.relative_to(root).as_posix() for p in side_effects)
            save()
        if failure is not None:
            restore_from_git(root, side_effects)
            raise failure
        try:
            after = load_receipt(candidate_root, name, strict_json)
        except (OSError, ValueError, EOFError):
            if mode != "path":
                restore()
            restore_from_git(root, side_effects)
            raise RuntimeError(f"{name} write did not produce a readable receipt")
        candidate_outputs = owned_outputs(candidate_root, name, after)
        # A digest of anything this producer just wrote is payload, not a pin.
        own = set(candidate_outputs) | {root / relative for relative in backed_up} | side_effects
        own = {path.resolve() for path in own}
        diff = differences(science_payload(root, name, before), science_payload(root, name, after),
                           file_digest=file_digest_for(own))
        with report_lock:
            result["after_sha256"] = digest(candidate)
            result["non_provenance_differences"] = diff
            if diff and args.float_tolerance is not None:
                # Classification for the report only; science compares exactly.
                result["differences_beyond_float_tolerance"] = differences(
                    science_payload(root, name, before), science_payload(root, name, after),
                    file_digest=file_digest_for(own), float_tolerance=args.float_tolerance)
            save()
        if diff:
            if mode != "path":
                restore()
            restore_from_git(root, side_effects)
            raise RuntimeError(f"{name} changed non-provenance payload; backups={args.output / 'before'}")
        relative_outputs = {p.relative_to(candidate_root).as_posix(): p for p in candidate_outputs}
        if relative_outputs.keys() != backed_up.keys():
            if mode != "path":
                restore()
            restore_from_git(root, side_effects)
            raise RuntimeError(f"{name} changed its output catalogue; review candidate outputs")

        def commit():
            if any(digest(root / relative) != expected for relative, expected in backed_up.items()):
                raise RuntimeError(f"{name} canonical output changed during replay")
            # Archives first and the receipt last. Dependent DAG nodes are
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

        # A checker given a path reads the candidate; a flag-only checker reads
        # the canonical receipt, so the candidate is committed first and the
        # backed-up outputs return if the check fails. Some checkers accept
        # only the repository-relative canonical path.
        canonical_check = interface["check"] in ("flag", "verify")
        if mode == "path" and canonical_check:
            commit()
        try:
            if interface["check"] == "path":
                run("check", prefix + ["--check", str(candidate) if mode == "path" else name])
            elif interface["check"] == "flag":
                run("check", prefix + ["--check"])
            elif interface["check"] == "verify":
                run("check", prefix + ["--verify"])
        except RuntimeError:
            if mode == "path" and canonical_check:
                restore()
            restore_from_git(root, side_effects)
            raise
        if mode == "path" and not canonical_check:
            commit()
        with current_lock:
            for relative in relative_outputs:
                path = (root / relative).resolve()
                holders.get(current.get(path), set()).discard(path)
                current[path] = digest(path)
                holders.setdefault(current[path], set()).add(path)
        return name

    pending = set(affected)
    complete = set(replay) - affected
    blocked = set()
    running = {}
    failures = []
    stopped = False
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        while pending or running:
            if args.stop_after is not None and time.monotonic() - started > 60 * args.stop_after:
                stopped = True
            if not stopped and (args.keep_going or not failures):
                ready = [name for name in order if name in pending and graph[name] <= complete]
                for name in ready[:args.jobs - len(running)]:
                    pending.remove(name)
                    running[pool.submit(execute_one, name)] = name
            if not running:
                if pending and not failures and not stopped:
                    failures.append({"receipt": None, "error": "No ready DAG nodes remain"})
                break
            finished, _ = wait(running, timeout=60, return_when=FIRST_COMPLETED)
            for future in finished:
                name = running.pop(future)
                try:
                    future.result()
                except Exception as error:
                    failures.append({"receipt": name, "error": str(error)})
                    if args.keep_going:
                        # Dependents of a failed receipt never start.
                        for other in order:
                            if other in pending and graph[other] & (blocked | {name}):
                                pending.remove(other)
                                blocked.add(other)
                else:
                    complete.add(name)
            with report_lock:
                report["failures"] = list(failures)
                report["completed"] = sorted(complete & affected)
                report["blocked_by_failure"] = sorted(blocked)
                report["not_scheduled"] = sorted(pending)
                report["running"] = sorted(running.values())
                report["stopped_by_time_limit"] = stopped
                save()
        # Exiting the executor waits for every already-running producer.
        # Without --keep-going no further node starts after the first failure.
    if failures:
        raise SystemExit("Stopped after awaiting running nodes: " + repr(failures))
    if stopped:
        raise SystemExit(f"Time limit reached; {len(pending)} affected receipts not started: {report_path}")
    # Fresh authoritative validation, not the local copy of discovery logic.
    try:
        from model import verify_current as authority
    except ImportError:
        authority = verify
    authority.verify_receipt_source_bindings()
    report["authoritative_final_source_check"] = "PASS"
    save()
    print(f"All {len(affected)} affected receipts replayed and checked: {report_path}", flush=True)


if __name__ == "__main__":
    main()
