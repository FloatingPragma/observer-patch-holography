"""Conservative observer regression selection; never a scientific replay cache.

Selection follows both old and current Python import edges and every declared
receipt input edge, including nested child receipts. A test without a receipt
owner, with dynamic imports, or with unresolved repository imports always runs.
Changes to configuration, module membership or a receipt's input-key inventory
force the original full suite. Unknown executable/data paths also force full.
The unchanged global source guards run separately before this plan is accepted.
Scheduled, manual and release runs keep the complete original suite.

This uses the repository's declared source contracts, not a claim that arbitrary
Python file access can be inferred statically. A receipt may own a test only
when its source inventory explicitly includes that test and it is a nearest
enclosing package receipt. Directory proximity alone is never ownership. Unowned tests are never silently omitted.
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict
import configparser
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import posixpath
from statistics import median
import math
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HASH = re.compile(r"[0-9a-f]{64}\Z")
SOURCE_MAPS = {"sources", "source_sha256", "input_sha256", "source_bindings", "source_hashes",
               "source_binding_sha256", "source_control_bindings", "source_pins"}
PACKAGE_LEVELS = {"implementation_sha256": 0, "source_interface_hashes": 1}
ENDPOINT_PINS = "root.retrospective_endpoint_audit.source.pins"
GLOBAL = {"pytest.ini", "pyproject.toml", "setup.cfg", "conftest.py", "requirements.txt",
          "model/requirements.txt", "model/__init__.py", "tools/__init__.py",
          "model/pytest_shards.py", "model/pytest_shard_costs.json"}
DOCUMENT_SUFFIXES = {".md", ".tex", ".bib", ".svg", ".png", ".pdf"}


ARCHIVE = "evidence/observer_dynamics_20260925/oph-physics-sim"
ARCHIVE_NAMES_SHA256 = "64b88e8de19eaa19a1923f85de303b15c6023bfbfc9a66f3fa4504254a3871cb"
ARCHIVE_CONTENT_SHA256 = "c6386d8adb6114b5221b677e8f6616bb92cf050e9fbe12da392d2ccbebb665fe"
LOADER_SHA256 = {
    "model/observer/p_response/experiment.py": "44049e5702af83ab8a62e40f8f7d3b25b088ed24e47078a0b412c3b2400ed61e",
    "model/common_source/slow_geometry.py": "91e808e06cc4fd883429405f3454dc53ee48eab7152eac36db93591749eebf38",
    "model/federated_observer/source.py": "4a91b7f254fcee89c93b1ce3e9972ec0c809290ec987dbc87fcb0392319e549d",
}


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def fixed_archive(root, expected_names=ARCHIVE_NAMES_SHA256, expected_content=ARCHIVE_CONTENT_SHA256):
    paths = sorted(p.relative_to(root).as_posix() for p in (root / ARCHIVE).rglob("*")
                   if p.is_file() and "__pycache__" not in p.parts and ".pytest_cache" not in p.parts
                   and p.suffix not in {".pyc", ".pyo"})
    if digest_json(paths) != expected_names:
        raise ValueError("reviewed archive file membership changed")
    bindings = {}
    for path in paths:
        target = root / path
        if not target.resolve().is_relative_to(root):
            raise ValueError("reviewed archive escapes checkout")
        raw = target.read_bytes()
        bindings[path] = hashlib.sha256(raw).hexdigest()
        if path.endswith(".py"):
            for node in ast.walk(ast.parse(raw, filename=path)):
                imports_here = ([a.name for a in node.names] if isinstance(node, ast.Import)
                                else [node.module or ""] if isinstance(node, ast.ImportFrom) and not node.level else [])
                if any(name.split(".")[0] in {"model", "tools"} for name in imports_here):
                    raise ValueError("archived library imports current RER code")
                if (isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant)
                        and isinstance(node.args[0].value, str)
                        and node.args[0].value.split(".")[0] in {"model", "tools"}
                        and ((isinstance(node.func, ast.Name) and node.func.id == "__import__")
                             or (isinstance(node.func, ast.Attribute) and node.func.attr == "import_module"))):
                    raise ValueError("archived library dynamically imports current RER code")
    if digest_json(bindings) != expected_content:
        raise ValueError("reviewed archive source bytes changed")
    return frozenset(paths)


def reviewed_loaders(root):
    """Three finite loader contracts; changed bytes require complete fallback."""
    for path, expected in LOADER_SHA256.items():
        if hashlib.sha256((root / path).read_bytes()).hexdigest() != expected:
            raise ValueError("fixed loader requires renewed review: " + path)
    archive = fixed_archive(root)
    def constants(path):
        result = {}
        def literal(node):
            if isinstance(node, ast.Name):
                return result[node.id]
            if isinstance(node, ast.Dict):
                return {literal(k): literal(v) for k, v in zip(node.keys, node.values)}
            return ast.literal_eval(node)
        for node in ast.parse((root / path).read_bytes()).body:
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                try:
                    result[node.targets[0].id] = literal(node.value)
                except (ValueError, TypeError, KeyError):
                    pass
        return result
    response, geometry, federation = LOADER_SHA256
    pins = constants(response)["PINS"]
    for name, expected in pins.items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError("fixed response input pin changed: " + name)
    for name, expected in constants(federation)["SOURCE_PINS"].items():
        if hashlib.sha256((root / ARCHIVE / name).read_bytes()).hexdigest() != expected:
            raise ValueError("fixed federation pin changed: " + name)
    spec = "evidence/native_geometric_source_20260925/spec.json"
    if hashlib.sha256((root / spec).read_bytes()).hexdigest() != "a6b14ad3c89ab17d3c213375022975e1a35239a82810406b73157a6b856cadcc":
        raise ValueError("fixed geometry source contract changed")
    for name, expected in json.loads((root / spec).read_bytes())["source_sha256"].items():
        if hashlib.sha256((root / ARCHIVE / name).read_bytes()).hexdigest() != expected:
            raise ValueError("fixed geometry pin changed: " + name)
    return {response: frozenset(pins), geometry: archive | {spec}, federation: archive}


def git(root, *args):
    return subprocess.check_output(["git", *args], cwd=root)


def source_paths(value, receipt="model/results.json", location="root", source_bases=None):
    """All preflight-recognized nested inputs, including package-relative maps."""
    found = set()
    source_bases = source_bases or {}
    def walk(value, location):
        if isinstance(value, dict):
            for key, child in value.items():
                here = location + "." + key
                recognized = key in SOURCE_MAPS or key in PACKAGE_LEVELS or here == ENDPOINT_PINS
                if recognized:
                    if isinstance(child, dict):
                        entries = child.items()
                    elif key == "source_bindings" and isinstance(child, list):
                        if not all(isinstance(row, dict) and set(row) == {"path", "sha256"} for row in child):
                            raise ValueError("malformed source-binding records")
                        entries = [(row["path"], row["sha256"]) for row in child]
                    elif key == "sources":
                        entries = []  # Mixed external metadata belongs to its owner.
                    else:
                        raise ValueError("malformed local source map")
                    seen = set()
                    for name, digest in entries:
                        if not isinstance(name, str) or not isinstance(digest, str) or not HASH.fullmatch(digest):
                            if key == "sources":
                                continue
                            raise ValueError("invalid local source binding")
                        if name in seen:
                            raise ValueError("duplicate source binding")
                        seen.add(name)
                        override = source_bases.get((receipt, here))
                        if override is None and key not in PACKAGE_LEVELS and name in found:
                            continue
                        base = (PurePosixPath(override) if override is not None else
                                PurePosixPath(receipt).parents[PACKAGE_LEVELS[key]] if key in PACKAGE_LEVELS else PurePosixPath("."))
                        path = PurePosixPath(posixpath.normpath(str(base / name)))
                        if path.is_absolute() or ".." in path.parts or "://" in name:
                            if key == "sources":
                                continue
                            raise ValueError("escaping declared input: " + name)
                        found.add(path.as_posix())
                if isinstance(child, (dict, list)):
                    walk(child, here)
        elif isinstance(value, list):
            for i, child in enumerate(value):
                if isinstance(child, (dict, list)):
                    walk(child, f"{location}[{i}]")
    walk(value, location)
    return found

def imports(path, text, names, loader_contracts=None):
    """Resolve relative imports and package initializers; ambiguity stays visible."""
    tree = ast.parse(text, filename=path)
    module = path[:-3].replace("/", ".")
    package = module.removesuffix(".__init__") if module.endswith(".__init__") else module.rpartition(".")[0]
    edges, ambiguous = set(), False
    loader_aliases = {alias.asname or alias.name.split(".")[0]
                      for node in ast.walk(tree) if isinstance(node, ast.Import)
                      for alias in node.names if alias.name.split(".")[0] in {"importlib", "builtins", "sys"}}

    parents = {id(child): node for node in ast.walk(tree) for child in ast.iter_child_nodes(node)}
    reflection_maps, reflected_values = set(), set()
    def reflection(node):
        return ((isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                 and node.func.id in {"vars", "globals", "locals"})
                or (isinstance(node, ast.Attribute) and node.attr == "__dict__")
                or (isinstance(node, ast.Name) and node.id in reflection_maps))
    assignments = [node for node in ast.walk(tree) if isinstance(node, (ast.Assign, ast.AnnAssign))]
    changed_aliases = True
    while changed_aliases:
        before = (len(reflection_maps), len(reflected_values))
        for node in assignments:
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            names_assigned = {t.id for t in targets if isinstance(t, ast.Name)}
            value = node.value
            if reflection(value):
                reflection_maps.update(names_assigned)
            if ((isinstance(value, ast.Subscript) and reflection(value.value))
                    or (isinstance(value, ast.Name) and value.id in reflected_values)):
                reflected_values.update(names_assigned)
        changed_aliases = before != (len(reflection_maps), len(reflected_values))

    def add(name, required=False):
        nonlocal ambiguous
        stem = name.replace(".", "/")
        candidates = {stem + ".py", stem + "/__init__.py"}
        hits = candidates & names
        edges.update(hits)
        parts = stem.split("/")
        edges.update("/".join(parts[:i]) + "/__init__.py" for i in range(1, len(parts))
                     if "/".join(parts[:i]) + "/__init__.py" in names)
        if required and name.startswith(("model.", "tools.")) and not hits:
            ambiguous = True

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] in {"runpy", "subprocess"} or (alias.name.startswith("importlib.") and alias.name != "importlib.util"):
                    ambiguous = True
                add(alias.name, required=True)
        elif isinstance(node, ast.ImportFrom):
            target = "." * node.level + (node.module or "")
            if node.level:
                target = importlib.util.resolve_name(target, package)
            if (target.split(".")[0] in {"runpy", "subprocess"}
                    or (target.startswith("importlib") and not
                        ((target == "importlib.util" and all(a.name == "resolve_name" for a in node.names))
                         or (target == "importlib" and all(a.name == "util" for a in node.names))))):
                ambiguous = True
            add(target, required=True)
            for alias in node.names:
                if alias.name == "*":
                    ambiguous = True
                else:
                    add(target + "." + alias.name)
        elif isinstance(node, ast.Attribute) and node.attr in {
                "import_module", "run_module", "run_path", "spec_from_file_location",
                "module_from_spec", "SourceFileLoader", "__globals__", "__builtins__"}:
            ambiguous = True
        elif (isinstance(node, ast.Attribute) and node.attr == "__dict__"
              and isinstance(node.value, ast.Name) and node.value.id in loader_aliases):
            ambiguous = True
        elif isinstance(node, ast.Subscript) and reflection(node.value):
            parent = parents.get(id(node))
            if ((isinstance(parent, ast.Call) and parent.func is node)
                    or (isinstance(node.slice, ast.Constant) and node.slice.value in
                        {"__import__", "eval", "exec", "compile", "import_module", "run_module"})):
                ambiguous = True
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id in {"__import__", "eval", "exec"} | reflected_values:
                ambiguous = True
            elif (isinstance(node.func, ast.Name) and node.func.id == "getattr" and node.args
                  and isinstance(node.args[0], ast.Name) and node.args[0].id in loader_aliases):
                ambiguous = True
            elif isinstance(node.func, ast.Attribute) and node.func.attr in {
                    "import_module", "run_module", "run_path", "spec_from_file_location"}:
                ambiguous = True
    if loader_contracts and path in loader_contracts:
        if hashlib.sha256(text.encode()).hexdigest() != LOADER_SHA256[path]:
            raise ValueError("fixed loader source changed during planning")
        edges.update(loader_contracts[path])
        ambiguous = False
    return edges, ambiguous


def closure(graph, seeds):
    pending, seen = list(seeds), set()
    while pending:
        path = pending.pop()
        if path not in seen:
            seen.add(path)
            pending.extend(graph.get(path, ()))
    return seen


def partition_modules(modules, count, costs):
    ordered = sorted(set(modules))
    if not ordered:
        return [[]]
    if not 1 <= count <= len(ordered):
        raise ValueError("partition count exceeds selected modules")
    default = median(costs.values()) if costs else 1.0
    loads, groups = [0.0] * count, [[] for _ in range(count)]
    for path in sorted(ordered, key=lambda p: -costs.get(p, default)):
        index = min(range(count), key=lambda i: loads[i])
        groups[index].append(path)
        loads[index] += costs.get(path, default)
    return [sorted(group) for group in groups]


def plan_from_files(files, changed, old=None, force_reason=None, source_bases=None, metadata_cache=None, loader_contracts=None, all_paths=None):
    """Pure plan builder. Old bytes are required for every changed existing file."""
    old = old or {}
    metadata_cache = {} if metadata_cache is None else metadata_cache
    names = set(files) | set(old) | set(all_paths or ())
    graph = defaultdict(set)
    executable_reverse = defaultdict(set)
    receipts, dynamic, ambiguous_receipts = {}, set(), set()
    ignored = {"model/completion/historical"}
    if "pytest.ini" in files:
        config = configparser.ConfigParser(interpolation=None)
        config.read_string(files["pytest.ini"].decode())
        tokens = shlex.split(config.get("pytest", "addopts", fallback=""), comments=True)
        for i, token in enumerate(tokens):
            if token.startswith("--ignore="):
                ignored.add(token.split("=", 1)[1].rstrip("/"))
            elif token == "--ignore" and i + 1 < len(tokens):
                ignored.add(tokens[i + 1].rstrip("/"))
    tests = sorted(p for p in files if p.startswith("model/")
                   and PurePosixPath(p).name.startswith("test_") and p.endswith(".py")
                   and not any(p == q or p.startswith(q + "/") for q in ignored))
    reasons = [] if force_reason is None else [force_reason]
    # Parse executable imports lazily from the test roots. A Python file held
    # only as documentary input is a dependency, not necessarily executable.
    pending, parsed = list(tests), set()
    while pending:
        path = pending.pop()
        if path in parsed:
            continue
        parsed.add(path)
        versions = [files[path]] if path in files else []
        if path in old and old[path] not in versions:
            versions.append(old[path])
        for raw in versions:
            try:
                edges, unknown = imports(path, raw.decode("utf-8"), names, loader_contracts)
                graph[path].update(edges)
                for dependency in edges:
                    executable_reverse[dependency].add(path)
                pending.extend(p for p in edges if p.endswith(".py") and p not in parsed
                               and not (loader_contracts and path in loader_contracts and p.startswith(ARCHIVE + "/")))
                if unknown:
                    dynamic.add(path)
            except (ValueError, SyntaxError, UnicodeError, ImportError) as error:
                reasons.append("unresolved executable dependency: " + path + ": " + str(error))
    for path in sorted(p for p in names if p.endswith(".json")):
        versions = [files[path]] if path in files else []
        if path in old and old[path] not in versions:
            versions.append(old[path])
        inventories = []
        for raw in versions:
            try:
                cache_key = (path, hashlib.sha256(raw).hexdigest())
                if cache_key not in metadata_cache:
                    metadata_cache[cache_key] = source_paths(json.loads(raw), path, source_bases=source_bases)
                inputs = metadata_cache[cache_key]
                inventories.append(inputs)
                if inputs:
                    graph[path].update(inputs)
                    receipts.setdefault(path, set()).update(inputs)
            except (ValueError, UnicodeError) as error:
                # An unchanged invalid-JSON mutation fixture is not a receipt.
                # Its owning test is still reached through imports/input maps.
                ambiguous_receipts.add(path)
                if path in changed:
                    reasons.append("unresolved dependency inventory: " + path + ": " + str(error))
        if len(inventories) == 2 and inventories[0] != inventories[1]:
            reasons.append("receipt input inventory changed: " + path)
    test_set = set(tests)
    candidates, owners = defaultdict(set), defaultdict(set)
    for receipt, inputs in receipts.items():
        directory = PurePosixPath(receipt).parent
        for test in inputs & test_set:
            if PurePosixPath(test).is_relative_to(directory):
                candidates[test].add(receipt)
    for test, options in candidates.items():
        depth = max(len(PurePosixPath(p).parent.parts) for p in options)
        owners[test].update(p for p in options if len(PurePosixPath(p).parent.parts) == depth)
        graph[test].update(owners[test])
    known = set(graph) | {p for deps in graph.values() for p in deps}
    for path in changed:
        if path in GLOBAL or path.startswith(".github/") or PurePosixPath(path).name == "conftest.py" or path.startswith("tools/ci_affected_model"):
            reasons.append("shared execution configuration changed: " + path)
        if path.endswith(".py") and ((path in files) != (path in old)):
            reasons.append("Python module membership changed: " + path)
        if path not in known and PurePosixPath(path).suffix not in DOCUMENT_SUFFIXES:
            reasons.append("changed path has no declared dependency owner: " + path)
    reverse = defaultdict(set)
    for consumer, dependencies in graph.items():
        for dependency in dependencies:
            reverse[dependency].add(consumer)
    impacted = closure(reverse, changed)
    dynamic_consumers = closure(executable_reverse, dynamic)
    metadata_consumers = closure(reverse, ambiguous_receipts)
    ambiguous_consumers = dynamic_consumers | metadata_consumers
    always = sorted(test for test in tests if not owners[test] or test in ambiguous_consumers)
    affected = sorted(set(tests) & impacted)
    selected = sorted(set(always) | set(affected))
    if not tests:
        reasons.append("no test inventory found")
    mode = "full" if reasons else "affected"
    if mode == "full":
        selected = tests
    # Keep one worker for global custody controls even when no test is affected.
    costs = {}
    if "model/pytest_shard_costs.json" in files:
        try:
            raw_costs = json.loads(files["model/pytest_shard_costs.json"])["seconds"]
            costs = {p: float(c) for p, c in raw_costs.items() if float(c) > 0 and math.isfinite(float(c))}
        except (ValueError, TypeError, KeyError):
            costs = {}  # Costs balance work; they never exclude tests.
    default = median(costs.values()) if costs else 1.0
    total_cost = sum(costs.get(p, default) for p in tests)
    selected_cost = sum(costs.get(p, default) for p in selected)
    count = 16 if mode == "full" else max(1, min(len(selected), 16, math.ceil(selected_cost / 600)))
    partitions = partition_modules(selected, min(count, len(selected)), costs) if selected else [[]]
    estimated_loads = [sum(costs.get(path, default) for path in group) for group in partitions]
    return {"schema": "oph.affected_model_plan.v1", "mode": mode,
            "changed": sorted(set(changed)), "full_reasons": sorted(set(reasons)),
            "all_test_modules": tests, "selected_test_modules": selected,
            "always_test_modules": always, "affected_dependencies": affected,
            "always_reason_counts": {"unowned": sum(not owners[t] for t in tests),
                                     "dynamic_imports": len(set(tests) & dynamic_consumers),
                                     "ambiguous_metadata": len(set(tests) & metadata_consumers)},
            "ambiguous_metadata_files": sorted(ambiguous_receipts),
            "dynamic_import_modules": sorted(dynamic),
            "test_module_count": len(tests), "selected_module_count": len(selected),
            "shard_count": count, "matrix": list(range(count)),
            "retained_cost_fraction": selected_cost / max(total_cost, 1),
            "collection_mode": "original-full" if mode == "full" else "assigned-modules",
            "module_partitions": [] if mode == "full" else partitions,
            "retained_estimated_selected_seconds": selected_cost,
            "retained_estimated_max_shard_seconds": max(estimated_loads, default=0),
            "largest_selected_module_costs": sorted(
                ((path, costs.get(path, default)) for path in selected),
                key=lambda row: (-row[1], row[0]))[:10],
            "full_source_guards_required": True,
            "scientific_comparisons_changed": False}


def build_plan(root, base=None, force_full=False, changed_override=None):
    root = root.resolve()
    names = git(root, "ls-files", "--cached", "--others", "--exclude-standard", "-z").decode().split("\0")
    if force_full:
        # Full coverage needs no affected-dependency discovery at all.
        fixtures = {p: b"" for p in names if p.startswith("model/")
                    and PurePosixPath(p).name.startswith("test_") and p.endswith(".py")}
        for path in ("pytest.ini", "model/pytest_shard_costs.json"):
            if (root / path).is_file():
                fixtures[path] = (root / path).read_bytes()
        result = plan_from_files(fixtures, [], force_reason="scheduled/manual/release complete coverage")
        result.update(base=base, head=git(root, "rev-parse", "HEAD").decode().strip())
        return result
    paths = [p for p in names if p and (root / p).is_file()
             and (p == "pytest.ini" or p.endswith(".py") or (p.startswith("model/") and p.endswith(".json")))]
    files = {p: (root / p).read_bytes() for p in paths}
    reason = "scheduled/manual/release complete coverage" if force_full else None
    source_bases, metadata_cache, loader_contracts = {}, {}, {}
    if all((root / path).is_file() for path in LOADER_SHA256):
        try:
            loader_contracts = reviewed_loaders(root)
        except (ValueError, OSError, SyntaxError, KeyError) as error:
            reason = "unresolved fixed-loader contract: " + str(error)
    # Reuse the exact, authenticated exception used by the independent guard.
    # Fixture repositories have no archived specification.
    if (root / "evidence/native_geometric_source_20260925/spec.json").is_file():
        try:
            from tools.audit_receipt_sources import archived_source_map_roots
            source_bases = {key: base.relative_to(root).as_posix()
                            for key, base in archived_source_map_roots(root).items()}
        except (ValueError, OSError, subprocess.CalledProcessError) as error:
            reason = "unresolved archive input convention: " + str(error)
    queue, seen_json = [p for p in files if p.endswith(".json")], set()
    while queue:
        path = queue.pop()
        if path in seen_json:
            continue
        seen_json.add(path)
        try:
            inputs = source_paths(json.loads(files[path]), path, source_bases=source_bases)
            metadata_cache[path, hashlib.sha256(files[path]).hexdigest()] = inputs
        except (ValueError, UnicodeError):
            continue  # The planner below conservatively marks affected consumers.
        for dependency in inputs:
            target = root / dependency
            if dependency.endswith(".json") and dependency not in files and target.is_file():
                if not target.resolve().is_relative_to(root):
                    reason = "escaping declared JSON input: " + dependency
                    continue
                files[dependency] = target.read_bytes()
                queue.append(dependency)
    old, changed = {}, []
    if changed_override is not None:
        changed = list(changed_override)
        old = {p: files[p] for p in changed if p in files}
    elif base:
        try:
            if not re.fullmatch(r"[0-9a-fA-F]{40,64}", base) or set(base) == {"0"}:
                raise ValueError("invalid base object")
            git(root, "cat-file", "-e", base + "^{commit}")
            changed = git(root, "diff", "--name-only", "--no-renames", "-z", base, "HEAD", "--").decode().split("\0")
            changed = [p for p in changed if p]
            for p in changed:
                if p.endswith(".py") or p.endswith(".json"):
                    result = subprocess.run(["git", "show", base + ":" + p], cwd=root, capture_output=True)
                    if result.returncode == 0:
                        old[p] = result.stdout
        except (ValueError, subprocess.CalledProcessError):
            reason = "Git base unavailable; complete coverage required"
    else:
        reason = reason or "Git base not supplied; complete coverage required"
    result = plan_from_files(files, changed, old, reason, source_bases, metadata_cache, loader_contracts, names)
    result["reviewed_loader_boundaries"] = {
        "loader_sha256": LOADER_SHA256 if loader_contracts else {},
        "archive_names_sha256": ARCHIVE_NAMES_SHA256 if loader_contracts else None,
        "archive_content_sha256": ARCHIVE_CONTENT_SHA256 if loader_contracts else None,
        "scope": "Only three fixed-source loaders; complete frozen archive is an opaque input, not generic dynamic-code exemption."}
    result["base"] = base
    result["head"] = git(root, "rev-parse", "HEAD").decode().strip()
    return result


FAST_TEST = "model/particle_identification/test_completion_certificate.py"
FAST_SOURCE_HASHES = {
    "model/native_effect_synthesis/synthesis.py": "559ea14f06b7d30ab017b2113feeeebc4d90a517bf85257fd3251d098cb4e372",
    "model/native_program/receipt.py": "9f507a9f2e14ee628fa3fa420e5b3ba393ccd8d9b5ca475314da0d12e3bc3cfd",
    "model/native_quantum_gauss/packet.py": "3a4b7e8b7e172d366fbe80e4db2bb6862fe00f478edcdef69a9f1ac9dd88d248",
    "model/native_sector_readout/packet.py": "914a3077bc6525e18740c247d9e5d50f63a85af8499b0133edfb68f44c22e4d1",
    "model/native_fair_local_preparation/packet.py": "880103c4330faf0813b82aa44c2391fab9a418ef2fb5183504effa2f137c71ac",
    "tools/immutable_walsh_replay.py": "69c0d8f565df654129e2d0aec648d62b7ecef5cc8cca3fea5cff4d878e55fa12",
    "tools/exact_walsh.py": "913ef7ab6d7bd9e567cfab5cc78dffbff36a5ea98f9e49f25e152252982b41f1",
    "tools/immutable_custody_replay.py": "19f030d5c246efeb8a4556154a3ac0789ea2974d1483813f52f8bfb6d9c2ad61",
    "tools/immutable_custody_pytest.py": "6867da943ec188b82eca69c60aba018f42ac2ed480068395ae7c8376a7f85ac5",
    "pytest.ini": "4f40f83f8d1457725aba1ca61b9d9f251ec33e561316395dfb793c7ef72a7448",
    FAST_TEST: "522b78c013be67a508aba3b48e5ee32fea8b0fc95b1120105fefd6147a1fad05",
}


def apply_fast_route(plan, root, requested):
    eligible = requested and FAST_TEST in plan["selected_test_modules"]
    if eligible:
        eligible = all((root / p).is_file() and hashlib.sha256((root / p).read_bytes()).hexdigest() == h
                       for p, h in FAST_SOURCE_HASHES.items())
        eligible = eligible and not any((root / p).exists() or (root / p).is_symlink() for p in
                     ("conftest.py", "model/conftest.py", "model/particle_identification/conftest.py"))
    plan["separate_test_modules"] = [FAST_TEST] if eligible else []
    plan["fast_identification"] = bool(eligible)
    plan["fast_route_source_sha256"] = FAST_SOURCE_HASHES if eligible else {}
    ordinary = sorted(set(plan["selected_test_modules"]) - set(plan["separate_test_modules"]))
    plan["ordinary_test_modules"] = ordinary
    if eligible and plan["collection_mode"] == "assigned-modules":
        raw = json.loads((root / "model/pytest_shard_costs.json").read_bytes())["seconds"]
        costs = {p: float(v) for p, v in raw.items() if float(v) > 0 and math.isfinite(float(v))}
        default = median(costs.values()) if costs else 1.0
        total = sum(costs.get(p, default) for p in ordinary)
        count = max(1, min(16, len(ordinary), math.ceil(total / 600)))
        plan["shard_count"], plan["matrix"] = count, list(range(count))
        plan["module_partitions"] = partition_modules(ordinary, count, costs)
        plan["retained_estimated_ordinary_seconds"] = total
        plan["retained_estimated_max_shard_seconds"] = max(
            (sum(costs.get(p, default) for p in g) for g in plan["module_partitions"]), default=0)
    return plan


def authenticate_plan(root, path, expected, base):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate plan key")
            result[key] = value
        return result
    def reject_constant(value):
        raise ValueError("non-finite plan number: " + value)
    plan = json.loads(path.read_text(), object_pairs_hook=unique_pairs, parse_constant=reject_constant)
    claimed = plan.pop("plan_sha256", None)
    actual = hashlib.sha256(json.dumps(plan, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if not expected or expected != actual or claimed != actual:
        raise ValueError("selection plan SHA256 differs from prerequisite job")
    if plan.get("head") != git(root, "rev-parse", "HEAD").decode().strip() or plan.get("base") != base:
        raise ValueError("selection plan Git head/base mismatch")
    if plan.get("selector_sha256") != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
        raise ValueError("selection plan selector source mismatch")
    if plan.get("mode") not in {"full", "affected"} or type(plan.get("shard_count")) is not int:
        raise ValueError("invalid selection mode or shard count")
    count = plan["shard_count"]
    if not 1 <= count <= 16 or plan.get("matrix") != list(range(count)):
        raise ValueError("invalid selection matrix")
    if plan["mode"] == "full" and count != 16:
        raise ValueError("full coverage must retain all sixteen shards")
    for key in ("all_test_modules", "selected_test_modules"):
        paths = plan.get(key)
        if not isinstance(paths, list) or not all(isinstance(p, str) and p.startswith("model/")
                and PurePosixPath(p).name.startswith("test_") and p.endswith(".py")
                and ".." not in PurePosixPath(p).parts for p in paths) or paths != sorted(set(paths)):
            raise ValueError("invalid test module inventory")
    if not set(plan["selected_test_modules"]) <= set(plan["all_test_modules"]):
        raise ValueError("selected tests escape complete inventory")
    separate = plan.get("separate_test_modules", [])
    ordinary = plan.get("ordinary_test_modules", plan["selected_test_modules"])
    if (separate not in ([], [FAST_TEST]) or not isinstance(ordinary, list)
            or any(not isinstance(p, str) for p in ordinary)
            or ordinary != sorted(set(ordinary)) or set(ordinary) & set(separate)
            or sorted(ordinary + separate) != plan["selected_test_modules"]):
        raise ValueError("ordinary and separate routes drop or duplicate selected tests")
    if separate and (not plan.get("fast_identification") or plan.get("fast_route_source_sha256") != FAST_SOURCE_HASHES
                     or any(not (root / p).is_file() or hashlib.sha256((root / p).read_bytes()).hexdigest() != h
                            for p, h in FAST_SOURCE_HASHES.items())):
        raise ValueError("exact-test route requires renewed source review")
    if plan.get("collection_mode") == "assigned-modules":
        partitions = plan.get("module_partitions")
        if not isinstance(partitions, list) or len(partitions) != count or not all(isinstance(g, list) for g in partitions):
            raise ValueError("invalid module partitions")
        flat = [p for group in partitions for p in group]
        if any(not isinstance(p, str) for p in flat) or len(flat) != len(set(flat)) or sorted(flat) != ordinary:
            raise ValueError("module partition duplicates or drops a selected test")
    elif plan.get("collection_mode") != "original-full" or plan["mode"] != "full":
        raise ValueError("invalid collection mode")
    plan["plan_sha256"] = claimed
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--base")
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--changed", action="append")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--fast-identification", action="store_true")
    parser.add_argument("--run-fast", action="store_true")
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--github-output", type=Path)
    parser.add_argument("--run-shard", type=int)
    parser.add_argument("--expect-plan-sha256")
    args = parser.parse_args()
    if args.plan:
        plan = authenticate_plan(args.root.resolve(), args.plan, args.expect_plan_sha256, args.base)
    else:
        plan = build_plan(args.root.resolve(), args.base, args.full, args.changed)
        plan = apply_fast_route(plan, args.root.resolve(), args.fast_identification)
        plan["selector_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        canonical = json.dumps(plan, sort_keys=True, separators=(",", ":"))
        plan["plan_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
        if args.expect_plan_sha256 and plan["plan_sha256"] != args.expect_plan_sha256:
            raise ValueError("selection plan differs from the prerequisite job")
    if args.output:
        args.output.write_text(json.dumps(plan, indent=2) + "\n")
    if args.github_output:
        with args.github_output.open("a") as out:
            out.write("matrix=" + json.dumps(plan["matrix"]) + "\n")
            out.write("mode=" + plan["mode"] + "\n")
            out.write("fast_identification=" + json.dumps(plan.get("fast_identification", False)) + "\n")
            out.write("shard_count=" + str(plan["shard_count"]) + "\n")
            out.write("plan_sha256=" + plan["plan_sha256"] + "\n")
    print(json.dumps({k: plan[k] for k in ("mode", "selected_module_count", "test_module_count", "shard_count", "full_reasons", "plan_sha256")}))
    if args.run_fast:
        if not args.plan or not args.audit or not plan.get("fast_identification"):
            parser.error("exact-test execution requires an authenticated selected plan and external audit")
        subprocess.run([sys.executable, "tools/immutable_walsh_replay.py", "--audit", str(args.audit)],
                       cwd=args.root, check=True)
    if args.run_shard is not None:
        if not 0 <= args.run_shard < plan["shard_count"]:
            parser.error("invalid shard")
        prepartitioned = plan["collection_mode"] == "assigned-modules"
        paths = plan["module_partitions"][args.run_shard] if prepartitioned else ["model"]
        if paths:
            # The unchanged plugin still installs authenticated documentary
            # fixtures. A one-part local shard keeps every assigned test node.
            index, count = (0, 1) if prepartitioned else (args.run_shard, plan["shard_count"])
            command = [sys.executable, "-m", "pytest", "-q", "-W", "error", "--durations=30",
                       "-p", "model.pytest_shards", "--model-shard-index", str(index),
                       "--model-shard-count", str(count),
                       "--deselect", "model/charged_lepton_response/test_packet.py::test_generated_results_and_runtime_custody", *paths]
            if plan.get("separate_test_modules") and not prepartitioned:
                command.insert(-1, "--ignore=" + FAST_TEST)
            subprocess.run(command, cwd=args.root, check=True)



if __name__ == "__main__":
    sys.path.insert(0, str(ROOT))
    main()
