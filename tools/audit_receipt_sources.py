"""Preflight retained receipt dependencies before their original expensive replay.

This read-only audit checks custody, not scientific correctness.  It discovers
CLI owners from literal OUTPUT declarations, augments the actual current replay
inventory, and visits every source map recursively.  A write preflight permits
only the explicitly selected roots to be stale; retained children stay strict.
Frozen captures are reported separately, never silently renewed.  Original
writers, cold comparisons and their documentary adapters remain mandatory.
"""
from __future__ import annotations

import argparse
import ast
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MAP_FIELDS = frozenset({"sources", "source_sha256", "source_bindings", "input_sha256",
    "source_binding_sha256", "source_control_bindings", "source_pins"})
PACKAGE_LEVELS = {"implementation_sha256": 0, "source_interface_hashes": 1}
ENDPOINT_PINS = "root.retrospective_endpoint_audit.source.pins"
HEX = re.compile(r"[0-9a-f]{64}\Z")
ARCHIVED_SPEC = 'evidence/native_geometric_source_20260925/spec.json'
ARCHIVED_SPEC_SHA256 = 'a6b14ad3c89ab17d3c213375022975e1a35239a82810406b73157a6b856cadcc'
ARCHIVED_OWNER = 'code/native_geometric_source/build.py'
ARCHIVED_OWNER_SHA256 = 'cb121da03e48b3968ec7c11e0b521faad9f1e07b5ec13c0c9a0130923db72671'
ARCHIVED_MANIFEST = 'evidence/observer_dynamics_20260925/manifest.json'
ARCHIVED_MANIFEST_SHA256 = '2b4b8090d84e7b4f646e0884c9db10ef9c174393269760c2f7a2552cba17c548'
ARCHIVED_COMMIT = '14d16994b12adb2320b994044bc75db02378b8b1'
ARCHIVED_VENDOR = 'evidence/observer_dynamics_20260925/oph-physics-sim'



@dataclass(frozen=True)
class Owner:
    module: str
    command: str | None
    receipt: str
    authority: str


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(path: Path) -> Any:
    raw = path.read_bytes()
    if path.name.endswith('.gz'):
        raw = gzip.decompress(raw)
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('nonfinite JSON value: ' + value)
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid)


def literal_path(node: ast.AST, values: dict[str, Any]) -> Any:
    """Evaluate only literal path construction; never import or execute a producer."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return values[node.id]
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return Path(literal_path(node.left, values)) / literal_path(node.right, values)
    if isinstance(node, ast.Attribute):
        value = literal_path(node.value, values)
        if isinstance(value, Path) and node.attr in ('parent', 'parents'):
            return value.parent if node.attr == 'parent' else tuple(value.parents)
    if isinstance(node, ast.Subscript):
        return literal_path(node.value, values)[literal_path(node.slice, values)]
    if isinstance(node, ast.Call) and not node.keywords:
        if isinstance(node.func, ast.Name) and node.func.id == 'Path' and len(node.args) == 1:
            return Path(literal_path(node.args[0], values))
        if isinstance(node.func, ast.Attribute):
            value = literal_path(node.func.value, values)
            if isinstance(value, Path):
                if node.func.attr == 'resolve' and not node.args:
                    return value.resolve()
                if node.func.attr == 'with_name' and len(node.args) == 1:
                    return value.with_name(literal_path(node.args[0], values))
    raise ValueError('not a literal path expression')


def discover_owners(root: Path, replay_rows=()) -> tuple[dict[str, Owner], dict[str, list[str]]]:
    owners = {}
    for module, command, receipt in replay_rows:
        owner = Owner(module, command, receipt, 'current replay inventory')
        if receipt in owners and owners[receipt] != owner:
            raise ValueError('conflicting current executable owners: ' + receipt)
        owners[receipt] = owner
    candidates: dict[str, list[str]] = {}
    for path in sorted((root / 'model').rglob('*.py')):
        if 'historical' in path.parts or 'predecessor' in path.parts:
            continue
        text = path.read_text()
        if not all(token in text for token in ('OUTPUT', '__name__', '--write', '--check')):
            continue
        try:
            tree = ast.parse(text)
        except SyntaxError:
            continue
        guards = [node for node in tree.body if isinstance(node, ast.If)
                  and isinstance(node.test, ast.Compare)
                  and isinstance(node.test.left, ast.Name) and node.test.left.id == '__name__'
                  and len(node.test.ops) == len(node.test.comparators) == 1
                  and isinstance(node.test.ops[0], ast.Eq)
                  and isinstance(node.test.comparators[0], ast.Constant)
                  and node.test.comparators[0].value == '__main__']
        writes_output = any(isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name) and node.func.value.id == 'OUTPUT'
            and node.func.attr in ('write_text', 'write_bytes') for node in ast.walk(tree))
        if not guards or not writes_output:
            continue
        values: dict[str, Any] = {'__file__': str(path)}
        for statement in tree.body:
            if not isinstance(statement, (ast.Assign, ast.AnnAssign)):
                continue
            targets = statement.targets if isinstance(statement, ast.Assign) else [statement.target]
            try:
                value = literal_path(statement.value, values)
            except (ValueError, KeyError, TypeError, IndexError):
                continue
            for target in targets:
                if isinstance(target, ast.Name):
                    values[target.id] = value
        output = values.get('OUTPUT')
        if not isinstance(output, Path):
            continue
        output = output.resolve()
        if not output.is_relative_to(root) or not output.name.endswith(('.json', '.json.gz')):
            continue
        relative = str(output.relative_to(root))
        module = '.'.join(path.relative_to(root).with_suffix('').parts)
        candidates.setdefault(relative, []).append(module)
    ambiguous = {}
    for receipt, modules in candidates.items():
        if receipt in owners:
            continue
        if len(modules) == 1:
            owners[receipt] = Owner(modules[0], None, receipt, 'literal OUTPUT and original write/check CLI')
        else:
            ambiguous[receipt] = modules
    return owners, ambiguous


def frozen_path(path: str) -> bool:
    """Recognize historical names for reporting; this grants no exemption."""
    parts = Path(path).parts
    return any(part in ('historical', 'predecessor') for part in parts)


class Inventory:
    def __init__(self, root: Path, source_map_roots=None):
        self.root = root.resolve()
        self.source_map_roots = source_map_roots or {}
        self.digests: dict[Path, str] = {}
        self.resolved: dict[tuple[str, str], tuple[Path | None, str | None]] = {}
        self.maps: dict[str, list[dict[str, Any]]] = {}
        self.receipt_hashes: dict[str, str] = {}
        self.parse_errors: dict[str, str] = {}

    def resolve(self, base: Path, relative: Any):
        key = (str(base), relative) if isinstance(relative, str) else None
        if key in self.resolved:
            return self.resolved[key]
        error = None
        target = None
        if not isinstance(relative, str) or '://' in relative or Path(relative).is_absolute():
            error = 'invalid local source path'
        else:
            candidate = (base / relative).resolve()
            if not candidate.is_relative_to(self.root) or not candidate.is_file():
                error = 'missing or escaping local source path'
            else:
                target = candidate
        result = (target, error)
        if key is not None:
            self.resolved[key] = result
        return result

    def digest(self, target: Path) -> str:
        if target not in self.digests:
            self.digests[target] = sha(target.read_bytes())
        return self.digests[target]

    def read(self, receipt: str) -> list[dict[str, Any]]:
        if receipt in self.maps:
            return self.maps[receipt]
        findings: list[dict[str, Any]] = []
        self.maps[receipt] = findings
        try:
            self.receipt_hashes[receipt] = sha((self.root / receipt).read_bytes())
            value = load_json(self.root / receipt)
        except (OSError, ValueError, EOFError) as error:
            self.parse_errors[receipt] = str(error)
            return findings
        def walk(value, location='root'):
            if isinstance(value, dict):
                for key, child in value.items():
                    here = location + '.' + key
                    recognized = key in MAP_FIELDS or key in PACKAGE_LEVELS or here == ENDPOINT_PINS
                    local = recognized and key != 'sources'
                    records = (key == 'source_bindings' and isinstance(child, list) and bool(child)
                        and all(isinstance(row, dict) and set(row) == {'path', 'sha256'}
                                and isinstance(row['path'], str) for row in child))
                    if records and len({row['path'] for row in child}) != len(child):
                        findings.append({'receipt': receipt, 'location': here, 'error': 'duplicate source-binding record'})
                    if local and not records and (not isinstance(child, dict) or not child):
                        findings.append({'receipt': receipt, 'location': here, 'error': 'empty or malformed local source map'})
                    if recognized and (isinstance(child, dict) or records):
                        base = ((self.root / receipt).parents[PACKAGE_LEVELS[key]]
                                if key in PACKAGE_LEVELS else self.root)
                        base = self.source_map_roots.get((receipt, here), base)
                        entries = child.items() if isinstance(child, dict) else [(row['path'], row['sha256']) for row in child]
                        for relative, expected in entries:
                            target, error = self.resolve(base, relative)
                            if not local and (error or not isinstance(expected, str) or not HEX.fullmatch(expected)):
                                continue  # Mixed external metadata uses its owning verifier.
                            row = {'receipt': receipt, 'location': here, 'source': relative, 'expected': expected}
                            if error:
                                row['error'] = error
                            elif not isinstance(expected, str) or not HEX.fullmatch(expected):
                                row['error'] = 'invalid SHA256 in local source map'
                            else:
                                row['source'] = str(target.relative_to(self.root))
                                row['actual'] = self.digest(target)
                                row['stale'] = row['actual'] != expected
                            findings.append(row)
                    walk(child, here)
            elif isinstance(value, list):
                for index, child in enumerate(value):
                    walk(child, f'{location}[{index}]')
        walk(value)
        return findings


def audit(root: Path, roots: list[str], owners: dict[str, Owner], *, mode='check',
          all_model=True, documentary_exemptions=frozenset(), frozen_captures=None, source_map_roots=None) -> dict[str, Any]:
    root = root.resolve()
    if not roots or mode not in ('check', 'write'):
        raise ValueError('preflight requires nonempty explicit roots and a valid mode')
    roots = list(dict.fromkeys(roots))
    for receipt in roots:
        path = Path(receipt)
        target = (root / path).resolve()
        if path.is_absolute() or not target.is_relative_to(root) or not target.is_file():
            raise ValueError('missing or escaping explicit receipt: ' + receipt)
    frozen_captures = frozen_captures or {}
    for receipt, expected in frozen_captures.items():
        path = (root / receipt).resolve()
        if not path.is_relative_to(root) or sha(path.read_bytes()) != expected:
            raise ValueError('fixed capture bytes changed: ' + receipt)
    inventory = Inventory(root, source_map_roots)
    if all_model:
        for path in sorted((root / 'model').rglob('*')):
            if path.is_file() and path.name.endswith(('.json', '.json.gz')):
                inventory.read(str(path.relative_to(root)))
    graph: dict[str, set[str]] = {}
    unknown = []
    def visit(receipt):
        if receipt in graph:
            return
        graph[receipt] = set()
        if receipt in roots and receipt not in owners:
            unknown.append({'receipt': receipt, 'error': 'selected live receipt has no authenticated owner'})
        for row in inventory.read(receipt):
            child = row.get('source', '')
            if not isinstance(child, str) or not child.endswith(('.json', '.json.gz')):
                continue
            if child in frozen_captures:
                continue
            # Every source-bearing dependency is checked, regardless of its
            # filename or whether a replayable producer has been identified.
            child_maps = inventory.read(child)
            if child in owners or child_maps or child in inventory.parse_errors:
                graph[receipt].add(child)
                visit(child)
    for receipt in roots:
        visit(receipt)
    order: list[str] = []
    visiting: set[str] = set()
    visited: set[str] = set()
    cycles = []
    def sort(receipt):
        if receipt in visited:
            return
        if receipt in visiting:
            cycles.append(receipt)
            return
        visiting.add(receipt)
        for child in sorted(graph[receipt]):
            sort(child)
        visiting.remove(receipt)
        visited.add(receipt)
        order.append(receipt)
    for receipt in sorted(graph):
        sort(receipt)
    failures = list(unknown)
    permitted = []
    adapter = []
    changed = set()
    refreshable = set(roots) - {child for children in graph.values() for child in children}
    for receipt in graph:
        if receipt in inventory.parse_errors:
            failures.append({'receipt': receipt, 'error': inventory.parse_errors[receipt]})
        for row in inventory.read(receipt):
            if not row.get('error') and not row.get('stale'):
                continue
            key = (receipt, row['location'], row.get('source'), row.get('expected'))
            if not row.get('error') and key in documentary_exemptions:
                adapter.append(row)
            elif mode == 'write' and receipt in refreshable and not row.get('error'):
                permitted.append(row)
                changed.add(receipt)
            else:
                failures.append(row)
                changed.add(receipt)
    affected = set(changed)
    for receipt in order:
        if graph[receipt] & affected:
            affected.add(receipt)
    outside = [row for receipt, rows in inventory.maps.items() if receipt not in graph
               for row in rows if row.get('stale') or row.get('error')]
    changed_during_audit = []
    for path, expected in list(inventory.digests.items()) + [(root / p, h) for p, h in inventory.receipt_hashes.items()]:
        try:
            actual = sha(path.read_bytes())
        except OSError:
            actual = None
        if actual != expected:
            changed_during_audit.append(str(path.relative_to(root)))
    return {'schema': 'oph.receipt_source_preflight.v1', 'mode': mode,
        'roots': roots, 'verdict': 'FAIL' if failures or cycles or changed_during_audit else 'PASS',
        'scientific_replay_performed': False,
        'original_writers_colds_and_documentary_adapters_required': True,
        'model_JSON_files_scanned': sum(p.startswith('model/') for p in inventory.maps),
        'all_JSON_files_scanned': len(inventory.maps), 'unique_source_files': len(inventory.digests),
        'live_dependency_receipts': sorted(graph),
        'source_bearing_inputs_without_executable_owner': sorted(set(graph) - set(owners)),
        'authenticated_fixed_captures': frozen_captures,
        'affected_dependencies_without_executable_owner': sorted(affected - set(owners)),
        'live_dependency_edges': {p: sorted(graph[p]) for p in sorted(graph)},
        'failures': failures, 'dependency_cycles': sorted(set(cycles)),
        'files_changed_during_audit': sorted(set(changed_during_audit)),
        'permitted_selected_root_refreshes': permitted,
        'required_original_documentary_adapters': adapter,
        'ordered_affected_owners': [asdict(owners[p]) for p in order if p in affected and p in owners],
        'background_findings': [dict(row, classification=('authenticated_fixed_capture' if row['receipt'] in frozen_captures else ('historical_path_outside_closure' if frozen_path(row['receipt']) else 'outside_live_dependency_closure'))) for row in outside],
        'outside_closure_parse_errors': {p: e for p, e in inventory.parse_errors.items() if p not in graph}}


def current_inventory(root: Path):
    if root != ROOT.resolve():
        raise ValueError('CLI must execute from the repository whose current routing is authenticated')
    from model.verify_current import current_replays
    return current_replays()


def documentary_exemptions(root: Path):
    """Authenticate the one existing retained document pin; do not replay science."""
    from model import native_control_custody as custody
    control = root / custody.CONTROL
    if sha(control.read_bytes()) != custody.CONTROL_SHA256:
        raise ValueError('original documentary control changed')
    value = load_json(control)
    expected = value['source_sha256'][custody.DOCUMENT]
    old = subprocess.check_output(['git', 'show', custody.COMMIT + ':' + custody.DOCUMENT], cwd=root)
    if sha(old) != expected:
        raise ValueError('original documentary input pin is unauthenticated')
    return frozenset({(custody.CONTROL, 'root.source_sha256', custody.DOCUMENT, expected)})


def archived_source_map_roots(root: Path):
    # The original owner resolves these seven source paths against VENDOR in
    # native_modules(), not against RER.  Pin that meaning and check the archived
    # manifest's declared native commit before applying the exact path base.
    for name, expected in ((ARCHIVED_SPEC, ARCHIVED_SPEC_SHA256),
            (ARCHIVED_OWNER, ARCHIVED_OWNER_SHA256),
            (ARCHIVED_MANIFEST, ARCHIVED_MANIFEST_SHA256)):
        if sha((root / name).read_bytes()) != expected:
            raise ValueError('archived source contract changed: ' + name)
    spec = load_json(root / ARCHIVED_SPEC)
    if spec.get('schema') != 'oph.native-geometric-source.spec.v1' or spec.get('native_commit') != ARCHIVED_COMMIT:
        raise ValueError('archived source schema or commit changed')
    manifest = load_json(root / ARCHIVED_MANIFEST)
    for path, digest in spec['source_sha256'].items():
        row = manifest['files']['oph-physics-sim/' + path]
        if row['origin_commit'] != ARCHIVED_COMMIT or row['origin_path'] != path or row['sha256'] != digest:
            raise ValueError('archived source manifest differs: ' + path)
    return {(ARCHIVED_SPEC, 'root.source_sha256'): root / ARCHIVED_VENDOR}


def fixed_captures(root: Path):
    contract_path = root / 'model/protected_readout_numerics/contract.json'
    contract = load_json(contract_path)
    if contract.get('schema') != 'oph.protected_readout_numerics.contract.v1':
        raise ValueError('unrecognized fixed-capture contract')
    return {contract['retained_receipt']: contract['retained_receipt_sha256']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', action='append', help='Explicit live root; repeat for more than one')
    parser.add_argument('--mode', choices=('check', 'write'), default='check')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and (args.output.resolve().is_relative_to(ROOT.resolve())
                        or args.output.is_symlink()
                        or (args.output.exists() and args.output.stat().st_nlink != 1)):
        parser.error('--output must be an external unaliased report path')
    rows = current_inventory(ROOT.resolve())
    owners, ambiguous = discover_owners(ROOT.resolve(), rows)
    roots = args.receipt or [receipt for _, _, receipt in rows]
    if args.mode == 'write' and not args.receipt:
        parser.error('write preflight requires an explicitly named root')
    report = audit(ROOT, roots, owners, mode=args.mode,
                   documentary_exemptions=documentary_exemptions(ROOT), frozen_captures=fixed_captures(ROOT),
                   source_map_roots=archived_source_map_roots(ROOT))
    report['ambiguous_literal_owners'] = ambiguous
    text = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text)
    print(json.dumps({key: report[key] for key in ('verdict', 'mode', 'roots',
        'model_JSON_files_scanned', 'unique_source_files', 'failures',
        'ordered_affected_owners')}, indent=2))
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.path.insert(0, str(ROOT))
    raise SystemExit(main())
