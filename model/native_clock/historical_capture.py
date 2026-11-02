"""Narrow compatibility for one retained host capture, never new timestamps."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

CAPTURE_SHA256 = '8a7590e9cd5aa13c7a703651c1851a46c51dfcffd49b125d072914434890900f'
MANIFEST_SHA256 = '2368b601d1244066d34da17ca5218c250463a91d4c7060ca5e9af4a0516152fb'
SOURCE = 'model/native_clock/host_clock.py'
ORIGINAL_SOURCE_SHA256 = '84906bc4afb3d70e026c6b7115e1169a901c4d185399a3f4e72576b18c47e23a'
REQUIREMENTS = 'model/requirements.txt'
ORIGINAL_REQUIREMENTS_SHA256 = '2c585acfcb6e6b57cf3e822eba9e7822b877ce373eaf598be3045669588447cd'
REVIEWED_REQUIREMENTS_SHA256 = '7abaa47730e8e286d366dfe826151677735f20e92beeb7f197dd30e39bd648d2'
SOURCE_SNAPSHOT = 'model/native_clock/historical/capture_source.py.txt'
REQUIREMENTS_SNAPSHOT = 'model/native_clock/historical/capture_requirements.txt'
ANALYSIS_INPUTS = (SOURCE_SNAPSHOT, REQUIREMENTS_SNAPSHOT,
                   'model/native_clock/historical_capture.py',
                   'model/native_clock/HISTORICAL_CAPTURE_COMPATIBILITY.md')


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _require(value, message):
    if not value:
        raise ValueError('capture source bindings: ' + message)


def timed_source_ast(text):
    """All original syntax except analysis and this exact compatibility import.

    The CLI, imports, policy, capture, timer loop and entire record replay stay
    in the comparison. AST locations/comments do not affect execution.
    """
    tree = ast.parse(text)
    allowed_import = ast.parse('from .historical_capture import source_compatibility, analysis_bindings').body[0]
    tree.body = [node for node in tree.body
                 if not (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == 'analyse')
                 and ast.dump(node, include_attributes=False) != ast.dump(allowed_import, include_attributes=False)]
    return ast.dump(tree, include_attributes=False)


def analysis_bindings(root):
    root = Path(root)
    return {path: _sha((root/path).read_bytes()) for path in ANALYSIS_INPUTS}


def source_compatibility(raw, current, root, canonical):
    """Preserve capture-time provenance and validate today's structural replay.

    Modified records with the same historical source declaration can exercise
    the existing countercontrols, but are explicitly not the retained measured
    capture. Their event/scheduler replay and metadata checks are never skipped.
    """
    root = Path(root)
    captured = raw['capture_sources']
    original_data = _sha(canonical(raw).encode()) == CAPTURE_SHA256
    if captured == current:
        mode = 'exact_current_source_declaration'
    else:
        _require(type(captured) is dict and _sha(json.dumps(captured, sort_keys=True,
                 separators=(',', ':')).encode()) == MANIFEST_SHA256,
                 'not the exact reviewed historical source manifest')
        _require(captured.keys() == current.keys(), 'source path catalogue changed')
        _require(_sha((root/SOURCE_SNAPSHOT).read_bytes()) == ORIGINAL_SOURCE_SHA256,
                 'original capture source snapshot changed')
        _require(_sha((root/REQUIREMENTS_SNAPSHOT).read_bytes()) == ORIGINAL_REQUIREMENTS_SHA256,
                 'original declared requirements snapshot changed')
        _require(captured[REQUIREMENTS] == ORIGINAL_REQUIREMENTS_SHA256
                 and current[REQUIREMENTS] == REVIEWED_REQUIREMENTS_SHA256,
                 'unreviewed requirements transition')
        for path, value in captured.items():
            if path not in (SOURCE, REQUIREMENTS):
                _require(current[path] == value, 'timed/native dependency changed: ' + path)
        _require(timed_source_ast((root/SOURCE).read_text()) ==
                 timed_source_ast((root/SOURCE_SNAPSHOT).read_text()),
                 'capture, scheduler, policy, timer or native replay syntax changed')
        mode = 'specific_historical_capture_with_unchanged_execution_syntax'
    return {'mode': mode, 'original_capture_data_match': original_data,
            'retained_original_capture_sha256': CAPTURE_SHA256,
            'capture_source_manifest_preserved': captured,
            'capture_declared_requirements_sha256': captured[REQUIREMENTS],
            'current_analysis_requirements_sha256': current[REQUIREMENTS],
            'current_timing_measurement_claimed': False,
            'timing_portability_claimed': False,
            'unrecorded_runtime_or_hardware_reconstructed': False,
            'modified_historical_input_is_a_measurement': False}
