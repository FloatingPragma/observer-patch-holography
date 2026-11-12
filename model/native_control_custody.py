"""Authenticate the native control's retained documentary input.

The six executable inputs and complete scientific result remain current.
Only the historical control's named documentary input is restored privately;
the current prose is separately consumed by the unchanged live producer.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile

from tools import strict_json
from tools import check_completion_current_code as views

CONTROL = "model/common_source/native_defect_decision_results.json"
LIVE = "model/common_source/native_defect_decision_live_results.json"
DOCUMENT = "docs/research/OPH_FIELD_RECONSTRUCTION.md"
COMMIT = "6d39bfd75dda28a7affb199f7cb5f70d39feedf5"
CONTROL_SHA256 = "fa53cbaef0f65752edfa385cbb7d15ef2bddada656a29d8ac774844208767bb6"
SOURCE_ERROR = f"Receipt source digest mismatch: {CONTROL} at root.source_sha256[{DOCUMENT}]"
REPLAY = ("model.common_source.native_defect_decision", None, CONTROL)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def authenticate(root, fresh):
    """Check immutable control bytes, every source, and the full live result."""
    root = Path(root).resolve()
    historical = subprocess.check_output(["git", "show", COMMIT + ":" + CONTROL], cwd=root)
    views.require(views.sha(historical) == CONTROL_SHA256
                  and (root / CONTROL).read_bytes() == historical,
                  "retained native control bytes changed")
    retained = strict_json.loads(historical.decode())
    documents = views.authenticated_documents(root, documentary_paths=(DOCUMENT,))
    source = retained["source_sha256"]
    views.require(len(source) == 7 and DOCUMENT in source,
                  "native control executable/documentary inventory changed")
    for name, expected in source.items():
        target = (root / name).resolve()
        views.require(target.is_relative_to(root) and target.is_file(),
                      "missing or escaping native source: " + name)
        data = documents[name] if name == DOCUMENT else target.read_bytes()
        views.require(views.sha(data) == expected, "native executable or documentary input changed: " + name)
    live = strict_json.load(root / LIVE)
    views.require(canonical(live) == canonical(fresh), "live native producer replay changed")
    science = lambda value: {key: child for key, child in value.items() if key != "source_sha256"}
    views.require(canonical(science(live)) == canonical(science(retained)),
                  "complete native scientific payload changed")
    views.require(set(live["source_sha256"]) == set(source), "live native source inventory changed")
    for name, expected in live["source_sha256"].items():
        views.require(views.sha((root / name).read_bytes()) == expected,
                      "live native source custody changed: " + name)
        if name != DOCUMENT:
            views.require(expected == source[name], "native executable source changed: " + name)
    return documents


def verify_in_documentary_view(retained):
    """Run the original artifact/source verifier with one authenticated input."""
    root = retained.ROOT.resolve()
    from model.common_source import native_defect_decision
    documents = authenticate(root, native_defect_decision.report())
    before = views.current_files(root)
    with tempfile.TemporaryDirectory(prefix="oph-native-control-custody-") as directory:
        view = Path(directory) / "view"
        views.install_current_inputs(root, view, before, documents, documentary_paths=(DOCUMENT,))
        try:
            retained.ROOT = view
            count = retained.verify_artifacts()
        finally:
            retained.ROOT = root
        views.verify_view(view, before, documentary_paths=(DOCUMENT,))
    views.require(views.current_files(root) == before, "current files changed during native custody check")
    print("Native documentary custody passed: one authenticated retained document; "
          "all other artifacts and sources current, complete native science replayed", flush=True)
    return count


def replay_live_control(root):
    """Actually recalculate both the live and retained scientific comparisons."""
    from model.common_source import native_defect_decision
    authenticate(root, native_defect_decision.report())
    print("Native control replay passed: unchanged executable sources and complete scientific payload; "
          "one authenticated retained documentary input, with current prose checked separately", flush=True)
