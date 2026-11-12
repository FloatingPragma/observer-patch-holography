"""Execute one cold current module with authenticated finite documentary inputs.

The parent owns and verifies the complete current-worktree snapshot. This
entry point checks every file read through the two scoped roots against the
current source, except their two exact authenticated documentary inputs.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import runpy
import sys

from tools import check_completion_current_code as views
from tools.completion_contract_inputs import bind_consumers, restore_consumers


def authenticate_view(root, view):
    root, view = Path(root).resolve(), Path(view).resolve()
    views.require(view != root and view.is_dir(), "separate finite input view required")
    documents = views.authenticated_documents(root)
    campaign = "model/completion/campaign.json"
    clock = "model/clock_defect_energy/campaign.json"
    control = "model/common_source/native_defect_decision_results.json"
    family = json.loads((root / campaign).read_text())
    diagnostic = json.loads((root / clock).read_text())
    native = json.loads((root / control).read_text())
    names = {campaign, clock, control}
    names.update(row["path"] for group in ("source_bindings", "retained_controls") for row in family[group])
    names.update(diagnostic["source_sha256"])
    names.update(native["source_sha256"])
    for name in names:
        current, copied = (root / name).resolve(), (view / name).resolve()
        views.require(current.is_relative_to(root) and copied.is_relative_to(view), "escaping scoped input")
        expected = documents[name] if name in documents else current.read_bytes()
        views.require(copied.read_bytes() == expected, "scoped scientific or documentary input changed: " + name)


def execute(module, arguments, view, *, root=views.ROOT):
    views.require(module.startswith("model.") and all(part.isidentifier() for part in module.split(".")),
                  "explicit current model module required")
    authenticate_view(root, view)
    if module == "model.replay_native_fixed":
        # This unchanged wrapper launches another cold Python process. Its
        # children must naturally read the same authenticated input view;
        # execve also clears all helper imports from their runtime inventory.
        scoped_view = Path(view).resolve()
        os.chdir(scoped_view)
        os.execve(sys.executable, [sys.executable, "-m", module, *arguments],
                  views.environment(scoped_view))
    from model.completion import family
    from model.clock_defect_energy import contract
    modules = (family, contract)
    scoped_view = Path(view).resolve()
    original = bind_consumers(modules, scoped_view)
    argv = sys.argv
    try:
        # Original validators execute with every scientific assertion active.
        family.load_contract()
        contract.load_contract()
        sys.argv = [module, *arguments]
        runpy.run_module(module, run_name="__main__", alter_sys=True)
    finally:
        sys.argv = argv
        restore_consumers(modules, original, scoped_view)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--view", type=Path, required=True)
    parser.add_argument("--module", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
    execute(args.module, arguments, args.view)


if __name__ == "__main__":
    main()
