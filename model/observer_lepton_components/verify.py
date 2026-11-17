#!/usr/bin/env python3
"""Verify the workspace component register; no scientific experiment is run."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from completion_replay import verify_packet

ROOT = Path(__file__).resolve().parents[2]
WORKSPACE = ROOT.parent
LEDGER = ROOT / "claims/openwave_rows_11_20_components.json"
OUTPUT = Path(__file__).with_name("results.json")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_path(relative):
    path = (WORKSPACE / relative).resolve()
    if not path.is_relative_to(WORKSPACE.resolve()) or not path.is_file():
        raise ValueError("Invalid or missing workspace artifact: " + relative)
    return path


def verify(ledger=LEDGER):
    data = json.loads(ledger.read_text())
    if data["physical_rows_established"] != [] or data["official_awards_added"] != 0:
        raise ValueError("Component register cannot award physical rows")
    leaves = 0
    for packet in data["sealed_packets"]:
        manifest = source_path(packet["manifest"])
        if digest(manifest) != packet["manifest_sha256"]:
            raise ValueError("Changed sealed manifest: " + packet["manifest"])
        for line in manifest.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            expected, name = line.split(None, 1)
            path = manifest.parent / name.strip().lstrip("*")
            if not path.is_file():
                path = source_path(name.strip().lstrip("*"))
            if not path.resolve().is_relative_to(WORKSPACE.resolve()) or digest(path) != expected:
                raise ValueError("Changed sealed artifact: " + str(path))
            leaves += 1
    snapshots = 0
    for packet in data["unsealed_work"]:
        pause = source_path(packet["snapshot"])
        if digest(pause) != packet["snapshot_sha256"]:
            raise ValueError("Changed draft snapshot: " + packet["snapshot"])
        snapshot = json.loads(pause.read_text())
        if snapshot["final_review_complete"] or snapshot["physical_row_established"]:
            raise ValueError("Unreviewed draft promoted in snapshot")
        for name, expected in snapshot["files_sha256"].items():
            if digest(pause.parent / name) != expected:
                raise ValueError("Changed paused draft: " + str(pause.parent / name))
            snapshots += 1
    reviews = data.get("component_reviews", [])
    if len(reviews) != 1:
        raise ValueError("Expected the separately scoped component review")
    review = reviews[0]
    if review["same_group"] is not True or review["full_row_solved"] is not False or review["G6_confirmed"] is not False:
        raise ValueError("Component review scope changed")
    portable = source_path(review["portable_manifest"])
    if digest(portable) != review["portable_manifest_sha256"]:
        raise ValueError("Portable review binding changed")
    completion = verify_packet(portable.parent)
    original = portable.parent / review["review_manifest_inside_portable"]
    if digest(original) != review["manifest_sha256"]:
        raise ValueError("Original review seal changed")
    for relative, expected in data.get("integration_sources", {}).items():
        if digest(source_path(relative)) != expected:
            raise ValueError("Changed integration source: " + relative)
    return {
        "schema": "oph.observer_lepton_components.custody.v2",
        "register_sha256": digest(ledger),
        "sealed_manifests_verified": len(data["sealed_packets"]),
        "sealed_file_hashes_verified": leaves,
        "unsealed_snapshots_verified": len(data["unsealed_work"]),
        "unsealed_file_hashes_verified": snapshots,
        "same_group_component_reviews_verified": len(reviews),
        "portable_completion": completion,
        "integration_sources_verified": len(data.get("integration_sources", {})),
        "scientific_experiments_executed_by_this_verifier": 0,
        "physical_rows_established": [],
        "official_awards_added": 0,
        "all_pass": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    text = json.dumps(result, indent=2) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != text:
            raise SystemExit("Component custody receipt is absent or stale")
    else:
        OUTPUT.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
