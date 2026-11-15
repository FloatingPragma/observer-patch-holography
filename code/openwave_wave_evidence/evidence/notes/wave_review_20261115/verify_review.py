#!/usr/bin/env python3
"""Check this review's receipts and comparator; optionally seal its artifacts."""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path

from reproduce import HERE, compare_json, save, sha, verify_manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seal", action="store_true")
    parser.add_argument("--check-seal", action="store_true")
    args = parser.parse_args()
    checks = []

    def check(condition, name):
        if not condition:
            raise ValueError(name)
        checks.append(name)

    for path in HERE.rglob("*.py"):
        ast.parse(path.read_text(encoding="utf-8"))
    checks.append("all review Python sources parse")
    base = {"integer": 2, "flag": True, "values": [1.5], "hash": "a" * 64}
    check(compare_json(base, base)["exact_mismatch_count"] == 0, "comparator identical baseline")
    for label, changed in [
        ("nonfinite NaN", dict(base, values=[math.nan])),
        ("nonfinite infinity", dict(base, values=[math.inf])),
        ("integer change", dict(base, integer=3)),
        ("boolean/integer type", dict(base, flag=1)),
        ("missing sequence element", dict(base, values=[])),
        ("extra field", dict(base, added=0)),
        ("hash change", dict(base, hash="b" * 64)),
    ]:
        check(compare_json(base, changed)["exact_mismatch_count"] > 0, "comparator rejects " + label)
    check(compare_json(base, dict(base, values=[1.6]))["above_1e_minus_10"] == 1,
          "comparator reports changed numerical value")
    independent = json.loads((HERE / "results/INDEPENDENT.json").read_bytes())
    check(independent["script_sha256"] == sha(HERE / "independent.py"), "independent source hash")
    check(independent["all_review_checks_pass"], "independent arithmetic receipt")
    check(independent["free_tapes"]["actual_current_records"] == 131072, "independent actual tape count")
    check(independent["free_tapes"]["branches"] == 37908, "independent branch count")
    check(independent["dispersion"]["sampled_cases"] == 2000, "dispersion sample count")
    profiles = json.loads((HERE / "results/PROFILE_DIAGNOSTICS.json").read_bytes())
    check(profiles["script_sha256"] == sha(HERE / "diagnose_profiles.py"), "profile diagnostic source hash")
    for label in ("windows", "linux_numpy235"):
        directory = HERE / "results" / label
        receipt_path = directory / "REPRODUCTION.json"
        receipt = json.loads(receipt_path.read_bytes())
        prefix = label + ": "
        check(receipt["subject_commit"] == independent["subject_commit"], prefix + "same commit")
        check(receipt["outer_manifest_entries"] == receipt["outer_manifest_entries_after"] == 420,
              prefix + "outer manifest coverage")
        check(receipt["copied_pins_verified"] == 409 and receipt["subject_unchanged"],
              prefix + "copied leaf coverage and unchanged subject")
        check(len(receipt["ancestor_manifests"]) == 6, prefix + "six ancestor manifests")
        check(receipt["custody"]["exit_code"] == 0, prefix + "exact-byte custody passed")
        check(receipt["original_replay"]["exit_code"] == 1, prefix + "original replay failure retained")
        check(receipt["optimized_original_replay"]["exit_code"] == 0, prefix + "optimization bypass reproduced")
        check(not all(receipt["optimized_original_replay"]["json"]["byte_replay"].values()),
              prefix + "optimized replay explicitly contains failed comparisons")
        check(receipt["origin_mutation_normal"]["exit_code"] == 1 and
              receipt["origin_mutation_optimized"]["exit_code"] == 0,
              prefix + "origin mutation rejected normally and accepted optimized")
        check(sha(directory / "REPRODUCE_SOURCE.py") == receipt["reviewer_script_sha256"],
              prefix + "exact executed reviewer source retained")
        commands = {name: receipt[name] for name in (
            "environment", "custody", "original_replay", "optimized_original_replay",
            "origin_mutation_normal", "origin_mutation_optimized")}
        for name, leaf in receipt["leaf_diagnostics"].items():
            check(leaf["exit_code"] == 0, prefix + name + " normal scientific leaf passed")
            commands["leaf_" + name] = leaf
        for name, command in commands.items():
            for stream in ("stdout", "stderr"):
                path = directory / (name + "." + stream + ".txt")
                check(sha(path) == command[stream + "_sha256"], prefix + name + " " + stream + " bound")
        diagnostic = profiles["environments"][label]
        check(diagnostic["reproduction_sha256"] == sha(receipt_path), prefix + "profile report binding")
        check(not diagnostic["all_large_differences_are_sign_only"], prefix + "magnitude changes retained")
        for profile in diagnostic["profiles"]:
            check(profile["expected_internal_reference"]["max_readout_error"] < 1e-12 and
                  profile["actual_internal_reference"]["max_readout_error"] < 1e-12,
                  prefix + "reference-port reconstruction " + str((profile["population"], profile["band"])))
    if args.check_seal:
        print(json.dumps(dict(review_receipt_checks=len(checks),
                              manifest_entries=verify_manifest(HERE, HERE / "MANIFEST.sha256"))))
        return
    save(HERE / "results/VALIDATION.json",
         dict(all_review_consistency_checks_pass=True, checks=checks,
              scope="Validation of review receipts/tools, not a passed original bundle or physical row."))
    if args.seal:
        files = sorted(p for p in HERE.rglob("*") if p.is_file()
                       and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        content = "".join(sha(p) + "  " + p.relative_to(HERE).as_posix() + "\n" for p in files)
        (HERE / "MANIFEST.sha256").write_bytes(content.encode())
        count = verify_manifest(HERE, HERE / "MANIFEST.sha256")
    else:
        count = None
    print(json.dumps(dict(review_receipt_checks=len(checks), sealed_files=count), indent=2))


if __name__ == "__main__":
    main()
