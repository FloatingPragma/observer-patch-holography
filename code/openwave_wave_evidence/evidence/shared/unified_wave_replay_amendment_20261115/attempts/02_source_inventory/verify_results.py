#!/usr/bin/env python3
"""Recheck saved replay artifacts and optionally compare multiple environments."""
import argparse
from pathlib import Path
import tempfile
import replay as r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    r.require(not args.output.exists(), "Choose a new verification output")
    registration = r.verify_registration()
    policy = r.strict_json((r.HERE / "POLICY.json").read_bytes())
    result = dict(status="pass", checked_utc=r.utc(), runs={}, cross_environment={})
    with tempfile.TemporaryDirectory(prefix="oph_wave_verify_") as directory:
        bundle = r.extract_git(Path(directory), policy)
        result["custody"] = r.custody(bundle, policy)
        for index, run in enumerate(args.runs):
            report = r.strict_json((run / "REPLAY.json").read_bytes())
            r.require(report["status"] == "pass" and report["source_unchanged"] is True, f"Failed replay: {run}")
            r.require(report["registration_sha256"] == registration and report["policy_sha256"] == r.sha(r.HERE / "POLICY.json"), "Replay input binding mismatch")
            r.require(set(report["leaves"]) == set(r.LEAVES), "Missing/extra replay leaf")
            r.require(report["total_native_calls_in_overlapping_lineage"] == 1009734, "Wrong total calls")
            summaries = {}
            for name, (relative, _, _, _, _, _) in r.LEAVES.items():
                leaf = run / name
                calculated = r.compare_artifacts(name, bundle / relative, leaf, policy)
                process = r.strict_json((leaf / "PROCESS.json").read_bytes())
                for stream in ("stdout", "stderr"):
                    r.require(r.sha(leaf / (stream + ".txt")) == process[stream + "_sha256"], "Log digest mismatch")
                r.require(process["exit_code"] == 0 and process["environment"]["optimize"] == 0, "Unverified/optimized leaf")
                r.require(r.strict_json((leaf / "stderr.txt").read_bytes().splitlines()[0]) == process["environment"], "Environment log mismatch")
                calculated["process"] = process
                r.require(calculated == report["leaves"][name], f"Saved comparison mismatch: {run}/{name}")
                summaries[name] = {k: v for k, v in calculated.items() if k not in ("amended_diagnostics", "process")}
                if index:
                    key = f"0_vs_{index}/{name}"
                    result["cross_environment"][key] = r.compare_artifacts(name, args.runs[0] / name, leaf, policy, amend_reference=False)
            result["runs"][str(index)] = dict(path=str(run), report_sha256=r.sha(run / "REPLAY.json"), leaves=summaries)
    r.save(args.output, result)
    print(f"Verified {len(args.runs)} replay(s), {len(result['cross_environment'])} cross-environment leaf comparisons")


if __name__ == "__main__":
    main()
