"""Assemble the paired receipt and byte inventory from completed runs."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

if __package__:
    from .experiment import ARCHIVE, HERE, ROOT, save_json, sha
else:
    from experiment import ARCHIVE, HERE, ROOT, save_json, sha


def wiring_census():
    rows = []
    for level in (3, 4, 5):
        g = dict(np.load(HERE / f"geometry/geometry_l{level}.npz"))
        c = len(g["faces"])
        base, centres = g["vertices"][:12], g["centres"]
        antipode = np.argmin(base@base.T, axis=1)
        reference = np.tile([0., 0., 1.], (c, 1))
        reference[abs(centres[:, 2]) > .9] = [1., 0., 0.]
        tangent = np.cross(reference, centres)
        tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
        frame = np.stack([tangent, np.cross(centres, tangent), centres], axis=2)
        incident = np.bincount(g["faces"].ravel(), minlength=len(g["vertices"]))
        defect = (incident[g["faces"]] == 5).any(axis=1)
        for wiring in ("w12", "w3"):
            w = g[wiring]
            degree = np.bincount(w[:, [0, 2]].ravel(), minlength=c)
            a = np.einsum("nij,nj->ni", frame[w[:, 0]], base[w[:, 1]])
            b = np.einsum("nij,nj->ni", frame[w[:, 2]], base[w[:, 3]])
            cosine = np.sum(a*b, axis=1)
            rows.append({"level": level, "wiring": wiring, "carriers": c, "glued_seams": len(w),
                         "degree_histogram": {str(k): v for k, v in sorted(Counter(map(int, degree)).items())},
                         "unused_ports": int(np.sum(12-degree)), "defect_cells": int(defect.sum()),
                         "unused_on_defect_cells": int(np.sum(12-degree[defect])),
                         "antipodal_label_pairs": int(np.sum(antipode[w[:, 1]] == w[:, 3])),
                         "global_direction_dot": {"min": float(cosine.min()), "mean": float(cosine.mean()), "max": float(cosine.max())}})
    return rows


def historical():
    directory = ROOT / "evidence/exact_federation_L6_canonical_20260909"
    def read(name):
        return json.loads((directory / name).read_text())
    kernel = read("kernel_readout.json")
    float_w3 = read("float_law_port_pair.json")
    float_isolated = read("float_law_isolated.json")
    integer = read("integer_law_port_pair.json")
    return {"archive": directory.relative_to(ROOT).as_posix(), "level": 6, "carriers": 81920,
            "kernel_sample_count": len(kernel["cells"]), "kernel_steps": kernel["steps"],
            "w3_kernel_summary": kernel["summary_glued"],
            "isolated_kernel_readout": kernel["per_cell"][0]["isolated"],
            "isolated_n300_max_deviation_from_4P": kernel["isolated_max_abs_deviation_from_4_P_slow_at_n_300"],
            "canonical_w3": {k: float_w3[k] for k in ("law", "sweeps", "attempts", "terminated", "budget", "max_abs_deviation_from_component_mean")},
            "canonical_isolated": {"law": float_isolated["law"], "schedules": len(float_isolated["schedules"]),
                                   "sweeps": float_isolated["sweeps"], "all_terminated": float_isolated["all_terminated"],
                                   "component_mean_within_1e-9_all": float_isolated["component_mean_within_1e-9_all"]},
            "historical_integer_control": {k: integer[k] for k in ("law", "sweeps", "unique_quotient_hash_count", "quotient_hash_equals_expected_all")},
            "sha256": {p.name: sha(p) for p in sorted(directory.iterdir()) if p.is_file()},
            "comparison_limit": "Historical L6 uses a different regulator, preparation, schedule and budget; matched L3-L5 controls are supplied separately. Integer nearest agreement is never the provenance law."}


def assemble(output, mirror=True):
    read = lambda name: json.loads((output/name).read_text())
    trace = read("trace/trace.json")
    receipt = {"schema": "oph.support_wiring.paired_receipt.v1", "issue": 776,
               "specification_sha256": sha(HERE/"SPECIFICATION.md"), "geometry": json.loads((HERE/"geometry/geometry.json").read_text()),
               "wiring": wiring_census(), "canonical_controls": read("controls.json"),
               "provenance_intervals": read("provenance.json"), "record_metric_q13": read("q13.json"),
               "record_metric_q13_controls": read("q13_controls.json"),
               "historical_L6": historical(),
               "execution": {"law": trace["binding"]["law"], "events": trace["events"], "phases": len(trace["chunks"]),
                             "mean_actions": sum(c["event_count"]//2 for c in trace["chunks"] if c["kind"] in ("intra", "glued")),
                             "glued_actions": sum(c["event_count"]//2 for c in trace["chunks"] if c["kind"] == "glued"),
                             "preparation_writes": sum(c["event_count"] for c in trace["chunks"] if c["kind"] == "prepare"),
                             "refinement_writes": sum(c["event_count"] for c in trace["chunks"] if c["kind"] == "copy"),
                             "final_denominator_exponent": trace["chunks"][-1]["exponent"], "final_chain": trace["final_chain"]},
               "claim_boundary": {"declared_wiring_port_assignment_schedule_loads_and_joins": True,
                                  "canonical_law_for_every_provenance_seam": True,
                                  "physical_clock": False, "continuum_limit": False, "source_selection": False,
                                  "M1_derived": False, "universal_finite_schedule_confluence": False,
                                  "dimension_is_acceptance_target": False,
                                  "port_label_antipodes_imposed_on_gluing": False,
                                  "join_is_commutative_cell_pullback_in_common_port_labels": True},
               "interpretation": {"confluence": "Each component has one fixed point at its preserved mean. Repetition of a finite sweep containing every seam converges to it. The four-sweep endpoints remain schedule dependent; no finite exact termination is claimed.",
                                  "kernels": "Expectation-operator all-port probes, not the response of the particular phase history. Four declared cells per glued wiring; extrema cover only those samples.",
                                  "positions": "Each interval has both coordinate placements. Coordinates select anchors and annotate events; the authenticated order and its counts do not change when a fixed interval is replotted.",
                                  "growth": "Raw interval counts against phase gaps and longest-chain height, plus successive ratios. Neither is a physical clock or a fitted continuum dimension.",
                                  "paired_gap": "The two routes use different laws, event granularities, preparations and clocks. The measured contrast does not isolate wiring as the sole cause, and is not a derivation of the record-metric read law.",
                                  "q13_boundary": "The lag-four q13 interval is clipped. Its approximately 4.15 estimate must not be labelled an interior-diamond result. Lags one through three are retained with their interior flags."}}
    save_json(output/"support_wiring_receipt.json", receipt)
    if mirror:
        destination = ROOT/"evidence/source_net_causal_poset/support_wiring_receipt.json"
        destination.write_bytes((output/"support_wiring_receipt.json").read_bytes())
        path = destination.parent/"archive_manifest.json"
        old = json.loads(path.read_text())
        old["inventory"] = [r for r in old["inventory"] if r["path"] != destination.name]
        old["inventory"].append({"path": destination.name, "bytes": destination.stat().st_size, "sha256": sha(destination)})
        old["inventory"].sort(key=lambda r: r["path"])
        old["curated_archive"]["file_count"] = len(old["inventory"])
        old["curated_archive"]["total_bytes"] = sum(row["bytes"] for row in old["inventory"])
        inventory_text = "".join(f"{row['sha256']}  {row['bytes']}  {row['path']}\n" for row in old["inventory"])
        old["curated_archive"]["inventory_sha256"] = hashlib.sha256(inventory_text.encode()).hexdigest()
        selection = old["curated_archive"]["selection"]
        suffix = "; paired canonical support-wiring diagnostic (#776)"
        if not selection.endswith(suffix):
            old["curated_archive"]["selection"] = selection+suffix
        old["result"] = {("central_diamond_myrheim_meyer_dimension_3d" if key == "interior_diamond_myrheim_meyer_dimension_3d" else key): value
                         for key, value in old["result"].items()}
        source_family = json.loads((destination.parent/"source_net_causal_limit_receipt.json").read_text())
        old["result"]["central_diamond_inside_cube_3d"] = [
            next(f for f in level["families"] if f["dimension"] == 3)["vertical_intervals"][-1]["continuum_diamond_inside_cube"]
            for level in source_family["levels"]]
        # Preserve every existing historical pin; only append the new attachment.
        path.write_text(json.dumps(old, indent=1)+"\n", encoding="utf-8", newline="\n")
    excluded = {"archive_manifest.json"}
    files = [p for p in output.rglob("*") if p.is_file() and p.name not in excluded]
    source_files = [p for p in HERE.rglob("*") if p.is_file() and "__pycache__" not in p.parts]
    workflow = ROOT/".github/workflows/support-wiring.yml"
    if workflow.exists():
        source_files.append(workflow)
    inventory = {"schema": "oph.support_wiring.archive.v1",
                 "artifacts": {p.relative_to(output).as_posix(): {"sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(files)},
                 "sources": {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(source_files)},
                 "mirror": "evidence/source_net_causal_poset/support_wiring_receipt.json"}
    save_json(output/"archive_manifest.json", inventory)
    print(f"Assembled paired receipt; {len(files)} evidence artifacts, {sum(p.stat().st_size for p in files):,} bytes", flush=True)
    return receipt


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=ARCHIVE)
    p.add_argument("--no-mirror", action="store_true")
    args = p.parse_args()
    assemble(args.output, not args.no_mirror)
