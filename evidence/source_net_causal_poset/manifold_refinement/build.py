#!/usr/bin/env python3
"""Reconcile retained finite-order diagnostics; never execute a source-net run."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ROOT = PACKAGE.parents[1]
INPUTS = (
    "evidence/source_net_causal_poset/source_net_manifold_sampled_q55_q89_2026-09-25.json",
    "evidence/source_net_causal_poset/source_net_manifold_sampled_q144_2026-09-25.json",
    "evidence/source_net_causal_poset/manifold_refinement/clock_control_receipt.json",
)
FIELDS = ["q", "spatial_dimension", "region", "direction", "layers", "event_count", "draws",
          "statistic", "estimate", "reference", "signed_deviation", "relative_deviation",
          "marginal_standard_error", "standard_error_over_reference", "reference_kind", "uncertainty"]
POLICY = {
    "verdict": "RETAINED_POINT_ESTIMATES_DO_NOT_UNIFORMLY_IMPROVE",
    "decision_scope": "Descriptive comparison of three retained finite levels, not statistical rejection of convergence.",
    "joint_chain_covariance": "not_retained",
    "weighted_mean_and_cdf_uncertainty": "not_retained",
    "raw_large_run_draws": "not_retained",
    "q144_chunk_parallel_driver": "not_retained; neighbour digest and recorded q21/q34 identities only",
    "physical_population_read_law_clock_selected": False,
    "continuum_manifoldlikeness_established": False,
    "new_sampling_performed": False,
    "future_uncertainty_contract": {
        "bounded_control": "Fix one region and one scalar statistic before drawing; validate on an exactly enumerable order.",
        "joint_chains": "Retain per-draw C2/C3/C4 summands, their cross-products, sample count and sampling design.",
        "weighted_spectrum": "Retain each weight and weighted value/grid indicator, or joint sufficient statistics including denominator cross-products.",
        "comparison": "Include sampled-reference uncertainty, cross-level coupling if any and multiplicity for the chosen decision family.",
        "custody": "Pin the actual driver, graph construction, seeds and all outcomes; declare precision and decision thresholds before sampling.",
    },
}


def chi(k, d):
    return (math.gamma(d + 1) / 2) ** (k - 1) * math.gamma(d / 2) * math.gamma(d) / (
        k * math.gamma(k * d / 2) * math.gamma((k + 1) * d / 2))


def metric(estimate, reference, error=None, kind="analytic_flat", uncertainty="not_retained"):
    return {"estimate": estimate, "reference": reference,
            "signed_deviation": estimate - reference,
            "relative_deviation": estimate / reference - 1 if reference else None,
            "marginal_standard_error": error,
            "standard_error_over_reference": error / reference if error is not None and reference else None,
            "reference_kind": kind, "uncertainty": uncertainty}


def region_key(row):
    return row["region"] + (":" + ",".join(map(str, row["direction"])) if row.get("direction") else "")


def build():
    documents = [json.loads((ROOT / name).read_text()) for name in INPUTS]
    regions = []
    for source, doc in zip(INPUTS[:2], documents[:2]):
        for level in doc["levels"]:
            for family in level["families"]:
                dim = family["dimension"]
                for location in ("regions", "homogeneity"):
                    originals = family["regions"] if location == "regions" else family["homogeneity"]["diamonds"]
                    for index, row in enumerate(originals):
                        n = row["event_count"]
                        stats = {}
                        for k in range(2, 5):
                            ch = row["chains"][str(k)]
                            falling = math.prod(range(n - k + 1, n + 1))
                            estimate, se = ch["count"]["estimate"], ch["count"]["standard_error"]
                            stats[f"C{k}_count"] = metric(estimate, falling * chi(k, dim + 1), se,
                                                          uncertainty="reported_marginal_not_reconstructed")
                            stats[f"C{k}_coefficient"] = metric(estimate / falling, chi(k, dim + 1), se / falling,
                                                                uncertainty="reported_marginal_not_reconstructed")
                            stats[f"C{k}_dimension"] = metric(ch["inverted_dimension"], dim + 1)
                        c2, c3, c4 = [row["chains"][str(k)]["count"]["estimate"] for k in (2, 3, 4)]
                        stats["C3_N_over_C2_squared"] = metric(c3 * n / c2**2, chi(3, dim + 1) / chi(2, dim + 1)**2)
                        stats["C4_N2_over_C2_cubed"] = metric(c4 * n**2 / c2**3, chi(4, dim + 1) / chi(2, dim + 1)**3)
                        stats["count_volume_coefficient"] = metric(row["count_volume_coefficient"], 1.0,
                            kind="supplied_continuum_volume", uncertainty="deterministic_finite_geometry_not_sampling_error")
                        spectrum = row["interval_spectrum"]
                        stats["interval_mean"] = metric(spectrum["weighted_mean"], chi(3, dim + 1) / chi(2, dim + 1))
                        stats["empty_interval_weight"] = metric(spectrum["empty_interval_weight"], 0.0,
                            kind="continuum_zero_lattice_sensitive_control")
                        cdf = {}
                        for ref_dim, ref in doc["continuum_references"].items():
                            signed = [a - b for a, b in zip(spectrum["cdf_on_grid"], ref["cdf_on_grid"])]
                            cdf[ref_dim] = signed
                            stats[f"cdf_grid_distance_d{ref_dim}"] = metric(max(map(abs, signed)), 0.0,
                                kind="sampled_flat_reference_25_point_grid", uncertainty="source_and_reference_joint_uncertainty_unavailable")
                            stats[f"interval_mean_sampled_reference_d{ref_dim}"] = metric(spectrum["weighted_mean"], ref["mean"],
                                kind="sampled_flat_reference", uncertainty="source_and_reference_joint_uncertainty_unavailable")
                            stats[f"interval_second_moment_d{ref_dim}"] = metric(spectrum["weighted_second_moment"], ref["second_moment"],
                                kind="sampled_flat_reference", uncertainty="source_and_reference_joint_uncertainty_unavailable")
                        regions.append({"q": level["q"], "spatial_dimension": dim, "region": row["region"],
                            "direction": row.get("direction"), "layers": row["layers"], "event_count": n, "draws": row["draws"],
                            "source_file": source, "source_location": location, "source_index": index,
                            "sampled_events_sha256": row["sampled_events_sha256"], "partners_sha256": row["partners_sha256"],
                            "statistics": stats, "cdf_signed_deviations_by_reference_dimension": cdf})
    regions.sort(key=lambda r: (r["spatial_dimension"], region_key(r), r["q"]))
    trends = []
    groups = {}
    for row in regions:
        groups.setdefault((row["spatial_dimension"], region_key(row)), []).append(row)
    for (dim, key), rows in sorted(groups.items()):
        # Raw Ck counts grow with volume; compare normalized deviations instead.
        names = [f"C{k}_{kind}" for k in (2, 3, 4) for kind in ("coefficient", "dimension")]
        names += ["C3_N_over_C2_squared", "C4_N2_over_C2_cubed", "count_volume_coefficient", "interval_mean",
                  f"cdf_grid_distance_d{dim + 1}"]
        for name in names:
            values = [row["statistics"][name]["signed_deviation"] for row in rows]
            absolute = list(map(abs, values))
            trends.append({"spatial_dimension": dim, "region_key": key, "statistic": name,
                "levels": [row["q"] for row in rows], "signed_deviations": values,
                "absolute_deviation_ratios_next_over_previous": [b / a if a else None for a, b in zip(absolute, absolute[1:])],
                "strictly_shrinks_at_both_steps": all(b < a for a, b in zip(absolute, absolute[1:])),
                "statistical_trend_decision": "not_identified_from_retained_summaries"})
    clock_fields = ("q", "law", "beta_target", "direction", "actual_beta", "status", "rank_relative_error",
                    "volume_relative_error", "corrected_volume_relative_error")
    clocks = [{key: row.get(key) for key in clock_fields} for row in documents[2]["rows"]]
    return {"schema": "oph.finite-manifold-refinement.v1", "inputs": {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in INPUTS},
        "interpretation": POLICY, "spectrum_grid": documents[0]["spectrum_grid"],
        "continuum_reference_inputs": documents[0]["continuum_references"],
        "regions": regions, "trends": trends, "clock_outcomes": clocks,
        "clock_boundary": "All scheduled outcomes from the distinct q8/q13/q21/q34 clock control are retained; these are not extra q55/q89/q144 manifold samples."}


def table(receipt):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in receipt["regions"]:
        metadata = {k: row[k] for k in FIELDS[:7]}
        metadata["direction"] = json.dumps(row["direction"], separators=(",", ":"))
        for name, stat in sorted(row["statistics"].items()):
            writer.writerow({**metadata, "statistic": name, **stat})
    return output.getvalue()


if __name__ == "__main__":
    receipt = build()
    (HERE / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (HERE / "refinement.csv").write_text(table(receipt))
    print(f"RECONCILED {len(receipt['regions'])} regions; {len(receipt['trends'])} descriptive trends; {len(receipt['clock_outcomes'])} clock outcomes")
