#!/usr/bin/env python3
"""Independent arithmetic, coverage and uncertainty check of the reconciliation.

No producer import, large event replay, covariance reconstruction or simulation.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE_PINS = {'evidence/source_net_causal_poset/source_net_manifold_sampled_q55_q89_2026-09-25.json': '8dba828403e6d6da9ba16dcec28e53a2271bb26ca68682a70dee17f5d6ae23e0', 'evidence/source_net_causal_poset/source_net_manifold_sampled_q144_2026-09-25.json': 'a3ca8a2b743631a75f9d6433ab4d12f0536763698d56ca3fb8e8ec505ba210b4', 'evidence/source_net_causal_poset/manifold_refinement/clock_control_receipt.json': '3c4a57ce11dcdc9b6244c1eede9e903044757c2faa1f55cf91af47c1bb3b6120'}


def require(value, message):
    if not value:
        raise ValueError(message)


def same(actual, expected, where):
    if isinstance(expected, float):
        require(type(actual) in (int, float) and math.isfinite(actual), where + ": finite number")
        require(math.isclose(actual, expected, rel_tol=8e-10, abs_tol=3e-15), where + ": arithmetic")
    else:
        require(actual == expected, where + ": value")


def coefficient(k, d):
    # Log-gamma evaluation independent of the generator's direct products.
    return math.exp(-math.log(k) + (k - 1) * (math.lgamma(d + 1) - math.log(2))
                    + math.lgamma(d / 2) + math.lgamma(d)
                    - math.lgamma(k * d / 2) - math.lgamma((k + 1) * d / 2))


def check_stat(stat, estimate, reference, error=None, kind="analytic_flat", uncertainty="not_retained"):
    expected = {"estimate": estimate, "reference": reference, "signed_deviation": estimate - reference,
                "relative_deviation": (estimate - reference) / reference if reference else None,
                "marginal_standard_error": error,
                "standard_error_over_reference": error / reference if error is not None and reference else None,
                "reference_kind": kind, "uncertainty": uncertainty}
    require(set(stat) == set(expected), "statistic fields")
    for key, value in expected.items():
        same(stat[key], value, key)


def verify(receipt, table_path=None):
    require(set(receipt) == {"schema", "inputs", "interpretation", "spectrum_grid", "continuum_reference_inputs",
                             "regions", "trends", "clock_outcomes", "clock_boundary"}, "receipt field inventory")
    require(receipt["schema"] == "oph.finite-manifold-refinement.v1", "schema")
    require(receipt["inputs"] == SOURCE_PINS, "frozen input pin set")
    sources = {}
    for path, digest in SOURCE_PINS.items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, "input digest " + path)
        sources[path] = json.loads((ROOT / path).read_text())
    expected_regions = {}
    sampled = [path for path in sources if "source_net_manifold_sampled_" in path]
    for path in sampled:
        doc = sources[path]
        require(doc["continuum_references"] == receipt["continuum_reference_inputs"], "all sampled reference inputs")
        require(doc["spectrum_grid"] == receipt["spectrum_grid"], "grid")
        for level in doc["levels"]:
            for family in level["families"]:
                for location in ("regions", "homogeneity"):
                    rows = family["regions"] if location == "regions" else family["homogeneity"]["diamonds"]
                    for index, row in enumerate(rows):
                        identity = (level["q"], family["dimension"], row["region"], tuple(row.get("direction", [])))
                        require(identity not in expected_regions, "unique source region")
                        expected_regions[identity] = (path, location, index, row, doc)
    actual = {}
    for row in receipt["regions"]:
        require(set(row) == {"q", "spatial_dimension", "region", "direction", "layers", "event_count", "draws",
                             "source_file", "source_location", "source_index", "sampled_events_sha256", "partners_sha256",
                             "statistics", "cdf_signed_deviations_by_reference_dimension"}, "region field inventory")
        key = (row["q"], row["spatial_dimension"], row["region"], tuple(row["direction"] or []))
        require(key not in actual, "no duplicate region")
        actual[key] = row
    require(set(actual) == set(expected_regions), "all levels, regions and dimensional controls retained")
    require(len(actual) == 96, "96 regions")
    for key, row in actual.items():
        path, location, index, original, doc = expected_regions[key]
        for name, value in {"source_file": path, "source_location": location, "source_index": index,
                            **{f: original[f] for f in ("layers", "event_count", "draws", "sampled_events_sha256", "partners_sha256")}}.items():
            same(row[name], value, "region " + name)
        dim, n, stats = key[1] + 1, original["event_count"], row["statistics"]
        wanted = set()

        def check(name, *args, **kwargs):
            wanted.add(name)
            require(name in stats, "missing statistic " + name)
            check_stat(stats[name], *args, **kwargs)

        for k in (2, 3, 4):
            count = original["chains"][str(k)]["count"]
            ff = math.factorial(n) // math.factorial(n - k) if n < 100 else math.prod(n - j for j in range(k))
            c, se = count["estimate"], count["standard_error"]
            require(se >= 0, "nonnegative marginal error")
            check(f"C{k}_count", c, ff * coefficient(k, dim), se, uncertainty="reported_marginal_not_reconstructed")
            check(f"C{k}_coefficient", c / ff, coefficient(k, dim), se / ff, uncertainty="reported_marginal_not_reconstructed")
            inverted = original["chains"][str(k)]["inverted_dimension"]
            same(coefficient(k, inverted), c / ff, "dimension inversion")
            check(f"C{k}_dimension", inverted, dim)
        c2, c3, c4 = [original["chains"][str(k)]["count"]["estimate"] for k in (2, 3, 4)]
        check("C3_N_over_C2_squared", n * c3 / c2 / c2, coefficient(3, dim) / coefficient(2, dim)**2)
        check("C4_N2_over_C2_cubed", n * n * c4 / c2 / c2 / c2, coefficient(4, dim) / coefficient(2, dim)**3)
        check("count_volume_coefficient", original["count_volume_coefficient"], 1.0,
              kind="supplied_continuum_volume", uncertainty="deterministic_finite_geometry_not_sampling_error")
        spectrum = original["interval_spectrum"]
        check("interval_mean", spectrum["weighted_mean"], coefficient(3, dim) / coefficient(2, dim))
        check("empty_interval_weight", spectrum["empty_interval_weight"], 0.0,
              kind="continuum_zero_lattice_sensitive_control")
        deviations = row["cdf_signed_deviations_by_reference_dimension"]
        require(set(deviations) == set(doc["continuum_references"]), "all five reference dimensions")
        for d, reference in doc["continuum_references"].items():
            difference = [spectrum["cdf_on_grid"][j] - reference["cdf_on_grid"][j] for j in range(25)]
            require(len(deviations[d]) == 25, "all CDF points")
            for a, e in zip(deviations[d], difference):
                same(a, e, "signed CDF deviation")
            check("cdf_grid_distance_d" + d, max(abs(x) for x in difference), 0.0,
                  kind="sampled_flat_reference_25_point_grid", uncertainty="source_and_reference_joint_uncertainty_unavailable")
            check("interval_mean_sampled_reference_d" + d, spectrum["weighted_mean"], reference["mean"],
                  kind="sampled_flat_reference", uncertainty="source_and_reference_joint_uncertainty_unavailable")
            check("interval_second_moment_d" + d, spectrum["weighted_second_moment"], reference["second_moment"],
                  kind="sampled_flat_reference", uncertainty="source_and_reference_joint_uncertainty_unavailable")
        require(set(stats) == wanted, "complete statistic inventory")
    expected_trends = {}
    for q, spatial, region, direction in actual:
        if q != 55:
            continue
        rkey = region + (":" + ",".join(map(str, direction)) if direction else "")
        metrics = [f"C{k}_{s}" for k in (2, 3, 4) for s in ("coefficient", "dimension")]
        metrics += ["C3_N_over_C2_squared", "C4_N2_over_C2_cubed", "count_volume_coefficient", "interval_mean", f"cdf_grid_distance_d{spatial + 1}"]
        for metric in metrics:
            expected_trends[(spatial, rkey, metric)] = [actual[(level, spatial, region, direction)]["statistics"][metric]["signed_deviation"] for level in (55, 89, 144)]
    seen = set()
    for trend in receipt["trends"]:
        require(set(trend) == {"spatial_dimension", "region_key", "statistic", "levels", "signed_deviations",
                               "absolute_deviation_ratios_next_over_previous", "strictly_shrinks_at_both_steps",
                               "statistical_trend_decision"}, "trend field inventory")
        key = (trend["spatial_dimension"], trend["region_key"], trend["statistic"])
        require(key not in seen and key in expected_trends, "unique expected trend")
        seen.add(key)
        deviations = expected_trends[key]
        same(trend["levels"], [55, 89, 144], "trend levels")
        same(trend["signed_deviations"], deviations, "trend deviations")
        ratios = [abs(b / a) if a else None for a, b in zip(deviations, deviations[1:])]
        same(trend["absolute_deviation_ratios_next_over_previous"], ratios, "trend ratios")
        same(trend["strictly_shrinks_at_both_steps"], abs(deviations[0]) > abs(deviations[1]) > abs(deviations[2]), "descriptive trend")
        require(trend["statistical_trend_decision"] == "not_identified_from_retained_summaries", "no invented trend significance")
    require(seen == set(expected_trends), "all region trends")
    require(any(not t["strictly_shrinks_at_both_steps"] for t in receipt["trends"]), "descriptive negative witness")
    policy = receipt["interpretation"]
    require(set(policy) == {"verdict", "decision_scope", "joint_chain_covariance", "weighted_mean_and_cdf_uncertainty",
                            "raw_large_run_draws", "q144_chunk_parallel_driver", "physical_population_read_law_clock_selected",
                            "continuum_manifoldlikeness_established", "new_sampling_performed", "future_uncertainty_contract"},
            "interpretation field inventory")
    for name in ("physical_population_read_law_clock_selected", "continuum_manifoldlikeness_established", "new_sampling_performed"):
        require(policy[name] is False, "no promotion " + name)
    require(policy["verdict"] == "RETAINED_POINT_ESTIMATES_DO_NOT_UNIFORMLY_IMPROVE", "bounded verdict")
    require(policy["decision_scope"] == "Descriptive comparison of three retained finite levels, not statistical rejection of convergence.", "descriptive scope")
    for name in ("joint_chain_covariance", "weighted_mean_and_cdf_uncertainty", "raw_large_run_draws"):
        require(policy[name] == "not_retained", "missing evidence " + name)
    require(policy["q144_chunk_parallel_driver"] == "not_retained; neighbour digest and recorded q21/q34 identities only", "driver custody")
    require(policy["future_uncertainty_contract"] == {
        "bounded_control": "Fix one region and one scalar statistic before drawing; validate on an exactly enumerable order.",
        "joint_chains": "Retain per-draw C2/C3/C4 summands, their cross-products, sample count and sampling design.",
        "weighted_spectrum": "Retain each weight and weighted value/grid indicator, or joint sufficient statistics including denominator cross-products.",
        "comparison": "Include sampled-reference uncertainty, cross-level coupling if any and multiplicity for the chosen decision family.",
        "custody": "Pin the actual driver, graph construction, seeds and all outcomes; declare precision and decision thresholds before sampling.",
    }, "future measurement contract")
    clock_path = next(path for path in sources if path.endswith("/clock_control_receipt.json"))
    clock_fields = ("q", "law", "beta_target", "direction", "actual_beta", "status", "rank_relative_error", "volume_relative_error", "corrected_volume_relative_error")
    clocks = [{field: r.get(field) for field in clock_fields} for r in sources[clock_path]["rows"]]
    require(receipt["clock_outcomes"] == clocks and len(clocks) == 120, "all 120 clock outcomes and failures")
    require(receipt["clock_boundary"] == "All scheduled outcomes from the distinct q8/q13/q21/q34 clock control are retained; these are not extra q55/q89/q144 manifold samples.", "clock sample scope")
    if table_path is not None:
        with Path(table_path).open(newline="") as stream:
            lines = list(csv.DictReader(stream))
        expected_table = {}
        for row in receipt["regions"]:
            for name, statistic in row["statistics"].items():
                entry = {field: row[field] for field in ("q", "spatial_dimension", "region", "direction", "layers", "event_count", "draws")}
                entry["direction"] = json.dumps(row["direction"], separators=(",", ":"))
                entry.update(statistic= name, **statistic)
                entry = {k: "" if v is None else str(v) for k, v in entry.items()}
                identity = tuple(entry[field] for field in ("q", "spatial_dimension", "region", "direction", "statistic"))
                expected_table[identity] = entry
        require(len(lines) == len(expected_table), "complete CSV table")
        found = set()
        for entry in lines:
            identity = tuple(entry[field] for field in ("q", "spatial_dimension", "region", "direction", "statistic"))
            require(identity not in found and entry == expected_table.get(identity), "CSV row identity and values")
            found.add(identity)
    return {"regions": len(actual), "trends": len(seen), "clock_outcomes": len(clocks)}


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "receipt.json"
    print("FINITE_MANIFOLD_REFINEMENT_VERIFIED", json.dumps(verify(json.loads(path.read_text()), HERE / "refinement.csv"), sort_keys=True))
