#!/usr/bin/env python3
"""Independently check selection accounting; do not import its producer.

Root certificates and the substitution scorecard remain imported evidence.
Their owning interval verifiers establish those numerical enclosures.  This
verifier checks custody, selection arithmetic, distinctions and non-promotion.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys

from mpmath.ctx_mp import MPContext

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pin_scope  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
BASE = "code/P_derivation/"
CERT = BASE + "runtime/p_interval_contraction_certificate_2026-07-14.json"
TRUNK = BASE + "runtime/p_closure_trunk_current.json"
TARGET = BASE + "codata_2022_alpha_fixture.json"
SCORECARD = BASE + "null_model_scorecard.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify(data: dict, root: Path = ROOT) -> None:
    require(set(data) == {"artifact", "promotion_allowed", "inputs", "comparison_target", "stage5_search", "historical_outputs", "substitution_accounting", "probability_nonidentification", "interpretation"}, "receipt schema changed")
    require(data["artifact"] == "oph.alpha_selection_accounting.v1" and data["promotion_allowed"] is False, "selection evidence promoted")
    names = set(pin_scope.bound([CERT, TRUNK, TARGET, SCORECARD, BASE + "paper_math.py"]))
    pins = pin_scope.bound_pins(data["inputs"])
    require(isinstance(pins, dict) and set(pins) == names, "source inventory changed")
    for name in names:
        require(pins[name] == hashlib.sha256((root / name).read_bytes()).hexdigest(), f"source digest changed: {name}")
    cert = json.loads((root / CERT).read_text(encoding="utf-8"))
    trunk = json.loads((root / TRUNK).read_text(encoding="utf-8"))
    target_fixture = json.loads((root / TARGET).read_text(encoding="utf-8"))
    target_text = target_fixture["inverse_fine_structure_constant"]["value"]
    require(target_fixture["claim_status"] == "compare_only_empirical_input", "target promoted")
    require(data["comparison_target"] == {"value": target_text, "role": "compare_only_empirical_input"}, "target changed")
    target = Fraction(target_text)
    mp = MPContext()
    mp.dps = 85
    scan = data["stage5_search"]
    require(set(scan) == {"fixed_mu_tau_exponents", "electron_candidates", "selected_vector", "measured_mass_read_by_this_search", "source_selected_physical_mass_claim"}, "search schema changed")
    require(scan["fixed_mu_tau_exponents"] == [4, 3] and scan["measured_mass_read_by_this_search"] is False and scan["source_selected_physical_mass_claim"] is False, "search premise or physical status changed")
    values = sorted((1 + mp.sqrt(2) * mp.cos(mp.mpf(2) / 9 + 2 * mp.pi * k / 3)) ** 2 for k in range(3))
    expected = {}
    for n in range(5, 12):
        vector = (n, 4, 3)
        expected[n] = max(abs(mp.log(values[i] / values[j]) + (vector[i] - vector[j]) * mp.log(6)) for i in range(3) for j in range(i + 1, 3))
    require([row["electron_exponent"] for row in scan["electron_candidates"]] == list(range(5, 12)), "electron search menu incomplete")
    for row in scan["electron_candidates"]:
        require(set(row) == {"electron_exponent", "max_log_ratio_residual"}, "search row schema changed")
        require(abs(mp.mpf(row["max_log_ratio_residual"]) - expected[row["electron_exponent"]]) < mp.mpf("1e-55"), "independent exponent residual disagrees")
    require(scan["selected_vector"] == [min(expected, key=expected.get), 4, 3] == [7, 4, 3], "wrong selected exponent")
    outputs = data["historical_outputs"]
    modes = ["thomson_structured_running", "thomson_structured_running_plus_gauge_width"]
    require([row["id"] for row in outputs] == modes + ["compressed_asymptotic_trunk"], "map inventory changed")
    require(set(cert["modes"]) == set(modes) and cert["promotion_allowed"] is False, "certificate scope changed")
    for row in outputs:
        point = Fraction(row["alpha_inv_point_display"])
        exact_ppm = abs(point - target) / target * 1000000
        require(abs(Fraction(row["distance_ppm_point_display"]) - exact_ppm) < Fraction("1e-60"), "point distance disagrees")
        if row["id"] in modes:
            require(set(row) == {"id", "support", "alpha_inv_point_display", "alpha_inv_enclosure", "distance_ppm_point_display", "su2_cutoff", "su3_cutoff", "tail_bounds_included"}, "certified row schema changed")
            block = cert["modes"][row["id"]]
            require(row["support"] == "interval_certified_fixed_point", "root support changed")
            require(row["alpha_inv_point_display"] == block["fixed_point_point_estimate_display_only"]["alpha_inv"], "root display changed")
            require(row["alpha_inv_enclosure"] == {k: block["certified_enclosure"]["alpha_inv"][k] for k in ("lo", "hi")}, "root enclosure changed")
            require(Fraction(row["alpha_inv_enclosure"]["lo"]) <= point <= Fraction(row["alpha_inv_enclosure"]["hi"]), "point not enclosed")
            require(all(block["banach"][k] is True for k in ("contraction", "existence", "uniqueness_in_interval", "g_maps_interval_into_interior")), "missing Banach certificate")
            require(row["su2_cutoff"] == block["su2_cutoff"] == 120 and row["su3_cutoff"] == block["su3_cutoff"] == 90, "cutoff changed")
            require(row["tail_bounds_included"] is True and block["edge_sum_tail_bounds"]["included"] is True, "tail control absent")
        else:
            require(set(row) == {"id", "support", "alpha_inv_point_display", "alpha_inv_enclosure", "distance_ppm_point_display", "relative_printed_pair_defect", "reported_alpha_fixed_point_residual"}, "trunk row schema changed")
            require(row["support"] == "noncertified_approximate_candidate" and row["alpha_inv_enclosure"] is None, "approximate trunk falsely certified")
            require(trunk["claim_status"] == "compressed_candidate_trunk_not_final_particle_root" and trunk["consumer_policy"]["may_feed_live_particle_predictions"] is False, "trunk promoted")
            fixed = trunk["fixed_point_candidate"]
            require(row["alpha_inv_point_display"] == fixed["alpha_inv"] and row["reported_alpha_fixed_point_residual"] == fixed["alpha_fixed_point_residual"], "trunk values changed")
            inverse = 1 / mp.mpf(fixed["alpha_inv"])
            reported_residual = mp.mpf(fixed["alpha_fixed_point_residual"])
            require(abs((inverse - mp.mpf(fixed["alpha"])) - reported_residual) < mp.mpf("1e-28"), "trunk reported residual inconsistent")
            for key, exact in (("phi", (1 + mp.sqrt(5)) / 2), ("sqrt_pi", mp.sqrt(mp.pi))):
                printed = trunk["closed_form_candidate"][key]
                require(re.fullmatch(r"[0-9]+\.[0-9]{28,}", printed) is not None, "trunk constant display precision changed")
                last_place = mp.power(10, -len(printed.split(".")[1]))
                require(abs(mp.mpf(printed) - exact) <= last_place, "trunk displayed constant incorrect")
            from_p = (mp.mpf(fixed["P"]) - (1 + mp.sqrt(5)) / 2) / mp.sqrt(mp.pi)
            defect = abs(from_p - inverse) / inverse
            require(defect > mp.mpf("1e-6") and abs(mp.mpf(row["relative_printed_pair_defect"]) - defect) < mp.mpf("1e-60"), "trunk pair discrepancy hidden")
    accounting = data["substitution_accounting"]
    require(set(accounting) == {"source_role", "rows", "label_count", "distinct_numeric_pair_count", "distinct_alternative_count", "relative_threshold", "certified_alternative_hits"}, "substitution schema changed")
    require(accounting["source_role"] == "imported_previously_certified_scorecard_rows_not_a_new_root_solve", "substitution evidence scope changed")
    section = (root / SCORECARD).read_text(encoding="utf-8").split("## W3a:", 1)[1].split("\n## ", 1)[0]
    expected_rows = []
    for line in section.splitlines():
        cells = [part.strip().strip("`") for part in line.split("|")[1:-1]]
        if len(cells) == 6 and cells[4].startswith("["):
            expected_rows.append({"c1": cells[0], "c2": cells[1], "relative_distance_enclosure": cells[4][1:-1].split(", ")})
    require(accounting["rows"] == expected_rows and len(expected_rows) == accounting["label_count"] == 48, "substitution rows incomplete")
    distinct = {}
    for row in expected_rows:
        key = ("8/5" if row["c1"] == "1.60" else row["c1"], row["c2"])
        lo, hi = map(Fraction, row["relative_distance_enclosure"])
        require(0 <= lo <= hi, "reversed distance enclosure")
        if key in distinct:
            require(distinct[key] == (lo, hi), "numeric aliases disagree")
        distinct[key] = (lo, hi)
    require(len(distinct) == accounting["distinct_numeric_pair_count"] == 42, "alias accounting changed")
    alternatives = {key: bounds for key, bounds in distinct.items() if key != ("φ", "π")}
    require(len(alternatives) == accounting["distinct_alternative_count"] == 41, "canonical pair included as alternative")
    require(accounting["relative_threshold"] == "2.5e-6", "retrospective threshold changed")
    threshold = Fraction(accounting["relative_threshold"])
    require(all(hi <= threshold or lo > threshold for lo, hi in alternatives.values()), "threshold-straddling interval")
    require(sum(hi <= threshold for lo, hi in alternatives.values()) == accounting["certified_alternative_hits"] == 0, "alternative hit count changed")
    counter = data["probability_nonidentification"]
    require(set(counter) == {"role", "menu", "hit_vector_at_relative_threshold_2_5e_6", "probability_models", "conclusion"}, "countermodel schema changed")
    require(counter["role"] == "countermodels_not_inferred_sampling_laws" and counter["menu"] == modes, "countermodels mistaken for a selected law")
    hits = []
    for row in outputs[:2]:
        lo, hi = (Fraction(row["alpha_inv_enclosure"][key]) for key in ("lo", "hi"))
        require(hi < target, "counterexample interval location changed")
        low_distance, high_distance = (target - hi) / target, (target - lo) / target
        require(high_distance <= threshold or low_distance > threshold, "counterexample hit unresolved")
        hits.append(high_distance <= threshold)
    require(counter["hit_vector_at_relative_threshold_2_5e_6"] == hits == [False, True], "counterexample hit set changed")
    probabilities = []
    require(len(counter["probability_models"]) == 2, "two countermodels required")
    for model in counter["probability_models"]:
        require(set(model) == {"weights", "hit_probability"} and len(model["weights"]) == len(hits), "countermodel weights missing")
        weights = list(map(Fraction, model["weights"]))
        require(all(weight > 0 for weight in weights) and sum(weights) == 1, "countermodel measure invalid")
        probability = sum(weight for weight, hit in zip(weights, hits) if hit)
        require(probability == Fraction(model["hit_probability"]), "countermodel probability incorrect")
        probabilities.append(probability)
    require(probabilities[0] != probabilities[1], "countermodels do not show nonidentification")
    require(counter["conclusion"] == "Identical output values, range, cardinality and hit set admit different hit probabilities; numerical coverage alone does not determine significance.", "countermodel scope changed")
    interpretation = data["interpretation"]
    require(set(interpretation) == {"recorded_map_outputs", "interval_certified_map_roots", "exhaustive_configuration_menu", "precomparison_menu", "configuration_sampling_measure", "global_chance_probability", "claimed_112_configuration_scan_reproduced", "counts_are_probabilities", "boundary"}, "interpretation schema changed")
    require(interpretation["recorded_map_outputs"] == 3 and interpretation["interval_certified_map_roots"] == 2, "map support count changed")
    for key in ("exhaustive_configuration_menu", "precomparison_menu", "claimed_112_configuration_scan_reproduced", "counts_are_probabilities"):
        require(interpretation[key] is False, f"unsupported {key}")
    require(interpretation["configuration_sampling_measure"] is None and interpretation["global_chance_probability"] is None, "unsupported null probability")
    require(interpretation["boundary"] == "Three historical outputs and a separate retrospective substitution grid do not define a map-form null distribution; the alleged 112-label product supplies neither a named complete menu nor a sampling measure.", "statistical boundary changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", type=Path, default=ROOT / BASE / "runtime/selection_accounting.json")
    args = parser.parse_args()
    verify(json.loads(args.path.read_text(encoding="utf-8")))
    print("SELECTION_ACCOUNTING_VALID: 7 exponent candidates; 2 certified roots + 1 approximate output; 48 substitution labels / 42 distinct pairs; no null probability")


if __name__ == "__main__":
    main()
