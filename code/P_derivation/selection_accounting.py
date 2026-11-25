#!/usr/bin/env python3
"""Reproduce finite selection accounting, without assigning a null probability.

The stored map outputs and substitution scorecard are inputs, not a new search
or an exhaustive menu.  This records their distinct numerical support classes.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import re
import sys

from paper_math import PaperMathContext, decimal_pi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pin_scope  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = "code/P_derivation/"
CERT = PACKAGE + "runtime/p_interval_contraction_certificate_2026-07-14.json"
TRUNK = PACKAGE + "runtime/p_closure_trunk_current.json"
TARGET = PACKAGE + "codata_2022_alpha_fixture.json"
SCORECARD = PACKAGE + "null_model_scorecard.md"
DEFAULT_OUT = ROOT / PACKAGE / "runtime/selection_accounting.json"
PARENTS = (CERT, TRUNK, TARGET, SCORECARD, PACKAGE + "paper_math.py")


def build(root: Path = ROOT) -> dict:
    cert = json.loads((root / CERT).read_text(encoding="utf-8"))
    trunk = json.loads((root / TRUNK).read_text(encoding="utf-8"))
    target = json.loads((root / TARGET).read_text(encoding="utf-8"))["inverse_fine_structure_constant"]["value"]
    ctx = PaperMathContext(precision=60, su2_cutoff=0, su3_cutoff=0)
    rows = []
    with localcontext() as work:
        work.prec = 70
        roots = ctx._koide_roots(ctx.stage5_vectors["delta"])
        for electron in range(5, 12):
            vector = (electron, 4, 3)
            residual = max(
                abs(((roots[i] / roots[j]) ** 2 /
                     ctx.stage5_vectors["epsilon"] ** (vector[i] - vector[j])).ln())
                for i in range(3) for j in range(i + 1, 3)
            )
            rows.append({"electron_exponent": electron, "max_log_ratio_residual": str(+residual)})
        outputs = []
        for mode, block in cert["modes"].items():
            point = block["fixed_point_point_estimate_display_only"]["alpha_inv"]
            outputs.append({
                "id": mode, "support": "interval_certified_fixed_point",
                "alpha_inv_point_display": point,
                "alpha_inv_enclosure": {k: block["certified_enclosure"]["alpha_inv"][k] for k in ("lo", "hi")},
                "distance_ppm_point_display": str(+(abs(Decimal(point) - Decimal(target)) / Decimal(target) * 1000000)),
                "su2_cutoff": block["su2_cutoff"], "su3_cutoff": block["su3_cutoff"],
                "tail_bounds_included": block["edge_sum_tail_bounds"]["included"],
            })
        fixed = trunk["fixed_point_candidate"]
        point = fixed["alpha_inv"]
        phi = (Decimal(1) + Decimal(5).sqrt()) / 2
        root_from_p = (Decimal(fixed["P"]) - phi) / decimal_pi(90).sqrt()
        inverse_readout = 1 / Decimal(point)
        outputs.append({
            "id": "compressed_asymptotic_trunk", "support": "noncertified_approximate_candidate",
            "alpha_inv_point_display": point, "alpha_inv_enclosure": None,
            "distance_ppm_point_display": str(+(abs(Decimal(point) - Decimal(target)) / Decimal(target) * 1000000)),
            "relative_printed_pair_defect": str(+(abs(root_from_p - inverse_readout) / inverse_readout)),
            "reported_alpha_fixed_point_residual": fixed["alpha_fixed_point_residual"],
        })
    text = (root / SCORECARD).read_text(encoding="utf-8").split("## W3a:", 1)[1].split("\n## ", 1)[0]
    substitutions = []
    for line in text.splitlines():
        fields = re.findall(r"`([^`]*)`", line)
        if len(fields) == 5 and fields[3].startswith("["):
            lo, hi = fields[4].strip("[]").split(", ")
            substitutions.append({"c1": fields[0], "c2": fields[1], "relative_distance_enclosure": [lo, hi]})
    aliases = {"8/5": "1.60"}
    unique = {(aliases.get(row["c1"], row["c1"]), row["c2"]): row for row in substitutions}
    alternatives = [row for key, row in unique.items() if key != ("φ", "π")]
    return {
        "artifact": "oph.alpha_selection_accounting.v1",
        "promotion_allowed": False,
        "inputs": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                   for name in pin_scope.bound(PARENTS)},
        "comparison_target": {"value": target, "role": "compare_only_empirical_input"},
        "stage5_search": {
            "fixed_mu_tau_exponents": [4, 3], "electron_candidates": rows,
            "selected_vector": list(ctx.stage5_vectors["n_e"]),
            "measured_mass_read_by_this_search": False,
            "source_selected_physical_mass_claim": False,
        },
        "historical_outputs": outputs,
        "substitution_accounting": {
            "source_role": "imported_previously_certified_scorecard_rows_not_a_new_root_solve",
            "rows": substitutions, "label_count": len(substitutions),
            "distinct_numeric_pair_count": len(unique), "distinct_alternative_count": len(alternatives),
            "relative_threshold": "2.5e-6",
            "certified_alternative_hits": sum(Decimal(row["relative_distance_enclosure"][1]) <= Decimal("2.5e-6") for row in alternatives),
        },
        "probability_nonidentification": {
            "role": "countermodels_not_inferred_sampling_laws",
            "menu": [row["id"] for row in outputs[:2]],
            "hit_vector_at_relative_threshold_2_5e_6": [False, True],
            "probability_models": [
                {"weights": ["1/2", "1/2"], "hit_probability": "1/2"},
                {"weights": ["1/4", "3/4"], "hit_probability": "3/4"},
            ],
            "conclusion": "Identical output values, range, cardinality and hit set admit different hit probabilities; numerical coverage alone does not determine significance.",
        },
        "interpretation": {
            "recorded_map_outputs": 3, "interval_certified_map_roots": 2,
            "exhaustive_configuration_menu": False, "precomparison_menu": False,
            "configuration_sampling_measure": None, "global_chance_probability": None,
            "claimed_112_configuration_scan_reproduced": False,
            "counts_are_probabilities": False,
            "boundary": "Three historical outputs and a separate retrospective substitution grid do not define a map-form null distribution; the alleged 112-label product supplies neither a named complete menu nor a sampling measure.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build(), indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
