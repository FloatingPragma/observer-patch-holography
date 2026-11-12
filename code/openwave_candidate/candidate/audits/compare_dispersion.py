#!/usr/bin/env python3
"""Audit two already-saved definition-only recomputations and a pinned reference.

This post-seal comparison is not a new blind calculation, external certification,
or physical-data test. It compares full exact polynomials and independently
re-evaluates the saved numerical rows with their analytic remainder bounds.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import mpmath as mp
import sympy as sp


PROJECT = Path(__file__).resolve().parents[2]
PUBLIC_BUNDLE = PROJECT.name == "openwave_candidate" and PROJECT.parent.name == "code"
REPOSITORY = PROJECT.parents[1] if PUBLIC_BUNDLE else PROJECT.parent
DEFAULT_RER_REPOSITORY = REPOSITORY if PUBLIC_BUNDLE else REPOSITORY/"reverse-engineering-reality"
RER_COMMIT = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"
RER_RECEIPT = "code/a5_fingerprint/runtime/spin_six_primitive_port_prediction_receipt.json"
x, y, z, a, t = sp.symbols("x y z a t", real=True)
XYZ = (x, y, z)
SYMBOLS = {"x": x, "y": y, "z": z, "a": a, "t": t, "sqrt": sp.sqrt}


def expr(value):
    return sp.sympify(value, locals=SYMBOLS)


def same(left, right):
    return sp.simplify(sp.expand(expr(left)-expr(right))) == 0


def relative(path):
    return os.path.relpath(Path(path).resolve(), REPOSITORY)


def digest(content):
    return hashlib.sha256(content).hexdigest()


def polynomial_from_map(coefficients):
    return sp.expand(sum(expr(value)*sp.prod(variable**power for variable, power
        in zip(XYZ, map(int, key.split(",")))) for key, value in coefficients.items()))


def compare_polynomial_maps(first, second, second_multiplier=1):
    keys = sorted(set(first) | set(second))
    differences = {key: str(sp.simplify(expr(first.get(key, "0"))
                    -second_multiplier*expr(second.get(key, "0")))) for key in keys}
    mismatches = {key: value for key, value in differences.items() if value != "0"}
    return {"passed": not mismatches, "monomial_count": len(keys),
            "mismatching_coefficient_differences": mismatches}


def exact_comparison(first, second):
    aa = first["outputs"]
    comparisons = {}
    maps = {}
    for degree in (2, 4, 6):
        label = str(degree)
        am = aa["moment_polynomials"][label]
        bm = second["moments"][label]
        ah = aa["omega_homogeneous_polynomials"][label]
        bh = second["homogeneous_polynomials"][label]
        scale = a**(degree-2)
        comparisons[f"moment_{degree}_expression"] = same(am["expression"], bm["polynomial"])
        maps[f"moment_{degree}_all_coefficients"] = compare_polynomial_maps(am["monomials"], bm["coefficients"])
        comparisons[f"moment_{degree}_A_map_consistency"] = same(am["expression"], polynomial_from_map(am["monomials"]))
        comparisons[f"moment_{degree}_B_map_consistency"] = same(bm["polynomial"], polynomial_from_map(bm["coefficients"]))
        comparisons[f"symbol_{degree}_expression"] = same(ah["expression"], scale*expr(bh["polynomial"]))
        maps[f"symbol_{degree}_all_coefficients"] = compare_polynomial_maps(ah["monomials"], bh["coefficients"], scale)
        comparisons[f"symbol_{degree}_A_map_consistency"] = same(ah["expression"], polynomial_from_map(ah["monomials"]))
        comparisons[f"symbol_{degree}_B_map_consistency"] = same(bh["polynomial"], polynomial_from_map(bh["coefficients"]))
        comparisons[f"rotational_part_{degree}"] = same(aa["omega_rotational_parts"][label], scale*expr(bh["rotationally_invariant_polynomial"]))
        multiplier = sp.Rational((-1)**(degree//2+1), 2*sp.factorial(degree))
        comparisons[f"spherical_mean_{degree}"] = same(multiplier*expr(aa["moment_sphere_averages"][label]), bh["rotationally_invariant_coefficient"])
        comparisons[f"nonconstant_part_{degree}"] = same(expr(ah["expression"])-expr(aa["omega_rotational_parts"][label]), scale*expr(bh["nonconstant_angular_polynomial"]))
    ai = aa["I6_homogeneous"]
    bi = second["angular_normalization"]
    comparisons["I6_expression"] = same(ai["expression"], bi["I6_homogeneous_polynomial"])
    maps["I6_all_coefficients"] = compare_polynomial_maps(ai["monomials"], bi["I6_homogeneous_coefficients"])
    comparisons["I6_A_map_consistency"] = same(ai["expression"], polynomial_from_map(ai["monomials"]))
    comparisons["I6_B_map_consistency"] = same(bi["I6_homogeneous_polynomial"], polynomial_from_map(bi["I6_homogeneous_coefficients"]))
    comparisons["P6_polynomial"] = same(expr(aa["legendre_P6"]).subs(t, x), bi["legendre_P6"])
    comparisons["legendre_sum_homogeneous"] = same(aa["legendre_sum_homogeneous"]["expression"], bi["homogeneous_sum_P6"])
    comparisons["vertex_normalization"] = same(aa["I6_vertex_denominator"], bi["sum_P6_at_vertex"])
    harmonic_b = expr(second["moments"]["6"]["polynomial"])-expr(aa["moment_sphere_averages"]["6"])*(x*x+y*y+z*z)**3
    comparisons["harmonic_sextic_moment"] = same(aa["harmonic_sextic_moment"]["expression"], harmonic_b)
    for source_key, target_key in (("c2", "c2"), ("c4", "c4"), ("c6_isotropic", "c6_isotropic"), ("c6_angular", "c6_I6")):
        comparisons[f"normalized_{source_key}"] = same(aa["coefficients"][source_key], second["normalized_coefficients"][target_key])
    for source_key, target_key in (
            ("c6_angular/c6_isotropic", "c6_I6_over_c6_isotropic"),
            ("c2*c6_isotropic/c4^2", "c2_times_c6_isotropic_over_c4_squared"),
            ("c2*c6_angular/c4^2", "c2_times_c6_I6_over_c4_squared")):
        comparisons[f"ratio_{source_key}"] = same(aa["scale_free_ratios"][source_key], second["scale_free_ratios"][target_key])
    for source_key, target_key in (
            ("quadratic_multiplier", "quadratic_term_multiplier"),
            ("quartic_multiplier", "quartic_term_multiplier"),
            ("sextic_multiplier", "sextic_isotropic_term_multiplier"),
            ("sextic_multiplier", "sextic_I6_term_multiplier"),
            ("remainder_bound_multiplier", "eighth_order_remainder_bound_multiplier")):
        comparisons[f"length_halving_{target_key}"] = same(aa["a_to_a_over_2_at_fixed_k_and_n"][source_key], second["length_halving"][target_key])
    tensor = sp.Matrix([[expr(value) for value in row] for row in second["broken_control"]["broken_quadratic_tensor"]])
    broken_quadratic = (sp.Matrix(XYZ).T*tensor*sp.Matrix(XYZ))[0]
    comparisons["removed_port_quadratic"] = same(aa["broken_control_quadratics"]["one_port_removed"]["expression"], broken_quadratic)
    comparisons.update({label: result["passed"] for label, result in maps.items()})
    return {"passed": all(comparisons.values()), "checks": comparisons,
            "coefficient_map_comparisons": maps,
            "total_compared_polynomial_coefficient_entries": sum(row["monomial_count"] for row in maps.values()),
            "convention_adjustment": "A includes a^(degree-2) in each omega polynomial; B factors it outside. The comparison restores this factor symbolically."}


def numerical_audit(first, second):
    """Recompute every saved row, rather than merely trust receipt check flags."""
    mp.mp.dps = 130
    phi = (1+mp.sqrt(5))/2
    radius = mp.sqrt(1+phi*phi)
    raw = [(0, i, j*phi) for i in (-1, 1) for j in (-1, 1)]
    raw += [(i, j*phi, 0) for i in (-1, 1) for j in (-1, 1)]
    raw += [(j*phi, 0, i) for i in (-1, 1) for j in (-1, 1)]
    raw_directions = {"axis": (1, 0, 0), "axis_x": (1, 0, 0),
        "diagonal": (1, 1, 1), "body_diagonal": (1, 1, 1),
        "generic": (1, 2, 3), "generic_1_2_3": (1, 2, 3),
        "port": (0, 1, phi), "defined_vertex": (0, 1, phi),
        "generic_2_minus3_5": (2, -3, 5)}
    polynomials = {degree: sp.lambdify(XYZ, expr(second["homogeneous_polynomials"][str(degree)]["polynomial"]), "mpmath") for degree in (2, 4, 6)}
    row_sets = {"A": first["outputs"]["direct_symbol"]["rows"],
                "B": second["direct_symbol_check"]["rows"]}
    audit = {}

    def near(saved, recomputed, digits):
        return abs(mp.mpf(saved)-recomputed) <= mp.mpf(10)**(-digits)*max(abs(recomputed), mp.mpf("1e-110"))

    for label, rows in row_sets.items():
        checks = {"saved_symbols_reproduced": True, "saved_truncations_reproduced": True,
            "saved_errors_reproduced": True, "saved_directional_bounds_reproduced": True,
            "analytic_absolute_bounds_hold": True, "cosine_sine_identity": True,
            "all_arguments_below_one": True, "saved_ratios_reproduced": True,
            "at_least_200fold_decay_on_halving_k": True}
        if label == "A":
            checks.update({"next_bounds_reproduced": True, "alternating_signed_interval_holds": True})
        else:
            checks.update({"uniform_bounds_reproduced": True, "directional_bound_below_uniform": True})
        groups = {}
        maximum_argument = mp.mpf(0)
        maximum_ratio = mp.mpf(0)
        minimum_decay = mp.inf
        for row in rows:
            direction = raw_directions[row["direction"]]
            norm = mp.sqrt(sum(v*v for v in direction))
            n = tuple(mp.mpf(v)/norm for v in direction)
            dots = [sum(u*q for u, q in zip(port, n))/radius for port in raw]
            length, momentum = mp.mpf(row["a"]), mp.mpf(row["k"])
            arguments = [length*momentum*value for value in dots]
            direct = sum(1-mp.cos(value) for value in arguments)/(2*length**2)
            sine = sum(mp.sin(value/2)**2 for value in arguments)/length**2
            truncation = sum(length**(degree-2)*momentum**degree*polynomials[degree](*n) for degree in (2, 4, 6))
            error = direct-truncation
            b8 = length**6*momentum**8*sum(value**8 for value in dots)/(2*mp.factorial(8))
            b10 = length**8*momentum**10*sum(value**10 for value in dots)/(2*mp.factorial(10))
            uniform = length**6*momentum**8/mp.mpf(6720)
            max_argument = max(abs(value) for value in arguments)
            maximum_argument = max(maximum_argument, max_argument)
            maximum_ratio = max(maximum_ratio, abs(error)/b8)
            digits = 32 if label == "A" else 42
            symbol_key = "symbol" if label == "A" else "direct_cosine_symbol"
            truncation_key = "truncation" if label == "A" else "degree_six_truncation"
            error_key = "signed_error" if label == "A" else "absolute_error"
            bound_key = "absolute_bound_B8" if label == "A" else "directional_remainder_bound"
            ratio_key = "abs_error_over_B8" if label == "A" else "error_over_directional_bound"
            checks["saved_symbols_reproduced"] &= near(row[symbol_key], direct, 60 if label == "A" else 42)
            checks["saved_truncations_reproduced"] &= near(row[truncation_key], truncation, 60 if label == "A" else 42)
            checks["saved_errors_reproduced"] &= near(row[error_key], error if label == "A" else abs(error), digits)
            checks["saved_directional_bounds_reproduced"] &= near(row[bound_key], b8, digits)
            checks["saved_ratios_reproduced"] &= near(row[ratio_key], abs(error)/b8, digits)
            checks["analytic_absolute_bounds_hold"] &= abs(error) <= b8+mp.mpf("1e-120")
            checks["cosine_sine_identity"] &= abs(direct-sine) < mp.mpf("1e-120")
            checks["all_arguments_below_one"] &= max_argument < 1
            if label == "A":
                checks["next_bounds_reproduced"] &= near(row["next_bound_B10"], b10, digits)
                checks["alternating_signed_interval_holds"] &= -b8-mp.mpf("1e-120") <= error <= -b8+b10+mp.mpf("1e-120")
            else:
                checks["uniform_bounds_reproduced"] &= near(row["uniform_remainder_bound"], uniform, digits)
                checks["directional_bound_below_uniform"] &= b8 <= uniform+mp.mpf("1e-120")
            groups.setdefault((row["direction"], row["a"]), []).append((momentum, abs(error)))
        for group in groups.values():
            group.sort(reverse=True)
            for previous, current in zip(group, group[1:]):
                if previous[0] == 2*current[0]:
                    decay = previous[1]/current[1]
                    minimum_decay = min(minimum_decay, decay)
                    checks["at_least_200fold_decay_on_halving_k"] &= decay > 200
        audit[label] = {"rows_rechecked": len(rows), "checks": {key: bool(value) for key, value in checks.items()},
            "maximum_absolute_argument": mp.nstr(maximum_argument, 30),
            "maximum_error_over_B8": mp.nstr(maximum_ratio, 30),
            "minimum_error_decay_factor_on_halving_k": mp.nstr(minimum_decay, 30),
            "passed": all(checks.values())}
    return {"decimal_precision": mp.mp.dps, "receipts": audit,
        "passed": all(item["passed"] for item in audit.values()),
        "absolute_bound_interpretation": "Taylor's degree-seven remainder gives |error| <= sum |a*k*u.n|^8/(2*a^2*8!) for all real arguments; |u.n| <= 1 yields a^6*k^8/6720.",
        "A_signed_refinement_interpretation": "For these arguments, whose magnitudes are below one, the alternating cosine series gives -B8 <= error <= -B8+B10. This extra refinement is not needed for B's valid absolute bound.",
        "precision_note": "Analytic bounds are mathematical statements; archived decimal values are finite-precision computations, independently reproduced here at 130 digits with tolerances respecting their stored significant digits."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=Path, default=PROJECT/"certification/recompute_a.json")
    parser.add_argument("--b", type=Path, default=PROJECT/"certification/recompute_b.json")
    parser.add_argument("--packet", type=Path, default=PROJECT/"certification/DEFINITION_ONLY.md")
    parser.add_argument("--rer-repo", type=Path, default=DEFAULT_RER_REPOSITORY)
    parser.add_argument("--rer-commit", default=RER_COMMIT)
    parser.add_argument("--output", type=Path, default=PROJECT/"certification/comparison.json")
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.rer_commit):
        parser.error("--rer-commit must be a full forty-character commit hash")
    files = {"recompute_a": args.a, "recompute_b": args.b, "definition_packet": args.packet,
             "implementation_a": PROJECT/"candidate/oph_openwave/dispersion.py",
             "implementation_b": PROJECT/"candidate/audits/dispersion_oracle.py",
             "comparison_implementation": Path(__file__).resolve()}
    contents = {key: path.read_bytes() for key, path in files.items()}
    first, second = json.loads(contents["recompute_a"]), json.loads(contents["recompute_b"])
    source_bytes = subprocess.check_output(["git", "-C", str(args.rer_repo), "show", f"{args.rer_commit}:{RER_RECEIPT}"])
    source = json.loads(source_bytes)
    source_blob = subprocess.check_output(["git", "-C", str(args.rer_repo), "rev-parse", f"{args.rer_commit}:{RER_RECEIPT}"], text=True).strip()
    exact = exact_comparison(first, second)
    numerical = numerical_audit(first, second)

    modified = copy.deepcopy(second)
    changed = modified["homogeneous_polynomials"]["6"]
    epsilon = sp.Rational(1, 1000000)
    changed["coefficients"]["4,2,0"] = str(expr(changed["coefficients"]["4,2,0"])+epsilon)
    changed["polynomial"] = str(expr(changed["polynomial"])+epsilon*x**4*y**2)
    modified_result = exact_comparison(first, modified)
    mutation_rejected = not modified_result["passed"] and not modified_result["checks"]["symbol_6_all_coefficients"]

    reference_checks = {}
    for source_key, our_key in (("C4_over_a2", "c4"), ("B0_over_a4", "c6_isotropic"), ("B6_over_a4", "c6_I6")):
        reference_checks[source_key] = same(source["exact_prediction"]["coefficients"][source_key], second["normalized_coefficients"][our_key])
    for source_key, our_key in (("B0_over_C4_squared", "c2_times_c6_isotropic_over_c4_squared"),
            ("B6_over_B0", "c6_I6_over_c6_isotropic"),
            ("B6_over_C4_squared", "c2_times_c6_I6_over_c4_squared")):
        reference_checks[source_key] = same(source["exact_prediction"]["scale_free_relations"][source_key], second["scale_free_ratios"][our_key])
    prefactor_match = re.search(r"I6\(n\)\s*=\s*\(([^)]+)\)", source["exact_prediction"]["template"])
    reference_checks["I6_normalization"] = bool(prefactor_match) and same(prefactor_match.group(1), 1/expr(second["angular_normalization"]["sum_P6_at_vertex"]))

    packet_sha = digest(contents["definition_packet"])
    provenance_checks = {
        "A_declares_definition_only_scientific_input": len(first["scientific_input_files_read"]) == 1 and first["scientific_input_files_read"][0].endswith("certification/DEFINITION_ONLY.md"),
        "A_declares_no_target_values": first["target_values_not_provided"] is True,
        "B_declares_no_target_values": second["provenance"]["target_values_not_provided"] is True,
        "B_packet_hash_matches": second["provenance"]["packet_sha256"] == packet_sha,
        "B_implementation_hash_matches": second["provenance"]["implementation_sha256"] == digest(contents["implementation_b"]),
        "receipt_provenance_paths_relative": all(not Path(path).is_absolute() for path in first["files_read"]+second["provenance"]["files_read"]),
        "A_internal_checks_pass": all(first["checks"].values()) and all(first["mutations"].values()),
        "B_internal_checks_pass": all(second["checks"].values()),
    }
    checks = {"exact_AB_math_agreement": exact["passed"], "all_saved_numerical_rows_reproduced_with_bounds": numerical["passed"],
        "modified_nonheadline_sextic_coefficient_rejected": mutation_rejected,
        "pinned_RER_coefficients_ratios_normalization_agree": all(reference_checks.values()),
        "source_physical_comparison_unarmed": source["exposure_and_custody_boundary"]["comparison_permitted"] is False,
        **provenance_checks}
    result = {
        "classification": "post-seal local cross-method mathematical audit",
        "status": "passed" if all(checks.values()) else "failed",
        "exact_math_agreement": exact["passed"],
        "full_criterion_met": False,
        "external_blind_certification": False,
        "empirical_prediction_validated": False,
        "packet_sha256": packet_sha,
        "inputs": {key: {"path": relative(path), "sha256": digest(contents[key])} for key, path in files.items()},
        "pinned_RER_input": {"repository": relative(args.rer_repo), "commit": args.rer_commit,
            "path": RER_RECEIPT, "git_blob_id": source_blob, "sha256": digest(source_bytes),
            "declared_receipt_sha256": source["receipt_sha256"]},
        "method_provenance": {
            "A": "Context-isolated from OPH sources, target values, and B during initial implementation. It directly expanded twelve-port moments over Q(sqrt(5)), used exact spherical monomial averages and a generated Legendre polynomial, then evaluated 40 original-symbol rows at 80 digits.",
            "B": "Context-isolated from OPH sources, target values, and A during initial implementation. It eliminated sign pairs into a binomial moment generator, used exact spherical monomial averages and homogenized Legendre vertex normalization, then evaluated 50 original-symbol rows at 110 digits.",
            "shared_dependencies": "Both are local agent implementations using SymPy and mpmath; this is computational cross-checking, not an organizationally independent external review.",
            "comparison_phase": "Only after both implementations and receipts were saved did the comparison agent read A and the pinned RER reference. This script is explicitly informed by those artifacts.",
            "B_provenance_revision": "Absolute task-content paths were replaced with repository-relative paths; equality of every non-provenance field before and after regeneration was checked.",
            "B_nonprovenance_content_sha256": digest(json.dumps({key: value for key, value in second.items() if key != "provenance"}, sort_keys=True, separators=(",", ":")).encode()),
        },
        "exact_AB_comparison": exact,
        "direct_symbol_audit": numerical,
        "modified_coefficient_control": {
            "mutation": "Add 1/1000000 to B's x^4*y^2 coefficient in H6, updating both its expression and coefficient map, while leaving the four headline coefficients unchanged.",
            "rejected": mutation_rejected,
            "failing_checks": [key for key, passed in modified_result["checks"].items() if not passed],
            "sextic_coefficient_comparison": modified_result["coefficient_map_comparisons"]["symbol_6_all_coefficients"],
        },
        "pinned_reference_math_comparison": {"passed": all(reference_checks.values()), "checks": reference_checks,
            "coverage": "The pinned RER receipt states scalar coefficients, scale-free ratios, and I6 normalization; full Cartesian-polynomial agreement is separately established between A and B."},
        "physical_scope": {
            "packet": "Defines an equal-weight twelve-direction scalar symbol and a positive mathematical length. It supplies no physical field identification or physical lattice-scale mapping.",
            "reference_type": source["prediction_scope"]["type"],
            "reference_branch": source["prediction_scope"]["name"],
            "reference_spin_language": source["prediction_scope"]["spin_language"],
            "reference_bridge_open": source["prediction_scope"]["stronger_derivation_open"],
            "reference_branch_premises": source["branch_premises"],
            "reference_comparison_state": source["exposure_and_custody_boundary"]["comparison_state"],
            "reference_comparison_permitted": source["exposure_and_custody_boundary"]["comparison_permitted"],
            "scope_of_possible_physical_failure": source["prospective_decision_rule"]["scope_of_failure"],
            "this_audit": "Establishes agreement for the declared scalar operator only. It supplies no physical-sector bridge, dataset-specific preregistration, external referee approval, or empirical result; it does not map a to an OpenWave geometric parameter.",
        },
        "checks": checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"output": relative(args.output), "status": result["status"], "checks": checks,
                      "coefficient_entries_compared": exact["total_compared_polynomial_coefficient_entries"]}, indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
