#!/usr/bin/env python3
"""Independent decimal replay of the fixed-purity full-entropy selector.

This module does not import the producer or SymPy. It reconstructs the six
icosahedral axes and the constrained entropy branches by separate formulas.
Finite high-precision checks audit the executable algebra; the global theorem
also uses the analytic stationary-point and boundary argument in the paper.
"""

from __future__ import annotations

import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE.parent / "runs/flavor/full_entropy_w5_selector.json"
D = Decimal

# The receipt is a scientific schema, not an arbitrary expression language.
# Fix the reviewed normal forms before doing independent numerical replay:
# finitely many samples alone accept interpolating or sub-tolerance changes.
CHARACTERISTIC = "lambda**3 - 4*lambda/5 - 16*sqrt(15)/225"
BRANCH_FORMS = (
    ("sqrt(15)/6", "-sqrt(15)/30", "sqrt(15)*r/36 + 1/6", "-sqrt(15)*r/180 + 1/6", "2*sqrt(15)"),
    ("sqrt(6)/6", "-sqrt(6)/12", "sqrt(6)*r/36 + 1/6", "-sqrt(6)*r/72 + 1/6", "2*sqrt(6)"),
    ("sqrt(3)/6", "-sqrt(3)/6", "sqrt(3)*r/36 + 1/6", "-sqrt(3)*r/36 + 1/6", "2*sqrt(3)"),
    ("sqrt(6)/12", "-sqrt(6)/6", "sqrt(6)*r/72 + 1/6", "-sqrt(6)*r/36 + 1/6", "sqrt(6)"),
    ("sqrt(15)/30", "-sqrt(15)/6", "sqrt(15)*r/180 + 1/6", "-sqrt(15)*r/36 + 1/6", "2*sqrt(15)/5"),
)
BOUNDARY_FORMS = (
    "epsilon*(2 - 3*epsilon)/4", "epsilon*(1 - 2*epsilon)/3",
    "epsilon*(2 - 5*epsilon)/8", "epsilon*(1 - 3*epsilon)/5",
)
MATRIX_FORMS = {
    "0", "-2*sqrt(15)/15", "-sqrt(3)/5 + sqrt(15)/15",
    "sqrt(15)/15 + sqrt(3)/5", "2*sqrt(3)/5", "-2*sqrt(3)/5",
}
AXIS_FORMS = {"0", "1", "-1", "1/2 + sqrt(5)/2", "-sqrt(5)/2 - 1/2"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def canonical_fields(actual: dict, expected: dict, context: str) -> None:
    """Require the reviewed scientific statements, independently of replay."""
    for key, value in expected.items():
        require(actual.get(key) == value, f"{context}: noncanonical {key}")


def expression(text: str, symbols: dict[str, Decimal] | None = None) -> Decimal:
    """Read the receipt's small algebraic expressions without eval or imports."""
    require(isinstance(text, str) and len(text) <= 200, "invalid expression")
    tree = ast.parse(text, mode="eval")

    def read(node: ast.AST, depth: int = 0) -> Decimal:
        require(depth <= 20, "expression too deep")
        if isinstance(node, ast.Name) and symbols and node.id in symbols:
            return symbols[node.id]
        if isinstance(node, ast.Constant) and type(node.value) is int:
            require(abs(node.value) <= 100000, "integer too large")
            return D(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = read(node.operand, depth + 1)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp):
            left, right = read(node.left, depth + 1), read(node.right, depth + 1)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                require(right == int(right) and abs(right) <= 6, "invalid power")
                return left ** int(right)
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "sqrt" and len(node.args) == 1 and not node.keywords):
            return read(node.args[0], depth + 1).sqrt()
        raise ValueError("unsupported algebraic expression")

    return read(tree.body)


def entropy(probabilities: list[Decimal]) -> Decimal:
    require(all(p >= 0 for p in probabilities), "negative probability")
    return sum((p * (6 * p).ln() for p in probabilities if p), D(0))


def moments(probabilities: list[Decimal]) -> tuple[Decimal, Decimal]:
    return sum(probabilities, D(0)), sum((p * p for p in probabilities), D(0))


def stationary_profile(purity: Decimal, high_count: int) -> list[Decimal]:
    """Solve normalization/purity directly, independently of amplitude r."""
    m, n = D(high_count), D(6 - high_count)
    delta = 6 * purity - 1
    high = (1 + (n * delta / m).sqrt()) / 6
    low = (1 - (m * delta / n).sqrt()) / 6
    return [high] * high_count + [low] * (6 - high_count)


def geometry() -> list[list[list[Decimal]]]:
    """Compute the forward quadrupole on each of six unoriented axes."""
    phi = (1 + D(5).sqrt()) / 2
    norm = (1 + phi * phi).sqrt()
    raw = [(0, 1, phi), (0, 1, -phi), (1, phi, 0),
           (1, -phi, 0), (phi, 0, 1), (phi, 0, -1)]
    axes = [[D(x) / norm for x in row] for row in raw]
    frame = [[sum((v[i] * v[j] for v in axes), D(0))
              for j in range(3)] for i in range(3)]
    tolerance = D("1e-65")
    require(max(abs(frame[i][j] - (2 if i == j else 0))
                for i in range(3) for j in range(3)) < tolerance, "frame identity")
    for i in range(6):
        for j in range(i):
            dot = sum((a * b for a, b in zip(axes[i], axes[j])), D(0))
            require(abs(dot * dot - D(1) / 5) < tolerance, "axis overlap")
    scale = D(60).sqrt()
    matrices = []
    for selected in range(6):
        weights = [D(5 if i == selected else -1) / scale for i in range(6)]
        q = [[2 * sum((a * (v[i] * v[j] - (D(1) / 3 if i == j else 0))
                       for a, v in zip(weights, axes)), D(0))
              for j in range(3)] for i in range(3)]
        # (Q + 4/sqrt60 I)(Q - 8/sqrt60 I)=0 fixes a double
        # eigenvalue together with trace zero; no eigensolver is shared.
        for i in range(3):
            for j in range(3):
                polynomial = sum((q[i][k] * q[k][j] for k in range(3)), D(0))
                polynomial -= 4 / scale * q[i][j]
                polynomial -= D(8) / 15 if i == j else D(0)
                require(abs(polynomial) < tolerance, "quadrupole minimal polynomial")
        require(abs(sum((q[i][i] for i in range(3)), D(0))) < tolerance, "trace")
        matrices.append(q)
    return matrices


def replay() -> dict:
    """Independent controls at fixed purities, including both boundary types."""
    with localcontext() as context:
        context.prec = 80
        tolerance = D("1e-65")
        branch_count = 0
        for purity in (D(1) / 5, D(1) / 4, D(1) / 3, D(1) / 2, D(3) / 4, D(99) / 100):
            best = stationary_profile(purity, 1)
            for count in range(1, 6):
                candidate = stationary_profile(purity, count)
                if min(candidate) < 0:
                    continue
                total, square = moments(candidate)
                require(abs(total - 1) < tolerance and abs(square - purity) < tolerance,
                        "stationary constraints")
                require(entropy(candidate) >= entropy(best) - tolerance, "wrong entropy minimum")
                if count > 1 and min(candidate) > 0:
                    high, low = candidate[0], candidate[-1]
                    curvature = 1 / high - (high / low).ln() / (high - low)
                    require(curvature < 0, "repeated-high saddle")
                branch_count += 1

        epsilon = D("1e-12")
        for support in range(2, 6):
            k = D(support)
            base = [1 / k] * support + [D(0)] * (6 - support)
            delta = (epsilon / k - epsilon * epsilon * (k + 1) / (2 * k)).sqrt()
            perturbed = [(1 - epsilon) / k] * support + [epsilon] + [D(0)] * (5 - support)
            perturbed[0] += delta
            perturbed[1] -= delta
            require(max(abs(a - b) for a, b in zip(moments(base), moments(perturbed))) < tolerance,
                    "equal-support constraints")
            require(entropy(perturbed) < entropy(base), "equal-support boundary descent")

        base = [D(1) / 2, D(1) / 3, D(1) / 6, D(0), D(0), D(0)]
        total = base[0] + base[1] - epsilon
        squares = base[0] ** 2 + base[1] ** 2 - epsilon ** 2
        gap = (2 * squares - total * total).sqrt()
        perturbed = [(total + gap) / 2, (total - gap) / 2, base[2], epsilon, D(0), D(0)]
        require(max(abs(a - b) for a, b in zip(moments(base), moments(perturbed))) < tolerance,
                "unequal-support constraints")
        require(entropy(perturbed) < entropy(base), "unequal-support boundary descent")
        matrices = geometry()
        return {"stationary_branches_checked": branch_count,
                "equal_support_boundary_controls": 4,
                "unequal_support_boundary_controls": 1,
                "independent_axis_quadrupoles": len(matrices),
                "numeric_checks_are_not_a_global_proof": True}


def verify_payload(payload: dict) -> dict:
    require(payload["schema"] == "oph.full_entropy_w5_selector.v1", "schema")
    require(payload["status"] ==
            "EXACT_FULL_ENTROPY_SELECTOR_CLASSIFICATION__SINGLE_QUADRUPOLE_HIERARCHY_EXCLUDED", "status")
    require(payload["claim_class"] == "conditional_mechanism_specific_analytic_theorem", "claim class")
    require(payload["proof_grade"] ==
            "analytic_global_argument_with_executable_exact_algebra_not_formal_machine_proof", "proof grade")
    scope = payload["scope"]
    require(scope["source_law"] ==
            "declared full-relative-entropy minimization at fixed norm on the antipodal W5 record",
            "source-law scope")
    for key in ("all_OPH_completions_excluded", "core_axiom_derivation_claimed",
                "empirical_target_consumed", "legacy_flavor_template_consumed",
                "numeric_P_consumed", "old_quartic_certificate_superseded",
                "source_selected_physical_readout"):
        require(scope[key] is False, f"invalid promotion: {key}")
    require(type(scope["extra_continuous_fitted_coefficients"]) is int
            and scope["extra_continuous_fitted_coefficients"] == 0, "fitted input")
    for key in ("absolute_masses_emitted", "quark_ratios_emitted", "physical_comparison_performed",
                "positive_physical_evidence_claimed", "prediction_freeze_claimed"):
        require(payload["physical_boundary"][key] is False, f"invalid physical claim: {key}")
    repo = HERE.parents[2]
    expected_paths = {"code/particles/flavor/full_entropy_w5_selector.py",
                      "code/particles/flavor/full_entropy_w5_selector_spec.json"}
    pins = payload["source_pins"]
    require(len(pins) == 2 and {pin["path"] for pin in pins} == expected_paths, "source pins")
    for pin in pins:
        raw = (repo / pin["path"]).read_bytes()
        require(pin["bytes"] == len(raw) and pin["sha256"] == hashlib.sha256(raw).hexdigest(), "source hash")
    spec = json.loads((HERE / "full_entropy_w5_selector_spec.json").read_text())
    require(payload["physical_boundary"] == spec["physical_boundary"], "boundary/spec mismatch")
    require(spec["selector"]["objective"] == "minimize D(p||uniform) = sum_i p_i*log(6*p_i)", "wrong objective")
    require(spec["selector"]["constraints"] ==
            ["p_i >= 0", "sum_i p_i = 1", "sum_i p_i^2 = 1/6+r^2/72"], "changed feasible set")
    require(spec["source_structure"]["pair_reference"] == "1/6", "changed reference")

    def check_flags(value: object) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if key in ("exact_checks", "endpoint_checks", "controls"):
                    require(isinstance(child, dict) and child and all(v is True for v in child.values()),
                            "failed algebraic check")
                else:
                    check_flags(child)
        elif isinstance(value, list):
            for child in value:
                check_flags(child)

    check_flags(payload)
    selector = payload["selector"]
    require(selector["pair_count"] == 6 and selector["amplitude_domain"] == "0 <= r <= sqrt(60)", "domain")
    canonical_fields(selector, {
        "purity": "r**2/72 + 1/6",
        "amplitude_map": "the result holds for every admissible r and therefore any scalar r(P); no such source map is supplied or needed for this exclusion",
        "radial_rescaling_separated": "fixing a nonzero source norm changes purity under a positive radial rescaling unless the multiplier is one; the selected shape is nonetheless degenerate",
    }, "selector")
    require(selector["zero_amplitude_boundary"] ==
            "at r=0 the probability record is uniform for every normalized a; no direction is selected and the actual quadrupole r*Q is zero",
            "zero-amplitude scope")
    winner = selector["minimizer"]
    require(winner["high_multiplicity"] == 1 and winner["low_multiplicity"] == 5, "winner multiplicities")
    canonical_fields(winner, {"high": BRANCH_FORMS[0][2], "low": BRANCH_FORMS[0][3]}, "winner")
    branches = selector["stationary_branches"]
    require(len(branches) == 5, "incomplete branch enumeration")
    canonical_fields(selector["differential_certificate"], {
        "stationarity_equation": "log(p_i)+1=lambda+2*mu*p_i (constant log(6) absorbed into lambda)",
        "stationarity_strict_concavity": "-1/p**2",
        "at_most_two_positive_roots_reason": "a strictly concave logarithm minus a line has at most two zeros",
        "repeated_high_second_variation": "2*(a - b + log((b/a)**b))/(b*(a - b))",
        "log_secant_gap": "log(t) - 1 + 1/t",
        "log_secant_gap_derivative": "(t - 1)/t**2",
        "strict_sign_domain": "a>0 and t=b/a>1",
    }, "differential certificate")
    canonical_fields(selector["boundary_certificate"], {
        "unequal_support_purity_residual": "-a**2 - b**2 + epsilon**2 + (a + u)**2 + (-b + epsilon + u)**2",
        "unequal_support_implicit_derivative_at_origin": "2*a - 2*b",
        "entropy_descent": "new epsilon log epsilon term dominates the O(epsilon) change of existing positive coordinates",
        "analytic_tools": "implicit function theorem for unequal support; Taylor expansion on positive support for both cases",
    }, "boundary certificate")
    with localcontext() as context:
        context.prec = 80
        tolerance = D("1e-65")
        for radius in (D(0), D(1), D(3), D(7), D(60).sqrt()):
            symbols = {"r": radius}
            high, low = expression(winner["high"], symbols), expression(winner["low"], symbols)
            require(abs(high - (1 + 5 * radius / D(60).sqrt()) / 6) < tolerance, "winner high")
            require(abs(low - (1 - radius / D(60).sqrt()) / 6) < tolerance, "winner low")
            require(abs(expression(selector["purity"], symbols) - (D(1) / 6 + radius ** 2 / 72)) < tolerance,
                    "purity formula")
        for count, row in enumerate(branches, 1):
            require(row["high_multiplicity"] == count and row["low_multiplicity"] == 6 - count, "branch multiplicity")
            canonical_fields(row, dict(zip(
                ("high_amplitude", "low_amplitude", "high_probability", "low_probability", "strict_support_upper_r"),
                BRANCH_FORMS[count - 1])), "stationary branch")
            high, low = expression(row["high_amplitude"]), expression(row["low_amplitude"])
            require(abs(count * high + (6 - count) * low) < tolerance, "centered branch")
            require(abs(2 * (count * high ** 2 + (6 - count) * low ** 2) - 1) < tolerance, "branch norm")
            require(high > low, "branch order")
            require(abs(expression(row["strict_support_upper_r"]) + 1 / low) < tolerance, "support endpoint")
            for radius in (D(0), D(1)):
                require(abs(expression(row["high_probability"], {"r": radius}) - (1 + radius * high) / 6) < tolerance,
                        "branch high probability")
                require(abs(expression(row["low_probability"], {"r": radius}) - (1 + radius * low) / 6) < tolerance,
                        "branch low probability")
            expected = ("unique_global_minimum_up_to_permutation" if count == 1
                        else "excluded_by_negative_constrained_second_variation")
            require(row["minimum_classification"] == expected, "stationary classification")
        boundary = selector["boundary_certificate"]["equal_support_rows"]
        require([row["positive_support_size"] for row in boundary] == [2, 3, 4, 5], "boundary supports")
        for row, form in zip(boundary, BOUNDARY_FORMS):
            require(row["delta_squared"] == form, "noncanonical boundary perturbation")
            k, epsilon = D(row["positive_support_size"]), D("0.000001")
            delta2 = expression(row["delta_squared"], {"epsilon": epsilon})
            require(abs(delta2 - epsilon / k + epsilon ** 2 * (k + 1) / (2 * k)) < tolerance, "boundary perturbation")

        shape = payload["geometry"]
        require(shape["axis_count"] == 6 and shape["minimum_adjacent_gap"] == "0", "geometry scope")
        canonical_fields(shape, {
            "construction_field": "Q(sqrt(5),sqrt(15))",
            "pair_record_dimension": 5,
            "characteristic_polynomial": CHARACTERISTIC,
            "eigenvalues": ["-2*sqrt(15)/15", "-2*sqrt(15)/15", "4*sqrt(15)/15"],
            "quadrupole_convention": "Q=2 sum_i a_i (n_i n_i^T-I/3); actual centered record r*a gives r*Q",
            "forward_map_readback": "n_i^T Q n_i=(8/5)*a_i; Q is not the inverse map with unit coefficient",
            "winning_axis_scope": "0<r<=sqrt(60); at r=0 no normalized direction is selected and r*Q=0",
            "unoriented_axis_stabilizer": "D5 of order 10",
            "single_spectral_readout_boundary": "any scalar spectral function of one Q retains its double eigenvalue; independent coupled tensors are outside this candidate",
            "quadrupole_gram_entries": [["2/3" if i == j else "-2/15" for j in range(6)] for i in range(6)],
        }, "quadrupole certificate")
        wanted = [-2 / D(15).sqrt(), -2 / D(15).sqrt(), 4 / D(15).sqrt()]
        got = [expression(value) for value in shape["eigenvalues"]]
        require(len(got) == 3 and max(abs(a - b) for a, b in zip(got, wanted)) < tolerance, "eigenvalues")
        matrices = geometry()
        rows = shape["winning_axis_rows"]
        require(len(rows) == 6 and {row["axis_index"] for row in rows} == set(range(6)), "axis enumeration")
        seen = set()
        for row in rows:
            entries = row["exact_Q_entries"]
            require(len(entries) == 3 and all(len(v) == 3 for v in entries), "matrix shape")
            require(all(value in MATRIX_FORMS for line in entries for value in line), "noncanonical matrix entry")
            require(all(value in AXIS_FORMS for value in row["unnormalized_axis"]), "noncanonical axis entry")
            require(row["characteristic_polynomial"] == CHARACTERISTIC, "noncanonical axis characteristic polynomial")
            q = [[expression(value) for value in line] for line in entries]
            matches = [index for index, expected in enumerate(matrices)
                       if max(abs(q[i][j] - expected[i][j]) for i in range(3) for j in range(3)) < tolerance]
            require(len(matches) == 1 and matches[0] not in seen, "axis quadrupole mismatch or duplicate")
            seen.add(matches[0])
            axis = [expression(value) for value in row["unnormalized_axis"]]
            require(len(axis) == 3, "axis shape")
            norm2 = sum((v * v for v in axis), D(0))
            require(norm2 > 0, "zero axis")
            for i in range(3):
                for j in range(3):
                    expected = 12 / D(60).sqrt() * axis[i] * axis[j] / norm2
                    expected -= 4 / D(60).sqrt() if i == j else 0
                    require(abs(q[i][j] - expected) < tolerance, "wrong forward normalization")
            require(row["discriminant"] == "0", "split eigenvalue")
            for x in (D(0), D(1), D(2), D(3)):
                expected = x ** 3 - 4 * x / 5 - 16 * D(15).sqrt() / 225
                require(abs(expression(row["characteristic_polynomial"].replace("lambda", "x"), {"x": x}) - expected) < tolerance,
                        "characteristic polynomial")
        for x in (D(0), D(1), D(2), D(3)):
            expected = x ** 3 - 4 * x / 5 - 16 * D(15).sqrt() / 225
            require(abs(expression(shape["characteristic_polynomial"].replace("lambda", "x"), {"x": x}) - expected) < tolerance,
                    "global characteristic polynomial")
    return {"receipt": "verified", **replay()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    result = verify_payload(json.loads(args.receipt.read_text()))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
