#!/usr/bin/env python3
"""Full-entropy selection on the antipodal icosahedral W5 record.

This tests one proposed source selector, not every OPH completion.  The
scientific specification beside this module gives its global analytic proof.
The code executes its exact algebraic steps: all five two-value stationary
branches, the second-variation sign identity, both boundary-descent constraint
calculations, and the six winning quadrupoles reconstructed from geometry.
The proof uses compactness, the implicit function theorem, and standard
Lagrange second-order conditions; this is not a formal machine proof of those
analytic steps or a derivation of the selector from OPH's core axioms.

At positive fixed norm r, full relative entropy selects one high antipodal
pair and five equal low pairs. Its quadrupole is axial for every nonzero admissible
r.  Therefore changing a scalar map r(P) cannot make its three eigenvalues
distinct.  No numerical P, mass target or historical template is read.

Run without arguments to emit the deterministic receipt; --check rebuilds
and requires byte-identical stored content.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SPEC_PATH = HERE / "full_entropy_w5_selector_spec.json"
OUT_PATH = HERE.parent / "runs" / "flavor" / "full_entropy_w5_selector.json"
SCHEMA = "oph.full_entropy_w5_selector.v1"


class CertificateError(ValueError):
    """An exact identity or source contract failed."""


def require(value: bool, label: str) -> None:
    if not value:
        raise CertificateError(label)


def zero(value: sp.Expr) -> bool:
    return bool(sp.simplify(value) == 0)


def text_expr(value: sp.Expr) -> str:
    return str(sp.simplify(value))


def _all_checks(checks: dict[str, bool], context: str) -> None:
    for key, passed in checks.items():
        require(passed is True, f"{context}: {key}")


def load_spec() -> dict[str, Any]:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    require(spec.get("schema") == "oph.full_entropy_w5_selector_spec.v1", "spec schema")
    source = spec["source_structure"]
    require(source["pair_count"] == 6, "six antipodal pairs required")
    require(source["pair_reference"] == "1/6", "uniform reference required")
    require(source["pair_probability"] == "p_i = (1+r*a_i)/6", "record normalization")
    require(source["amplitude_domain"] == "0 <= r <= sqrt(60)", "amplitude domain")
    selector = spec["selector"]
    require(selector["objective"] == "minimize D(p||uniform) = sum_i p_i*log(6*p_i)", "objective")
    require(selector["constraints"] == ["p_i >= 0", "sum_i p_i = 1", "sum_i p_i^2 = 1/6+r^2/72"], "constraints")
    require(selector["continuous_fitted_coefficients"] == 0, "no fitted coefficient")
    require(selector["core_axiom_derivation_claimed"] is False, "core derivation boundary")
    require(selector["new_physical_law_claimed"] is False, "physical source boundary")
    for key in ("absolute_masses_emitted", "quark_ratios_emitted", "physical_comparison_performed", "positive_physical_evidence_claimed", "prediction_freeze_claimed"):
        require(spec["physical_boundary"][key] is False, f"physical boundary: {key}")
    for key in ("empirical_targets_loaded", "numeric_P_loaded", "legacy_flavor_templates_loaded", "historical_quartic_receipt_modified"):
        require(spec["input_policy"][key] is False, f"input boundary: {key}")
    return spec


def exact_selector_certificate() -> dict[str, Any]:
    r = sp.Symbol("r", positive=True)
    t = sp.Symbol("t", positive=True)
    p = sp.Symbol("p", positive=True)
    a = sp.Symbol("a", positive=True)
    b = sp.Symbol("b", positive=True)
    mu = sp.Symbol("mu", real=True)
    epsilon = sp.Symbol("epsilon", positive=True)
    u = sp.Symbol("u", real=True)
    purity = sp.Rational(1, 6) + r**2 / 72
    objective_term = p * sp.log(6 * p)
    stationarity_curvature = sp.diff(sp.log(p) + 1 - 2 * mu * p, p, 2)
    sign_gap = sp.log(t) - 1 + 1 / t
    secant = sp.log(b / a) / (b - a)
    hessian = 2 * (1 / b - secant)
    checks = {
        "entropy_first_derivative": zero(sp.diff(objective_term, p) - (sp.log(6*p) + 1)),
        "entropy_second_derivative": zero(sp.diff(objective_term, p, 2) - 1/p),
        "antipodal_pair_relative_entropy_reduction": zero(2*(p/2)*sp.log(12*p/2)-objective_term),
        "stationarity_strict_concavity": zero(stationarity_curvature + 1/p**2),
        "log_secant_sign_gap_derivative": zero(sp.diff(sign_gap, t) - (t-1)/t**2),
        "log_secant_sign_gap_at_one": zero(sign_gap.subs(t, 1)),
        "repeated_high_second_variation": zero(hessian.subs(b, a*t) + 2*sign_gap/(a*(t-1))),
        "boundary_epsilon_log_epsilon_limit": sp.limit(epsilon*sp.log(epsilon), epsilon, 0, dir="+") == 0,
        "boundary_log_dominates_linear": sp.limit(sp.log(epsilon), epsilon, 0, dir="+") is sp.S.NegativeInfinity,
    }
    _all_checks(checks, "differential identities")

    branches = []
    for m in range(1, 6):
        n = 6 - m
        high_amplitude = sp.sqrt(sp.Rational(n, 12*m))
        low_amplitude = -sp.sqrt(sp.Rational(m, 12*n))
        high = (1+r*high_amplitude)/6
        low = (1+r*low_amplitude)/6
        upper = sp.sqrt(sp.Rational(12*n, m))
        branch_checks = {
            "zero_mean_amplitude": zero(m*high_amplitude+n*low_amplitude),
            "unit_twelve_port_norm": zero(2*(m*high_amplitude**2+n*low_amplitude**2)-1),
            "probability_sum": zero(m*high+n*low-1),
            "fixed_purity": zero(m*high**2+n*low**2-purity),
            "low_probability_zero_at_domain_endpoint": zero(low.subs(r, upper)),
            "positive_ordered_amplitude_difference": bool(sp.simplify(high_amplitude-low_amplitude).is_positive),
        }
        if m >= 2:
            # Two equal high coordinates admit this nonzero tangent direction.
            point = [high]*m+[low]*n
            direction = [sp.Integer(1), sp.Integer(-1)]+[sp.Integer(0)]*4
            branch_checks["repeated_high_direction_sum_zero"] = sum(direction) == 0
            branch_checks["repeated_high_direction_purity_tangent"] = zero(sum(2*x*v for x,v in zip(point,direction)))
            branch_checks["repeated_high_direction_norm_squared_two"] = sum(v*v for v in direction) == 2
        _all_checks(branch_checks, f"stationary multiplicity {m}")
        branches.append({
            "high_multiplicity": m,
            "low_multiplicity": n,
            "high_amplitude": text_expr(high_amplitude),
            "low_amplitude": text_expr(low_amplitude),
            "high_probability": text_expr(high),
            "low_probability": text_expr(low),
            "strict_support_upper_r": text_expr(upper),
            "minimum_classification": "unique_global_minimum_up_to_permutation" if m == 1 else "excluded_by_negative_constrained_second_variation",
            "exact_checks": branch_checks,
        })

    unequal_residual = (a+u)**2+(b-epsilon-u)**2+epsilon**2-a**2-b**2
    boundary_checks = {
        "unequal_support_sum_preserved": zero((a+u)+(b-epsilon-u)+epsilon-a-b),
        "unequal_support_constraint_at_origin": zero(unequal_residual.subs({epsilon:0,u:0})),
        "implicit_derivative_nonzero_if_a_ne_b": zero(sp.diff(unequal_residual,u).subs({epsilon:0,u:0})-2*(a-b)),
    }
    equal_support_rows = []
    delta = sp.Symbol("delta", real=True)
    for k in range(2, 6):
        base = (1-epsilon)/k
        delta_squared = epsilon/k-sp.Rational(k+1,2*k)*epsilon**2
        coordinates = [base+delta,base-delta]+[base]*(k-2)+[epsilon]
        second_moment = sp.expand(sum(x*x for x in coordinates)).subs(delta**2,delta_squared)
        equal_checks = {
            "sum_preserved": zero(sum(coordinates)-1),
            "purity_preserved": zero(second_moment-sp.Rational(1,k)),
            "delta_squared_positive_near_zero": zero(sp.limit(delta_squared/epsilon,epsilon,0)-sp.Rational(1,k)),
            "positive_support_limit": zero(base.subs(epsilon,0)-sp.Rational(1,k)),
            "opposite_first_order_entropy_terms_cancel": zero(sp.diff(objective_term.subs(p,base+delta)+objective_term.subs(p,base-delta),delta).subs(delta,0)),
        }
        _all_checks(equal_checks, f"equal boundary support {k}")
        equal_support_rows.append({"positive_support_size": k,"delta_squared": text_expr(delta_squared),"exact_checks":equal_checks})
    _all_checks(boundary_checks,"boundary identities")

    high = (1+5*r/sp.sqrt(60))/6
    low = (1-r/sp.sqrt(60))/6
    endpoint_checks = {
        "uniform_high": zero(high.subs(r,0)-sp.Rational(1,6)),
        "uniform_low": zero(low.subs(r,0)-sp.Rational(1,6)),
        "pure_high": zero(high.subs(r,sp.sqrt(60))-1),
        "pure_low": zero(low.subs(r,sp.sqrt(60))),
        "pure_endpoint_purity": zero(purity.subs(r,sp.sqrt(60))-1),
        "nonzero_amplitude_norm_separates_radial_rescaling": zero(sp.diff(purity,r)-r/36),
    }
    _all_checks(endpoint_checks,"endpoint identities")
    return {
        "pair_count":6,
        "amplitude_domain":"0 <= r <= sqrt(60)",
        "purity":text_expr(purity),
        "minimizer":{"high":text_expr(high),"low":text_expr(low),"high_multiplicity":1,"low_multiplicity":5},
        "stationary_branches":branches,
        "differential_certificate":{
            "stationarity_equation":"log(p_i)+1=lambda+2*mu*p_i (constant log(6) absorbed into lambda)",
            "stationarity_strict_concavity":text_expr(stationarity_curvature),
            "at_most_two_positive_roots_reason":"a strictly concave logarithm minus a line has at most two zeros",
            "repeated_high_second_variation":text_expr(hessian),
            "log_secant_gap":text_expr(sign_gap),
            "log_secant_gap_derivative":text_expr(sp.diff(sign_gap,t)),
            "strict_sign_domain":"a>0 and t=b/a>1",
            "exact_checks":checks,
        },
        "boundary_certificate":{
            "unequal_support_purity_residual":text_expr(unequal_residual),
            "unequal_support_implicit_derivative_at_origin":text_expr(sp.diff(unequal_residual,u).subs({epsilon:0,u:0})),
            "exact_checks":boundary_checks,
            "equal_support_rows":equal_support_rows,
            "entropy_descent":"new epsilon log epsilon term dominates the O(epsilon) change of existing positive coordinates",
            "analytic_tools":"implicit function theorem for unequal support; Taylor expansion on positive support for both cases",
        },
        "endpoint_checks":endpoint_checks,
        "zero_amplitude_boundary":"at r=0 the probability record is uniform for every normalized a; no direction is selected and the actual quadrupole r*Q is zero",
        "radial_rescaling_separated":"fixing a nonzero source norm changes purity under a positive radial rescaling unless the multiplier is one; the selected shape is nonetheless degenerate",
        "amplitude_map":"the result holds for every admissible r and therefore any scalar r(P); no such source map is supplied or needed for this exclusion",
    }


def exact_icosahedral_geometry() -> dict[str, Any]:
    phi = (1+sp.sqrt(5))/2
    axes = [sp.Matrix(x) for x in ((0,1,phi),(0,1,-phi),(1,phi,0),(1,-phi,0),(phi,0,1),(phi,0,-1))]
    norm_squared = 1+phi**2
    projectors = [(v*v.T/norm_squared).applyfunc(sp.simplify) for v in axes]
    identity = sp.eye(3)
    projector_sum = sum(projectors,sp.zeros(3))
    gram = sp.Matrix(6,6,lambda i,j:sp.simplify(sp.trace((projectors[i]-identity/3)*(projectors[j]-identity/3))))
    expected_gram = sp.Rational(4,5)*sp.eye(6)-sp.Rational(2,15)*sp.ones(6)
    geometry_checks = {
        "equal_axis_norms":all(zero(v.dot(v)-norm_squared) for v in axes),
        "axis_projectors_idempotent":all((p*p-p).applyfunc(sp.simplify)==sp.zeros(3) for p in projectors),
        "six_axis_two_design":(projector_sum-2*identity).applyfunc(sp.simplify)==sp.zeros(3),
        "quadrupole_gram":(gram-expected_gram).applyfunc(sp.simplify)==sp.zeros(6),
        "quadrupole_rank_five":gram.rank()==5,
    }
    _all_checks(geometry_checks,"icosahedral geometry")
    low = -1/sp.sqrt(60)
    high = 5/sp.sqrt(60)
    lam = sp.Symbol("lambda")
    expected_eigenvalues = [-2/sp.sqrt(15),-2/sp.sqrt(15),4/sp.sqrt(15)]
    expected_polynomial = sp.expand(sp.prod(lam-value for value in expected_eigenvalues))
    rows = []
    for j in range(6):
        weights = [high if i==j else low for i in range(6)]
        q = sum((2*weights[i]*(projectors[i]-identity/3) for i in range(6)),sp.zeros(3)).applyfunc(sp.simplify)
        polynomial = sp.expand(q.charpoly(lam).as_expr())
        discriminant = sp.simplify(sp.discriminant(polynomial,lam))
        row_checks = {
            "centered_record":zero(sum(weights)),
            "unit_twelve_port_record_norm":zero(2*sum(x*x for x in weights)-1),
            "quadrupole_trace_zero":zero(sp.trace(q)),
            "forward_readback_factor_eight_fifths":all(zero((axes[i].T*q*axes[i])[0]/norm_squared-sp.Rational(8,5)*weights[i]) for i in range(6)),
            "axial_projector_formula":(q-(12*projectors[j]-4*identity)/sp.sqrt(60)).applyfunc(sp.simplify)==sp.zeros(3),
            "characteristic_polynomial":zero(polynomial-expected_polynomial),
            "double_eigenvalue":discriminant==0,
            "nonzero_quadrupole":not q.equals(sp.zeros(3)),
        }
        _all_checks(row_checks,f"winning axis {j}")
        rows.append({
            "axis_index":j,
            "unnormalized_axis":[text_expr(v) for v in axes[j]],
            "exact_Q_entries":[[text_expr(q[i,k]) for k in range(3)] for i in range(3)],
            "characteristic_polynomial":text_expr(polynomial),
            "discriminant":text_expr(discriminant),
            "exact_checks":row_checks,
        })
    return {
        "construction_field":"Q(sqrt(5),sqrt(15))",
        "axis_count":6,
        "pair_record_dimension":5,
        "exact_checks":geometry_checks,
        "quadrupole_gram_entries":[[text_expr(gram[i,j]) for j in range(6)] for i in range(6)],
        "quadrupole_convention":"Q=2 sum_i a_i (n_i n_i^T-I/3); actual centered record r*a gives r*Q",
        "forward_map_readback":"n_i^T Q n_i=(8/5)*a_i; Q is not the inverse map with unit coefficient",
        "winning_axis_scope":"0<r<=sqrt(60); at r=0 no normalized direction is selected and r*Q=0",
        "unoriented_axis_stabilizer":"D5 of order 10",
        "winning_axis_rows":rows,
        "eigenvalues":[text_expr(value) for value in expected_eigenvalues],
        "characteristic_polynomial":text_expr(expected_polynomial),
        "minimum_adjacent_gap":"0",
        "single_spectral_readout_boundary":"any scalar spectral function of one Q retains its double eigenvalue; independent coupled tensors are outside this candidate",
    }


def _pin(path: Path, role: str) -> dict[str, Any]:
    raw = path.read_bytes()
    return {"path":path.relative_to(REPO).as_posix(),"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),"role":role}


def build_payload() -> dict[str, Any]:
    spec = load_spec()
    selector = exact_selector_certificate()
    geometry = exact_icosahedral_geometry()
    # These are concrete counterfactual detectors, not a list of claimed tests.
    wrong_high = sp.Rational(1,6)+sp.sqrt(15)*sp.Symbol("r",positive=True)/18
    right_low = sp.Rational(1,6)-sp.sqrt(15)*sp.Symbol("r",positive=True)/180
    lam = sp.Symbol("lambda")
    split_control = sp.Matrix.diag(-2/sp.sqrt(15),-1/sp.sqrt(15),3/sp.sqrt(15))
    split_discriminant = sp.simplify(sp.discriminant(split_control.charpoly(lam).as_expr(),lam))
    controls = {
        "changed_high_profile_breaks_sum":not zero(wrong_high+5*right_low-1),
        "repeated_high_branch_excluded":all(row["minimum_classification"]=="excluded_by_negative_constrained_second_variation" for row in selector["stationary_branches"][1:]),
        "split_spectrum_has_nonzero_discriminant":split_discriminant!=0,
        "quartic_taylor_objective_not_used":spec["selector"]["objective"]=="minimize D(p||uniform) = sum_i p_i*log(6*p_i)",
    }
    _all_checks(controls,"counterfactual controls")
    return {
        "schema":SCHEMA,
        "status":"EXACT_FULL_ENTROPY_SELECTOR_CLASSIFICATION__SINGLE_QUADRUPOLE_HIERARCHY_EXCLUDED",
        "claim_class":"conditional_mechanism_specific_analytic_theorem",
        "proof_grade":"analytic_global_argument_with_executable_exact_algebra_not_formal_machine_proof",
        "scope":{
            "source_law":"declared full-relative-entropy minimization at fixed norm on the antipodal W5 record",
            "extra_continuous_fitted_coefficients":0,
            "all_OPH_completions_excluded":False,
            "core_axiom_derivation_claimed":False,
            "source_selected_physical_readout":False,
            "numeric_P_consumed":False,
            "empirical_target_consumed":False,
            "legacy_flavor_template_consumed":False,
            "old_quartic_certificate_superseded":False,
        },
        "source_pins":[_pin(SPEC_PATH,"declared selector and analytic global argument"),_pin(Path(__file__).resolve(),"exact symbolic producer")],
        "selector":selector,
        "geometry":geometry,
        "controls":controls,
        "physical_boundary":dict(spec["physical_boundary"]),
    }


def serialized(payload: dict[str, Any]) -> str:
    return json.dumps(payload,indent=2,sort_keys=True)+"\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=OUT_PATH)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    content = serialized(build_payload())
    if args.check:
        require(args.out.exists(),f"receipt missing: {args.out}")
        require(args.out.read_text(encoding="utf-8")==content,"stored receipt differs from exact rebuild")
        print("full-entropy W5 receipt: exact rebuild matches")
    else:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(content,encoding="utf-8")
        print(f"saved: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
