#!/usr/bin/env python3
"""Definition-only scalar-symbol oracle, independently derived.

The mathematical input is certification/DEFINITION_ONLY.md.  This script does
not import the candidate implementation or read other recomputations.  Exact
coefficients come from a sign-paired binomial generating polynomial and exact
spherical monomial averages; the numerical path evaluates the original cosine
symbol with independently constructed mpmath directions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
PACKET = ROOT / "certification" / "DEFINITION_ONLY.md"
x, y, z, a, k = sp.symbols("x y z a k", real=True)
XYZ = (x, y, z)
R2 = x*x + y*y + z*z
PHI = (1 + sp.sqrt(5)) / 2
NORM2 = 1 + PHI**2


def exact(value):
    return sp.simplify(sp.expand(value))


def serial(value):
    return str(exact(value))


def polynomial_coefficients(expression):
    return {
        ",".join(map(str, powers)): serial(coefficient)
        for powers, coefficient in sp.Poly(sp.expand(expression), XYZ).terms()
    }


def paired_moment(order):
    """Eliminate signs first, then expand the three coordinate-plane orbits."""
    assert order % 2 == 0
    expression = 0
    for first, second in ((x, y), (y, z), (z, x)):
        expression += sum(
            sp.binomial(order, even) * first**(order-even)
            * PHI**even * second**even
            for even in range(0, order+1, 2)
        )
    return exact(4 * expression / NORM2**(order//2))


def sphere_mean(expression):
    """Exact normalized S^2 integral using Cartesian monomial moments."""
    result = 0
    for powers, coefficient in sp.Poly(sp.expand(expression), XYZ).terms():
        if any(power % 2 for power in powers):
            continue
        numerator = sp.prod(sp.factorial2(power-1) for power in powers)
        denominator = sp.factorial2(sum(powers)+1)
        result += coefficient * numerator / denominator
    return exact(result)


def make_exact_ports():
    ports = []
    norm = sp.sqrt(NORM2)
    for first_sign in (-1, 1):
        for second_sign in (-1, 1):
            ports.extend([
                sp.Matrix([0, first_sign, second_sign*PHI])/norm,
                sp.Matrix([first_sign, second_sign*PHI, 0])/norm,
                sp.Matrix([second_sign*PHI, 0, first_sign])/norm,
            ])
    return ports


def matrix_serial(matrix):
    return [[serial(value) for value in matrix.row(index)]
            for index in range(matrix.rows)]


def quadratic_tensor(ports):
    return sum((port*port.T/4 for port in ports), sp.zeros(3)).applyfunc(exact)


def tensor_isotropy(tensor):
    trace_free = (tensor-sp.trace(tensor)*sp.eye(3)/3).applyfunc(exact)
    norm_squared = exact(sp.trace(trace_free.T*trace_free))
    return trace_free == sp.zeros(3), trace_free, norm_squared


def numeric_checks(homogeneous_terms, angular_polynomial, normalization):
    mp.mp.dps = 110
    phi = (1+mp.sqrt(5))/2
    norm = mp.sqrt(1+phi*phi)
    ports = []
    # Construct original directions independently of the exact polynomial path.
    for s in (-1, 1):
        for t in (-1, 1):
            ports += [
                (mp.mpf(0), s/norm, t*phi/norm),
                (s/norm, t*phi/norm, mp.mpf(0)),
                (t*phi/norm, mp.mpf(0), s/norm),
            ]

    def unit(raw):
        magnitude = mp.sqrt(sum(value*value for value in raw))
        return tuple(mp.mpf(value)/magnitude for value in raw)

    directions = {
        "axis_x": unit((1, 0, 0)),
        "body_diagonal": unit((1, 1, 1)),
        "defined_vertex": unit((0, 1, phi)),
        "generic_1_2_3": unit((1, 2, 3)),
        "generic_2_minus3_5": unit((2, -3, 5)),
    }
    evaluate_terms = {degree: sp.lambdify(XYZ, term, "mpmath")
                      for degree, term in homogeneous_terms.items()}
    evaluate_i6 = sp.lambdify(XYZ, angular_polynomial, "mpmath")
    norm_number = mp.mpf(str(sp.N(normalization, 115)))
    rows = []
    angular_rows = []
    passed_remainder = True
    passed_identity = True
    passed_angular = True
    passed_asymptotic = True
    base_length = mp.mpf(7)/5
    k_values = [mp.mpf(1)/denominator for denominator in (5, 10, 20, 40, 80)]
    rounding_guard = mp.mpf("1e-100")

    def number(value):
        return mp.nstr(value, 45)

    for name, direction in directions.items():
        projections = [sum(u_i*n_i for u_i, n_i in zip(port, direction))
                       for port in ports]
        p6_sum = sum((231*t**6-315*t**4+105*t*t-5)/16
                     for t in projections)
        angular_from_polynomial = evaluate_i6(*direction)
        angular_error = abs(p6_sum/norm_number-angular_from_polynomial)
        angular_ok = angular_error < rounding_guard
        passed_angular &= angular_ok
        angular_rows.append({
            "direction": name,
            "unit_vector": [number(value) for value in direction],
            "I6": number(angular_from_polynomial),
            "direct_legendre_error": number(angular_error),
            "passed": bool(angular_ok),
        })
        for length_label, length in (("a", base_length), ("a_over_2", base_length/2)):
            prior_error = None
            for momentum in k_values:
                arguments = [length*momentum*t for t in projections]
                symbol = sum(1-mp.cos(value) for value in arguments)/(2*length**2)
                stable_symbol = sum(mp.sin(value/2)**2 for value in arguments)/length**2
                truncation = sum(
                    length**(degree-2)*momentum**degree*evaluate_terms[degree](*direction)
                    for degree in (2, 4, 6)
                )
                error = abs(symbol-truncation)
                # Taylor's theorem at degree seven: |cos(t)-T_6(t)| <= |t|^8/8!.
                directional_bound = sum(abs(value)**8 for value in arguments)/(2*length**2*mp.factorial(8))
                uniform_bound = len(ports)*length**6*momentum**8/(2*mp.factorial(8))
                remainder_ok = error <= directional_bound+rounding_guard
                equivalent_ok = abs(symbol-stable_symbol) < rounding_guard
                shrinking_ok = prior_error is None or error < prior_error
                passed_remainder &= remainder_ok and directional_bound <= uniform_bound+rounding_guard
                passed_identity &= equivalent_ok
                passed_asymptotic &= shrinking_ok
                rows.append({
                    "direction": name,
                    "length_case": length_label,
                    "a": number(length),
                    "k": number(momentum),
                    "direct_cosine_symbol": number(symbol),
                    "degree_six_truncation": number(truncation),
                    "absolute_error": number(error),
                    "directional_remainder_bound": number(directional_bound),
                    "uniform_remainder_bound": number(uniform_bound),
                    "error_over_directional_bound": number(error/directional_bound),
                    "error_ratio_on_halving_k": None if prior_error is None else number(prior_error/error),
                    "cosine_vs_sine_identity_error": number(abs(symbol-stable_symbol)),
                    "remainder_bound_passed": bool(remainder_ok),
                    "shrinking_error_passed": bool(shrinking_ok),
                })
                prior_error = error
    return {
        "decimal_precision": mp.mp.dps,
        "roundoff_comparison_guard": number(rounding_guard),
        "symbol_evaluation": "Original sum(1-cos(a*k*u.n))/(2*a^2); independent sin-squared identity cross-check.",
        "directional_bound_formula": "sum_u |a*k*(u.n)|^8 / (2*a^2*8!)",
        "uniform_bound_formula": "12*a^6*k^8/(2*8!) = a^6*k^8/6720",
        "bound_justification": "Taylor degree seven with bounded eighth cosine derivative; its degree-seven polynomial equals the degree-six truncation.",
        "directions": angular_rows,
        "rows": rows,
    }, {
        "direct_symbol_all_remainder_bounds": bool(passed_remainder),
        "direct_cosine_sine_identity": bool(passed_identity),
        "numeric_legendre_polynomial_identity": bool(passed_angular),
        "errors_shrink_with_decreasing_k": bool(passed_asymptotic),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "certification" / "recompute_b.json")
    args = parser.parse_args()
    packet_bytes = PACKET.read_bytes()
    script_path = Path(__file__).resolve()
    script_bytes = script_path.read_bytes()

    ports = make_exact_ports()
    moments = {degree: paired_moment(degree) for degree in (2, 4, 6)}
    # The Taylor factors follow from (1-cos(t))/(2*a^2).
    terms = {degree: exact((-1)**(degree//2+1)*moments[degree]/(2*sp.factorial(degree)))
             for degree in (2, 4, 6)}
    means = {degree: sphere_mean(term) for degree, term in terms.items()}
    isotropic_terms = {degree: means[degree]*R2**(degree//2) for degree in terms}
    anisotropic_terms = {degree: exact(terms[degree]-isotropic_terms[degree]) for degree in terms}

    # Homogenize the packet's P6 sum before normalizing at one declared vertex.
    p6_sum = exact((231*moments[6]-315*R2*moments[4]+105*R2**2*moments[2]-5*len(ports)*R2**3)/16)
    vertex = ports[0]
    normalization = exact(p6_sum.subs(dict(zip(XYZ, vertex))))
    i6_homogeneous = exact(p6_sum/normalization)
    # Determine the angular coefficient by polynomial coefficient matching.
    probe_powers, probe_value = sp.Poly(i6_homogeneous, XYZ).terms()[0]
    angular_coefficient = exact(sp.Poly(anisotropic_terms[6], XYZ).coeff_monomial(probe_powers)/probe_value)

    direct_moment_checks = {}
    for degree in moments:
        direct_moment = sum((port.dot(sp.Matrix(XYZ)))**degree for port in ports)
        direct_moment_checks[f"paired_vs_direct_moment_{degree}"] = exact(direct_moment-moments[degree]) == 0
    original_tensor = quadratic_tensor(ports)
    original_isotropic, original_traceless, original_norm = tensor_isotropy(original_tensor)
    broken_ports = ports[:-1]
    broken_tensor = quadratic_tensor(broken_ports)
    broken_isotropic, broken_traceless, broken_norm = tensor_isotropy(broken_tensor)
    renormalized_tensor = broken_tensor*sp.Rational(len(ports), len(broken_ports))
    renormalized_isotropic, _, _ = tensor_isotropy(renormalized_tensor)

    checks = {
        "twelve_ports": len(ports) == 12,
        "ports_unit_length": all(exact(port.dot(port)-1) == 0 for port in ports),
        "positive_nonzero_legendre_normalization": bool(normalization > 0),
        "quadratic_rotationally_invariant": anisotropic_terms[2] == 0,
        "quartic_rotationally_invariant": anisotropic_terms[4] == 0,
        "sextic_first_nonconstant_angular_term": anisotropic_terms[6] != 0,
        "angular_mean_zero": sphere_mean(i6_homogeneous) == 0,
        "angular_solid_harmonic": exact(sum(sp.diff(i6_homogeneous, variable, 2) for variable in XYZ)) == 0,
        "all_vertex_normalizations_one": all(exact(i6_homogeneous.subs(dict(zip(XYZ, port)))-1) == 0 for port in ports),
        "sextic_harmonic_decomposition_exact": exact(terms[6]-isotropic_terms[6]-angular_coefficient*i6_homogeneous) == 0,
        "original_second_order_tensor_isotropic": bool(original_isotropic),
        "removed_port_rejected": bool(not broken_isotropic and broken_norm > 0),
        "removed_port_rejected_after_scalar_weight_renormalization": bool(not renormalized_isotropic),
        "a_halving_quartic_factor": exact((a*a*terms[4]).subs(a, a/2)-a*a*terms[4]/4) == 0,
        "a_halving_sextic_factor": exact((a**4*terms[6]).subs(a, a/2)-a**4*terms[6]/16) == 0,
        **direct_moment_checks,
    }
    numerical, numerical_checks = numeric_checks(terms, i6_homogeneous, normalization)
    checks.update(numerical_checks)
    result = {
        "status": "passed" if all(checks.values()) else "failed",
        "scope": "Local definition-only mathematical recomputation; neither external certification nor empirical prediction.",
        "method": "Sign-paired binomial generating moments; exact spherical monomial averaging; homogenized Legendre harmonic and vertex normalization; independent 110-digit original-symbol evaluation with analytic remainder bounds.",
        "provenance": {
            "targetvalues_not_provided": True,
            "target_values_not_provided": True,
            "files_read": [str(PACKET.relative_to(ROOT.parent)), str(script_path.relative_to(ROOT.parent))],
            "files_read_scope": "Task-content files; Python and installed numerical-library implementation files are not mathematical inputs.",
            "provenance_revision": {
                "change": "Recorded task-content paths relative to the repository root instead of machine-local absolute paths.",
                "mathematical_values_changed": False,
            },
            "packet_sha256": hashlib.sha256(packet_bytes).hexdigest(),
            "implementation_sha256": hashlib.sha256(script_bytes).hexdigest(),
            "candidate_dispersion_imported_or_read": False,
            "other_recomputation_read": False,
            "web_or_oph_sources_read": False,
        },
        "conventions": {
            "momentum_vector": "p = k*n = (x,y,z), with |n|=1",
            "homogeneous_terms": "omega_squared = H2(p) + a^2*H4(p) + a^4*H6(p) + O(a^6*|p|^8)",
            "I6_homogeneous": "|p|^6*I6(p/|p|); its restriction to the unit sphere is the packet's I6.",
            "polynomial_coefficient_keys": "powers of x,y,z",
        },
        "moments": {str(degree): {"polynomial": serial(moment), "coefficients": polynomial_coefficients(moment)}
                    for degree, moment in moments.items()},
        "homogeneous_polynomials": {str(degree): {
            "polynomial": serial(term),
            "coefficients": polynomial_coefficients(term),
            "a_power": degree-2,
            "rotationally_invariant_coefficient": serial(means[degree]),
            "rotationally_invariant_polynomial": serial(isotropic_terms[degree]),
            "nonconstant_angular_polynomial": serial(anisotropic_terms[degree]),
        } for degree, term in terms.items()},
        "angular_normalization": {
            "legendre_P6": serial(sp.legendre(6, x)),
            "sum_P6_at_vertex": serial(normalization),
            "homogeneous_sum_P6": serial(p6_sum),
            "I6_homogeneous_polynomial": serial(i6_homogeneous),
            "I6_homogeneous_coefficients": polynomial_coefficients(i6_homogeneous),
            "sphere_mean": serial(sphere_mean(i6_homogeneous)),
            "reference_vertex": [serial(value) for value in vertex],
        },
        "normalized_coefficients": {
            "equation": "omega_squared/k^2 = c2 + c4*(a*k)^2 + (c6_isotropic + c6_I6*I6(n))*(a*k)^4 + O((a*k)^6)",
            "c2": serial(means[2]),
            "c4": serial(means[4]),
            "c6_isotropic": serial(means[6]),
            "c6_I6": serial(angular_coefficient),
        },
        "scale_free_ratios": {
            "c6_I6_over_c6_isotropic": serial(angular_coefficient/means[6]),
            "c2_times_c6_isotropic_over_c4_squared": serial(means[2]*means[6]/means[4]**2),
            "c2_times_c6_I6_over_c4_squared": serial(means[2]*angular_coefficient/means[4]**2),
            "note": "These are invariant under a rescaling and an overall multiplicative normalization of omega_squared. Angular normalization is fixed by the packet.",
        },
        "length_halving": {
            "quadratic_term_multiplier": "1",
            "quartic_term_multiplier": "1/4",
            "sextic_isotropic_term_multiplier": "1/16",
            "sextic_I6_term_multiplier": "1/16",
            "eighth_order_remainder_bound_multiplier": "1/64",
            "scale_free_ratios_unchanged": True,
            "comparison_held_fixed": "physical k and direction n",
        },
        "broken_control": {
            "mutation": "Remove the final unit direction exactly once, leaving eleven equally weighted ports.",
            "removed_port": [serial(value) for value in ports[-1]],
            "test": "For H2=p^T*T*p, require T-tr(T)*Identity/3 to vanish exactly.",
            "original_quadratic_tensor": matrix_serial(original_tensor),
            "original_traceless_tensor": matrix_serial(original_traceless),
            "original_traceless_frobenius_squared": serial(original_norm),
            "broken_quadratic_tensor": matrix_serial(broken_tensor),
            "broken_traceless_tensor": matrix_serial(broken_traceless),
            "broken_traceless_frobenius_squared": serial(broken_norm),
            "broken_isotropic_coefficient": serial(sp.trace(broken_tensor)/3),
            "isotropy_test_rejects": bool(not broken_isotropic),
            "scalar_total_weight_renormalization_still_rejected": bool(not renormalized_isotropic),
        },
        "direct_symbol_check": numerical,
        "checks": checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output.resolve()), "status": result["status"],
                      "checks": checks, "normalized_coefficients": result["normalized_coefficients"],
                      "normalization": serial(normalization), "ratios": result["scale_free_ratios"]}, indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
