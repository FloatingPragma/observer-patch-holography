"""Definition-only local mathematical recomputation of a scalar cosine symbol.

Scientific input: certification/DEFINITION_ONLY.md only. This is neither
external certification nor an empirical prediction or physical identification.
Run this file to print the complete exact and numerical result as JSON.
"""
from __future__ import annotations

import json
from math import factorial

import mpmath as mp
import sympy as s

X, Y, Z, A = s.symbols("x y z a", real=True)
XYZ = (X, Y, Z)
R2 = X**2 + Y**2 + Z**2
PHI = (1 + s.sqrt(5)) / 2
RAW = [(0, i, j * PHI) for i in (-1, 1) for j in (-1, 1)]
RAW += [(i, j * PHI, 0) for i in (-1, 1) for j in (-1, 1)]
RAW += [(j * PHI, 0, i) for i in (-1, 1) for j in (-1, 1)]
NORM2 = s.simplify(1 + PHI**2)


def clean(expr):
    """Reduce each polynomial coefficient in Q(sqrt(5)) exactly."""
    return s.Poly.from_dict(
        {powers: s.simplify(c) for powers, c in s.Poly(s.expand(expr), XYZ).terms()},
        XYZ,
    ).as_expr()


def zero(expr):
    return clean(expr) == 0


def moment(degree, ports=RAW, weights=None):
    weights = [s.Integer(1)] * len(ports) if weights is None else weights
    return clean(sum(w * sum(v * n for v, n in zip(u, XYZ))**degree
                     for u, w in zip(ports, weights)) / NORM2**(degree // 2))


def sphere_average(poly):
    """Normalized S^2 average using exact even-coordinate moments."""
    result = 0
    for powers, coefficient in s.Poly(poly, XYZ).terms():
        if any(power % 2 for power in powers):
            continue
        numerator = s.prod(s.factorial2(power - 1) for power in powers)
        result += coefficient * numerator / s.factorial2(sum(powers) + 1)
    return s.simplify(result)


def isotropic(poly, degree):
    return zero(poly - sphere_average(poly) * R2**(degree // 2))


def polynomial_record(poly):
    return {"expression": str(poly), "monomials": {
        ",".join(map(str, powers)): str(c)
        for powers, c in s.Poly(poly, XYZ).terms()}}


def exact_calculation():
    moments = {degree: moment(degree) for degree in (2, 4, 6)}
    averages = {degree: sphere_average(poly) for degree, poly in moments.items()}
    harmonic = clean(moments[6] - averages[6] * R2**3)
    # Homogenize the generated Legendre polynomial before summing over ports.
    t = s.Symbol("t")
    legendre = s.Poly(s.legendre(6, t), t)
    legendre_sum = clean(sum(coefficient * (moments[power[0]] if power[0] else 12)
                            * R2**((6 - power[0]) // 2)
                            for power, coefficient in legendre.terms()))
    normalizations = [s.simplify(legendre_sum.subs(dict(zip(XYZ, u))) / NORM2**3)
                      for u in RAW]
    normalization = normalizations[0]
    i6 = clean(legendre_sum / normalization)
    angular_multiplier = s.simplify(s.Poly(harmonic, XYZ).coeff_monomial(X**6)
                                    / s.Poly(i6, XYZ).coeff_monomial(X**6))
    prefactors = {degree: s.Rational((-1)**(degree // 2 + 1), 2 * factorial(degree))
                  for degree in moments}
    omega = {degree: clean(prefactors[degree] * A**(degree - 2) * poly)
             for degree, poly in moments.items()}
    c2 = prefactors[2] * averages[2]
    c4 = prefactors[4] * averages[4]
    c6 = prefactors[6] * averages[6]
    c6a = s.simplify(prefactors[6] * angular_multiplier)
    checks = {
        "twelve_unique_ports": len(RAW) == len(set(RAW)) == 12,
        "all_ports_unit": all(s.simplify(sum(v*v for v in u) / NORM2) == 1 for u in RAW),
        "quadratic_isotropic": isotropic(moments[2], 2),
        "quartic_isotropic": isotropic(moments[4], 4),
        "sextic_nonconstant": not isotropic(moments[6], 6),
        "sextic_remainder_harmonic": zero(sum(s.diff(harmonic, v, 2) for v in XYZ)),
        "legendre_matches_harmonic": zero(legendre_sum - s.Rational(231, 16)*harmonic),
        "normalization_nonzero": normalization != 0,
        "all_vertex_normalizations_equal": all(v == normalization for v in normalizations),
        "i6_all_vertices_one": all(s.simplify(i6.subs(dict(zip(XYZ, u))) / NORM2**3) == 1 for u in RAW),
        "i6_spherical_mean_zero": sphere_average(i6) == 0,
        "exact_angular_reconstruction": zero(harmonic - angular_multiplier*i6),
        "length_halving_identity": all(zero(poly.subs(A, A/2) - poly/2**(degree-2))
                                       for degree, poly in omega.items()),
    }
    broken = {"one_port_removed": moment(2, RAW[:-1]),
              "one_weight_doubled": moment(2, weights=[2] + [1]*11)}
    mutations = {name + "_rejected_by_quadratic_isotropy": not isotropic(poly, 2)
                 for name, poly in broken.items()}
    outputs = {
        "legendre_P6": str(legendre.as_expr()),
        "moment_polynomials": {str(d): polynomial_record(p) for d, p in moments.items()},
        "omega_homogeneous_polynomials": {str(d): polynomial_record(p) for d, p in omega.items()},
        "moment_sphere_averages": {str(d): str(v) for d, v in averages.items()},
        "omega_rotational_parts": {str(d): str(s.expand(prefactors[d]*A**(d-2)*averages[d]*R2**(d//2)))
                                   for d in moments},
        "first_nonconstant_angular_degree": 6,
        "harmonic_sextic_moment": polynomial_record(harmonic),
        "legendre_sum_homogeneous": polynomial_record(legendre_sum),
        "I6_vertex_denominator": str(normalization),
        "I6_homogeneous": polynomial_record(i6),
        "I6_usage": "Evaluate the homogeneous polynomial on n with n dot n = 1.",
        "expansion": "omega_squared = c2*k^2 + c4*a^2*k^4 + (c6_isotropic + c6_angular*I6(n))*a^4*k^6 + O(a^6*k^8)",
        "coefficients": dict(zip(("c2", "c4", "c6_isotropic", "c6_angular"),
                                 map(str, (c2, c4, c6, c6a)))),
        "scale_free_ratios": {"c6_angular/c6_isotropic": str(s.cancel(c6a/c6)),
                              "c2*c6_isotropic/c4^2": str(s.cancel(c2*c6/c4**2)),
                              "c2*c6_angular/c4^2": str(s.cancel(c2*c6a/c4**2))},
        "a_to_a_over_2_at_fixed_k_and_n": {"quadratic_multiplier": "1", "quartic_multiplier": "1/4",
                                           "sextic_multiplier": "1/16", "remainder_bound_multiplier": "1/64",
                                           "scale_free_ratios": "unchanged"},
        "broken_control_quadratics": {name: polynomial_record(p/4) for name, p in broken.items()},
    }
    return outputs, checks, mutations, i6, (c2, c4, c6, c6a)


def direct_check(i6, coefficients):
    """Independent 80-digit cosine evaluations with rigorous analytic bounds."""
    rows = []
    with mp.workdps(80):
        phi = (1 + mp.sqrt(5)) / 2
        radius = mp.sqrt(1 + phi**2)
        raw = [(0, i, j*phi) for i in (-1, 1) for j in (-1, 1)]
        raw += [(i, j*phi, 0) for i in (-1, 1) for j in (-1, 1)]
        raw += [(j*phi, 0, i) for i in (-1, 1) for j in (-1, 1)]
        directions = {"axis": (1, 0, 0), "diagonal": (1, 1, 1),
                      "generic": (1, 2, 3), "port": (0, 1, phi)}
        i6fun = s.lambdify(XYZ, i6, "mpmath")
        c2, c4, c6, c6a = [mp.mpf(str(s.N(v, 85))) for v in coefficients]
        valid = stable = decay = True
        for label, direction in directions.items():
            length = mp.sqrt(sum(v*v for v in direction))
            n = [mp.mpf(v)/length for v in direction]
            dots = [sum(v*q for v, q in zip(u, n))/radius for u in raw]
            for a in (mp.mpf(1), mp.mpf("0.5")):
                previous = None
                for k in map(mp.mpf, ("0.4", "0.2", "0.1", "0.05", "0.025")):
                    direct = sum(1-mp.cos(a*k*d) for d in dots)/(2*a*a)
                    sine = sum(mp.sin(a*k*d/2)**2 for d in dots)/(a*a)
                    truncated = c2*k*k + c4*a*a*k**4 + (c6+c6a*i6fun(*n))*a**4*k**6
                    error = direct-truncated
                    b8 = a**6*k**8*sum(d**8 for d in dots)/(2*factorial(8))
                    b10 = a**8*k**10*sum(d**10 for d in dots)/(2*factorial(10))
                    tolerance = mp.mpf("1e-70")
                    # Taylor theorem gives |R| <= b8 at any k. The alternating
                    # refinement used here holds because all |a*k*d| <= 0.4.
                    valid &= bool(max(abs(a*k*d) for d in dots) < 1
                                  and abs(error) <= b8+tolerance
                                  and -b8-tolerance <= error <= -b8+b10+tolerance)
                    stable &= bool(abs(direct-sine) < tolerance)
                    if previous is not None:
                        decay &= bool(abs(error) < previous/200)
                    previous = abs(error)
                    rows.append({"direction": label, "a": str(a), "k": str(k),
                                 "symbol": mp.nstr(direct, 65), "truncation": mp.nstr(truncated, 65),
                                 "signed_error": mp.nstr(error, 35), "absolute_bound_B8": mp.nstr(b8, 35),
                                 "next_bound_B10": mp.nstr(b10, 35), "abs_error_over_B8": mp.nstr(abs(error)/b8, 35)})
    return {"decimal_precision": 80,
            "analytic_absolute_bound": "B8 = a^6*k^8*sum_u(u dot n)^8/(2*8!)",
            "analytic_interval": "-B8 <= omega_squared - truncation <= -B8+B10 when max|a*k*u dot n| <= 1",
            "B10": "a^8*k^10*sum_u(u dot n)^10/(2*10!)", "rows": rows}, {
                "direct_symbol_remainder_bounds": valid,
                "direct_cosine_matches_stable_sine": stable,
                "remainder_decays_at_least_200fold_per_k_halving": decay}


def run():
    outputs, checks, mutations, i6, coefficients = exact_calculation()
    outputs["direct_symbol"], numerical = direct_check(i6, coefficients)
    checks.update(numerical)
    payload = {
        "classification": "definition-only local mathematical recompute",
        "target_values_not_provided": True,
        "scientific_input_files_read": ["certification/DEFINITION_ONLY.md"],
        "files_read": [
            "certification/DEFINITION_ONLY.md", "candidate/oph_openwave/dispersion.py",
            "candidate/tests/test_dispersion.py", "candidate/methods/dispersion.md",
            "certification/recompute_a.json"],
        "file_audit_scope": "Original implementation-phase read declaration in repository-relative paths, not runtime filesystem reads. Scientific packet and own artifacts were read during implementation and final review; standard installed Python libraries are excluded.",
        "checks": checks, "mutations": mutations, "outputs": outputs,
        "full_criterion_met": False,
        "scope": "No physical identification of the scalar variable, field, or length a is supplied or inferred. No empirical prediction, external certification, experimental comparison, or physical lattice-scale claim is made.",
    }
    if not all(checks.values()) or not all(mutations.values()):
        raise AssertionError("Internal mathematical verification failed")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
