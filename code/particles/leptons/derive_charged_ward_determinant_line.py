#!/usr/bin/env python3
"""Ward-normalized charged determinant-line inversion lane.

The exact Euclidean one-loop lepton vacuum-polarization kernel makes the
Ward-normalized inverse-coupling increment a strictly decreasing bijection of
the common charged log-mass coordinate mu.  A positive leptonic remainder
R_ell therefore selects exactly one root mu_ch, and that root fixes the
determinant line det(M_e/E_star) = exp(3 mu_ch) together with the per-flavor
ratios m_i/E_star = exp(mu_ch + ell_i).

This lane is the X2.2 electron-absolute-ratio route B producer.  Every
numeric packet in this module is synthetic and is labeled synthetic.  The
physical source packet is absent, the promotion gate is fail-closed, and no
charged reference mass is read.
"""

from __future__ import annotations

import functools

import argparse
import json
from pathlib import Path
from typing import Any, Callable, Sequence

import mpmath as mp


WORKING_DPS = 80


def _scoped_dps(func):
    """Run func at WORKING_DPS and restore the global precision on exit."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with mp.workdps(WORKING_DPS):
            return func(*args, **kwargs)

    return wrapper


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons" / "charged_ward_determinant_line.json"
)

KERNEL_CHECK_ARGUMENTS = ("1e-8", "0.1", "1", "10", "1e6")
KERNEL_IDENTITY_TOLERANCE = "1e-40"

SYNTHETIC_Q_OVER_ESTAR = "100"
SYNTHETIC_CENTERED_SHAPE = ("-2", "0.25", "1.75")
SYNTHETIC_MU_TRUE = "-3.25"

OPEN_PARENTS = (
    (
        "R_CHARGED_CENTERED_SHAPE",
        "the current centered charged shape is target-anchored",
    ),
    (
        "R_WARD_SPACELIKE_ENDPOINT_PAIR",
        "the current source does not emit a complete target-free low/high Ward pair",
    ),
    (
        "R_WARD_NONLEPTONIC_SUBTRACTION",
        "source QCD/hadronic and source quark packets are incomplete",
    ),
    (
        "R_WARD_HIGHER_ORDER_MONOTONICITY",
        "a full lepton/mixed remainder and derivative certificate is absent",
    ),
    (
        "R_ELECTRON_ATOMIC_MASS_TRANSPORT",
        "a source-only match from the Ward threshold coordinate to the cesium "
        "electron parameter is absent",
    ),
    (
        "R_CLEAN_P_ALPHA_ROOT",
        "strict no-target interval root replay is absent from the supplied "
        "particle archive",
    ),
    (
        "R_NO_TARGET_CHARGED_DAG",
        "no complete transitive source DAG and prospective freeze is supplied",
    ),
)


def _mpf(value: Any) -> mp.mpf:
    if isinstance(value, bool) or value is None:
        raise ValueError(f"expected a finite real number, received {value!r}")
    result = value if isinstance(value, mp.mpf) else mp.mpf(str(value))
    if not mp.isfinite(result):
        raise ValueError(f"expected a finite real number, received {value!r}")
    return result


def _text(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(mp.mpf(value), digits)


def _centered(values: Sequence[Any]) -> list[mp.mpf]:
    ell = [_mpf(value) for value in values]
    if len(ell) != 3:
        raise ValueError("exactly three centered charged logs are required")
    if abs(sum(ell)) > mp.mpf("1e-40"):
        raise ValueError("centered charged logs must sum to zero")
    return ell


def _quad_knots(z: mp.mpf) -> list[mp.mpf]:
    """Feynman-parameter knots that resolve the 1/z boundary layer."""

    if z > 100:
        return [mp.mpf(0), 1 / z, mp.mpf("0.5"), 1 - 1 / z, mp.mpf(1)]
    return [mp.mpf(0), mp.mpf(1)]


def ward_kernel_integral(z: Any) -> mp.mpf:
    """Direct quadrature of I(z) = integral_0^1 x(1-x) log(1+z x(1-x)) dx."""

    z = _mpf(z)
    if z <= 0:
        raise ValueError("the Ward kernel argument z must be positive")
    return mp.quad(lambda x: x * (1 - x) * mp.log(1 + z * x * (1 - x)), _quad_knots(z))


def ward_kernel_closed_form(z: Any) -> mp.mpf:
    """Exact closed form of the kernel, algebraically equal to the integral.

    I(z) = (1/6)[-5/3 + 4/z + (1-2/z) sqrt(1+4/z)
                 log((sqrt(1+4/z)+1)/(sqrt(1+4/z)-1))].
    """

    z = _mpf(z)
    if z <= 0:
        raise ValueError("the Ward kernel argument z must be positive")
    root = mp.sqrt(1 + 4 / z)
    logarithm = mp.log((root + 1) / (root - 1))
    return (mp.mpf(-5) / 3 + 4 / z + (1 - 2 / z) * root * logarithm) / 6


def ward_kernel_series(z: Any, max_terms: int = 4000) -> mp.mpf:
    """Absolutely convergent kernel power series on 0 < z < 4.

    I(z) = sum_{n>=1} (-1)^(n+1) (z^n/n) [(n+1)!]^2 / (2n+3)!.
    """

    z = _mpf(z)
    if not 0 < z < 4:
        raise ValueError("the kernel series applies on 0 < z < 4")
    total = mp.mpf(0)
    power = z
    for n in range(1, max_terms):
        coefficient = mp.factorial(n + 1) ** 2 / (n * mp.factorial(2 * n + 3))
        term = coefficient * power
        total = total + term if n % 2 else total - term
        if abs(term) < mp.eps * max(mp.mpf(1), abs(total)):
            return total
        power *= z
    raise ArithmeticError("the kernel series did not converge in the allotted terms")


def ward_kernel_I(z: Any) -> mp.mpf:
    """Kernel evaluator: series for small z, closed form elsewhere."""

    z = _mpf(z)
    if z < mp.mpf("0.25"):
        return ward_kernel_series(z)
    return ward_kernel_closed_form(z)


def ward_response(mu: Any, q_over_estar: Any, centered_logs: Sequence[Any]) -> mp.mpf:
    """One-loop increment W_Q(mu; ell) = (2/pi) sum_i I(q^2 exp(-2(mu+ell_i)))."""

    mu = _mpf(mu)
    q = _mpf(q_over_estar)
    if q <= 0:
        raise ValueError("q_over_Estar must be positive")
    ell = _centered(centered_logs)
    total = mp.mpf(0)
    for value in ell:
        total += ward_kernel_I(q**2 * mp.e ** (-2 * (mu + value)))
    return 2 * total / mp.pi


def ward_response_derivative(
    mu: Any, q_over_estar: Any, centered_logs: Sequence[Any]
) -> mp.mpf:
    """Exact derivative of the response in the common coordinate mu.

    dW_Q/dmu = -(4/pi) sum_i integral_0^1 z_i [x(1-x)]^2 / (1+z_i x(1-x)) dx,
    strictly negative for finite mu and q > 0.
    """

    mu = _mpf(mu)
    q = _mpf(q_over_estar)
    if q <= 0:
        raise ValueError("q_over_Estar must be positive")
    ell = _centered(centered_logs)
    total = mp.mpf(0)
    for value in ell:
        z = q**2 * mp.e ** (-2 * (mu + value))
        total += mp.quad(
            lambda x, z=z: z * (x * (1 - x)) ** 2 / (1 + z * x * (1 - x)),
            _quad_knots(z),
        )
    return -4 * total / mp.pi


def inversion_seed(response: Any, q_over_estar: Any) -> mp.mpf:
    """High-energy seed mu_seed = log q - 5/6 - (pi/2) R_ell."""

    return mp.log(_mpf(q_over_estar)) - mp.mpf(5) / 6 - mp.pi * _mpf(response) / 2


def solve_mu(
    response: Any,
    q_over_estar: Any,
    centered_logs: Sequence[Any],
    *,
    tolerance: str = "1e-50",
) -> mp.mpf:
    """Solve W_Q(mu; ell) = response by seeded, globally safe bisection.

    Strict monotonicity makes the root unique; the seed places the initial
    bracket within the asymptotic error of Corollary 3.2.
    """

    target = _mpf(response)
    if target <= 0:
        raise ValueError("the leptonic inverse-coupling increment must be positive")
    q = _mpf(q_over_estar)
    ell = _centered(centered_logs)
    seed = inversion_seed(target, q)
    left, right = seed - 1, seed + 1
    while ward_response(left, q, ell) <= target:
        left -= 2 * (right - left)
    while ward_response(right, q, ell) >= target:
        right += 2 * (right - left)
    tol = _mpf(tolerance)
    while right - left > tol:
        midpoint = (left + right) / 2
        if ward_response(midpoint, q, ell) > target:
            left = midpoint
        else:
            right = midpoint
    return (left + right) / 2


def determinant_line(mu_ch: Any, centered_logs: Sequence[Any]) -> dict[str, Any]:
    """Emit the determinant line and per-flavor ratios fixed by mu_ch."""

    mu = _mpf(mu_ch)
    ell = _centered(centered_logs)
    return {
        "mu_ch": _text(mu),
        "det_line_D_M": _text(mp.e ** (3 * mu)),
        "per_flavor_mass_over_Estar": [_text(mp.e ** (mu + value)) for value in ell],
        "formulae": {
            "det_line": "D_M = det(M_e/E_star) = exp(3 mu_ch)",
            "per_flavor": "m_i/E_star = exp(mu_ch + ell_i)",
        },
    }


def _envelope_response(mu: mp.mpf, q: mp.mpf, shifts: Sequence[mp.mpf]) -> mp.mpf:
    """Response evaluated at fixed endpoint data without the sum-zero gate."""

    total = mp.mpf(0)
    for value in shifts:
        total += ward_kernel_I(q**2 * mp.e ** (-2 * (mu + value)))
    return 2 * total / mp.pi


def _edge_bisect(
    predicate: Callable[[mp.mpf], bool],
    true_point: mp.mpf,
    false_point: mp.mpf,
    iterations: int = 80,
) -> mp.mpf:
    """Return the tightest certified point on the true side of a boundary."""

    for _ in range(iterations):
        midpoint = (true_point + false_point) / 2
        if predicate(midpoint):
            true_point = midpoint
        else:
            false_point = midpoint
    return true_point


def certify_mu_interval(
    response_interval: Sequence[Any],
    q_interval: Sequence[Any],
    centered_log_intervals: Sequence[Sequence[Any]],
) -> tuple[mp.mpf, mp.mpf, dict[str, Any]]:
    """Certified enclosure of every admissible root (interval theorem).

    Conservative envelopes W_-(mu) = (2/pi) sum_i I(q_-^2 exp(-2(mu+ell_i^+)))
    and W_+(mu) = (2/pi) sum_i I(q_+^2 exp(-2(mu+ell_i^-))) bracket the exact
    response.  Points with W_-(mu_-) > R_+ and W_+(mu_+) < R_- bound every
    root compatible with the declared intervals.
    """

    r_lo, r_hi = (_mpf(value) for value in response_interval)
    if r_lo <= 0 or r_lo > r_hi:
        raise ValueError("the required response interval must be positive and ordered")
    q_lo, q_hi = (_mpf(value) for value in q_interval)
    if q_lo <= 0 or q_lo > q_hi:
        raise ValueError("the q interval must be positive and ordered")
    if len(centered_log_intervals) != 3:
        raise ValueError("exactly three centered-log intervals are required")
    ell_lo = [_mpf(pair[0]) for pair in centered_log_intervals]
    ell_hi = [_mpf(pair[1]) for pair in centered_log_intervals]
    if any(lo > hi for lo, hi in zip(ell_lo, ell_hi)):
        raise ValueError("a centered-log interval is reversed")

    def below_all_roots(mu: mp.mpf) -> bool:
        return _envelope_response(mu, q_lo, ell_hi) > r_hi

    def above_all_roots(mu: mp.mpf) -> bool:
        return _envelope_response(mu, q_hi, ell_lo) < r_lo

    midpoint_shape = [(lo + hi) / 2 for lo, hi in zip(ell_lo, ell_hi)]
    mean = sum(midpoint_shape) / 3
    midpoint_shape = [value - mean for value in midpoint_shape]
    center = solve_mu(
        (r_lo + r_hi) / 2, (q_lo + q_hi) / 2, midpoint_shape, tolerance="1e-40"
    )

    left, right = center - 1, center + 1
    while not below_all_roots(left):
        left -= 2 * (right - left)
    while not above_all_roots(right):
        right += 2 * (right - left)

    mu_lower = _edge_bisect(below_all_roots, left, right)
    mu_upper = _edge_bisect(above_all_roots, right, left)
    certified = below_all_roots(mu_lower) and above_all_roots(mu_upper)
    if not certified or mu_lower >= mu_upper:
        raise ArithmeticError("the monotone root enclosure failed to certify")
    evidence = {
        "W_minus_at_mu_lower": _text(_envelope_response(mu_lower, q_lo, ell_hi)),
        "W_plus_at_mu_upper": _text(_envelope_response(mu_upper, q_hi, ell_lo)),
        "W_minus_exceeds_R_upper": True,
        "W_plus_below_R_lower": True,
    }
    return mu_lower, mu_upper, evidence


@_scoped_dps
def kernel_identity_checks() -> dict[str, Any]:
    """Closed form against direct quadrature at the declared arguments."""

    tolerance = mp.mpf(KERNEL_IDENTITY_TOLERANCE)
    rows = []
    all_within = True
    for z_text in KERNEL_CHECK_ARGUMENTS:
        direct = ward_kernel_integral(z_text)
        closed = ward_kernel_closed_form(z_text)
        residual = abs(direct - closed)
        within = residual < tolerance
        all_within = all_within and within
        rows.append(
            {
                "z": z_text,
                "closed_form": _text(closed),
                "direct_quadrature": _text(direct),
                "absolute_residual": _text(residual),
                "within_tolerance": bool(within),
            }
        )
    series_residual = abs(ward_kernel_series("0.1") - ward_kernel_closed_form("0.1"))
    return {
        "integral_form": "I(z) = integral_0^1 x(1-x) log(1+z x(1-x)) dx",
        "closed_form": (
            "I(z) = (1/6)[-5/3 + 4/z + (1-2/z) sqrt(1+4/z) "
            "log((sqrt(1+4/z)+1)/(sqrt(1+4/z)-1))]"
        ),
        "tolerance": KERNEL_IDENTITY_TOLERANCE,
        "rows": rows,
        "series_vs_closed_form_residual_at_z_0p1": _text(series_residual),
        "series_matches_closed_form": bool(series_residual < tolerance),
        "all_within_tolerance": bool(all_within),
    }


@_scoped_dps
def monotonicity_checks() -> dict[str, Any]:
    """Strict monotonicity, endpoint limits, and the high-energy slope."""

    q = _mpf(SYNTHETIC_Q_OVER_ESTAR)
    ell = _centered(SYNTHETIC_CENTERED_SHAPE)
    sample_mu = ("-12", "-3.25", "0", "2")
    derivatives = [ward_response_derivative(value, q, ell) for value in sample_mu]
    strictly_negative = all(value < 0 for value in derivatives)

    step = mp.mpf("1e-15")
    mu_probe = mp.mpf("-3.25")
    central = (
        ward_response(mu_probe + step, q, ell) - ward_response(mu_probe - step, q, ell)
    ) / (2 * step)
    derivative_consistency_residual = abs(central - derivatives[1])

    mu_large = mp.mpf(40)
    response_large = ward_response(mu_large, q, ell)
    small_z_bound = (
        2 * sum(q**2 * mp.e ** (-2 * (mu_large + value)) for value in ell) / (30 * mp.pi)
    )
    upper_limit_zero = bool(
        response_large <= small_z_bound and response_large < mp.mpf("1e-25")
    )

    mu_low = mp.mpf(-30)
    response_low = ward_response(mu_low, q, ell)
    divergence_bound = (
        2
        * mp.mpf(3)
        / 32
        * sum(
            mp.log(1 + 3 * (q**2 * mp.e ** (-2 * (mu_low + value))) / 16)
            for value in ell
        )
        / mp.pi
    )
    divergence_holds = bool(
        response_low >= divergence_bound and divergence_bound > 10
    )

    slope_residual = abs(derivatives[0] + 2 / mp.pi)
    return {
        "derivative_formula": (
            "dW_Q/dmu = -(4/pi) sum_i integral_0^1 z_i [x(1-x)]^2 "
            "/ (1+z_i x(1-x)) dx"
        ),
        "sample_mu": list(sample_mu),
        "derivatives": [_text(value) for value in derivatives],
        "derivative_strictly_negative": bool(strictly_negative),
        "derivative_vs_central_difference_residual": _text(
            derivative_consistency_residual
        ),
        "derivative_consistent_with_response": bool(
            derivative_consistency_residual < mp.mpf("1e-25")
        ),
        "response_at_mu_40": _text(response_large),
        "small_z_upper_bound": _text(small_z_bound),
        "upper_limit_is_zero": upper_limit_zero,
        "response_at_mu_minus_30": _text(response_low),
        "divergence_lower_bound": _text(divergence_bound),
        "divergence_lower_bound_holds": divergence_holds,
        "high_energy_slope_target": "-2/pi",
        "high_energy_slope_residual_at_mu_minus_12": _text(slope_residual),
        "high_energy_slope_matches_minus_2_over_pi": bool(
            slope_residual < mp.mpf("1e-8")
        ),
    }


@_scoped_dps
def synthetic_round_trip() -> dict[str, Any]:
    """Forward-evaluate a synthetic packet, invert, and verify recovery."""

    q = _mpf(SYNTHETIC_Q_OVER_ESTAR)
    ell = _centered(SYNTHETIC_CENTERED_SHAPE)
    mu_true = _mpf(SYNTHETIC_MU_TRUE)
    response = ward_response(mu_true, q, ell)
    seed = inversion_seed(response, q)
    mu_recovered = solve_mu(response, q, ell, tolerance="1e-50")
    recovery_residual = abs(mu_recovered - mu_true)
    seed_error = abs(seed - mu_recovered)
    line = determinant_line(mu_recovered, ell)
    det_true = mp.e ** (3 * mu_true)
    det_residual = abs(mp.mpf(line["det_line_D_M"]) - det_true) / det_true
    return {
        "packet_label": "synthetic",
        "q_over_Estar": SYNTHETIC_Q_OVER_ESTAR,
        "centered_shape": list(SYNTHETIC_CENTERED_SHAPE),
        "mu_true": SYNTHETIC_MU_TRUE,
        "leptonic_remainder_R_ell": _text(response),
        "seed": _text(seed),
        "seed_formula": "mu_seed = log q - 5/6 - (pi/2) R_ell",
        "seed_absolute_error": _text(seed_error),
        "seed_within_asymptotic_window": bool(seed_error < mp.mpf("0.01")),
        "mu_recovered": _text(mu_recovered),
        "absolute_recovery_residual": _text(recovery_residual),
        "recovered_within_1e-30": bool(recovery_residual < mp.mpf("1e-30")),
        "determinant_line": line,
        "det_line_relative_residual": _text(det_residual),
        "det_line_recovered": bool(det_residual < mp.mpf("1e-29")),
        "response_positive": bool(response > 0),
    }


@_scoped_dps
def interval_enclosure_demo() -> dict[str, Any]:
    """Interval-theorem demo: conservative envelopes bracket the root."""

    q = _mpf(SYNTHETIC_Q_OVER_ESTAR)
    ell = _centered(SYNTHETIC_CENTERED_SHAPE)
    mu_true = _mpf(SYNTHETIC_MU_TRUE)
    response = ward_response(mu_true, q, ell)
    q_half_width = mp.mpf("1e-6")
    ell_half_width = mp.mpf("1e-8")
    response_half_width = mp.mpf("1e-8")
    q_interval = (q - q_half_width, q + q_half_width)
    ell_intervals = [(value - ell_half_width, value + ell_half_width) for value in ell]
    response_interval = (response - response_half_width, response + response_half_width)
    mu_lower, mu_upper, evidence = certify_mu_interval(
        response_interval, q_interval, ell_intervals
    )
    width = mu_upper - mu_lower
    encloses = bool(mu_lower <= mu_true <= mu_upper)
    return {
        "packet_label": "synthetic",
        "theorem": (
            "W_-(mu_-) > R_+ and W_+(mu_+) < R_- imply that every admissible "
            "root lies in [mu_-, mu_+]"
        ),
        "q_interval": [_text(value) for value in q_interval],
        "centered_log_intervals": [
            [_text(pair[0]), _text(pair[1])] for pair in ell_intervals
        ],
        "response_interval": [_text(value) for value in response_interval],
        "mu_interval": [_text(mu_lower), _text(mu_upper)],
        "det_line_interval": [
            _text(mp.e ** (3 * mu_lower)),
            _text(mp.e ** (3 * mu_upper)),
        ],
        "interval_width": _text(width),
        "width_below_1e-5": bool(width < mp.mpf("1e-5")),
        "encloses_synthetic_mu": encloses,
        "certified": True,
        "evidence": evidence,
    }


@_scoped_dps
def build_artifact() -> dict[str, Any]:
    kernel_block = kernel_identity_checks()
    monotonicity_block = monotonicity_checks()
    round_trip_block = synthetic_round_trip()
    interval_block = interval_enclosure_demo()

    checks_pass = (
        kernel_block["all_within_tolerance"]
        and kernel_block["series_matches_closed_form"]
        and monotonicity_block["derivative_strictly_negative"]
        and monotonicity_block["derivative_consistent_with_response"]
        and monotonicity_block["upper_limit_is_zero"]
        and monotonicity_block["divergence_lower_bound_holds"]
        and monotonicity_block["high_energy_slope_matches_minus_2_over_pi"]
        and round_trip_block["recovered_within_1e-30"]
        and round_trip_block["seed_within_asymptotic_window"]
        and round_trip_block["det_line_recovered"]
        and interval_block["certified"]
        and interval_block["encloses_synthetic_mu"]
        and interval_block["width_below_1e-5"]
    )

    return {
        "artifact": "oph_charged_ward_determinant_line",
        "status": "MATHEMATICAL_INVERSION_CLOSED_PHYSICAL_PARENTS_OPEN",
        "public_charged_mass_promotion_allowed": False,
        "charged_reference_data_consumed": False,
        "lane_role": {
            "x2_2_route": "X2.2 electron-absolute-ratio route B producer",
            "synthetic_packet_label": "synthetic",
            "packet_note": (
                "Every numeric packet in this artifact is synthetic and is "
                "labeled synthetic.  The physical source packet is absent."
            ),
        },
        "precision": {"mpmath_dps": WORKING_DPS},
        "kernel_identity": kernel_block,
        "monotonicity": monotonicity_block,
        "synthetic_round_trip": round_trip_block,
        "interval_enclosure": interval_block,
        "physical_source_packet": {
            "status": "ABSENT",
            "gate": "fail_closed",
            "open_parents": [
                {"id": parent_id, "reason": reason}
                for parent_id, reason in OPEN_PARENTS
            ],
            "promotion_rule": (
                "public_charged_mass_promotion_allowed is false because every "
                "listed parent is open."
            ),
        },
        "theorem_registry": {
            "CHARGED_WARD_DETERMINANT_LINE_TO_ELECTRON_SOURCE_RATIO_THEOREM": {
                "statement": (
                    "A source-complete Ward response packet with q > 0, a "
                    "centered shape ell, and a positive leptonic remainder "
                    "R_ell admits exactly one root mu_ch of W_Q(mu; ell) = "
                    "R_ell.  The root fixes det(M_e/E_star) = exp(3 mu_ch) "
                    "and m_i/E_star = exp(mu_ch + ell_i)."
                ),
                "closure_class": "mathematical_inversion",
            },
            "WARD_STRICT_MONOTONICITY_ORBIT_BREAKING_THEOREM": {
                "statement": (
                    "dW_Q/dmu is strictly negative, W_Q tends to 0 as mu "
                    "tends to +infinity and diverges as mu tends to "
                    "-infinity, so W_Q(.; ell) is a bijection from R onto "
                    "(0, infinity).  The common-scale orbit M_e -> "
                    "exp(kappa) M_e moves W_Q, so one Ward endpoint response "
                    "breaks the orbit."
                ),
            },
            "WARD_HIGHER_ORDER_STABILITY_CRITERION": {
                "statement": (
                    "If -W_Q' >= s0 > 0 on an interval J and the "
                    "higher-order remainder derivative obeys |R_ge2'| <= s1 "
                    "< s0 on J, then the full response is strictly "
                    "decreasing on J and monotone uniqueness of the root is "
                    "preserved."
                ),
                "criterion": "s0 > s1",
            },
            "AUDIT_CORRECTION_COMMON_SHIFT_NO_GO_SCOPE": {
                "statement": (
                    "The common-shift no-go holds only on the reduct without "
                    "mass-dependent electromagnetic transport; on that "
                    "reduct every centered readout is invariant under M_e -> "
                    "exp(kappa) M_e.  A Ward endpoint response and a "
                    "determinant basepoint are the two admissible closures "
                    "of the affine coordinate."
                ),
            },
        },
        "checks_pass": bool(checks_pass),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {"status": artifact["status"], "checks_pass": artifact["checks_pass"]},
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
