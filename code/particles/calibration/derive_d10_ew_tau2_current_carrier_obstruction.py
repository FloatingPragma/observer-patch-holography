#!/usr/bin/env python3
"""Certify the current-carrier obstruction beneath direct D10 exact-W/Z emission.

Chain role: prove when the current one-variable D10 carrier cannot emit exact
`W` and `Z` simultaneously through a single `tau2_tree_exact` scalar.

Mathematics: fiberwise neutral-leg law, a current-point affine-germ diagnostic,
and a closed-form nonlinear point test.  The germ sign mismatch is not used as
a global no-go: on the physical domain ``tau2 > -1`` the exact W formula is
strictly increasing, so the W target fixes a unique ``tau2``.  Substitution of
that value into the exact fiberwise Z formula then decides the declared
central-value pair.

OPH-derived inputs: the selected D10 population point and the emitted
fiberwise population tree law on the current carrier.

Output: the smallest obstruction artifact opening the next neutral residual
`delta_n_tree_exact`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import TypeAlias


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_POPULATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"
DEFAULT_FIBERWISE_TREE_LAW = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"


RationalInterval: TypeAlias = tuple[Fraction, Fraction]
PI_MACHIN_TERMS = 40
INPUT_COMPATIBILITY_ABS_TOLERANCE = Fraction(1, 10**16)

EXPECTED_SOURCE_PAIR_ARTIFACT = "oph_d10_ew_source_transport_pair"
EXPECTED_POPULATION_ARTIFACT = "oph_d10_ew_population_evaluator"
EXPECTED_FIBERWISE_ARTIFACT = (
    "oph_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity"
)
EXPECTED_FIBERWISE_FORMULAS = {
    "coordinate_symbol": "tau2_tree_exact",
    "eta_source_formula": "alpha_u_from_seed * beta_EW",
    "fiber_population_functional_formula": (
        "J_pop_EW(tauY,tau2) = tau2^2 + (tauY*tau2)^2 + "
        "(0.5*(tauY + tau2) + eta_source)^2"
    ),
    "fiber_stationarity_formula": (
        "(1 + 4*tau2_tree_exact^2)*tau_Y + tau2_tree_exact + 2*eta_source = 0"
    ),
    "fiber_second_derivative_formula": "1/2 + 2*tau2_tree_exact^2",
    "tauY_formula": (
        "-(tau2_tree_exact + 2*eta_source) / (1 + 4*tau2_tree_exact^2)"
    ),
    "n_EW_formula": (
        "1 + (alphaY_mz * tau_Y + alpha2_mz * tau2_tree_exact) / "
        "(alphaY_mz + alpha2_mz)"
    ),
    "u_EW_formula": "1 + tau2_tree_exact",
    "MW_formula": "v_inherited * sqrt(pi * alpha2_mz * (1 + tau2_tree_exact))",
    "MZ_formula": (
        "v_inherited * sqrt(pi * (alphaY_mz * (1 + tau_Y) + "
        "alpha2_mz * (1 + tau2_tree_exact)))"
    ),
}
EXPECTED_SOURCE_PAIR_FIELDS = {
    "source_pair_symbol": "Tau_EW_D10 = (tau_Y, tau_2)",
    "two_scalar_population_status": "closed_current_carrier",
    "population_selector_status": "closed",
    "population_selector_formula": (
        "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)"
    ),
    "eta_source_formula": "alpha_u_from_seed * beta_EW",
    "predictive_population_closed": True,
}
EXPECTED_POPULATION_FIELDS = {
    "object_id": "EWGaugeSourceTransportPairPopulationEvaluator_D10",
    "population_functional_status": "closed",
    "selector_formula": (
        "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)"
    ),
    "population_selector_rule": (
        "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)"
    ),
    "eta_source_formula": "alpha_u_from_seed * beta_EW",
    "population_functional_formula_sigma_eta": (
        "J_pop_EW(sigma_EW,eta_EW) = (sigma_EW + eta_EW)^2 + "
        "(sigma_EW^2 - eta_EW^2)^2 + (sigma_EW + eta_source)^2"
    ),
    "predictive_promotion_allowed": False,
}
EXPECTED_POPULATION_BASIS_FORMULAS = {
    "beta_EW_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
    "u_EW_formula": "1 + sigma_EW + eta_EW",
    "n_EW_formula": "1 + sigma_EW + beta_EW * eta_EW",
}
EXPECTED_POPULATION_MASS_FORMULAS = {
    "mW_formula": "v_inherited * sqrt(pi * alpha2_mz * u_EW)",
    "mZ_formula": (
        "v_inherited * sqrt(pi * (alphaY_mz + alpha2_mz) * n_EW)"
    ),
}


def _q(value: object) -> Fraction:
    """Interpret an input's serialized decimal value as an exact rational."""

    return Fraction(str(value))


def _payload_digest(payload: object) -> str:
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _require_literal(actual: object, expected: object, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label} must be {expected!r}, got {actual!r}")


def _require_expected_fields(
    payload: dict,
    expected_fields: dict[str, object],
    prefix: str,
) -> None:
    for field, expected in expected_fields.items():
        _require_literal(payload.get(field), expected, f"{prefix}.{field}")


def _require_decimal_close(
    actual: object,
    expected: object,
    label: str,
    *,
    tolerance: Fraction = INPUT_COMPATIBILITY_ABS_TOLERANCE,
) -> None:
    actual_q = _q(actual)
    expected_q = _q(expected)
    if abs(actual_q - expected_q) > tolerance:
        raise ValueError(
            f"{label} differs by {abs(actual_q - expected_q)}; "
            f"maximum allowed serialization drift is {tolerance}"
        )


def _validate_input_contracts(
    source_pair: dict,
    population: dict,
    fiberwise_tree_law: dict,
) -> dict[str, object]:
    """Fail closed unless the three receipts describe one current carrier.

    The obstruction code evaluates a machine-readable specialization of the
    formula strings emitted by the fiber-law producer. Consequently every
    formula, scalar, and anchor used by that specialization is checked here
    before the hard-coded arithmetic implementation is allowed to run. Tiny
    decimal differences created by upstream JSON float serialization are
    accepted only within the declared absolute tolerance; the exact interval
    proof itself uses the fiber law's serialized ``eta_source``.
    """

    _require_literal(
        source_pair.get("artifact"),
        EXPECTED_SOURCE_PAIR_ARTIFACT,
        "source_pair.artifact",
    )
    _require_literal(
        source_pair.get("status"),
        "selected_two_scalar_family",
        "source_pair.status",
    )
    _require_expected_fields(
        source_pair,
        EXPECTED_SOURCE_PAIR_FIELDS,
        "source_pair",
    )
    _require_literal(
        population.get("artifact"),
        EXPECTED_POPULATION_ARTIFACT,
        "population.artifact",
    )
    _require_literal(
        population.get("status"),
        "closed_current_carrier",
        "population.status",
    )
    _require_literal(
        population.get("proof_status"),
        "population_functional_closed_on_current_carrier",
        "population.proof_status",
    )
    _require_literal(
        population.get("source_transport_pair_artifact"),
        EXPECTED_SOURCE_PAIR_ARTIFACT,
        "population.source_transport_pair_artifact",
    )
    _require_expected_fields(
        population,
        EXPECTED_POPULATION_FIELDS,
        "population",
    )
    _require_literal(
        fiberwise_tree_law.get("artifact"),
        EXPECTED_FIBERWISE_ARTIFACT,
        "fiberwise_tree_law.artifact",
    )
    _require_literal(
        fiberwise_tree_law.get("status"),
        "closed_smaller_primitive",
        "fiberwise_tree_law.status",
    )
    _require_literal(
        fiberwise_tree_law.get("proof_status"),
        "fiberwise_unique_J_pop_minimizer_on_fixed_tau2",
        "fiberwise_tree_law.proof_status",
    )
    _require_expected_fields(
        fiberwise_tree_law,
        EXPECTED_FIBERWISE_FORMULAS,
        "fiberwise_tree_law",
    )

    source_basis = dict(source_pair["population_basis"])
    population_basis = dict(population["population_basis"])
    source_mass_formulas = dict(source_pair["population_atomic_quartet"])
    population_mass_formulas = dict(population["population_atomic_quartet"])
    for prefix, payload, expected in (
        (
            "source_pair.population_basis",
            source_basis,
            EXPECTED_POPULATION_BASIS_FORMULAS,
        ),
        (
            "population.population_basis",
            population_basis,
            EXPECTED_POPULATION_BASIS_FORMULAS,
        ),
        (
            "source_pair.population_atomic_quartet",
            source_mass_formulas,
            EXPECTED_POPULATION_MASS_FORMULAS,
        ),
        (
            "population.population_atomic_quartet",
            population_mass_formulas,
            EXPECTED_POPULATION_MASS_FORMULAS,
        ),
    ):
        _require_expected_fields(payload, expected, prefix)

    source_slots = dict(source_pair["source_pair"])
    duplicate_source_slots = dict(source_pair["source_slots"])
    fiber_slots = dict(fiberwise_tree_law["carrier_basis_scalar"])
    for field in ("alphaY_mz", "alpha2_mz", "v_inherited"):
        if _q(source_slots[field]) <= 0:
            raise ValueError(f"source_pair.source_pair.{field} must be positive")
        _require_decimal_close(
            duplicate_source_slots[field],
            source_slots[field],
            f"source_pair.source_slots.{field}",
        )
        _require_decimal_close(
            fiber_slots[field],
            source_slots[field],
            f"fiberwise_tree_law.carrier_basis_scalar.{field}",
        )

    alpha_y = _q(source_slots["alphaY_mz"])
    alpha2 = _q(source_slots["alpha2_mz"])
    beta_ew = (alpha2 - alpha_y) / (alpha2 + alpha_y)
    _require_decimal_close(
        source_basis["beta_EW"],
        beta_ew,
        "source_pair.population_basis.beta_EW",
    )
    _require_decimal_close(
        population_basis["beta_EW"],
        beta_ew,
        "population.population_basis.beta_EW",
    )
    _require_decimal_close(
        fiber_slots["beta_EW"],
        beta_ew,
        "fiberwise_tree_law.carrier_basis_scalar.beta_EW",
    )

    eta_fiber = _q(fiberwise_tree_law["eta_source"])
    anchor = dict(fiberwise_tree_law["anchor_point"])
    _require_decimal_close(
        anchor["eta_EW"],
        eta_fiber,
        "fiberwise_tree_law.anchor_point.eta_EW",
    )
    _require_literal(
        _q(anchor["tau2_tree_exact"]),
        Fraction(0),
        "fiberwise_tree_law.anchor_point.tau2_tree_exact",
    )
    _require_decimal_close(
        anchor["sigma_EW"],
        -eta_fiber,
        "fiberwise_tree_law.anchor_point.sigma_EW",
    )
    _require_decimal_close(
        anchor["tau_Y"],
        -2 * eta_fiber,
        "fiberwise_tree_law.anchor_point.tau_Y",
    )

    selected_point = dict(population["selected_population_point"])
    _require_decimal_close(
        population["eta_source"],
        eta_fiber,
        "population.eta_source",
    )
    _require_decimal_close(
        source_pair["eta_source"],
        eta_fiber,
        "source_pair.eta_source",
    )
    for field, expected in (
        ("eta_EW", eta_fiber),
        ("sigma_EW", -eta_fiber),
        ("tau_2", Fraction(0)),
        ("tau_Y", -2 * eta_fiber),
    ):
        _require_decimal_close(
            selected_point[field],
            expected,
            f"population.selected_population_point.{field}",
        )

    selected_basis = dict(population["selected_population_basis_point"])
    expected_n = Fraction(1) - (Fraction(1) - beta_ew) * eta_fiber
    _require_decimal_close(
        selected_basis["u_EW"],
        Fraction(1),
        "population.selected_population_basis_point.u_EW",
    )
    _require_decimal_close(
        selected_basis["n_EW"],
        expected_n,
        "population.selected_population_basis_point.n_EW",
    )

    derived_from = fiberwise_tree_law.get("derived_from_artifacts", [])
    if EXPECTED_POPULATION_ARTIFACT not in derived_from:
        raise ValueError(
            "fiberwise_tree_law.derived_from_artifacts must include the "
            "population evaluator"
        )

    return {
        "status": "PASS",
        "binding": "fail_closed_formula_scalar_anchor_validation",
        "canonical_eta_source": "fiberwise_tree_law.eta_source",
        "serialization_compatibility_absolute_tolerance": _decimal(
            INPUT_COMPATIBILITY_ABS_TOLERANCE
        ),
        "validated_artifacts": [
            EXPECTED_SOURCE_PAIR_ARTIFACT,
            EXPECTED_POPULATION_ARTIFACT,
            EXPECTED_FIBERWISE_ARTIFACT,
        ],
        "validated_formula_fields": sorted(EXPECTED_FIBERWISE_FORMULAS),
        "validated_formula_contracts": {
            "source_pair": sorted(EXPECTED_SOURCE_PAIR_FIELDS),
            "population": sorted(EXPECTED_POPULATION_FIELDS),
            "population_basis": sorted(EXPECTED_POPULATION_BASIS_FORMULAS),
            "population_mass_map": sorted(EXPECTED_POPULATION_MASS_FORMULAS),
        },
    }


def _declared_path(path: Path) -> str:
    """Render an input path relative to the scientific repo when possible,
    always with forward slashes so the artifact is platform independent."""

    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def _point(value: Fraction) -> RationalInterval:
    return value, value


def _interval_add(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    return left[0] + right[0], left[1] + right[1]


def _interval_neg(value: RationalInterval) -> RationalInterval:
    return -value[1], -value[0]


def _interval_sub(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    return _interval_add(left, _interval_neg(right))


def _interval_mul(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    products = tuple(a * b for a in left for b in right)
    return min(products), max(products)


def _interval_reciprocal(value: RationalInterval) -> RationalInterval:
    if value[0] <= 0 <= value[1]:
        raise ValueError("cannot invert an interval containing zero")
    endpoints = Fraction(1, value[0]), Fraction(1, value[1])
    return min(endpoints), max(endpoints)


def _interval_div(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    return _interval_mul(left, _interval_reciprocal(right))


def _arctan_bounds(unit_denominator: int, terms: int) -> RationalInterval:
    """Exact alternating-series enclosure of arctan(1/unit_denominator)."""

    x = Fraction(1, unit_denominator)
    partial = Fraction(0)
    for k in range(terms):
        term = x ** (2 * k + 1) / (2 * k + 1)
        partial += term if k % 2 == 0 else -term
    next_term = x ** (2 * terms + 1) / (2 * terms + 1)
    next_signed = next_term if terms % 2 == 0 else -next_term
    other_endpoint = partial + next_signed
    return min(partial, other_endpoint), max(partial, other_endpoint)


def certified_pi_bounds(terms: int = PI_MACHIN_TERMS) -> RationalInterval:
    """Enclose pi using Machin's identity and exact rational series bounds."""

    if terms < 2:
        raise ValueError("at least two Machin-series terms are required")
    atan_fifth = _arctan_bounds(5, terms)
    atan_239th = _arctan_bounds(239, terms)
    # Machin identity: pi = 16*atan(1/5) - 4*atan(1/239).
    return (
        16 * atan_fifth[0] - 4 * atan_239th[1],
        16 * atan_fifth[1] - 4 * atan_239th[0],
    )


def _decimal(value: Fraction, precision: int = 70) -> str:
    with localcontext() as context:
        context.prec = precision
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def _fraction_payload(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal_approximation": _decimal(value),
    }


def _interval_payload(value: RationalInterval) -> dict[str, object]:
    return {
        "lower": _fraction_payload(value[0]),
        "upper": _fraction_payload(value[1]),
    }


def _exact_nonlinear_interval_certificate(
    *,
    alpha_y: Fraction,
    alpha2: Fraction,
    v_value: Fraction,
    eta_source: Fraction,
    mw_target: Fraction,
    mz_target: Fraction,
    pi_bounds: RationalInterval,
) -> dict[str, object]:
    """Bound MZ(tau_W)^2-MZ_target^2 with exact rational arithmetic."""

    if not (alpha_y > 0 and alpha2 > 0 and v_value > 0):
        raise ValueError("electroweak couplings and inherited scale must be positive")
    if not (mw_target > 0 and mz_target > 0):
        raise ValueError("mass targets must be positive")
    if not (pi_bounds[0] > 0 and pi_bounds[0] <= pi_bounds[1]):
        raise ValueError("pi bounds must be ordered and positive")

    one = _point(Fraction(1))
    two = _point(Fraction(2))
    four = _point(Fraction(4))
    alpha_y_i = _point(alpha_y)
    alpha2_i = _point(alpha2)
    v_i = _point(v_value)
    eta_i = _point(eta_source)
    mw_i = _point(mw_target)
    mz_i = _point(mz_target)

    v_squared = _interval_mul(v_i, v_i)
    tau_w = _interval_sub(
        _interval_div(
            _interval_mul(mw_i, mw_i),
            _interval_mul(_interval_mul(v_squared, pi_bounds), alpha2_i),
        ),
        one,
    )
    if tau_w[0] <= -1:
        raise ValueError("certified W preimage leaves the tau2 > -1 domain")
    tau_w_squared = _interval_mul(tau_w, tau_w)
    tau_y = _interval_neg(
        _interval_div(
            _interval_add(tau_w, _interval_mul(two, eta_i)),
            _interval_add(one, _interval_mul(four, tau_w_squared)),
        )
    )
    n_fiber = _interval_add(
        one,
        _interval_div(
            _interval_add(
                _interval_mul(alpha_y_i, tau_y),
                _interval_mul(alpha2_i, tau_w),
            ),
            _interval_add(alpha_y_i, alpha2_i),
        ),
    )
    if n_fiber[0] <= 0:
        raise ValueError("interval cannot certify the positive MZ fiber domain")
    mz_squared = _interval_mul(
        _interval_mul(
            _interval_mul(v_squared, pi_bounds),
            _interval_add(alpha_y_i, alpha2_i),
        ),
        n_fiber,
    )
    residual_squared = _interval_sub(mz_squared, _interval_mul(mz_i, mz_i))
    if residual_squared[1] < 0:
        sign = "strictly_negative"
    elif residual_squared[0] > 0:
        sign = "strictly_positive"
    else:
        sign = "undetermined_interval_contains_zero"

    return {
        "arithmetic": "exact_rational_interval_arithmetic",
        "input_interpretation": (
            "each serialized decimal input is interpreted as the exact rational "
            "with that finite decimal expansion"
        ),
        "tested_expression": "MZ_fiber(tau2_from_exact_W)^2 - MZ_target^2",
        "tau2_from_exact_W_interval": _interval_payload(tau_w),
        "tauY_fiber_at_exact_W_interval": _interval_payload(tau_y),
        "n_EW_fiber_at_exact_W_interval": _interval_payload(n_fiber),
        "MZ_squared_at_exact_W_interval_gev2": _interval_payload(mz_squared),
        "MZ_squared_residual_interval_gev2": _interval_payload(residual_squared),
        "residual_squared_sign": sign,
        "zero_excluded": sign in {"strictly_negative", "strictly_positive"},
    }


def _fiber_tau_y(tau2: float, eta_source: float) -> float:
    return -(tau2 + 2.0 * eta_source) / (1.0 + 4.0 * tau2 * tau2)


def _fiber_n(
    tau2: float,
    eta_source: float,
    alpha_y: float,
    alpha2: float,
) -> float:
    tau_y = _fiber_tau_y(tau2, eta_source)
    return 1.0 + (alpha_y * tau_y + alpha2 * tau2) / (alpha_y + alpha2)


def _mw(tau2: float, v_value: float, alpha2: float) -> float:
    if tau2 <= -1.0:
        raise ValueError("MW current-carrier domain requires tau2 > -1")
    return v_value * math.sqrt(math.pi * alpha2 * (1.0 + tau2))


def _mz_fiber(
    tau2: float,
    eta_source: float,
    alpha_y: float,
    alpha2: float,
    v_value: float,
) -> float:
    n_value = _fiber_n(tau2, eta_source, alpha_y, alpha2)
    if n_value <= 0.0:
        raise ValueError("MZ current-carrier domain requires n_EW_fiber > 0")
    return v_value * math.sqrt(math.pi * (alpha_y + alpha2) * n_value)


def build_artifact(
    source_pair: dict,
    population: dict,
    fiberwise_tree_law: dict,
    references: dict,
    *,
    pi_bounds: RationalInterval | None = None,
    source_pair_path: Path = DEFAULT_SOURCE_PAIR,
    population_path: Path = DEFAULT_POPULATION,
    fiberwise_tree_law_path: Path = DEFAULT_FIBERWISE_TREE_LAW,
    reference_path: Path = REFERENCE_JSON,
) -> dict:
    input_contract_validation = _validate_input_contracts(
        source_pair,
        population,
        fiberwise_tree_law,
    )
    selected_point = dict(population.get("selected_population_point", {}))
    if not selected_point:
        raise ValueError("selected population point is required")
    source_slots = dict(source_pair["source_pair"])
    alpha_y = float(source_slots["alphaY_mz"])
    alpha2 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    eta_source = float(fiberwise_tree_law["eta_source"])
    beta_ew = (alpha2 - alpha_y) / (alpha2 + alpha_y)
    n0 = 1.0 - (1.0 - beta_ew) * eta_source

    mw_current = _mw(0.0, v_value, alpha2)
    mz_current = _mz_fiber(0.0, eta_source, alpha_y, alpha2, v_value)
    mw_exact = float(references["w_boson"]["value_gev"])
    mz_exact = float(references["z_boson"]["value_gev"])
    mw_sigma = float(references["w_boson"]["error_plus_gev"])
    mz_sigma = float(references["z_boson"]["error_plus_gev"])

    # First-order germ coefficients at the current point and the tau2 value each
    # leg would need on its own. A single tau2 serves both legs at first order
    # only when the two required values share a sign.
    germ_w = 0.5
    germ_z = beta_ew / (2.0 * n0)
    w_offset = mw_current - mw_exact
    z_offset = mz_current - mz_exact
    tau2_required_w = (mw_exact / mw_current - 1.0) / germ_w
    tau2_required_z = (mz_exact / mz_current - 1.0) / germ_z
    single_tau2_possible_at_first_order = (
        (tau2_required_w > 0.0) == (tau2_required_z > 0.0)
    )

    # Exact nonlinear point test.  MW(tau2) is strictly increasing for
    # tau2 > -1, hence its central target has the unique closed-form preimage
    # below.  Testing MZ at that same coordinate is sufficient to decide
    # simultaneous attainment of the two declared central values on this
    # current carrier.  This does not constrain a modified fiber law or an
    # additional neutral-leg coordinate.
    tau2_from_exact_w = (mw_exact / mw_current) ** 2 - 1.0
    tau_y_at_exact_w = _fiber_tau_y(tau2_from_exact_w, eta_source)
    n_at_exact_w = _fiber_n(tau2_from_exact_w, eta_source, alpha_y, alpha2)
    mw_at_exact_w = _mw(tau2_from_exact_w, v_value, alpha2)
    mz_at_exact_w = _mz_fiber(
        tau2_from_exact_w,
        eta_source,
        alpha_y,
        alpha2,
        v_value,
    )
    mw_residual_at_exact_w = mw_at_exact_w - mw_exact
    mz_residual_at_exact_w = mz_at_exact_w - mz_exact
    machin_pi_bounds = certified_pi_bounds()
    selected_pi_bounds = machin_pi_bounds if pi_bounds is None else pi_bounds
    if not (
        selected_pi_bounds[0] <= machin_pi_bounds[0]
        and selected_pi_bounds[1] >= machin_pi_bounds[1]
    ):
        raise ValueError(
            "supplied pi bounds must contain the certified Machin-series enclosure"
        )
    exact_interval = _exact_nonlinear_interval_certificate(
        alpha_y=_q(source_slots["alphaY_mz"]),
        alpha2=_q(source_slots["alpha2_mz"]),
        v_value=_q(source_slots["v_inherited"]),
        eta_source=_q(fiberwise_tree_law["eta_source"]),
        mw_target=_q(references["w_boson"]["value_gev"]),
        mz_target=_q(references["z_boson"]["value_gev"]),
        pi_bounds=selected_pi_bounds,
    )
    obstruction_established = bool(exact_interval["zero_excluded"])
    simultaneous_exact_pair_possible: bool | None = (
        False if obstruction_established else None
    )

    status = (
        "closed_smaller_primitive"
        if obstruction_established
        else "obstruction_not_established_interval_contains_zero"
    )
    proof_status = (
        "exact_interval_excludes_single_tau2_central_WZ_pair_on_current_carrier"
        if obstruction_established
        else "exact_interval_does_not_exclude_single_tau2_central_WZ_pair"
    )

    return {
        "artifact": "oph_d10_ew_tau2_current_carrier_obstruction",
        "object_id": "EWCurrentCarrierTau2Obstruction_D10",
        "status": status,
        "proof_status": proof_status,
        "strictly_smaller_than": "tau2_tree_exact",
        "diagnostic_only": True,
        "input_contract_validation": input_contract_validation,
        "input_provenance": {
            "digest_scope": (
                "sha256 of canonical sorted compact JSON for the exact "
                "in-memory payload consumed by this producer"
            ),
            "source_pair": {
                "declared_path": _declared_path(source_pair_path),
                "payload_sha256": _payload_digest(source_pair),
            },
            "population": {
                "declared_path": _declared_path(population_path),
                "payload_sha256": _payload_digest(population),
            },
            "fiberwise_tree_law": {
                "declared_path": _declared_path(fiberwise_tree_law_path),
                "payload_sha256": _payload_digest(fiberwise_tree_law),
            },
            "reference_entries": {
                "declared_path": _declared_path(reference_path),
                "subobject": "entries",
                "payload_sha256": _payload_digest(references),
            },
        },
        "population_evaluator_artifact": population.get("artifact"),
        "fiberwise_population_tree_law_artifact": fiberwise_tree_law.get("artifact"),
        "coordinate_symbol": "tau2_tree_exact",
        "eta_source": eta_source,
        "beta_EW": beta_ew,
        "n0": n0,
        "fiberwise_tauY_formula": "-(tau2_tree_exact + 2*eta_source) / (1 + 4*tau2_tree_exact^2)",
        "n_EW_fiber_formula": "1 + (alphaY_mz * fiberwise_tauY + alpha2_mz * tau2_tree_exact) / (alphaY_mz + alpha2_mz)",
        "MW_formula": "v_inherited * sqrt(pi * alpha2_mz * (1 + tau2_tree_exact))",
        "MZ_fiber_formula": "v_inherited * sqrt(pi * (alphaY_mz + alpha2_mz) * n_EW_fiber)",
        "local_affine_germ": {
            "n0_formula": "1 - (1 - beta_EW) * eta_source",
            "MW_relative_formula": "1 + 0.5 * tau2_tree_exact + O(tau2_tree_exact^2)",
            "MZ_relative_formula": "1 + (beta_EW / (2 * n0)) * tau2_tree_exact + O(tau2_tree_exact^2)",
        },
        "current_point": {
            "MW_pole": mw_current,
            "MZ_pole": mz_current,
        },
        "reference_targets": {
            "MW_pole": mw_exact,
            "MZ_pole": mz_exact,
        },
        "proof_scope": (
            "exact_nonlinear_current_carrier_formulas_against_current_reference_"
            "central_values; no uncertainty-band, modified-fiber-law, or future-"
            "extension impossibility is claimed"
        ),
        "exactness_convention": {
            "status": "conditional_exact_diagnostic",
            "upstream_decimal_interpretation": (
                "the finite-decimal fields present in the consumed payloads "
                "are treated as exact rationals; this does not assert that "
                "the underlying measured or derived physical quantities are exact"
            ),
            "scientific_scope": (
                "diagnostic obstruction for the declared central pair and "
                "current carrier formulas only"
            ),
        },
        "closed_form_nonlinear_point_test": {
            "domain": "tau2_tree_exact > -1 and n_EW_fiber(tau2_tree_exact) > 0",
            "uniqueness_argument": (
                "MW(tau2) = MW(0)*sqrt(1+tau2) is strictly increasing on "
                "tau2 > -1, so exact MW fixes the displayed unique tau2"
            ),
            "analytic_W_preimage_status": "exact_unique_on_tau2_tree_exact_greater_than_minus_one",
            "tau2_from_exact_W": tau2_from_exact_w,
            "tauY_fiber_at_exact_W": tau_y_at_exact_w,
            "n_EW_fiber_at_exact_W": n_at_exact_w,
            "MW_at_unique_W_coordinate_gev": mw_at_exact_w,
            "MZ_at_unique_W_coordinate_gev": mz_at_exact_w,
            "MW_residual_gev": mw_residual_at_exact_w,
            "MZ_residual_gev": mz_residual_at_exact_w,
            "floating_point_values_role": (
                "display-only central-value diagnostics; neither the binary-float "
                "residual nor a numerical tolerance is used in the proof"
            ),
            "pi_enclosure": {
                "method": (
                    "Machin identity pi=16*atan(1/5)-4*atan(1/239), with "
                    "alternating-series remainder bounds evaluated as exact rationals"
                ),
                "terms_per_arctan_series": PI_MACHIN_TERMS,
                "canonical_machin_bounds": _interval_payload(machin_pi_bounds),
                "bounds_used": _interval_payload(selected_pi_bounds),
            },
            "exact_interval_certificate": exact_interval,
            "simultaneous_exact_central_pair_possible": simultaneous_exact_pair_possible,
            "obstruction_established": obstruction_established,
            "nonlinear_routes_left_open": [
                "an additional neutral-leg scalar delta_n_tree_exact",
                "a source-derived modification of the fiberwise tauY law",
                "a different carrier or selector",
                "fits to reference uncertainty regions rather than the two central values",
            ],
        },
        "direction_obstruction": {
            "scope": "first_order_affine_germ_at_tau2_tree_exact=0_only",
            "W_current_minus_exact_sign": "positive" if w_offset > 0 else "negative" if w_offset < 0 else "zero",
            "Z_current_minus_exact_sign": "positive" if z_offset > 0 else "negative" if z_offset < 0 else "zero",
            "germ_coefficient_W": germ_w,
            "germ_coefficient_Z": germ_z,
            "tau2_required_for_W_first_order": tau2_required_w,
            "tau2_required_for_Z_first_order": tau2_required_z,
            "single_tau2_possible_at_first_order": single_tau2_possible_at_first_order,
            "global_no_go_inferred_from_germ": False,
            "evaluation": (
                "computed from the two offset signs and germ coefficient signs "
                "at the current point; the separate closed-form nonlinear point "
                "test, not this germ, establishes the central-value obstruction"
            ),
        },
        "reference_distance": {
            "W_offset_gev": w_offset,
            "Z_offset_gev": z_offset,
            "W_reference_sigma_gev": mw_sigma,
            "Z_reference_sigma_gev": mz_sigma,
            "W_offset_sigma": w_offset / mw_sigma,
            "Z_offset_sigma": z_offset / mz_sigma,
        },
        "minimal_extra_scalar_or_invariant": {
            "symbol": "delta_n_tree_exact",
            "equivalent_source_scalar": "delta_tauY_tree_exact",
            "definition_formula": "n_EW_exact = n_EW_fiber + delta_n_tree_exact",
            "equivalent_tauY_formula": "tauY_exact = fiberwise_tauY + ((alphaY_mz + alpha2_mz) / alphaY_mz) * delta_n_tree_exact",
        },
        "next_single_residual_object": (
            "delta_n_tree_exact" if obstruction_established else None
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 current-carrier tau2 obstruction artifact.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--fiberwise-tree-law", default=str(DEFAULT_FIBERWISE_TREE_LAW))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    source_pair_path = Path(args.source_pair)
    population_path = Path(args.population)
    fiberwise_tree_law_path = Path(args.fiberwise_tree_law)
    source_pair = json.loads(source_pair_path.read_text(encoding="utf-8"))
    population = json.loads(population_path.read_text(encoding="utf-8"))
    fiberwise_tree_law = json.loads(fiberwise_tree_law_path.read_text(encoding="utf-8"))
    artifact = build_artifact(
        source_pair,
        population,
        fiberwise_tree_law,
        references,
        source_pair_path=source_pair_path,
        population_path=population_path,
        fiberwise_tree_law_path=fiberwise_tree_law_path,
        reference_path=REFERENCE_JSON,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
