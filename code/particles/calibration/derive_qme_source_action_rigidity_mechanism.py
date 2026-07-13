#!/usr/bin/env python3
"""Quotient-MaxEnt source-action rigidity mechanism receipt.

An explicit two-regulator source packet (Q0 = {0,1}, Q1 = {00,01,10,11},
coarse map (x,y) -> x, uniform reference measures, Bernoulli feature on the
coarse level) carries the finite quotient-MaxEnt selection mechanism at the
audited pixel branch moment p = (P - phi)/sqrt(pi).  The builder verifies the
unique MaxEnt laws, the exponential-family multiplier, exact projective
compatibility under the coarse map, the Pythagorean identity, a strictly
positive Pinsker selector gap, and the Legendre Hessian 1/(p(1-p)).

This artifact sharpens the Local MaxEnt axiom into a quantitative selection
principle and is the mechanism receipt for the D10 discrete law decision.
The physical Standard-Model moment packet is absent and the promotion gate
is fail-closed.
"""

from __future__ import annotations

import functools

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp


WORKING_DPS = 60


def _scoped_dps(func):
    """Run func at WORKING_DPS and restore the global precision on exit."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with mp.workdps(WORKING_DPS):
            return func(*args, **kwargs)

    return wrapper


ROOT = Path(__file__).resolve().parents[2]
PIXEL_CERTIFICATE = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "qme_source_action_rigidity_mechanism.json"
)

EXPECTED_P_BERNOULLI = "0.007299545168876648"
P_AGREEMENT_TOLERANCE = "1e-15"
RESIDUAL_TOLERANCE = "1e-40"

R1_STATES = ("00", "01", "10", "11")
R1_FEATURE_X = (0, 0, 1, 1)
R1_FEATURE_Y = (0, 1, 0, 1)
COARSE_MAP = {"00": "0", "01": "0", "10": "1", "11": "1"}

OPEN_PROOFS = (
    "operator completeness modulo BRST",
    "c_r = c_r(P_star, N_star) source emission",
    "reflection positivity and continuum nontriviality",
    "source-complete matter spectral cuts",
    "physical 1PI two-point blocks with analytic continuation",
    "source clock closure",
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


def load_certificate(path: Path = PIXEL_CERTIFICATE) -> dict[str, Any]:
    """Load the pixel certificate with float literals kept as strings."""

    return json.loads(path.read_text(encoding="utf-8"), parse_float=str)


def golden_ratio() -> mp.mpf:
    return (1 + mp.sqrt(5)) / 2


def bernoulli_moment(p_cand: Any) -> mp.mpf:
    """Source moment p = (P - phi)/sqrt(pi) from the audited branch value P."""

    p = (_mpf(p_cand) - golden_ratio()) / mp.sqrt(mp.pi)
    if not 0 < p < 1:
        raise ValueError("the Bernoulli source moment must lie in (0, 1)")
    return p


def relative_entropy(law: Sequence[Any], reference: Sequence[Any]) -> mp.mpf:
    """D(law || reference) in natural logarithm units."""

    if len(law) != len(reference):
        raise ValueError("relative entropy requires equal-length laws")
    total = mp.mpf(0)
    for a, b in zip(law, reference):
        a, b = _mpf(a), _mpf(b)
        if a < 0 or b <= 0:
            raise ValueError("relative entropy requires a >= 0 and b > 0")
        if a > 0:
            total += a * mp.log(a / b)
    return total


def l1_distance(left: Sequence[Any], right: Sequence[Any]) -> mp.mpf:
    return sum(abs(_mpf(a) - _mpf(b)) for a, b in zip(left, right))


def maxent_law_r0(p: Any) -> tuple[mp.mpf, mp.mpf]:
    """Unique MaxEnt law on Q0 = {0,1} with uniform reference and E[X] = p."""

    p = _mpf(p)
    return (1 - p, p)


def maxent_law_r1(p: Any) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """Unique MaxEnt law on Q1 with E[X] = p and E[Y] = 1/2."""

    p = _mpf(p)
    return ((1 - p) / 2, (1 - p) / 2, p / 2, p / 2)


def multiplier_theta(p: Any) -> mp.mpf:
    """Exponential-family multiplier theta_star = log((1-p)/p)."""

    p = _mpf(p)
    return mp.log((1 - p) / p)


def gibbs_law(
    reference: Sequence[Any],
    features: Sequence[Sequence[Any]],
    theta: Sequence[Any],
) -> tuple[mp.mpf, ...]:
    """Normalized Gibbs law mu(q) = m(q) exp(-theta . F(q)) / Z."""

    weights = []
    for m_value, feature_row in zip(reference, features):
        exponent = -sum(_mpf(t) * _mpf(f) for t, f in zip(theta, feature_row))
        weights.append(_mpf(m_value) * mp.e**exponent)
    partition = sum(weights)
    return tuple(weight / partition for weight in weights)


def coarse_pushforward(mu1: Sequence[Any]) -> tuple[mp.mpf, mp.mpf]:
    """Pushforward of a Q1 law under the coarse map (x,y) -> x."""

    values = [_mpf(value) for value in mu1]
    return (values[0] + values[1], values[2] + values[3])


def legendre_hessian(p: Any) -> mp.mpf:
    """Coarse-level Legendre Hessian 1/(p(1-p)) = Var(X)^-1."""

    p = _mpf(p)
    return 1 / (p * (1 - p))


@_scoped_dps
def mechanism_checks(p: Any) -> dict[str, Any]:
    """All rigidity-mechanism computations at the supplied moment p."""

    p = _mpf(p)
    tolerance = mp.mpf(RESIDUAL_TOLERANCE)
    half = mp.mpf("0.5")
    mu0 = maxent_law_r0(p)
    mu1 = maxent_law_r1(p)
    m0 = (half, half)
    m1 = (mp.mpf("0.25"),) * 4

    mu0_sum_residual = abs(sum(mu0) - 1)
    mu1_sum_residual = abs(sum(mu1) - 1)
    r0_moment_residual = abs(mu0[1] - p)
    r1_x_moment_residual = abs((mu1[2] + mu1[3]) - p)
    r1_y_moment_residual = abs((mu1[1] + mu1[3]) - half)

    theta = multiplier_theta(p)
    gibbs_r0 = gibbs_law(m0, ((0,), (1,)), (theta,))
    gibbs_r1 = gibbs_law(
        m1, tuple(zip(R1_FEATURE_X, R1_FEATURE_Y)), (theta, mp.mpf(0))
    )
    gibbs_r0_residual = max(abs(a - b) for a, b in zip(gibbs_r0, mu0))
    gibbs_r1_residual = max(abs(a - b) for a, b in zip(gibbs_r1, mu1))

    push = coarse_pushforward(mu1)
    projective_residual = max(abs(a - b) for a, b in zip(push, mu0))

    displacement = min(p, 1 - p) / 8
    nu1 = (
        mu1[0] + displacement,
        mu1[1] - displacement,
        mu1[2] - displacement,
        mu1[3] + displacement,
    )
    nu_full_support = bool(min(nu1) > 0)
    nu_x_moment_residual = abs((nu1[2] + nu1[3]) - p)
    nu_y_moment_residual = abs((nu1[1] + nu1[3]) - half)
    lhs = relative_entropy(nu1, m1) - relative_entropy(mu1, m1)
    rhs = relative_entropy(nu1, mu1)
    pythagorean_residual = abs(lhs - rhs)

    delta = l1_distance(nu1, mu1)
    selector_gap = rhs
    pinsker_bound = delta**2 / 2
    gap_exceeds_bound = bool(selector_gap >= pinsker_bound and pinsker_bound > 0)

    hessian = legendre_hessian(p)

    return {
        "p": _text(p),
        "theta_star": _text(theta),
        "theta_star_formula": "theta_star = log((1-p)/p)",
        "mu0_star": [_text(value) for value in mu0],
        "mu1_star": [_text(value) for value in mu1],
        "mu0_sum_residual": _text(mu0_sum_residual),
        "mu1_sum_residual": _text(mu1_sum_residual),
        "laws_normalized": bool(
            mu0_sum_residual < tolerance and mu1_sum_residual < tolerance
        ),
        "moment_residuals": {
            "r0_E_X": _text(r0_moment_residual),
            "r1_E_X": _text(r1_x_moment_residual),
            "r1_E_Y": _text(r1_y_moment_residual),
        },
        "moment_constraints_hold": bool(
            r0_moment_residual < tolerance
            and r1_x_moment_residual < tolerance
            and r1_y_moment_residual < tolerance
        ),
        "gibbs_reconstruction_residuals": {
            "r0": _text(gibbs_r0_residual),
            "r1": _text(gibbs_r1_residual),
        },
        "gibbs_reconstruction_matches": bool(
            gibbs_r0_residual < tolerance and gibbs_r1_residual < tolerance
        ),
        "projective_residual": _text(projective_residual),
        "projective_compatibility_exact": bool(
            projective_residual == 0 or projective_residual < mp.mpf("1e-30")
        ),
        "perturbation_displacement": _text(displacement),
        "perturbed_feasible_law_nu1": [_text(value) for value in nu1],
        "perturbation_feasible_full_support": bool(
            nu_full_support
            and nu_x_moment_residual < tolerance
            and nu_y_moment_residual < tolerance
        ),
        "pythagorean_residual": _text(pythagorean_residual),
        "pythagorean_identity_holds": bool(pythagorean_residual < mp.mpf("1e-30")),
        "l1_distance_delta": _text(delta),
        "selector_gap": _text(selector_gap),
        "pinsker_lower_bound": _text(pinsker_bound),
        "pinsker_bound_formula": "Delta(delta) >= delta^2/2",
        "selector_gap_exceeds_pinsker_bound": gap_exceeds_bound,
        "Gamma_hessian": _text(hessian),
        "Gamma_hessian_formula": "1/(p(1-p))",
    }


@_scoped_dps
def build_artifact(certificate: Mapping[str, Any]) -> dict[str, Any]:
    p_cand_raw = certificate["P_cand"]
    phi = golden_ratio()
    p = bernoulli_moment(p_cand_raw)
    p_agreement_residual = abs(p - mp.mpf(EXPECTED_P_BERNOULLI))
    p_agrees = bool(p_agreement_residual < mp.mpf(P_AGREEMENT_TOLERANCE))

    checks_block = mechanism_checks(p)

    checks_pass = (
        p_agrees
        and checks_block["laws_normalized"]
        and checks_block["moment_constraints_hold"]
        and checks_block["gibbs_reconstruction_matches"]
        and checks_block["projective_compatibility_exact"]
        and checks_block["perturbation_feasible_full_support"]
        and checks_block["pythagorean_identity_holds"]
        and checks_block["selector_gap_exceeds_pinsker_bound"]
    )

    return {
        "artifact": "oph_qme_source_action_rigidity_mechanism",
        "status": "SELECTION_MECHANISM_CLOSED_PHYSICAL_MOMENT_VECTOR_OPEN",
        "promotion_allowed": False,
        "axiom_role": {
            "sharpens": "Local MaxEnt axiom",
            "into": (
                "a quantitative selection principle with a strictly positive "
                "Pinsker selector gap"
            ),
            "d10_role": "mechanism receipt for the D10 discrete law decision",
        },
        "branch": {
            "P_cand": str(p_cand_raw),
            "P_cand_source": str(PIXEL_CERTIFICATE.relative_to(ROOT)),
            "phi": _text(phi),
            "p_definition": "(P - phi)/sqrt(pi)",
            "p": _text(p),
            "expected_p": EXPECTED_P_BERNOULLI,
            "p_agreement_residual": _text(p_agreement_residual),
            "p_agrees_to_1e-15": p_agrees,
        },
        "source_packet": {
            "regulator_r0": {
                "Q": ["0", "1"],
                "reference_measure": ["0.5", "0.5"],
                "features": {"X": [0, 1]},
                "moments": {"E[X]": checks_block["p"]},
                "mu_star": checks_block["mu0_star"],
                "theta": [checks_block["theta_star"]],
                "Gamma_hessian": checks_block["Gamma_hessian"],
            },
            "regulator_r1": {
                "Q": list(R1_STATES),
                "coarse_map_to_r0": dict(COARSE_MAP),
                "reference_measure": ["0.25", "0.25", "0.25", "0.25"],
                "features": {"X": list(R1_FEATURE_X), "Y": list(R1_FEATURE_Y)},
                "moments": {"E[X]": checks_block["p"], "E[Y]": "0.5"},
                "mu_star": checks_block["mu1_star"],
                "theta": [checks_block["theta_star"], "0.0"],
            },
        },
        "mechanism_checks": checks_block,
        "physical_standard_model_moment_packet": {
            "status": "ABSENT",
            "gate": "fail_closed",
            "open_proofs": list(OPEN_PROOFS),
            "promotion_rule": (
                "promotion_allowed is false because the physical moment vector "
                "and its completeness proofs are open."
            ),
        },
        "theorem_registry": {
            "SOURCE_ACTION_RIGIDITY_THEOREM": {
                "statement": (
                    "On a finite quotient with full-support reference measure "
                    "and interior feasible moments, the relative-entropy "
                    "minimizer is a unique full-support Gibbs law, and the "
                    "source action theta_star . F is unique modulo additive "
                    "constants, exact feature redundancies, and BRST-exact "
                    "terms."
                ),
            },
            "PINSKER_SELECTOR_GAP": {
                "statement": (
                    "Every feasible law nu obeys D(nu||m) - D(mu_star||m) = "
                    "D(nu||mu_star) >= ||nu - mu_star||_1^2 / 2, so the "
                    "selector gap outside an L1 ball of radius delta is at "
                    "least delta^2/2.  Uniqueness is quantitative."
                ),
            },
            "REFINEMENT_PUSHFORWARD_THEOREM": {
                "statement": (
                    "For a surjective coarse map with compatible features, "
                    "compatible reference pushforward, and entropy-maximal "
                    "fiber conditionals, the fine MaxEnt law pushes forward "
                    "exactly to the coarse MaxEnt law."
                ),
            },
            "BRST_CLASS_UNIQUENESS": {
                "statement": (
                    "On a gauge-fixed lift with nilpotent BRST differential "
                    "and BRST-invariant measure, expectations of BRST-closed "
                    "observables are independent of the gauge-fixing fermion, "
                    "so the source action is unique as the BRST equivalence "
                    "class [S_star + s Psi]."
                ),
            },
        },
        "precision": {"mpmath_dps": WORKING_DPS},
        "checks_pass": bool(checks_pass),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=PIXEL_CERTIFICATE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact(load_certificate(args.certificate))
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
