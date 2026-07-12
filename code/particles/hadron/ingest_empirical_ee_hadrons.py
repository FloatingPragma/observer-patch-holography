#!/usr/bin/env python3
"""Build the empirical e+e- hadronic spectral measure payload.

Populates code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json
per the schema in this directory and the source registry policy. Row class is
oph_plus_empirical_hadron_closure: external cross-section data is used, the
payload is usable for public final values on the empirical closure surface,
and it is never promotable as an OPH source theorem.

Compilation. The R(s) input is a documented piecewise compilation built from
PDG 2025 resonance parameters and perturbative QCD, in the style of the
classic dispersion evaluations:

  * two-pion channel 0.32-1.05 GeV via vector-meson dominance with the
    s-dependent p-wave rho(770) width (PDG mass/width), R = (beta^3/4)|F_pi|^2;
  * omega(782) and phi(1020) as narrow resonances via their PDG electronic
    widths and hadronic branching fractions;
  * inclusive continuum regions with flat or linear R and declared
    normalization budgets (1.05-2.0 rising 1.4->2.2 at 15 percent; 2.0-3.7 at
    R=2.15, 6 percent; 3.7-5.0 charm threshold at R=3.6, 15 percent; 5.0-10.5
    at R=3.6, 6 percent; 10.5-11.2 at R=3.66, 8 percent);
  * J/psi, psi(2S), Upsilon(1S,2S,3S) as narrow resonances;
  * a five-flavor perturbative QCD tail above 11.2 GeV with
    R = (11/3)(1 + a + 1.409 a^2 - 12.767 a^3), a = alpha_s(s)/pi at one loop,
    2 percent budget, integrated numerically to 200 GeV and analytically
    beyond.

Kernel. The payload integral is the subtracted vacuum-polarization dispersion
integral in the on-shell convention,

  Delta_alpha_had^(5)(M_Z^2)
    = (alpha M_Z^2 / 3 pi) P.V. Int_{s_th}^inf ds R(s) / (s (M_Z^2 - s)),

with the principal value across s = M_Z^2 handled by subtracting R(M_Z^2)
inside the perturbative tail. Narrow resonances contribute

  delta = (3 M Gamma_ee B_had / alpha) * M_Z^2 / (M^2 (M_Z^2 - M^2)).

Uncertainty. Region normalization budgets combine in quadrature across
regions; electronic-width errors enter the resonance terms; a global
compilation-coarseness term of 2 percent is added in quadrature. The known
data-driven compilations (KNT19: 0.02761(11)) serve as an external
cross-check row, not as an input.

Run:
    python3 code/particles/hadron/ingest_empirical_ee_hadrons.py
writes code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json.
"""

from __future__ import annotations

import json
import math
import pathlib
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
OUT_PATH = HERE.parent / "runs" / "hadron" / "empirical_ee_hadronic_spectral_measure.json"

ALPHA0 = 1.0 / 137.035999177  # CODATA 2022, used only to normalize measured widths
M_Z = 91.1876  # GeV, PDG
M_Z2 = M_Z * M_Z
M_PI = 0.13957039  # GeV, charged pion

# rho(770) vector-meson-dominance parameters (PDG 2025)
M_RHO, G_RHO = 0.77526, 0.1474

# narrow resonances: (name, mass GeV, Gamma_ee GeV, B_had, rel. width error)
NARROW = [
    ("omega(782)", 0.78266, 0.625e-6, 0.916, 0.02),
    ("phi(1020)", 1.019461, 1.266e-6, 0.99, 0.02),
    ("J/psi(1S)", 3.096900, 5.53e-6, 0.877, 0.02),
    ("psi(2S)", 3.686097, 2.33e-6, 0.982, 0.03),
    ("Upsilon(1S)", 9.46040, 1.340e-6, 0.94, 0.02),
    ("Upsilon(2S)", 10.02326, 0.612e-6, 0.94, 0.03),
    ("Upsilon(3S)", 10.35520, 0.443e-6, 0.94, 0.03),
]

# continuum regions: (label, sqrt_s_lo, sqrt_s_hi, R_lo, R_hi, norm_budget)
REGIONS = [
    ("continuum_1.05_2.0", 1.05, 2.00, 1.4, 2.2, 0.15),
    ("uds_2.0_3.7", 2.00, 3.70, 2.15, 2.15, 0.06),
    ("charm_threshold_3.7_5.0", 3.70, 5.00, 3.6, 3.6, 0.15),
    ("udsc_5.0_10.5", 5.00, 10.50, 3.6, 3.6, 0.06),
    ("threshold_b_10.5_11.2", 10.50, 11.20, 3.66, 3.66, 0.08),
]

TWO_PI_BUDGET = 0.03
PQCD_BUDGET = 0.02
COARSENESS_BUDGET = 0.02
PQCD_START = 11.20  # GeV
PQCD_NUMERIC_END = 200.0  # GeV
LAMBDA_QCD5 = 0.210  # GeV, one-loop five-flavor scale


def kernel(s: float) -> float:
    """Dispersion kernel M_Z^2 / (s (M_Z^2 - s)) without the alpha/3pi factor."""
    return M_Z2 / (s * (M_Z2 - s))


def rho_pion_form_factor_sq(s: float) -> float:
    """|F_pi(s)|^2 from vector-meson dominance with the p-wave running width."""
    if s <= 4.0 * M_PI * M_PI:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * M_PI * M_PI / s)
    beta_rho = math.sqrt(1.0 - 4.0 * M_PI * M_PI / (M_RHO * M_RHO))
    gamma_s = G_RHO * (s / (M_RHO * M_RHO)) * (beta / beta_rho) ** 3
    denom = (M_RHO * M_RHO - s) ** 2 + (M_RHO * gamma_s) ** 2
    return M_RHO ** 4 / denom


def r_two_pion(s: float) -> float:
    if s <= 4.0 * M_PI * M_PI:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * M_PI * M_PI / s)
    return 0.25 * beta ** 3 * rho_pion_form_factor_sq(s)


def alpha_s_one_loop(s: float) -> float:
    return 4.0 * math.pi / (23.0 / 3.0 * math.log(s / (LAMBDA_QCD5 * LAMBDA_QCD5)))


def r_pqcd_nf5(s: float) -> float:
    a = alpha_s_one_loop(s) / math.pi
    return (11.0 / 3.0) * (1.0 + a + 1.409 * a * a - 12.767 * a ** 3)


def integral_two_pion() -> float:
    lo, hi, n = 0.32, 1.05, 2000
    total = 0.0
    for i in range(n):
        rs = lo + (hi - lo) * (i + 0.5) / n
        s = rs * rs
        total += r_two_pion(s) * kernel(s) * 2.0 * rs * (hi - lo) / n
    return total


def integral_region(lo: float, hi: float, r_lo: float, r_hi: float,
                    n: int = 2000) -> float:
    total = 0.0
    for i in range(n):
        rs = lo + (hi - lo) * (i + 0.5) / n
        r_val = r_lo + (r_hi - r_lo) * (rs - lo) / (hi - lo)
        s = rs * rs
        total += r_val * kernel(s) * 2.0 * rs * (hi - lo) / n
    return total


def integral_narrow(mass: float, g_ee: float, b_had: float) -> float:
    """Narrow-width contribution to Int R(s) kernel(s) ds."""
    return (9.0 * math.pi * mass * g_ee * b_had / (ALPHA0 * ALPHA0)) * kernel(mass * mass)


def integral_pqcd() -> float:
    """Perturbative tail with principal-value subtraction across s = M_Z^2."""
    r_mz = r_pqcd_nf5(M_Z2)
    s_lo, s_hi, n = PQCD_START ** 2, PQCD_NUMERIC_END ** 2, 40000
    total = 0.0
    for i in range(n):
        s = s_lo + (s_hi - s_lo) * (i + 0.5) / n
        total += (r_pqcd_nf5(s) - r_mz) * kernel(s) * (s_hi - s_lo) / n
    # analytic principal value of the constant piece over [s_lo, s_hi]:
    # P.V. Int ds M_Z^2/(s(M_Z^2-s)) = ln(s/|M_Z^2-s|) evaluated at the ends
    total += r_mz * (math.log(s_hi / abs(M_Z2 - s_hi)) - math.log(s_lo / abs(M_Z2 - s_lo)))
    # constant-R remainder above s_hi (kernel negative there):
    # Int_{s_hi}^inf ds M_Z^2/(s(M_Z^2-s)) = -ln(1/(1 - M_Z^2/s_hi)) for s_hi > M_Z^2
    a_tail = alpha_s_one_loop(s_hi) / math.pi
    r_tail = (11.0 / 3.0) * (1.0 + a_tail)
    total += -r_tail * math.log(1.0 / (1.0 - M_Z2 / s_hi))
    return total


def build_payload() -> dict:
    prefactor = ALPHA0 / (3.0 * math.pi)

    contributions = {}
    budgets = {}

    two_pi = prefactor * integral_two_pion()
    contributions["two_pion_rho"] = two_pi
    budgets["two_pion_rho"] = abs(two_pi) * TWO_PI_BUDGET

    for label, lo, hi, r_lo, r_hi, budget in REGIONS:
        val = prefactor * integral_region(lo, hi, r_lo, r_hi)
        contributions[label] = val
        budgets[label] = abs(val) * budget

    for name, mass, g_ee, b_had, err in NARROW:
        val = prefactor * integral_narrow(mass, g_ee, b_had)
        contributions[name] = val
        budgets[name] = abs(val) * err

    pqcd = prefactor * integral_pqcd()
    contributions["pqcd_tail_nf5"] = pqcd
    budgets["pqcd_tail_nf5"] = abs(pqcd) * PQCD_BUDGET

    value = sum(contributions.values())
    stat_sys = math.sqrt(sum(b * b for b in budgets.values()))
    coarse = abs(value) * COARSENESS_BUDGET
    uncertainty = math.sqrt(stat_sys * stat_sys + coarse * coarse)

    # export grid for the continuum pieces (schema requires explicit arrays)
    grid, r_vals = [], []
    for rs_int in range(32, 105, 4):
        rs = rs_int / 100.0
        grid.append(rs)
        r_vals.append(r_two_pion(rs * rs))
    for label, lo, hi, r_lo, r_hi, _b in REGIONS:
        for k in range(8):
            rs = lo + (hi - lo) * k / 8.0
            grid.append(round(rs, 6))
            r_vals.append(round(r_lo + (r_hi - r_lo) * k / 8.0, 6))
    for rs in (11.2, 15.0, 25.0, 40.0, 91.1876, 150.0, 200.0):
        grid.append(rs)
        r_vals.append(round(r_pqcd_nf5(rs * rs), 6))

    return {
        "artifact": "oph_empirical_ee_hadronic_spectral_measure",
        "format_version": 1,
        "source_compilation": {
            "id": "oph_documented_piecewise_compilation_v1",
            "name": "Documented piecewise R(s) compilation from PDG 2025 resonance "
                    "parameters plus perturbative QCD",
            "url": "https://pdg.lbl.gov/",
            "citation": "PDG 2025 resonance masses, electronic widths, and hadronic "
                        "branching fractions; five-flavor perturbative QCD tail; "
                        "compilation structure follows the standard data-driven "
                        "dispersion evaluations (KNT19 arXiv:1911.00367 as external "
                        "cross-check only)",
        },
        "data_release": {
            "release_id": "pdg2025_piecewise_v1",
            "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        },
        "row_class": "oph_plus_empirical_hadron_closure",
        "guards": {
            "source_only": False,
            "empirical_hadron_closure": True,
            "external_cross_section_data_used": True,
            "promotable_as_oph_source_theorem": False,
            "usable_for_public_final_values": True,
        },
        "energy_grid": {"variable": "sqrt_s", "unit": "GeV", "values": grid},
        "R_values": {
            "definition": "R(s) = sigma(e+e- -> hadrons) / sigma(e+e- -> mu+mu-), "
                          "bare, five flavors; narrow resonances carried "
                          "analytically, not on the grid",
            "values": r_vals,
        },
        "covariance": {
            "statistical": {
                "policy": "statistical errors subsumed into region normalization "
                          "budgets for this compilation class",
            },
            "systematic": {
                "policy": "region-wise fully correlated normalization budgets, "
                          "quadrature across regions, plus a 2 percent global "
                          "coarseness term",
                "region_budgets": {k: round(v, 8) for k, v in budgets.items()},
            },
        },
        "correction_policy": {
            "radiative_corrections": "resonance electronic widths are PDG physical "
                                     "values; no additional ISR correction applied "
                                     "at this compilation granularity",
            "vacuum_polarization_undressing": "not applied; the induced bias is "
                                              "within the declared region budgets",
            "final_state_radiation": "not separated; within declared budgets",
        },
        "region_policy": {
            "exclusive_channel_sum": "two-pion channel modeled via vector-meson "
                                     "dominance with s-dependent p-wave rho width; "
                                     "omega/phi carried as narrow resonances",
            "inclusive_region": "flat or linear R with declared normalization "
                                "budgets between 1.05 and 11.2 GeV",
            "pQCD_tail": "five-flavor massless pQCD above 11.2 GeV with one-loop "
                         "alpha_s, principal value across M_Z^2, analytic remainder "
                         "above 200 GeV",
        },
        "kernel": {
            "name": "subtracted_vacuum_polarization_dispersion",
            "formula": "Delta_alpha_had^(5)(M_Z^2) = (alpha M_Z^2 / 3 pi) P.V. "
                       "Int_{s_th}^inf ds R(s) / (s (M_Z^2 - s))",
            "target": "Delta_alpha_had_5_MZ",
        },
        "integral": {
            "value": value,
            "uncertainty": uncertainty,
            "unit": "dimensionless",
            "method": "piecewise midpoint quadrature (2000 points per region, "
                      "40000 in the perturbative tail), narrow-width resonance "
                      "terms, principal-value subtraction at M_Z^2, analytic "
                      "asymptotic remainder",
            "contributions": {k: round(v, 8) for k, v in contributions.items()},
            "external_cross_check": {
                "id": "knt19",
                "value": 0.02761,
                "uncertainty": 0.00011,
                "role": "compare_only_not_an_input",
            },
        },
    }


def main() -> int:
    payload = build_payload()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1)
        f.write("\n")
    integ = payload["integral"]
    print(f"Delta_alpha_had^(5)(M_Z) = {integ['value']:.6f} +- {integ['uncertainty']:.6f}")
    for k, v in integ["contributions"].items():
        print(f"  {k:28s} {v:+.6f}")
    print(f"external cross-check KNT19 = 0.02761 +- 0.00011 (compare only)")
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
