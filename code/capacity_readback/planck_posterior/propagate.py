#!/usr/bin/env python3
"""Joint Planck posterior propagated through the Lambda->N map (GAP-A5 pending item).

Chain:  Lambda = 3*Omega_L*H0^2/c^2 ;  N_Lambda = 3*pi/(Lambda*l_P^2) = pi*c^5/(Omega_L*H0^2*hbar*G).
Sufficient statistic: Omega_L*h^2 with H0 = 100h km/s/Mpc, so
    ln N_Lambda = ln(pi*c^5/(hbar*G)) - ln(Omega_L*h^2) - 2*ln(100 km/s/Mpc).
The bridge side N_CRC is certified at relative width 1.6e-25 (delta function here).

Uncertainty model: Gaussian in ln(Omega_L h^2). Three independent routes to the
relative sigma of Omega_L h^2, each using published Planck 2018 marginals plus one
degeneracy assumption, are computed and compared; the spread across routes and
across the correlation sweep is the model uncertainty of the propagation.

Everything is conditional on CP-1..CP-3 and the readback map F (CL-7): this is a
propagation of the comparison's measurement-side uncertainty, not a promotion of
the row to an unconditional test.
"""
import json
import math

# --- constants (SI; c, hbar exact by SI definition; G CODATA 2022) ---
c = 299_792_458.0
hbar = 1.054_571_817e-34
G = 6.674_30e-11
sigma_G_rel = 2.2e-5  # CODATA 2022 relative; negligible but carried
MPC = 3.085_677_581_491_367e22  # m
H100 = 100.0e3 / MPC  # 100 km/s/Mpc in 1/s

# --- bridge side (certified, G2-GAP-1 conditional fixed point) ---
N_CRC = 3.5321315434189358e122
N_CRC_rel_width = 1.58e-25

# --- Planck 2018 VI base-LCDM marginals (Table 2) ---
COMBOS = {
    "TTTEEE_lowE_lensing": {
        "H0": (67.36, 0.54),
        "OmegaL": (0.6847, 0.0073),
        "Omh2": (0.14240, 0.00087),
        "note": "column consumed implicitly by the corpus (N_Lambda = 3.313e122 back-solves to these centrals)",
    },
    "TTTEEE_lowE_lensing_BAO": {
        "H0": (67.66, 0.42),
        "OmegaL": (0.6889, 0.0056),
        "Omh2": (0.14100, 0.00087),
        "note": "sensitivity case; not the combination the ledger consumed",
    },
}


def n_lambda(omega_l, h0_kms):
    lam = 3.0 * omega_l * (h0_kms * 1e3 / MPC) ** 2 / c**2
    lp2 = hbar * G / c**3
    return 3.0 * math.pi / (lam * lp2)


def rel_sigma_route1(combo, rho):
    """Route 1: Omega_L*h^2 = Omega_L * h^2 from the (Omega_L, H0) pair, correlation rho."""
    (oL, s_oL) = combo["OmegaL"]
    (h0, s_h0) = combo["H0"]
    a = s_oL / oL
    b = 2.0 * s_h0 / h0
    v = a * a + b * b + 2.0 * rho * a * b
    return math.sqrt(max(v, 0.0))


def rel_sigma_route2(combo, rho):
    """Route 2: Omega_L*h^2 = h^2 - Om*h^2 from (H0, Om h^2), correlation rho (expected negative)."""
    (h0, s_h0) = combo["H0"]
    (omh2, s_omh2) = combo["Omh2"]
    h = h0 / 100.0
    olh2 = h * h - omh2
    s_h2 = 2.0 * h * (s_h0 / 100.0)
    v = s_h2**2 + s_omh2**2 - 2.0 * rho * s_h2 * s_omh2
    return math.sqrt(max(v, 0.0)) / olh2


def rel_sigma_route3(combo):
    """Route 3: use the tight Planck degeneracy Om*h^3 ~= const to eliminate Om h^2.

    Omega_L h^2 = h^2 - (Om h^3)/h ; treat Om h^3 as fixed, propagate h only.
    d(OLh2)/dh = 2h + (Om h^3)/h^2 = 2h + (Om h^2)/h.
    """
    (h0, s_h0) = combo["H0"]
    (omh2, _) = combo["Omh2"]
    h = h0 / 100.0
    olh2 = h * h - omh2
    d = 2.0 * h + omh2 / h
    return d * (s_h0 / 100.0) / olh2


results = {"inputs": {"N_CRC": N_CRC, "N_CRC_rel_width": N_CRC_rel_width,
                      "G_rel_sigma": sigma_G_rel, "combos": COMBOS},
           "combos": {}}

for name, combo in COMBOS.items():
    oL, h0 = combo["OmegaL"][0], combo["H0"][0]
    NL = n_lambda(oL, h0)
    gap_ln = math.log(N_CRC / NL)
    rows = []
    # correlation sweep for routes 1 and 2
    for rho in (0.90, 0.95, 0.99, 1.00):
        rows.append(("route1_OmegaL_H0", rho, rel_sigma_route1(combo, rho)))
    for rho in (-0.90, -0.95, -0.99, -1.00):
        rows.append(("route2_H0_Omh2", rho, rel_sigma_route2(combo, rho)))
    rows.append(("route3_Omh3_degeneracy", None, rel_sigma_route3(combo)))
    zrows = []
    for route, rho, rs in rows:
        # ln-space sigma: sigma(ln N) = rel sigma(Omega_L h^2) (+ G in quadrature)
        s_ln = math.sqrt(rs * rs + sigma_G_rel * sigma_G_rel)
        zrows.append({"route": route, "rho": rho, "rel_sigma_OLh2": rs,
                      "sigma_lnN": s_ln, "z": gap_ln / s_ln})
    zs = [r["z"] for r in zrows]
    results["combos"][name] = {
        "note": combo["note"],
        "N_Lambda": NL,
        "gap_ratio": N_CRC / NL,
        "gap_ln": gap_ln,
        "rows": zrows,
        "z_min": min(zs), "z_max": max(zs),
    }

print(f"{'combo':34s} {'N_Lambda':>12s} {'gap':>7s} {'z range':>18s}")
for name, r in results["combos"].items():
    print(f"{name:34s} {r['N_Lambda']:.4e} {r['gap_ratio']-1:6.2%} "
          f"{r['z_min']:7.2f} .. {r['z_max']:5.2f}")
    for row in r["rows"]:
        rho = "  --" if row["rho"] is None else f"{row['rho']:+.2f}"
        print(f"    {row['route']:26s} rho={rho}  rel_sigma={row['rel_sigma_OLh2']:.4%}"
              f"  z={row['z']:.3f}")

with open(__file__.replace("propagate.py", "planck_lambda_to_N_propagation.json"), "w") as f:
    json.dump(results, f, indent=2)
print("\nartifact written: planck_lambda_to_N_propagation.json")
