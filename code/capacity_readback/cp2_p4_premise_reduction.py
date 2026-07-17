#!/usr/bin/env python3
"""Machine checks for the CP-2 (GAP-A3) and P4 (GAP-A8) premise reductions.

Certifies, with mpmath interval arithmetic, the two reduction theorems recorded
in CP2_INVERSION_FORM_REDUCTION_2026-07-17.md and
P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md.

Part A (CP-2, GAP-A3). Every load-mediated readback candidate is
F(N) = pi*exp(phi(X_read(N))) for a strictly increasing load gauge phi
(gauge parameterization lemma; exact, since the load chart X(N) = log(N/pi)
is a proven bijection). Checks:

  A1. Read consistency (RC: the reconstructed capacity carries the load that
      was read, X(G(x)) = x) holds for the inversion gauge phi = id, with the
      residual certified below 1e-40 on the probe grid.
  A2. RC is independent of the remaining class axioms: deformation families
      (power, affine-in-capacity, shift, balance-pinned sine) either violate a
      named axiom (seed) or satisfy every remaining axiom while violating RC
      with a certified residual. The surviving freedom of the class without RC
      is exactly the gauge phi.
  A3. P4 forcing of the fixed-point location: under the CP-4 carrier
      (count density stationary exactly at load balance), the coupled fixed
      point of every gauge with phi(x_EW) = x_EW sits in the bridge box, and
      every gauge with phi(x_EW) != x_EW has a certified off-balance fixed
      point. The gauge freedom is harmless for the location, load-bearing for
      the contraction rate.

Part B (P4, GAP-A8). Checks:

  B1. Aggregate reading of membership clause 4: support collapse. The coupled
      readback residual is certified nonzero off balance and encloses zero at
      balance, so the clause-4-exact support is the fixed-point box and the
      MAR argmax over the support is the fixed point (degenerate coherence).
  B2. Per-configuration reading, positive control: for a count density
      satisfying CP-4 the derivative-sign certificate holds and the argmax box
      agrees with the fixed-point box.
  B3. Negative control, location clause dropped: stationarity off balance
      gives a certified divergence of the two representations.
  B4. Negative control, sign-carrier clause dropped: a remote maximum above
      the balance value gives a certified divergence while the balance point
      stays stationary to below 1e-80.
  B5. Recorded-lattice control: the interior stationary points of the
      2026-07-14 declared-structure counts sit at c/(1-rho) <= 18 nats,
      certified at least 121 decimal orders below the bridge capacity
      (matches p4_coherent_rows = 0 of the recorded run).

Blindness status, recorded before any comparison: the evaluation cone
contains the bridge exponent 6*pi/(P*alpha_U) by construction. These are
theorem-level checks with the same standing as F_candidate_coupled.py; no
blind CL-7 landing is claimed. Moves no ledger row.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_CANON = _HERE
if not (_CANON / "toy_readback.py").exists():
    _CANON = (
        _HERE.parents[2]
        / "reverse-engineering-reality"
        / "code"
        / "capacity_readback"
    )
    sys.path.insert(0, str(_CANON))

from mpmath import iv, mp, mpf

from toy_readback import _contained, _endpoints, _interval_json
from F_candidate_capL import P_HI, P_LO
from F_candidate_coupled import ALPHA_U_HI, ALPHA_U_LO

ARTIFACT_NAME = "oph_cp2_p4_premise_reduction_check_2026-07-17"
PRECISION = 60
LAMBDA = "0.5"


def _thin(x) -> object:
    return iv.mpf(mp.nstr(mpf(x), mp.dps))


def _abs_lo(box) -> mpf:
    """Certified lower bound of |box|."""
    lo, hi = _endpoints(box)
    if lo > 0:
        return lo
    if hi < 0:
        return -hi
    return mpf(0)


def _abs_hi(box) -> mpf:
    lo, hi = _endpoints(box)
    return max(abs(lo), abs(hi))


def _banach_box(C, Cp, seed, half_width: str, iterations: int = 400) -> dict:
    """Iterate y <- C(y) from seed, then certify the box y* +- half_width."""
    point = mpf(seed)
    for _ in range(iterations):
        point = mpf(C(_thin(point)).mid)
    h = mpf(half_width)
    box = iv.mpf([mp.nstr(point - h, mp.dps), mp.nstr(point + h, mp.dps)])
    m = _thin(mpf(box.mid))
    image = C(m) + Cp(box) * (box - m)
    ok = _contained(image, box)
    dlo, dhi = _endpoints(Cp(box))
    L = max(abs(dlo), abs(dhi))
    return {
        "box": _interval_json(box),
        "image": _interval_json(image),
        "self_map_pass": bool(ok),
        "lipschitz_bound_L": mp.nstr(L, 30),
        "banach_pass": bool(ok and L < 1),
        "enclosure": _interval_json(image),
    }


def build() -> dict:
    iv.dps = PRECISION
    mp.dps = PRECISION

    P = iv.mpf([P_LO, P_HI])
    alpha_u = iv.mpf([ALPHA_U_LO, ALPHA_U_HI])
    lam = iv.mpf(LAMBDA)
    one = iv.mpf(1)
    x_ew = 6 * iv.pi / (P * alpha_u)
    x_lo, x_hi = _endpoints(x_ew)
    n_crc = iv.pi * iv.exp(x_ew)

    def x_read(y):
        return (one - lam) * y + lam * x_ew

    # ------------------------------------------------------------------
    # A1. Read consistency holds exactly for the inversion gauge.
    # ------------------------------------------------------------------
    rc_probes = ["0", "0.5", "1", "5", "50", mp.nstr(mpf(x_ew.mid), 40)]
    rc_residuals = []
    rc_pass = True
    for probe in rc_probes:
        x0 = iv.mpf(probe)
        residual = iv.log(iv.pi * iv.exp(x0) / iv.pi) - x0
        bound = _abs_hi(residual)
        rc_residuals.append({"x": probe, "abs_residual_hi": mp.nstr(bound, 8)})
        if bound > mpf("1e-40"):
            rc_pass = False
    a1 = {
        "statement": "RC residual X(G(x)) - x for G(x) = pi*exp(x), certified below 1e-40 on the probe grid",
        "probes": rc_residuals,
        "pass": bool(rc_pass),
    }

    # ------------------------------------------------------------------
    # A2/A3. Gauge deformation families.
    # Class axioms checked per family: seed (phi(0) = 0), monotone gauge,
    # coupled Banach contraction, RC, coupled fixed-point location vs balance.
    # ------------------------------------------------------------------
    families = []

    def _record_family(
        name,
        gauge,
        seed_residual,
        monotone_note,
        coupled_cert,
        separation_box,
        rc_probe,
        rc_residual,
        classification,
    ):
        sep_lo = _abs_lo(separation_box)
        rc_bound = _abs_lo(rc_residual)
        families.append(
            {
                "family": name,
                "gauge": gauge,
                "seed_residual_abs_hi": mp.nstr(_abs_hi(seed_residual), 8),
                "seed_pass": bool(_abs_hi(seed_residual) < mpf("1e-40")),
                "monotone": monotone_note,
                "coupled_banach_pass": bool(coupled_cert["banach_pass"]),
                "coupled_fixed_point_load": coupled_cert["box"],
                "rc_probe_x": rc_probe,
                "rc_residual_abs_lo": mp.nstr(rc_bound, 8),
                "rc_pass": bool(rc_bound == 0),
                "separation_from_balance_lo": mp.nstr(sep_lo, 8),
                "p4_cp4_coherent": bool(sep_lo < mpf("1e-8")),
                "classification": classification,
            }
        )

    # Inversion member phi = id (reference row).
    def C_id(y):
        return x_read(y)

    def Cp_id(y):
        return (one - lam) + 0 * y

    cert_id = _banach_box(C_id, Cp_id, x_ew.mid, "1e-10")
    box_id = iv.mpf([cert_id["box"]["lo"], cert_id["box"]["hi"]])
    _record_family(
        "inversion (phi = id)",
        "phi(x) = x",
        iv.mpf(0),
        "phi' = 1 exactly",
        cert_id,
        box_id - x_ew,
        "1",
        iv.mpf(0) * 0,
        "the CP-2 form; RC holds (A1); coupled fixed point in the bridge box",
    )
    # Correction: RC residual for the identity row is zero; overwrite verdicts.
    families[-1]["rc_residual_abs_lo"] = "0.0"
    families[-1]["rc_pass"] = True

    # Power gauges phi_s(x) = s*x: the CAP-P degradation branches.
    s_values = [
        ("s = 1/2", iv.mpf("0.5")),
        ("s = exp(-P/24) (poisson reserve)", iv.exp(-P / 24)),
        ("s = 1 - P/24 (presence reserve)", one - P / 24),
        ("s = 0.99", iv.mpf("0.99")),
    ]
    for label, s in s_values:
        def C_s(y, s=s):
            return s * x_read(y)

        def Cp_s(y, s=s):
            return s * (one - lam) + 0 * y

        closed_form = s * lam * x_ew / (one - s * (one - lam))
        cert = _banach_box(C_s, Cp_s, closed_form.mid, "1e-10")
        box = iv.mpf([cert["box"]["lo"], cert["box"]["hi"]])
        # decoupled member: C0(y) = s*y fixes 0, so N = pi exactly.
        decoupled = "fixed point pi exact (s*y = y forces y = 0)"
        rc_res = s * iv.mpf(1) - iv.mpf(1)
        _record_family(
            f"power gauge, {label}",
            "phi(x) = s*x",
            iv.mpf(0),
            "phi' = s > 0",
            cert,
            box - x_ew,
            "1",
            rc_res,
            "live member of the load-gauge freedom: satisfies P1-P3/P5 and the "
            "seed, violates RC; coupled fixed point off balance, so P4 plus "
            "CP-4 excludes it; " + decoupled,
        )

    # Affine-in-capacity gauges F0(N) = a*N + pi*(1-a):
    # phi_a(x) = log(a*e^x + 1 - a); seed holds by construction.
    for a_str in ["0.3", "0.7"]:
        a = iv.mpf(a_str)

        def phi_a(x, a=a):
            return iv.log(a * iv.exp(x) + one - a)

        def C_a(y, a=a):
            return phi_a(x_read(y))

        def Cp_a(y, a=a):
            u = x_read(y)
            e = a * iv.exp(u)
            return (one - lam) * e / (e + one - a)

        cert = _banach_box(C_a, Cp_a, mpf(x_ew.mid) + 2 * mp.log(mpf(a_str)), "1e-10")
        box = iv.mpf([cert["box"]["lo"], cert["box"]["hi"]])
        rc_res = phi_a(iv.mpf(1)) - iv.mpf(1)
        _record_family(
            f"affine capacity gauge, a = {a_str}",
            "phi(x) = log(a*e^x + 1 - a), from F0(N) = a*N + pi*(1-a)",
            phi_a(iv.mpf(0)),
            "phi' = a*e^x/(a*e^x + 1 - a) > 0",
            cert,
            box - x_ew,
            "1",
            rc_res,
            "live member of the load-gauge freedom (seed holds, monotone, "
            "contraction); violates RC; coupled fixed point sits near "
            "x_EW + 2*log(a), certified off balance; P4 plus CP-4 excludes it",
        )

    # Shift gauge phi(x) = x + 1/10: violates the seed normalization.
    shift = iv.mpf("0.1")

    def C_shift(y):
        return x_read(y) + shift

    def Cp_shift(y):
        return (one - lam) + 0 * y

    cert = _banach_box(C_shift, Cp_shift, mpf(x_ew.mid) + mpf("0.2"), "1e-10")
    box = iv.mpf([cert["box"]["lo"], cert["box"]["hi"]])
    _record_family(
        "shift gauge",
        "phi(x) = x + 1/10",
        shift,
        "phi' = 1",
        cert,
        box - x_ew,
        "1",
        shift,
        "violates the seed normalization phi(0) = 0 (named principle: the "
        "certified decoupled fixed point pi); decoupled map N*e^(1/10) has no "
        "fixed point; coupled fixed point also off balance",
    )

    # Balance-pinned sine gauge: phi(x) = x + eps*(x_EW/pi)*sin(pi*x/x_EW).
    # phi(0) = 0 and phi(x_EW) = x_EW exactly; phi != id everywhere else.
    eps = iv.mpf("0.1")
    amp = eps * x_ew / iv.pi

    def phi_pin(x):
        return x + amp * iv.sin(iv.pi * x / x_ew)

    def C_pin(y):
        return phi_pin(x_read(y))

    def Cp_pin(y):
        u = x_read(y)
        return (one - lam) * (one + eps * iv.cos(iv.pi * u / x_ew))

    cert = _banach_box(C_pin, Cp_pin, x_ew.mid, "1e-8")
    box = iv.mpf([cert["box"]["lo"], cert["box"]["hi"]])
    rc_probe = x_ew / 2
    rc_res = phi_pin(rc_probe) - rc_probe
    _record_family(
        "balance-pinned sine gauge",
        "phi(x) = x + (1/10)*(x_EW/pi)*sin(pi*x/x_EW)",
        phi_pin(iv.mpf(0)),
        "phi' = 1 + (1/10)*cos(pi*x/x_EW) in [9/10, 11/10] > 0",
        cert,
        box - x_ew,
        "x_EW/2",
        rc_res,
        "harmlessness witness: phi != id (RC violated at x_EW/2 by about "
        "8.9 nats) yet phi(x_EW) = x_EW, and the coupled fixed point is "
        "certified inside the bridge box; the gauge freedom cannot move the "
        "fixed-point location once P4 plus CP-4 pin it at balance",
    )

    a2_a3 = {
        "statement": "gauge parameterization: every load-mediated candidate is "
        "pi*exp(phi(X_read)); the class axioms pin phi only at 0 (seed) and, "
        "under P4 plus CP-4, at x_EW; RC (phi = id) is independent of the rest",
        "families": families,
    }

    # ------------------------------------------------------------------
    # B1. Aggregate reading: support collapse.
    # ------------------------------------------------------------------
    collapse_probes = []
    collapse_pass = True
    for delta in ["-30", "-10", "-1", "-0.1", "0.1", "1", "10", "30"]:
        y = x_ew + iv.mpf(delta)
        residual = C_id(y) - y  # = lam*(x_EW - y)
        lo_bound = _abs_lo(residual)
        expected = mpf(LAMBDA) * abs(mpf(delta))
        ok = bool(lo_bound > mpf("0.9") * expected)
        collapse_pass = collapse_pass and ok
        collapse_probes.append(
            {
                "load_offset": delta,
                "abs_readback_residual_lo": mp.nstr(lo_bound, 8),
                "nonzero_certified": ok,
            }
        )
    at_balance = C_id(x_ew) - x_ew
    balance_width = _abs_hi(at_balance)
    b1 = {
        "statement": "aggregate reading of membership clause 4: the readback "
        "residual is certified nonzero at every off-balance probe and "
        "encloses zero at balance, so the clause-4-exact support is the "
        "fixed-point box and the MAR argmax over the support equals the fixed "
        "point (degenerate coherence, the count representation carries no "
        "independent information)",
        "off_balance_probes": collapse_probes,
        "balance_residual_abs_hi": mp.nstr(balance_width, 8),
        "pass": bool(collapse_pass and balance_width < mpf("1e-20")),
    }

    # ------------------------------------------------------------------
    # B2. Per-configuration positive control: CP-4 carrier
    #     l(y) = -(y - x_EW)^2/2, H proportional to -(y - x_EW).
    # ------------------------------------------------------------------
    sign_probes = []
    sign_pass = True
    for delta in ["-30", "-5", "-0.5", "0.5", "5", "30"]:
        y = x_ew + iv.mpf(delta)
        H = -(y - x_ew)  # sign of dl/dN on the load coordinate
        lo, hi = _endpoints(H)
        want_positive = mpf(delta) < 0
        ok = bool(lo > 0) if want_positive else bool(hi < 0)
        sign_pass = sign_pass and ok
        sign_probes.append(
            {"load_offset": delta, "H_sign_certified": ok}
        )
    argmax_box = x_ew  # unique zero of -(y - x_EW)
    agreement = argmax_box - box_id  # box_id: coupled fixed point load box
    b2 = {
        "statement": "CP-4 carrier l = -(X - x_EW)^2/2: derivative-sign "
        "certificate holds (H > 0 below balance, H < 0 above, unique zero at "
        "balance) and the argmax box agrees with the coupled fixed-point box",
        "derivative_sign_probes": sign_probes,
        "argmax_load": _interval_json(argmax_box),
        "fixed_point_load": _interval_json(box_id),
        "agreement_abs_hi": mp.nstr(_abs_hi(agreement), 8),
        "pass": bool(sign_pass and _abs_hi(agreement) < mpf("1e-8")),
    }

    # ------------------------------------------------------------------
    # B3. Negative control: stationarity off balance.
    # ------------------------------------------------------------------
    b3_rows = []
    for a_label, a_box in [
        ("x_EW - 5", x_ew - 5),
        ("0.85 * x_EW", iv.mpf("0.85") * x_ew),
    ]:
        separation = _abs_lo(a_box - box_id)
        b3_rows.append(
            {
                "stationary_point": a_label,
                "separation_from_fixed_point_lo": mp.nstr(separation, 8),
                "diverges": bool(separation > mpf("1")),
            }
        )
    b3 = {
        "statement": "location clause of CP-4 dropped: the argmax moves with "
        "the stationary point while the readback fixed point stays at "
        "balance; the two representations diverge",
        "rows": b3_rows,
        "pass": bool(all(r["diverges"] for r in b3_rows)),
    }

    # ------------------------------------------------------------------
    # B4. Negative control: sign carrier dropped (remote maximum).
    #     l(y) = -(y - x_EW)^2/2 + A*exp(-(y - y0)^2/2), A = 300, y0 = x_EW - 20.
    # ------------------------------------------------------------------
    # Evaluated in the load-offset coordinate delta = y - x_EW (balance at
    # delta = 0 exactly; avoids the spurious interval width of x_EW - x_EW).
    A = iv.mpf(300)

    def ell_bump_offset(delta):
        return -(delta**2) / 2 + A * iv.exp(-((delta + 20) ** 2) / 2)

    value_gap = ell_bump_offset(iv.mpf(-20)) - ell_bump_offset(iv.mpf(0))
    zero = iv.mpf(0)
    H_at_balance = -zero - A * (zero + 20) * iv.exp(-((zero + 20) ** 2) / 2)
    b4 = {
        "statement": "sign-carrier clause of CP-4 dropped: a remote maximum "
        "20 nats below balance tops the balance value by about 100 while the "
        "balance point stays stationary to below 1e-80; the MAR argmax leaves "
        "the fixed point",
        "value_gap_lo": mp.nstr(_abs_lo(value_gap), 8),
        "H_at_balance_abs_hi": mp.nstr(_abs_hi(H_at_balance), 8),
        "pass": bool(
            _abs_lo(value_gap) > mpf(90) and _abs_hi(H_at_balance) < mpf("1e-80")
        ),
    }

    # ------------------------------------------------------------------
    # B5. Recorded-lattice control (2026-07-14 run).
    # ------------------------------------------------------------------
    reserve = iv.log(one - P / 24)
    rho_axis = [
        ("R1 poisson per-cell", one - iv.mpf(1) / 6),
        ("R2 presence per-cell", one + (4 / P) * reserve),
        ("R4 poisson per-edge", iv.mpf("0.5")),
        ("R5 presence per-edge", one + (12 / P) * reserve),
    ]
    lattice_rows = []
    min_log10_sep = None
    for label, rho in rho_axis:
        for c in (1, 2, 3):
            n_star = c / (one - rho)
            log10_sep = (iv.log(n_crc) - iv.log(n_star)) / iv.log(iv.mpf(10))
            sep_lo = mpf(_endpoints(log10_sep)[0])
            if min_log10_sep is None or sep_lo < min_log10_sep:
                min_log10_sep = sep_lo
            lattice_rows.append(
                {
                    "branch": label,
                    "c": c,
                    "N_star_nats": _interval_json(n_star),
                    "log10_separation_lo": mp.nstr(sep_lo, 8),
                }
            )
    b5 = {
        "statement": "every interior stationary point of the recorded "
        "declared-structure counts sits at c/(1-rho) <= 18 nats, certified at "
        "least 121 decimal orders below the bridge capacity; restates "
        "p4_coherent_rows = 0 of the 2026-07-14 run",
        "rows": lattice_rows,
        "min_log10_separation": mp.nstr(min_log10_sep, 8),
        "pass": bool(min_log10_sep > mpf(121)),
    }

    return {
        "artifact": ARTIFACT_NAME,
        "status": "reductions_certified_no_premise_discharged",
        "derivations": [
            "CP2_INVERSION_FORM_REDUCTION_2026-07-17.md",
            "P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md",
        ],
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": PRECISION,
            "rounding": "mpmath_interval_outward",
            "promotion_backend_required": "arb_or_mpfi_directed_outward",
        },
        "inputs": {
            "P_certified_enclosure": {"lo": P_LO, "hi": P_HI},
            "alpha_U_certified_enclosure": {"lo": ALPHA_U_LO, "hi": ALPHA_U_HI},
            "source": "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json",
            "lambda": LAMBDA,
        },
        "x_ew_enclosure": _interval_json(x_ew),
        "part_A_cp2": {"A1_read_consistency": a1, "A2_A3_gauge_class": a2_a3},
        "part_B_p4": {
            "B1_aggregate_support_collapse": b1,
            "B2_positive_control_cp4": b2,
            "B3_negative_control_location": b3,
            "B4_negative_control_sign_carrier": b4,
            "B5_recorded_lattice_control": b5,
        },
        "blindness": {
            "cone_contains_cl3_bridge_expression": True,
            "v08_note": "theorem-level checks; not a blind CL-7 candidate; "
            "same standing as F_candidate_coupled.py",
        },
        "gap_status": {
            "GAP-A3_cp2": "open; reduced to load mediation (LM) plus read "
            "consistency (RC); surviving freedom named: the load gauge phi, "
            "pinned at 0 by the seed and at x_EW by P4 plus CP-4",
            "GAP-A8_p4": "open; per-configuration reading reduced to CP-4 "
            "(stationarity at balance with the derivative-sign carrier); "
            "aggregate reading closes it degenerately by support collapse",
        },
        "moves_cl7": False,
        "cl7_status": "open",
    }


def main() -> int:
    certificate = build()
    out = _CANON / "runtime" / "cp2_p4_premise_reduction_check.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "A1_rc_identity": certificate["part_A_cp2"]["A1_read_consistency"]["pass"],
        "A2_families": len(certificate["part_A_cp2"]["A2_A3_gauge_class"]["families"]),
        "B1_support_collapse": certificate["part_B_p4"]["B1_aggregate_support_collapse"]["pass"],
        "B2_positive": certificate["part_B_p4"]["B2_positive_control_cp4"]["pass"],
        "B3_negative_location": certificate["part_B_p4"]["B3_negative_control_location"]["pass"],
        "B4_negative_sign": certificate["part_B_p4"]["B4_negative_control_sign_carrier"]["pass"],
        "B5_lattice": certificate["part_B_p4"]["B5_recorded_lattice_control"]["pass"],
        "status": certificate["status"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
