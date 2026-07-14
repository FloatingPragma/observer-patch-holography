#!/usr/bin/env python3
"""CAP-L candidate family for the capacity readback map F (construction run 2026-07-14).

Log-count readback: F(N) = log|Omega^sc_N| with the count assembled from declared
structure only (derivation: F_CONSTRUCTION_2026-07-14.md, Steps 2-7):

    log|Omega^sc_N| = rho*N + c*log N + d

with the (rho, c, d) branch lattice

    rho : BR-1 x BR-2 reserve semantics/attachment (5 readings)
    c,d : BR-3 readback-record effect (4) x BR-4 observer marking (3)
          x BR-5 symmetry factor (3)

180 rows. Per row: existence analysis of the stable fixed point of
N = rho*N + c*log N + d, interval Banach certificate (self-map + derivative
enclosure, mpmath.iv, outward rounding), certified fixed-point enclosure, P4
stationarity record (N_star = c/(1-rho) where interior), blindness record.

Blindness: the only non-integer input is the certified forward-closure enclosure of
P from code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json.
No measured Lambda, no SL-4 estimate, no CL-3 bridge value, no alpha_U.
"""

from __future__ import annotations

import json
from pathlib import Path

from mpmath import iv, mp, mpf

from toy_readback import _contained, _endpoints, _interval_json

ARTIFACT_NAME = "oph_capacity_readback_candidate_capL_2026-07-14"
PRECISION = 40

# Certified forward-closure enclosure of P (P_derivation interval contraction
# certificate, 2026-07-14). The only admissible P input per the construction rule:
# computed fixed point only, never a measured value.
P_LO = "1.6309720958588973769645139031692398220074613523151957129273490076"
P_HI = "1.6309720958588973769645139038446714349507911573627096569324070023"

ENCLOSURE_HALF_WIDTH = "1e-25"


def build_axes():
    P = iv.mpf([P_LO, P_HI])
    one = iv.mpf(1)
    ln_pi = iv.log(iv.pi)
    ln_4_over_P = iv.log(4 / P)
    ln_60 = iv.log(iv.mpf(60))
    reserve = iv.log(one - P / 24)  # log(1 - P/24) < 0

    rho_axis = [
        # (id, description, rho interval)
        ("R1", "poisson reserve per cell: rho = 1 - (4/P)(P/24) = 5/6", iv.mpf(5) / 6),
        ("R2", "presence reserve per cell: rho = 1 + (4/P)ln(1-P/24)", one + (4 / P) * reserve),
        ("R3", "no reserve exclusion: rho = 1", one),
        ("R4", "poisson reserve per shared edge: rho = 1 - (12/P)(P/24) = 1/2", iv.mpf(1) / 2),
        ("R5", "presence reserve per shared edge: rho = 1 + (12/P)ln(1-P/24)", one + (12 / P) * reserve),
    ]
    br3_axis = [
        ("Xp", "port-orbit multiplicity e^X = N/pi", 1, -ln_pi),
        ("X0", "readback as delta constraint on defect data", 0, iv.mpf(0)),
        ("Xm", "readback record cost X nats in bulk", -1, ln_pi),
        ("Xq", "write+check record cost 2X on oriented register", -2, 2 * ln_pi),
    ]
    br4_axis = [
        ("K0", "unmarked at-least-one observer chain", 0, iv.mpf(0)),
        ("K1", "one marked host cell: factor K = 4N/P", 1, ln_4_over_P),
        ("K2", "marked host + checkpoint cells: factor K^2", 2, 2 * ln_4_over_P),
    ]
    br5_axis = [
        ("Sm", "A5 quotient: divide by 60", -ln_60),
        ("S0", "nf already quotients: no factor", iv.mpf(0)),
        ("Sp", "sixty face-corner flag anchorings: multiply by 60", ln_60),
    ]
    return rho_axis, br3_axis, br4_axis, br5_axis


def _thin(x: mpf):
    return iv.mpf(mp.nstr(x, mp.dps))


def centered_image(F, Fp, interval):
    """Mean-value form enclosure of F(interval): F(m) + F'(I)*(I - m).

    Valid image enclosure by the mean value theorem; avoids the dependency
    blow-up of the naive interval evaluation when F' has mixed-sign terms.
    """
    m = _thin(mpf(interval.mid))
    return F(m) + Fp(interval) * (interval - m)


def contraction_certificate_centered(F, Fp, interval) -> dict:
    """Stage-2 Banach check with a centered image enclosure and sup|F'| bound.

    Mirrors toy_readback.contraction_certificate; the image uses the mean-value
    form and the Lipschitz bound uses max(|lo|, |hi|) of the derivative enclosure."""
    image = centered_image(F, Fp, interval)
    self_map_pass = _contained(image, interval)
    derivative_enclosure = Fp(interval)
    d_lo, d_hi = _endpoints(derivative_enclosure)
    abs_L = max(abs(d_lo), abs(d_hi))
    return {
        "interval": _interval_json(interval),
        "image_form": "centered mean-value: F(m) + F'(I)*(I - m)",
        "image": _interval_json(image),
        "self_map_pass": bool(self_map_pass),
        "derivative_enclosure": _interval_json(derivative_enclosure),
        "lipschitz_bound_L": mp.nstr(abs_L, 30),
        "lipschitz_pass": bool(abs_L < 1),
        "monotone_nonnegative_pass": bool(d_lo >= 0),
        "banach_pass": bool(self_map_pass and abs_L < 1),
    }


def fixed_point_enclosure_centered(
    F, Fp, interval, half_width: str, max_iterations: int = 5000, lipschitz_L: mpf | None = None
) -> dict:
    """Banach iteration to the fixed point, then a centered-form box enclosure.

    The stop tolerance is scaled by (1 - L) so that the remaining distance to the
    fixed point, bounded by step/(1 - L), stays below half_width/100."""
    lo, hi = _endpoints(interval)
    point = (lo + hi) / 2
    slack = (1 - lipschitz_L) if (lipschitz_L is not None and lipschitz_L < 1) else mpf(1)
    tolerance = mpf(half_width) * slack / 100
    iterations = 0
    for _ in range(max_iterations):
        next_point = mpf(F(_thin(point)).mid)
        iterations += 1
        if abs(next_point - point) < tolerance:
            point = next_point
            break
        point = next_point
    h = mpf(half_width)
    box = iv.mpf([mp.nstr(point - h, mp.dps), mp.nstr(point + h, mp.dps)])
    box_image = centered_image(F, Fp, box)
    box_self_map_pass = _contained(box_image, box)
    image_lo, image_hi = _endpoints(box_image)
    return {
        "located": bool(box_self_map_pass),
        "iterations": iterations,
        "candidate_box": _interval_json(box),
        "box_self_map_pass": bool(box_self_map_pass),
        "enclosure": _interval_json(box_image),
        "width": mp.nstr(image_hi - image_lo, 8),
    }


def make_map(rho, c, d):
    def F(x):
        return rho * x + c * iv.log(x) + d

    def Fprime(x):
        return rho + c / x

    return F, Fprime


def locate_stable_root(rho_mid: mpf, c: int, d_mid: mpf):
    """Scalar analysis of g(N) = (rho-1)N + c*log N + d on (0, inf).

    Returns (status, root_or_None, note). status in
    {candidate, no_fixed_point, rejected_no_contraction, rejected_trivial}.
    """
    one = mpf(1)
    if rho_mid == one:
        if c > 0:
            return ("rejected_no_contraction", None,
                    "rho = 1, c > 0: F' = 1 + c/N > 1 everywhere; no stable fixed point")
        if c == 0:
            if d_mid == 0:
                return ("rejected_trivial", None, "identity map; excluded by P5")
            return ("no_fixed_point", None, "translation by d != 0; no fixed point")
        root = mp.e ** (d_mid / abs(c))
        return ("candidate", root, "rho = 1, c < 0: root exp(d/|c|), F' = 1 - |c|/N < 1")

    slope_gap = one - rho_mid  # > 0 on all executed rho < 1 rows

    def g(n):
        return (rho_mid - one) * n + c * mp.log(n) + d_mid

    if c > 0:
        n_peak = mpf(c) / slope_gap
        g_peak = c * (mp.log(n_peak) - 1) + d_mid
        if g_peak <= 0:
            return ("no_fixed_point", None,
                    "sup_N [F(N) - N] = c(log(c/(1-rho)) - 1) + d <= 0: F(N) < N on (0, inf)")
        hi = n_peak * 2
        while g(hi) > 0:
            hi *= 2
        lo = n_peak
        for _ in range(300):
            mid = (lo + hi) / 2
            if g(mid) > 0:
                lo = mid
            else:
                hi = mid
        return ("candidate", (lo + hi) / 2,
                "stable upper root (F' < 1 there); lower root is the unstable companion")
    if c == 0:
        if d_mid <= 0:
            return ("no_fixed_point", None, "F(N) = rho*N + d with d <= 0: no positive fixed point")
        return ("candidate", d_mid / slope_gap, "affine map: unique fixed point d/(1-rho)")
    # c < 0: g strictly decreasing, g(0+) = +inf, unique root.
    lo = mpf("1e-30")
    hi = mpf(2)
    while g(hi) > 0:
        hi *= 2
    for _ in range(300):
        mid = (lo + hi) / 2
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    return ("candidate", (lo + hi) / 2, "c < 0: unique root of strictly decreasing g")


def certify_row(rho, c, d, root: mpf) -> dict:
    F, Fp = make_map(rho, c, d)
    result = None
    for eps in ("0.05", "0.02", "0.01", "0.005", "0.002", "0.001"):
        e = mpf(eps)
        interval = iv.mpf([mp.nstr(root * (1 - e), mp.dps), mp.nstr(root * (1 + e), mp.dps)])
        cert = contraction_certificate_centered(F, Fp, interval)
        if cert["banach_pass"]:
            result = (eps, interval, cert)
            break
    if result is None:
        deriv_at_root = mpf(Fp(iv.mpf(mp.nstr(root, mp.dps))).mid)
        unstable = abs(deriv_at_root) >= 1
        return {
            "banach_certified": False,
            "fixed_point_unstable": bool(unstable),
            "derivative_at_root": mp.nstr(deriv_at_root, 20),
            "reason": (
                "fixed point exists with |F'| >= 1: unstable, fails P2/A4, no Banach certificate"
                if unstable else
                "no certificate interval found near located root"
            ),
        }
    eps, interval, cert = result
    lipschitz_L = mpf(cert["lipschitz_bound_L"])
    enclosure = None
    for half_width in (ENCLOSURE_HALF_WIDTH, "1e-22", "1e-18", "1e-12"):
        enclosure = fixed_point_enclosure_centered(F, Fp, interval, half_width, lipschitz_L=lipschitz_L)
        if enclosure["located"]:
            break
    return {
        "banach_certified": bool(cert["banach_pass"] and enclosure["located"]),
        "certificate_interval_relative_half_width": eps,
        "contraction_certificate": cert,
        "fixed_point": enclosure,
    }


def p4_record(rho, c, d, row_status: str, fp_mid: mpf | None) -> dict:
    rho_hi = _endpoints(rho)[1]
    if rho_hi < 1 and c > 0:
        n_star = c / (1 - rho)
        rec = {
            "interior_stationary_point": True,
            "n_star_enclosure": _interval_json(n_star),
            "ell_second_derivative": "l'' = -c/N^2 < 0 (concave; derivative-sign certificate holds)",
        }
        if fp_mid is not None:
            n_star_mid = mpf(n_star.mid)
            rec["n_crc_equals_n_star"] = bool(abs(fp_mid - n_star_mid) < mpf("1e-12"))
            rec["discrepancy_n_crc_minus_n_star"] = mp.nstr(fp_mid - n_star_mid, 12)
            rec["verdict"] = "registered_discrepancy" if not rec["n_crc_equals_n_star"] else "coherent"
        else:
            rec["verdict"] = "selector exists; no readback fixed point on this row"
        return rec
    return {
        "interior_stationary_point": False,
        "reason": "l'(N) = (rho-1) + c/N has no interior zero for this (rho, c); MAR argmax sits on a boundary",
        "verdict": "p4_fails_structurally",
    }


def build() -> dict:
    iv.dps = PRECISION
    mp.dps = PRECISION
    rho_axis, br3_axis, br4_axis, br5_axis = build_axes()

    rows = []
    for r_id, r_desc, rho in rho_axis:
        for x_id, x_desc, c3, d3 in br3_axis:
            for k_id, k_desc, ck, dk in br4_axis:
                for s_id, s_desc, ds in br5_axis:
                    c = c3 + ck
                    d = d3 + dk + ds
                    branch_id = f"capL.{r_id}.{x_id}.{k_id}.{s_id}"
                    rho_mid, d_mid = mpf(rho.mid), mpf(d.mid)
                    status, root, note = locate_stable_root(rho_mid, c, d_mid)
                    row = {
                        "branch": branch_id,
                        "readings": [r_desc, x_desc, k_desc, s_desc],
                        "rho": _interval_json(rho),
                        "c": c,
                        "d": _interval_json(d),
                        "existence_note": note,
                    }
                    fp_mid = None
                    if status == "candidate":
                        cert = certify_row(rho, c, d, root)
                        row.update(cert)
                        if cert.get("banach_certified"):
                            row["status"] = "fixed_point_certified"
                            lo = mpf(cert["fixed_point"]["enclosure"]["lo"])
                            hi = mpf(cert["fixed_point"]["enclosure"]["hi"])
                            fp_mid = (lo + hi) / 2
                            row["fixed_point_nats"] = mp.nstr(fp_mid, 25)
                        elif cert.get("fixed_point_unstable"):
                            row["status"] = "fixed_point_unstable_rejected"
                        else:
                            row["status"] = "certificate_failed"
                    else:
                        row["status"] = status
                    row["count_density_p4"] = p4_record(rho, c, d, row["status"], fp_mid)
                    rows.append(row)

    by_status: dict[str, int] = {}
    fps = []
    for row in rows:
        by_status[row["status"]] = by_status.get(row["status"], 0) + 1
        if row["status"] == "fixed_point_certified":
            fps.append(mpf(row["fixed_point_nats"]))
    summary = {
        "rows": len(rows),
        "by_status": by_status,
        "certified_fixed_point_min_nats": mp.nstr(min(fps), 20) if fps else None,
        "certified_fixed_point_max_nats": mp.nstr(max(fps), 20) if fps else None,
        "p4_coherent_rows": sum(
            1 for row in rows
            if row["count_density_p4"].get("verdict") == "coherent"
        ),
    }

    return {
        "artifact": ARTIFACT_NAME,
        "specification": "F_READBACK_SPEC.md",
        "derivation": "F_CONSTRUCTION_2026-07-14.md",
        "family": "CAP-L log-count readback: F(N) = rho*N + c*log N + d, unit nat normalization",
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": PRECISION,
            "rounding": "mpmath_interval_outward",
            "promotion_backend_required": "arb_or_mpfi_directed_outward",
        },
        "inputs": {
            "P_certified_enclosure": {"lo": P_LO, "hi": P_HI},
            "P_source": "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json",
            "structure_integers": [12, 24, 6, 20, 60],
        },
        "blindness": {
            "inputs": [
                "certified forward-closure P enclosure (source-side computed fixed point)",
                "declared structure integers 12/24/6/20/60 and pi",
            ],
            "reads_measured_lambda": False,
            "reads_sl4_estimate": False,
            "reads_cl3_bridge_value": False,
            "reads_alpha_U": False,
            "dependency_cone": ["mpmath", "toy_readback interval machinery", "P_derivation certified enclosure"],
        },
        "summary": summary,
        "rows": rows,
        "moves_cl7": False,
        "cl7_status": "open",
    }


def main() -> int:
    certificate = build()
    out = Path(__file__).resolve().parent / "runtime" / "F_candidate_capL_certificates.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    out.write_text(text, encoding="utf-8")
    print(json.dumps(certificate["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
