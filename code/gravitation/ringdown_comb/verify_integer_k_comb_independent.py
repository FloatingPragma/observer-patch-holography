"""Independent verifier for the integer-k comb template receipt.

Re-derives every number in
`runtime/integer_k_comb_template_receipt.json` from scratch on a code
path disjoint from `integer_k_comb_template.py` and byte-compares the
canonical JSON serializations. The producer is not imported.

Independent numeric primitives used here, against the producer's:

- pi from the Euler decomposition pi/4 = arctan(1/2) + arctan(1/3)
  (producer: Machin 16*arctan(1/5) - 4*arctan(1/239));
- natural logarithms from the atanh series
  ln(k) = 2*sum_{i>=0} z^(2i+1)/(2i+1) with z = (k-1)/(k+1)
  (producer: the decimal library ln primitive);
- square roots by Newton iteration (producer: the decimal library sqrt
  primitive);
- g(chi) through the surface-gravity identity g = 4*G*M*kappa/c^3
  (producer: the closed form 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)));
- tooth spacing through base = kappa/(4*pi^2)
  (producer: c^3*g(chi)/(16*pi^2*G*M));
- working precision 60 significant digits (producer: 50), both rendered
  to 40 significant digits, so agreement requires both computations to
  be correct at the rendered precision.

Exit status 0 on byte equality, 1 on any difference. The SHA-256 of both
serializations is printed.

What is not proved here. Byte agreement certifies that two disjoint
implementations of the frozen formulas produce the same rendered
numbers; it does not certify the physical reading of the template, does
not register a prediction, and touches no event data or comparison
dataset.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from decimal import Decimal, getcontext, localcontext

V_PRECISION = 60
V_SIG_DIGITS = 40

V_C = "299792458"
V_GM_SUN = "1.3271244E20"
V_P0 = "2E-4"
V_A_RANGE = ("1", "10")
V_K_MIN = 2
V_K_MAX = 12
V_MASS_SOLAR = "62"
V_CHI = "0.67"
V_M_AZ = 2

V_RECEIPT_BASENAME = "integer_k_comb_template_receipt.json"


def v_arctan_inv(x: int) -> Decimal:
    """arctan(1/x) by the alternating Taylor series (guarded context).

    Terminates when the power magnitude falls below one part in
    10^(prec + 5) of unity."""
    limit = Decimal(1).scaleb(-(getcontext().prec + 5))
    inv = Decimal(1) / Decimal(x)
    x2 = Decimal(x) * Decimal(x)
    term = inv
    total = Decimal(0)
    sign = 1
    n = 1
    while term.copy_abs() > limit:
        total += sign * term / Decimal(n)
        term = term / x2
        sign = -sign
        n += 2
    return total


def v_pi() -> Decimal:
    """pi from pi/4 = arctan(1/2) + arctan(1/3)."""
    with localcontext() as ctx:
        ctx.prec = V_PRECISION + 12
        quarter = v_arctan_inv(2) + v_arctan_inv(3)
        four_pi = 4 * quarter
    return +four_pi


def v_ln_nat(k: int) -> Decimal:
    """ln(k) for natural k >= 1 by the atanh series, no library ln."""
    if k == 1:
        return Decimal(0)
    with localcontext() as ctx:
        ctx.prec = V_PRECISION + 12
        limit = Decimal(1).scaleb(-(ctx.prec + 5))
        z = Decimal(k - 1) / Decimal(k + 1)
        z2 = z * z
        power = z
        total = Decimal(0)
        n = 1
        while power > limit:
            total += power / Decimal(n)
            power = power * z2
            n += 2
        result = 2 * total
    return +result


def v_sqrt(x: Decimal) -> Decimal:
    """Square root by Newton iteration, no library sqrt.

    Newton iteration doubles the digit count per step; eight steps from
    a sixteen-digit float seed exceed the guarded working precision by a
    wide margin, so a fixed iteration count is used and no convergence
    test is needed."""
    if x == 0:
        return Decimal(0)
    with localcontext() as ctx:
        ctx.prec = V_PRECISION + 12
        y = Decimal(float(x) ** 0.5)
        for _ in range(8):
            y = (y + x / y) / 2
        result = y
    return +result


def v_sig40(x: Decimal) -> str:
    return format(x, ".%dE" % (V_SIG_DIGITS - 1))


def v_build_receipt() -> dict:
    getcontext().prec = V_PRECISION
    pi = v_pi()
    c = Decimal(V_C)
    gm = Decimal(V_MASS_SOLAR) * Decimal(V_GM_SUN)
    chi = Decimal(V_CHI)
    p0 = Decimal(V_P0)
    ks = list(range(V_K_MIN, V_K_MAX + 1))
    ln = {k: v_ln_nat(k) for k in ks + [2]}

    s = v_sqrt(Decimal(1) - chi * chi)
    rhat = Decimal(1) + s
    kappa = c ** 3 * s / (2 * gm * rhat)
    # Independent route: surface-gravity identity, not the closed form.
    g_chi = 4 * gm * kappa / (c ** 3)
    omega_h = (c ** 3 / (2 * gm)) * (chi / rhat)
    r_plus_m = gm / (c * c) * rhat
    rot = Decimal(V_M_AZ) * omega_h / (2 * pi)
    # Independent route: base spacing through kappa/(4*pi^2).
    base = kappa / (4 * pi * pi)

    display: dict[str, float] = {}
    ladder = []
    for k in ks:
        xk = ln[k] / (8 * pi)
        rk = ln[k] / ln[2]
        wk = Decimal(k - 1) / Decimal(k)
        key = "k%02d" % k
        ladder.append(
            {
                "k": k,
                "x_exact": "ln(%d)/(8*pi)" % k,
                "x_sig40": v_sig40(xk),
                "ratio_to_k2_exact": "ln(%d)/ln(2)" % k,
                "ratio_to_k2_sig40": v_sig40(rk),
                "kms_weight_exact": "%d/%d" % (k - 1, k),
                "kms_weight_sig40": v_sig40(wk),
            }
        )
        display["universal_ladder.%s.x" % key] = float(xk)
        display["universal_ladder.%s.ratio_to_k2" % key] = float(rk)
        display["universal_ladder.%s.kms_weight" % key] = float(wk)

    display["reference.omega_h_rad_per_s"] = float(omega_h)
    display["reference.kappa_per_s"] = float(kappa)
    display["reference.g_chi"] = float(g_chi)
    display["reference.r_plus_hat"] = float(rhat)
    display["reference.r_plus_m"] = float(r_plus_m)
    display["reference.rotation_line_hz"] = float(rot)
    display["reference.base_spacing_hz_per_nat"] = float(base)

    teeth = []
    a_lo = Decimal(V_A_RANGE[0])
    a_hi = Decimal(V_A_RANGE[1])
    for k in ks:
        dfk = base * ln[k]
        fk = rot + dfk
        lw_a10 = 64 * pi * pi * p0 / (a_hi * ln[k])
        lw_a1 = 64 * pi * pi * p0 / (a_lo * ln[k])
        key = "k%02d" % k
        teeth.append(
            {
                "k": k,
                "delta_f_hz_sig40": v_sig40(dfk),
                "f_hz_sig40": v_sig40(fk),
                "linewidth_fraction_a10_sig40": v_sig40(lw_a10),
                "linewidth_fraction_a1_sig40": v_sig40(lw_a1),
            }
        )
        display["reference.teeth.%s.delta_f_hz" % key] = float(dfk)
        display["reference.teeth.%s.f_hz" % key] = float(fk)
        display["reference.teeth.%s.linewidth_fraction_a10" % key] = float(lw_a10)
        display["reference.teeth.%s.linewidth_fraction_a1" % key] = float(lw_a1)

    receipt = {
        "schema": "oph.ringdown.integer_k_comb_template.v1",
        "status": (
            "build-stage instrument, target-blind draft; not a registered, "
            "frozen, or scored prediction"
        ),
        "frozen_law": {
            "source_of_record": (
                "falsification/frozen_targets/fz01_2026-07-17/"
                "frozen_target_integer_k_comb_2026-07-17.md"
            ),
            "companion_statement": (
                "proof/epic_wins/ringdown_comb/INTEGER_K_COMB_STATEMENT.md"
            ),
            "erratum": (
                "falsification/frozen_targets/fz01_2026-07-17/"
                "SCIENTIFIC_STATUS_ERRATUM_2026-07-29.md"
            ),
            "ratio_law": (
                "(f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) "
                "= ln(k_a)/ln(k_b), integers k >= 2"
            ),
            "universal_coordinate": (
                "x = (G*M/(c^3*g(chi))) * (omega - m*Omega_H); "
                "x_k = ln(k)/(8*pi)"
            ),
            "g_of_chi": (
                "g(chi) = 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)) "
                "= 4*G*M*kappa(M,chi)/c^3; statement-pinned via the KMS "
                "reading omega - m*Omega_H = (kappa/(2*pi))*ln(k)"
            ),
            "tooth_offset": "Delta_f_k = c^3*g(chi)*ln(k)/(16*pi^2*G*M)",
            "kms_weight": "(k-1)/k",
            "linewidth_fraction": "64*pi^2*p_0/(a*ln(k))",
            "kerr_functions": {
                "r_plus": "(G*M/c^2)*(1 + sqrt(1-chi^2))",
                "omega_h": "c^3*chi/(2*G*M*(1 + sqrt(1-chi^2)))",
                "kappa": "c^3*sqrt(1-chi^2)/(2*G*M*(1 + sqrt(1-chi^2)))",
            },
        },
        "constants_exact": {
            "c_m_per_s": V_C,
            "c_provenance": "SI defined value, exact",
            "gm_sun_nominal_m3_per_s2": V_GM_SUN,
            "gm_sun_provenance": (
                "IAU 2015 Resolution B3 nominal solar mass parameter, "
                "exact nominal value"
            ),
        },
        "declared_selections": {
            "a_range": {
                "value": list(V_A_RANGE),
                "status": "declared; frozen-statement nuisance interval",
            },
            "a_display_endpoints": {
                "value": list(V_A_RANGE),
                "status": "declared; display endpoints of the a range",
            },
            "p_0": {
                "value": V_P0,
                "status": (
                    "declared; Page emission coefficient as pinned in the "
                    "companion statement"
                ),
            },
            "k_range": {
                "value": [V_K_MIN, V_K_MAX],
                "status": (
                    "declared; matches the frozen target KILL-condition "
                    "ladder set {2, ..., 12}"
                ),
            },
            "mass_parameterization": {
                "value": "G*M = mass_solar * (GM)_sun nominal",
                "status": (
                    "declared; nominal solar mass parameter in place of a "
                    "separate measured G"
                ),
            },
            "g_of_chi_provenance": {
                "value": "statement-pinned, derived from the KMS reading",
                "status": (
                    "not a free selection; the identity "
                    "g = 4*G*M*kappa/c^3 is recorded in frozen_law.g_of_chi"
                ),
            },
            "reference_point": {
                "value": {
                    "mass_nominal_solar": V_MASS_SOLAR,
                    "chi": V_CHI,
                    "m_azimuthal": V_M_AZ,
                },
                "status": (
                    "declared; synthetic-reference, not an event fit; "
                    "matches no published remnant posterior"
                ),
            },
        },
        "universal_ladder": ladder,
        "reference_point_synthetic": {
            "label": "synthetic-reference, not an event fit",
            "mass_nominal_solar": V_MASS_SOLAR,
            "chi": V_CHI,
            "m_azimuthal": V_M_AZ,
            "frame": "source frame; detector frame carries 1/(1+z)",
            "omega_h_rad_per_s_sig40": v_sig40(omega_h),
            "kappa_per_s_sig40": v_sig40(kappa),
            "g_chi_sig40": v_sig40(g_chi),
            "r_plus_hat_sig40": v_sig40(rhat),
            "r_plus_m_sig40": v_sig40(r_plus_m),
            "rotation_line_hz_sig40": v_sig40(rot),
            "base_spacing_hz_per_nat_sig40": v_sig40(base),
            "teeth": teeth,
        },
        "derived_for_display": display,
        "numerics": {
            "working_precision_decimal_digits": 50,
            "rendered_significant_digits": 40,
            "pi_method": "Machin: 16*arctan(1/5) - 4*arctan(1/239)",
        },
        "boundary": (
            "Target-blind build-stage instrument. No gravitational-wave "
            "event data, remnant posterior, detector likelihood, or "
            "comparison dataset was read, fetched, or evaluated. The "
            "physical reading of the template scale as a Kerr remnant "
            "quantity is the frozen target's declared identification. "
            "Registration is open pending the owner's freeze."
        ),
    }
    return receipt


def v_canonical_bytes(receipt: dict) -> bytes:
    return (
        json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode("ascii")
        + b"\n"
    )


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "runtime", V_RECEIPT_BASENAME)
    with open(path, "rb") as fh:
        on_disk = fh.read()
    rebuilt = v_canonical_bytes(v_build_receipt())
    print("on-disk  sha256 %s" % hashlib.sha256(on_disk).hexdigest())
    print("rebuilt  sha256 %s" % hashlib.sha256(rebuilt).hexdigest())
    if on_disk == rebuilt:
        print("VERIFIED: byte-identical canonical receipts")
        return 0
    print("MISMATCH: canonical receipts differ")
    return 1


if __name__ == "__main__":
    sys.exit(main())
