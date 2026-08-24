"""Integer-k Kerr ringdown comb template. Build-stage numeric instrument.

Source of record (not modified here):
`falsification/frozen_targets/fz01_2026-07-17/frozen_target_integer_k_comb_2026-07-17.md`
with companion statement
`proof/epic_wins/ringdown_comb/INTEGER_K_COMB_STATEMENT.md`, under the
scientific-status erratum
`falsification/frozen_targets/fz01_2026-07-17/SCIENTIFIC_STATUS_ERRATUM_2026-07-29.md`.

Frozen ratio law. Candidate spectral features above the rotation line
satisfy, under the comb hypothesis,

    (f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) = ln(k_a)/ln(k_b)

for integers k >= 2, equivalently universal-coordinate positions
x_k = ln(k)/(8*pi) with x = (G*M/(c^3*g(chi))) * (omega - m*Omega_H).
Secondary structure: the (k-1)/k KMS weight hierarchy and the
mass-independent fractional linewidth 64*pi^2*p_0/(a*ln(k)) with declared
a in [1, 10].

Kerr horizon functions (geometric input, SI output). With
s(chi) = sqrt(1 - chi^2):

    r_plus      = (G*M/c^2) * (1 + s(chi))
    Omega_H     = c^3 * chi     / (2*G*M*(1 + s(chi)))      [rad/s]
    kappa       = c^3 * s(chi)  / (2*G*M*(1 + s(chi)))      [1/s]

Derivation of g(chi) from the frozen law's universal coordinate (the
KMS/temperature reading pinned by the companion statement). The first-law
discreteness condition of the statement is

    hbar*(omega - m*Omega_H) = k_B*T_H * ln(k),

and with the Hawking temperature k_B*T_H = hbar*kappa/(2*pi) (kappa in
1/s as above) this reads

    omega - m*Omega_H = (kappa/(2*pi)) * ln(k).

The frozen universal coordinate x = (G*M/(c^3*g(chi)))*(omega - m*Omega_H)
takes the value ln(k)/(8*pi) on these lines exactly when

    kappa/(2*pi) = (c^3*g(chi)/(G*M)) / (8*pi),

that is g(chi) = 4*G*M*kappa/c^3 = 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)),
which is the explicit g(chi) of the companion statement. g(chi) is
therefore statement-pinned, not a free selection; the receipt records the
identity g(chi) = 4*G*M*kappa/c^3 as its provenance. Note that the
alternative normalization g = G*M*kappa/c^3 is inconsistent with
x_k = ln(k)/(8*pi) by a factor of 4 and is not used.

Factor-of-2*pi bookkeeping for tooth frequencies in Hz. With
omega = 2*pi*f,

    omega_k - m*Omega_H = (c^3*g(chi)/(G*M)) * x_k
                        = (c^3*g(chi)/(G*M)) * ln(k)/(8*pi),
    Delta_f_k = (omega_k - m*Omega_H)/(2*pi)
              = c^3*g(chi)*ln(k) / (16*pi^2*G*M),
    f_{k,m}   = m*Omega_H/(2*pi) + Delta_f_k.

Constants. Only exact definitional constants enter: the SI defined
c = 299792458 m/s and the IAU 2015 Resolution B3 nominal solar mass
parameter (GM)_sun = 1.3271244e20 m^3/s^2. Masses are parameterized in
nominal solar masses so that G*M = mass_solar * (GM)_sun with no separate
measured G; this parameterization is a declared selection recorded in the
receipt.

Declared selections (all recorded in the receipt): the linewidth nuisance
range a in [1, 10] and the display endpoints a in {1, 10}; the Page
emission coefficient p_0 = 2e-4 as pinned in the companion statement; the
tooth range k in {2, ..., 12} matching the frozen target's KILL-condition
ladder set; the mass parameterization above; and the synthetic reference
point M = 62 nominal solar masses, chi = 0.67, m = 2, which is a
synthetic reference, not an event fit, and matches no published remnant
posterior.

Numeric discipline. All receipt numbers are computed in decimal
arithmetic at 50 significant digits; pi is computed from the Machin
formula pi = 16*arctan(1/5) - 4*arctan(1/239); logarithms and square
roots use the decimal library primitives. Rendered strings carry 40
significant digits; float renderings live under `derived_for_display`.
The canonical serialization is sorted-key, separator-minimal JSON with a
trailing newline, no timestamps, and no machine paths, so the receipt is
checksum-stable. The independent verifier
`verify_integer_k_comb_independent.py` re-derives every number through
different numeric primitives and byte-compares the canonical JSON.

What is not proved here. This module is a build-stage template
instrument, target-blind by construction: no gravitational-wave event
data, no remnant posterior, no detector likelihood, and no comparison
dataset is read, fetched, or evaluated anywhere in this directory. The
physical reading of the template scale (that G*M/(c^3*g(chi)) belongs to
a Kerr remnant) is the frozen target's declared identification, not
derived here. Nothing in this module or its receipt is a registered,
frozen, or scored prediction; the registration contract in this
directory is a draft pending the owner's freeze, and the frozen strain
likelihood, prior normalization, event selection, and trials accounting
demanded by the frozen target are open.
"""

from __future__ import annotations

import hashlib
import json
import os
from decimal import Decimal, getcontext, localcontext

# Working precision (significant decimal digits) for all internal
# arithmetic. Rendered strings carry SIG_DIGITS significant digits.
WORKING_PRECISION = 50
SIG_DIGITS = 40

# Exact definitional constants (strings; parsed into Decimal).
C_LIGHT_M_PER_S = "299792458"  # SI defined value, exact.
GM_SUN_NOMINAL_M3_PER_S2 = "1.3271244E20"  # IAU 2015 B3 nominal, exact.

# Declared selections.
DECLARED_A_RANGE = ("1", "10")  # linewidth nuisance a in [1, 10], declared.
DECLARED_A_DISPLAY = ("1", "10")  # display endpoints, declared.
DECLARED_P0 = "2E-4"  # Page emission coefficient, statement-pinned, declared.
DECLARED_K_MIN = 2
DECLARED_K_MAX = 12  # KILL-condition ladder set {2, ..., 12}, declared.
DECLARED_REFERENCE_MASS_SOLAR = "62"  # synthetic reference, declared.
DECLARED_REFERENCE_CHI = "0.67"  # synthetic reference, declared.
DECLARED_REFERENCE_M_AZIMUTHAL = 2  # synthetic reference, declared.

RECEIPT_BASENAME = "integer_k_comb_template_receipt.json"


def _arctan_inv(x: int) -> Decimal:
    """arctan(1/x) for integer x >= 2 by the alternating Taylor series.

    Terminates when the term magnitude falls below one part in
    10^(prec + 5) of unity; the truncation error is then far below the
    rendered precision."""
    limit = Decimal(1).scaleb(-(getcontext().prec + 5))
    one_over_x = Decimal(1) / Decimal(x)
    x2 = Decimal(x) * Decimal(x)
    term = one_over_x
    total = Decimal(0)
    n = 0
    while term.copy_abs() > limit:
        total += term if n % 2 == 0 else -term
        n += 1
        term = term / x2 * Decimal(2 * n - 1) / Decimal(2 * n + 1)
    return total


def compute_pi() -> Decimal:
    """pi from the Machin formula 16*arctan(1/5) - 4*arctan(1/239)."""
    with localcontext() as ctx:
        ctx.prec = WORKING_PRECISION + 10
        pi_guard = 16 * _arctan_inv(5) - 4 * _arctan_inv(239)
    return +pi_guard


def dec(value: str | int) -> Decimal:
    return Decimal(value)


def sqrt_one_minus_chi_squared(chi: Decimal) -> Decimal:
    """s(chi) = sqrt(1 - chi^2), the Kerr root factor."""
    return (Decimal(1) - chi * chi).sqrt()


def r_plus_hat(chi: Decimal) -> Decimal:
    """Outer horizon radius in units of G*M/c^2: 1 + sqrt(1 - chi^2)."""
    return Decimal(1) + sqrt_one_minus_chi_squared(chi)


def gm_si(mass_solar: Decimal) -> Decimal:
    """G*M in m^3/s^2 from the nominal solar mass parameter."""
    return mass_solar * dec(GM_SUN_NOMINAL_M3_PER_S2)


def r_plus_si(mass_solar: Decimal, chi: Decimal) -> Decimal:
    """Outer horizon radius in meters: (G*M/c^2)*(1 + sqrt(1 - chi^2))."""
    c = dec(C_LIGHT_M_PER_S)
    return gm_si(mass_solar) / (c * c) * r_plus_hat(chi)


def omega_h_si(mass_solar: Decimal, chi: Decimal) -> Decimal:
    """Horizon angular frequency in rad/s:
    c^3*chi / (2*G*M*(1 + sqrt(1 - chi^2)))."""
    c = dec(C_LIGHT_M_PER_S)
    return c ** 3 * chi / (2 * gm_si(mass_solar) * r_plus_hat(chi))


def kappa_si(mass_solar: Decimal, chi: Decimal) -> Decimal:
    """Surface gravity in 1/s:
    c^3*sqrt(1 - chi^2) / (2*G*M*(1 + sqrt(1 - chi^2)))."""
    c = dec(C_LIGHT_M_PER_S)
    s = sqrt_one_minus_chi_squared(chi)
    return c ** 3 * s / (2 * gm_si(mass_solar) * r_plus_hat(chi))


def g_of_chi(chi: Decimal) -> Decimal:
    """Statement-pinned spin factor
    g(chi) = 2*sqrt(1 - chi^2)/(1 + sqrt(1 - chi^2)) = 4*G*M*kappa/c^3."""
    s = sqrt_one_minus_chi_squared(chi)
    return 2 * s / (Decimal(1) + s)


def base_spacing_hz_per_nat(mass_solar: Decimal, chi: Decimal, pi: Decimal) -> Decimal:
    """Tooth spacing per nat of ln(k): c^3*g(chi) / (16*pi^2*G*M) in Hz."""
    c = dec(C_LIGHT_M_PER_S)
    return c ** 3 * g_of_chi(chi) / (16 * pi * pi * gm_si(mass_solar))


def rotation_line_hz(mass_solar: Decimal, chi: Decimal, m: int, pi: Decimal) -> Decimal:
    """Rotation line m*Omega_H/(2*pi) in Hz."""
    return Decimal(m) * omega_h_si(mass_solar, chi) / (2 * pi)


def universal_position(k: int, pi: Decimal) -> Decimal:
    """Frozen universal-coordinate tooth position x_k = ln(k)/(8*pi)."""
    return Decimal(k).ln() / (8 * pi)


def ladder_ratio(k: int) -> Decimal:
    """Offset-subtracted ratio against the k = 2 tooth: ln(k)/ln(2)."""
    return Decimal(k).ln() / Decimal(2).ln()


def kms_weight(k: int) -> Decimal:
    """KMS detailed-balance weight (k-1)/k."""
    return Decimal(k - 1) / Decimal(k)


def tooth_offset_hz(mass_solar: Decimal, chi: Decimal, k: int, pi: Decimal) -> Decimal:
    """Delta_f_k = c^3*g(chi)*ln(k) / (16*pi^2*G*M) in Hz."""
    return base_spacing_hz_per_nat(mass_solar, chi, pi) * Decimal(k).ln()


def tooth_frequency_hz(
    mass_solar: Decimal, chi: Decimal, m: int, k: int, pi: Decimal
) -> Decimal:
    """f_{k,m} = m*Omega_H/(2*pi) + Delta_f_k in Hz (source frame)."""
    return rotation_line_hz(mass_solar, chi, m, pi) + tooth_offset_hz(
        mass_solar, chi, k, pi
    )


def linewidth_fraction(a: Decimal, k: int, pi: Decimal) -> Decimal:
    """Mass-independent fractional linewidth 64*pi^2*p_0/(a*ln(k)).

    Takes no mass argument: the fraction is mass-independent at every
    (a, k) by the frozen statement."""
    p0 = dec(DECLARED_P0)
    return 64 * pi * pi * p0 / (a * Decimal(k).ln())


def sig40(x: Decimal) -> str:
    """Render exactly SIG_DIGITS significant digits in scientific form."""
    return format(x, ".%dE" % (SIG_DIGITS - 1))


def build_receipt() -> dict:
    """Assemble the full receipt dictionary (pure; no I/O)."""
    getcontext().prec = WORKING_PRECISION
    pi = compute_pi()

    mass = dec(DECLARED_REFERENCE_MASS_SOLAR)
    chi = dec(DECLARED_REFERENCE_CHI)
    m_az = DECLARED_REFERENCE_M_AZIMUTHAL
    ks = list(range(DECLARED_K_MIN, DECLARED_K_MAX + 1))

    display: dict[str, float] = {}

    ladder = []
    for k in ks:
        xk = universal_position(k, pi)
        rk = ladder_ratio(k)
        wk = kms_weight(k)
        key = "k%02d" % k
        ladder.append(
            {
                "k": k,
                "x_exact": "ln(%d)/(8*pi)" % k,
                "x_sig40": sig40(xk),
                "ratio_to_k2_exact": "ln(%d)/ln(2)" % k,
                "ratio_to_k2_sig40": sig40(rk),
                "kms_weight_exact": "%d/%d" % (k - 1, k),
                "kms_weight_sig40": sig40(wk),
            }
        )
        display["universal_ladder.%s.x" % key] = float(xk)
        display["universal_ladder.%s.ratio_to_k2" % key] = float(rk)
        display["universal_ladder.%s.kms_weight" % key] = float(wk)

    omega_h = omega_h_si(mass, chi)
    kappa = kappa_si(mass, chi)
    g_chi = g_of_chi(chi)
    rhat = r_plus_hat(chi)
    r_plus_m = r_plus_si(mass, chi)
    rot = rotation_line_hz(mass, chi, m_az, pi)
    base = base_spacing_hz_per_nat(mass, chi, pi)

    display["reference.omega_h_rad_per_s"] = float(omega_h)
    display["reference.kappa_per_s"] = float(kappa)
    display["reference.g_chi"] = float(g_chi)
    display["reference.r_plus_hat"] = float(rhat)
    display["reference.r_plus_m"] = float(r_plus_m)
    display["reference.rotation_line_hz"] = float(rot)
    display["reference.base_spacing_hz_per_nat"] = float(base)

    teeth = []
    a_lo = dec(DECLARED_A_DISPLAY[0])
    a_hi = dec(DECLARED_A_DISPLAY[1])
    for k in ks:
        dfk = tooth_offset_hz(mass, chi, k, pi)
        fk = rot + dfk
        lw_lo = linewidth_fraction(a_hi, k, pi)  # a = 10: narrow end
        lw_hi = linewidth_fraction(a_lo, k, pi)  # a = 1: wide end
        key = "k%02d" % k
        teeth.append(
            {
                "k": k,
                "delta_f_hz_sig40": sig40(dfk),
                "f_hz_sig40": sig40(fk),
                "linewidth_fraction_a10_sig40": sig40(lw_lo),
                "linewidth_fraction_a1_sig40": sig40(lw_hi),
            }
        )
        display["reference.teeth.%s.delta_f_hz" % key] = float(dfk)
        display["reference.teeth.%s.f_hz" % key] = float(fk)
        display["reference.teeth.%s.linewidth_fraction_a10" % key] = float(lw_lo)
        display["reference.teeth.%s.linewidth_fraction_a1" % key] = float(lw_hi)

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
            "c_m_per_s": C_LIGHT_M_PER_S,
            "c_provenance": "SI defined value, exact",
            "gm_sun_nominal_m3_per_s2": GM_SUN_NOMINAL_M3_PER_S2,
            "gm_sun_provenance": (
                "IAU 2015 Resolution B3 nominal solar mass parameter, "
                "exact nominal value"
            ),
        },
        "declared_selections": {
            "a_range": {
                "value": list(DECLARED_A_RANGE),
                "status": "declared; frozen-statement nuisance interval",
            },
            "a_display_endpoints": {
                "value": list(DECLARED_A_DISPLAY),
                "status": "declared; display endpoints of the a range",
            },
            "p_0": {
                "value": DECLARED_P0,
                "status": (
                    "declared; Page emission coefficient as pinned in the "
                    "companion statement"
                ),
            },
            "k_range": {
                "value": [DECLARED_K_MIN, DECLARED_K_MAX],
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
                    "mass_nominal_solar": DECLARED_REFERENCE_MASS_SOLAR,
                    "chi": DECLARED_REFERENCE_CHI,
                    "m_azimuthal": DECLARED_REFERENCE_M_AZIMUTHAL,
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
            "mass_nominal_solar": DECLARED_REFERENCE_MASS_SOLAR,
            "chi": DECLARED_REFERENCE_CHI,
            "m_azimuthal": DECLARED_REFERENCE_M_AZIMUTHAL,
            "frame": "source frame; detector frame carries 1/(1+z)",
            "omega_h_rad_per_s_sig40": sig40(omega_h),
            "kappa_per_s_sig40": sig40(kappa),
            "g_chi_sig40": sig40(g_chi),
            "r_plus_hat_sig40": sig40(rhat),
            "r_plus_m_sig40": sig40(r_plus_m),
            "rotation_line_hz_sig40": sig40(rot),
            "base_spacing_hz_per_nat_sig40": sig40(base),
            "teeth": teeth,
        },
        "derived_for_display": display,
        "numerics": {
            "working_precision_decimal_digits": WORKING_PRECISION,
            "rendered_significant_digits": SIG_DIGITS,
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


def canonical_bytes(receipt: dict) -> bytes:
    """Canonical serialization: sorted keys, minimal separators, ASCII,
    trailing newline. No timestamps and no machine paths appear."""
    return (
        json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode("ascii")
        + b"\n"
    )


def receipt_path() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "runtime", RECEIPT_BASENAME)


def main() -> int:
    receipt = build_receipt()
    payload = canonical_bytes(receipt)
    path = receipt_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(payload)
    digest = hashlib.sha256(payload).hexdigest()
    print("wrote %s" % os.path.relpath(path, os.path.dirname(os.path.abspath(__file__))))
    print("sha256 %s" % digest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
