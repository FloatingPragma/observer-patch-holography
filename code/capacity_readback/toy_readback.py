#!/usr/bin/env python3
"""Toy schema demonstration of the OPH capacity readback map F = Cap o Obs o nf.

This is a schema demonstration on a toy model. It carries no physical content
and does not move CL-7. The toy universe, the toy repair surface, the toy
observer sector, and the toy readback normalization are all declared here and
have no relation to the OPH grammar. The point of the file is the pipeline
shape required by F_READBACK_SPEC.md: an explicit finite state space, an
explicit terminating/confluent nf, an explicit Obs selection, an explicit
Cap_read, a real-interval extension, a stage-2 Banach certificate (self-map
enclosure plus interval derivative bound L < 1), a certified fixed-point
enclosure, and a blindness record. A deliberately non-contracting variant
exercises the rejection path.

Toy model
---------
- Universe candidate U_n at integer capacity n: all strings of length n over
  the alphabet {0, 1, 2} (3^n raw states; fully enumerable for small n).
- nf: sort the string (repeated adjacent swaps of out-of-order pairs; a
  terminating, confluent rewriting whose normal forms are the non-decreasing
  strings, i.e. multisets).
- Obs: keep the normal forms containing at least one '0' (free record room)
  and at least one '1' (an observer marker cell). Closed-form sector count:
  |Omega_n| = n(n-1)/2.
- Cap_read: KAPPA * log|sector| + MU with declared toy constants KAPPA, MU.

Real extension: F_toy(N) = MU + KAPPA*log(N(N-1)/2) on N > 1, using the
closed-form sector count. The enumeration check verifies the closed form
against the explicit pipeline on a range of small n.

Negative control: the "non_contracting" variant replaces Cap_read by
MU + 2*sqrt(|sector|); its interval certificate fails both the self-map and
the Lipschitz check and the builder reports status "rejected".
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from mpmath import iv, mp, mpf

ARTIFACT_NAME = "oph_capacity_readback_toy_schema"
DISCLAIMER = (
    "Schema demonstration on a toy model. No physical content. "
    "Does not move CL-7. The physical readback map F is not constructed."
)

ALPHABET = ("0", "1", "2")

# Declared toy readback normalization (toy constants, not derived quantities).
KAPPA = 6
MU = 4

DEFAULT_INTERVAL = ("40", "50")
DEFAULT_PRECISION = 40
DEFAULT_ENUMERATION_RANGE = (3, 10)
DEFAULT_ENCLOSURE_HALF_WIDTH = "1e-25"


# ---------------------------------------------------------------------------
# Finite toy pipeline: explicit state space, nf, Obs, Cap_read.
# ---------------------------------------------------------------------------

def toy_universe(n: int):
    """All raw boundary states at integer capacity n."""
    return ("".join(cells) for cells in itertools.product(ALPHABET, repeat=n))


def toy_nf(state: str) -> str:
    """Quotient normal form: the sorted string.

    The rewriting 'swap any adjacent out-of-order pair' terminates (each step
    lowers the inversion count) and is confluent (all maximal rewrites reach
    the sorted string), so nf is a function.
    """
    return "".join(sorted(state))


def toy_obs(normal_forms: set[str]) -> set[str]:
    """Stable observer sector: normal forms with >= one '0' and >= one '1'."""
    return {s for s in normal_forms if "0" in s and "1" in s}


def toy_sector_count_closed_form(n: int) -> int:
    """|Omega_n| = n(n-1)/2 (multisets over 3 letters with >=1 zero, >=1 one)."""
    return n * (n - 1) // 2


def enumerate_sector_count(n: int) -> int:
    """Run the explicit finite pipeline nf -> Obs and count the sector."""
    normal_forms = {toy_nf(state) for state in toy_universe(n)}
    return len(toy_obs(normal_forms))


def enumeration_check(n_min: int, n_max: int) -> dict:
    counts = {}
    all_match = True
    for n in range(n_min, n_max + 1):
        enumerated = enumerate_sector_count(n)
        closed = toy_sector_count_closed_form(n)
        counts[str(n)] = {"enumerated": enumerated, "closed_form": closed}
        if enumerated != closed:
            all_match = False
    return {
        "range": {"n_min": n_min, "n_max": n_max},
        "raw_state_space_size": "3^n",
        "counts": counts,
        "closed_form": "n*(n-1)/2",
        "closed_form_matches": all_match,
    }


# ---------------------------------------------------------------------------
# Real-interval extension of the toy readback map.
# ---------------------------------------------------------------------------

def toy_readback(x):
    """F_toy(N) = MU + KAPPA*log(N(N-1)/2); valid for interval or scalar N > 1."""
    return MU + KAPPA * iv.log(x * (x - 1) / 2)


def toy_readback_derivative(x):
    """F_toy'(N) = KAPPA*(2N-1)/(N(N-1)); interval-evaluated on the certificate interval."""
    return KAPPA * (2 * x - 1) / (x * (x - 1))


def toy_readback_non_contracting(x):
    """Negative-control readback: MU + 2*sqrt(N(N-1)/2). Slope > 1 for N >= 3."""
    return MU + 2 * iv.sqrt(x * (x - 1) / 2)


def toy_readback_non_contracting_derivative(x):
    """Derivative of the negative-control map: (2N-1)/(2*sqrt(N(N-1)/2))."""
    return (2 * x - 1) / (2 * iv.sqrt(x * (x - 1) / 2))


VARIANTS = {
    "toy": {
        "map": toy_readback,
        "derivative": toy_readback_derivative,
        "cap_read": f"{MU} + {KAPPA}*log|sector|",
    },
    "non_contracting": {
        "map": toy_readback_non_contracting,
        "derivative": toy_readback_non_contracting_derivative,
        "cap_read": f"{MU} + 2*sqrt(|sector|)",
    },
}


# ---------------------------------------------------------------------------
# Interval certificate machinery.
# ---------------------------------------------------------------------------

def _endpoints(x) -> tuple[mpf, mpf]:
    return mpf(x.a), mpf(x.b)


def _interval_json(x) -> dict[str, str]:
    lo, hi = _endpoints(x)
    return {"lo": mp.nstr(lo, 30), "hi": mp.nstr(hi, 30)}


def _contained(inner, outer) -> bool:
    inner_lo, inner_hi = _endpoints(inner)
    outer_lo, outer_hi = _endpoints(outer)
    return outer_lo <= inner_lo and inner_hi <= outer_hi


def contraction_certificate(readback, derivative, interval) -> dict:
    """Stage-2 Banach check: readback(I) subset of I and sup|F'| <= L < 1 on I."""
    image = readback(interval)
    self_map_pass = _contained(image, interval)
    derivative_enclosure = derivative(interval)
    _, derivative_hi = _endpoints(derivative_enclosure)
    derivative_lo, _ = _endpoints(derivative_enclosure)
    lipschitz_pass = derivative_hi < 1
    monotone_pass = derivative_lo >= 0
    return {
        "interval": _interval_json(interval),
        "image": _interval_json(image),
        "self_map_pass": self_map_pass,
        "derivative_enclosure": _interval_json(derivative_enclosure),
        "lipschitz_bound_L": mp.nstr(derivative_hi, 30),
        "lipschitz_pass": bool(lipschitz_pass),
        "monotone_nonnegative_pass": bool(monotone_pass),
        "banach_pass": bool(self_map_pass and lipschitz_pass),
    }


def fixed_point_enclosure(readback, interval, half_width: str, max_iterations: int = 400) -> dict:
    """Locate the fixed point by Banach iteration, then certify a tight enclosure.

    The point iteration runs at scalar precision from the interval midpoint.
    The enclosure step re-runs the interval self-map check on the small box
    J = [n - h, n + h]; readback(J) subset of J plus the outer certificate's
    contraction bound places the unique fixed point inside readback(J).
    """
    lo, hi = _endpoints(interval)
    point = (lo + hi) / 2
    tolerance = mpf(half_width) / 100
    iterations = 0
    for _ in range(max_iterations):
        next_point = mpf(readback(iv.mpf(mp.nstr(point, mp.dps))).mid)
        iterations += 1
        if abs(next_point - point) < tolerance:
            point = next_point
            break
        point = next_point
    h = mpf(half_width)
    box = iv.mpf([mp.nstr(point - h, mp.dps), mp.nstr(point + h, mp.dps)])
    box_image = readback(box)
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


def count_density_record(n_min: int, n_max: int) -> dict:
    """Report the toy count density ell(n) = log|Omega_n| - n on the enumerated range.

    The toy sector count grows polynomially, so ell is strictly decreasing on
    the toy range and has no interior stationary point; the count-density
    selector is reported as a schema field with no verdict. Coherence of the
    readback fixed point with the count-density stationary point is a property
    required of the physical F (spec property P4), not enforced by the toy.
    """
    values = {}
    for n in range(n_min, n_max + 1):
        count = toy_sector_count_closed_form(n)
        ell = mp.log(count) - n
        values[str(n)] = mp.nstr(ell, 12)
    return {
        "definition": "ell(n) = log|Omega_n| - n",
        "values": values,
        "stationary_point_in_toy_range": None,
        "note": (
            "Toy count is polynomial; ell has no interior stationary point. "
            "Spec property P4 (coherence with the readback fixed point) binds "
            "the physical F only."
        ),
    }


def build_certificate(
    variant: str,
    interval_lo: str,
    interval_hi: str,
    precision: int,
    enumeration_min: int,
    enumeration_max: int,
    enclosure_half_width: str,
) -> dict:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant: {variant}")
    iv.dps = precision
    mp.dps = precision

    spec = VARIANTS[variant]
    interval = iv.mpf([interval_lo, interval_hi])
    contraction = contraction_certificate(spec["map"], spec["derivative"], interval)

    if contraction["banach_pass"]:
        fixed_point = fixed_point_enclosure(spec["map"], interval, enclosure_half_width)
        status = "pass" if fixed_point["located"] else "rejected"
    else:
        fixed_point = {
            "located": False,
            "reason": "no stage-2 contraction certificate on the stated interval",
        }
        status = "rejected"

    return {
        "artifact": ARTIFACT_NAME,
        "disclaimer": DISCLAIMER,
        "physical_content": False,
        "moves_cl7": False,
        "cl7_status": "open",
        "specification": "F_READBACK_SPEC.md",
        "variant": variant,
        "status": status,
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": precision,
            "rounding": "mpmath_interval_outward",
            "promotion_backend_required": "arb_or_mpfi_directed_outward",
        },
        "toy_model": {
            "alphabet": list(ALPHABET),
            "nf": "sort string (adjacent-swap rewriting; terminating, confluent)",
            "obs": "normal forms containing >= one '0' and >= one '1'",
            "cap_read": spec["cap_read"],
            "kappa": KAPPA,
            "mu": MU,
            "sector_count_closed_form": "n*(n-1)/2",
        },
        "enumeration_check": enumeration_check(enumeration_min, enumeration_max),
        "contraction_certificate": contraction,
        "fixed_point": fixed_point,
        "count_density": count_density_record(enumeration_min, enumeration_max),
        "blindness": {
            "inputs": [
                "toy alphabet {0,1,2}",
                "toy constants kappa, mu",
                "certificate interval endpoints",
            ],
            "reads_measured_lambda": False,
            "reads_sl4_estimate": False,
            "dependency_cone": ["itertools", "mpmath"],
            "note": "V-08 shape only; the toy has no physical comparison basin.",
        },
        "promotion_allowed": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the toy capacity-readback schema end to end and emit a "
            "certificate JSON. Toy model only; no physical content."
        )
    )
    parser.add_argument(
        "--variant",
        choices=sorted(VARIANTS),
        default="toy",
        help="Toy readback variant; 'non_contracting' is the negative control.",
    )
    parser.add_argument("--interval-lo", default=DEFAULT_INTERVAL[0], help="Certificate interval lower endpoint.")
    parser.add_argument("--interval-hi", default=DEFAULT_INTERVAL[1], help="Certificate interval upper endpoint.")
    parser.add_argument("--precision", type=int, default=DEFAULT_PRECISION, help="Decimal precision for interval arithmetic.")
    parser.add_argument("--enumeration-min", type=int, default=DEFAULT_ENUMERATION_RANGE[0], help="Smallest enumerated toy capacity.")
    parser.add_argument("--enumeration-max", type=int, default=DEFAULT_ENUMERATION_RANGE[1], help="Largest enumerated toy capacity.")
    parser.add_argument(
        "--enclosure-half-width",
        default=DEFAULT_ENCLOSURE_HALF_WIDTH,
        help="Half-width of the fixed-point enclosure box.",
    )
    parser.add_argument(
        "--output",
        default=str(Path(__file__).resolve().parent / "runtime" / "toy_readback_certificate.json"),
        help="Path for the JSON certificate.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    certificate = build_certificate(
        variant=args.variant,
        interval_lo=args.interval_lo,
        interval_hi=args.interval_hi,
        precision=args.precision,
        enumeration_min=args.enumeration_min,
        enumeration_max=args.enumeration_max,
        enclosure_half_width=args.enclosure_half_width,
    )
    text = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
