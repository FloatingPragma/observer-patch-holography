#!/usr/bin/env python3
"""Extended angular fingerprint: invariant-level support and odd blindness.

This producer extends the equal-port angular comb of the issue #643 sprint
receipt from level fourteen to level forty in exact quadratic-field
arithmetic and certifies two structural facts against the icosahedral
invariant tower:

* **Support alignment.** The icosahedral rotation group has harmonic
  invariants exactly at the levels of the generating function
  ``(1 + t^15) / ((1 - t^6)(1 - t^10))``: even invariant levels are the
  representable sums ``l = 6m + 10n`` and odd invariant levels are
  ``l = 15 + 6m + 10n``. On the certified window ``l <= 40`` the
  equal-port comb ``I_l`` is nonzero at an even level exactly when the
  invariant dimension is positive; the even zeros are exactly
  ``{2, 4, 8, 14}``. Every nonzero weight is an exact rational, pinned
  with zero free parameters.
* **Odd blindness.** Every odd-level comb value is exactly zero, while
  the odd invariant tower is nonempty from level fifteen upward. The
  twelve-port equal-weight readback is therefore exactly blind to the
  entire parity-odd invariant tower, so any certified parity-odd
  response through this channel falsifies the readback branch outright.

Both facts convert to one-sided conditional kill rules that require no
tuned parameter: certified power at an even non-invariant level in the
declared channel kills the carrier-measure branch, and certified
parity-odd response kills the readback branch. Both stay conditional on
the missing screen-to-sky template map recorded by the parent receipt,
and no comparison surface is opened here. Completions that share the
icosahedral support (any icosahedrally symmetric carrier) leave the
weights free; the equal-port comb pins every weight exactly, which is
the discriminating content of this extension.

The window boundary is typed: the support alignment is certified by
exhaustion for ``l <= 40`` and claimed for no higher level.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "angular_fingerprint_receipt.json"
PARENT_RECEIPT_PATH = RUNTIME / "angular_template_receipt.json"

SCHEMA = "oph.angular_fingerprint_receipt.v1"
STATUS = "EXACT_INVARIANT_SUPPORT_FINGERPRINT__CONDITIONAL_KILL_RULES_TYPED"
WINDOW = 40
PARENT_MAX_LEVEL = 14


class FingerprintError(ValueError):
    """The fingerprint certificate refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FingerprintError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Exact Q(sqrt5) arithmetic
# ---------------------------------------------------------------------------

Q5 = tuple[Fraction, Fraction]


def q5(a, b=0) -> Q5:
    return (Fraction(a), Fraction(b))


def q5_add(x: Q5, y: Q5) -> Q5:
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x: Q5, y: Q5) -> Q5:
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x: Q5, y: Q5) -> Q5:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_scale(x: Q5, factor: Fraction) -> Q5:
    return (x[0] * factor, x[1] * factor)


def legendre_values(t: Q5, max_level: int) -> list[Q5]:
    values = [q5(1), t]
    for level in range(1, max_level):
        term = q5_scale(q5_mul(t, values[level]), Fraction(2 * level + 1))
        term = q5_sub(term, q5_scale(values[level - 1], Fraction(level)))
        values.append(q5_scale(term, Fraction(1, level + 1)))
    return values[: max_level + 1]


def equal_port_comb(max_level: int) -> list[Fraction]:
    """Exact ``I_l`` for ``l <= max_level``; fails if any value leaves Q."""

    plus = legendre_values(q5(0, Fraction(1, 5)), max_level)
    minus = legendre_values(q5(0, Fraction(-1, 5)), max_level)
    sequence = []
    for level in range(max_level + 1):
        value = q5_scale(
            q5_add(
                q5(1 + (-1) ** level),
                q5_scale(q5_add(plus[level], minus[level]), Fraction(5)),
            ),
            Fraction(1, 12),
        )
        require(value[1] == 0, f"comb value irrational at level {level}")
        sequence.append(value[0])
    return sequence


# ---------------------------------------------------------------------------
# Icosahedral invariant dimensions by exhaustion
# ---------------------------------------------------------------------------


def even_invariant_dimension(level: int) -> int:
    """Number of pairs ``(m, n) >= 0`` with ``6m + 10n = level``."""

    return sum(
        1
        for m in range(level // 6 + 1)
        for n in range(level // 10 + 1)
        if 6 * m + 10 * n == level
    )


def odd_invariant_dimension(level: int) -> int:
    """Number of pairs ``(m, n) >= 0`` with ``15 + 6m + 10n = level``."""

    if level < 15:
        return 0
    return even_invariant_dimension(level - 15)


# ---------------------------------------------------------------------------
# Certificates
# ---------------------------------------------------------------------------


def parent_pin() -> dict[str, Any]:
    payload = PARENT_RECEIPT_PATH.read_bytes()
    return {
        "path": "code/angular_sprint/runtime/angular_template_receipt.json",
        "bytes": len(payload),
        "sha256": tagged_sha256(payload),
    }


def build_certificates() -> dict[str, Any]:
    sequence = equal_port_comb(WINDOW)

    parent = json.loads(PARENT_RECEIPT_PATH.read_text(encoding="utf-8"))
    parent_sequence = [
        Fraction(value) for value in parent["equal_port_certificate"]["sequence"]
    ]
    require(
        sequence[: PARENT_MAX_LEVEL + 1] == parent_sequence,
        "extension disagrees with the pinned parent comb",
    )

    odd_all_zero = all(
        sequence[level] == 0 for level in range(1, WINDOW + 1, 2)
    )
    require(odd_all_zero, "an odd comb level is nonzero")

    support_rows = []
    even_zero_levels = []
    alignment = True
    for level in range(0, WINDOW + 1, 2):
        dimension = even_invariant_dimension(level)
        nonzero = sequence[level] != 0
        aligned = nonzero == (dimension >= 1)
        alignment = alignment and aligned
        if not nonzero:
            even_zero_levels.append(level)
        support_rows.append(
            {
                "level": level,
                "weight": str(sequence[level]),
                "invariant_dimension": dimension,
                "aligned": aligned,
            }
        )
    require(alignment, "comb support misaligns with the invariant levels")
    require(
        even_zero_levels == [2, 4, 8, 14],
        "even zero set differs from {2, 4, 8, 14} on the window",
    )

    odd_tower = [
        level
        for level in range(WINDOW + 1)
        if level % 2 == 1 and odd_invariant_dimension(level) >= 1
    ]
    require(odd_tower[:1] == [15], "odd invariant tower does not start at 15")

    return {
        "window": WINDOW,
        "generating_function": "(1 + t^15) / ((1 - t^6)(1 - t^10))",
        "support_alignment": {
            "rows": support_rows,
            "even_zero_levels": even_zero_levels,
            "aligned_on_window": alignment,
            "boundary": (
                "certified by exhaustion for levels at most forty; no "
                "claim at higher levels"
            ),
        },
        "odd_blindness": {
            "all_odd_comb_levels_zero": odd_all_zero,
            "odd_invariant_levels_on_window": odd_tower,
            "statement": (
                "the equal-weight twelve-port readback carries exactly "
                "zero response on every parity-odd level while the odd "
                "invariant tower is nonempty from level fifteen"
            ),
        },
        "kill_rules": {
            "even_non_invariant_power": (
                "certified power at an even level in {2, 4, 8, 14} in the "
                "declared channel falsifies the equal-port carrier-measure "
                "branch"
            ),
            "parity_odd_response": (
                "certified parity-odd response through the port readback "
                "falsifies the readback branch"
            ),
            "typing": (
                "both rules are one-sided, parameter-free, and conditional "
                "on the screen-to-sky template map recorded open by the "
                "parent receipt; no comparison surface is opened here"
            ),
        },
        "discrimination_note": (
            "an icosahedrally symmetric carrier fixes the support and "
            "leaves every weight free; the equal-port comb pins all "
            "weights exactly with zero free parameters, and an isotropic "
            "carrier has no comb"
        ),
    }


def build_receipt() -> dict[str, Any]:
    certificates = build_certificates()
    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 643,
        "parent_pins": [parent_pin()],
        "certificates": certificates,
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
    }
    receipt["receipt_sha256"] = tagged_sha256(
        canonical_json_bytes(
            {k: v for k, v in receipt.items() if k != "receipt_sha256"}
        )
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["certificates"]["support_alignment"], indent=1)[-400:])
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
