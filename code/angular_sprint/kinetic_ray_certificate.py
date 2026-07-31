#!/usr/bin/env python3
"""Issue #646 fast-falsification lane: port-coordinate kinetic-ray arithmetic.

The registered port-current pairing carries the exact Hilbert--Schmidt band
coefficients of the pinned port receipt: unit ``1/4``, frame ``5+sqrt5``,
kernel ``5-sqrt5``, quintet ``3+sqrt5``. Per the pinned block realization,
the even block is ``u(3)`` carrying the unit, frame, and quintet bands, and
the kernel block is ``so(3)``: the current algebra splits as
``u(1) + su(2) + su(3)`` with the unit band on ``u(1)``, the kernel band on
``su(2)``, and the frame and quintet bands jointly on the eight-dimensional
``su(3)``.

Typing, superseding the v1 wording of this file: the pairing
``B(f,h) = -Re tr(K(f) K(h))`` is ad-invariant by trace cyclicity, verified
exactly on the full structure-constant table by
``kinetic_form_selection_certificate.py``. The distinct frame and quintet
band coefficients are eigenvalue data of the pullback relative to the port
Euclidean metric; they witness that the port coordinate metric and the
representation trace metric are not isometrically identified, not an
invariance failure on the simple ``su(3)`` ideal. This producer therefore
enumerates port-coordinate ray arithmetic:

* the two raw block rays read the band coefficients per port-normalized
  band, one branch per ``su(3)`` block;
* the dimension-weighted band average
  ``(3 k_frame + 5 k_quintet)/8 = (15+4 sqrt5)/4`` is retained as
  port-metric arithmetic; no projection metric or repair operation selects
  it, and it is not the canonical invariant coefficient (the Killing-relative
  coefficients live in the selection receipt);
* the tested reference ray is the representation-index ray ``(5/3, 1, 1)``.

The exact tests run without any measured value: whether any single overall
scale carries a port-coordinate ray onto the reference ray, and whether the
conditional quadratic-commutant relation ``k1 = 3 k2 - 2 k3`` holds on any
port-coordinate ray. Neither test refutes the representation-index ray: that
ray arises from the rank-fifteen matter trace form, a different quadratic
form whose Killing-relative arithmetic and frozen renormalization-line
statistic are recorded in ``kinetic_form_selection_certificate.py``. The
statistic's beta column is frozen there with exact values; this file records
the same frozen columns and defers ownership to that receipt.

No public measurement is read and no comparison is opened.
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
RECEIPT_PATH = RUNTIME / "kinetic_ray_receipt.json"
PORT_RECEIPT_PATH = (
    REPO_ROOT
    / "code"
    / "a5_closure"
    / "receipts"
    / "port_current_inner_reference.receipt.json"
)

SCHEMA = "oph.kinetic_ray_receipt.v2"
STATUS = (
    "EXACT_PORT_COORDINATE_RAY_ARITHMETIC__"
    "REFERENCE_RAY_DISTINCT_FORM_NOT_REFUTED"
)


class KineticError(ValueError):
    """The kinetic ray certificate refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise KineticError(message)


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


Q5 = tuple[Fraction, Fraction]


def q5(a, b=0) -> Q5:
    return (Fraction(a), Fraction(b))


def q5_add(x: Q5, y: Q5) -> Q5:
    return (x[0] + y[0], x[1] + y[1])


def q5_mul(x: Q5, y: Q5) -> Q5:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_scale(x: Q5, factor: Fraction) -> Q5:
    return (x[0] * factor, x[1] * factor)


def q5_div(x: Q5, y: Q5) -> Q5:
    norm = y[0] * y[0] - 5 * y[1] * y[1]
    numerator = q5_mul(x, (y[0], -y[1]))
    return (numerator[0] / norm, numerator[1] / norm)


def q5_str(x: Q5) -> str:
    return f"{x[0]}+{x[1]}*sqrt5"


def parse_q5(text: str) -> Q5:
    cleaned = text.replace(" ", "").replace("sqrt(5)", "sqrt5")
    head, _, _ = cleaned.partition("*sqrt5")
    if "*sqrt5" not in cleaned:
        return (Fraction(cleaned), Fraction(0))
    a_part, _, b_part = head.rpartition("+")
    if a_part == "":
        a_part, _, b_part = head.rpartition("-")
        if a_part != "":
            b_part = "-" + b_part
        else:
            a_part, b_part = "0", head
    return (Fraction(a_part), Fraction(b_part))


def load_pinned_bands() -> dict[str, Q5]:
    receipt = json.loads(PORT_RECEIPT_PATH.read_text(encoding="utf-8"))
    coefficients = receipt["compactness"][
        "hilbert_schmidt_pullback_band_coefficients"
    ]
    bands = {name: parse_q5(text) for name, text in coefficients.items()}
    require(
        bands["unit_band"] == q5(Fraction(1, 4))
        and bands["frame_band"] == q5(5, 1)
        and bands["kernel_band"] == q5(5, -1)
        and bands["quintet_band"] == q5(3, 1),
        "pinned band coefficient drift",
    )
    return bands


def ideal_decomposition(bands: dict[str, Q5]) -> dict[str, Any]:
    frame = bands["frame_band"]
    quintet = bands["quintet_band"]
    average_su3 = q5_scale(
        q5_add(q5_scale(frame, Fraction(3)), q5_scale(quintet, Fraction(5))),
        Fraction(1, 8),
    )
    require(average_su3 == q5(Fraction(15, 4), Fraction(1)), "su3 average drift")
    non_isometric = frame != quintet
    return {
        "ideal_dimensions": {"u1": 1, "su2": 3, "su3": 8},
        "band_to_ideal": {
            "unit_band": "u1 (even block center)",
            "kernel_band": "su2 (the so(3) kernel block)",
            "frame_band": "su3 (three-dimensional even-block component)",
            "quintet_band": "su3 (five-dimensional even-block component)",
        },
        "block_source": (
            "port_to_generator_map.band_realization in the pinned port "
            "receipt: the even block is u(3) on the unit, frame, and "
            "quintet bands, and the kernel block is so(3)"
        ),
        "su3_blocks_disagree": non_isometric,
        "su3_block_values": [q5_str(frame), q5_str(quintet)],
        "su3_dimension_weighted_average": q5_str(average_su3),
        "invariance_statement": (
            "the pairing is ad-invariant by trace cyclicity; the two "
            "distinct su(3) block coefficients are eigenvalue data of the "
            "pullback relative to the port Euclidean metric and witness the "
            "non-isometry between port coordinates and the representation "
            "trace metric; the dimension-weighted average of the blocks is "
            "port-metric arithmetic, and the canonical invariant "
            "coefficients are the Killing-relative values in the "
            "kinetic-form selection receipt"
        ),
    }


def candidate_rays(bands: dict[str, Q5], decomposition: dict[str, Any]) -> list[dict[str, Any]]:
    average_su3 = parse_q5(decomposition["su3_dimension_weighted_average"])
    return [
        {
            "ray_id": "raw-block-ray-frame-branch",
            "components": [
                q5_str(bands["unit_band"]),
                q5_str(bands["kernel_band"]),
                q5_str(bands["frame_band"]),
            ],
            "note": (
                "reads the three-dimensional su(3) block; declared as one "
                "branch of the non-invariant raw pairing"
            ),
        },
        {
            "ray_id": "raw-block-ray-quintet-branch",
            "components": [
                q5_str(bands["unit_band"]),
                q5_str(bands["kernel_band"]),
                q5_str(bands["quintet_band"]),
            ],
            "note": (
                "reads the five-dimensional su(3) block; the second branch "
                "of the non-invariant raw pairing"
            ),
        },
        {
            "ray_id": "dimension-weighted-average-ray",
            "components": [
                q5_str(bands["unit_band"]),
                q5_str(bands["kernel_band"]),
                q5_str(average_su3),
            ],
            "note": (
                "port-metric arithmetic: unit and kernel band coefficients "
                "with the dimension-weighted su(3) block average; no "
                "projection metric or repair operation selects this "
                "average, and the canonical Killing-relative coefficients "
                "live in the kinetic-form selection receipt"
            ),
        },
    ]


def ray_tests(rays: list[dict[str, Any]]) -> dict[str, Any]:
    reference = [q5(Fraction(5, 3)), q5(1), q5(1)]
    results = []
    for ray in rays:
        components = [parse_q5(text) for text in ray["components"]]
        proportional = True
        scale = q5_div(components[0], reference[0])
        for component, target in zip(components[1:], reference[1:]):
            if q5_div(component, target) != scale:
                proportional = False
        k1, k2, k3 = components
        commutant = k1 == q5_add(
            q5_scale(k2, Fraction(3)), q5_scale(k3, Fraction(-2))
        )
        results.append(
            {
                "ray_id": ray["ray_id"],
                "proportional_to_reference_5_3_1_1": proportional,
                "quadratic_commutant_k1_eq_3k2_minus_2k3": commutant,
                "commutant_left": q5_str(k1),
                "commutant_right": q5_str(
                    q5_add(q5_scale(k2, Fraction(3)), q5_scale(k3, Fraction(-2)))
                ),
            }
        )
    return {
        "reference_ray": ["5/3", "1", "1"],
        "reference_normalization": (
            "components stated per unit-index generator basis; the single "
            "overall ray scale is the only quotiented freedom"
        ),
        "results": results,
        "reference_ray_hit": any(
            row["proportional_to_reference_5_3_1_1"] for row in results
        ),
        "commutant_relation_hit": any(
            row["quadratic_commutant_k1_eq_3k2_minus_2k3"] for row in results
        ),
        "scope": (
            "the port-coordinate rays are distinct from the reference ray; "
            "this does not refute the representation-index ray, which "
            "arises from the rank-fifteen matter trace form recorded in "
            "the kinetic-form selection receipt"
        ),
    }


def frozen_rg_statistic(decomposition: dict[str, Any]) -> dict[str, Any]:
    return {
        "statistic": "det(alpha_inverse, k, b) = 0",
        "definition": (
            "the three-by-three determinant of the inverse-coupling column, "
            "the kinetic column, and the one-loop beta column vanishes "
            "exactly when the inverse-coupling vector lies in the plane "
            "spanned by the kinetic and beta columns"
        ),
        "owner": (
            "kinetic_form_selection_certificate.py freezes the matter-trace "
            "branch of this statistic with exact columns and cofactors; "
            "this receipt records the same frozen values"
        ),
        "kinetic_column": ["10/3", "2", "2"],
        "kinetic_column_premise": (
            "the rank-fifteen matter-trace branch of the open kinetic-form "
            "selection premise; the port-coordinate rays of this receipt "
            "are not kinetic columns"
        ),
        "beta_column": ["41/6", "-19/6", "-7"],
        "beta_column_premise": (
            "one-loop imported QFT law at the declared (nG, nH) = (3, 1) "
            "completion in the census hypercharge normalization"
        ),
        "exact_cofactors": ["-23/3", "37", "-218/9"],
        "integer_zero_locus": "69 x1 - 333 x2 + 218 x3 = 0",
        "alpha_column": (
            "sealed: measured inverse couplings enter only through the "
            "issue-639 custody surface at its single comparison"
        ),
        "scheme_and_threshold_budget": (
            "one-loop, single-threshold, no extra fields; any change after "
            "scoring voids the statistic rather than repairing it"
        ),
        "frozen_before_comparison": True,
    }


def build_receipt() -> dict[str, Any]:
    bands = load_pinned_bands()
    decomposition = ideal_decomposition(bands)
    rays = candidate_rays(bands, decomposition)
    tests = ray_tests(rays)
    require(
        tests["reference_ray_hit"] is False,
        "reference ray unexpectedly hit; recheck normalization",
    )
    statistic = frozen_rg_statistic(decomposition)
    payload = PORT_RECEIPT_PATH.read_bytes()
    receipt = {
        "schema": SCHEMA,
        "issue": 646,
        "status": STATUS,
        "pinned_band_coefficients": {
            name: q5_str(value) for name, value in bands.items()
        },
        "ideal_decomposition": decomposition,
        "candidate_rays": rays,
        "ray_tests": tests,
        "frozen_rg_statistic": statistic,
        "parent_pins": [
            {
                "path": PORT_RECEIPT_PATH.relative_to(REPO_ROOT).as_posix(),
                "bytes": len(payload),
                "sha256": tagged_sha256(payload),
            }
        ],
        "kinetic_action_bridge": (
            "OPEN: the identity of any finite invariant form with the "
            "physical continuum kinetic action is unproved; per the frozen "
            "stop rule, three independent invariant kinetic coefficients "
            "remain the default freedom, the selection between the "
            "port-response pullback and the matter trace form is a named "
            "open premise, and any ray becomes predictive only after that "
            "bridge is proved or independently tested"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a proved or independently tested kinetic-action bridge and a "
            "resolved kinetic-form selection premise, at which point the "
            "frozen matter-branch determinant statistic becomes the "
            "issue-639 candidate under its declared scheme and threshold "
            "budget"
        ),
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def write_runtime() -> Path:
    RUNTIME.mkdir(exist_ok=True)
    RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    return RECEIPT_PATH


def verify_runtime() -> None:
    if RECEIPT_PATH.read_bytes() != canonical_json_bytes(build_receipt()):
        raise SystemExit("kinetic ray receipt is stale")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        print(write_runtime())
    if args.verify:
        verify_runtime()
        print("KINETIC_RAY_VALID")
    if not args.write and not args.verify:
        receipt = build_receipt()
        print(
            json.dumps(
                {
                    "status": receipt["status"],
                    "reference_ray_hit": receipt["ray_tests"]["reference_ray_hit"],
                    "commutant_hit": receipt["ray_tests"]["commutant_relation_hit"],
                    "su3_average": receipt["ideal_decomposition"][
                        "su3_dimension_weighted_average"
                    ],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
