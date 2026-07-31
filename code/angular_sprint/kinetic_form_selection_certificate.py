#!/usr/bin/env python3
"""Issue #646: exact kinetic-form selection certificate.

The port-current pullback form ``B(f,h) = -Re tr(K(f) K(h))`` on the declared
charged-double-triplet current fixture is ad-invariant: trace cyclicity gives
``B([x,y],z) + B(y,[x,z]) = 0`` once the bracket is transported through the
faithful representation ``K``. This producer verifies that identity exactly on
the full twelve-dimensional structure-constant table. The two distinct
Hilbert--Schmidt band coefficients on the frame and quintet bands therefore do
not witness an invariance failure on the simple ``su(3)`` ideal; they witness
that the port Euclidean metric and the representation trace metric are not
isometrically identified. The earlier ``su(3)``-non-invariance wording and the
promotion of the dimension-weighted band average to a canonical invariant
coefficient are superseded by this receipt; the average is retained as
port-metric arithmetic only.

On each simple ideal every invariant symmetric form is one multiple of the
Killing form, so the canonical content of ``B`` is one Killing-relative
coefficient per simple ideal. This producer computes the Killing Gram matrix
from the exact structure constants, verifies the proportionality on complete
ideal bases, and records the coefficients

* ``c2 = B / (-kappa)`` on the ``su(2)`` ideal,
* ``c3 = B / (-kappa)`` on the ``su(3)`` ideal,

together with the dimensionless ratio ``rho_B = c2/c3``. The rank-fifteen
matter trace form supplies the second candidate kinetic form through the
pinned per-copy representation indices ``(T1, T2, T3) = (10/3, 2, 2)`` with
``c = T / h_dual`` on the simple ideals, hence ``rho_matter = 3/2``. The two
candidate kinetic forms are exactly distinct. Which form, combination, or
family the repair dynamics selects as the physical kinetic action is a named
open source premise; no selection theorem is claimed here.

The frozen matter-branch renormalization-line statistic is recorded with its
exact cofactor coefficients: with kinetic column ``(10/3, 2, 2)`` and one-loop
beta column ``(41/6, -19/6, -7)`` at the declared ``(nG, nH) = (3, 1)``
completion, ``det(x, k, b) = (-23/3) x1 + 37 x2 + (-218/9) x3``, equivalently
``69 x1 - 333 x2 + 218 x3 = 0`` on the zero locus. The inverse-coupling column
``x`` stays sealed with the issue-639 comparison surface. Per-factor coupling
renormalization rescales determinant rows, so the zero locus is
normalization-covariant. The ``u(1)`` coefficient has no Killing
normalization; its scale is typed open pending the primitive-period selection
owned by issue #567.

No public measurement is read and no comparison is opened.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "kinetic_form_selection_receipt.json"
A5_CLOSURE = REPO_ROOT / "code" / "a5_closure"
PORT_RECEIPT_PATH = A5_CLOSURE / "receipts" / "port_current_inner_reference.receipt.json"
MATTER_INDEX_PIN_PATH = REPO_ROOT / "Lean" / "Screen" / "RGRepresentationFrontier.lean"

sys.path.insert(0, str(A5_CLOSURE))

import port_current_inner_certificate as pcc  # noqa: E402

SCHEMA = "oph.kinetic_form_selection_receipt.v1"
STATUS = "EXACT_KINETIC_FORM_DICHOTOMY__AD_INVARIANCE_RESTORED__SELECTION_PREMISE_OPEN"


class KineticSelectionError(ValueError):
    """The kinetic-form selection certificate refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise KineticSelectionError(message)


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


def f5_str(value: Any) -> str:
    return value.text()


def rebuild_current() -> dict[str, Any]:
    """Rebuild the declared current fixture exactly as the pinned receipt does."""

    manifest = pcc.load_json(
        pcc.MODULE_DIR / "manifests" / "port_current_response_reference.json"
    )
    params = pcc.validate_manifest(manifest, pcc.MODULE_DIR)
    carrier, _group_row, plus, carrier_manifest = pcc.load_carrier(
        manifest, pcc.MODULE_DIR
    )
    verts = pcc.standard_vertices()
    matched = pcc.orientation_matched_assignments(carrier, verts)
    require(params["artifact_ref"] is not None, "semantic artifact reference missing")
    artifact = pcc.load_semantic_artifact(params["artifact_ref"], pcc.MODULE_DIR, False)
    binding = pcc.bind_semantic_artifact(
        artifact, carrier, carrier_manifest, verts, matched, plus, params
    )
    frame = pcc.FrameRealization(carrier, binding["psi"], verts)
    model = pcc.MODELS[params["model"]](frame, params)
    generators = [model.generator(field) for field in pcc.BASIS_FIELDS]
    flat = [pcc.flatten(blocks) for blocks in generators]
    require(pcc.rank([row[:] for row in flat]) == 12, "current image rank is not twelve")
    return {
        "frame": frame,
        "model": model,
        "generators": generators,
        "flat": flat,
    }


def hs_gram(generators: list[Any]) -> list[list[Any]]:
    def pairing(x: Any, y: Any) -> Any:
        total = pcc.ZERO
        for bx, by in zip(x, y, strict=True):
            total = total + pcc.ctrace(pcc.cmul(bx, by)).re
        return -total

    gram = [[pairing(generators[p], generators[q]) for q in range(12)] for p in range(12)]
    for p in range(12):
        for q in range(12):
            require(
                (gram[p][q] - gram[q][p]).is_zero(),
                "pullback form is not symmetric",
            )
    return gram


def structure_constants(
    generators: list[Any], flat: list[list[Any]]
) -> list[list[list[Any]]]:
    """Exact bracket table c[i][j] with [e_i, e_j] = sum_k c[i][j][k] e_k."""

    table: list[list[list[Any]]] = [
        [[pcc.ZERO] * 12 for _ in range(12)] for _ in range(12)
    ]
    for i in range(12):
        for j in range(i + 1, 12):
            bracket = tuple(
                pcc.commutator(generators[i][b], generators[j][b])
                for b in range(len(generators[i]))
            )
            coeffs = pcc.solve_in_span(flat, pcc.flatten(bracket))
            table[i][j] = list(coeffs)
            table[j][i] = [-c for c in coeffs]
    return table


def check_ad_invariance(
    gram: list[list[Any]], table: list[list[list[Any]]]
) -> int:
    """Verify B([e_i,e_j],e_k) + B(e_j,[e_i,e_k]) = 0 on all basis triples."""

    checked = 0
    for i in range(12):
        for j in range(12):
            for k in range(12):
                left = pcc.ZERO
                for m in range(12):
                    left = left + table[i][j][m] * gram[m][k]
                right = pcc.ZERO
                for m in range(12):
                    right = right + table[i][k][m] * gram[j][m]
                require(
                    (left + right).is_zero(),
                    f"ad-invariance fails on basis triple ({i},{j},{k})",
                )
                checked += 1
    return checked


def band_field_bases(frame: Any) -> dict[str, list[list[Any]]]:
    """Extract exact port-field bases of the four response bands."""

    projectors = pcc.band_projectors(frame)
    expected = {"unit_band": 1, "kernel_band": 3, "frame_band": 3, "quintet_band": 5}
    bases: dict[str, list[list[Any]]] = {}
    for name, projector in projectors.items():
        columns = [[projector[i][j] for i in range(12)] for j in range(12)]
        basis: list[list[Any]] = []
        for column in columns:
            candidate = basis + [column]
            if pcc.rank([row[:] for row in candidate]) == len(candidate):
                basis.append(column)
        require(
            len(basis) == expected[name],
            f"{name} basis dimension {len(basis)} != {expected[name]}",
        )
        bases[name] = basis
    return bases


def bracket_vector(
    table: list[list[list[Any]]], x: list[Any], y: list[Any]
) -> list[Any]:
    out = [pcc.ZERO] * 12
    for i in range(12):
        if x[i].is_zero():
            continue
        for j in range(12):
            if y[j].is_zero():
                continue
            factor = x[i] * y[j]
            row = table[i][j]
            for k in range(12):
                out[k] = out[k] + factor * row[k]
    return out


def is_zero_vector(vector: list[Any]) -> bool:
    return all(entry.is_zero() for entry in vector)


def in_span(space: list[list[Any]], vector: list[Any]) -> bool:
    if is_zero_vector(vector):
        return True
    stacked = [row[:] for row in space]
    base_rank = pcc.rank([row[:] for row in stacked])
    stacked.append(vector[:])
    return pcc.rank(stacked) == base_rank


def verify_ideal_structure(
    table: list[list[list[Any]]], bases: dict[str, list[list[Any]]]
) -> dict[str, Any]:
    center = bases["unit_band"][0]
    su2 = [column[:] for column in bases["kernel_band"]]
    su3 = [column[:] for column in bases["frame_band"]] + [
        column[:] for column in bases["quintet_band"]
    ]
    require(pcc.rank([row[:] for row in su3]) == 8, "su3 candidate span is not eight-dimensional")
    for basis in (su2, su3):
        for vector in basis:
            require(
                is_zero_vector(bracket_vector(table, center, vector)),
                "center fails to commute with an ideal basis vector",
            )
    for x in su2:
        for y in su3:
            require(
                is_zero_vector(bracket_vector(table, x, y)),
                "su2 and su3 candidate ideals fail to commute",
            )
    for name, basis in (("su2", su2), ("su3", su3)):
        for x in basis:
            for y in basis:
                require(
                    in_span(basis, bracket_vector(table, x, y)),
                    f"{name} candidate span is not bracket-closed",
                )
    return {"center": center, "su2": su2, "su3": su3}


def killing_gram_on(
    table: list[list[list[Any]]], vectors: list[list[Any]]
) -> list[list[Any]]:
    def ad_matrix(x: list[Any]) -> list[list[Any]]:
        matrix = [[pcc.ZERO] * 12 for _ in range(12)]
        for i in range(12):
            if x[i].is_zero():
                continue
            for j in range(12):
                row = table[i][j]
                for k in range(12):
                    matrix[k][j] = matrix[k][j] + x[i] * row[k]
        return matrix

    ads = [ad_matrix(x) for x in vectors]
    size = len(vectors)
    gram = [[pcc.ZERO] * size for _ in range(size)]
    for a in range(size):
        for b in range(size):
            product = pcc.rmul(ads[a], ads[b])
            total = pcc.ZERO
            for d in range(12):
                total = total + product[d][d]
            gram[a][b] = total
    return gram


def form_gram_on(
    gram12: list[list[Any]], vectors: list[list[Any]]
) -> list[list[Any]]:
    size = len(vectors)
    out = [[pcc.ZERO] * size for _ in range(size)]
    for a in range(size):
        for b in range(size):
            total = pcc.ZERO
            for i in range(12):
                if vectors[a][i].is_zero():
                    continue
                for j in range(12):
                    total = total + vectors[a][i] * vectors[b][j] * gram12[i][j]
            out[a][b] = total
    return out


def killing_relative_coefficient(
    b_gram: list[list[Any]], k_gram: list[list[Any]], name: str
) -> Any:
    size = len(b_gram)
    pivot = None
    for d in range(size):
        if not k_gram[d][d].is_zero():
            pivot = d
            break
    require(pivot is not None, f"{name}: Killing Gram has no nonzero diagonal entry")
    coefficient = b_gram[pivot][pivot] / (-k_gram[pivot][pivot])
    for a in range(size):
        for b in range(size):
            expected = coefficient * (-k_gram[a][b])
            require(
                (b_gram[a][b] - expected).is_zero(),
                f"{name}: pullback form is not proportional to the Killing form",
            )
    return coefficient


def port_metric_band_data(gram12: list[list[Any]], frame: Any) -> dict[str, Any]:
    projectors = pcc.band_projectors(frame)
    coefficients: dict[str, Any] = {}
    for name, projector in projectors.items():
        product = pcc.rmul(gram12, projector)
        trace_bp = pcc.ZERO
        trace_p = pcc.ZERO
        for d in range(12):
            trace_bp = trace_bp + product[d][d]
            trace_p = trace_p + projector[d][d]
        coefficients[name] = trace_bp / trace_p
    frame_c = coefficients["frame_band"]
    quintet_c = coefficients["quintet_band"]
    average = (frame_c * pcc.F5(3) + quintet_c * pcc.F5(5)) * pcc.F5(Fraction(1, 8))
    return {
        "band_coefficients": {name: f5_str(value) for name, value in coefficients.items()},
        "su3_dimension_weighted_average": f5_str(average),
        "typing": (
            "eigenvalue data of the Hilbert--Schmidt pullback relative to the "
            "port Euclidean metric: the frame/quintet spread measures the "
            "non-isometry between port coordinates and the representation "
            "trace metric, and the dimension-weighted average is port-metric "
            "arithmetic, not a canonical invariant coefficient"
        ),
    }


def matter_branch() -> dict[str, Any]:
    t1 = Fraction(10, 3)
    t2 = Fraction(2)
    t3 = Fraction(2)
    c2 = t2 / 2
    c3 = t3 / 3
    rho = c2 / c3
    b_y = Fraction(41, 6)
    b_2 = Fraction(-19, 6)
    b_3 = Fraction(-7)
    k = (t1, t2, t3)
    cof1 = k[1] * b_3 - k[2] * b_2
    cof2 = -(k[0] * b_3 - k[2] * b_y)
    cof3 = k[0] * b_2 - k[1] * b_y
    require(cof1 == Fraction(-23, 3), "matter cofactor one drift")
    require(cof2 == Fraction(37), "matter cofactor two drift")
    require(cof3 == Fraction(-218, 9), "matter cofactor three drift")
    scaled = tuple(value * Fraction(-9) for value in (cof1, cof2, cof3))
    require(scaled == (Fraction(69), Fraction(-333), Fraction(218)), "integer form drift")
    return {
        "per_copy_weyl_indices": {"u1": str(t1), "su2": str(t2), "su3": str(t3)},
        "index_pin": (
            "Lean/Screen/RGRepresentationFrontier.lean, theorem "
            "representation_indices, on the registered rank-fifteen census"
        ),
        "killing_relative": {"su2": str(c2), "su3": str(c3)},
        "ratio_su2_over_su3": str(rho),
        "frozen_rg_statistic": {
            "statistic": "det(alpha_inverse, k, b) = 0",
            "kinetic_column_k": [str(v) for v in k],
            "beta_column_b": [str(b_y), str(b_2), str(b_3)],
            "beta_premises": (
                "one-loop imported QFT law at the declared (nG, nH) = (3, 1) "
                "completion in the census hypercharge normalization; "
                "single-threshold, no extra fields; any change after scoring "
                "voids the statistic"
            ),
            "exact_cofactors": [str(cof1), str(cof2), str(cof3)],
            "integer_zero_locus": "69 x1 - 333 x2 + 218 x3 = 0",
            "alpha_column": (
                "sealed: measured inverse couplings enter only through the "
                "issue-639 custody surface at its single comparison"
            ),
            "normalization_covariance": (
                "per-factor coupling renormalization rescales one determinant "
                "row per gauge factor, so the zero locus is invariant under "
                "nonzero per-factor rescaling"
            ),
            "baseline_note": (
                "minimal Standard Model plus general relativity carries no "
                "concurrency constraint on the three renormalization lines: "
                "the three one-loop lines are free to be non-concurrent, so "
                "the vanishing determinant is a genuine frozen discriminator "
                "of the matter-trace branch"
            ),
        },
    }


def build_receipt() -> dict[str, Any]:
    current = rebuild_current()
    gram12 = hs_gram(current["generators"])
    table = structure_constants(current["generators"], current["flat"])
    triples = check_ad_invariance(gram12, table)
    bases = band_field_bases(current["frame"])
    ideals = verify_ideal_structure(table, bases)

    su2_b = form_gram_on(gram12, ideals["su2"])
    su2_k = killing_gram_on(table, ideals["su2"])
    c2 = killing_relative_coefficient(su2_b, su2_k, "su2")

    su3_b = form_gram_on(gram12, ideals["su3"])
    su3_k = killing_gram_on(table, ideals["su3"])
    c3 = killing_relative_coefficient(su3_b, su3_k, "su3")

    rho_port = c2 / c3

    center_b = form_gram_on(gram12, [ideals["center"]])[0][0]
    center_k = killing_gram_on(table, [ideals["center"]])[0][0]
    require(center_k.is_zero(), "center has a nonzero Killing norm")

    matter = matter_branch()
    port_ratio = rho_port
    matter_ratio_fraction = Fraction(3, 2)
    port_is_matter = (
        port_ratio - pcc.F5(matter_ratio_fraction)
    ).is_zero()

    payload = PORT_RECEIPT_PATH.read_bytes()
    matter_pin = MATTER_INDEX_PIN_PATH.read_bytes()
    receipt = {
        "schema": SCHEMA,
        "issue": 646,
        "status": STATUS,
        "current_fixture": (
            "declared charged-double-triplet current on the certified "
            "carrier lineage; the fixture is conditional input, not "
            "source-selected, per the pinned port receipt"
        ),
        "ad_invariance": {
            "identity": "B([x,y],z) + B(y,[x,z]) = 0",
            "verified_basis_triples": triples,
            "reason": (
                "trace cyclicity of the pullback form under the bracket "
                "transported through the faithful representation"
            ),
        },
        "superseded_typing": {
            "withdrawn": (
                "the kinetic_ray_receipt v1 wording that the raw pairing is "
                "not an ad-invariant kinetic form on su(3), and the "
                "promotion of the dimension-weighted band average to the "
                "unique ad-invariant projection"
            ),
            "replacement": (
                "the pairing is ad-invariant; on each simple ideal it is one "
                "Killing-form multiple; the frame/quintet coefficient spread "
                "and the dimension-weighted average are port-metric data of "
                "the embedding, retained below at that typing"
            ),
        },
        "port_metric_band_data": port_metric_band_data(gram12, current["frame"]),
        "killing_relative_coefficients": {
            "su2": f5_str(c2),
            "su3": f5_str(c3),
            "ratio_su2_over_su3": f5_str(rho_port),
            "normalization": (
                "c = B / (-kappa) with kappa the Killing form computed from "
                "the exact structure constants; basis-independent on each "
                "simple ideal"
            ),
        },
        "u1_coefficient": {
            "pullback_on_constant_field": f5_str(center_b),
            "typing": (
                "the abelian ideal has no Killing normalization; the u(1) "
                "kinetic scale is open pending the primitive-period "
                "selection owned by issue #567"
            ),
        },
        "matter_trace_branch": matter,
        "dichotomy": {
            "port_response_ratio_su2_over_su3": f5_str(rho_port),
            "matter_trace_ratio_su2_over_su3": str(matter_ratio_fraction),
            "branches_exactly_distinct": not port_is_matter,
            "selection_premise": (
                "which invariant form the repair dynamics selects as the "
                "physical kinetic action -- the port-response pullback, the "
                "rank-fifteen matter trace, a derived combination, or a "
                "nonunique family -- is a named open source premise; no "
                "selection theorem is claimed"
            ),
        },
        "consumers": {
            "registry": (
                "code/invariant_mining: the kinetic candidates consume this "
                "receipt's typing; the frozen matter-branch determinant is "
                "the branch-conditional issue-639 candidate"
            ),
            "lean": "Lean/Screen/KineticFormDichotomy.lean checks the branch arithmetic",
        },
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "parent_pins": [
            {
                "path": PORT_RECEIPT_PATH.relative_to(REPO_ROOT).as_posix(),
                "bytes": len(payload),
                "sha256": tagged_sha256(payload),
            },
            {
                "path": MATTER_INDEX_PIN_PATH.relative_to(REPO_ROOT).as_posix(),
                "bytes": len(matter_pin),
                "sha256": tagged_sha256(matter_pin),
            },
        ],
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def write_runtime() -> Path:
    RUNTIME.mkdir(exist_ok=True)
    RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    return RECEIPT_PATH


def verify_runtime() -> None:
    if RECEIPT_PATH.read_bytes() != canonical_json_bytes(build_receipt()):
        raise SystemExit("kinetic form selection receipt is stale")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        print(write_runtime())
    if args.verify:
        verify_runtime()
        print("KINETIC_FORM_SELECTION_VALID")
    if not args.write and not args.verify:
        receipt = build_receipt()
        print(
            json.dumps(
                {
                    "status": receipt["status"],
                    "c2": receipt["killing_relative_coefficients"]["su2"],
                    "c3": receipt["killing_relative_coefficients"]["su3"],
                    "port_ratio": receipt["dichotomy"][
                        "port_response_ratio_su2_over_su3"
                    ],
                    "matter_ratio": receipt["dichotomy"][
                        "matter_trace_ratio_su2_over_su3"
                    ],
                    "distinct": receipt["dichotomy"]["branches_exactly_distinct"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
