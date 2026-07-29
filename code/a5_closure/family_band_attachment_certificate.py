"""Family band attachment certificate (issue #569, Lane 2).

Exact selection theorem for the family multiplicity object inside the
source-visible screen coefficient space.  Everything is computed in exact
rational or Q(sqrt5) arithmetic on the pinned twelve-port carrier; no float
enters the payload.

WHAT IS PROVED (exact, on the pinned carrier):

* The complexified screen coefficient space decomposes under the listed
  sixty-rotation action into four isotypic bands 1 + 3 + 3' + 5, exhibited
  by four exact spectral projectors of the port adjacency (eigenvalues
  5, sqrt5, -sqrt5, -1), each idempotent, mutually orthogonal, complete,
  and equivariant under all sixty rotations.  The ordered port pairs fall
  into exactly four orbits and the four distance matrices commute
  pairwise, so the commutant is the four-dimensional commutative
  Bose-Mesner algebra and the module is multiplicity free: every band
  embedding is unique up to scalar.
* Admissibility filters drawn from clauses in the record, namely the
  single-complete-object clause quoted in the #569 body (A3 acts on one
  fixed complete ambient multiplicity object), faithful family exchange
  (the band action has trivial kernel), and the exact physical window
  3 <= N_g <= 5 pinned from the #617 receipt, leave exactly three
  candidates: the 3 band, the 3' band, and the 5 band.
* The operational comparison order of the #625 receipt (quadratic
  seam-mismatch readback, the graph Laplacian 5I - A, per unit norm)
  evaluates on the three candidates to 5 - sqrt5, 5 + sqrt5, and 6.  The
  order is strict, so the comparison selects a unique minimizer: the
  3 band, of complex dimension exactly three.
* All three actual carrier refinements preserve the four exact projectors
  and satisfy their declared cocycle.  Their tensor transports U x I_15
  preserve the rank-forty-five projector P_3 x I_15 exactly.
* The realized attachment object, band tensor generation, has complex rank
  exactly 3 x 15 = 45.  The fifteen-state generation is recomputed from the
  pinned block charges, every three-copy anomaly form is zero, and exhaustive
  replay of the thirty-six center candidates leaves the same diagonal Z6
  kernel after triplication.
* Every character in the pinned #627 Z2, Z3, and Z6 seam menus transports
  uniformly because carrier refinement acts on the family factor and the
  conditional character acts on the matter factor.  This classifies the
  finite conditional menu without selecting a physical seam action.

THE NAMED INTERFACE (this is the boundary; read it before citing):

The selection binds only under `screen_realized_multiplicity_object`,
which carries TWO clauses, each with a control proving it load-bearing:

  (R) realization: the physical pole-residue multiplicity object is a
      single complete subobject of the source-visible screen coefficient
      space (without it, the #617 copy-count invisibility applies and
      nothing is selected; control `external_copy_reduct`);
  (S) selection: the attachment is compared by the #625 operational
      cost order (with the excluded form 6I + A the minimizer flips to
      the Galois partner; control `excluded_cone`).

Neither clause is derived from A1-A3.  The #599 simulator response
artifact realizes clause S on the declared finite channel: the per-band
adjacency channels are executable readbacks, the operational cost evaluated
on them has the frame triplet as strict minimizer, and conjugation swaps the
frame and kernel bands while the readback order separates them.  Clause R is
realized for the response resolvent of the declared Laplacian generator by the
simulator pole-residue artifact.  The propagated dynamics has exactly four pole
clusters at the band costs, the residue at the minimal positive pole is
the rank-three frame projector (faithful, equivariant, Galois partner
at the maximal pole), and the finite screen assembly, simulator-read
multiplicity times the pinned fifteen-state generation, has complex
rank exactly forty-five with the generation factor stated as an
import.  The matter-pole identification, chirality and spin data,
Spin/locality, physical seam selection, and laboratory-current receipts stay
open on issue #569.  The #617 invisibility theorem for external C^n
completions is re-verified and holds unchanged.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import port_current_inner_certificate as p566  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

F5 = p566.F5

SCHEMA = "oph.family_band_attachment_certificate.v6"
MANIFEST_PATH = MODULE_DIR / "manifests" / "family_band_attachment_reference.json"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"
CARRIER_RECEIPT_NAME = "echosahedral_federation_reference.receipt.json"
CURRENT_MANIFEST_NAME = "port_current_response_reference.json"
CURRENT_RECEIPT_NAME = "port_current_inner_reference.receipt.json"
WINDOW_MANIFEST_NAME = "multiplicity_window_reference.json"
READBACK_MANIFEST_NAME = "load_fiber_readback_reference.json"
MATTER_MANIFEST_NAME = "super_tannakian_matter_reference.json"
MATTER_RECEIPT_NAME = "super_tannakian_matter_reference.receipt.json"
GLOBAL_FORM_MANIFEST_NAME = "axis_center_descent_reference.json"
GLOBAL_FORM_RECEIPT_NAME = "axis_center_descent_reference.receipt.json"
SEAM_CLASSIFICATION_NAME = "seam_grammar_matter_classification_reference.json"
RESPONSE_ARTIFACT_NAME = "charged_response_semantic_artifact.json"
POLE_RESIDUE_ARTIFACT_NAME = "charged_response_pole_residue_artifact.json"
MATTER_ATTACHMENT_RECEIPT_NAME = "matter_attachment_receipt.json"

# Measured artifact band labels against this certificate's band names.
ARTIFACT_BAND_MAP = {
    "unit_band": "1",
    "frame_band": "3",
    "kernel_band": "3p",
    "quintet_band": "5",
}

ISSUE = 569
PORTS = 12

SQRT5 = F5(0, 1)
ONE = F5(1, 0)
ZERO = F5(0, 0)


# ---------------------------------------------------------------------------
# Exact F5 matrix helpers on the twelve ports
# ---------------------------------------------------------------------------


Matrix = list  # list[list[F5]], 12 x 12


def mat(fill: F5) -> Matrix:
    return [[fill for _ in range(PORTS)] for _ in range(PORTS)]


def identity() -> Matrix:
    out = mat(ZERO)
    for i in range(PORTS):
        out[i][i] = ONE
    return out


def mat_add(x: Matrix, y: Matrix) -> Matrix:
    return [[x[i][j] + y[i][j] for j in range(PORTS)] for i in range(PORTS)]


def mat_sub(x: Matrix, y: Matrix) -> Matrix:
    return [[x[i][j] - y[i][j] for j in range(PORTS)] for i in range(PORTS)]


def mat_scale(c: F5, x: Matrix) -> Matrix:
    return [[c * x[i][j] for j in range(PORTS)] for i in range(PORTS)]


def mat_mul(x: Matrix, y: Matrix) -> Matrix:
    out = mat(ZERO)
    for i in range(PORTS):
        for k in range(PORTS):
            xik = x[i][k]
            if xik.is_zero():
                continue
            row = y[k]
            oi = out[i]
            for j in range(PORTS):
                oi[j] = oi[j] + xik * row[j]
    return out


def mat_eq(x: Matrix, y: Matrix) -> bool:
    return all(x[i][j] == y[i][j] for i in range(PORTS) for j in range(PORTS))


def mat_trace(x: Matrix) -> F5:
    total = ZERO
    for i in range(PORTS):
        total = total + x[i][i]
    return total


def mat_conj(x: Matrix) -> Matrix:
    return [[x[i][j].conj() for j in range(PORTS)] for i in range(PORTS)]


def mat_is_zero(x: Matrix) -> bool:
    return all(x[i][j].is_zero() for i in range(PORTS) for j in range(PORTS))


def f5_str(value: F5) -> str:
    return f"{value.a}+{value.b}*sqrt5"


def f5_lt(x: F5, y: F5) -> bool:
    return (y - x).is_positive()


def f5_sorted(names: Sequence[str], key: Callable[[str], F5]) -> list[str]:
    """Insertion sort by the exact Q(sqrt5) order (lexicographic pair order
    is NOT the numeric order and must never be used here)."""

    ordered: list[str] = []
    for name in names:
        position = 0
        while position < len(ordered) and f5_lt(key(ordered[position]), key(name)):
            position += 1
        ordered.insert(position, name)
    return ordered


def compose_permutations(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    """Composition `left` after `right` in the carrier convention."""

    require(
        len(left) == len(right),
        "PERMUTATION_SIZE",
        "permutation composition requires equal finite carriers",
    )
    return tuple(int(left[int(right[index])]) for index in range(len(right)))


def tensor_identity_permutation(
    port_permutation: Sequence[int], matter_dimension: int
) -> tuple[int, ...]:
    """The exact basis permutation U tensor I on port-major tensor indices."""

    require(
        sorted(int(value) for value in port_permutation) == list(range(PORTS)),
        "REFINEMENT_BIJECTION",
        "a carrier transport must be a twelve-port permutation",
    )
    require(
        matter_dimension > 0,
        "MATTER_DIMENSION",
        "the tensor identity factor must have positive dimension",
    )
    return tuple(
        int(port_permutation[port]) * matter_dimension + matter
        for port in range(PORTS)
        for matter in range(matter_dimension)
    )


# ---------------------------------------------------------------------------
# Carrier, rotations, and orbit structure
# ---------------------------------------------------------------------------


def load_carrier() -> tuple[Any, list[tuple[int, ...]], dict[str, Any]]:
    manifest = load_json(MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME)
    carrier = e565.validate_carrier(manifest)
    _, plus_group, _ = e565.group_certificate(carrier)
    require(len(plus_group) == 60, "ROTATION_COUNT", "expected sixty listed rotations")
    pin = {
        "path": f"manifests/{CARRIER_MANIFEST_NAME}",
        "sha256": sha256_json(manifest),
        "issue": 565,
    }
    return carrier, [tuple(g) for g in plus_group], pin


def adjacency_int(carrier: Any) -> list[list[int]]:
    return [
        [1 if j in carrier.adjacency[i] else 0 for j in range(PORTS)]
        for i in range(PORTS)
    ]


def distance_matrices(carrier: Any) -> dict[int, list[list[int]]]:
    out: dict[int, list[list[int]]] = {}
    for d in range(4):
        out[d] = [
            [1 if carrier.distances[i][j] == d else 0 for j in range(PORTS)]
            for i in range(PORTS)
        ]
    return out


def lift(x: Sequence[Sequence[int]]) -> Matrix:
    return [[F5(x[i][j], 0) for j in range(PORTS)] for i in range(PORTS)]


def verify_pair_orbits(
    carrier: Any, rotations: Sequence[tuple[int, ...]]
) -> dict[str, Any]:
    """The listed action has exactly four orbits on ordered port pairs, and
    they coincide with the four distance classes.

    This makes the commutant dimension four self-contained here (an
    equivariant matrix is constant on pair orbits), instead of citing the
    Lean module for it.
    """

    seen: set[tuple[int, int]] = set()
    orbits: list[set[tuple[int, int]]] = []
    for i in range(PORTS):
        for j in range(PORTS):
            if (i, j) in seen:
                continue
            orbit = {(g[i], g[j]) for g in rotations}
            require((i, j) in orbit, "ORBIT_IDENTITY", "orbit must contain its seed")
            orbits.append(orbit)
            seen |= orbit
    require(len(orbits) == 4, "PAIR_ORBIT_COUNT", "expected exactly four ordered-pair orbits")
    for orbit in orbits:
        distances = {carrier.distances[i][j] for (i, j) in orbit}
        require(
            len(distances) == 1,
            "ORBIT_DISTANCE_MIXED",
            "a pair orbit must sit inside one distance class",
        )
    return {
        "ordered_pair_orbits": 4,
        "orbits_are_distance_classes": True,
        "consequence": "the commutant of the listed action is the span of the four distance matrices",
    }


# ---------------------------------------------------------------------------
# Spectral resolution of the adjacency over Q(sqrt5)
# ---------------------------------------------------------------------------


def spectral_projectors(adjacency: Matrix) -> dict[str, Matrix]:
    """The four exact eigenprojectors of the port adjacency.

    Eigenvalues: 5 on the invariant line, sqrt5 on the 3 band, -sqrt5 on
    the 3' band, -1 on the 5 band.  Each projector is a cubic polynomial
    in the adjacency divided by the exact eigenvalue-gap product.
    """

    ident = identity()

    def poly(shifts: Sequence[F5]) -> Matrix:
        out = ident
        for shift in shifts:
            out = mat_mul(out, mat_sub(adjacency, mat_scale(shift, ident)))
        return out

    five = F5(5, 0)
    minus_one = F5(-1, 0)

    # Gap products: at 5, (5-sqrt5)(5+sqrt5)(5+1) = 120; at sqrt5,
    # (sqrt5-5)(2 sqrt5)(sqrt5+1) = -40; at -sqrt5 the Galois mirror -40;
    # at -1, (-6)(-1-sqrt5)(-1+sqrt5) = 24.
    p1 = mat_scale(F5(Fraction(1, 120), 0), poly([SQRT5, -SQRT5, minus_one]))
    p3 = mat_scale(F5(Fraction(-1, 40), 0), poly([five, -SQRT5, minus_one]))
    p3p = mat_scale(F5(Fraction(-1, 40), 0), poly([five, SQRT5, minus_one]))
    p5 = mat_scale(F5(Fraction(1, 24), 0), poly([five, SQRT5, -SQRT5]))

    return {"1": p1, "3": p3, "3p": p3p, "5": p5}


def verify_spectral_resolution(
    adjacency: Matrix, projectors: Mapping[str, Matrix]
) -> dict[str, Any]:
    ident = identity()
    names = ["1", "3", "3p", "5"]
    eigen = {"1": F5(5, 0), "3": SQRT5, "3p": -SQRT5, "5": F5(-1, 0)}
    dims = {"1": 1, "3": 3, "3p": 3, "5": 5}

    total = mat(ZERO)
    for name in names:
        p = projectors[name]
        require(mat_eq(mat_mul(p, p), p), "PROJ_IDEMPOTENT", f"P{name} not idempotent")
        require(
            mat_eq(mat_mul(adjacency, p), mat_scale(eigen[name], p)),
            "PROJ_EIGENBAND",
            f"P{name} is not the {f5_str(eigen[name])} eigenband",
        )
        require(
            mat_trace(p) == F5(dims[name], 0),
            "PROJ_TRACE",
            f"trace of P{name} must equal {dims[name]}",
        )
        for other in names:
            if other != name:
                require(
                    mat_is_zero(mat_mul(p, projectors[other])),
                    "PROJ_ORTHOGONAL",
                    f"P{name} P{other} must vanish",
                )
        total = mat_add(total, p)
    require(mat_eq(total, ident), "PROJ_COMPLETE", "projectors must sum to the identity")
    require(
        mat_eq(mat_conj(projectors["3"]), projectors["3p"]),
        "PROJ_GALOIS_PAIR",
        "the Galois conjugate of P3 must equal P3'",
    )
    return {
        "eigenvalues": {name: f5_str(eigen[name]) for name in names},
        "band_dimensions": dims,
        "idempotent_orthogonal_complete": True,
        "galois_pair": "sqrt5 -> -sqrt5 exchanges P3 and P3'",
    }


def verify_equivariance(
    projectors: Mapping[str, Matrix], rotations: Sequence[tuple[int, ...]]
) -> dict[str, Any]:
    for name, p in projectors.items():
        for g in rotations:
            for i in range(PORTS):
                gi = g[i]
                row_g = p[gi]
                row = p[i]
                for j in range(PORTS):
                    if row_g[g[j]] != row[j]:
                        raise CertificateError(
                            "PROJ_EQUIVARIANCE",
                            f"P{name} is not invariant under a listed rotation",
                        )
    return {"rotations_checked": len(rotations), "all_projectors_equivariant": True}


def verify_multiplicity_free(carrier: Any) -> dict[str, Any]:
    """The four distance matrices commute pairwise (Bose-Mesner algebra).

    Together with the four-orbit count this makes the commutant a
    four-dimensional commutative algebra, so the permutation module is
    multiplicity free and each band embedding is unique up to scalar.
    """

    dist = distance_matrices(carrier)
    lifted = {d: lift(m) for d, m in dist.items()}
    for a in range(4):
        for b in range(a + 1, 4):
            require(
                mat_eq(
                    mat_mul(lifted[a], lifted[b]), mat_mul(lifted[b], lifted[a])
                ),
                "COMMUTANT_NOT_COMMUTATIVE",
                f"distance matrices {a} and {b} must commute",
            )
    return {
        "distance_matrices": 4,
        "pairwise_commuting": True,
        "consequence": "multiplicity-free permutation module; every band embedding is unique up to scalar",
        "lean_companion": "Screen/A5Commutant.lean carries the same commutant independently",
    }


# ---------------------------------------------------------------------------
# Pinned upstream receipts
# ---------------------------------------------------------------------------


def _pin_named(
    directory: str, name: str, issue: int, role: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    path = MODULE_DIR / directory / name
    payload = load_json(path)
    return payload, {
        "path": f"{directory}/{name}",
        "sha256": sha256_json(payload),
        "issue": issue,
        "role": role,
    }


def _verify_self_hash(payload: Mapping[str, Any], code: str) -> None:
    body = {key: value for key, value in payload.items() if key != "manifest_sha256"}
    require(
        payload.get("manifest_sha256") == "sha256:" + sha256_json(body),
        code,
        "the generated parent manifest does not match its self-hash",
    )


def pin_structural_chain(
    carrier_manifest_sha256: str,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """Resolve #565/#566/#314/#567/#627 as one exact finite chain."""

    carrier_receipt, carrier_receipt_pin = _pin_named(
        "receipts",
        CARRIER_RECEIPT_NAME,
        565,
        "carrier_refinement_receipt",
    )
    current_manifest, current_manifest_pin = _pin_named(
        "manifests",
        CURRENT_MANIFEST_NAME,
        566,
        "current_manifest",
    )
    current_receipt, current_receipt_pin = _pin_named(
        "receipts",
        CURRENT_RECEIPT_NAME,
        566,
        "current_refinement_receipt",
    )
    matter_manifest, matter_manifest_pin = _pin_named(
        "manifests",
        MATTER_MANIFEST_NAME,
        314,
        "matter_manifest",
    )
    matter_receipt, matter_receipt_pin = _pin_named(
        "receipts",
        MATTER_RECEIPT_NAME,
        314,
        "matter_refinement_receipt",
    )
    global_manifest, global_manifest_pin = _pin_named(
        "manifests",
        GLOBAL_FORM_MANIFEST_NAME,
        567,
        "global_form_manifest",
    )
    global_receipt, global_receipt_pin = _pin_named(
        "receipts",
        GLOBAL_FORM_RECEIPT_NAME,
        567,
        "global_form_receipt",
    )
    seam_manifest, seam_manifest_pin = _pin_named(
        "manifests",
        SEAM_CLASSIFICATION_NAME,
        627,
        "conditional_seam_menu",
    )

    require(
        carrier_receipt.get("schema")
        == "oph.echosahedral_selector_receipt.v1"
        and carrier_receipt.get("issue") == 565
        and carrier_receipt.get("manifest_sha256")
        == carrier_manifest_sha256,
        "CARRIER_RECEIPT_CHAIN",
        "the #565 receipt is not bound to the selected carrier manifest",
    )
    require(
        current_manifest.get("schema") == "oph.port_current_response_manifest.v5"
        and current_manifest.get("carrier_manifest_sha256")
        == carrier_manifest_sha256,
        "CURRENT_MANIFEST_CHAIN",
        "the #566 current manifest is not bound to the selected carrier",
    )
    require(
        current_receipt.get("schema") == "oph.port_current_inner_receipt.v5"
        and current_receipt.get("issue") == 566
        and current_receipt.get("manifest_sha256")
        == sha256_json(current_manifest)
        and current_receipt.get("carrier_manifest_sha256")
        == carrier_manifest_sha256,
        "CURRENT_RECEIPT_CHAIN",
        "the #566 receipt does not resolve to its manifest and carrier",
    )

    expected_maps = [
        {"source": "r0", "target": "r1", "intertwined": True},
        {"source": "r1", "target": "r2", "intertwined": True},
        {"source": "r0", "target": "r2", "intertwined": True},
    ]
    expected_physical_maps = [
        {
            "source_level": 0,
            "target_level": 1,
            "origin": "defect_port_persistence_in_geodesic_tower",
            "intertwined": True,
        },
        {
            "source_level": 0,
            "target_level": 2,
            "origin": "defect_port_persistence_in_geodesic_tower",
            "intertwined": True,
        },
    ]
    current_refinement = current_receipt.get("refinement", {})
    require(
        current_refinement.get("natural") is True
        and current_refinement.get("naturality") == expected_maps
        and current_refinement.get("physical_naturality")
        == expected_physical_maps
        and current_refinement.get("carrier_tower", {}).get(
            "checked_cocycle_triangles"
        )
        == 1,
        "CURRENT_REFINEMENT_CHAIN",
        "the #566 current packet lost exact three-map refinement naturality",
    )

    require(
        matter_manifest.get("schema")
        == "oph.super_tannakian_matter_manifest.v5"
        and matter_manifest.get("current_manifest_sha256")
        == sha256_json(current_manifest)
        and matter_manifest.get("current_receipt_sha256")
        == sha256_json(current_receipt),
        "MATTER_MANIFEST_CHAIN",
        "the #314 matter manifest is not bound to the #566 current packet",
    )
    require(
        matter_receipt.get("schema")
        == "oph.super_tannakian_matter_receipt.v5"
        and matter_receipt.get("issue") == 314
        and matter_receipt.get("manifest_sha256")
        == sha256_json(matter_manifest)
        and matter_receipt.get("refinement", {}).get("natural") is True
        and matter_receipt.get("refinement", {}).get("maps") == expected_maps
        and matter_receipt.get("refinement", {}).get("physical_maps")
        == expected_physical_maps,
        "MATTER_RECEIPT_CHAIN",
        "the #314 matter receipt lost its exact #566 refinement binding",
    )

    require(
        global_manifest.get("schema") == "oph.axis_center_descent_manifest.v4"
        and global_manifest.get("matter_receipt_sha256")
        == sha256_json(matter_receipt),
        "GLOBAL_FORM_MANIFEST_CHAIN",
        "the #567 global-form manifest is not bound to the matter receipt",
    )
    kernel = global_receipt.get("kernel_on_realized_tensors", {})
    expected_kernel = [
        [0, 0, 0],
        [0, 1, 3],
        [1, 0, 4],
        [1, 1, 1],
        [2, 0, 2],
        [2, 1, 5],
    ]
    require(
        global_receipt.get("schema") == "oph.axis_center_descent_receipt.v4"
        and global_receipt.get("issue") == 567
        and global_receipt.get("manifest_sha256")
        == sha256_json(global_manifest)
        and global_receipt.get("matter_receipt_sha256")
        == sha256_json(matter_receipt)
        and kernel.get("candidates_enumerated") == 36
        and kernel.get("kernel_order") == 6
        and kernel.get("cyclic_generator") == [1, 1, 1]
        and kernel.get("kernel_elements") == expected_kernel
        and kernel.get("matches_emitted_kernel_data") is True
        and global_receipt.get("weight_level_refinement_invariance", {}).get(
            "physical_loop_or_bundle_refinement_naturality_derived"
        )
        is True
        and global_receipt.get("sector_transport_consistency", {}).get(
            "unique_menu_matching_form"
        )
        == "z6_quotient",
        "GLOBAL_FORM_RECEIPT_CHAIN",
        "the #567 receipt does not resolve to the matter packet and Z6 kernel",
    )

    _verify_self_hash(seam_manifest, "SEAM_MANIFEST_HASH")
    seam_pins = seam_manifest.get("upstream_pins", {})
    require(
        seam_manifest.get("schema")
        == "oph.seam_grammar_matter_classification_certificate.v2"
        and seam_manifest.get("issue") == 627
        and seam_pins.get("matter_receipt", {}).get("sha256")
        == sha256_json(matter_receipt)
        and seam_pins.get("diagonal_global_form_receipt", {}).get("sha256")
        == sha256_json(global_receipt),
        "SEAM_MANIFEST_CHAIN",
        "the #627 seam menu is not bound to #314 and #567",
    )
    interface = seam_manifest.get("matter_action_interface", {})
    diagonal = seam_manifest.get("diagonal_kernel_action", {})
    character_groups = seam_manifest.get(
        "hypercharge_character_menu", {}
    ).get("groups", [])
    require(
        interface.get("class") == "conditional_open_interface"
        and interface.get("owner_issue") == 569
        and "selection requires" in interface.get("statement", ""),
        "SEAM_SELECTION_BOUNDARY",
        "the #627 menu must leave physical seam action selection open on #569",
    )
    require(
        diagonal.get("generator_color_weak_hypercharge") == [1, 1, 1]
        and diagonal.get("module_dimension") == 15
        and diagonal.get("fixed_subspace_dimension") == 15
        and all(
            row.get("phase_sixths") == 0
            for row in diagonal.get("fields", [])
        )
        and len(diagonal.get("fields", [])) == 5
        and [row.get("group_order") for row in character_groups]
        == [2, 3, 6]
        and [len(row.get("characters", [])) for row in character_groups]
        == [2, 3, 6],
        "SEAM_MENU_CONTENT",
        "the #627 diagonal action or complete character menu has drifted",
    )

    payloads = {
        "current_manifest": current_manifest,
        "current_receipt": current_receipt,
        "matter_manifest": matter_manifest,
        "matter_receipt": matter_receipt,
        "global_form_manifest": global_manifest,
        "global_form_receipt": global_receipt,
        "seam_manifest": seam_manifest,
    }
    pins = {
        "carrier_refinement_receipt": carrier_receipt_pin,
        "current_manifest": current_manifest_pin,
        "current_receipt": current_receipt_pin,
        "matter_manifest": matter_manifest_pin,
        "matter_receipt": matter_receipt_pin,
        "global_form_manifest": global_manifest_pin,
        "global_form_receipt": global_receipt_pin,
        "seam_classification": seam_manifest_pin,
    }
    return payloads, pins


def pin_window() -> tuple[int, int, dict[str, Any]]:
    manifest = load_json(MODULE_DIR / "manifests" / WINDOW_MANIFEST_NAME)
    window = manifest["family_multiplicity_window"]
    lower = int(window["cp_capability_lower_edge"]["lower_edge"])
    upper = int(window["su2_ultraviolet_upper_edge"]["upper_edge"])
    require((lower, upper) == (3, 5), "WINDOW_MISMATCH", "expected the exact window [3, 5]")
    require(
        window["in_window_non_selection"]["count_inside_window_source_selected"] is False,
        "WINDOW_SELECTION_DRIFT",
        "the pinned receipt must record in-window non-selection",
    )
    pin = {
        "path": f"manifests/{WINDOW_MANIFEST_NAME}",
        "sha256": sha256_json(manifest),
        "issue": 617,
    }
    return lower, upper, pin


def pin_cost_cone() -> tuple[tuple[int, int], tuple[int, int], dict[str, Any]]:
    manifest = load_json(MODULE_DIR / "manifests" / READBACK_MANIFEST_NAME)
    cone = manifest["operational_cost_cone"]
    grammar = cone["a2_comparison_grammar"]
    require(
        "5I - A" in str(grammar["seam_translation_access"]),
        "CONE_LAPLACIAN_MISSING",
        "the pinned cone must generate the seam Laplacian 5I - A",
    )
    excluded = cone["candidate_6I_plus_A"]
    require(
        excluded["classification"] == "excluded_from_operational_comparison_cone",
        "CONE_EXCLUSION_DRIFT",
        "the pinned receipt must exclude 6I + A",
    )
    pin = {
        "path": f"manifests/{READBACK_MANIFEST_NAME}",
        "sha256": sha256_json(manifest),
        "issue": 625,
    }
    return (5, -1), (6, 1), pin


def pin_matter() -> tuple[dict[str, Fraction], dict[str, Any]]:
    manifest = load_json(MODULE_DIR / "manifests" / MATTER_MANIFEST_NAME)
    charges_raw = manifest["exterior_matter_contract"]["block_trace_charges"]
    charges = {
        "color_block": Fraction(str(charges_raw["color_block"])),
        "weak_block": Fraction(str(charges_raw["weak_block"])),
    }
    require(
        charges == {"color_block": Fraction(-1, 3), "weak_block": Fraction(1, 2)},
        "MATTER_CHARGES_DRIFT",
        "the pinned block charges must be (-1/3, 1/2)",
    )
    pin = {
        "path": f"manifests/{MATTER_MANIFEST_NAME}",
        "sha256": sha256_json(manifest),
        "issue": 314,
    }
    return charges, pin


# ---------------------------------------------------------------------------
# The measured response artifact (issue #599) and the clause receipts
# ---------------------------------------------------------------------------


def parse_channel(text: str) -> F5:
    """Parse an exact channel string of the artifact.

    Accepted forms: an integer or fraction ('5', '-1'), or
    'a + b*sqrt(5)' with integer or fraction parts.
    """

    raw = str(text).strip()
    if "sqrt" not in raw:
        return F5(Fraction(raw), 0)
    left, right = raw.split("+")
    rational = Fraction(left.strip())
    coeff = Fraction(right.strip().split("*")[0].strip())
    return F5(rational, coeff)


def pin_response_artifact(carrier_pin: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    artifact = load_json(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    require(
        artifact["carrier_binding"]["carrier_manifest_sha256"]
        == carrier_pin["sha256"],
        "ARTIFACT_CARRIER_MISMATCH",
        "the measured artifact must bind the same carrier manifest",
    )
    pin = {
        "path": f"manifests/{RESPONSE_ARTIFACT_NAME}",
        "sha256": sha256_json(artifact),
        "issue": 599,
    }
    return artifact, pin


def pin_pole_residue_artifact(
    carrier_pin: Mapping[str, Any], response_artifact: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    artifact = load_json(MODULE_DIR / "manifests" / POLE_RESIDUE_ARTIFACT_NAME)
    binding = artifact["carrier_binding"]
    require(
        binding["carrier_manifest_sha256"] == carrier_pin["sha256"],
        "POLE_ARTIFACT_CARRIER_MISMATCH",
        "the pole-residue artifact must bind the same carrier manifest",
    )
    require(
        binding["parent_artifact_sha256"] == response_artifact["artifact_sha256"],
        "POLE_ARTIFACT_CHAIN",
        "the pole-residue artifact must pin the response artifact it extends",
    )
    pin = {
        "path": f"manifests/{POLE_RESIDUE_ARTIFACT_NAME}",
        "sha256": sha256_json(artifact),
        "issue": 569,
    }
    return artifact, pin


def pin_matter_attachment_receipt(
    response_artifact: Mapping[str, Any],
    pole_artifact: Mapping[str, Any],
    receipt_override: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Pin the finite Spin/locality attachment receipt of the local domain.

    The receipt is produced on the issue-634 typed domain and must bind
    the same response and pole-residue artifacts this certificate pins,
    carry the attained verdict with every clause true, and record the
    rank forty-five attachment with the measured rank-three band.  The
    override parameter exists for the mutation control only."""

    receipt = (
        json.loads(json.dumps(dict(receipt_override)))
        if receipt_override is not None
        else load_json(MODULE_DIR / "manifests" / MATTER_ATTACHMENT_RECEIPT_NAME)
    )
    require(
        receipt["schema"] == "oph.local-domain-matter-attachment.v1"
        and receipt["issue"] == 569
        and receipt["physical_promotion_allowed"] is False,
        "MATTER_ATTACHMENT_SCHEMA",
        "the attachment receipt must carry the declared schema and issue",
    )
    pins = receipt["upstream_pins"]
    require(
        pins["response_artifact_sha256"] == response_artifact["artifact_sha256"]
        and pins["pole_residue_artifact_sha256"]
        == pole_artifact["artifact_sha256"],
        "MATTER_ATTACHMENT_CHAIN",
        "the attachment receipt must bind the pinned response and "
        "pole-residue artifacts",
    )
    require(
        receipt["attachment"]["complex_rank"] == 45
        and receipt["attachment"]["band_rank_measured"] == 3,
        "MATTER_ATTACHMENT_RANK",
        "the attachment receipt must carry the rank forty-five object "
        "from the measured rank-three band",
    )
    require(
        receipt["verdict"] == "ATTAINED"
        and bool(receipt["clause_verdicts"])
        and all(receipt["clause_verdicts"].values())
        and receipt["blockers"] == [],
        "MATTER_ATTACHMENT_VERDICT",
        "the attachment receipt must be attained with every clause true",
    )
    pin = {
        "path": f"manifests/{MATTER_ATTACHMENT_RECEIPT_NAME}",
        "sha256": sha256_json(receipt),
        "issue": 569,
        "scope": "finite Spin/locality layer on the issue-634 local domain",
    }
    return receipt, pin


def pole_residue_receipt(
    pole_artifact: Mapping[str, Any], generation: Mapping[str, Any]
) -> dict[str, Any]:
    """The simulator pole-residue realization of the multiplicity object.

    The simulator artifact carries the recorded pole clusters and exact
    reconstructed residues of the response resolvent of the propagated
    dynamics. This receipt verifies, against the pinned data, that the four
    recorded poles are the Laplacian band costs with multiplicities
    1, 3, 3, 5; the residue at the minimal positive pole is the rank-three
    frame projector, faithful, equivariant, with its Galois partner at the
    maximal pole; and the finite screen assembly, simulator-read multiplicity
    times the pinned fifteen-state generation, has complex rank exactly
    forty-five.  The generation factor is the pinned exact packet, stated
    as an import. The simulator reads the multiplicity factor from its
    declared generator; this is not a laboratory measurement.
    """

    readback = pole_artifact["pole_residue_readback"]
    poles = readback["measured_poles"]
    expected = {
        "unit": (F5(0, 0), 1),
        "frame": (F5(5, -1), 3),
        "quintet": (F5(6, 0), 5),
        "kernel": (F5(5, 1), 3),
    }
    for band, (value, multiplicity) in expected.items():
        require(
            parse_channel(poles[band]["pole"]) == value,
            "POLE_TABLE",
            f"the simulator-read {band} pole must reconstruct to its band cost",
        )
        require(
            int(poles[band]["multiplicity"]) == multiplicity,
            "POLE_TABLE",
            f"the simulator-read {band} pole multiplicity must equal {multiplicity}",
        )
    residue = readback["family_band_residue"]
    require(
        residue
        == {
            "band": "frame",
            "measured_rank": 3,
            "equals_exact_frame_projector": True,
            "lowest_positive_generator_frequency": True,
            "unitary_mode_norms_conserved": True,
            "faithful_kernel_order": 1,
            "equivariant_under_all_automorphisms": True,
            "galois_partner_at_maximal_pole": True,
        },
        "FAMILY_BAND_RESIDUE",
        "the simulator family-band residue receipt does not match its required form",
    )
    rank = int(residue["measured_rank"]) * int(generation["weyl_state_count"])
    require(rank == 45, "POLE_RANK_45", "the realized attachment rank must be forty-five")
    return {
        "measured_poles": {
            band: dict(poles[band]) for band in ("unit", "frame", "quintet", "kernel")
        },
        "family_band_residue": dict(residue),
        "attachment_rank": {
            "measured_multiplicity_rank": int(residue["measured_rank"]),
            "generation_weyl_states_imported": int(generation["weyl_state_count"]),
            "complex_rank": rank,
        },
        "realization_scope": (
            "the pole-residue object of the response resolvent of the "
            "declared unitary screen Laplacian, realized inside the screen "
            "coefficient space; this finite simulator receipt carries no "
            "relaxation rate, matter-pole identification, chirality or spin "
            "data, or laboratory attachment"
        ),
    }


def control_matter_attachment_mutation(
    response_artifact: Mapping[str, Any],
    pole_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    """A doctored attachment rank must be refused by the receipt pin."""

    doctored = json.loads(
        json.dumps(
            load_json(MODULE_DIR / "manifests" / MATTER_ATTACHMENT_RECEIPT_NAME)
        )
    )
    doctored["attachment"]["complex_rank"] = 60
    try:
        pin_matter_attachment_receipt(
            response_artifact, pole_artifact, receipt_override=doctored
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "MATTER_ATTACHMENT_RANK",
            "meaning": (
                "the pin path reads the receipt rank; a doctored rank is "
                "refused by the same loader the certificate uses"
            ),
        }
    return {"expected_failure": True, "failed": False}


def control_pole_table_mutation(pole_artifact: Mapping[str, Any]) -> dict[str, Any]:
    """A doctored pole table (frame pole moved onto the quintet cost) must be
    refused by the pole gate."""

    doctored = json.loads(json.dumps(dict(pole_artifact)))
    doctored["pole_residue_readback"]["measured_poles"]["frame"]["pole"] = "6"
    generation_stub = {"weyl_state_count": 15}
    try:
        pole_residue_receipt(doctored, generation_stub)
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "POLE_TABLE",
            "meaning": "the receipt reads the simulator pole values; a moved pole is refused",
        }
    return {"expected_failure": True, "failed": False}


def measured_band_receipt(
    carrier: Any,
    adjacency: Matrix,
    projectors: Mapping[str, Matrix],
    artifact: Mapping[str, Any],
    operational: tuple[int, int],
    kernels: Mapping[str, int],
) -> dict[str, Any]:
    """Clause receipts read from the #599 simulator response artifact.

    The artifact carries, as exact data: the per-band adjacency channel
    values, the sector dimensions, the Galois pairing, the response band
    scales, and the antipode polynomial of the recorded response
    operator. This receipt binds those simulator readbacks to the selection:

    * the simulator channel values reproduce the exact band spectrum, so
      the operational cost evaluated on those channels gives the band costs
      and their strict order (clause S realized on the declared simulator);
    * the simulator-read frame band is the strict minimizer among the faithful
      bands, and it is a subobject of the screen coefficient space
      exhibited by the recorded response basis;
    * the recorded antipode polynomial (A^3 - 4A^2 - 5A + 10I)/10 equals
      the carrier antipode permutation exactly, and the real part of the
      scaled band projector equals ten times identity-minus-antipode, so
      the simulator response algebra and the spectral selection algebra
      are one object;
    * conjugation swaps the recorded frame and kernel bands while the
      recorded cost order separates them, so the Galois resolution is
      explicit in the simulator artifact.
    """

    basis = artifact["response_basis"]
    channels = {
        name: parse_channel(value)
        for name, value in basis["adjacency_channel_values"].items()
    }
    dims = {name: int(value) for name, value in basis["sector_dimensions"].items()}
    require(
        set(channels) == set(ARTIFACT_BAND_MAP) and set(dims) == set(ARTIFACT_BAND_MAP),
        "ARTIFACT_BANDS",
        "the artifact must carry exactly the four named bands",
    )

    eigen = {"1": F5(5, 0), "3": SQRT5, "3p": -SQRT5, "5": F5(-1, 0)}
    for artifact_name, band in ARTIFACT_BAND_MAP.items():
        require(
            channels[artifact_name] == eigen[band],
            "MEASURED_CHANNEL_MISMATCH",
            f"measured channel of {artifact_name} must equal the {band} band eigenvalue",
        )
        require(
            dims[artifact_name] == BAND_DIMS[band],
            "MEASURED_DIMENSION_MISMATCH",
            f"measured dimension of {artifact_name} must equal {BAND_DIMS[band]}",
        )

    a, b = operational
    measured_costs = {
        ARTIFACT_BAND_MAP[name]: F5(a, 0) + F5(b, 0) * value
        for name, value in channels.items()
    }
    faithful = [band for band, count in kernels.items() if count == 1]
    ordered = f5_sorted(faithful, lambda name: measured_costs[name])
    for left, right in zip(ordered, ordered[1:]):
        require(
            f5_lt(measured_costs[left], measured_costs[right]),
            "MEASURED_COST_ORDER",
            "the measured band costs must be strictly ordered",
        )
    require(
        ordered[0] == "3" and ARTIFACT_BAND_MAP["frame_band"] == "3",
        "MEASURED_MINIMIZER",
        "the measured frame band must be the strict cost minimizer",
    )

    pairing = basis["galois_pairing"]
    require(
        pairing["frame_and_kernel_swapped_by_conjugation"] is True
        and pairing["unit_and_quintet_galois_stable"] is True,
        "MEASURED_GALOIS_PAIRING",
        "the artifact must record the measured Galois pairing",
    )

    scales = artifact["derived"]["response_band_scales"]
    require(
        scales == {"frame_band": "1", "kernel_band": "1", "quintet_band": "-1", "unit_band": "-1"},
        "MEASURED_RESPONSE_SCALES",
        "the measured response must scale the double triplet by one and the complement by minus one",
    )

    # Bind the measured antipode polynomial to the carrier and to the
    # spectral selection: 10*antipode = A^3 - 4A^2 - 5A + 10I, and the
    # real part of the scaled 3-band projector is 10*(I - antipode).
    ident = identity()
    a2 = mat_mul(adjacency, adjacency)
    a3 = mat_mul(a2, adjacency)
    poly = mat_add(
        mat_sub(a3, mat_scale(F5(4, 0), a2)),
        mat_add(mat_scale(F5(-5, 0), adjacency), mat_scale(F5(10, 0), ident)),
    )
    antipode = mat(ZERO)
    for i in range(PORTS):
        antipode[i][carrier.antipode[i]] = ONE
    require(
        mat_eq(poly, mat_scale(F5(10, 0), antipode)),
        "ANTIPODE_POLYNOMIAL",
        "the measured antipode polynomial must equal the carrier antipode",
    )
    x_real = mat_scale(F5(10, 0), mat_sub(ident, antipode))
    p3_scaled = mat_scale(F5(40, 0), projectors["3"])
    require(
        all(
            p3_scaled[i][j].a == x_real[i][j].a and x_real[i][j].b == 0
            for i in range(PORTS)
            for j in range(PORTS)
        ),
        "BAND_RESPONSE_IDENTITY",
        "the real part of the scaled 3-band projector must equal ten times identity minus antipode",
    )

    return {
        "artifact_issue": 599,
        "measured_channels": {name: f5_str(value) for name, value in channels.items()},
        "measured_band_costs": {name: f5_str(value) for name, value in measured_costs.items()},
        "measured_cost_order": [
            {"object": name, "cost": f5_str(measured_costs[name])} for name in ordered
        ],
        "measured_minimizer": "frame_band (the 3 band)",
        "galois_pairing_measured": True,
        "response_double_triplet_scales": scales,
        "antipode_polynomial_bound": "10*antipode = A^3 - 4A^2 - 5A + 10I on the pinned carrier",
        "band_response_identity": "Re(40 P3) = 10 (I - antipode)",
        "clause_S": "simulator_realized_on_declared_channel",
        "clause_R": (
            "simulator_realized_for_declared_response_resolvent__"
            "matter_pole_identification_open"
        ),
        "open_receipt": (
            "the physical matter-pole identification of the finite rank-forty-five "
            "screen assembly; its gates are specified on issue #569"
        ),
    }


def control_measured_channel_swap(
    artifact: Mapping[str, Any], operational: tuple[int, int], kernels: Mapping[str, int]
) -> dict[str, Any]:
    """Swapping the measured frame and kernel channel values must flip the
    measured minimizer to the Galois partner, so the receipt consumes the
    measured values and not the band labels."""

    basis = artifact["response_basis"]
    swapped = dict(basis["adjacency_channel_values"])
    swapped["frame_band"], swapped["kernel_band"] = (
        swapped["kernel_band"],
        swapped["frame_band"],
    )
    a, b = operational
    costs = {
        ARTIFACT_BAND_MAP[name]: F5(a, 0) + F5(b, 0) * parse_channel(value)
        for name, value in swapped.items()
    }
    faithful = [band for band, count in kernels.items() if count == 1]
    ordered = f5_sorted(faithful, lambda name: costs[name])
    flipped = ordered[0]
    try:
        require(
            flipped == "3",
            "MEASURED_SWAP_DETECTED",
            "the swapped channel table selects the wrong band",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "MEASURED_SWAP_DETECTED",
            "swapped_minimizer": flipped,
            "meaning": "the measured receipt reads the channel values, so a frame/kernel value swap is detected as the Galois partner",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# The fifteen-state generation from the pinned block charges
# ---------------------------------------------------------------------------


def generation_certificate(charges: Mapping[str, Fraction]) -> dict[str, Any]:
    """Branch Lambda^2 V + Lambda^4 V for V = C + W, dims 3 + 2, exactly.

    Lambda^2 gives u_c (Lambda^2 C), Q (C x W), e_c (Lambda^2 W); Lambda^4
    of the five-dimensional V is det(V) tensor V*, giving d_c and L with
    hypercharges det - y_block.  Anomaly forms and the doublet count are
    recomputed from these rows; nothing is imported as a number.
    """

    yc = charges["color_block"]
    yw = charges["weak_block"]
    det = 3 * yc + 2 * yw
    states = [
        {"label": "u_c", "color": 3, "weak": 1, "y": 2 * yc, "states": 3},
        {"label": "Q", "color": 3, "weak": 2, "y": yc + yw, "states": 6},
        {"label": "e_c", "color": 1, "weak": 1, "y": 2 * yw, "states": 1},
        {"label": "d_c", "color": 3, "weak": 1, "y": det - yc, "states": 3},
        {"label": "L", "color": 1, "weak": 2, "y": det - yw, "states": 2},
    ]
    total = sum(int(s["states"]) for s in states)
    require(total == 15, "GENERATION_COUNT", "the exterior branch must carry fifteen Weyl states")
    for s in states:
        require(
            int(s["states"]) == int(s["color"]) * int(s["weak"]),
            "GENERATION_STATE_COUNT",
            "each row's Weyl count must equal color times weak dimension",
        )

    u1_cubed = Fraction(0)
    su3_sq_u1 = Fraction(0)
    su2_sq_u1 = Fraction(0)
    grav_u1 = Fraction(0)
    for s in states:
        y = Fraction(s["y"])
        color = int(s["color"])
        weak = int(s["weak"])
        count = color * weak
        u1_cubed += count * y**3
        grav_u1 += count * y
        if color == 3:
            su3_sq_u1 += weak * y
        if weak == 2:
            su2_sq_u1 += color * y
    forms = {
        "u1_cubed": u1_cubed,
        "su3_sq_u1": su3_sq_u1,
        "su2_sq_u1": su2_sq_u1,
        "grav_u1": grav_u1,
    }
    require(
        all(value == 0 for value in forms.values()),
        "GENERATION_ANOMALY",
        "the fifteen-state generation must cancel every listed anomaly form",
    )
    doublets = sum(int(s["color"]) for s in states if int(s["weak"]) == 2)
    require(doublets == 4, "GENERATION_DOUBLETS", "one generation must carry four weak doublets")
    return {
        "states": [
            {
                "label": s["label"],
                "color": s["color"],
                "weak": s["weak"],
                "hypercharge": str(Fraction(s["y"])),
                "weyl_states": s["states"],
            }
            for s in states
        ],
        "weyl_state_count": total,
        "per_family_anomaly_forms": {k: str(v) for k, v in forms.items()},
        "weak_doublets_per_family": doublets,
    }


# ---------------------------------------------------------------------------
# Exact refinement, rank-45 tensor, anomaly, and seam-menu transport
# ---------------------------------------------------------------------------


def projector_fixed_by(
    projector: Matrix, permutation: Sequence[int]
) -> bool:
    return all(
        projector[int(permutation[i])][int(permutation[j])] == projector[i][j]
        for i in range(PORTS)
        for j in range(PORTS)
    )


def verify_refinement_transport(
    carrier: Any,
    rotations: Sequence[tuple[int, ...]],
    projectors: Mapping[str, Matrix],
    matter_dimension: int,
    matter_receipt: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[tuple[str, str], tuple[int, ...]]]:
    """Replay all three #565 maps on P_band and P3 tensor I_15."""

    carrier_manifest = load_json(
        MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME
    )
    refinement_summary = e565.validate_refinement(
        carrier_manifest,
        carrier,
        rotations,
        e565.gram_matrix(carrier),
    )
    require(
        refinement_summary.get("declared_nonidentity_maps") == 3
        and refinement_summary.get("checked_cocycle_triangles") == 1,
        "REFINEMENT_TOWER",
        "the carrier must expose the three-map refinement triangle",
    )

    rotation_set = set(rotations)
    maps: dict[tuple[str, str], tuple[int, ...]] = {}
    rows = []
    tensor_dimension = PORTS * matter_dimension
    for item in carrier_manifest["refinement_tower"]["maps"]:
        source = str(item["source"])
        target = str(item["target"])
        permutation = e565.parse_port_permutation(item["port_map"], carrier)
        require(
            permutation in rotation_set,
            "REFINEMENT_ROTATION",
            "each actual refinement must lie in the listed order-sixty action",
        )
        require(
            all(
                projector_fixed_by(projector, permutation)
                for projector in projectors.values()
            ),
            "REFINEMENT_PROJECTOR",
            "each actual refinement must preserve all four exact projectors",
        )
        tensor_permutation = tensor_identity_permutation(
            permutation, matter_dimension
        )
        require(
            sorted(tensor_permutation) == list(range(tensor_dimension)),
            "TENSOR_TRANSPORT_BIJECTION",
            "U tensor I_15 must be a permutation of the tensor basis",
        )

        p3 = projectors["3"]
        tensor_checks = 0
        for i in range(PORTS):
            for j in range(PORTS):
                for matter_left in range(matter_dimension):
                    for matter_right in range(matter_dimension):
                        lhs = (
                            p3[permutation[i]][permutation[j]]
                            if matter_left == matter_right
                            else ZERO
                        )
                        rhs = (
                            p3[i][j]
                            if matter_left == matter_right
                            else ZERO
                        )
                        require(
                            lhs == rhs,
                            "RANK45_PROJECTOR_TRANSPORT",
                            "U tensor I_15 must intertwine P3 tensor I_15",
                        )
                        tensor_checks += 1

        maps[(source, target)] = permutation
        rows.append(
            {
                "source": source,
                "target": target,
                "port_permutation": list(permutation),
                "all_four_projectors_natural": True,
                "tensor_formula": (
                    "U_refinement tensor I_15 in refinement-transported "
                    "coordinates on M15"
                ),
                "tensor_dimension": tensor_dimension,
                "tensor_permutation": list(tensor_permutation),
                "tensor_permutation_sha256": sha256_json(
                    list(tensor_permutation)
                ),
                "rank45_projector_intertwined": True,
                "rank45_entry_checks": tensor_checks,
            }
        )

    expected_keys = {("r0", "r1"), ("r1", "r2"), ("r0", "r2")}
    require(
        set(maps) == expected_keys,
        "REFINEMENT_MAP_SET",
        "the exact r0/r1/r2 refinement triangle must be present",
    )
    require(
        compose_permutations(maps[("r1", "r2")], maps[("r0", "r1")])
        == maps[("r0", "r2")],
        "REFINEMENT_COCYCLE",
        "U_12 after U_01 must equal the direct U_02 transport",
    )
    tensor_maps = {
        key: tensor_identity_permutation(value, matter_dimension)
        for key, value in maps.items()
    }
    require(
        compose_permutations(
            tensor_maps[("r1", "r2")], tensor_maps[("r0", "r1")]
        )
        == tensor_maps[("r0", "r2")],
        "TENSOR_REFINEMENT_COCYCLE",
        "the U tensor I_15 transports must satisfy the same cocycle",
    )

    p3_trace = mat_trace(projectors["3"])
    require(
        p3_trace == F5(3, 0),
        "FAMILY_PROJECTOR_RANK",
        "the selected family projector must have exact rank three",
    )
    tensor_rank = int(p3_trace.a) * matter_dimension
    require(
        tensor_rank == 45,
        "RANK_45_TRANSPORT",
        "P3 tensor I_15 must have exact rank forty-five",
    )
    fock_dimension = int(
        matter_receipt.get("auxiliary_car_fock", {}).get("dimension", 0)
    )
    require(
        fock_dimension == 32,
        "FOCK_DIMENSION",
        "the #314 ambient CAR Fock carrier must have dimension thirty-two",
    )
    return (
        {
            "carrier_levels": ["r0", "r1", "r2"],
            "actual_map_count": len(rows),
            "maps": rows,
            "projectors_checked": ["1", "3", "3p", "5"],
            "all_projectors_refinement_natural": True,
            "port_cocycle": "U_12 after U_01 = U_02",
            "port_cocycle_checked": True,
            "invariant_transport": (
                "mapwise U_3 tensor (gamma_refinement restricted to M15)"
            ),
            "ambient_projector": "P3 tensor Q15",
            "coordinate_transport": (
                "U_refinement tensor I_15 on the abstract fifteen-state "
                "label-copy carrier"
            ),
            "coordinate_trivialization_scope": (
                "finite label bookkeeping only; it does not replace the "
                "nontrivial #314 gamma intertwiner"
            ),
            "restricted_label_carrier_dimension": tensor_dimension,
            "full_fock_ambient_dimension": PORTS * fock_dimension,
            "restricted_tensor_projector": "P3 tensor I_M15",
            "tensor_projector_rank": tensor_rank,
            "tensor_projector_refinement_natural": True,
            "tensor_cocycle_checked": True,
            "matter_intertwiner_source": (
                "hash-pinned #314 gamma refinement on the unordered "
                "conjugate rank-fifteen projector pair"
            ),
            "matter_maps_intertwined": (
                matter_receipt["refinement"]["maps"]
            ),
            "physical_persistence_maps_intertwined": (
                matter_receipt["refinement"]["physical_maps"]
            ),
            "actual_gamma_strict_triangle_claimed": False,
            "binary_lift_centre_cocycle_retained": True,
            "one_conjugate_projector_selected": False,
            "physical_identification_promoted": False,
        },
        maps,
    )


def physical_persistence_transport(
    artifact: Mapping[str, Any],
    matter_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    rows = artifact["physical_refinement_maps"]["port_persistence_maps"]
    require(
        len(rows) == 2,
        "PHYSICAL_REFINEMENT_MAPS",
        "the #599 artifact must carry both physical persistence maps",
    )
    emitted = []
    for row in rows:
        body = {key: value for key, value in row.items() if key != "map_hash"}
        require(
            row["port_map"] == list(range(PORTS))
            and row["map_hash"] == "sha256:" + sha256_json(body),
            "PHYSICAL_REFINEMENT_HASH",
            "a physical persistence map or its hash has drifted",
        )
        emitted.append(dict(row))
    require(
        matter_receipt["refinement"]["physical_maps"]
        == [
            {
                "source_level": row["source_level"],
                "target_level": row["target_level"],
                "origin": row["origin"],
                "intertwined": True,
            }
            for row in rows
        ],
        "PHYSICAL_MATTER_REFINEMENT",
        "the #314 physical matter maps do not bind the #599 persistence maps",
    )
    return {
        "namespace": "artifact_physical_persistence",
        "maps": emitted,
        "distinct_from_algebraic_r_tower": True,
        "matter_gamma_intertwined_mapwise": True,
        "laboratory_current_identified": False,
    }


def verify_tripled_anomaly_and_z6(
    generation: Mapping[str, Any],
    matter_receipt: Mapping[str, Any],
    global_receipt: Mapping[str, Any],
    seam_manifest: Mapping[str, Any],
    families: int,
) -> dict[str, Any]:
    """Recompute three-copy anomaly forms and the complete Z6 kernel."""

    recorded_anomalies = matter_receipt["anomalies"]
    require(
        recorded_anomalies["traces"]
        == {
            "SU3_cubed": "0",
            "SU3_squared_U1": "0",
            "SU2_squared_U1": "0",
            "U1_cubed": "0",
            "gravity_squared_U1": "0",
        }
        and recorded_anomalies["witten_parity"]["weak_doublets"] == 4
        and recorded_anomalies["witten_parity"]["even"] is True,
        "MATTER_ANOMALY_PACKET",
        "the #314 five-form anomaly and Witten packet has drifted",
    )
    three_copy_forms = {
        key: str(families * Fraction(value))
        for key, value in recorded_anomalies["traces"].items()
    }
    require(
        all(Fraction(value) == 0 for value in three_copy_forms.values()),
        "THREE_COPY_ANOMALY",
        "every listed anomaly form must remain zero under triplication",
    )

    weights = global_receipt["realized_weight_table"]
    state_rows = generation["states"]
    labels = {str(row["label"]) for row in state_rows}
    require(
        labels == {"Q", "u_c", "e_c", "d_c", "L"},
        "Z6_GENERATION_LABELS",
        "the generation and #567 weight table must name the same five fields",
    )

    kernel = []
    for color in range(3):
        for weak in range(2):
            for hypercharge in range(6):
                fixed = True
                for label in labels:
                    weight = weights[label]
                    phase = (
                        2 * color * int(weight["triality"])
                        + 3 * weak * int(weight["duality"])
                        + hypercharge * int(weight["q"])
                    ) % 6
                    fixed = fixed and phase == 0
                if fixed:
                    kernel.append((color, weak, hypercharge))

    recorded_kernel = sorted(
        tuple(int(value) for value in row)
        for row in global_receipt["kernel_on_realized_tensors"][
            "kernel_elements"
        ]
    )
    require(
        sorted(kernel) == recorded_kernel
        and len(kernel) == 6,
        "Z6_KERNEL_REPLAY",
        "the exhaustive thirty-six-candidate replay must recover the #567 Z6",
    )
    generator = tuple(
        int(value)
        for value in global_receipt["kernel_on_realized_tensors"][
            "cyclic_generator"
        ]
    )
    require(
        generator == (1, 1, 1),
        "Z6_GENERATOR",
        "the pinned diagonal kernel generator must remain (1,1,1)",
    )

    tripled_rows = []
    fixed_states = 0
    for row in state_rows:
        label = str(row["label"])
        weight = weights[label]
        phase = (
            2 * generator[0] * int(weight["triality"])
            + 3 * generator[1] * int(weight["duality"])
            + generator[2] * int(weight["q"])
        ) % 6
        require(
            phase == 0,
            "Z6_TRIPLED_ACTION",
            "the diagonal Z6 generator must fix each tripled matter field",
        )
        multiplicity = families * int(row["weyl_states"])
        fixed_states += multiplicity
        tripled_rows.append(
            {
                "field": label,
                "phase_sixths": phase,
                "one_family_states": int(row["weyl_states"]),
                "three_copy_states": multiplicity,
            }
        )
    require(
        fixed_states == 45,
        "Z6_TRIPLED_RANK",
        "the diagonal kernel must fix all forty-five tensor states",
    )

    diagonal = seam_manifest["diagonal_kernel_action"]
    require(
        tuple(diagonal["generator_color_weak_hypercharge"]) == generator
        and diagonal["module_dimension"] == 15
        and diagonal["fixed_subspace_dimension"] == 15,
        "Z6_SEAM_CONSISTENCY",
        "the #627 diagonal action must agree with the #567 kernel replay",
    )
    weak_doublets = families * int(generation["weak_doublets_per_family"])
    require(
        weak_doublets == 12 and weak_doublets % 2 == 0,
        "THREE_COPY_WEAK_PARITY",
        "three generations must carry twelve weak doublets",
    )
    return {
        "family_copies": families,
        "three_copy_anomaly_forms": three_copy_forms,
        "all_listed_anomalies_zero": True,
        "three_copy_weak_doublets": weak_doublets,
        "weak_parity_even": True,
        "center_candidates_replayed": 36,
        "kernel_elements": [list(row) for row in sorted(kernel)],
        "kernel_order": len(kernel),
        "kernel_generator": list(generator),
        "kernel_after_triplication": "same diagonal Z6",
        "triplication_does_not_change_kernel": True,
        "tripled_field_action": tripled_rows,
        "fixed_subspace_dimension": fixed_states,
    }


def uniform_seam_menu_transport(
    seam_manifest: Mapping[str, Any],
    refinement: Mapping[str, Any],
    families: int,
) -> dict[str, Any]:
    """Transport every #627 conditional character, without choosing one."""

    menu = seam_manifest["hypercharge_character_menu"]["groups"]
    require(
        [int(row["group_order"]) for row in menu] == [2, 3, 6],
        "SEAM_MENU_GROUPS",
        "the conditional seam menu must contain exactly Z2, Z3, and Z6",
    )
    map_count = int(refinement["actual_map_count"])
    require(
        map_count == 3
        and refinement["tensor_projector_refinement_natural"] is True,
        "SEAM_REFINEMENT_INPUT",
        "the seam menu requires the exact three-map rank-45 transport",
    )

    rows = []
    character_count = 0
    for group in menu:
        characters = []
        for character in group["characters"]:
            fixed = int(character["fixed_subspace_dimension"])
            require(
                0 <= fixed <= 15,
                "SEAM_FIXED_DIMENSION",
                "a conditional character has an invalid fixed-space dimension",
            )
            characters.append(
                {
                    "character_exponent": int(
                        character["character_exponent"]
                    ),
                    "faithful_on_module": bool(
                        character["faithful_on_module"]
                    ),
                    "matter_fixed_dimension": fixed,
                    "rank3_family_tensor_fixed_dimension": families * fixed,
                    "refinement_maps_checked": map_count,
                    "label_transport_commutes_with_U_tensor_I15": True,
                    "full_gamma_transport_status": (
                        "conditional_on_supplied_character_compatibility"
                    ),
                }
            )
            character_count += 1
        rows.append(
            {
                "group_order": int(group["group_order"]),
                "characters": characters,
                "uniform_transport": True,
            }
        )
    require(
        character_count == 11,
        "SEAM_MENU_SIZE",
        "the complete conditional Z2/Z3/Z6 character menu has eleven rows",
    )
    return {
        "source_issue": 627,
        "classification": "uniform_conditional_tensor_transport",
        "factorization": (
            "carrier refinement acts as U on the family factor; each "
            "conditional seam character acts on the fifteen-state matter "
            "factor, so U tensor I_15 commutes with I_family tensor chi"
        ),
        "groups": rows,
        "characters_checked": character_count,
        "refinement_maps_checked_per_character": map_count,
        "physical_seam_action_selected": False,
        "selected_character": None,
        "selected_two_representation": None,
        "selection_interface": "physical_sector_mechanism_selection",
        "selection_status": "open",
    }


def control_refinement_projector_mutation(
    maps: Mapping[tuple[str, str], tuple[int, ...]],
    projectors: Mapping[str, Matrix],
) -> dict[str, Any]:
    doctored = list(maps[("r0", "r1")])
    doctored[0], doctored[2] = doctored[2], doctored[0]
    try:
        require(
            all(
                projector_fixed_by(projector, doctored)
                for projector in projectors.values()
            ),
            "REFINEMENT_PROJECTOR",
            "a doctored port map must not preserve the exact projector packet",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "REFINEMENT_PROJECTOR",
        }
    return {"expected_failure": True, "failed": False}


def control_refinement_cocycle_mutation(
    maps: Mapping[tuple[str, str], tuple[int, ...]],
) -> dict[str, Any]:
    doctored_direct = tuple(range(PORTS))
    try:
        require(
            compose_permutations(
                maps[("r1", "r2")], maps[("r0", "r1")]
            )
            == doctored_direct,
            "REFINEMENT_COCYCLE",
            "the direct map may not be replaced by the identity",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "REFINEMENT_COCYCLE",
        }
    return {"expected_failure": True, "failed": False}


def control_matter_factor_mutation(
    maps: Mapping[tuple[str, str], tuple[int, ...]],
    matter_dimension: int,
) -> dict[str, Any]:
    port_map = maps[("r0", "r1")]
    expected = tensor_identity_permutation(port_map, matter_dimension)
    doctored = tuple(
        int(port_map[port]) * matter_dimension
        + ((matter + 1) % matter_dimension)
        for port in range(PORTS)
        for matter in range(matter_dimension)
    )
    try:
        require(
            doctored == expected,
            "LABEL_FACTOR_IDENTITY",
            "the declared label-bookkeeping map must retain I_15",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "LABEL_FACTOR_IDENTITY",
        }
    return {"expected_failure": True, "failed": False}


def control_z6_generator_mutation(
    generation: Mapping[str, Any],
    global_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    weights = global_receipt["realized_weight_table"]
    mutated_generator = (1, 1, 0)
    try:
        for row in generation["states"]:
            weight = weights[row["label"]]
            phase = (
                2 * mutated_generator[0] * int(weight["triality"])
                + 3 * mutated_generator[1] * int(weight["duality"])
                + mutated_generator[2] * int(weight["q"])
            ) % 6
            require(
                phase == 0,
                "Z6_GENERATOR_MUTATION",
                "a mutated center generator must fail on realized matter",
            )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "Z6_GENERATOR_MUTATION",
        }
    return {"expected_failure": True, "failed": False}


def control_seam_selection_promotion(
    seam_transport: Mapping[str, Any],
) -> dict[str, Any]:
    try:
        require(
            seam_transport["physical_seam_action_selected"] is True,
            "PHYSICAL_SEAM_SELECTION_OPEN",
            "uniform conditional transport does not select a physical seam action",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "PHYSICAL_SEAM_SELECTION_OPEN",
        }
    return {"expected_failure": True, "failed": False}


def control_strict_matter_cocycle_promotion(
    refinement: Mapping[str, Any],
) -> dict[str, Any]:
    try:
        require(
            refinement["actual_gamma_strict_triangle_claimed"] is True,
            "STRICT_MATTER_COCYCLE_NOT_DERIVED",
            "mapwise #314 intertwiners do not supply a strict spin-lift section",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "STRICT_MATTER_COCYCLE_NOT_DERIVED",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Candidate enumeration and the selection theorem
# ---------------------------------------------------------------------------


BAND_ORDER = ["1", "3", "3p", "5"]
BAND_DIMS = {"1": 1, "3": 3, "3p": 3, "5": 5}


def band_costs(coefficients: tuple[int, int]) -> dict[str, F5]:
    """Quadratic readback value per unit norm on each adjacency eigenband.

    A form aI + bA restricted to the eigenband of eigenvalue lambda acts
    as the scalar a + b*lambda; on the Laplacian convention (5, -1) this
    is the exact per-unit-norm seam-mismatch cost of the band.
    """

    a, b = coefficients
    eigen = {"1": F5(5, 0), "3": SQRT5, "3p": -SQRT5, "5": F5(-1, 0)}
    return {name: F5(a, 0) + F5(b, 0) * value for name, value in eigen.items()}


def band_action_kernels(
    projectors: Mapping[str, Matrix], rotations: Sequence[tuple[int, ...]]
) -> dict[str, int]:
    """For each band, the number of listed rotations acting as the identity.

    The columns of the symmetric projector P span the band, and a
    permutation g fixes every band vector exactly when P[g(i)][j] equals
    P[i][j] for all ports i, j.  Counting over the listed group is
    inverse-symmetric, so the count is the kernel order of the band
    action.
    """

    kernels: dict[str, int] = {}
    for name, p in projectors.items():
        count = 0
        for g in rotations:
            gp = [[p[g[i]][j] for j in range(PORTS)] for i in range(PORTS)]
            if mat_eq(gp, p):
                count += 1
        kernels[name] = count
    return kernels


def enumerate_candidates(
    window: tuple[int, int],
    kernels: Mapping[str, int],
    costs: Mapping[str, F5],
) -> dict[str, Any]:
    lower, upper = window
    rows: list[dict[str, Any]] = []
    for mask in range(1, 16):
        parts = [BAND_ORDER[k] for k in range(4) if mask & (1 << k)]
        dim = sum(BAND_DIMS[p] for p in parts)
        label = "+".join(parts)
        row: dict[str, Any] = {"object": label, "dimension": dim}
        if len(parts) > 1:
            row["excluded_by"] = "single_complete_object_clause"
            row["reason"] = "a source-visible proper splitting projector exists"
        elif kernels[parts[0]] == 60:
            row["excluded_by"] = "faithful_family_exchange"
            row["reason"] = "every listed rotation acts as the identity"
        elif not (lower <= dim <= upper):
            row["excluded_by"] = "physical_window"
            row["reason"] = f"dimension outside [{lower}, {upper}]"
        else:
            row["admissible"] = True
            row["cost_per_unit_norm"] = f5_str(costs[parts[0]])
        rows.append(row)
    admissible = [r["object"] for r in rows if r.get("admissible")]
    require(
        sorted(admissible) == ["3", "3p", "5"],
        "CANDIDATE_SET",
        "the admissible candidates must be exactly the 3, 3', and 5 bands",
    )
    return {"rows": rows, "admissible": admissible}


def strict_minimizer(
    admissible: Sequence[str], costs: Mapping[str, F5]
) -> dict[str, Any]:
    ordered = f5_sorted(admissible, lambda name: costs[name])
    for left, right in zip(ordered, ordered[1:]):
        require(
            f5_lt(costs[left], costs[right]),
            "COST_ORDER_NOT_STRICT",
            "the band costs must be strictly totally ordered",
        )
    winner = ordered[0]
    return {
        "order": [{"object": name, "cost": f5_str(costs[name])} for name in ordered],
        "strict": True,
        "minimizer": winner,
    }


# ---------------------------------------------------------------------------
# Controls (every control must fail closed)
# ---------------------------------------------------------------------------


def control_external_copy_reduct(generation: Mapping[str, Any]) -> dict[str, Any]:
    """External C^n completions stay reduct-indistinguishable (#617 intact).

    The per-copy reduct data is derived from the actual generation
    certificate, not from literals: n copies scale every per-family
    anomaly form (all exactly zero) and the doublet parity (4n mod 2).
    Attempting to select a copy count from that reduct must refuse.
    """

    per_copy = {}
    for n in (3, 4):
        per_copy[n] = {
            "anomaly_forms_scaled": {
                key: str(Fraction(value) * n)
                for key, value in generation["per_family_anomaly_forms"].items()
            },
            "weak_parity": (int(generation["weak_doublets_per_family"]) * n) % 2,
            "kernel": "Z6",
        }
    indistinguishable = per_copy[3] == per_copy[4]
    try:
        require(
            not indistinguishable,
            "REDUCT_COUNT_NOT_SELECTED",
            "the family-free reduct cannot distinguish in-window copy counts",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "REDUCT_COUNT_NOT_SELECTED",
            "reduct_rows": {str(n): row for n, row in per_copy.items()},
            "meaning": "the #617 invisibility theorem is preserved; without clause (R) nothing is selected",
        }
    return {"expected_failure": True, "failed": False}


def control_reducible_object(window: tuple[int, int]) -> dict[str, Any]:
    """The reducible 1+3 object (dimension four, inside the window) must be
    rejected by the single-complete-object clause, not by the window."""

    lower, upper = window
    dim = BAND_DIMS["1"] + BAND_DIMS["3"]
    inside = lower <= dim <= upper
    try:
        require(
            not inside,
            "NOT_SINGLE_COMPLETE_OBJECT",
            "1+3 sits inside the window and must be excluded by the object clause",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "NOT_SINGLE_COMPLETE_OBJECT",
            "meaning": "the window alone does not exclude reducible objects; the single-complete-object clause is load-bearing",
        }
    return {"expected_failure": True, "failed": False}


def control_excluded_cone(
    admissible: Sequence[str], excluded_coefficients: tuple[int, int]
) -> dict[str, Any]:
    """Selecting with the excluded comparison readback 6I + A must flip the
    minimizer, so clause (S) is load-bearing: the certificate refuses the
    imported-only form."""

    wrong = band_costs(excluded_coefficients)
    ordered = f5_sorted(admissible, lambda name: wrong[name])
    flipped = ordered[0]
    try:
        require(
            flipped == "3",
            "COST_CONE_VIOLATION",
            "the excluded readback selects the wrong band",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "COST_CONE_VIOLATION",
            "excluded_readback_minimizer": flipped,
            "excluded_readback_costs": {name: f5_str(wrong[name]) for name in admissible},
            "meaning": "the selection genuinely consumes the #625 operational cone; the imported-only form 6I + A picks the Galois partner instead",
        }
    return {"expected_failure": True, "failed": False}


def control_galois_transport(
    projectors: Mapping[str, Matrix],
    rotations: Sequence[tuple[int, ...]],
    costs: Mapping[str, F5],
) -> dict[str, Any]:
    """The Galois automorphism swaps the bands and reverses the cost order,
    and no listed transport realizes it, so it must be refused as a
    source transport."""

    p3, p3p = projectors["3"], projectors["3p"]
    realized = False
    for g in rotations:
        gp = [[p3[g[i]][g[j]] for j in range(PORTS)] for i in range(PORTS)]
        if mat_eq(gp, p3p):
            realized = True
            break
    order_reversed = f5_lt(costs["3"], costs["3p"]) and f5_lt(
        costs["3p"].conj(), costs["3"].conj()
    )
    try:
        require(
            realized or not order_reversed,
            "GALOIS_NOT_A_TRANSPORT",
            "the Galois swap is not induced by any listed rotation and reverses the measured cost order",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "GALOIS_NOT_A_TRANSPORT",
            "meaning": "sqrt5 -> -sqrt5 exchanges the bands abstractly, but it is not an A2 transport and does not preserve the measured cost order, so the 3 versus 3' ambiguity is operationally resolved",
        }
    return {"expected_failure": True, "failed": False}


def control_block_swap(costs: Mapping[str, F5]) -> dict[str, Any]:
    """The unitary block swap between the two isometric band embeddings
    (the 2026-07-20 reopening witness) changes the exact cost value, so it
    is excluded from the family-relabelling groupoid."""

    swapped_cost = costs["3p"]
    try:
        require(
            swapped_cost == costs["3"],
            "BLOCK_SWAP_COST_DETECTED",
            "the block swap is not cost-preserving",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "BLOCK_SWAP_COST_DETECTED",
            "cost_before": f5_str(costs["3"]),
            "cost_after": f5_str(swapped_cost),
            "meaning": "the refinement witness that defeated the 2026-07-20 closure attempt is detected by the cost readback; in-band relabellings preserve the cost exactly",
        }
    return {"expected_failure": True, "failed": False}


def control_dropped_faithfulness(costs: Mapping[str, F5]) -> dict[str, Any]:
    """Without the faithfulness clause the minimizer over all four bands is
    the trivial band at cost zero, so cost minimization alone must never
    be claimed to force three families."""

    ordered = f5_sorted(BAND_ORDER, lambda name: costs[name])
    unguarded = ordered[0]
    try:
        require(
            unguarded == "3",
            "TRIVIAL_BAND_WITHOUT_FAITHFULNESS",
            "cost minimization alone selects the trivial band",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "TRIVIAL_BAND_WITHOUT_FAITHFULNESS",
            "unguarded_minimizer": unguarded,
            "unguarded_cost": f5_str(costs[unguarded]),
            "meaning": "the faithfulness clause is load-bearing; the selection is cost minimization among faithful in-window single objects, not cost minimization alone",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Payload assembly
# ---------------------------------------------------------------------------


def require_no_floats(value: Any, path: str = "$") -> None:
    require(
        not isinstance(value, float),
        "FLOAT_FORBIDDEN",
        f"a float appears in the payload at {path}",
    )
    if isinstance(value, Mapping):
        for key, item in value.items():
            require_no_floats(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            require_no_floats(item, f"{path}[{index}]")


def build_payload() -> dict[str, Any]:
    carrier, rotations, carrier_pin = load_carrier()
    adjacency = lift(adjacency_int(carrier))
    structural, structural_pins = pin_structural_chain(carrier_pin["sha256"])

    projectors = spectral_projectors(adjacency)
    spectral = verify_spectral_resolution(adjacency, projectors)
    equivariance = verify_equivariance(projectors, rotations)
    pair_orbits = verify_pair_orbits(carrier, rotations)
    multiplicity_free = verify_multiplicity_free(carrier)

    lower, upper, window_pin = pin_window()
    operational, excluded, cone_pin = pin_cost_cone()
    charges, matter_pin = pin_matter()
    require(
        matter_pin["sha256"] == structural_pins["matter_manifest"]["sha256"],
        "MATTER_PIN_CHAIN",
        "the generation charge packet and structural matter packet differ",
    )
    artifact, response_pin = pin_response_artifact(carrier_pin)
    pole_artifact, pole_pin = pin_pole_residue_artifact(carrier_pin, artifact)
    attachment_receipt, attachment_pin = pin_matter_attachment_receipt(
        artifact, pole_artifact
    )

    costs = band_costs(operational)
    kernels = band_action_kernels(projectors, rotations)
    require(
        kernels == {"1": 60, "3": 1, "3p": 1, "5": 1},
        "BAND_KERNELS",
        "the trivial band must absorb all sixty rotations and every other band action must be faithful",
    )
    measured = measured_band_receipt(
        carrier, adjacency, projectors, artifact, operational, kernels
    )

    candidates = enumerate_candidates((lower, upper), kernels, costs)
    minimizer = strict_minimizer(candidates["admissible"], costs)
    require(minimizer["minimizer"] == "3", "SELECTED_BAND", "the strict minimizer must be the 3 band")

    generation = generation_certificate(charges)
    pole_receipt = pole_residue_receipt(pole_artifact, generation)
    families = BAND_DIMS["3"]
    refinement_transport, refinement_maps = verify_refinement_transport(
        carrier,
        rotations,
        projectors,
        int(generation["weyl_state_count"]),
        structural["matter_receipt"],
    )
    physical_transport = physical_persistence_transport(
        artifact, structural["matter_receipt"]
    )
    three_copy_consistency = verify_tripled_anomaly_and_z6(
        generation,
        structural["matter_receipt"],
        structural["global_form_receipt"],
        structural["seam_manifest"],
        families,
    )
    seam_transport = uniform_seam_menu_transport(
        structural["seam_manifest"],
        refinement_transport,
        families,
    )
    attachment = {
        "family_object": "the 3 band of the screen coefficient space",
        "family_dimension": families,
        "generation_weyl_states": generation["weyl_state_count"],
        "complex_rank": families * int(generation["weyl_state_count"]),
        "invariant_projector": "P3 tensor Q15",
        "coordinate_projector": (
            "P3 tensor I_M15 in a #314 refinement-transported basis"
        ),
        "one_conjugate_matter_projector_selected": False,
        "three_family_anomaly_forms": three_copy_consistency[
            "three_copy_anomaly_forms"
        ],
        "three_family_weak_doublets": three_copy_consistency[
            "three_copy_weak_doublets"
        ],
        "weak_parity_even": three_copy_consistency["weak_parity_even"],
        "common_kernel": "Z6, unchanged under family triplication",
        "common_kernel_recomputed": True,
    }
    require(
        attachment["complex_rank"] == 45,
        "RANK_45",
        "the realized attachment must have complex rank forty-five",
    )

    uniqueness = {
        "embedding": "multiplicity one: the equivariant embedding of the 3 band is unique up to a scalar",
        "relabellings": "the induced family-exchange image is the faithful icosahedral rotation image on the band; in-band relabellings preserve the cost exactly",
        "galois_branch": "the 3' band is the Galois partner; it is operationally separated by the strict cost order and by the absence of any listed transport realizing the swap",
    }

    controls = {
        "external_copy_reduct": control_external_copy_reduct(generation),
        "reducible_object": control_reducible_object((lower, upper)),
        "excluded_cone": control_excluded_cone(candidates["admissible"], excluded),
        "galois_transport": control_galois_transport(projectors, rotations, costs),
        "block_swap_refinement": control_block_swap(costs),
        "dropped_faithfulness": control_dropped_faithfulness(costs),
        "measured_channel_swap": control_measured_channel_swap(
            artifact, operational, kernels
        ),
        "pole_table_mutation": control_pole_table_mutation(pole_artifact),
        "matter_attachment_mutation": control_matter_attachment_mutation(
            artifact, pole_artifact
        ),
        "refinement_projector_mutation": (
            control_refinement_projector_mutation(
                refinement_maps, projectors
            )
        ),
        "refinement_cocycle_mutation": control_refinement_cocycle_mutation(
            refinement_maps
        ),
        "matter_factor_mutation": control_matter_factor_mutation(
            refinement_maps, int(generation["weyl_state_count"])
        ),
        "z6_generator_mutation": control_z6_generator_mutation(
            generation, structural["global_form_receipt"]
        ),
        "seam_selection_promotion": control_seam_selection_promotion(
            seam_transport
        ),
        "strict_matter_cocycle_promotion": (
            control_strict_matter_cocycle_promotion(refinement_transport)
        ),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "claim_boundary": (
            "Exact selection inside the source-visible screen: among single "
            "complete faithful in-window multiplicity objects the #625 "
            "operational comparison order has the 3 band as unique strict "
            "minimizer, fixing N_g = 3 with attachment rank forty-five. "
            "The #617 copy-count invisibility for external completions is "
            "preserved. Clause S is simulator-realized by the #599 response "
            "artifact; clause R is simulator-realized for the response "
            "resolvent of the declared Laplacian generator by the pole-residue "
            "artifact. Its rank-three frame residue sits at the lowest positive "
            "generator frequency with the generation factor imported. All "
            "three actual refinements preserve the finite projector and its "
            "rank-forty-five tensor transport. The #627 seam character menu "
            "transports conditionally without selecting an action. Matter-pole "
            "identification, continuum Spin/locality, physical seam selection, "
            "and laboratory-current attachment stay open."
        ),
        "named_interface": {
            "id": "screen_realized_multiplicity_object",
            "class": "conditional_open_interface",
            "clauses": {
                "R_realization": (
                    "the physical pole-residue multiplicity object is realized "
                    "as a single complete subobject of the source-visible "
                    "screen coefficient space"
                ),
                "S_selection": (
                    "the attachment is compared by the #625 operational cost "
                    "order (quadratic seam-mismatch readback per unit norm)"
                ),
            },
            "clause_controls": {
                "R_realization": "external_copy_reduct",
                "S_selection": "excluded_cone",
            },
            "clause_status": {
                "R_realization": (
                    "simulator_realized_for_declared_response_resolvent__"
                    "matter_pole_identification_open"
                ),
                "S_selection": measured["clause_S"],
            },
            "open_receipts": [
                "matter-pole identification",
                "continuum Spin/locality receipt",
                "physical seam action selection",
                "laboratory current identification",
            ],
        },
        "upstream_pins": {
            "carrier": carrier_pin,
            **structural_pins,
            "multiplicity_window": window_pin,
            "operational_cost_cone": cone_pin,
            "matter_packet": matter_pin,
            "measured_response_artifact": response_pin,
            "measured_pole_residue_artifact": pole_pin,
            "finite_spin_locality_receipt": attachment_pin,
        },
        "measured_receipt": measured,
        "pole_residue_receipt": pole_receipt,
        "finite_spin_locality_attachment": {
            "receipt": attachment_pin,
            "domain": (
                "issue-634 typed local domain, bound by capture hash to "
                "the frozen stage receipts"
            ),
            "layer_clauses_all_true": True,
            "generation_recomputed_states": attachment_receipt[
                "generation_certificate"
            ]["weyl_state_count"],
            "z6_all_states_fixed": attachment_receipt[
                "z6_kernel_certificate"
            ]["all_states_fixed"],
            "chirality_nondegenerate": attachment_receipt[
                "chirality_certificate"
            ]["chirality_nondegenerate"],
            "gap_inherited_exact": attachment_receipt[
                "gap_inheritance_certificate"
            ]["inherited"],
            "lift_ambiguity_rank": attachment_receipt["spin_layer"][
                "lift_ambiguity_rank"
            ],
            "continuum_gate_unchanged": True,
            "scope_note": (
                "the finite layer types chirality, spin data, locality, "
                "and refinement on the local domain; the continuum "
                "Spin/locality receipt, matter-pole identification, "
                "physical seam selection, and laboratory identification "
                "stay open"
            ),
        },
        "spectral_resolution": spectral,
        "equivariance": equivariance,
        "refinement_transport": refinement_transport,
        "physical_persistence_transport": physical_transport,
        "three_copy_consistency": three_copy_consistency,
        "conditional_seam_menu_transport": seam_transport,
        "pair_orbits": pair_orbits,
        "multiplicity_free": multiplicity_free,
        "band_action_kernels": kernels,
        "physical_window": {"lower": lower, "upper": upper, "source": "pinned #617 receipt"},
        "operational_cost": {
            "form": "5I - A per unit norm on each band",
            "coefficients": list(operational),
            "excluded_comparison_form": list(excluded),
            "band_costs": {name: f5_str(value) for name, value in costs.items()},
        },
        "candidate_enumeration": candidates,
        "selection": minimizer,
        "generation": generation,
        "attachment": attachment,
        "uniqueness": uniqueness,
        "controls": controls,
        "promotion": {
            "matter_pole_identified": False,
            "continuum_spin_locality_derived": False,
            "physical_seam_action_selected": False,
            "laboratory_current_identified": False,
            "promotion_allowed": False,
        },
        "open_gates": [
            "matter_pole_identification",
            "continuum_Spin_locality",
            "physical_seam_action_selection",
            "laboratory_current_identification",
        ],
        "invisibility_preserved": True,
        "bounded_exit": "exact_finite_refinement_transport",
        "lean_spine": [
            "Screen/A5PortAction.lean (sixty listed rotations, kernel-decided)",
            "Screen/A5Commutant.lean (four-dimensional commutant)",
            (
                "Screen/A5FamilyBand.lean (spectral split, actual refinement "
                "cocycle, rank-45 tensor transport, anomaly and Z6 checks)"
            ),
        ],
    }
    require_no_floats(payload)
    return payload


def build_manifest() -> dict[str, Any]:
    payload = build_payload()
    manifest = dict(payload)
    manifest["manifest_sha256"] = "sha256:" + sha256_json(payload)
    return manifest


def verify_stored() -> dict[str, Any]:
    stored = load_json(MANIFEST_PATH)
    body = {key: value for key, value in stored.items() if key != "manifest_sha256"}
    require(
        stored.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "stored manifest hash does not match its body",
    )
    rebuilt = build_payload()
    require(
        body == rebuilt,
        "MANIFEST_DRIFT",
        "stored manifest does not match a deterministic rebuild",
    )
    return {"status": "PASS", "manifest": str(MANIFEST_PATH)}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Family band attachment certificate for issue #569")
    parser.add_argument("--verify", action="store_true", help="compare the stored manifest with a rebuild")
    parser.add_argument("--output", type=Path, default=MANIFEST_PATH)
    args = parser.parse_args(argv)
    if args.verify:
        print(json.dumps(verify_stored(), indent=2))
        return 0
    manifest = build_manifest()
    write_json(args.output, manifest)
    print(json.dumps({"status": "WROTE", "manifest": str(args.output), "manifest_sha256": manifest["manifest_sha256"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
