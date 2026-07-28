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
* The realized attachment object, band tensor generation, has complex
  rank exactly 3 x 15 = 45; the fifteen-state generation is recomputed
  from the pinned block charges, every per-family anomaly form is zero,
  three-family weak parity is even, and the common kernel is unchanged.

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

Neither clause is derived from A1-A3.  The measured #599 response
artifact carries clause S as a measurement: the per-band adjacency
channels are measured data, the operational cost evaluated on them has
the frame triplet as strict minimizer, and conjugation swaps the
measured frame and kernel bands while the measured order separates
them.  The carrier component of clause R is measured the same way (the
frame band is a subobject of the measured response basis); the
pole-residue factorization with complex rank forty-five, Spin/locality,
refinement, and laboratory receipts stay open on issue #569.  The #617
invisibility theorem for external C^n completions is re-verified and
holds unchanged.
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

SCHEMA = "oph.family_band_attachment_certificate.v2"
MANIFEST_PATH = MODULE_DIR / "manifests" / "family_band_attachment_reference.json"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"
WINDOW_MANIFEST_NAME = "multiplicity_window_reference.json"
READBACK_MANIFEST_NAME = "load_fiber_readback_reference.json"
MATTER_MANIFEST_NAME = "super_tannakian_matter_reference.json"
RESPONSE_ARTIFACT_NAME = "charged_response_semantic_artifact.json"

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


def measured_band_receipt(
    carrier: Any,
    adjacency: Matrix,
    projectors: Mapping[str, Matrix],
    artifact: Mapping[str, Any],
    operational: tuple[int, int],
    kernels: Mapping[str, int],
) -> dict[str, Any]:
    """Clause receipts read from the measured #599 response artifact.

    The artifact carries, as exact data: the per-band adjacency channel
    values, the sector dimensions, the Galois pairing, the response band
    scales, and the antipode polynomial of the measured response
    operator.  This receipt binds those measurements to the selection:

    * the measured channel values reproduce the exact band spectrum, so
      the operational cost evaluated on measured channels gives the
      measured band costs, and their strict order is a measured fact
      (clause S realized by measurement);
    * the measured frame band is the strict minimizer among the faithful
      bands, and it is a subobject of the screen coefficient space
      exhibited by the measured response basis (the carrier component of
      clause R measured; the pole-residue factorization receipt with
      complex rank forty-five stays open);
    * the measured antipode polynomial (A^3 - 4A^2 - 5A + 10I)/10 equals
      the carrier antipode permutation exactly, and the real part of the
      scaled band projector equals ten times identity-minus-antipode, so
      the measured response algebra and the spectral selection algebra
      are one object;
    * conjugation swaps the measured frame and kernel bands while the
      measured cost order separates them, so the Galois resolution is
      itself measured.
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
        "clause_S": "measured_realized",
        "clause_R": "carrier_component_measured__pole_residue_receipt_open",
        "open_receipt": (
            "the pole-residue factorization with complex rank forty-five "
            "through band x generation; its gates are specified on issue #569"
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

    projectors = spectral_projectors(adjacency)
    spectral = verify_spectral_resolution(adjacency, projectors)
    equivariance = verify_equivariance(projectors, rotations)
    pair_orbits = verify_pair_orbits(carrier, rotations)
    multiplicity_free = verify_multiplicity_free(carrier)

    lower, upper, window_pin = pin_window()
    operational, excluded, cone_pin = pin_cost_cone()
    charges, matter_pin = pin_matter()
    artifact, response_pin = pin_response_artifact(carrier_pin)

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
    families = BAND_DIMS["3"]
    three_family_forms = {
        key: str(Fraction(value) * families)
        for key, value in generation["per_family_anomaly_forms"].items()
    }
    require(
        all(Fraction(value) == 0 for value in three_family_forms.values()),
        "THREE_FAMILY_ANOMALY",
        "the realized three-family tensors must cancel every listed anomaly form",
    )
    doublets_total = int(generation["weak_doublets_per_family"]) * families
    attachment = {
        "family_object": "the 3 band of the screen coefficient space",
        "family_dimension": families,
        "generation_weyl_states": generation["weyl_state_count"],
        "complex_rank": families * int(generation["weyl_state_count"]),
        "three_family_anomaly_forms": three_family_forms,
        "three_family_weak_doublets": doublets_total,
        "weak_parity_even": doublets_total % 2 == 0,
        "common_kernel": "Z6, unchanged under family triplication (per-family forms are zero and the kernel is charge-determined)",
    }
    require(attachment["complex_rank"] == 45, "RANK_45", "the realized attachment must have complex rank forty-five")

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
            "preserved. Clause S is realized by the measured #599 response "
            "artifact and the carrier component of clause R is measured; "
            "the rank-forty-five pole-residue factorization is the open "
            "receipt."
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
                "R_realization": measured["clause_R"],
                "S_selection": measured["clause_S"],
            },
            "open_receipts": [
                "complex rank-45 pole-residue factorization",
                "Spin/locality receipt",
                "refinement receipt",
                "laboratory current identification",
            ],
        },
        "upstream_pins": {
            "carrier": carrier_pin,
            "multiplicity_window": window_pin,
            "operational_cost_cone": cone_pin,
            "matter_packet": matter_pin,
            "measured_response_artifact": response_pin,
        },
        "measured_receipt": measured,
        "spectral_resolution": spectral,
        "equivariance": equivariance,
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
        "invisibility_preserved": True,
        "bounded_exit": "exact_named_realization",
        "lean_spine": [
            "Screen/A5PortAction.lean (sixty listed rotations, kernel-decided)",
            "Screen/A5Commutant.lean (four-dimensional commutant)",
            "Screen/A5FamilyBand.lean (spectral split, strict cost order, minimizer)",
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
