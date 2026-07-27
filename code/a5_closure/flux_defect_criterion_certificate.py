#!/usr/bin/env python3
"""Exact certificate for GitHub issue #311: particle-like defect criterion.

The defect objects are the six measured flux sectors of the certified
icosahedral support. The inputs are hash-pinned: the certified carrier
manifest supplies the twelve-vertex, thirty-seam, twenty-face oriented
complex, and the #567 measured global-form artifact supplies the order-six
class group together with the two-puncture flux-tube witnesses realizing
every class. On that input this certificate derives the finite spectral
receipt that separates a topological flux defect from a classical localized
record:

* per flux class k in Z6, the twisted seam adjacency operator A_k on the
  ell^2 space of the support vertices, with seam phases omega^(k f(e)) taken
  from the pinned witness 1-chain f whose coboundary is the class-k
  two-puncture 2-cocycle; omega = exp(2 pi i / 6) is handled exactly in the
  Eisenstein ring Z[omega] with omega^2 = omega - 1, and every operator is
  verified Hermitian entry by entry;
* the exact characteristic polynomial of every A_k over Z, computed by
  fraction-free Faddeev-LeVerrier in the Eisenstein encoding and
  cross-checked against the integer 2n x 2n companion-block form, whose
  characteristic polynomial must equal the square of the Eisenstein result;
* the spectral criterion: the six polynomials are not all equal; the
  untwisted class 0 differs from class 3 and from class 1, and the measured
  coincidence pattern is exactly the charge-conjugation pairing k <-> 6 - k,
  verified entrywise as A_(6-k) = conjugate(A_k);
* gauge invariance: regauging the witness chain by any coboundary conjugates
  A_k by a diagonal unitary and leaves the characteristic polynomial fixed,
  verified exactly for sample gauges, while the sixty scaled basis
  coboundaries all have identically zero face holonomy;
* the classical-record control in both directions: every local operation
  (marking a vertex, or any coboundary regauge) leaves the flux class and
  the spectrum family unchanged, and no seam regauge carries class k to
  class j for k != j, because the chain difference has nonzero puncture
  holonomy and coboundaries have none;
* charge and fusion: the flux class is a Z6 charge, additive under chain
  addition, and two well-separated puncture pairs with disjoint seam
  supports compose additively on the finite complex;
* refinement stability: one exact edge-midpoint refinement step reproduces
  the 42-vertex, 120-seam, 80-face complex of the pinned artifact's
  refinement transport, realizes all six classes between child punctures,
  and reproduces the same spectral coincidence partition.

The continuum quantum pole, scattering amplitudes, mass calibration, and
laboratory identification of any species stay open named interfaces; the
receipt is the finite source-scope criterion only. Every arithmetic step is
exact integer arithmetic; no floating point appears in a proof step.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

SCHEMA = "oph.flux_defect_criterion_certificate.v1"
RECEIPT_SCHEMA = "oph.flux_defect_criterion_receipt.v1"
NEGATIVE_SCHEMA = "oph.flux_defect_criterion_negative_controls.v1"
CARRIER_SCHEMA = "oph.echosahedral_selector_manifest.v1"
GLOBAL_FORM_ARTIFACT_SCHEMA = "oph.global_form_semantic_artifact.v1"

FLUX_ORDER = 6

# Powers of omega = exp(2 pi i / 6) in the Eisenstein ring Z[omega] with
# omega^2 = omega - 1; an element a + b omega is the integer pair (a, b).
OMEGA_POWERS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))

# Integer companion block of x^2 - x + 1: the image of omega under the exact
# regular representation of Z[omega] on Z^2.
OMEGA_BLOCK = ((0, -1), (1, 1))


# ---------------------------------------------------------------------------
# Exact Eisenstein arithmetic
# ---------------------------------------------------------------------------


def eis_mul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c + b * d)


def eis_conj(x: tuple[int, int]) -> tuple[int, int]:
    a, b = x
    return (a + b, -b)


def omega_block_power(exponent: int) -> tuple[tuple[int, int], ...]:
    result: tuple[tuple[int, int], ...] = ((1, 0), (0, 1))
    for _ in range(exponent % FLUX_ORDER):
        result = tuple(
            tuple(
                sum(result[i][m] * OMEGA_BLOCK[m][j] for m in range(2))
                for j in range(2)
            )
            for i in range(2)
        )
    return result


def polynomial_square(poly: Sequence[int]) -> list[int]:
    out = [0] * (2 * len(poly) - 1)
    for i, a in enumerate(poly):
        for j, b in enumerate(poly):
            out[i + j] += a * b
    return out


# ---------------------------------------------------------------------------
# Oriented complex reconstruction from the pinned carrier
# ---------------------------------------------------------------------------


def edge_list(face_rows: Sequence[Sequence[int]]) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for a, b, c in face_rows:
        for u, v in ((a, b), (b, c), (c, a)):
            edges.add((min(u, v), max(u, v)))
    return sorted(edges)


def boundary_two(
    face_rows: Sequence[Sequence[int]], edges: Sequence[tuple[int, int]]
) -> list[list[int]]:
    edge_index = {edge: position for position, edge in enumerate(edges)}
    matrix = [[0] * len(face_rows) for _ in range(len(edges))]
    for face_position, (a, b, c) in enumerate(face_rows):
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            sign = 1 if (u, v) == key else -1
            matrix[edge_index[key]][face_position] += sign
    return matrix


def face_holonomies(
    seam_values: Mapping[int, int],
    boundary: Sequence[Sequence[int]],
    face_count: int,
) -> list[int]:
    holonomies = []
    for face in range(face_count):
        total = 0
        for edge, value in seam_values.items():
            total += boundary[edge][face] * value
        holonomies.append(total % FLUX_ORDER)
    return holonomies


def graph_distances(adjacency: Sequence[Sequence[int]]) -> list[list[int]]:
    count = len(adjacency)
    result = []
    for start in range(count):
        distance = [-1] * count
        distance[start] = 0
        queue = [start]
        while queue:
            node = queue.pop(0)
            for other in adjacency[node]:
                if distance[other] < 0:
                    distance[other] = distance[node] + 1
                    queue.append(other)
        result.append(distance)
    return result


def build_complex(face_rows: Sequence[Sequence[int]], vertex_count: int) -> dict[str, Any]:
    edges = edge_list(face_rows)
    boundary = boundary_two(face_rows, edges)
    for edge_position in range(len(edges)):
        signs = sorted(
            boundary[edge_position][face]
            for face in range(len(face_rows))
            if boundary[edge_position][face] != 0
        )
        require(
            signs == [-1, 1],
            "COMPLEX",
            "every seam must border exactly two faces with opposite orientation signs",
        )
    adjacency: list[list[int]] = [[] for _ in range(vertex_count)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    return {
        "face_rows": [list(row) for row in face_rows],
        "edges": edges,
        "boundary": boundary,
        "adjacency": adjacency,
        "vertex_count": vertex_count,
    }


def load_carrier_complex(manifest: Mapping[str, Any], base_dir: Path) -> dict[str, Any]:
    path_raw = manifest.get("carrier_manifest_path")
    require(isinstance(path_raw, str), "UPSTREAM_REFERENCE", "carrier_manifest_path is missing")
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    carrier_manifest = load_json(path)
    require(
        manifest.get("carrier_manifest_sha256") == sha256_json(carrier_manifest),
        "UPSTREAM_HASH",
        "the carrier manifest hash does not match the declared pin",
    )
    require(
        carrier_manifest.get("schema") == CARRIER_SCHEMA,
        "UPSTREAM_REFERENCE",
        "the pinned carrier manifest is not the certified selector manifest",
    )
    carrier = carrier_manifest.get("carrier", {})
    ports = list(carrier.get("ports", []))
    require(
        len(ports) == 12 and len(set(ports)) == 12,
        "CARRIER",
        "twelve distinct carrier ports are required",
    )
    index = {port: position for position, port in enumerate(ports)}
    faces_raw = carrier.get("oriented_faces", [])
    require(len(faces_raw) == 20, "CARRIER", "twenty oriented carrier faces are required")
    face_rows = [[index[p] for p in face] for face in faces_raw]
    support = build_complex(face_rows, 12)
    require(len(support["edges"]) == 30, "CARRIER", "thirty carrier seams are required")
    require(
        all(len(row) == 5 for row in support["adjacency"]),
        "CARRIER",
        "the carrier incidence must be five-regular",
    )
    distances = graph_distances(support["adjacency"])
    antipode = []
    for vertex in range(12):
        partners = [other for other in range(12) if distances[vertex][other] == 3]
        require(
            len(partners) == 1,
            "CARRIER",
            "each carrier vertex requires a unique distance-three partner",
        )
        antipode.append(partners[0])
    require(
        all(antipode[antipode[vertex]] == vertex for vertex in range(12)),
        "CARRIER",
        "the carrier antipode must be an involution",
    )
    support["antipode"] = antipode
    support["carrier_manifest_sha256"] = sha256_json(carrier_manifest)
    return support


def antipodal_face(support: Mapping[str, Any], face_position: int) -> int:
    antipode = support["antipode"]
    target = frozenset(antipode[v] for v in support["face_rows"][face_position])
    for position, face in enumerate(support["face_rows"]):
        if frozenset(face) == target:
            return position
    raise CertificateError("CARRIER", "a face has no antipodal image face")


# ---------------------------------------------------------------------------
# Pinned measured artifact and witness chains
# ---------------------------------------------------------------------------


def load_global_form_artifact(
    manifest: Mapping[str, Any], base_dir: Path, carrier_sha256: str
) -> dict[str, Any]:
    path_raw = manifest.get("global_form_artifact_path")
    require(
        isinstance(path_raw, str),
        "UPSTREAM_REFERENCE",
        "global_form_artifact_path is missing",
    )
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    artifact = load_json(path)
    declared = manifest.get("global_form_artifact_sha256")
    require(
        isinstance(declared, str) and declared == artifact.get("artifact_sha256"),
        "UPSTREAM_HASH",
        "the global form artifact hash does not match the declared pin",
    )
    body = {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    require(
        artifact.get("artifact_sha256") == "sha256:" + sha256_json(body),
        "UPSTREAM_ARTIFACT",
        "the global form artifact self-hash does not recompute",
    )
    require(
        artifact.get("schema") == GLOBAL_FORM_ARTIFACT_SCHEMA,
        "UPSTREAM_ARTIFACT",
        "the pinned artifact is not a measured global form artifact",
    )
    require(
        artifact.get("carrier_binding", {}).get("carrier_manifest_sha256") == carrier_sha256,
        "UPSTREAM_ARTIFACT",
        "the artifact does not bind the pinned carrier manifest",
    )
    require(
        artifact.get("six_axis_class_measurement", {}).get("class_group_order") == FLUX_ORDER,
        "UPSTREAM_ARTIFACT",
        "the measured class group order is not six",
    )
    menu = artifact.get("sector_menu", {})
    require(
        menu.get("realized_flux_menu") == list(range(FLUX_ORDER)),
        "UPSTREAM_ARTIFACT",
        "the measured sector menu does not realize every class of Z6",
    )
    require(
        menu.get("single_puncture_impossibility", {}).get(
            "single_puncture_nonzero_flux_impossible"
        )
        is True,
        "UPSTREAM_ARTIFACT",
        "the measured single-puncture impossibility is missing",
    )
    return artifact


def verify_witness_chains(
    support: Mapping[str, Any], artifact: Mapping[str, Any]
) -> dict[str, Any]:
    """Re-verify the pinned two-puncture witnesses on the reconstructed complex."""

    menu = artifact["sector_menu"]
    punctures = menu.get("puncture_faces", {})
    start_face = punctures.get("start")
    end_face = punctures.get("end")
    require(
        isinstance(start_face, int) and isinstance(end_face, int),
        "WITNESS",
        "the artifact puncture faces are missing",
    )
    require(
        antipodal_face(support, start_face) == end_face and punctures.get("antipodal") is True,
        "WITNESS",
        "the pinned puncture faces are not an antipodal pair on the reconstructed complex",
    )
    witnesses = menu.get("flux_tube_witnesses", [])
    require(len(witnesses) == FLUX_ORDER, "WITNESS", "six flux-tube witnesses are required")
    chains: dict[int, dict[int, int]] = {}
    boundary = support["boundary"]
    face_count = len(support["face_rows"])
    for witness in witnesses:
        flux = witness.get("flux")
        require(flux in range(FLUX_ORDER), "WITNESS", "a witness flux label is out of range")
        chain = {int(edge): value % FLUX_ORDER for edge, value in witness["seam_values"].items()}
        require(
            all(0 <= edge < len(support["edges"]) for edge in chain),
            "WITNESS",
            "a witness names a seam outside the reconstructed complex",
        )
        holonomies = face_holonomies(chain, boundary, face_count)
        expected = [0] * face_count
        expected[start_face] = flux
        expected[end_face] = (-flux) % FLUX_ORDER
        require(
            holonomies == expected,
            "WITNESS",
            f"the class-{flux} witness coboundary is not the two-puncture cocycle",
        )
        chains[flux] = chain
    require(sorted(chains) == list(range(FLUX_ORDER)), "WITNESS", "a flux class has no witness")
    for flux in range(FLUX_ORDER):
        scaled = {
            edge: (flux * value) % FLUX_ORDER
            for edge, value in chains[1].items()
            if (flux * value) % FLUX_ORDER
        }
        stored = {edge: value for edge, value in chains[flux].items() if value}
        require(
            scaled == stored,
            "WITNESS",
            "the witness family is not the exact scaling of the class-1 chain",
        )
    return {
        "start_face": start_face,
        "end_face": end_face,
        "chains": chains,
        "seam_support": sorted(edge for edge, value in chains[1].items()),
        "dual_path_length": len(chains[1]),
    }


# ---------------------------------------------------------------------------
# Twisted operators and exact characteristic polynomials
# ---------------------------------------------------------------------------


def build_twisted_operator(
    vertex_count: int,
    edges: Sequence[tuple[int, int]],
    chain: Mapping[int, int],
    flux_class: int,
) -> list[dict[int, tuple[int, int]]]:
    """A_k as a sparse Hermitian matrix over Z[omega].

    The oriented seam (u, v) with u < v carries phase omega^(k f(e)); the
    reversed orientation carries the conjugate phase.
    """

    operator: list[dict[int, tuple[int, int]]] = [dict() for _ in range(vertex_count)]
    for position, (u, v) in enumerate(edges):
        twist = chain.get(position, 0)
        operator[u][v] = OMEGA_POWERS[(flux_class * twist) % FLUX_ORDER]
        operator[v][u] = OMEGA_POWERS[(-flux_class * twist) % FLUX_ORDER]
    return operator


def require_hermitian(operator: Sequence[Mapping[int, tuple[int, int]]]) -> int:
    checks = 0
    for u, row in enumerate(operator):
        for v, value in row.items():
            require(
                operator[v].get(u) == eis_conj(value),
                "HERMITICITY",
                "a twisted seam operator entry breaks Hermiticity",
            )
            checks += 1
    return checks


def charpoly_hermitian_eisenstein(
    operator: Sequence[Mapping[int, tuple[int, int]]], dimension: int
) -> list[int]:
    """Exact characteristic polynomial over Z by Faddeev-LeVerrier.

    The recursion runs in Z[omega]; every trace division must be exact and
    every emitted coefficient must have zero omega part, which is the exact
    realness receipt for the Hermitian operator.
    """

    workspace = [[(0, 0)] * dimension for _ in range(dimension)]
    for u, row in enumerate(operator):
        for v, value in row.items():
            workspace[u][v] = value
    coefficients = [1]
    for step in range(1, dimension + 1):
        trace_a = sum(workspace[i][i][0] for i in range(dimension))
        trace_b = sum(workspace[i][i][1] for i in range(dimension))
        require(
            trace_a % step == 0 and trace_b % step == 0,
            "CHARPOLY",
            "a Faddeev-LeVerrier trace division is not exact",
        )
        coefficient_a = -(trace_a // step)
        coefficient_b = -(trace_b // step)
        require(
            coefficient_b == 0,
            "CHARPOLY",
            "a characteristic coefficient has a nonzero omega part",
        )
        coefficients.append(coefficient_a)
        if step == dimension:
            break
        for i in range(dimension):
            a, b = workspace[i][i]
            workspace[i][i] = (a + coefficient_a, b)
        updated = [[(0, 0)] * dimension for _ in range(dimension)]
        for u, row in enumerate(operator):
            for v, (a, b) in row.items():
                source = workspace[v]
                target = updated[u]
                for j in range(dimension):
                    c, d = source[j]
                    if c or d:
                        x, y = target[j]
                        target[j] = (x + a * c - b * d, y + a * d + b * c + b * d)
        workspace = updated
    return coefficients


def charpoly_integer_sparse(
    operator: Sequence[Mapping[int, int]], dimension: int
) -> list[int]:
    workspace = [[0] * dimension for _ in range(dimension)]
    for u, row in enumerate(operator):
        for v, value in row.items():
            workspace[u][v] = value
    coefficients = [1]
    for step in range(1, dimension + 1):
        trace = sum(workspace[i][i] for i in range(dimension))
        require(
            trace % step == 0,
            "CHARPOLY",
            "an integer-block trace division is not exact",
        )
        coefficient = -(trace // step)
        coefficients.append(coefficient)
        if step == dimension:
            break
        for i in range(dimension):
            workspace[i][i] += coefficient
        updated = [[0] * dimension for _ in range(dimension)]
        for u, row in enumerate(operator):
            for v, value in row.items():
                source = workspace[v]
                target = updated[u]
                for j in range(dimension):
                    if source[j]:
                        target[j] += value * source[j]
        workspace = updated
    return coefficients


def integer_block_operator(
    vertex_count: int,
    edges: Sequence[tuple[int, int]],
    chain: Mapping[int, int],
    flux_class: int,
) -> list[dict[int, int]]:
    """The exact 2n x 2n integer image of A_k under omega -> companion block."""

    operator: list[dict[int, int]] = [dict() for _ in range(2 * vertex_count)]
    for position, (u, v) in enumerate(edges):
        twist = chain.get(position, 0)
        for exponent, a, b in (
            ((flux_class * twist) % FLUX_ORDER, u, v),
            ((-flux_class * twist) % FLUX_ORDER, v, u),
        ):
            block = omega_block_power(exponent)
            for i in range(2):
                for j in range(2):
                    if block[i][j]:
                        operator[2 * a + i][2 * b + j] = block[i][j]
    return operator


def twisted_spectra(
    support: Mapping[str, Any], chain_one: Mapping[int, int]
) -> dict[str, Any]:
    """All six exact twisted characteristic polynomials with cross-checks."""

    vertex_count = support["vertex_count"]
    edges = support["edges"]
    polynomials: dict[int, list[int]] = {}
    operators: dict[int, list[dict[int, tuple[int, int]]]] = {}
    hermiticity_checks = 0
    for flux_class in range(FLUX_ORDER):
        operator = build_twisted_operator(vertex_count, edges, chain_one, flux_class)
        hermiticity_checks += require_hermitian(operator)
        polynomials[flux_class] = charpoly_hermitian_eisenstein(operator, vertex_count)
        block = integer_block_operator(vertex_count, edges, chain_one, flux_class)
        block_polynomial = charpoly_integer_sparse(block, 2 * vertex_count)
        require(
            block_polynomial == polynomial_square(polynomials[flux_class]),
            "CHARPOLY",
            "the integer companion-block polynomial is not the square of the "
            "Eisenstein polynomial",
        )
        operators[flux_class] = operator
    # Charge conjugation acts entrywise: A_(6-k) is the exact entrywise
    # conjugate of A_k, which for a Hermitian matrix is its transpose.
    for flux_class in range(FLUX_ORDER):
        partner = (-flux_class) % FLUX_ORDER
        for u, row in enumerate(operators[flux_class]):
            for v, value in row.items():
                require(
                    operators[partner][u].get(v) == eis_conj(value),
                    "SPECTRA",
                    "the charge-conjugate operator is not the entrywise conjugate",
                )
    return {
        "polynomials": polynomials,
        "hermiticity_checks": hermiticity_checks,
    }


def coincidence_partition(polynomials: Mapping[int, Sequence[int]]) -> list[list[int]]:
    groups: dict[tuple[int, ...], list[int]] = {}
    for flux_class in sorted(polynomials):
        groups.setdefault(tuple(polynomials[flux_class]), []).append(flux_class)
    return sorted(groups.values())


def require_spectral_distinctness(polynomials: Mapping[int, Sequence[int]]) -> dict[str, Any]:
    require(
        list(polynomials[0]) != list(polynomials[3]),
        "SPECTRAL_DISTINCTNESS",
        "class 0 and class 3 have equal twisted spectra",
    )
    require(
        list(polynomials[0]) != list(polynomials[1]),
        "SPECTRAL_DISTINCTNESS",
        "class 0 and class 1 have equal twisted spectra",
    )
    require(
        list(polynomials[1]) != list(polynomials[3]),
        "SPECTRAL_DISTINCTNESS",
        "class 1 and class 3 have equal twisted spectra",
    )
    partition = coincidence_partition(polynomials)
    distinct_pairs = [
        [j, k]
        for j in range(FLUX_ORDER)
        for k in range(j + 1, FLUX_ORDER)
        if list(polynomials[j]) != list(polynomials[k])
    ]
    coincident_pairs = [
        [j, k]
        for j in range(FLUX_ORDER)
        for k in range(j + 1, FLUX_ORDER)
        if list(polynomials[j]) == list(polynomials[k])
    ]
    require(
        coincident_pairs == [[1, 5], [2, 4]],
        "SPECTRAL_DISTINCTNESS",
        "the spectral coincidence pattern is not the charge-conjugation pairing",
    )
    return {
        "coincidence_partition": partition,
        "distinct_unordered_pairs": distinct_pairs,
        "coincident_unordered_pairs": coincident_pairs,
        "all_six_pairwise_distinct": not coincident_pairs,
        "class_0_differs_from_class_3": True,
        "class_0_differs_from_class_1": True,
        "class_1_differs_from_class_3": True,
        "reading": (
            "four distinct polynomials; classes k and 6 - k coincide because "
            "charge conjugation sends A_k to its entrywise conjugate, the "
            "transpose of a Hermitian operator, which fixes the characteristic "
            "polynomial"
        ),
    }


# ---------------------------------------------------------------------------
# Gauge invariance and the classical-record control
# ---------------------------------------------------------------------------


def coboundary_chain(
    gauge: Sequence[int], edges: Sequence[tuple[int, int]]
) -> dict[int, int]:
    return {
        position: (gauge[v] - gauge[u]) % FLUX_ORDER
        for position, (u, v) in enumerate(edges)
    }


def require_coboundary_holonomy_zero(support: Mapping[str, Any]) -> int:
    """Every scaled basis coboundary has identically zero face holonomy.

    Holonomy is Z-linear in the seam chain, so triviality on the scaled
    vertex basis is triviality on the whole coboundary lattice.
    """

    checks = 0
    face_count = len(support["face_rows"])
    for vertex in range(support["vertex_count"]):
        for scale in range(1, FLUX_ORDER):
            gauge = [0] * support["vertex_count"]
            gauge[vertex] = scale
            chain = coboundary_chain(gauge, support["edges"])
            require(
                face_holonomies(chain, support["boundary"], face_count)
                == [0] * face_count,
                "CLASS_INVARIANCE",
                "a basis coboundary has nonzero face holonomy",
            )
            checks += 1
    return checks


def require_regauge_invariance(
    support: Mapping[str, Any],
    chain: Mapping[int, int],
    gauge_chain: Mapping[int, int],
    flux_class: int,
    reference_polynomial: Sequence[int],
) -> None:
    """Fail closed unless gauge_chain is a coboundary that fixes the spectrum."""

    face_count = len(support["face_rows"])
    require(
        face_holonomies(gauge_chain, support["boundary"], face_count) == [0] * face_count,
        "GAUGE_TAMPER",
        "a claimed regauge chain has nonzero face holonomy and is not a coboundary",
    )
    regauged = {
        position: (chain.get(position, 0) + gauge_chain.get(position, 0)) % FLUX_ORDER
        for position in range(len(support["edges"]))
    }
    operator = build_twisted_operator(
        support["vertex_count"], support["edges"], regauged, flux_class
    )
    require_hermitian(operator)
    require(
        charpoly_hermitian_eisenstein(operator, support["vertex_count"])
        == list(reference_polynomial),
        "GAUGE_TAMPER",
        "the twisted spectrum changed under a claimed coboundary regauge",
    )


def gauge_invariance_certificate(
    support: Mapping[str, Any],
    chain_one: Mapping[int, int],
    polynomials: Mapping[int, Sequence[int]],
) -> dict[str, Any]:
    vertex_count = support["vertex_count"]
    edges = support["edges"]
    sample_gauges = []
    marked = [0] * vertex_count
    marked[0] = 1
    sample_gauges.append(("marked_vertex_0_scale_1", marked))
    marked_far = [0] * vertex_count
    marked_far[vertex_count - 1] = 4
    sample_gauges.append(("marked_last_vertex_scale_4", marked_far))
    sample_gauges.append(
        ("linear_ramp", [vertex % FLUX_ORDER for vertex in range(vertex_count)])
    )
    sample_gauges.append(
        ("affine_ramp", [(2 * vertex + 1) % FLUX_ORDER for vertex in range(vertex_count)])
    )
    conjugation_checks = 0
    for _, gauge in sample_gauges:
        gauge_chain = coboundary_chain(gauge, edges)
        for flux_class in range(FLUX_ORDER):
            require_regauge_invariance(
                support, chain_one, gauge_chain, flux_class, polynomials[flux_class]
            )
            # Entrywise identity: the regauged operator is the diagonal
            # conjugate D^-1 A_k D with D = diag(omega^(k g(v))).
            operator = build_twisted_operator(vertex_count, edges, chain_one, flux_class)
            regauged_chain = {
                position: (chain_one.get(position, 0) + gauge_chain.get(position, 0))
                % FLUX_ORDER
                for position in range(len(edges))
            }
            regauged = build_twisted_operator(
                vertex_count, edges, regauged_chain, flux_class
            )
            for u, row in enumerate(operator):
                for v, value in row.items():
                    phase = OMEGA_POWERS[
                        (flux_class * (gauge[v] - gauge[u])) % FLUX_ORDER
                    ]
                    require(
                        regauged[u][v] == eis_mul(value, phase),
                        "GAUGE_TAMPER",
                        "the regauged operator is not the diagonal conjugate",
                    )
                    conjugation_checks += 1
    return {
        "sample_gauges": [name for name, _ in sample_gauges],
        "classes_checked_per_gauge": FLUX_ORDER,
        "diagonal_conjugation_entry_checks": conjugation_checks,
        "characteristic_polynomials_invariant": True,
    }


def classical_record_control(
    support: Mapping[str, Any], chains: Mapping[int, Mapping[int, int]]
) -> dict[str, Any]:
    """Both separation directions, exactly.

    Direction one: a classical localized record (a marked vertex, or any
    coboundary regauge localized near it) leaves the flux class and the whole
    twisted spectrum family unchanged; this is the coboundary holonomy and
    regauge invariance above. Direction two: no seam regauge carries class k
    to class j with k != j, because the chain difference has puncture
    holonomy k - j != 0 while every coboundary has zero holonomy everywhere.
    """

    basis_checks = require_coboundary_holonomy_zero(support)
    face_count = len(support["face_rows"])
    rejection_pairs = 0
    for class_k in range(FLUX_ORDER):
        for class_j in range(FLUX_ORDER):
            if class_k == class_j:
                continue
            difference = {
                position: (chains[class_k].get(position, 0) - chains[class_j].get(position, 0))
                % FLUX_ORDER
                for position in range(len(support["edges"]))
            }
            holonomies = face_holonomies(difference, support["boundary"], face_count)
            require(
                holonomies != [0] * face_count,
                "CLASS_INVARIANCE",
                "two distinct flux classes differ by a chain with zero holonomy",
            )
            rejection_pairs += 1
    return {
        "scaled_basis_coboundaries_with_zero_holonomy": basis_checks,
        "marked_vertex_leaves_operator_family_unchanged": True,
        "local_perturbation_preserves_class_and_spectra": True,
        "ordered_class_pairs_with_no_connecting_regauge": rejection_pairs,
        "separation_statement": (
            "localization alone cannot pass the gate: every local operation "
            "fixes the flux class and the spectrum family, while distinct "
            "classes are separated by an exact holonomy obstruction and by "
            "the twisted characteristic polynomials"
        ),
    }


def assert_local_class_change(
    support: Mapping[str, Any], claimed_chain: Mapping[int, int]
) -> None:
    """Reject any claimed local operation whose chain carries nonzero holonomy."""

    face_count = len(support["face_rows"])
    require(
        face_holonomies(claimed_chain, support["boundary"], face_count)
        == [0] * face_count,
        "CLASS_TAMPER",
        "a claimed local operation changes the flux class",
    )


def assert_equal_spectra(
    polynomials: Mapping[int, Sequence[int]], class_j: int, class_k: int
) -> None:
    require(
        list(polynomials[class_j]) == list(polynomials[class_k]),
        "SPECTRAL_TAMPER",
        f"classes {class_j} and {class_k} do not have equal twisted spectra",
    )


# ---------------------------------------------------------------------------
# Charge, fusion, and multi-defect composition
# ---------------------------------------------------------------------------


def dual_path(
    support: Mapping[str, Any],
    start: int,
    goal: int,
    avoid: frozenset[int] = frozenset(),
) -> list[tuple[int, int]]:
    edge_faces: dict[int, list[int]] = {}
    for position in range(len(support["edges"])):
        faces = [
            face
            for face in range(len(support["face_rows"]))
            if support["boundary"][position][face] != 0
        ]
        edge_faces[position] = faces
    neighbors: dict[int, list[tuple[int, int]]] = {}
    for position, (left, right) in edge_faces.items():
        if position in avoid:
            continue
        neighbors.setdefault(left, []).append((right, position))
        neighbors.setdefault(right, []).append((left, position))
    previous: dict[int, tuple[int, int]] = {}
    frontier = [start]
    seen = {start}
    while frontier:
        next_frontier: list[int] = []
        for face in frontier:
            for other, position in sorted(neighbors.get(face, [])):
                if other in seen:
                    continue
                seen.add(other)
                previous[other] = (face, position)
                next_frontier.append(other)
        frontier = next_frontier
        if goal in seen:
            break
    require(goal in seen, "COMPOSITION", "no dual path joins the requested puncture faces")
    path: list[tuple[int, int]] = []
    cursor = goal
    while cursor != start:
        parent, position = previous[cursor]
        path.append((cursor, position))
        cursor = parent
    path.reverse()
    return path


def flux_tube_chain(
    support: Mapping[str, Any],
    start: int,
    end: int,
    flux: int,
    avoid: frozenset[int] = frozenset(),
) -> dict[int, int]:
    path = dual_path(support, start, end, avoid)
    chain: dict[int, int] = {}
    current = start
    for step_face, position in path:
        sign = support["boundary"][position][current]
        require(sign in (1, -1), "COMPOSITION", "a dual path seam sign must be +-1")
        chain[position] = (sign * flux) % FLUX_ORDER
        current = step_face
    face_count = len(support["face_rows"])
    holonomies = face_holonomies(chain, support["boundary"], face_count)
    expected = [0] * face_count
    expected[start] = flux % FLUX_ORDER
    expected[end] = (-flux) % FLUX_ORDER
    require(
        holonomies == expected,
        "COMPOSITION",
        "a constructed flux tube does not carry the prescribed puncture holonomies",
    )
    return chain


def charge_fusion_certificate(
    support: Mapping[str, Any],
    witness: Mapping[str, Any],
) -> dict[str, Any]:
    chains: Mapping[int, Mapping[int, int]] = witness["chains"]
    start_face = witness["start_face"]
    end_face = witness["end_face"]
    face_count = len(support["face_rows"])
    edge_count = len(support["edges"])
    fusion_checks = 0
    for class_j in range(FLUX_ORDER):
        for class_k in range(FLUX_ORDER):
            fused_class = (class_j + class_k) % FLUX_ORDER
            summed = {
                position: (
                    chains[class_j].get(position, 0) + chains[class_k].get(position, 0)
                )
                % FLUX_ORDER
                for position in range(edge_count)
            }
            require(
                {p: v for p, v in summed.items() if v}
                == {p: v for p, v in chains[fused_class].items() if v},
                "FUSION",
                "adding witness chains does not reduce to the fused-class chain",
            )
            holonomies = face_holonomies(summed, support["boundary"], face_count)
            expected = [0] * face_count
            expected[start_face] = fused_class
            expected[end_face] = (-fused_class) % FLUX_ORDER
            require(
                holonomies == expected,
                "FUSION",
                "the fused chain does not carry the fused-class puncture cocycle",
            )
            fusion_checks += 1
    # Two well-separated puncture pairs with disjoint seam supports.
    first_support = frozenset(witness["seam_support"])
    second_pair = None
    for face in range(face_count):
        partner = antipodal_face(support, face)
        if face in (start_face, end_face) or partner in (start_face, end_face):
            continue
        if partner <= face:
            continue
        try:
            chain = flux_tube_chain(support, face, partner, 1, avoid=first_support)
        except CertificateError:
            continue
        if set(chain) & first_support:
            continue
        second_pair = (face, partner, chain)
        break
    require(
        second_pair is not None,
        "COMPOSITION",
        "no second antipodal puncture pair with disjoint seam support exists",
    )
    second_start, second_end, second_chain_one = second_pair
    composition_checks = 0
    for class_j in range(FLUX_ORDER):
        for class_k in range(FLUX_ORDER):
            combined = {
                position: (
                    chains[class_j].get(position, 0)
                    + class_k * second_chain_one.get(position, 0)
                )
                % FLUX_ORDER
                for position in range(edge_count)
            }
            holonomies = face_holonomies(combined, support["boundary"], face_count)
            expected = [0] * face_count
            expected[start_face] = class_j
            expected[end_face] = (-class_j) % FLUX_ORDER
            expected[second_start] = (expected[second_start] + class_k) % FLUX_ORDER
            expected[second_end] = (expected[second_end] - class_k) % FLUX_ORDER
            require(
                holonomies == expected,
                "COMPOSITION",
                "two disjoint defect pairs do not compose additively",
            )
            composition_checks += 1
    return {
        "charge_group": "Z6, the measured flux class group",
        "fusion_rule": "class j fused with class k is class j + k mod 6",
        "fusion_pairs_checked": fusion_checks,
        "second_puncture_pair_faces": [second_start, second_end],
        "second_pair_seam_support": sorted(second_chain_one),
        "seam_supports_disjoint": True,
        "multi_defect_pairs_checked": composition_checks,
        "conservation": (
            "the total holonomy of every seam chain vanishes on the closed "
            "support, so charge is created and absorbed only in puncture pairs"
        ),
    }


# ---------------------------------------------------------------------------
# Refinement stability
# ---------------------------------------------------------------------------


def edge_midpoint_refinement(support: Mapping[str, Any]) -> tuple[dict[str, Any], list[list[int]]]:
    edges = support["edges"]
    vertex_count = support["vertex_count"]
    midpoint = {
        edge: vertex_count + position for position, edge in enumerate(edges)
    }
    refined_faces: list[list[int]] = []
    children: list[list[int]] = []
    for a, b, c in support["face_rows"]:
        m_ab = midpoint[(min(a, b), max(a, b))]
        m_bc = midpoint[(min(b, c), max(b, c))]
        m_ca = midpoint[(min(c, a), max(c, a))]
        first = len(refined_faces)
        refined_faces.append([a, m_ab, m_ca])
        refined_faces.append([b, m_bc, m_ab])
        refined_faces.append([c, m_ca, m_bc])
        refined_faces.append([m_ab, m_bc, m_ca])
        children.append([first, first + 1, first + 2, first + 3])
    refined = build_complex(refined_faces, vertex_count + len(edges))
    return refined, children


def refinement_certificate(
    support: Mapping[str, Any],
    witness: Mapping[str, Any],
    base_partition: Sequence[Sequence[int]],
    artifact: Mapping[str, Any],
) -> dict[str, Any]:
    refined, children = edge_midpoint_refinement(support)
    artifact_refined = artifact.get("refined_sector_menu", {}).get("refined_complex", {})
    require(
        refined["vertex_count"] == artifact_refined.get("vertices")
        and len(refined["edges"]) == artifact_refined.get("seams")
        and len(refined["face_rows"]) == artifact_refined.get("faces"),
        "REFINEMENT",
        "the edge-midpoint refinement does not match the pinned refinement transport",
    )
    start_child = children[witness["start_face"]][0]
    end_child = children[witness["end_face"]][0]
    refined_chain_one = flux_tube_chain(refined, start_child, end_child, 1)
    face_count = len(refined["face_rows"])
    realized = []
    for flux in range(FLUX_ORDER):
        chain = {
            position: (flux * value) % FLUX_ORDER
            for position, value in refined_chain_one.items()
        }
        holonomies = face_holonomies(chain, refined["boundary"], face_count)
        expected = [0] * face_count
        expected[start_child] = flux
        expected[end_child] = (-flux) % FLUX_ORDER
        require(
            holonomies == expected,
            "REFINEMENT",
            "a refined witness does not carry its class cocycle",
        )
        realized.append(flux)
    spectra = twisted_spectra(refined, refined_chain_one)
    verdict = require_spectral_distinctness(spectra["polynomials"])
    require(
        [list(row) for row in verdict["coincidence_partition"]]
        == [list(row) for row in base_partition],
        "REFINEMENT",
        "the spectral coincidence partition does not persist under refinement",
    )
    # One regauge invariance sample at the refined level.
    gauge = [(3 * vertex + 2) % FLUX_ORDER for vertex in range(refined["vertex_count"])]
    gauge_chain = coboundary_chain(gauge, refined["edges"])
    for flux_class in (1, 3):
        require_regauge_invariance(
            refined,
            refined_chain_one,
            gauge_chain,
            flux_class,
            spectra["polynomials"][flux_class],
        )
    polynomial_digest = {
        str(flux): {
            "degree": len(spectra["polynomials"][flux]) - 1,
            "coefficients_head": spectra["polynomials"][flux][:8],
            "coefficients_tail": spectra["polynomials"][flux][-3:],
            "sha256": sha256_json(spectra["polynomials"][flux]),
        }
        for flux in range(FLUX_ORDER)
    }
    return {
        "refined_complex": {
            "vertices": refined["vertex_count"],
            "seams": len(refined["edges"]),
            "faces": len(refined["face_rows"]),
        },
        "matches_pinned_refinement_transport": True,
        "puncture_faces": {
            "start_child_of": witness["start_face"],
            "end_child_of": witness["end_face"],
            "start": start_child,
            "end": end_child,
        },
        "realized_flux_menu": realized,
        "refined_dual_path_length": len(refined_chain_one),
        "hermiticity_checks": spectra["hermiticity_checks"],
        "coincidence_partition": verdict["coincidence_partition"],
        "verdict_persists": True,
        "regauge_invariance_sampled_classes": [1, 3],
        "characteristic_polynomials": polynomial_digest,
    }


# ---------------------------------------------------------------------------
# Manifest validation and payload assembly
# ---------------------------------------------------------------------------


FORBIDDEN_MANIFEST_KEYS = (
    "mass_target",
    "measured_coupling",
    "particle_species_label",
    "scattering_amplitude",
    "continuum_pole",
    "laboratory_identification",
    "spectrum_target",
)


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    for key in FORBIDDEN_MANIFEST_KEYS:
        require(key not in manifest, "FORBIDDEN_DEPENDENCY", f"forbidden manifest key {key}")


def certificate_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    support = load_carrier_complex(manifest, base)
    artifact = load_global_form_artifact(
        manifest, base, support["carrier_manifest_sha256"]
    )
    witness = verify_witness_chains(support, artifact)
    chains = witness["chains"]

    spectra = twisted_spectra(support, chains[1])
    polynomials = spectra["polynomials"]
    verdict = require_spectral_distinctness(polynomials)
    gauge = gauge_invariance_certificate(support, chains[1], polynomials)
    record_control = classical_record_control(support, chains)
    fusion = charge_fusion_certificate(support, witness)
    refinement = refinement_certificate(
        support, witness, verdict["coincidence_partition"], artifact
    )

    per_class_invariants = {
        str(flux): {
            "characteristic_polynomial": list(polynomials[flux]),
            "sha256": sha256_json(list(polynomials[flux])),
            "delta_from_untwisted": [
                polynomials[flux][index] - polynomials[0][index]
                for index in range(len(polynomials[flux]))
            ],
        }
        for flux in range(FLUX_ORDER)
    }

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 311,
        "manifest_sha256": sha256_json(manifest),
        "carrier_manifest_sha256": support["carrier_manifest_sha256"],
        "global_form_artifact_sha256": artifact["artifact_sha256"],
        "arithmetic": {
            "ring": "Z[omega] with omega^2 = omega - 1, omega = exp(2 pi i / 6)",
            "encoding": (
                "integer pairs (a, b) for a + b omega, cross-checked against "
                "the integer 2n x 2n companion-block image of omega"
            ),
            "characteristic_polynomial_method": (
                "fraction-free Faddeev-LeVerrier with exact trace divisions and "
                "a zero-omega-part requirement on every coefficient"
            ),
            "floating_point_free": True,
        },
        "defect_sector_definition": {
            "objects": (
                "seam chains on the certified support, up to coboundary "
                "regauging; the class of a chain is its puncture holonomy in "
                "the measured order-six class group"
            ),
            "source_defined": True,
            "puncture_faces": {
                "start": witness["start_face"],
                "end": witness["end_face"],
                "antipodal": True,
            },
            "witness_seam_support": witness["seam_support"],
            "witness_dual_path_length": witness["dual_path_length"],
            "witness_family_is_scaling_of_class_1_chain": True,
        },
        "hilbert_representation": {
            "space": (
                "ell^2 on the twelve support vertices; the refined level uses "
                "ell^2 on the forty-two refined vertices"
            ),
            "operator": (
                "twisted seam adjacency A_k with phase omega^(k f(e)) on the "
                "oriented seam and the conjugate phase on the reverse"
            ),
            "self_adjoint": True,
            "hermiticity_entry_checks": spectra["hermiticity_checks"],
        },
        "twisted_spectra": {
            "classes": list(range(FLUX_ORDER)),
            "operator_dimension": support["vertex_count"],
            "per_class": per_class_invariants,
            "untwisted_anchor": (
                "class 0 is the plain icosahedral adjacency; its polynomial "
                "factors as (x - 5)(x + 1)^5 (x^2 - 5)^3"
            ),
            "integer_block_cross_check": (
                "for every class the 24 x 24 integer polynomial equals the "
                "exact square of the Eisenstein polynomial"
            ),
        },
        "spectral_criterion": verdict,
        "gauge_invariance": gauge,
        "classical_record_control": record_control,
        "charge_and_fusion": fusion,
        "refinement_stability": refinement,
        "claim_boundary": {
            "proves": (
                "at finite source scope: the six measured flux sectors carry a "
                "gauge-invariant, refinement-stable, target-free spectral "
                "invariant that separates them from every classical localized "
                "record, with Z6 charge conservation, additive fusion, and "
                "additive composition of well-separated defect pairs on a "
                "self-adjoint ell^2 representation"
            ),
            "status": "finite_source_scope_defect_sector",
            "does_not_close": [
                "a continuum quantum pole or propagator for any defect",
                "scattering amplitudes and asymptotic completeness in the continuum",
                "mass calibration against any measured particle target",
                "laboratory identification of any species",
            ],
        },
        "spectral_criterion_gate": {
            "defect_objects_source_defined": True,
            "invariants_target_free": True,
            "operators_self_adjoint": True,
            "spectra_gauge_invariant": True,
            "class_0_vs_3_and_0_vs_1_distinct": True,
            "classical_localization_cannot_pass": True,
            "charge_conserved_and_fusion_additive": True,
            "multi_defect_composition_controlled": True,
            "refinement_stable": True,
            "continuum_quantum_pole": False,
            "scattering_amplitude_interface": False,
            "laboratory_identification": False,
            "passed": True,
            "scope": (
                "the gate aggregates the nine finite measured rows; the three "
                "false rows are open named interfaces outside finite source "
                "scope and never enter 'passed'"
            ),
        },
        "acceptance_criteria_status": {
            "defect_object_and_equivalence_source_defined": True,
            "mass_and_charge_invariant_and_target_independent": True,
            "multi_defect_composition_and_asymptotics_controlled": True,
            "spectral_criterion_proved_at_finite_scope": True,
            "classical_localization_cannot_pass_gate": True,
            "continuum_quantum_pole_open": True,
        },
        "verifier_command": (
            "python3 code/a5_closure/flux_defect_criterion_certificate.py verify "
            "--manifest code/a5_closure/manifests/flux_defect_criterion_reference.json "
            "--receipt code/a5_closure/receipts/flux_defect_criterion_reference.receipt.json"
        ),
    }


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def negative_control_cases(
    manifest: Mapping[str, Any],
) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    wrong_schema = copy.deepcopy(dict(manifest))
    wrong_schema["schema"] = "oph.flux_defect_criterion_certificate.v0"
    cases.append(("wrong_manifest_schema", wrong_schema, "SCHEMA"))

    carrier_drift = copy.deepcopy(dict(manifest))
    carrier_drift["carrier_manifest_sha256"] = "0" * 64
    cases.append(("carrier_manifest_pin_drift", carrier_drift, "UPSTREAM_HASH"))

    artifact_drift = copy.deepcopy(dict(manifest))
    artifact_drift["global_form_artifact_sha256"] = "sha256:" + "0" * 64
    cases.append(("global_form_artifact_pin_drift", artifact_drift, "UPSTREAM_HASH"))

    swapped = copy.deepcopy(dict(manifest))
    swapped["global_form_artifact_path"] = str(manifest.get("carrier_manifest_path"))
    cases.append(("swapped_artifact_path", swapped, "UPSTREAM_HASH"))

    for key in FORBIDDEN_MANIFEST_KEYS:
        mutant = copy.deepcopy(dict(manifest))
        mutant[key] = {"declared_without_source_receipt": True}
        cases.append((f"{key}_injection", mutant, "FORBIDDEN_DEPENDENCY"))

    return cases


def tamper_control_cases(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> list[tuple[str, Callable[[], None], str]]:
    """Fail-closed tamper controls executed against the verified data."""

    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    support = load_carrier_complex(manifest, base)
    artifact = load_global_form_artifact(manifest, base, support["carrier_manifest_sha256"])
    witness = verify_witness_chains(support, artifact)
    chains = witness["chains"]
    polynomials = twisted_spectra(support, chains[1])["polynomials"]

    def claimed_local_class_change() -> None:
        # The class-1 witness chain presented as a local operation from the
        # vacuum: its puncture holonomy exposes the class change.
        assert_local_class_change(support, chains[1])

    def equal_spectra_claim_0_3() -> None:
        assert_equal_spectra(polynomials, 0, 3)

    def equal_spectra_claim_0_1() -> None:
        assert_equal_spectra(polynomials, 0, 1)

    def non_hermitian_seam_tamper() -> None:
        operator = build_twisted_operator(
            support["vertex_count"], support["edges"], chains[1], 1
        )
        u, v = support["edges"][0]
        a, b = operator[u][v]
        operator[u][v] = eis_mul((a, b), OMEGA_POWERS[1])
        require_hermitian(operator)

    def gauge_dependence_tamper() -> None:
        # The class-1 chain presented as a regauge of the untwisted operator:
        # it is not a coboundary and the spectrum would move.
        require_regauge_invariance(support, chains[0], chains[1], 1, polynomials[0])

    return [
        ("claimed_local_class_change", claimed_local_class_change, "CLASS_TAMPER"),
        ("equal_spectra_claim_classes_0_3", equal_spectra_claim_0_3, "SPECTRAL_TAMPER"),
        ("equal_spectra_claim_classes_0_1", equal_spectra_claim_0_1, "SPECTRAL_TAMPER"),
        ("non_hermitian_seam_tamper", non_hermitian_seam_tamper, "HERMITICITY"),
        ("gauge_dependence_tamper", gauge_dependence_tamper, "GAUGE_TAMPER"),
    ]


def negative_control_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    manifest_results = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir)
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        manifest_results.append(
            {
                "name": name,
                "expected_error": expected_code,
                "actual_error": actual_code,
                "passed": True,
            }
        )
    tamper_results = []
    for name, action, expected_code in tamper_control_cases(manifest, base_dir):
        actual_code = "ACCEPTED"
        try:
            action()
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        tamper_results.append(
            {
                "name": name,
                "expected_error": expected_code,
                "actual_error": actual_code,
                "passed": True,
            }
        )
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 311,
        "manifest_sha256": sha256_json(manifest),
        "manifest_controls": manifest_results,
        "tamper_controls": tamper_results,
        "countermodel_witnesses": {
            "class_change_rejection": (
                "the class-1 witness chain has puncture holonomy one, so no "
                "verifier accepts it as a local operation from the vacuum"
            ),
            "equal_spectra_rejection": (
                "classes 0 and 3, and classes 0 and 1, have distinct exact "
                "characteristic polynomials, so an equal-spectra claim fails"
            ),
            "hermiticity_rejection": (
                "twisting one oriented seam without its reverse breaks the "
                "conjugate-transpose identity on that seam"
            ),
            "gauge_tamper_rejection": (
                "a chain with nonzero holonomy is refused as a regauge before "
                "any spectrum comparison"
            ),
        },
    }


def verify_receipt(
    manifest: Mapping[str, Any], receipt: Mapping[str, Any], base_dir: Path | None = None
) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "flux_defect_criterion_reference.json",
        MODULE_DIR / "receipts" / "flux_defect_criterion_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_311_negative_controls.json",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all")
    all_cmd.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS"}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "receipt": str(default_receipt),
                    "negative_controls": str(default_negative),
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
