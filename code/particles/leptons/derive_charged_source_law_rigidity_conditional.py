#!/usr/bin/env python3
"""Verify the conditional CFQ packet-weight rigidity theorem.

The submitted CFQ construction stipulates eight finite response registers and
eight primitive path classes.  Conditional on that entire packet, normalized
rank-one traces fix the displayed scalar coefficients and exclude independent
source multipliers.  This module verifies that implication exactly and checks
that it reproduces the already frozen, historically target-informed affine
map.

It does *not* promote the packet to an OPH consequence.  In particular, the
physical response modules, record-algebra bridge, clock, path exhaustion,
signs, refinement intertwiners, and no-target ancestry have no supplied
evidence.  Every physical CFQ gate therefore remains false.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import mpmath as mp


mp.mp.dps = 100

HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
FACE_RECEIPT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_icosahedral_face_carrier_frontier.json"
)
DECLARED_MAP_RECEIPT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_face_incidence_conditional_theorem.json"
)
DEFAULT_OUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_source_law_rigidity_conditional.json"
)

RIGIDITY_ARCHIVE_SHA256 = (
    "792dad5c1524984c52d7582e751e1123099978164def80fad34310a34e314dff"
)
RIGIDITY_NARRATIVE_SHA256 = (
    "879f5afbdb704d8894fe826544e26dff90651040df45954485be422ea2c0779e"
)
RIGIDITY_SUBMITTED_RECEIPT_SHA256 = (
    "d0acb1242714a38ea7ccc3b9f335bc9accf8e8bbdd9d6959d7e7e5dcf93be4d3"
)


@dataclass(frozen=True)
class Register:
    name: str
    dimension: int
    event_rank: int = 1

    @property
    def trace_weight(self) -> Fraction:
        return Fraction(self.event_rank, self.dimension)


@dataclass(frozen=True)
class PathRule:
    name: str
    block: str
    monomial: str
    sign: int
    registers: tuple[str, ...]
    multiplicity: int = 1


Graph = dict[str, set[str]]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def mp_text(value: mp.mpf, digits: int = 80) -> str:
    return mp.nstr(value, digits)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def add_node(graph: Graph, node: str) -> None:
    graph.setdefault(node, set())


def add_edge(graph: Graph, left: str, right: str) -> None:
    if left == right:
        return
    add_node(graph, left)
    add_node(graph, right)
    graph[left].add(right)
    graph[right].add(left)


def graph_connected(graph: Graph) -> bool:
    if not graph:
        return False
    pending = [next(iter(graph))]
    seen: set[str] = set()
    while pending:
        node = pending.pop()
        if node in seen:
            continue
        seen.add(node)
        pending.extend(graph[node] - seen)
    return len(seen) == len(graph)


def graph_summary(graph: Graph, expected_nodes: int) -> dict[str, Any]:
    return {
        "nodes": len(graph),
        "edges": sum(len(neighbors) for neighbors in graph.values()) // 2,
        "connected": graph_connected(graph),
        "expected_nodes": expected_nodes,
        "dimension_matches": len(graph) == expected_nodes,
    }


def star_graph(items: Iterable[str]) -> Graph:
    graph: Graph = {}
    add_node(graph, "identity")
    for item in items:
        add_edge(graph, "identity", item)
    return graph


def cartesian_product(left: Graph, right: Graph) -> Graph:
    product: Graph = {}
    for left_node in left:
        for right_node in right:
            add_node(product, f"{left_node}|{right_node}")
    for left_node, neighbors in left.items():
        for left_neighbor in neighbors:
            if left_node < left_neighbor:
                for right_node in right:
                    add_edge(
                        product,
                        f"{left_node}|{right_node}",
                        f"{left_neighbor}|{right_node}",
                    )
    for right_node, neighbors in right.items():
        for right_neighbor in neighbors:
            if right_node < right_neighbor:
                for left_node in left:
                    add_edge(
                        product,
                        f"{left_node}|{right_node}",
                        f"{left_node}|{right_neighbor}",
                    )
    return product


def build_register_graphs(
    face_vertex_indices: list[list[int]],
    register_dimensions: dict[str, int],
) -> tuple[dict[str, Graph], dict[str, Any]]:
    """Construct the bundle's finite graph models without external graph code.

    Connectivity shows that each *declared* graph can generate a full matrix
    algebra from local matrix units.  It is not evidence that the physical OPH
    response realizes that graph or those matrix units.
    """

    faces = [tuple(sorted(face)) for face in face_vertex_indices]
    vertices = sorted({vertex for face in faces for vertex in face})
    edges = sorted(
        {
            tuple(sorted((face[left], face[right])))
            for face in faces
            for left, right in ((0, 1), (1, 2), (0, 2))
        }
    )
    adjacency: dict[int, set[int]] = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)

    chart_face = faces[0]
    chart_face_edges = {
        tuple(sorted((chart_face[0], chart_face[1]))),
        tuple(sorted((chart_face[1], chart_face[2]))),
        tuple(sorted((chart_face[0], chart_face[2]))),
    }

    kappa_0: Graph = {}
    for edge in edges:
        add_node(kappa_0, f"edge:{edge[0]}-{edge[1]}")
    for index, face in enumerate(faces):
        face_node = f"face:{index}"
        add_node(kappa_0, face_node)
        for edge in edges:
            if set(edge).issubset(face):
                add_edge(kappa_0, face_node, f"edge:{edge[0]}-{edge[1]}")

    transverse_vertices = [vertex for vertex in vertices if vertex not in chart_face]
    kappa_1 = star_graph(f"edge:{left}-{right}" for left, right in edges)
    kappa_2 = star_graph(f"vertex:{vertex}" for vertex in transverse_vertices)

    chi_0: Graph = {}
    binary_dimension = len(transverse_vertices)
    for record in range(2**binary_dimension):
        add_node(chi_0, str(record))
    for record in range(2**binary_dimension):
        for bit in range(binary_dimension):
            neighbor = record ^ (1 << bit)
            if record < neighbor:
                add_edge(chi_0, str(record), str(neighbor))

    anchor = chart_face[0]
    collar_vertices = sorted(adjacency[anchor])
    local7: Graph = {}
    for vertex in collar_vertices:
        add_node(local7, f"collar:{vertex}")
    for left in collar_vertices:
        for right in adjacency[left]:
            if right in collar_vertices and left < right:
                add_edge(local7, f"collar:{left}", f"collar:{right}")
    for closure in range(2):
        for vertex in collar_vertices:
            add_edge(local7, f"closure:{closure}", f"collar:{vertex}")

    compare11: Graph = {}
    nonanchor_vertices = [vertex for vertex in vertices if vertex != anchor]
    for vertex in nonanchor_vertices:
        add_node(compare11, f"vertex:{vertex}")
    for left, right in edges:
        if anchor not in (left, right):
            add_edge(compare11, f"vertex:{left}", f"vertex:{right}")
    chi_1 = cartesian_product(local7, compare11)

    zeta_0 = star_graph(f"face:{index}" for index in range(len(faces)))
    nonface_edges = [edge for edge in edges if edge not in chart_face_edges]
    zeta_1: Graph = {}
    for edge in nonface_edges:
        add_node(zeta_1, f"edge:{edge[0]}-{edge[1]}")
    for index, left in enumerate(nonface_edges):
        for right in nonface_edges[index + 1 :]:
            if set(left) & set(right):
                add_edge(
                    zeta_1,
                    f"edge:{left[0]}-{left[1]}",
                    f"edge:{right[0]}-{right[1]}",
                )

    zeta_2: Graph = {}
    for vertex in collar_vertices:
        add_node(zeta_2, f"collar:{vertex}")
    for left in collar_vertices:
        for right in adjacency[left]:
            if right in collar_vertices and left < right:
                add_edge(zeta_2, f"collar:{left}", f"collar:{right}")

    graphs = {
        "kappa_0": kappa_0,
        "kappa_1": kappa_1,
        "kappa_2": kappa_2,
        "chi_0": chi_0,
        "chi_1": chi_1,
        "zeta_0": zeta_0,
        "zeta_1": zeta_1,
        "zeta_2": zeta_2,
    }
    geometry = {
        "selected_face_for_audit_only": list(chart_face),
        "selected_anchor_for_audit_only": anchor,
        "face_choice_is_physical": False,
        "icosahedral_face_count": len(faces),
        "transverse_vertex_count": len(transverse_vertices),
        "nonface_edge_count": len(nonface_edges),
        "collar_vertex_count": len(collar_vertices),
    }
    summaries = {
        name: graph_summary(graph, register_dimensions[name])
        for name, graph in graphs.items()
    }
    return graphs, {"geometry": geometry, "graphs": summaries}


def packet_registers() -> dict[str, Register]:
    return {
        "kappa_0": Register("kappa_0", 50),
        "kappa_1": Register("kappa_1", 31),
        "kappa_2": Register("kappa_2", 10),
        "chi_0": Register("chi_0", 512),
        "chi_1": Register("chi_1", 77),
        "zeta_0": Register("zeta_0", 21),
        "zeta_1": Register("zeta_1", 27),
        "zeta_2": Register("zeta_2", 5),
    }


def packet_paths() -> tuple[PathRule, ...]:
    return (
        PathRule("kappa_write", "kappa", "alpha_U", +1, ("kappa_0",)),
        PathRule(
            "kappa_return_1", "kappa", "alpha_U*kappa", -1, ("kappa_1",)
        ),
        PathRule(
            "kappa_return_2",
            "kappa",
            "alpha_U^2*kappa",
            -1,
            ("kappa_1", "kappa_2"),
        ),
        PathRule("chi_write", "chi", "alpha_U^2", +1, ("chi_0",)),
        PathRule("chi_return_1", "chi", "alpha_U*chi", +1, ("chi_1",)),
        PathRule("zeta_write", "zeta", "alpha_P^2", +1, ("zeta_0",)),
        PathRule("zeta_return_1", "zeta", "alpha_U*zeta", +1, ("zeta_1",)),
        PathRule(
            "zeta_return_2",
            "zeta",
            "alpha_U^2*zeta",
            +1,
            ("zeta_1", "zeta_2"),
        ),
    )


def path_weight(path: PathRule, registers: dict[str, Register]) -> Fraction:
    weight = Fraction(path.sign * path.multiplicity, 1)
    for register_name in path.registers:
        weight *= registers[register_name].trace_weight
    return weight


def all_path_weights(
    paths: tuple[PathRule, ...], registers: dict[str, Register]
) -> dict[str, Fraction]:
    return {path.name: path_weight(path, registers) for path in paths}


def multiplier_admissible_inside_cfq(multiplier: Fraction) -> bool:
    """Return admissibility only for the fixed rank-one, unit-clock packet."""

    return multiplier == 1


def build_artifact(
    face_receipt: dict[str, Any],
    face_receipt_sha256: str,
    declared_map_receipt: dict[str, Any],
    declared_map_receipt_sha256: str,
) -> dict[str, Any]:
    if face_receipt.get("checks_pass") is not True:
        raise ValueError("the exact icosahedral face receipt is required")
    if face_receipt.get("public_charged_mass_promotion_allowed") is not False:
        raise ValueError("the face receipt must retain its promotion boundary")
    if declared_map_receipt.get("checks_pass") is not True:
        raise ValueError("the declared-map conditional theorem is required")
    if declared_map_receipt.get("historical_charged_target_informed") is not True:
        raise ValueError("historical target ancestry must not be erased")
    if declared_map_receipt.get("branch_tuple_coherent") is not False:
        raise ValueError("the hybrid-branch warning must remain active")
    if declared_map_receipt.get("public_prediction_promotion_allowed") is not False:
        raise ValueError("the declared map must remain unpromoted")

    incidence = face_receipt["incidence"]
    vertices = int(incidence["vertices"])
    edges = int(incidence["edges"])
    faces = int(incidence["faces"])
    euler = int(incidence["euler_characteristic"])
    vertex_degree = int(incidence["vertex_degrees"][0])
    face_degree = int(incidence["face_sizes"][0])
    n_c = n_g = 3

    structural_dimensions = {
        "kappa_0": edges + faces,
        "kappa_1": edges + 1,
        "kappa_2": vertices - n_g + 1,
        "chi_0": 2 ** (vertices - n_g),
        "chi_1": (vertex_degree + euler) * (vertices - 1),
        "zeta_0": faces + 1,
        "zeta_1": edges - n_g,
        "zeta_2": vertex_degree,
    }
    registers = packet_registers()
    register_dimensions = {
        name: register.dimension for name, register in registers.items()
    }
    paths = packet_paths()
    weights = all_path_weights(paths, registers)
    expected_weights = {
        "kappa_write": Fraction(1, 50),
        "kappa_return_1": Fraction(-1, 31),
        "kappa_return_2": Fraction(-1, 310),
        "chi_write": Fraction(1, 512),
        "chi_return_1": Fraction(1, 77),
        "zeta_write": Fraction(1, 21),
        "zeta_return_1": Fraction(1, 27),
        "zeta_return_2": Fraction(1, 135),
    }

    _, graph_certificate = build_register_graphs(
        incidence["face_vertex_indices"], register_dimensions
    )
    graph_checks = graph_certificate["graphs"]

    declared_inputs = declared_map_receipt["declared_hybrid_inputs"]
    alpha_u = mp.mpf(declared_inputs["alpha_U"])
    p_value = mp.mpf(declared_inputs["P"])
    alpha_p = (p_value - (1 + mp.sqrt(5)) / 2) / mp.sqrt(mp.pi)
    calculated_coefficients = {
        "s_kappa": alpha_u / 50,
        "q_kappa": alpha_u / 31 + alpha_u**2 / 310,
        "s_chi": alpha_u**2 / 512,
        "q_chi": alpha_u / 77,
        "s_zeta": alpha_p**2 / 21,
        "q_zeta": alpha_u / 27 + alpha_u**2 / 135,
    }
    declared_coefficients = {
        **{
            name: mp.mpf(value)
            for name, value in declared_map_receipt["repair_map"]["sources"].items()
        },
        **{
            name: mp.mpf(value)
            for name, value in declared_map_receipt["repair_map"]["feedback"].items()
        },
    }
    coefficient_residuals = {
        name: calculated_coefficients[name] - declared_coefficients[name]
        for name in calculated_coefficients
    }

    physical_gates = {
        "CFQ1_three_block_response_exhaustion": False,
        "CFQ2_quotient_homogeneity_and_invariant_state": False,
        "CFQ3_unit_primitive_checkpoint_clock": False,
        "CFQ4_complete_rank_one_path_table": False,
        "CFQ5_exact_sign_character": False,
        "CFQ6_primitive_depth_two_locality": False,
        "CFQ7_path_exhaustion_and_multiplicity_one": False,
        "CFQ8_refinement_implementation_and_ancestry_naturality": False,
    }
    evidence_hashes = {
        "physical_quotient_schema": None,
        "central_block_projectors": None,
        "normalized_trace_matrices": None,
        "rank_one_projector_checks": None,
        "record_algebra_bridge": None,
        "path_lists": None,
        "path_canonicalizer": None,
        "clock_calibration": None,
        "sign_matrix": None,
        "coarse_graining_intertwiners": None,
        "dependency_dag": None,
    }

    checks = {
        "icosahedral_counts_are_12_30_20": (vertices, edges, faces) == (12, 30, 20),
        "euler_vertex_and_face_degrees_are_2_5_3": (
            euler,
            vertex_degree,
            face_degree,
        )
        == (2, 5, 3),
        "declared_register_dimensions_match_structural_arithmetic": (
            register_dimensions == structural_dimensions
        ),
        "all_declared_graph_models_are_connected": all(
            entry["connected"] for entry in graph_checks.values()
        ),
        "all_declared_graph_models_match_dimensions": all(
            entry["dimension_matches"] for entry in graph_checks.values()
        ),
        "eight_exact_path_weights_match": weights == expected_weights,
        "all_path_multiplicities_are_one": all(path.multiplicity == 1 for path in paths),
        "declared_affine_map_coefficients_reproduced_to_emitted_precision": max(
            abs(value) for value in coefficient_residuals.values()
        )
        < mp.mpf("1e-75"),
        "unit_multiplier_is_only_multiplier_admissible_inside_fixed_cfq_packet": (
            multiplier_admissible_inside_cfq(Fraction(1, 1))
            and not multiplier_admissible_inside_cfq(Fraction(1, 2))
            and not multiplier_admissible_inside_cfq(Fraction(2, 1))
        ),
        "all_physical_cfq_gates_remain_unpassed": not any(physical_gates.values()),
        "all_physical_evidence_hashes_remain_absent": all(
            value is None for value in evidence_hashes.values()
        ),
        "historical_target_ancestry_retained": declared_map_receipt[
            "historical_charged_target_informed"
        ],
        "hybrid_branch_warning_retained": not declared_map_receipt[
            "branch_tuple_coherent"
        ],
        "public_promotion_remains_forbidden": not declared_map_receipt[
            "public_prediction_promotion_allowed"
        ],
    }

    return {
        "artifact": "oph_charged_source_law_rigidity_conditional",
        "schema_version": 1,
        "status": (
            "CONDITIONAL_CFQ_PACKET_WEIGHT_RIGIDITY_CLOSED_"
            "PHYSICAL_CFQ_REALIZATION_SOURCE_SELECTION_AND_ATTACHMENTS_OPEN"
        ),
        "runtime_charged_reference_consumed": False,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "mass_scheme_certified": False,
        "public_prediction_promotion_allowed": False,
        "runtime_dependency": {"mpmath": mp.__version__, "external_graph_library": None},
        "provenance": {
            "geometric_face_receipt": str(FACE_RECEIPT.relative_to(CODE_ROOT)),
            "geometric_face_receipt_sha256": face_receipt_sha256,
            "declared_map_receipt": str(DECLARED_MAP_RECEIPT.relative_to(CODE_ROOT)),
            "declared_map_receipt_sha256": declared_map_receipt_sha256,
            "submitted_rigidity_archive_sha256": RIGIDITY_ARCHIVE_SHA256,
            "submitted_narrative_attachment_sha256": RIGIDITY_NARRATIVE_SHA256,
            "submitted_receipt_sha256": RIGIDITY_SUBMITTED_RECEIPT_SHA256,
        },
        "theorem_scope": {
            "proved_implication": (
                "If CFQ1-CFQ8 physically hold exactly, normalized rank-one trace "
                "weights and the frozen path catalogue uniquely reproduce the declared "
                "three-coordinate scalar affine map."
            ),
            "rigidity_scope": (
                "Uniqueness holds inside the stipulated CFQ admissible class. It is not "
                "a proof that current OPH axioms select that class or the physical repair operator."
            ),
            "former_multiplier_witness_status": (
                "The source-multiplier witnesses remain countermodels to the weaker current "
                "OPH premises. They are excluded only after CFQ fixes rank, trace, clock, "
                "path catalogue, signs, exhaustion, and multiplicity one."
            ),
            "definitional_content": (
                "CFQ2-CFQ7 stipulate the coefficient-bearing dimensions, event ranks, "
                "clock, paths, monomials, signs, primitive depth, exhaustion, and multiplicities."
            ),
        },
        "structural_data": {
            "V": vertices,
            "E": edges,
            "F": faces,
            "euler_characteristic": euler,
            "vertex_valence": vertex_degree,
            "face_size": face_degree,
            "N_c": n_c,
            "N_g": n_g,
        },
        "declared_finite_graph_models": {
            **graph_certificate,
            "physical_evidence": False,
            "qualification": (
                "The models verify finite counts and connectivity only. They do not emit "
                "physical response modules, matrix units, recovery maps, clocks, path "
                "exhaustion, refinement intertwiners, or ancestry evidence."
            ),
        },
        "conditional_cfq_packet": {
            "registers": {
                name: {
                    "dimension": register.dimension,
                    "event_rank": register.event_rank,
                    "normalized_trace_weight": fraction_text(register.trace_weight),
                }
                for name, register in registers.items()
            },
            "paths": [
                {
                    "name": path.name,
                    "block": path.block,
                    "monomial": path.monomial,
                    "sign": path.sign,
                    "registers": list(path.registers),
                    "multiplicity": path.multiplicity,
                    "exact_trace_coefficient": fraction_text(weights[path.name]),
                }
                for path in paths
            ],
            "physical_gates": physical_gates,
            "evidence_hashes": evidence_hashes,
        },
        "declared_map_agreement": {
            "formula": "T(k,c,z)=(s_k-q_k*k, s_c+q_c*c, s_z+q_z*z)",
            "coefficient_values": {
                name: mp_text(value) for name, value in calculated_coefficients.items()
            },
            "coefficient_residuals_against_prior_receipt": {
                name: mp_text(value) for name, value in coefficient_residuals.items()
            },
            "inputs_are_historical_hybrid_tuple": True,
            "downstream_mass_values_recomputed_here": False,
        },
        "record_and_maxent_boundary": {
            "central_record_bridge_open": True,
            "statement": (
                "A rank-one projector in B(C^d) is noncentral for d>1, whereas OPH "
                "completed records are normally central. A process/readout time-slice or "
                "diagonal record-algebra theorem must connect the CFQ event projectors to "
                "the OPH record notion."
            ),
            "maxent_boundary": (
                "Current OPH uses constrained Gibbs MaxEnt. The maximally mixed state I/d "
                "follows here only from CFQ2's additional quotient-homogeneous, unconstrained "
                "finite-register premise."
            ),
            "projector_scaling_boundary": (
                "The elementary cP idempotence lemma forbids rescaling one fixed exact "
                "projector. It does not by itself exclude a different state, event rank, "
                "path multiplicity, clock, action, or source construction."
            ),
            "continuous_exact_record_counterexample": (
                "For an exact rank-one P in d dimensions, rho_t=tP+(1-t)(I-P)/(d-1) "
                "is a valid state with Tr(rho_t P)=t. Exact projectivity alone therefore "
                "does not select the weight 1/d."
            ),
        },
        "separate_conditional_corollaries": {
            "one_bit_balance": (
                "Requires an additional physical two-state orientation register plus a "
                "probability-to-Hilbert-Schmidt-power response isometry."
            ),
            "base_phase_2_over_9": (
                "Requires an additional eighteen-slot rank-four response register, action-angle "
                "normalization, and physical charged connection."
            ),
            "determinant_factor_6_to_minus_14": (
                "Requires fourteen independent physical rank-one Z6 records and their "
                "attachment to the charged determinant line."
            ),
            "closed_by_core_cfq_theorem": False,
        },
        "open": [
            "derive CFQ1-CFQ8 from OPH dynamics with frozen no-target ancestry",
            "bridge central completed OPH records to the CFQ rank-one event projectors",
            "derive the register dimensions and eight-path catalogue rather than stipulating them",
            "prove physical path exhaustion, signs, multiplicity one, clock, and refinement naturality",
            "attach the A5/C3 face-corner bundle equivariantly to physical charged-family lines",
            "derive the bare-to-endpoint block-balance response bridge and charged phase connection",
            "attach the primitive determinant character to the physical charged determinant line",
            "freeze one coherent source branch before comparison",
            "derive a declared running coordinate and controlled QED/electroweak pole-scheme map",
        ],
        "checks": checks,
        "checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--face-receipt", type=Path, default=FACE_RECEIPT)
    parser.add_argument("--declared-map-receipt", type=Path, default=DECLARED_MAP_RECEIPT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    face_raw = args.face_receipt.read_bytes()
    declared_raw = args.declared_map_receipt.read_bytes()
    artifact = build_artifact(
        json.loads(face_raw),
        sha256(face_raw),
        json.loads(declared_raw),
        sha256(declared_raw),
    )
    encoded = (json.dumps(artifact, indent=2, sort_keys=True) + "\n").encode()

    if args.check:
        actual = args.out.read_bytes() if args.out.exists() else None
        ok = actual == encoded
        print(json.dumps({"status": "OK" if ok else "DRIFT"}, indent=2))
        return 0 if ok else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "physical_cfq_gates_passed": sum(
                    artifact["conditional_cfq_packet"]["physical_gates"].values()
                ),
                "public_prediction_promotion_allowed": artifact[
                    "public_prediction_promotion_allowed"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
