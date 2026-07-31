from __future__ import annotations

from fractions import Fraction

import pytest

import angular_interpolant_certificate as ac


NARROWED_STATUS = (
    "EXACT_SOURCE_TEMPLATE__STATIC_BASE_PORT_TRANSFER_NONIDENTIFIABLE__"
    "REFINEMENT_REPAIR_SKY_TRANSFER_OPEN"
)
RETIRED_STATUS = "EXACT_SOURCE_TEMPLATE__TRANSFER_NONIDENTIFIABLE"


def test_receipt_builds_with_expected_status() -> None:
    receipt = ac.build_receipt()
    assert receipt["status"] == NARROWED_STATUS
    assert receipt["projector_certificate"]["mutually_orthogonal"] is True
    assert receipt["projector_certificate"]["resolution_of_identity"] is True
    assert receipt["projector_certificate"]["ranks"] == ["1", "3", "5", "3"]
    assert receipt["equivariance_certificate"]["all_equivariant"] is True
    assert receipt["band_binding_certificate"]["all_bound"] is True
    assert receipt["parity_response_certificate"]["all_exact"] is True
    assert receipt["comparison_boundary"]["comparison_permitted"] is False


def test_equal_port_vector_matches_plan_values() -> None:
    receipt = ac.build_receipt()
    vector = receipt["equal_port_certificate"]["even_initial_vector"]
    assert vector == {
        "2": "0",
        "4": "0",
        "6": "11/25",
        "8": "0",
        "10": "247/1875",
        "12": "1071/3125",
        "14": "0",
    }
    assert receipt["equal_port_certificate"]["all_odd_levels_zero"] is True


def test_transfer_decision_disagreement_is_exact() -> None:
    receipt = ac.build_receipt()
    decision = receipt["transfer_decision"]
    assert decision["verdict"] == "STATIC_BASE_PORT_UNDERDETERMINATION"
    assert decision["normalized_statistics_disagree_exactly"] is True
    assert decision["static_base_port_underdetermination"] is True
    assert decision["completion_a"]["normalized_power_at_level_6"] == "0"
    assert decision["completion_b"]["normalized_power_at_level_6"] == "11/25"
    assert decision["completion_a"]["object_type"] == (
        "smooth band-limited field of degree at most 3"
    )
    assert decision["completion_b"]["object_type"] == (
        "discrete equal-port carrier measure"
    )
    search = decision["geometry_imprint_search"]
    assert search["sky_field_emission_found"] is False
    assert search["second_independent_channel_found"] is False
    assert len(decision["distinct_source_objects"]) == 3


def test_stop_rule_keeps_transfer_row_open() -> None:
    receipt = ac.build_receipt()
    decision = receipt["transfer_decision"]
    assert "does not close" in decision["stop_rule"]
    assert "stays open" in decision["stop_rule"]
    assert receipt["row_state"]["refinement_repair_sky_transfer"] == "OPEN"
    assert receipt["row_state"]["static_base_port_underdetermination"] == (
        "CERTIFIED"
    )
    assert "intertwine refinement" in receipt["row_state"]["close_conditions"]


def test_refinement_frontier_is_recorded() -> None:
    receipt = ac.build_receipt()
    frontier = receipt["transfer_decision"]["refinement_frontier"]
    assert frontier["producer"] == (
        "code/angular_sprint/refinement_transfer_certificate.py"
    )
    assert frontier["repair_law"] == (
        "oph-physics-sim docs/CANONICAL_REPAIR_LAW.md"
    )
    assert "level-one refinement vertex set" in frontier["separation_witness"]
    assert "vanishing on the twelve base ports" in (
        frontier["separation_witness"]
    )
    assert "intertwine" in frontier["future_verdict_requirement"]
    assert "repair semigroup" in frontier["future_verdict_requirement"]


def test_narrowed_verdict_string_is_pinned() -> None:
    assert ac.STATUS == NARROWED_STATUS
    receipt = ac.build_receipt()
    assert receipt["status"] == NARROWED_STATUS


def test_receipt_drops_general_nonidentifiability_status() -> None:
    committed = ac.RECEIPT_PATH.read_bytes()
    assert RETIRED_STATUS.encode("ascii") not in committed
    assert NARROWED_STATUS.encode("ascii") in committed
    assert b"reopen_condition" not in committed


def test_committed_receipt_is_byte_exact() -> None:
    committed = ac.RECEIPT_PATH.read_bytes()
    assert committed == ac.canonical_json_bytes(ac.build_receipt())


def test_projectors_match_independent_lagrange_construction() -> None:
    """Independent check: build the band projectors from the adjacency
    matrix by exact Lagrange interpolation over its four eigenvalues and
    compare with the scaled Legendre Gram kernels."""

    data = ac.build_kernels()
    adjacency = ac.int_to_q5(data["components"]["adjacency"])
    identity = ac.int_to_q5(data["components"]["identity"])
    eigenvalues = [ac.q5(5), ac.q5(0, 1), ac.q5(-1), ac.q5(0, -1)]
    scales = [Fraction(12), Fraction(4), Fraction(12, 5), Fraction(4)]

    def q5_div(x, y):
        norm = y[0] * y[0] - 5 * y[1] * y[1]
        conj = (y[0], -y[1])
        num = ac.q5_mul(x, conj)
        return (num[0] / norm, num[1] / norm)

    for target in range(4):
        projector = identity
        for other in range(4):
            if other == target:
                continue
            shift = [
                [
                    ac.q5_sub(
                        adjacency[i][j],
                        ac.q5_mul(eigenvalues[other], identity[i][j]),
                    )
                    for j in range(12)
                ]
                for i in range(12)
            ]
            projector = ac.mat_mul(projector, shift)
            denominator = ac.q5_sub(eigenvalues[target], eigenvalues[other])
            projector = [
                [q5_div(value, denominator) for value in row]
                for row in projector
            ]
        expected = ac.mat_scale(data["kernels"][target], Fraction(1) / scales[target])
        assert ac.mat_is_zero(ac.mat_sub(projector, expected)), target


def test_legendre_values_against_closed_forms() -> None:
    plus = ac.legendre_at(ac.INV_SQRT5, 6)
    assert plus[2] == ac.q5(Fraction(-1, 5))
    assert plus[4] == ac.q5(Fraction(-1, 5))
    assert plus[6] == ac.q5(Fraction(41, 125))
    assert plus[1] == ac.INV_SQRT5
    assert plus[3][0] == 0


def test_kernel_tamper_fails_closed() -> None:
    data = ac.build_kernels()
    kernels = data["kernels"]
    kernels[1][0][1] = ac.q5_add(kernels[1][0][1], ac.q5(1))
    certificate = ac.projector_certificate(data)
    assert certificate["scaled_projectors"] is False


def test_antipode_construction_is_validated() -> None:
    ports, adjacency, antipode, second = ac.ports_and_structure()
    assert len(ports) == 12
    assert sum(sum(row) for row in adjacency) == 60
    assert sum(sum(row) for row in second) == 60
    assert sum(sum(row) for row in antipode) == 12
