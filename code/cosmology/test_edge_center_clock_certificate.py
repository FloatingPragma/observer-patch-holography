from __future__ import annotations

import itertools
import json
import math

from mpmath import iv, mp, mpf
import pytest

from edge_center_clock_certificate import (
    EXACT_DEFECT_WIDTH_BOUND,
    P_CERTIFICATE_PATH,
    SELECTED_PRIMARY_BRANCH,
    CertificateError,
    build,
    build_check_primary_branch,
    collar_tower_records,
    load_certified_p_interval,
    load_repair_round_invariant,
    reject_forbidden_inputs,
    survival_family_records,
)


@pytest.fixture(scope="module")
def payload() -> dict:
    return build()


# ---------------------------------------------------------------------------
# independent brute-force icosahedral edge-midpoint refinement
# ---------------------------------------------------------------------------


def _icosahedron_faces() -> list[tuple[int, int, int]]:
    golden = (1.0 + math.sqrt(5.0)) / 2.0
    coordinates = []
    for a, b in itertools.product((-1.0, 1.0), repeat=2):
        coordinates.append((0.0, a, b * golden))
        coordinates.append((a, b * golden, 0.0))
        coordinates.append((a * golden, 0.0, b))
    def squared_distance(u, v):
        return sum((x - y) ** 2 for x, y in zip(u, v))

    edges = {
        (i, j)
        for i in range(12)
        for j in range(i + 1, 12)
        if abs(squared_distance(coordinates[i], coordinates[j]) - 4.0) < 1e-9
    }
    faces = [
        (i, j, k)
        for i in range(12)
        for j in range(i + 1, 12)
        for k in range(j + 1, 12)
        if (i, j) in edges and (i, k) in edges and (j, k) in edges
    ]
    return faces


def _refine_and_count(depth: int) -> dict[str, int]:
    """Brute-force edge-midpoint refinement with skeleton tracking."""
    faces = _icosahedron_faces()
    skeleton = {}
    on_skeleton: dict[int, frozenset[int]] = {}
    next_vertex = 12
    for face in faces:
        for pair in itertools.combinations(face, 2):
            key = frozenset(pair)
            if key not in skeleton:
                skeleton[key] = len(skeleton)
    for vertex in range(12):
        on_skeleton[vertex] = frozenset(
            edge_id for key, edge_id in skeleton.items() if vertex in key
        )

    midpoint_memo: dict[frozenset[int], int] = {}

    def midpoint(u: int, w: int) -> int:
        nonlocal next_vertex
        key = frozenset((u, w))
        if key not in midpoint_memo:
            midpoint_memo[key] = next_vertex
            on_skeleton[next_vertex] = on_skeleton[u] & on_skeleton[w]
            next_vertex += 1
        return midpoint_memo[key]

    current = faces
    for _ in range(depth):
        refined = []
        for a, b, c in current:
            mab, mbc, mca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(
                [(a, mab, mca), (b, mbc, mab), (c, mca, mbc), (mab, mbc, mca)]
            )
        current = refined
        midpoint_memo.clear()

    mesh_edges = set()
    for face in current:
        for pair in itertools.combinations(face, 2):
            mesh_edges.add(frozenset(pair))

    def on_skeleton_edge(edge: frozenset[int]) -> bool:
        u, w = tuple(edge)
        return bool(on_skeleton[u] & on_skeleton[w])

    skeleton_sub_edges = sum(1 for edge in mesh_edges if on_skeleton_edge(edge))
    collar_faces = 0
    for face in current:
        if any(
            on_skeleton_edge(frozenset(pair))
            for pair in itertools.combinations(face, 2)
        ):
            collar_faces += 1
    vertices = {vertex for face in current for vertex in face}
    return {
        "faces": len(current),
        "edges": len(mesh_edges),
        "vertices": len(vertices),
        "skeleton_sub_edges": skeleton_sub_edges,
        "collar_cells": collar_faces,
    }


@pytest.mark.parametrize("depth", [0, 1, 2])
def test_tower_records_match_brute_force_mesh(depth: int) -> None:
    brute = _refine_and_count(depth)
    record = collar_tower_records(max(depth, 1))[depth]
    assert record["faces"] == brute["faces"]
    assert record["edges"] == brute["edges"]
    assert record["vertices"] == brute["vertices"]
    assert record["skeleton_sub_edges"] == brute["skeleton_sub_edges"]
    assert record["collar_cells"] == brute["collar_cells"]


def test_tower_exact_invariants() -> None:
    records = collar_tower_records(12)
    for record in records:
        assert record["euler_characteristic"] == 2
        assert (
            record["collar_cells_route_a_per_face"]
            == record["collar_cells_route_b_per_sub_edge"]
        )
        assert record["collar_cells"] + record["full_cells"] == record["faces"]
        assert record["orientation_balance_defect"] == 0
        involution = record["orientation_reversal_involution"]
        assert involution["fixed_point_free"] is True
        assert involution["half_over_full_slot_ratio"] == "1/2"
        assert 2 * involution["half_collar_slot_count"] == record["oriented_collar_slots"]
    # Exact closed forms at a few depths.
    assert records[1]["collar_cells"] == 60
    assert records[1]["full_cells"] == 20
    assert records[2]["collar_cells"] == 180
    assert records[2]["full_cells"] == 140


# ---------------------------------------------------------------------------
# generator, half-collar identity, and defect gates
# ---------------------------------------------------------------------------


def test_p_interval_is_loaded_from_certified_artifact() -> None:
    record = load_certified_p_interval()
    raw = json.loads(P_CERTIFICATE_PATH.read_text(encoding="utf-8"))
    enclosure = raw["modes"]["thomson_structured_running"]["certified_enclosure"]["P"]
    assert record["lo"] == enclosure["lo"]
    assert record["hi"] == enclosure["hi"]
    assert record["exact_alpha_promoted"] is False


def test_repair_round_invariant_is_loaded_and_checked(payload: dict) -> None:
    record = load_repair_round_invariant()
    assert record["m_rep"] == 24
    assert record["exponent_denominator"] == 48
    antecedent = payload["antecedents"]["repair_round_invariant"]
    assert antecedent["m_rep"] == 24
    assert antecedent["specialized_exponent"] == "-1/48"


def test_generator_is_p_over_24_interval(payload: dict) -> None:
    mp.dps = 60
    record = load_certified_p_interval()
    lo = mpf(record["lo"]) / 24
    hi = mpf(record["hi"]) / 24
    emitted = payload["generator"]["full_collar_derivative"]
    assert mpf(emitted["lo"]) <= lo <= hi <= mpf(emitted["hi"]) * (1 + mpf("1e-38"))
    assert abs(mpf(emitted["lo"]) - lo) < mpf("1e-38")


def test_half_collar_identity_intervals(payload: dict) -> None:
    mp.dps = 60
    half = payload["half_collar_identity"]
    theta_lo, theta_hi = mpf(half["theta"]["lo"]), mpf(half["theta"]["hi"])
    generator_lo = mpf(payload["generator"]["full_collar_derivative"]["lo"])
    generator_hi = mpf(payload["generator"]["full_collar_derivative"]["hi"])
    # theta = (1/2) * generator within display rounding.
    assert abs(theta_lo - generator_lo / 2) < mpf("1e-38")
    assert abs(theta_hi - generator_hi / 2) < mpf("1e-38")
    # n_s = 1 - theta within display rounding.
    assert abs(mpf(half["n_s"]["lo"]) - (1 - theta_hi)) < mpf("1e-38")
    assert abs(mpf(half["n_s"]["hi"]) - (1 - theta_lo)) < mpf("1e-38")
    assert half["half_over_full_slot_ratio"] == "1/2"


def test_kappa_rep_interval_and_issue_reconciliation(payload: dict) -> None:
    mp.dps = 60
    half = payload["half_collar_identity"]
    kappa_lo = mpf(half["kappa_rep"]["lo"])
    kappa_hi = mpf(half["kappa_rep"]["hi"])
    assert kappa_lo < kappa_hi
    assert mpf("2.6262") < kappa_lo and kappa_hi < mpf("2.6263")
    reconciliation = half["kappa_rep_issue_body_reconciliation"]
    assert reconciliation["inside_certified_interval"] is False
    assert reconciliation["issue_body_value"] == "2.627023712627471"


def test_survival_family_defect_gates(payload: dict) -> None:
    mp.dps = 60
    records = payload["survival_family"]["records"]
    assert len(records) == 13
    previous = None
    for record in records:
        defects = record["defects"]
        for key in ("semigroup_on_grid", "sub_step_derivative"):
            value = defects[key]["value"]
            assert mpf(value["lo"]) <= 0 <= mpf(value["hi"])
            assert mpf(value["hi"]) - mpf(value["lo"]) < mpf(EXACT_DEFECT_WIDTH_BOUND)
        for key in (
            "limit_family_derivative",
            "orientation_balance",
            "refinement_to_poisson_floor",
        ):
            value = defects[key]["value"]
            bound = mpf(defects[key]["bound"])
            assert mpf(value["hi"]) <= bound
        assert mpf(defects["refinement_to_poisson_floor"]["value"]["lo"]) > 0
        one_tick = mpf(record["one_tick_survival"]["lo"])
        if previous is not None:
            assert one_tick > previous
        previous = mpf(record["one_tick_survival"]["hi"])


def test_refinement_defect_shrinks_with_depth(payload: dict) -> None:
    mp.dps = 60
    records = payload["survival_family"]["records"]
    highs = [
        mpf(record["defects"]["refinement_to_poisson_floor"]["value"]["hi"])
        for record in records
    ]
    assert all(later < earlier for earlier, later in zip(highs, highs[1:]))


def test_presence_one_step_matches_depth_zero_family(payload: dict) -> None:
    mp.dps = 60
    presence = payload["generator"]["finite_one_step_presence_value"]
    depth_zero = payload["survival_family"]["records"][0]["one_tick_survival"]
    assert abs(mpf(presence["lo"]) - mpf(depth_zero["lo"])) < mpf("1e-38")
    assert abs(mpf(presence["hi"]) - mpf(depth_zero["hi"])) < mpf("1e-38")


def test_survival_family_direct_construction() -> None:
    iv.dps = 60
    mp.dps = 60
    record = load_certified_p_interval()
    P = iv.mpf([record["lo"], record["hi"]])
    family = survival_family_records(P, 4)
    assert [entry["sub_slots_per_tick"] for entry in family] == [1, 2, 4, 8, 16]


# ---------------------------------------------------------------------------
# diagnostics and rejected branches
# ---------------------------------------------------------------------------


def test_finite_transition_exponent_is_diagnostic_and_disjoint(payload: dict) -> None:
    mp.dps = 60
    diagnostic = payload["diagnostics"]["finite_transition_exponent"]
    assert diagnostic["diagnostic_only"] is True
    theta_hi = mpf(payload["half_collar_identity"]["theta"]["hi"])
    assert mpf(diagnostic["value"]["lo"]) > theta_hi


def test_e_branch_is_diagnostic_only(payload: dict) -> None:
    diagnostic = payload["diagnostics"]["e_branch"]
    assert diagnostic["diagnostic_only"] is True
    assert payload["selected_primary_branch"] == SELECTED_PRIMARY_BRANCH


def test_wrong_orientation_branch_is_flagged_rejected(payload: dict) -> None:
    mp.dps = 60
    branch = payload["rejected_branches"]["wrong_orientation_factor_1"]
    assert branch["status"] == "rejected_orientation_branch"
    generator = payload["generator"]["full_collar_derivative"]
    assert branch["theta"] == generator
    theta_hi = mpf(payload["half_collar_identity"]["theta"]["hi"])
    assert mpf(branch["theta"]["lo"]) > theta_hi


# ---------------------------------------------------------------------------
# controls: fail-closed rejections
# ---------------------------------------------------------------------------


def test_injected_declared_step_time_is_rejected() -> None:
    with pytest.raises(CertificateError, match="forbidden input key"):
        build(injected_inputs={"step_time_seconds": "5.39e-44"})


def test_injected_measured_tilt_target_is_rejected() -> None:
    with pytest.raises(CertificateError, match="measured-tilt-shaped"):
        build(injected_inputs={"primordial_amplitude_hint": 0.9649})


def test_injected_named_target_key_is_rejected() -> None:
    with pytest.raises(CertificateError, match="forbidden input key"):
        build(injected_inputs={"n_s_target": 0.965})


def test_e_branch_as_primary_fails_closed() -> None:
    with pytest.raises(CertificateError, match="diagnostics only"):
        build(primary_branch="repair_clock_e")
    with pytest.raises(CertificateError, match="diagnostics only"):
        build_check_primary_branch("repair_clock_e")


def test_forbidden_input_scan_accepts_clean_inputs() -> None:
    reject_forbidden_inputs(
        {"structure_integers": {"m_rep_round_count": 24}, "P_lo": "1.63"}
    )


def test_controls_are_recorded(payload: dict) -> None:
    controls = payload["controls"]
    for name in (
        "injected_declared_step_time",
        "injected_measured_tilt_target",
        "e_branch_selected_as_primary",
        "wrong_orientation_factor_1",
    ):
        assert controls[name]["rejected"] is True


# ---------------------------------------------------------------------------
# clock binding and freeze policy
# ---------------------------------------------------------------------------


def test_clock_binding_has_no_step_time(payload: dict) -> None:
    binding = payload["clock_binding"]
    assert binding["declared_command_line_step_time"] is None
    assert binding["clock_unit"] == "one accepted oriented repair event on the collar"
    source = binding["clock_normalization_source"]
    assert "m_rep = 24" in source["invariant"]
    assert source["artifacts"]
    hygiene = payload["input_hygiene"]
    assert hygiene["declared_step_time_present"] is False
    assert hygiene["time_dimension_constants_present"] is False
    assert hygiene["measured_alpha_present"] is False


def test_no_sky_comparison_and_not_frozen(payload: dict) -> None:
    assert payload["sky_comparison"]["performed"] is False
    assert payload["freeze_status"] == "not_frozen_here"
    text = json.dumps(payload).lower()
    for token in ("planck 2018", "likelihood_value", "sky_value", "cobe", "wmap"):
        assert token not in text


def test_acceptance_mapping_is_honest(payload: dict) -> None:
    mapping = payload["acceptance_mapping"]
    assert mapping["clock_normalization_source_identified"]["discharged_here"] is True
    assert (
        mapping["full_collar_derivative_and_half_collar_identity"]["discharged_here"]
        is True
    )
    assert (
        mapping["simulator_emits_p_over_48_only_from_receipt"]["discharged_here"]
        is False
    )


def test_payload_is_deterministic(payload: dict) -> None:
    again = build()
    assert json.dumps(payload, sort_keys=True) == json.dumps(again, sort_keys=True)


def test_schema_and_issue_fields(payload: dict) -> None:
    assert payload["schema"] == "oph.edge_center_clock_certificate.v1"
    assert payload["artifact"] == "oph_edge_center_clock_certificate"
    assert payload["github_issue"] == 522
    assert payload["status"] == "edge_center_generator_and_clock_certificate_emitted"
