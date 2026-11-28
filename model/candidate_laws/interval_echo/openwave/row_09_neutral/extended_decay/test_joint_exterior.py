"""Independent checks of the actual one-resident stock and native stationary edge."""

import json
from hashlib import sha256
import numpy as np
from scipy.linalg import eigh
from model.candidate_laws.interval_echo import candidate
from . import joint_exterior as joint


def records():
    return json.loads((joint.HERE / "joint_exterior.json").read_text())


def test_full_physical_one_resident_matrix_and_complete_spectrum_on_smallest_mount():
    r = records()["boards"]["20"]
    board = joint.exterior.federation(0)
    sites = sorted(
        {i for i in range(board.n) if i // 12 not in r["reserved_owners"]}
        | set(r["path"][1:-1])
    )
    index = {p: i for i, p in enumerate(sites)}
    actual, _ = joint.single(board, sites)
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    words = [law.SRC.index(w) for w in [(1, 0), (0, 1)]]
    pair = h[np.ix_(words, words)]
    oracle = np.diag((candidate.KAPPA * (6 - board.degree[sites])).astype(float))
    for p, q in board.seams:
        if p in index and q in index:
            local = [index[p], index[q]]
            oracle[np.ix_(local, local)] += pair
    np.testing.assert_allclose(actual.toarray(), oracle, atol=2e-13)
    e, v = eigh(oracle)
    rec = r["exterior"]["with_approach"]
    np.testing.assert_allclose(rec["first_gaps"], (e - e[0])[:8], atol=5e-11)
    lo, hi = rec["recoil_interval"]
    assert np.count_nonzero((e - e[0] >= lo) & (e - e[0] <= hi)) == 0
    bitphase = np.ones(len(sites))
    bitphase[index[r["path"][1]]] = -1
    phase = v.T @ (bitphase * v[:, 0])
    first = abs(e - e[1]) < 1e-7
    np.testing.assert_allclose(
        rec["opening"]["phase_first_excited_probability"],
        sum(phase[first] ** 2),
        atol=2e-12,
    )


def test_actual_source_provenance_and_ground_bounds_on_the_held_sizes():
    data = records()
    for key, protocol in (
        ("protocol_sha256", joint.DESIGN),
        ("stationary_protocol_sha256", "STATIC_EDGE.md"),
    ):
        assert data[key] == sha256((joint.HERE / protocol).read_bytes()).hexdigest()
    for level in range(3):
        b = joint.exterior.federation(level)
        rec = data["boards"][str(b.cells)]
        d, heads = rec["dimers"], rec["heads"]
        assert len(d) == 11 and len(heads) == 30
        assert all(v in b.neighbours[u] for u, v in d)
        assert all(v in b.neighbours[u] for u, v, _ in heads)
        assert rec["manufacture"]["total_load"] == rec["manufacture_ports"] - 11
        path = rec["path"]
        assert all(v in b.neighbours[u] for u, v in zip(path, path[1:]))
        assert path[0] in d[-1] and path[-1] // 12 not in rec["reserved_owners"]
        assert (
            abs(rec["source"]["native_lower_core_and_second"] - 0.09589617515795)
            < 5e-12
        )
        assert rec["source"]["source_cut_lower_core_and_second"] < 0.002
        assert rec["source"]["native_idle_lower_core_and_second"] < 2e-10
        bare_counts = rec["bare_interval_inertia"]
        matches = rec["bare_gaps_in_interval"]
        assert len(matches) == bare_counts[1] - bare_counts[0]
        assert len(matches) == (2 if b.cells == 320 else 0)
        for phase in rec["exterior"].values():
            left, right = phase["target_inertia"]
            assert left == right
            assert abs(phase["contrast_removed"]["ground"]) < 1e-11
        approach = rec["exterior"]["with_approach"]
        cert = approach["approach_certificate"]
        assert cert["first_gap_lower"] > rec["source"]["budget"]
        assert cert["first_gap_lower"] <= approach["first_gaps"][1]
        assert (
            1 - approach["ground_corridor_probability"]
            <= cert["ground_outside_probability_upper"]
        )
        assert approach["opening"]["idle_first_excited_probability"] == 0
        if b.cells > 20:
            np.testing.assert_array_equal(
                cert["approach_block"], [[16.5, -0.5], [-0.5, 17]]
            )
            assert approach["opening"]["phase_first_excited_probability"] > 0.75
            assert cert["first_gap_lower"] > 1.10


def test_all_stationary_native_edge_energies_and_sharp_pins():
    full = records()["stationary_boundary"]
    assert full == joint.classify()
    for pin in range(3):
        left, right = (full[f"{side}_resident_pin_{pin}"] for side in ("left", "right"))
        assert left == right
        shifts = {e1 - e0 for e0, e1, _ in left}
        assert shifts == ({0, -32, -33} if pin == 1 else {0})
        assert sum(d for _, _, d in left) == 8
    # The open, blank edge does supply +8.5, with nonzero physical exchange.
    h = candidate.native_law.seam_term(candidate.KAPPA)
    words = candidate.native_law.SRC
    p, q = words.index((1, 0)), words.index((0, 1))
    assert h[p, p] == 8.5 and h[p, q] == -0.5


def test_larger_heavy_remnant_exhausts_downward_strand_transitions():
    data = json.loads((joint.HERE / "larger_transition.json").read_text())
    assert (
        data["joint_archive_sha256"]
        == sha256((joint.HERE / "joint_exterior.json").read_bytes()).hexdigest()
    )
    assert (
        data["protocol_sha256"]
        == sha256((joint.HERE / "LARGER_TRANSITION.md").read_bytes()).hexdigest()
    )
    levels = data["strand_levels"]
    assert len(levels) == 15 and levels[-1] < data["larger_parent_loss"]
    assert data["boards"]["20"]["exterior_gaps"] == []
    for n in ("80", "320"):
        row = data["boards"][n]
        residuals = [
            levels[i] - levels[k] - row["exterior_gaps"][0]
            for i in range(1, len(levels))
            for k in range(i)
        ]
        assert len(row["tried"]) == len(residuals) == 105
        np.testing.assert_allclose([item[-1] for item in row["tried"]], residuals)
        assert (
            row["closest"][:3] == [4, 0, 0] and min(abs(x) for x in residuals) > 0.0503
        )
    rates = data["source_lower_remnant_and_level"]
    assert rates["native_pulse"][4] > 0.011
    assert rates["source_cut_pulse"][4] < 0.000043
    assert rates["native_idle"][4] < 1e-10
    measured = json.loads((joint.HERE / "strand_lines.json").read_text())["compared"][
        "8"
    ][3]
    assert abs(measured["gap"] - levels[4]) < measured["radius"]
