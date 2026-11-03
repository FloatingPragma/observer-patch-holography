"""Hostile checks of the live-source relative attachment audit."""

import copy
import fnmatch
import json
from pathlib import Path
import re
import sys

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import port_support_source_audit as producer
import verify_port_support_source_audit as independent


@pytest.fixture(scope="module")
def packet():
    return json.loads(producer.RECEIPT.read_text())


def test_independent_live_source_reconstruction(packet):
    assert independent.verify(packet)["raw_maps"] == 60


def test_ci_push_and_pr_triggers_cover_the_whole_source_closure():
    workflow = (producer.ROOT / ".github/workflows/arithmon-interface-contracts.yml").read_text()
    sections = [workflow.split("  push:", 1)[1].split("  pull_request:", 1)[0],
                workflow.split("  pull_request:", 1)[1].split("jobs:", 1)[0]]
    for section in sections:
        patterns = re.findall(r"- '([^']+)'", section)
        assert all(any(fnmatch.fnmatch(path, pattern) for pattern in patterns) for path in producer.SOURCES)


def test_every_seed_relabeling_has_typed_action_without_selecting(packet):
    group = [tuple(p) for p in packet["seed"]["raw_maps"]]
    identity = tuple(range(12))
    a = group[1]
    assert a != identity
    for b in group:
        # Local-only: h -> h a^-1. Support-only: h -> b h.
        assert producer.compose(b, producer.inverse(a)) in group
        assert producer.compose(a, b) in group
        assert producer.compose(producer.compose(a, b), producer.inverse(a)) in group
    assert producer.compose(identity, producer.inverse(a)) != identity
    assert producer.compose(a, identity) != identity
    assert producer.compose(producer.compose(a, identity), producer.inverse(a)) == identity
    assert packet["seed"]["fixed_support_left_orbit_sizes"] == [5]*12
    assert packet["seed"]["diagonal_conjugacy_orbit_sizes"] == [1, 12, 12, 15, 20]
    # The diagonal identity is only distinguished after a cross-object chart
    # matching has been supplied. It is not a source invariant raw selector.


def test_oriented_integral_chain_cannot_be_replaced_by_unsigned_faces():
    faces = producer.local_faces()
    one_reversed = [(faces[0][0], faces[0][2], faces[0][1]), *faces[1:]]
    assert {frozenset(f) for f in faces} == {frozenset(f) for f in one_reversed}
    assert producer.chain(faces) != producer.chain(one_reversed)
    all_reversed = [(a, c, b) for a, b, c in faces]
    assert producer.chain(faces) != producer.chain(all_reversed)
    assert producer.automorphisms(faces) == producer.automorphisms(all_reversed)


def test_exact_marked_generator_midpoint_falsifier():
    rows = producer.source_rows()
    p = (0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11)
    fine = producer.lift(p, 12, rows[0]["support"]["faces"])
    c = rows[1]["support"]["coarsen"]
    assert fine[18] == 17 and c[17] == 2 and p[c[18]] == 1
    # Transporting the declared c is covariant by definition; recomputing max
    # on the relabeled labels gives a different convention.
    transported = [None]*len(c)
    for v in range(len(c)):
        transported[fine[v]] = p[c[v]]
    assert transported != c
    assert all(transported[fine[v]] == p[c[v]] for v in range(len(c)))


def test_natural_geometric_midpoint_survives_every_seed_action(packet):
    for row in packet["native_support_refinement"]:
        assert row["subdivision_chain_symmetries"] == 60
        assert row["coarsening_symmetries"] == 1


def test_no_equivariant_endpoint_rounding_in_the_declared_class(packet):
    witnesses = packet["seed"]["endpoint_choice_obstruction"]["edge_swap_witnesses"]
    assert len(witnesses) == 30
    for witness in witnesses:
        a, b = witness["edge"]
        swap = witness["swap"]
        assert swap[a] == b and swap[b] == a
        # The lifted action fixes this edge's midpoint, but fixes neither
        # admissible endpoint. Either possible endpoint choice fails.
        assert all(swap[endpoint] != endpoint for endpoint in (a, b))


def test_cyclic_normalizers_separate_the_two_domain_poles(packet):
    for row in packet["seed"]["cyclic_axis_normalizers"]:
        assert row["normalizer_size"] == {2: 4, 3: 6, 5: 10}[row["prime"]]
        left, right = row["domain_axis_poles"]
        swap = row["pole_swap_involution"]
        assert {swap[v] for v in left} == set(right)
        assert {swap[v] for v in right} == set(left)


def test_corrupt_midpoint_and_skipped_level_fail():
    rows = [r["support"] for r in producer.source_rows()]
    group = producer.automorphisms(producer.local_faces())
    with pytest.raises(ValueError, match="consecutive"):
        producer.refinement_audit(group, [rows[0], rows[2], rows[3]])
    bad = copy.deepcopy(rows)
    bad[1]["coarsen"][18] = 11  # outside the parent edge {1,3}
    with pytest.raises(ValueError):
        producer.refinement_audit(group, bad)


@pytest.mark.parametrize("mutation", [
    "missing_map", "reverse_map", "support_count", "C5", "physical_orbits",
    "coarsening_natural", "midpoint_square", "skip_level", "unsigned_degree",
    "invented_local_tower", "invented_common_sphere", "physical_selection",
    "invented_drift", "source_hash", "implementation_hash", "target_selected",
])
def test_forged_scientific_receipts_fail(packet, mutation):
    bad = copy.deepcopy(packet)
    if mutation == "missing_map":
        bad["seed"]["raw_maps"].pop()
    elif mutation == "reverse_map":
        bad["seed"]["raw_maps"][0][1:3] = [2, 1]
    elif mutation == "support_count":
        bad["support_tower"][0]["nerve_carriers"] = 12
    elif mutation == "C5":
        bad["seed"]["fixed_support_stabilizer"].pop()
    elif mutation == "physical_orbits":
        bad["seed"]["fixed_support_left_orbit_sizes"] = [1]*12
    elif mutation == "coarsening_natural":
        bad["native_support_refinement"][0]["coarsening_symmetries"] = 60
    elif mutation == "midpoint_square":
        bad["native_support_refinement"][0]["first_counterexample"]["relabel_after_coarsen"] = 2
    elif mutation == "skip_level":
        bad["native_support_refinement"].pop(1)
    elif mutation == "unsigned_degree":
        bad["support_tower"][0]["bridge_chain_degree"] = -1
    elif mutation == "invented_local_tower":
        bad["local_refinement"] = "RE FINEMENT NATURAL".replace(" ", "_")
    elif mutation == "invented_common_sphere":
        bad["common_sphere"] = "CERTIFIED"
    elif mutation == "physical_selection":
        bad["physical_selection"] = True
    elif mutation == "invented_drift":
        # Every pin now matches, so the disclosure is empty; inventing a drift
        # row is the lie the verifier has to reject.
        bad["historical_drift"] = {"code/source_selection_model/DERIVATION.md": {
            "classification": "BASELINE_PIN_DRIFT",
            "frozen_sha256": "0" * 64, "live_sha256": "0" * 64}}
    elif mutation == "source_hash":
        bad["sources"][producer.GEOMETRY]["sha256"] = "0"*64
    elif mutation == "implementation_hash":
        bad["implementation_sha256"]["code/a5_closure/port_support_source_audit.py"] = "0"*64
    else:
        bad["verdicts"] = ["G_PLUS_SELECTED"]
    with pytest.raises(ValueError):
        independent.verify(bad)
