"""Semantic attacks on primal/dual certificates, native words and attempt costs."""
from copy import deepcopy
from fractions import Fraction as F
import os
from pathlib import Path
import re
import subprocess
import sys
from itertools import combinations, product

import pytest

from . import build, codec, verify, tomography, check_tomography


@pytest.fixture(scope="module")
def packet():
    return codec.load(codec.HERE/"controls.json")


def test_full_independent_replay(packet):
    codec.equal(verify.verify(packet), codec.load(codec.HERE/"receipt.json"), "receipt")


@pytest.mark.parametrize("rows,target", [
    ([], [1, 0, 0]),
    ([[0, 0, 0]], [0, 0, 0]),
    ([[1, 1, 0], [0, 1, 1]], [1, 0, -1]),
    ([[1, 1, 0], [0, 1, 1]], [1, 0, 0]),
    ([[0, 1, 2], [1, 0, 3], [1, 1, 5]], [1, -1, 1]),
    ([[0, 1, 2], [1, 0, 3], [1, 1, 5]], [1, -1, 0]),
    ([[1, 2, 3, 4], [2, 4, 6, 8]], [2, 4, 6, 8]),
    ([[1, 2, 3, 4], [2, 4, 6, 8]], [2, 4, 6, 7]),
])
def test_general_primal_or_dual_certificates(rows, target):
    rows = [list(map(F, row)) for row in rows]
    cert = build.certificate(rows, target)
    verify.check_certificate(rows, target, cert)
    forged = deepcopy(cert)
    if cert["kind"] == "decoder":
        # A decoder for one meaning cannot attest a different meaning.
        target = list(target)
        target[0] += 1
    else:
        forged["perturbation"] = ["0"]*len(target)
    with pytest.raises(ValueError):
        verify.check_certificate(rows, target, forged)


def test_incompatible_observation_is_not_an_ambiguity():
    with pytest.raises(ValueError, match="visible ambiguity"):
        verify.check_certificate([[F(1), F(0)]], [1, 0],
                                 {"kind": "ambiguity", "perturbation": ["1", "0"]})


@pytest.mark.parametrize("rows,target,cert", [
    ([[1, 42]], [1], {"kind": "decoder", "coefficients": ["1"], "sample_error_gain": "1"}),
    ([[0, 42]], [1], {"kind": "ambiguity", "perturbation": ["1"]}),
    ([[]], [1], {"kind": "ambiguity", "perturbation": ["1"]}),
    ([[1.0]], [1], {"kind": "decoder", "coefficients": ["1"], "sample_error_gain": "1"}),
    ([[True]], [1], {"kind": "decoder", "coefficients": ["1"], "sample_error_gain": "1"}),
    ([[1]], [True], {"kind": "decoder", "coefficients": ["1"], "sample_error_gain": "1"}),
])
def test_certificate_rejects_ragged_or_inexact_linear_problem(rows, target, cert):
    with pytest.raises(ValueError, match="exact"):
        verify.check_certificate(rows, target, cert)


def test_zero_dimensional_meaning_is_valid():
    assert verify.check_certificate([[]], [], {
        "kind": "decoder", "coefficients": ["0"], "sample_error_gain": "0"})


def test_all_small_two_column_problems_agree_with_independent_rank():
    for entries in product((-1, 0, 1), repeat=4):
        rows = [list(map(F, entries[:2])), list(map(F, entries[2:]))]
        for target in product((-1, 0, 1), repeat=2):
            cert = build.certificate(rows, target)
            assert verify.check_certificate(rows, target, cert) == verify.identifies_two(rows, target)


def test_aggregate_recovery_does_not_certify_individual_reads(packet):
    certs = packet["captured"][0]["prefix_certificates"][-1]
    assert [c["kind"] for c in certs] == ["ambiguity", "ambiguity", "decoder", "ambiguity"]
    forged = deepcopy(packet["captured"])
    forged[0]["prefix_certificates"][-1][0] = deepcopy(certs[2])
    with pytest.raises(ValueError, match="does not reconstruct"):
        verify.check_captured(forged)


@pytest.mark.parametrize("attack", ["coefficient", "gain", "future", "prefix", "word", "sample",
                                  "writer", "cost", "preparation", "missing", "extra"])
def test_captured_semantic_mutations(packet, attack):
    controls = deepcopy(packet["captured"])
    c = controls[1]
    decoder = c["prefix_certificates"][-1][0]
    if attack == "coefficient":
        decoder["coefficients"][-1] = str(F(decoder["coefficients"][-1])+1)
    elif attack == "gain":
        decoder["sample_error_gain"] = "0"
    elif attack == "future":
        c["prefix_certificates"][0][0] = deepcopy(decoder)
    elif attack == "prefix":
        c["prefix_certificates"].pop(2)
    elif attack == "word":
        c["word"].reverse()
    elif attack == "sample":
        c["cases"][0]["samples"][-1] = "42"
    elif attack == "writer":
        c["cases"][0]["tape_sha256"] = "0"*64
    elif attack == "cost":
        c["means_per_case"] -= 1
    elif attack == "preparation":
        c["full_support_preparation_writes"] = 5
    elif attack == "missing":
        controls.pop()
    else:
        c["m1_derived"] = True
    with pytest.raises(ValueError):
        verify.check_captured(controls)


def test_stopped_work_includes_abort_mass(packet):
    for row in packet["toy"]:
        for t in range(4):
            assert sum(row["first_counts"][t])+row["abort_counts"][t] == row["attempt_words"]
            exact = sum(k*n for k, n in enumerate(row["first_counts"][t]))
            exact += row["horizon"]*row["abort_counts"][t]
            assert row["stopped_mean_totals"][t] == exact
            assert F(row["expected_samples"][t]) == 1+F(exact, row["attempt_words"])
    assert packet["toy"][-1]["abort_counts"][0] > 0


def test_noise_gain_bounds_adversarial_sample_errors(packet):
    c = packet["captured"][1]
    coefficients = list(map(F, c["prefix_certificates"][-1][0]["coefficients"]))
    eps = F(1, 4096)
    noise = [eps if a > 0 else -eps if a < 0 else 0 for a in coefficients]
    gain = F(c["prefix_certificates"][-1][0]["sample_error_gain"])
    assert sum(a*e for a, e in zip(coefficients, noise)) == gain*eps
    for case in c["cases"]:
        samples = [F(v)-3+e for v, e in zip(case["samples"], noise)]
        error = abs(sum(a*v for a, v in zip(coefficients, samples))-case["payload"][0])
        assert error == gain*eps


@pytest.mark.parametrize("attack", ["scope", "pins", "abort_cost", "guard_cost", "erasure", "null", "duplicate"])
def test_real_cli_rejects_resealed_bad_packets(packet, tmp_path, attack):
    forged = deepcopy(packet)
    if attack == "scope":
        forged["m1_derived"] = True
    elif attack == "pins":
        forged["pins"] = {}
    elif attack == "abort_cost":
        forged["toy"][-1]["stopped_mean_totals"][0] -= 8*forged["toy"][-1]["abort_counts"][0]
    elif attack == "guard_cost":
        forged["guarded"]["horizons"][-1]["rejected_proposal_total"] = 0
    elif attack == "erasure":
        forged["erasure"]["sound_individual_abort_lower_bound"] = "0"
    path = tmp_path/"controls.json"
    raw = codec.canonical(forged)
    if attack == "null":
        raw = b"null\n"
    if attack == "duplicate":
        raw = raw.replace(b'{', b'{"m1_derived":false,', 1)
    path.write_bytes(raw)
    result = subprocess.run([sys.executable, "-m", "source_temporal_acceptance.verify",
                             "--controls", str(path)], cwd=codec.ROOT,
                            env=dict(os.environ, PYTHONPATH=str(codec.ROOT/"code")),
                            capture_output=True, text=True, timeout=120)
    assert result.returncode != 0, result.stdout
    assert "ValueError" in result.stderr, result.stderr


def test_complete_lean_gate_and_ci_entry():
    base = codec.ROOT/"Lean/Geometry"
    gate = (base/"SourceTemporalAcceptanceAxiomAudit.lean").read_text(encoding="utf-8")
    for path in base.glob("SourceTemporal*.lean"):
        if "AxiomAudit" in path.name:
            continue
        for name in re.findall(r"^theorem (\w+)", path.read_text(encoding="utf-8"), re.M):
            assert f"audit_reusable_bus_axioms OPH.{path.stem}.{name}" in gate
    workflow = (codec.ROOT/".github/workflows/lean.yml").read_text(encoding="utf-8")
    executable = "\n".join(line for line in workflow.splitlines() if not line.lstrip().startswith("#"))
    assert '"Geometry.SourceTemporalAcceptanceAxiomAudit"' in executable


@pytest.mark.parametrize("attack", ["missing", "duplicate", "remote", "unseen_relay", "unsupported"])
def test_topology_plan_semantics(attack):
    vertices = set(range(4))
    edges = {frozenset((0, 1)), frozenset((1, 2)), frozenset((2, 3))}
    paths = [[2, 3], [1, 2, 3], [0, 1, 2, 3]]
    if attack == "missing":
        paths.pop()
    elif attack == "duplicate":
        paths[2] = [0, 1, 2, 1, 2, 3]
    elif attack == "remote":
        paths[2] = [0, 1]
    elif attack == "unseen_relay":
        paths[0] = [1, 2, 3]
    else:
        paths[2] = [0, 2, 3]
    with pytest.raises(ValueError):
        check_tomography.check_paths(vertices, edges, 3, paths)


@pytest.mark.parametrize("edges,root", [
    ([(0, 1), (0, 2), (0, 3), (0, 4)], 0),
    ([(0, 1), (1, 2), (1, 3), (3, 4)], 2),
    ([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], 4),
])
def test_branching_native_completion_all_basis_vectors(edges, root):
    graph = {v: set() for v in range(5)}
    for a, b in edges:
        graph[a].add(b)
        graph[b].add(a)
    _, _, paths = tomography.plan(graph, root)
    check_tomography.check_paths(set(graph), {frozenset(e) for e in edges}, root, paths)
    for basis in range(5):
        initial = {v: F(v == basis) for v in graph}
        samples, final, _ = tomography.scalar_run(paths, root, initial)
        recovered, calibrated = check_tomography.recover(root, paths, samples)
        assert recovered == initial and calibrated == final


def test_initial_receiver_sample_is_a_necessary_guard_hypothesis():
    # The whole initial two-port state distinguishes the sources, but a mean
    # erases their difference unless the receiver's initial sample is retained.
    columns = []
    for basis in range(2):
        initial = {v: F(v == basis) for v in range(2)}
        samples, _, _ = tomography.scalar_run([[0, 1]], 1, initial)
        columns.append(samples)
        recovered, _ = check_tomography.recover(1, [[0, 1]], samples)
        assert recovered == initial
    rows = list(zip(*columns))
    trial_state = [rows[-1], rows[-1]]
    for target in ((1, 0), (0, 1)):
        assert verify.identifies_two(rows, target)
        assert not verify.identifies_two(trial_state, target)
        assert verify.identifies_two([rows[0]]+trial_state, target)


def test_every_connected_four_port_graph_all_receivers_and_basis_states():
    vertices = set(range(4))
    possible = list(combinations(range(4), 2))
    connected = 0
    for included in product((False, True), repeat=len(possible)):
        edges = [e for e, yes in zip(possible, included) if yes]
        graph = {v: set() for v in vertices}
        for a, b in edges:
            graph[a].add(b)
            graph[b].add(a)
        # Independent finite reachability check, separate from the planner.
        reached = {0}
        for _ in vertices:
            reached |= {b for a in reached for b in graph[a]}
        if reached != vertices:
            with pytest.raises(ValueError, match="disconnected"):
                tomography.plan(graph, 0)
            continue
        connected += 1
        for root in vertices:
            _, _, paths = tomography.plan(graph, root)
            work = check_tomography.check_paths(vertices, {frozenset(e) for e in edges}, root, paths)
            assert work <= 6 and len(paths) == 3
            for basis in vertices:
                initial = {v: F(v == basis) for v in vertices}
                samples, final, _ = tomography.scalar_run(paths, root, initial)
                recovered, calibrated = check_tomography.recover(root, paths, samples)
                assert recovered == initial and calibrated == final
    assert connected == 38


def test_single_port_completion_needs_only_the_initial_sample():
    _, _, paths = tomography.plan({0: set()}, 0)
    assert paths == []
    assert check_tomography.check_paths({0}, set(), 0, paths) == 0
    initial = {0: F(-7, 3)}
    samples, final, _ = tomography.scalar_run(paths, 0, initial)
    assert samples == [initial[0]]
    assert check_tomography.recover(0, paths, samples) == (initial, final)


@pytest.mark.parametrize("attack", ["sample", "tree", "cost", "precision", "radius", "replay"])
def test_completion_evidence_mutations(packet, attack):
    control = deepcopy(packet["tomography"])
    if attack == "sample":
        control["path_control"]["cases"][0]["samples"][-1] = "9"
    elif attack == "tree":
        control["captured_plan"]["tree_sha256"] = "0"*64
    elif attack == "cost":
        control["captured_plan"]["native_means"] -= 1
    elif attack == "precision":
        control["physical_precision_derived"] = True
    elif attack == "radius":
        control["metric_radius_selected"] = True
    else:
        control["captured_plan"]["payload_replay"] = True
    with pytest.raises(ValueError):
        check_tomography.check(control, verify.adjacency())


def test_continuation_is_stronger_than_waiting_or_current_read(packet):
    for row in packet["toy"]:
        for j in range(4):
            assert row["accepted_counts"][j] <= row["continuable_counts"][j]
            assert row["pending_completable_counts"][j] == (
                row["continuable_counts"][j]-row["accepted_counts"][j])
            assert row["irreversible_counts"][j]+row["continuable_counts"][j] == row["attempt_words"]
    first = packet["toy"][1]
    assert first["irreversible_counts"] == [1, 1, 0, 1]
    assert first["pending_completable_counts"] == [2, 2, 3, 2]


@pytest.mark.parametrize("bad", ["1e999999999999999999999", "NaN", "1/0", "2/4", "-0", "+1", "01"])
def test_reject_noncanonical_or_explosive_rational(bad):
    with pytest.raises(ValueError):
        verify.rational(bad)


def test_guard_rejects_erasure_and_charges_rejected_proposals(packet):
    idle, rescued, completion = packet["guarded"]["controls"]
    assert idle["admitted"] == [False]*6
    assert idle["means"] == 0 and idle["stopped_proposals"] == 6
    assert rescued["admitted"] == [False, True, True, True, True, True]
    assert rescued["first_complete"] == completion["first_complete"] == 6
    for row in packet["guarded"]["horizons"]:
        assert row["complete_count"]+row["pending_count"] == row["attempt_words"]
        assert row["executed_mean_total"]+row["rejected_proposal_total"] == row["stopped_proposal_total"]
        assert row["sample_total"] == row["attempt_words"]+row["executed_mean_total"]
        assert row["irreversible_count"] == 0


def test_universal_word_completes_after_every_small_guarded_prefix():
    completion = (2, 1, 2, 0, 1, 2)
    for length in range(5):
        for prefix in product(range(3), repeat=length):
            word = prefix+completion
            trace = verify.guarded_trace(word)
            assert trace["first_complete"] is not None
            assert trace["first_complete"] <= length+6
            assert build.guarded_word(word) == trace


@pytest.mark.parametrize("attack", ["allow_loss", "future_read", "missing_rejection", "row"])
def test_guard_control_forgeries(packet, attack):
    forged = deepcopy(packet["guarded"])
    control = forged["controls"][1]
    if attack == "allow_loss":
        control["admitted"][0] = True
    elif attack == "future_read":
        control["first_complete"] = 1
    elif attack == "missing_rejection":
        control["stopped_proposals"] -= 1
    else:
        control["receiver_rows"][-1] = ["1", "0"]
    with pytest.raises(ValueError, match="guard control"):
        verify.check_guarded(forged)
