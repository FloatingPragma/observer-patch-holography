import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from source_selection_model import channels, finite_model, geometry, grammar, verify_response

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


@pytest.fixture(scope="module")
def response():
    return verify_response.strict_load(HERE/"response.json")


@pytest.fixture(scope="module")
def records():
    return verify_response.strict_load(HERE/"records.json")


def test_response_is_complete_and_operational(response):
    result = verify_response.verify(response)
    assert result == {"response_rank": 12, "observable_response_rank": 12,
                      "block_diagonal_observable_control_rank": 11,
                      "ordered_brackets": 66, "cayley_factors": 20,
                      "closed_path_actions": 60, "response_naturality_pairs": 720,
                      "path_compositions": 3600}


@pytest.mark.parametrize("field", ["generators", "matrix_unit_response_probes",
                                     "ordered_mixed_responses", "cayley_factors", "closed_paths"])
def test_response_omissions_fail(response, field):
    bad = copy.deepcopy(response)
    bad[field].pop()
    with pytest.raises(ValueError):
        verify_response.verify(bad)


@pytest.mark.parametrize("target", ["generator", "tomography", "jet", "bracket", "control", "output", "path", "port"])
def test_response_forgery_is_not_a_passing_certificate(response, target):
    bad = copy.deepcopy(response)
    if target == "generator":
        bad["generators"][0][0][0][0][2] = "123"
    elif target == "tomography":
        bad["matrix_unit_response_probes"][0][-1][0] = "1"
    elif target == "jet":
        bad["ordered_mixed_responses"][0]["mixed_coefficient"][0][0][0][0] = "123"
    elif target == "bracket":
        bad["ordered_mixed_responses"][0]["in_port_basis"][0][0] = "123"
    elif target == "control":
        bad["cayley_factors"][0]["control"][0][0] = "123"
    elif target == "output":
        bad["cayley_factors"][0]["observed_response"][0][0][0][0] = "123"
    elif target == "path":
        bad["closed_paths"][-1]["factor_indices"] = []
    else:
        bad["closed_paths"][-1]["port_action"] = list(range(12))
    with pytest.raises(ValueError):
        verify_response.verify(bad)


@pytest.mark.parametrize("value", [True, 0.5, "01", "1/0", "NaN", "1"*128])
def test_noncanonical_or_unbounded_numbers_fail(response, value):
    bad = copy.deepcopy(response)
    bad["generators"][0][0][0][0][0] = value
    with pytest.raises((ValueError, ZeroDivisionError)):
        verify_response.verify(bad)


def test_central_record_interventions_and_selected_weights(records):
    assert finite_model.record_controls(records)["local_events"] == 6912
    assert finite_model.algebra_controls()["reference_eigenvalue"] == "1/72"
    assert finite_model.record_controls(records)["uniform_guess_failure"] == "11/12"


@pytest.mark.parametrize("field", range(7))
def test_every_local_record_field_is_checked(records, field):
    bad = copy.deepcopy(records)
    bad["tapes"][0][0][field] = "123" if field in (3, 4) else 123
    with pytest.raises(ValueError):
        finite_model.record_controls(bad)


@pytest.mark.parametrize("target", ["missing_case", "duplicate_case", "wrong_tape", "missing_event", "extra_metadata", "bool"])
def test_no_missing_intervention_or_unaccounted_channel(records, target):
    bad = copy.deepcopy(records)
    if target == "missing_case":
        bad["cases"].pop()
    elif target == "duplicate_case":
        bad["cases"][-1] = copy.deepcopy(bad["cases"][0])
    elif target == "wrong_tape":
        bad["cases"][1]["tapes"][1] = bad["cases"][0]["tapes"][0]
    elif target == "missing_event":
        bad["tapes"][0].pop()
    elif target == "extra_metadata":
        bad["tapes"][0][0].append("global-state-root")
    else:
        bad["cases"][0]["labels"][0] = False
    with pytest.raises(ValueError):
        finite_model.record_controls(bad)


def test_support_chain_maps_routes_and_mixture_weights():
    rows = geometry.tower()
    result = geometry.check_tower(rows)
    assert [r["carriers"] for r in result] == [12, 42, 162, 642]
    assert [r["tagged_scalar_seam_routes"] for r in result] == [0, 382, 1810, 7152]
    assert finite_model.refinement_controls(rows) == {
        "one_step_maps": 3, "composite_fibers": 12,
        "nonuniform_composite_fibers": 11, "selective_posterior": ["1/4", "3/4"]}


def test_all_omitted_state_coordinates_have_positive_hidden_interventions():
    assert grammar.controls()["omission_witnesses"] == 431
    rows = grammar.coordinates()
    for omitted in range(len(rows)):
        with pytest.raises(ValueError, match="grammar count"):
            grammar.verify_coordinates(rows[:omitted]+rows[omitted+1:])


@pytest.mark.parametrize("target", ["duplicate", "bool", "imaginary_diagonal", "foreign_block"])
def test_false_constraint_basis_fails(target):
    rows = grammar.coordinates()
    if target == "duplicate":
        rows[-1] = rows[0]
    elif target == "bool":
        rows[0][0] = False
    elif target == "imaginary_diagonal":
        rows[0][3] = 1
    else:
        rows[0][0] = 12
    with pytest.raises(ValueError):
        grammar.verify_coordinates(rows)


def test_channel_invertibility_is_not_exact_readability():
    result = channels.controls()
    assert result["identity"]["decoder"] == list(range(12))
    assert result["redundant_outputs"]["decoder"] == [x//2 for x in range(24)]
    assert result["full_rank_noisy"]["decoder"] is None
    assert result["full_rank_noisy"]["signed_inverse_is_stochastic"] is False


@pytest.mark.parametrize("rows", [[], [[]], [[True]], [["01"]], [["-1", "2"]],
                                  [["1/2"]], [["1"], ["1/2", "1/2"]]])
def test_bad_channel_probabilities_fail(rows):
    with pytest.raises(ValueError):
        channels.classify(rows)


@pytest.mark.parametrize("target", ["missing_face", "reverse_face", "collapse", "nonlocal", "population"])
def test_bad_geometric_bridge_rejected(target):
    rows = geometry.tower()
    if target == "missing_face":
        rows[1]["faces"].pop()
    elif target == "reverse_face":
        rows[1]["faces"][0] = tuple(reversed(rows[1]["faces"][0]))
    elif target == "collapse":
        rows[1]["coarsen"] = [0]*42
    elif target == "nonlocal":
        rows[1]["coarsen"][-1] = 0
    else:
        rows[1]["vertices"] = 43
    with pytest.raises(ValueError):
        geometry.check_tower(rows)


def test_nonicosahedral_sphere_with_correct_census_is_rejected():
    # An edge flip preserves V/E/F, spherical topology and every old chain
    # and routing check. Its four affected seed valences are 4,4,6,6.
    faces = geometry.seed()
    faces.remove((0, 1, 2))
    faces.remove((0, 3, 1))
    faces.extend([(2, 0, 3), (2, 3, 1)])
    count = 12
    rows = [{"level": 0, "vertices": count, "faces": faces, "coarsen": None}]
    for level in range(1, 4):
        count, faces, coarse = geometry.subdivide(count, faces)
        rows.append({"level": level, "vertices": count, "faces": faces, "coarsen": coarse})
    with pytest.raises(ValueError, match="icosahedral seed"):
        geometry.check_tower(rows)


def test_fine_edge_flip_cannot_substitute_another_spherical_refinement():
    rows = geometry.tower()
    faces = rows[1]["faces"]
    a, b, c = faces[0]
    j = next(i for i, f in enumerate(faces) if i != 0 and b in f and c in f)
    other = faces[j]
    d = next(v for v in other if v not in (b, c))
    # The common oriented edge is b -> c in the first face.
    faces[0], faces[j] = (a, b, d), (a, d, c)
    with pytest.raises(ValueError, match="midpoint"):
        geometry.check_tower(rows)


def test_checker_does_not_import_response_or_record_producers():
    program = '''
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name in {"source_selection_model.response", "source_selection_model.records"}:
            raise RuntimeError("producer import forbidden")
sys.meta_path.insert(0, Block())
from source_selection_model.verify import evidence
assert evidence()["response"]["observable_response_rank"] == 12
'''
    env = dict(os.environ, PYTHONPATH=str(ROOT/"code"))
    result = subprocess.run([sys.executable, "-c", program], cwd=ROOT, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr


def test_actual_cli_and_resealed_false_receipt(tmp_path):
    env = dict(os.environ, PYTHONPATH=str(ROOT/"code"))
    cmd = [sys.executable, "-m", "source_selection_model.verify"]
    good = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)
    assert good.returncode == 0, good.stdout+good.stderr
    bad = json.loads((HERE/"receipt.json").read_text(encoding="utf-8"))
    bad["response"]["observable_response_rank"] = 11
    # There is no trusted envelope hash that can make this semantic change valid.
    path = tmp_path/"forged.json"
    path.write_text(json.dumps(bad), encoding="utf-8")
    result = subprocess.run(cmd+["--receipt", str(path)], cwd=ROOT, env=env,
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "disagrees with independent replay" in result.stderr


def test_duplicate_json_keys_rejected(tmp_path):
    path = tmp_path/"duplicate.json"
    path.write_text('{"schema":1,"schema":2}', encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON"):
        verify_response.strict_load(path)


def test_capture_producers_reproduce_retained_packets():
    from source_selection_model import records, response
    assert json.loads(json.dumps(records.capture())) == verify_response.strict_load(HERE/"records.json")
    assert response.source_packet() == verify_response.strict_load(HERE/"response.json")


def test_source_edit_invalidates_actual_cli_receipt(tmp_path):
    import shutil
    from source_selection_model.verify import pins
    clone = tmp_path/"source"
    for relative in list(pins()) + ["code/source_selection_model/receipt.json"]:
        target = clone/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/relative, target)
    changed = clone/"code/source_selection_model/DERIVATION.md"
    changed.write_bytes(changed.read_bytes()+b"\nChanged model membership.\n")
    env = dict(os.environ, PYTHONPATH=str(clone/"code"))
    result = subprocess.run([sys.executable, "-m", "source_selection_model.verify"],
                            cwd=clone, env=env, capture_output=True, text=True)
    assert result.returncode != 0
    assert "disagrees with independent replay" in result.stderr


def test_lean_gate_covers_every_declaration():
    import re
    source = (ROOT/"Lean/Geometry/SourceSelectionLocality.lean").read_text(encoding="utf-8")
    gate = (ROOT/"Lean/Geometry/SourceSelectionLocalityAxiomAudit.lean").read_text(encoding="utf-8")
    names = re.findall(r"^theorem (\w+)", source, re.M)
    assert len(names) == 15
    assert {"OPH.SourceSelectionLocality."+n for n in names} == set(
        re.findall(r"^audit_reusable_bus_axioms (\S+)$", gate, re.M))
    ci = (ROOT/".github/workflows/lean.yml").read_text(encoding="utf-8")
    assert '"Geometry.SourceSelectionLocalityAxiomAudit"' in ci
