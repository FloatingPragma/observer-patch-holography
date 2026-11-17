"""The current suite cannot silently omit or reassign a completion test file."""
from pathlib import Path

import pytest

from tools import completion_test_scope as scope


def test_all_current_completion_files_have_one_domain():
    selected = scope.paths()
    assert 'model/completion/test_assumption_closure.py' in selected
    assert 'model/completion/test_theorem_map.py' in selected
    assert 'model/completion/test_acceptance_audit.py' in selected
    assert 'model/completion/test_observer.py' in selected
    assert 'model/completion/test_modular_access.py' in selected
    assert 'model/completion/test_particle_premise_joins.py' in selected
    assert 'model/completion/test_full_field_joins.py' in selected
    assert 'model/completion/test_four_lane_joins.py' in selected
    assert 'model/completion/test_field_current_measurement_consumers.py' in selected
    assert 'model/completion/test_program_class_v2.py::test_exact_comparison_rule_retains_inconclusive_and_calibrated_cases' in selected
    assert 'model/completion/test_program_class_v2.py::test_actual_continuing_observer_counts_fit_declared_finite_request_budget' in selected
    assert not any(name.removeprefix('model/completion/') in scope.COMPONENT_ROUTING
                   for name in selected)


def test_new_test_requires_explicit_source_domain(tmp_path):
    directory = tmp_path / 'model/completion'
    directory.mkdir(parents=True)
    for path in (scope.ROOT / 'model/completion').glob('test_*.py'):
        (directory / path.name).touch()
    (directory / 'test_new_operator.py').touch()
    with pytest.raises(ValueError, match='explicit source domain'):
        scope.paths(tmp_path)


def test_duplicate_domain_is_rejected(monkeypatch):
    monkeypatch.setattr(scope, 'CURRENT', (*scope.CURRENT, scope.COMPONENTS[0]))
    with pytest.raises(ValueError, match='explicit source domain'):
        scope.paths()


def test_model_suite_covers_every_model_collection_unit():
    selected = scope.model_paths()
    for path in (scope.ROOT / 'model').rglob('test_*.py'):
        relative = path.relative_to(scope.ROOT / 'model')
        if relative.parts[0] == 'completion' or any(
                part.startswith('.') or part == '__pycache__' for part in relative.parts):
            continue
        assert 'model/' + relative.parts[0] in selected, str(path)
    assert len(selected) == len(set(selected))
    assert 'model/ccf_lt2' in selected
    assert 'model/completion/test_assumption_closure.py' in selected
    assert 'model/completion' not in selected


def test_new_model_package_and_root_test_are_automatically_covered(tmp_path):
    completion = tmp_path / 'model/completion'
    completion.mkdir(parents=True)
    for path in (scope.ROOT / 'model/completion').glob('test_*.py'):
        (completion / path.name).touch()
    novel = tmp_path / 'model/new_quantum_operator/nested'
    novel.mkdir(parents=True)
    (novel / 'test_exact_current.py').touch()
    (tmp_path / 'model/test_new_gauss.py').touch()
    custom = tmp_path / 'model/custom_collection'
    custom.mkdir()
    (custom / 'conftest.py').touch()
    selected = scope.model_paths(tmp_path)
    assert 'model/new_quantum_operator' in selected
    assert 'model/test_new_gauss.py' in selected
    assert 'model/custom_collection' in selected


def test_only_explicitly_isolated_custody_node_can_be_deselected():
    command = scope.command(model_suite=True, shard_index=3, shard_count=16,
                            deselect=[scope.ISOLATED_CUSTODY_TEST])
    assert command.count('--deselect') == 1
    assert command[command.index('--deselect')+1] == scope.ISOLATED_CUSTODY_TEST
    assert command[command.index('--model-shard-index')+1] == '3'
    assert command[command.index('--model-shard-count')+1] == '16'
    assert 'model/charged_lepton_response' in command
    assert 'model/ccf_lt2' in command
    with pytest.raises(ValueError, match='separately executed'):
        scope.command(model_suite=True, deselect=['model/ccf_lt2/test_integrity.py::test_source'])
    with pytest.raises(ValueError, match='separately executed'):
        scope.command(deselect=[scope.ISOLATED_CUSTODY_TEST])
