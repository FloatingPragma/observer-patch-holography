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
