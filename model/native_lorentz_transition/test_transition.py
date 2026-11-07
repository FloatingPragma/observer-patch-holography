from fractions import Fraction

import pytest

from . import transition as t
from .certificate import _encode, build


def test_existing_native_kernel_replays():
    result = t.source_contract()
    assert result['exact_design_rank'] == 9
    assert result['exact_kernel_dimension'] == 1
    assert result['registered_binding']['registered_gram_entries_checked_exactly'] == 144


def test_cached_source_cannot_hide_changed_retained_bytes(tmp_path, monkeypatch):
    import shutil
    t.source_contract()
    for name in ('oph_exact/carrier.py', 'oph_fpe/dynamics/self_readback_repair_closure.py'):
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(t.ARCHIVE / name, target)
    target = tmp_path / 'oph_exact/carrier.py'
    target.write_text(target.read_text() + '\n# altered retained source\n')
    monkeypatch.setattr(t, 'ARCHIVE', tmp_path)
    with pytest.raises(ValueError, match='pinned native geometry'):
        t.source_contract()


@pytest.mark.parametrize('name', ['identity', 'boost', 'rotation', 'native_cyclic_rotation', 'dilation'])
def test_positive_controls(name):
    rows = t.controls()[name]
    result = t.proper_future(rows)
    assert result['future_preserving'] and result['proper']
    assert all(value == t.q(0) for value in t.null_charges(rows))
    assert t.pullback(rows) == t.scaled(result['metric_factor'], t.metric())
    assert result['metric_factor'] == t.q(Fraction(9, 4) if name == 'dilation' else 1)


def test_parity_satisfies_time_oriented_lean_contract_without_properness():
    rows = t.controls()['parity']
    result = t.classify(rows)
    assert result['future_preserving'] and not result['proper']
    assert result['determinant'] == t.q(-1)
    with pytest.raises(ValueError, match='spacetime orientation'):
        t.proper_future(rows)


def test_time_reversal_passes_null_equations_but_fails_future_test():
    rows = t.controls()['time_reversal']
    assert all(value == t.q(0) for value in t.null_charges(rows))
    assert not t.classify(rows)['future_preserving']
    with pytest.raises(ValueError, match='future cone'):
        t.proper_future(rows)


def test_singular_null_collapse_requires_explicit_invertibility():
    rows = t.controls()['null_collapse']
    assert all(value == t.q(0) for value in t.null_charges(rows))
    assert t.determinant(rows) == t.q(0)
    with pytest.raises(ValueError, match='invertible'):
        t.classify(rows)


def test_nonconformal_perturbation_fails_without_fitting():
    rows = [list(row) for row in t.boost()]
    rows[1][2] = t.q(Fraction(1, 1000))
    with pytest.raises(ValueError, match='null-direction'):
        t.classify(rows)
    result = t.metric_residual(rows)
    assert result['maximum_metric_residual'] != t.q(0)
    assert t.source().q5_sign(t.source().q5_sub(result['bound'], result['maximum_metric_residual'])) >= 0


def test_noncommuting_composition_retains_multiplicative_scale():
    examples = t.controls()
    first = t.scaled(Fraction(3, 2), examples['boost'])
    second = t.scaled(Fraction(2, 3), examples['rotation'])
    assert t.multiply(first, second) != t.multiply(second, first)
    result = t.compose(first, second)
    assert result['metric_factor'] == t.q(1)
    assert result['matrix'] == t.multiply(second, first)
    assert result['future_preserving'] and result['proper']


def test_orientation_signs_compose_separately():
    examples = t.controls()
    assert t.compose(examples['parity'], examples['parity'])['proper']
    assert t.compose(examples['time_reversal'], examples['time_reversal'])['future_preserving']
    both = t.compose(examples['parity'], examples['time_reversal'])
    assert both['proper'] and not both['future_preserving']


@pytest.mark.parametrize('factor', [0, -1, Fraction(3, 7)])
def test_general_form_kernel_does_not_infer_positive_scale(factor):
    assert t.conformal_form(t.scaled(factor, t.metric())) == t.q(factor)


def test_nonsymmetric_form_is_rejected():
    rows = [list(row) for row in t.metric()]
    rows[1][2] = t.q(1)
    with pytest.raises(ValueError, match='symmetric'):
        t.conformal_form(rows)


@pytest.mark.parametrize('value', [True, 0.1, (1.0, 0), (1, 2, 3), '1'])
def test_inexact_or_malformed_scalar_rejected(value):
    with pytest.raises(ValueError):
        t.q(value)


def test_bad_matrix_shape_rejected():
    with pytest.raises(ValueError, match='four by four'):
        t.classify([[1, 0], [0, 1]])


def test_boost_parameter_boundary_rejected():
    with pytest.raises(ValueError, match='future boost'):
        t.boost(1)


def test_certificate_is_json_stable_and_explicitly_conditional():
    import json
    result = build()
    assert json.loads(json.dumps(result)) == result
    assert sorted(result['native_cyclic_port_permutation']) == list(range(12))
    assert result['scope']['linear_transition_supplied']
    assert not result['scope']['source_selected_boost']
    assert result['all_pass']
    assert _encode({'a': (Fraction(1, 2), Fraction(0))}) == {'a': ['1/2', '0']}
