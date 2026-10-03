"""Independent replay and adversarial controls for the conditional top result."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import conditional_top_fixed_p as producer
import verify_conditional_top_fixed_p as verifier


@pytest.fixture(scope='module')
def receipt():
    return verifier.strict_load(producer.OUT_PATH)


def test_retained_three_truncations_replay_independently(receipt):
    checked = verifier.verify_payload(receipt)
    assert checked['branches_checked'] == 3
    assert checked['physical_matching_or_theory_uncertainty_certified'] is False


def test_forward_rebuild_does_not_consume_mass_benchmarks(receipt, monkeypatch):
    # Benchmark endpoints and the archived flavor-dependent root are deliberately
    # unusable. The forward run must use its independently supplied alpha alone.
    monkeypatch.setattr(producer.rge, 'COUPLINGS_AT_MT', (float('nan'),)*5)
    monkeypatch.setattr(producer.rge, 'BENCHMARK_AT_MPL', (float('nan'),)*5)
    monkeypatch.setattr(producer.rge, 'MT', float('nan'))
    monkeypatch.setattr(producer.rge, 'MPL', float('nan'))
    monkeypatch.setattr(producer.rge, 'validate', producer.forbidden)
    monkeypatch.setattr(producer.legacy, 'ARCHIVED_PIN', {})
    monkeypatch.setattr(producer.legacy, 'P_FALLBACK', float('nan'))
    monkeypatch.setattr(producer.legacy, 'ALPHA_U_FALLBACK', float('nan'))
    for method in producer.FORBIDDEN_METHODS:
        if hasattr(producer.PaperMathContext, method):
            monkeypatch.setattr(producer.PaperMathContext, method, producer.forbidden)
    rebuilt = producer.build_payload()
    # Provenance and the physical contract are exact. Floating-point ODE
    # outputs may vary slightly across SciPy/platform builds; their replay
    # tolerance is still much tighter than the independent mass check.
    for name in ('schema', 'specification', 'scope', 'source_pins',
                 'P_from_independent_alpha', 'D10'):
        assert rebuilt[name] == receipt[name]
    assert_rebuild_matches(rebuilt, receipt)
    assert verifier.verify_payload(rebuilt)['verified']


def assert_rebuild_matches(actual, expected, path='receipt'):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys(), path
        for key in expected:
            assert_rebuild_matches(actual[key], expected[key], f'{path}.{key}')
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for index, (value, reference) in enumerate(zip(actual, expected)):
            assert_rebuild_matches(value, reference, f'{path}[{index}]')
    elif isinstance(expected, float):
        assert not isinstance(actual, bool), path
        assert actual == pytest.approx(expected, abs=2e-8, rel=2e-11), path
    else:
        assert actual == expected, path


@pytest.mark.parametrize('path,value', [
    (('scope','physical_pole_prediction_claimed'), True),
    (('scope','all_six_quark_masses_claimed'), True),
    (('scope','source_selected_boundary_claimed'), True),
    (('scope','source_selected_criticality_claimed'), True),
    (('scope','external_Estar_scale_calibration'), False),
    (('scope','mass_targets_consumed'), True),
    (('scope','quark_template_consumed'), True),
    (('scope','numerical_enclosure_certified'), True),
    (('scope','theory_uncertainty_quantified'), True),
    (('scope','blind_prediction_claimed'), True),
    (('specification','boundary','root'), 'any positive solution'),
    (('specification','boundary','global_polynomial_uniqueness'), True),
    (('specification','transport','hypercharge'), 'g1_GUT=gY'),
    (('specification','transport','potential'), 'V=lambda*(Hdagger*H)^2/2'),
    (('specification','transport','threshold_complete'), True),
    (('specification','readout','physical_MSbar_identification'), 'already matched exactly'),
    (('specification','readout','n_l'), 6),
    (('specification','readout','QCD_order'), 4),
    (('specification','readout','electroweak_pole_matching'), True),
    (('specification','readout','v_chart'), 'measured Fermi constant input'),
    (('P_from_independent_alpha',), '1.630972095858897'),
    (('D10','alpha3_mz'), '0.2'),
    (('source_scales_GeV','E_star'), 1.0),
    (('source_scales_GeV','v_transmutation_gev'), 174.0),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'top_running_mass_coordinate_GeV'), 172.0),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'top_QCD_pole_proxy_GeV'), 162.059468),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'yukawa_self_scale'), 0.92932086*2**.5),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'boundary_top_yukawa'), 5.64560),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'critical_polynomial_roots_yt_squared'), [0.15150913813]),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'beta_lambda_boundary'), 1e-3),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'low_anchor_gauge_residual'), [0.,0.,1e-3]),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'top_fixed_point_residual_GeV'), 1e-3),
    (('fully_coupled_two_loop_same_low_gauge_anchor',1,'top_running_coordinate_over_Estar'), 1e-10),
    (('numerical_diagnostics','interpretation'), 'The numerical interval proves agreement with nature.'),
])
def test_corrupted_physics_numerics_and_scope_rejected(receipt, path, value):
    mutant = copy.deepcopy(receipt)
    cursor = mutant
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises((ValueError, KeyError)):
        verifier.verify_payload(mutant)


def test_nonperturbative_root_cannot_be_relabelled_as_selected(receipt):
    mutant = copy.deepcopy(receipt)
    row = mutant['fully_coupled_two_loop_same_low_gauge_anchor'][1]
    row['boundary_top_yukawa'] = row['critical_polynomial_roots_yt_squared'][2]**.5
    row['beta_lambda_boundary'] = 0.0
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


def test_all_branch_receipts_are_required(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['hybrid_rows'].pop('1')
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


def test_numerical_certainty_cannot_be_added_as_extra_metadata(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['physical_mass_prediction'] = {'certified': True}
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


def test_stale_source_pin_is_rejected(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['source_pins'][producer.PIN_PATHS[0]] = '0'*64
    with pytest.raises(ValueError, match='source pin'):
        verifier.verify_payload(mutant)


@pytest.mark.parametrize('bad', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}'])
def test_strict_json_rejects_duplicates_and_nonfinite(tmp_path, bad):
    path = tmp_path/'bad.json'
    path.write_text(bad)
    with pytest.raises(ValueError):
        verifier.strict_load(path)


def test_cli_verifier_runs_with_all_producer_imports_blocked():
    program = '''
import builtins, runpy, sys
original = builtins.__import__
def guard(name, *args, **kwargs):
    if name.split('.')[-1] in {'conditional_top_fixed_p','paper_math','derive_d11_criticality_boundary_scan','sm_two_loop_rge_engine'}:
        raise AssertionError('producer import attempted: '+name)
    return original(name,*args,**kwargs)
builtins.__import__ = guard
sys.argv = [sys.argv[1]]
runpy.run_path(sys.argv[0],run_name='__main__')
'''
    result = subprocess.run([sys.executable, '-c', program, str(HERE/'verify_conditional_top_fixed_p.py')],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)['verified'] is True
