"""Mutation and independence checks for the positive-capacity top range."""
import copy
import importlib.util
import json
from pathlib import Path
import sys

import pytest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import top_positive_capacity_band as producer
import verify_top_positive_capacity_band as verifier


@pytest.fixture
def receipt():
    return verifier.strict_load(verifier.OUT)


def test_independent_full_receipt_and_symbolic_identities(receipt):
    result=verifier.verify_payload(receipt)
    assert result['verified']
    assert result['one_loop_controls']==result['two_loop_controls']==5
    assert not result['formal_theorem_proof_claimed']
    assert not result['physical_matching_or_two_loop_interval_certified']


def test_no_equal_capacity_or_two_loop_band_promotion(receipt):
    assert not receipt['scope']['equal_capacities_assumed']
    assert not receipt['scope']['two_loop_continuum_enclosure_claimed']
    assert receipt['one_loop_controls'][0]['weight_E_over_sum']==0
    assert receipt['one_loop_controls'][-1]['weight_E_over_sum']==1
    assert receipt['one_loop_monotone_band_endpoints']==[
        receipt['one_loop_controls'][0]['mass_coordinate'],
        receipt['one_loop_controls'][-1]['mass_coordinate']]


@pytest.mark.parametrize('field',[
    'equal_capacities_assumed','mass_targets_consumed','quark_template_consumed',
    'controls_used_for_model_selection','numerical_endpoint_enclosure_certified',
    'two_loop_continuum_enclosure_claimed','theory_error_bound_claimed',
    'physical_pole_prediction_claimed','source_selected_criticality_claimed',
    'source_selected_anchor_records_claimed','source_only_alpha_claimed',
    'blind_prediction_claimed'])
def test_scope_promotions_rejected(receipt,field):
    receipt['scope'][field]=True
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('field',['P','E_star','alpha_U','v','mz','mu_U','E_cell'])
def test_mixed_source_packet_rejected(receipt,field):
    receipt['packet'][field]*=1.01
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('field',[
    'positive_inverse_gauge_minimum','positive_sensitivity_coefficient',
    'negative_Cprime_coefficient','uniform_lower_y_at_50','required_y_at_50',
    'uniform_upper_y_at_mu_U','required_y_at_mu_U'])
def test_false_analytic_margin_rejected(receipt,field):
    receipt['proof_checks'][field]*=-1
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('field',['c_over_b_on_anchor_endpoints','c_over_a_on_domain_endpoints'])
def test_endpoint_hypothesis_failure_rejected(receipt,field):
    receipt['proof_checks'][field][-1]=.49
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('collection',['one_loop_controls','coupled_two_loop_controls'])
def test_controls_cannot_be_dropped(receipt,collection):
    receipt[collection].pop(0)
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('collection',['one_loop_controls','coupled_two_loop_controls'])
def test_capacity_choice_cannot_be_tuned(receipt,collection):
    receipt[collection][1]['weight_E_over_sum']=.3
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('endpoint',[0,1])
def test_band_cannot_be_tuned(receipt,endpoint):
    receipt['one_loop_monotone_band_endpoints'][endpoint]+=1
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('field',['mass_coordinate','QCD_pole_proxy','Higgs_tree_proxy'])
def test_two_loop_mass_control_tampering_rejected(receipt,field):
    receipt['coupled_two_loop_controls'][2]['rows'][1][field]+=1
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


def test_quadrature_ode_disagreement_rejected(receipt):
    receipt['one_loop_controls'][1]['independent_ODE_mass']+=.001
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('field',['gauge_anchor_max_residual','criticality_residual','self_scale_residual'])
def test_bad_coupled_residual_rejected(receipt,field):
    receipt['coupled_two_loop_controls'][1]['rows'][0][field]=.1
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


@pytest.mark.parametrize('value',[float('nan'),float('inf'),True,'154.8'])
def test_non_numeric_or_nonfinite_mass_rejected(receipt,value):
    receipt['one_loop_controls'][0]['mass_coordinate']=value
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


def test_platform_tolerant_float_replay(receipt):
    row=receipt['one_loop_controls'][2]
    row['mass_coordinate']+=1e-9
    row['independent_ODE_mass']+=1e-9
    row['ODE_difference']=row['independent_ODE_mass']-row['mass_coordinate']
    assert verifier.verify_payload(receipt)['verified']


def test_hidden_mass_target_input_rejected_by_spec_pin(tmp_path,monkeypatch):
    spec=verifier.reviewed_spec();spec['measured_top_mass_GeV']=172.5
    path=tmp_path/'poisoned.json';path.write_text(json.dumps(spec))
    monkeypatch.setattr(producer,'SPEC',path)
    with pytest.raises(ValueError,match='unreviewed specification'):producer.reviewed_spec()


def test_changed_capacity_menu_rejected_by_spec_pin(tmp_path,monkeypatch):
    spec=verifier.reviewed_spec();spec['anchors']['control_fractions']=[.5]
    path=tmp_path/'poisoned.json';path.write_text(json.dumps(spec))
    monkeypatch.setattr(producer,'SPEC',path)
    with pytest.raises(ValueError,match='unreviewed specification'):producer.reviewed_spec()


def test_changed_alpha_calibration_rejected(tmp_path,monkeypatch):
    fixture={'inverse_fine_structure_constant':{'value':'137.0','standard_uncertainty':'0.1'}}
    path=tmp_path/'poisoned.json';path.write_text(json.dumps(fixture))
    monkeypatch.setattr(producer,'ALPHA',path)
    with pytest.raises(ValueError,match='alpha fixture'):producer.source_packet(producer.reviewed_spec())


@pytest.mark.parametrize('method',[
    'diagonal_quark_masses','charged_lepton_masses','alpha_external_from_d10',
    'structured_thomson_running','structured_thomson_running_asymptotic','solve_closure'])
def test_gauge_context_forbids_quark_and_legacy_closures(method):
    ctx=object.__new__(producer.GaugeOnlyContext)
    with pytest.raises(ValueError,match='forbidden'):getattr(ctx,method)()


def test_flavor_initialization_suppressed():
    ctx=object.__new__(producer.GaugeOnlyContext)
    assert ctx._derive_stage5_integer_vectors()=={}
    verifier.source_ancestry_check()


def test_only_alpha_and_spec_parsed_and_target_fixture_never_read(monkeypatch):
    """Poison every other parsed text input during gauge-only construction."""
    original=Path.read_text
    allowed={producer.SPEC.resolve(),producer.ALPHA.resolve()}
    def read(path,*args,**kwargs):
        if path.resolve() not in allowed:
            raise AssertionError('Undeclared numeric input: '+str(path))
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'read_text',read)
    # Do not rerun the slow high-precision root: mock only its final output
    # while checking the actual fixed-P argument and constructor suppression.
    expected=json.loads(original(verifier.OUT))['packet']
    from types import SimpleNamespace
    def build(ctx,p):
        assert abs(float(p)-expected['P'])<1e-14
        return SimpleNamespace(alpha_u=expected['alpha_U'],v=expected['v']/expected['E_star'],
            mz_run=expected['mz']/expected['E_star'],alpha_y_mz=expected['low_g_squared'][0]/(4*producer.pi),
            alpha2_mz=expected['low_g_squared'][1]/(4*producer.pi),alpha3_mz=expected['low_g_squared'][2]/(4*producer.pi))
    monkeypatch.setattr(producer.GaugeOnlyContext,'build_d10_from_p',build)
    def poison(*args,**kwargs):raise AssertionError('Quark template or closure was consumed')
    for method in ['diagonal_quark_masses','charged_lepton_masses','alpha_external_from_d10',
                   'structured_thomson_running','structured_thomson_running_asymptotic','_derive_stage5_integer_vectors']:
        monkeypatch.setattr(producer.PaperMathContext,method,poison)
    result=producer.source_packet(producer.reviewed_spec())
    assert abs(result['v']-expected['v'])<1e-10


def test_extra_target_field_or_missing_guard_rejected(receipt):
    receipt['packet']['observed_top']=172.5
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


def test_historical_or_comparison_source_pin_rejected(receipt):
    receipt['source_pins']['code/particles/flavor/compare_conditional_quark_masses.py']='0'*64
    with pytest.raises(ValueError):verifier.verify_payload(receipt)


def test_verifier_does_not_import_producer_or_beta_engine():
    tree=verifier.ast.parse(Path(verifier.__file__).read_text())
    imports=[]
    for node in verifier.ast.walk(tree):
        if isinstance(node,verifier.ast.Import):imports.extend(a.name for a in node.names)
        elif isinstance(node,verifier.ast.ImportFrom):imports.append(node.module)
    assert not set(imports)&{'top_positive_capacity_band','conditional_top_fixed_p','paper_math',
                             'sm_two_loop_rge_engine','derive_d11_criticality_boundary_scan'}
