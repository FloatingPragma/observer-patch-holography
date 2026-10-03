"""Target isolation, finite matching conventions, and numerical guard tests."""
from copy import deepcopy
import json
from pathlib import Path, PureWindowsPath
import sys
import math

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parent))
import source_ew_vev_matching as producer
import verify_source_ew_top_pole as top

@pytest.fixture(scope='module')
def receipt():return json.loads(producer.RECEIPT.read_text())


@pytest.mark.parametrize('implementation', ['producer', 'independent_verifier'])
def test_windows_relative_paths_preserve_receipt_keys(receipt, monkeypatch, implementation):
    import verify_source_ew_vev_matching as independent
    original = Path.relative_to
    # Exercise native Windows separators while reading the same real files.
    def windows_relative(path, *args, **kwargs):
        return PureWindowsPath(*original(path, *args, **kwargs).parts)
    with monkeypatch.context() as patch:
        patch.setattr(Path, 'relative_to', windows_relative)
        if implementation == 'producer':
            assert producer.source_pins() == receipt['source_pins']
        else:
            assert independent.verify(receipt)['matching_rows_checked'] == 9


def test_top_and_tadpole_reconstructed_independently(receipt):
    assert top.verify(receipt)


def test_target_and_old_coordinate_poison_does_not_change_inputs():
    source=json.loads((producer.ROOT/producer.load_spec()['source_receipt']).read_text())
    expected=producer.build_inputs(source)
    keep={'boundary_gauges','boundary_top_yukawa'}
    for branch in source['fully_coupled_two_loop_same_low_gauge_anchor']:
        for key in list(branch):
            if key not in keep:branch[key]='POISON_OBSERVED_OR_DERIVED_MASS_COORDINATE'
    source['hybrid_rows']='POISON';source['numerical_diagnostics']='POISON'
    source['experimental_quark_masses']={'t':1e99,'u':-1e99}
    assert producer.build_inputs(source)==expected


@pytest.mark.parametrize('field',[
    'quark_mass_targets_consumed','electroweak_mass_or_GF_targets_consumed',
    'blind_prediction_claimed','source_derived_vev_attachment_claimed',
    'complete_flavor_model_claimed','physical_nonzero_light_quark_masses_derived',
    'finite_muon_phase_space_term_retained','numerical_enclosure_certified',
    'theory_uncertainty_quantified','scale_controls_are_confidence_interval'])
def test_promotion_fails(receipt,field):
    changed=deepcopy(receipt);changed['scope'][field]=True
    with pytest.raises(ValueError):producer.validate_cached(changed)


@pytest.mark.parametrize('mutate',[
    lambda d:d['forward']['rows'].pop(),
    lambda d:d['forward']['rows'][0]['working'].__setitem__('v',250),
    lambda d:d['forward']['rows'][2]['top_QCD4_EW2'].__setitem__('mass_GeV',172),
    lambda d:d['forward']['rows'][2].__setitem__('GF_local_zero_muon_GeVm2',1.16e-5),
    lambda d:d['inputs'].__setitem__('Q0_GeV',91.1876),
    lambda d:d['inputs'].__setitem__('v_Q0_GeV',246.22),
    lambda d:d['source_pins'].pop(next(iter(d['source_pins']))),
])
def test_receipt_or_source_boundary_mutation_fails(receipt,mutate):
    changed=deepcopy(receipt);mutate(changed)
    with pytest.raises(ValueError):producer.validate_cached(changed)


@pytest.mark.parametrize('path,value,changed',[
    ('root/GF_local_zero_muon_GeVm2',1.18e-5,1.180001e-5),
    ('root/GF_removed_finite_muon_term_GeVm2',6e-12,6.1e-12),
    ('root/gY',.35,.3500001),
    ('root/mass_GeV',170.,170.001),
    ('root/matching_order',1,True),
    ('root/value',1.,True),
    ('root/matching_order',1,1.000000001),
    ('root/value',1.,float('nan')),
])
def test_replay_comparison_dimension_and_type_guards(path,value,changed):
    with pytest.raises(ValueError):producer.compare(value,changed,path)


def test_replay_allows_platform_roundoff():
    producer.compare(170.,170.+1e-7,'root/mass_GeV')
    producer.compare(1.18e-5,1.18e-5+1e-15,'root/GF_local_zero_muon_GeVm2')


def test_independent_top_rejects_wrong_pole_mass(receipt):
    changed=deepcopy(receipt);changed['forward']['rows'][1]['top_method0']['mass_GeV']+=.1
    with pytest.raises(ValueError):top.verify(changed)


def test_independent_top_rejects_width_sign(receipt):
    changed=deepcopy(receipt);changed['forward']['rows'][1]['top_method0']['width_GeV']*=-1
    with pytest.raises(ValueError):top.verify(changed)


def test_independent_tadpole_rejects_tree_minimum_at_one_loop(receipt):
    changed=deepcopy(receipt);m=changed['forward']['m2_minimum_controls'];m['1']=m['0']
    with pytest.raises(ValueError):top.verify(changed)


def test_feynman_integral_absorptive_part():
    # For two massless legs, every interior parameter is timelike.
    assert top.B(0,0,10,2).imag==pytest.approx(math.pi,abs=1e-14)
    assert top.B(20,30,1,2).imag==0


def test_no_empirical_loader_or_alt_evaluator_called():
    source=producer.C_SOURCE.read_text()
    for forbidden in ('SMDR_Read_','SMDR_Fit','SMDR_Eval_yt','SMDR_Eval_lambda','_1loopALT','SMDR_Eval_GFermi_DGG'):
        assert forbidden not in source
    assert 'source_ew_poison.h' in source
    assert '0.00000051862L/(sqrtl(2)*SMDR_v*SMDR_v)' in source
    assert 'SMDR_RGeval_SM(q,2)' in source


@pytest.mark.parametrize('row_index,channel,linux_width', [
    (5, 'top_method1', 1.3079294776870154),
    (8, 'top_method0', 1.1822835492905803),
    (8, 'top_QCD4_EW2', 1.3287015496076546),
])
def test_two_loop_top_width_platform_replay(receipt, row_index, channel, linux_width):
    # Independent x86_64 Linux/Clang execution of the same pinned SMDR source.
    expected = receipt['forward']['rows'][row_index]
    actual = deepcopy(expected)
    actual[channel]['width_GeV'] = linux_width
    producer.compare(expected, actual, f'root/forward/rows/{row_index}')


@pytest.mark.parametrize('order,channel,field,expected,delta', [
    (2, 'top_method1', 'width_GeV', 1.3, 2e-6),
    (1, 'top_method1', 'width_GeV', 1.3, 6e-7),
    (0, 'top_method1', 'width_GeV', 0.0, 6e-7),
    (2, 'W', 'width_GeV', 2.1, 6e-7),
    (2, 'top_method1', 'mass_GeV', 170.0, 4e-6),
    (2, 'working', 'gY', .35, 1e-7),
])
def test_top_width_budget_does_not_relax_other_controls(order, channel, field, expected, delta):
    old = {'matching_order': order, channel: {field: expected}}
    new = deepcopy(old)
    new[channel][field] += delta
    with pytest.raises(ValueError, match='numerical mismatch'):
        producer.compare(old, new, 'root/forward/rows/0')
