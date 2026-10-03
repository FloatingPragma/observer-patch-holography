"""Independent transport, charged-current and pole-convention regressions."""
import copy
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_source_ew_vev_matching as verifier


@pytest.fixture(scope='module')
def receipt():
    return verifier.strict_load(verifier.OUT)


def test_independent_ew_reconstruction(receipt, monkeypatch):
    import builtins
    original = builtins.__import__

    def checked(name, *args, **kwargs):
        if name.split('.')[-1] in {'source_ew_vev_matching', 'conditional_quark_mass_replay'}:
            raise AssertionError('producer/target import in independent verifier')
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', checked)
    result = verifier.verify(receipt)
    assert result['matching_rows_checked'] == 9
    assert result['independent_GF_orders'] == [0, 1]
    assert result['full_multiloop_replay_required']
    assert not result['independent_two_loop_matching_claimed']


@pytest.mark.parametrize('row,key,delta', [
    (0, 'v', .01), (3, 'v', 1.), (6, 'v', 2.),
    (1, 'Q', 1.), (4, 'gY', .0001), (4, 'g2', .0001),
    (4, 'g3', .0001), (7, 'yt', .0001), (7, 'lambda', .0001),
    (4, 'mt_MSbar', .01),
])
def test_running_parameter_mutations_fail_without_using_payload_digest(receipt, row, key, delta):
    forward = copy.deepcopy(receipt['forward'])
    forward['rows'][row]['working'][key] += delta
    with pytest.raises(ValueError):
        verifier.verify_gf_and_running(receipt['inputs'], forward)


def test_constant_vev_substitution_rejected(receipt):
    forward = copy.deepcopy(receipt['forward'])
    for row in forward['rows']:
        row['working']['v'] = receipt['inputs']['v_Q0_GeV']
    with pytest.raises(ValueError):
        verifier.verify_gf_and_running(receipt['inputs'], forward)


@pytest.mark.parametrize('row', [0, 1, 3, 4, 6, 7])
def test_coherent_GF_shift_rejected_by_independent_formula(receipt, row):
    forward = copy.deepcopy(receipt['forward'])
    for key in ('GF_raw_GeVm2', 'GF_local_zero_muon_GeVm2'):
        forward['rows'][row][key] += 1e-10
    with pytest.raises(ValueError):
        verifier.verify_gf_and_running(receipt['inputs'], forward)


@pytest.mark.parametrize('row', [0, 4, 8])
def test_empirical_muon_correction_cannot_survive(receipt, row):
    forward = copy.deepcopy(receipt['forward'])
    forward['rows'][row]['GF_local_zero_muon_GeVm2'] = forward['rows'][row]['GF_raw_GeVm2']
    with pytest.raises(ValueError):
        verifier.verify_gf_and_running(receipt['inputs'], forward)


@pytest.mark.parametrize('particle', ['W', 'Z'])
def test_complex_pole_is_not_BW_mass(receipt, particle):
    forward = copy.deepcopy(receipt['forward'])
    forward['rows'][2][particle]['BW_mass_GeV'] = forward['rows'][2][particle]['mass_GeV']
    with pytest.raises(ValueError):
        verifier.verify_gf_and_running(receipt['inputs'], forward)


def test_two_loop_tampering_requires_full_replay_not_one_loop_label(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['forward']['rows'][2]['top_QCD4_EW2']['mass_GeV'] += .1
    with pytest.raises(ValueError, match='multiloop payload'):
        verifier.verify(mutant)


@pytest.mark.parametrize('key', ['source_derived_vev_attachment_claimed',
                                 'complete_flavor_model_claimed',
                                 'theory_uncertainty_quantified'])
def test_physical_promotion_rejected(receipt, key):
    mutant = copy.deepcopy(receipt)
    mutant['scope'][key] = True
    with pytest.raises(ValueError):
        verifier.verify(mutant)


def test_no_scale_control_may_be_dropped(receipt):
    forward = copy.deepcopy(receipt['forward'])
    forward['rows'].pop()
    with pytest.raises(ValueError):
        verifier.verify_gf_and_running(receipt['inputs'], forward)


@pytest.mark.parametrize('text', ['{"x":1,"x":2}', '{"x":NaN}'])
def test_ambiguous_json_rejected(tmp_path, text):
    path = tmp_path/'bad.json'
    path.write_text(text)
    with pytest.raises(ValueError):
        verifier.strict_load(path)
