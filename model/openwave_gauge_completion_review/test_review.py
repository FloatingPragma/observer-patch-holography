from copy import deepcopy
import json
import pytest
from .review import ROOT, car_preparation, pi_certificate, complete_work_account

def test_complete_native_instrument_retains_rejected_unprojected_pump():
    r = car_preparation()
    assert r['exact_complete_instrument_Gram']
    assert r['unprojected_pump_squared_gain_exact'] == '4/3'

def test_owner_pi_range_reduction_constant_is_proved_by_machin_bounds():
    assert pi_certificate()['independent_rational_enclosure_verifies_constant']

def test_all_work_records_retain_bare_field_energy_and_detector_charge():
    receipt = json.loads((ROOT/'model/openwave_transverse_completion/results.json').read_text())
    assert complete_work_account(receipt)['complete_executions_checked'] == 30
    wrong = deepcopy(receipt)
    wrong['complete_native_executions'][0]['preparation_mean_source_work_electron_units_diagnostic'] += 1/16
    with pytest.raises(AssertionError):
        complete_work_account(wrong)
