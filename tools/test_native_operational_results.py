"""Exact observable roles and rejected physical promotion in the public ledger."""
import copy,json
from pathlib import Path
import pytest
from tools.native_operational_results import NAMES,project,render
from tools.native_physical_results import interval,decimal
ROOT=Path(__file__).resolve().parents[1]
@pytest.fixture(scope='module')
def records():return {n:json.loads((ROOT/'model'/('openwave_'+n)/'results.json').read_text())for n in NAMES}
def test_complete_numeric_source_domains(records):
 p=project(records)
 assert p['clock']['complete_control_cases']==6 and p['fock']['modes']==192
 assert p['gravity']['observed_G_calibration'] and not p['lambda']['energy_is_rest_mass']
 assert len(p['gauge'])==2 and p['new_independent_nature_validations']==0
 text='\n'.join(render(p,interval,decimal))
 assert '2.418046697462' in text and '0.734179279095' in text
 assert 'residual adiabatic coordinate mode' in text and 'not its rest mass' in text
@pytest.mark.parametrize('name,path',[
 ('finite_phase_clock',('scope','literal_openwave_row3_validated')),
 ('lambda_spectral_measure',('on_shell_native_weak_density_proved',)),
 ('nuclear_em_transfer',('new_mass_fit',)),
 ('energy_schur',('nontrivial_continuum_or_Clay_result',)),
 ('gravity_noise_normalization',('scope','independent_Newton_prediction')),
 ('capacity_normalization_transfer',('mode_and_observation_scope','metric_ratio_is_parent_finite_k_matter_growth_ratio'))])
def test_no_promotion_by_receipt_label(records,name,path):
 d=copy.deepcopy(records);x=d[name]
 for key in path[:-1]:x=x[key]
 x[path[-1]]=True
 with pytest.raises(ValueError):project(d)
