from copy import deepcopy
from pathlib import Path
import json,pytest
from tools.coherent_integration_results import NAMES,project,render
@pytest.fixture(scope='module')
def data():
 root=Path(__file__).resolve().parents[1]
 return {n:json.loads((root/f'model/openwave_{n}/results.json').read_text())for n in NAMES}
def test_same_hardware_retains_distinct_programs_and_no_physical_award(data):
 p=project(data);assert len(p['source_packages'])==7 and len(p['named_programs'])==4
 assert not p['joint_spatial_charge_clock_interacting_field_implemented']
 assert p['uncompensated_clock_counterexample']['coherence_multiplier_real_exact']=='85/96'
 assert p['exact_interfaces']['compensated_charge_clock_branch_identity']
 assert not p['new_independent_GeV_rest_masses']and p['official_OpenWave_untested_rows']==31
 assert 'zero independent Nature' in '\n'.join(render(p))
@pytest.mark.parametrize('package,path',[
 ('native_coherent_clock',('scope','literal_OpenWave_row3_complete')),
 ('native_calibrated_energy',('energy_attachment','physical_rest_pole_or_particle_identification')),
 ('native_charge_readout',('matching_clock_recovery_join','uncompensated_clock_theorem_transferred')),
 ('active_native_charge',('scope','physical_electric_charge_identified')),
 ('native_euclidean_transfer',('scope','four_dimensional_field_limit'))])
def test_interface_and_physics_promotions_rejected(data,package,path):
 d=deepcopy(data);v=d[package]
 for k in path[:-1]:v=v[k]
 v[path[-1]]=True
 with pytest.raises(ValueError):project(d)

def test_distinct_hardware_cannot_be_joined_by_program_labels(data):
 d=deepcopy(data);d['native_charge_readout']['shared_construction']={'different':'hardware'}
 with pytest.raises(ValueError,match='different hardware'):project(d)
