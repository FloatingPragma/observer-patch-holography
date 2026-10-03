"""Comparison conventions, complete controls and target-isolation guards."""
from copy import deepcopy
from decimal import Decimal
from pathlib import Path
import sys
import pytest
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import compare_source_ew_vev_matching as consumer

@pytest.fixture(scope='module')
def packet():return consumer.strict_load(consumer.PACKET)
@pytest.fixture(scope='module')
def targets():return consumer.read_targets()
@pytest.fixture(scope='module')
def comparison():return consumer.build()

def test_full_verified_receipt_rebuild(comparison):
 assert consumer.serialized(comparison)==consumer.serialized(consumer.strict_load(consumer.OUTPUT))
 assert len(comparison['retained_forward_output']['rows'])==9
 assert len(comparison['comparisons'])==57
 assert [(r['scale_multiplier'],r['matching_order']) for r in comparison['retained_forward_output']['rows']]==[(s,o) for s in (1,2,4) for o in (0,1,2)]
 assert all(r['scale_multiplier']==1 and r['matching_order']==2 for r in comparison['primary'])
 assert next(r for r in comparison['primary'] if r['quantity']=='top')['method']=='top_QCD4_EW2'
 assert comparison['joint_significance'] is None and comparison['theory_uncertainty'] is None

def test_bw_masses_and_zero_muon_coefficient_used(packet,comparison):
 source=packet['forward']['rows'][2];primary={r['quantity']:r for r in comparison['primary']}
 for q in ('W','Z'):
  assert Decimal(primary[q]['computed'])==Decimal(str(source[q]['BW_mass_GeV']))
  assert Decimal(primary[q]['computed'])!=Decimal(str(source[q]['mass_GeV']))
 assert primary['GF']['unit']=='GeV^-2'
 assert Decimal(primary['GF']['computed'])==Decimal(str(source['GF_local_zero_muon_GeVm2']))
 assert Decimal(primary['GF']['computed'])!=Decimal(str(source['GF_raw_GeVm2']))

def test_descriptive_residuals_have_correct_units(comparison):
 for r in comparison['comparisons']:
  computed=Decimal(r['computed']);reference=Decimal(r['reference']['central']);delta=computed-reference
  assert Decimal(r['signed_difference'])==delta
  assert Decimal(r['absolute_difference'])==abs(delta)
  assert abs(Decimal(r['relative_percent'])-100*delta/reference)<Decimal('1e-25')
 assert all(r['signed_difference']!='0' for r in comparison['primary'])

@pytest.mark.parametrize('mutate',[
 lambda p:p['forward']['rows'].pop(),
 lambda p:p['forward']['rows'].__setitem__(0,p['forward']['rows'][2]),
 lambda p:p['forward']['rows'][2]['W'].__setitem__('BW_mass_GeV',p['forward']['rows'][2]['W']['mass_GeV']),
 lambda p:p['forward']['rows'][2]['Z'].__setitem__('BW_mass_GeV',p['forward']['rows'][2]['Z']['mass_GeV']),
 lambda p:p['forward']['rows'][2].__setitem__('GF_local_zero_muon_GeVm2',1e6*p['forward']['rows'][2]['GF_local_zero_muon_GeVm2']),
 lambda p:p['forward']['rows'][2].pop('top_QCD4_EW2'),
 lambda p:p['forward']['rows'][2]['top_QCD4_EW2'].__setitem__('mass_GeV',172.4),
 lambda p:p['specification']['calculation'].__setitem__('primary_scale_multiplier',4),
])
def test_modified_or_dropped_forward_controls_rejected(packet,targets,mutate):
 changed=deepcopy(packet);mutate(changed)
 with pytest.raises(ValueError):consumer.compare(changed,targets)

@pytest.mark.parametrize('field,value',[
 ('unit','GeV^2'),('unit','MeV^-2'),('central','0.00001166378500000000001'),
 ('comparison_chart','operational_muon_decay_prediction'),
])
def test_GF_target_units_or_meaning_cannot_change(packet,targets,field,value):
 changed=deepcopy(targets);next(r for r in changed['targets'] if r['quantity']=='GF')[field]=value
 with pytest.raises(ValueError):consumer.compare(packet,changed)

@pytest.mark.parametrize('quantity,chart',[
 ('W','complex_pole_real_mass'),('Z','complex_pole_real_mass'),('top','Monte_Carlo_mass'),
])
def test_incompatible_mass_charts_rejected(packet,targets,quantity,chart):
 changed=deepcopy(targets);next(r for r in changed['targets'] if r['quantity']==quantity)['comparison_chart']=chart
 with pytest.raises(ValueError):consumer.compare(packet,changed)

@pytest.mark.parametrize('field',[
 'retuned_parameters','scale_selected_by_agreement','joint_significance_claimed','confidence_interval_claimed',
 'theory_uncertainty_supplied','complete_flavor_model_claimed','local_GF_is_operational_muon_lifetime',
 'blind_prediction_claimed','successful_postdiction_promoted'])
def test_comparison_promotion_rejected(packet,targets,field):
 changed=deepcopy(targets);changed['interpretation'][field]=True
 with pytest.raises(ValueError):consumer.compare(packet,changed)

def test_unfavorable_comparison_cannot_be_removed(comparison):
 changed=deepcopy(comparison);changed['comparisons'].pop()
 with pytest.raises(ValueError):consumer.verify_payload(changed)

def test_absolute_residual_tampering_rejected(comparison):
 changed=deepcopy(comparison);changed['primary'][0]['absolute_difference']='0'
 with pytest.raises(ValueError):consumer.verify_payload(changed)

def test_target_file_is_not_a_forward_input(monkeypatch):
 import source_ew_vev_matching as forward
 original=Path.read_text
 def guarded(path,*a,**kw):
  if path in (consumer.TARGETS,consumer.OUTPUT,consumer.PRIOR_TOP):raise AssertionError('empirical comparison input entered forward producer')
  return original(path,*a,**kw)
 monkeypatch.setattr(Path,'read_text',guarded)
 inputs=forward.build_inputs()
 assert inputs['Q0_GeV']>0

@pytest.mark.parametrize('text',['{"a":1,"a":2}','{"a":NaN}','{"a":Infinity}'])
def test_invalid_json_rejected(tmp_path,text):
 p=tmp_path/'invalid.json';p.write_text(text)
 with pytest.raises(ValueError):consumer.strict_load(p)
