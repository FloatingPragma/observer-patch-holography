from fractions import Fraction as F
import numpy as np
from .source import effects,exact_readout_control,interval_joint,infer_gap

def test_native_readout_retains_every_raw_slot_and_complement():
 d=exact_readout_control();assert len(d['fine_records'])==25
 for j in range(2):assert sum(F(x['input_probabilities_exact'][j])for x in d['fine_records'])==1
 assert any(x['input_probabilities_exact']==['0','0']for x in d['fine_records'])
 for spectator in range(3):
  for axis in ('X','Y'):
   P=effects(spectator,axis);assert np.linalg.norm(sum(P)-np.eye(36))<2e-15
   for a in P:assert np.linalg.norm(a@a-a)<2e-15

def test_complete_two_time_ideal_normalization_and_record_inverse():
 from model.characteristic_constraints.bounds import trig_interval
 for gap in (F(1,256),-F(1,1024),-F(1,3072)):
  for axis in ('X','Y'):
   b=interval_joint(gap,F(64),axis,F(1,10000));assert sum(x[0]for row in b for x in row)<=1<=sum(x[1]for row in b for x in row)
  l,u=infer_gap(trig_interval(gap*64),F(64));assert l<=gap<=u

def test_replay_has_two_blocks_no_reset_and_no_mass_promotion():
 import json
 from pathlib import Path
 p=Path(__file__).with_name('results.json')
 if not p.exists():return
 d=json.loads(p.read_text());assert len(d['repeated_records'])==36 and len(d['twirled_controls'])==3
 assert d['resource_contract']['two_blocks_without_repreparation']
 assert not d['scope']['literal_OpenWave_row3_complete']
 assert not d['scope']['three_spectator_values_are_three_lepton_species']
