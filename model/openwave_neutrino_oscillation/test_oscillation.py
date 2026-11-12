import numpy as np
import pytest
from .oscillation import construction,records,certificate

def test_actual_current_frame_and_complete_records():
 d=construction();assert d['ward']==12
 for a in range(3):
  row=records(d['C'],d['masses'],1.,10.,a);p=row['complete_original_probabilities']
  assert len(p)==5 and abs(sum(p)-1)<1e-12 and p[3]>0 and p[4]>0
  assert row['complete_law_TV_from_ideal']<=row['proved_uniform_TV_upper']+1e-12

def test_invalid_or_incomplete_instrument_rejected():
 d=construction()
 for C,p in((2*d['C'],1),(d['C'],0),(d['C'][:2],1)):
  with pytest.raises(ValueError):records(C,d['masses'],p,1,0)

def test_rephasing_and_no_mixing_controls():
 r=certificate();assert r['zero_port_phase_control_has_no_family_appearance']
 assert r['Majorana_column_phase_control_error']<1e-12
 assert r['maximum_executed_neutral_antineutral_probability_contrast']>1e-9
 assert r['imported_PMNS_angles_mass_splittings_or_neutrino_GeV_scale']is False
