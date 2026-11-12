from fractions import Fraction as F
from .calculation import certificate

def test_same_source_binding_and_complete_spectrum():
 d=certificate();assert F(d['binding_lower'])>1
 assert len(d['complete_triplet_spectrum']['all_eigenvalue_intervals'])==27
 assert len(d['all_other_fragment_thresholds'])==9

def test_complete_exchanged_partner_exclusion():
 d=certificate();assert len(set(d['exact_spin_weak_exchange']['full_sector_permutation']))==11700
 for row in d['complete_spin_zero_isotriplet_partners']:
  assert len(row['complete_spectrum_intervals'])==27
  assert F(row['ground_above_same_charge_two_baryon_threshold_lower'])>F(2,5)

def test_complete_original_records_and_held_out_failure():
 d=certificate();r=d['complete_original_records']
 assert abs(sum(r['complete_success_probabilities'])+sum(r['complete_failure_probabilities'])-1)<3e-12
 assert not d['physical_deuteron_binding_reproduced']
