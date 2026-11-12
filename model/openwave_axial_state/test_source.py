from fractions import Fraction as F
import numpy as np
from .source import parameters,sector
from .certificate import ground_certificate
from .fragments import fragment_controls

def test_exact_source_complements_and_native_colour_edge():
 p=parameters();assert p['mbar']>p['delta']/2>0
 for kind,n,k in(('baryon',2024,120),('meson',576,64)):
  d=sector(kind);assert d['H'].shape==(n,n)and d['h'].shape==(k,k)
  assert (d['parts']['Cg48']@d['Binteger']).nnz==0
  assert np.linalg.norm(d['H']@d['B']-d['B']@d['h'])<2e-8

def test_certified_native_state_and_complete_original_records():
 d=ground_certificate();lo,hi=map(F,d['native_ground_axial_interval_rational']);assert lo<F('1.2753')<hi and hi-lo<F('1e-8')
 r=d['complete_trial_original_records'];assert sum(map(F,r['success_native_occupation_probabilities']))+sum(map(F,r['failure_native_occupation_probabilities']))==1
 assert d['colocated_ground_probability_interval'][1]<.036

def test_retained_failed_copy_control_and_zero_same_edge_claim():
 d=fragment_controls();assert all(x['threshold_minus_baryon']<0 for x in d['copy_comparisons']if x['edge_colour_coefficient']==0)
 assert all(x['threshold_minus_baryon']>2 for x in d['copy_comparisons']if x['edge_colour_coefficient']==1)
 assert len(d['complete_fragment_spectra'])==10
