from fractions import Fraction as F
import numpy as np
from .neutral import construction,complete_vertex,exact_spin_contraction,free_source_check,flight,ward

def test_complete_massive_source_and_current():
 for epsilon in(2**-32,.25):
  d=construction(epsilon);assert free_source_check(d)['full_CAR_modes']==288
  labels,W=complete_vertex(d);assert W.shape==(144,2)
  assert np.linalg.norm(W.conj().T@W-np.eye(2))<3e-13
  assert len(exact_spin_contraction())==32

def test_complete_flight_failure_and_work():
 d=construction();labels,W=complete_vertex(d);v=flight(d,labels,W)
 assert 0<v['production_selected_probability']<=1+2e-13
 for row in v['complete_records']:
  assert abs(row['original_record_probability_sum']-1)<2e-13
  assert 0<row['production_filter_failure_probability']<1e-12
  assert row['original_mean_detector_work']<0

def test_native_gauge_and_unfiltered_bound():
 assert ward()['actual_native_weak_Higgs_intertwiner_identities']==3
 G=F(1,64);t=F(1,10**8);r=(G+G*G*41952**2)*t*t/2
 assert 0<r<G*t and (G*t-r)**2>0
