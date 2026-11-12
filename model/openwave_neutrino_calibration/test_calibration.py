from fractions import Fraction as F
import numpy as np
from .calibration import calibration,catalogue,rigid_obstruction
from .neutral import construction,complete_vertex,free_source_check,complete_flight,vacuum_kernel

def test_scale_obstruction_and_exact_calibration():
 o=rigid_obstruction();assert not o['one_scale_fits_both_mass_splittings']
 assert F(o['normal_discrepancy_factor_lower'])>7
 d=calibration();assert d['source_recoil_q_squared_exact']>F(49,100)**2;assert max(d['source_mass_squared_exact'])<F(5,10**10)**2
 sq=d['mass_squared_eV2_exact'];assert sq[1]-sq[0]==F('0.0000749');assert sq[2]-sq[0]==F('0.002513')
 for row,x in zip(catalogue()['rows'],sq):
  lo,hi=map(F,row['mass_eV_interval']);assert lo*lo<=x<=hi*hi
 assert np.linalg.norm(d['R'].conj().T@d['R']-np.eye(3))<1e-14

def test_actual_native_source_current_and_complete_flight():
 d=construction();assert free_source_check(d)['full_CAR_modes']==288
 labels,W=complete_vertex(d);r=complete_flight(d,labels,W)
 assert F(r['uniform_exact_complete_record_TV_bound'])<F(41,10**10)
 for row in r['complete_records']:
  assert abs(row['original_probability_total']-1)<2e-13
  assert row['production_filter_failure']>0
  assert row['original_mean_detector_work']<0
  assert row['complete_probability_TV_difference_diagnostic']<1e-12

def test_oscillation_controls_and_CP_conjugation():
 d=construction();assert np.linalg.norm(vacuum_kernel(d,0)-np.eye(3))<1e-14
 a=vacuum_kernel(d,np.pi);b=vacuum_kernel(d,np.pi,True);assert a[0,1]-b[0,1]>.01
 conjugate={**d,'R':d['R'].conj()};assert np.max(abs(vacuum_kernel(conjugate,np.pi)-b))<1e-14
 idle={**d,'R':np.eye(3)};assert np.max(abs(vacuum_kernel(idle,123)-np.eye(3)))<1e-14


def test_native_antiparent_current_and_CP_records():
 from .antiparticle import exact_antiparticle_spin_contraction,complete_antiparticle_vertex
 assert len(exact_antiparticle_spin_contraction())==32
 d=construction();labels,W=complete_antiparticle_vertex(d);r=complete_flight(d,labels,W,(0,1))
 for row in r['complete_records']:
  assert abs(row['original_probability_total']-1)<2e-13
  assert row['production_filter_failure']>0
  assert row['complete_probability_TV_difference_diagnostic']<1e-12
 assert abs(r['complete_records'][2]['antineutrino_e_detection_probability']-.0346502689391265)<1e-14
