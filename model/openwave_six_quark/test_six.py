from fractions import Fraction as F
from .source import reduction,native_Fierz
from .certificate import certificate

def test_actual_native_current_and_complete_exact_projector():
 assert native_Fierz()['exact_actual_native_Fierz_entries']==81
 d=reduction();assert d['full_dimension']==11700 and d['zero_Cartan_dimension']==1800
 assert sum(d['Pnum'].diagonal())==27*d['Pden'] and d['Gram'].det()>0

def test_same_source_binding_and_unbound_architecture_control():
 d=certificate();assert F(d['strict_binding_below_two_baryons_lower'])>F(6,25)
 assert F(d['unbalanced_S1_I0_ground_above_two_baryons_lower'])>F(9,10000)
 assert len(d['complete_balanced_spectrum']['all_eigenvalue_intervals'])==27

def test_complete_records_and_positive_intrinsic_quadrupole():
 d=certificate()['positive_bound_state']
 assert abs(sum(d['complete_success_occupation_probabilities'])+sum(d['complete_failure_occupation_probabilities'])-1)<1e-11
 assert F(d['ground_quadrupole_interval'][0])>0 and F(d['all_ground_vectors_quadrupole_lower_from_binding'])>0
 assert F(d['ground_line_trace_distance_upper'])<F(1,10**7)
 assert len(d['complete_success_occupation_probabilities'])==11700
