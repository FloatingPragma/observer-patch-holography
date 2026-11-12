from fractions import Fraction as F
from .source import spectral_certificate,exact_record_witness
from .calculation import certificate

def test_exact_full_spectrum_band_and_complete_records():
 s=spectral_certificate();rows=[[F(v) for v in r] for r in s['all_eigenvalue_intervals_rational']]
 assert len(rows)==120 and max(v[1] for v in rows[:4])<rows[4][0]
 for q in (0,1):
  a=exact_record_witness(q);p=list(map(F,a['complete_unconditional_native_occupation_probabilities']))
  assert len(p)==2024 and min(p)>=0
  assert sum(p)+F(a['exact_filter_failure_probability'])==1
  assert F(a['trace_distance_to_exact_charge_spin_ground_line_upper'])<F(1,10**7)

def test_signed_intrinsic_charge_and_beta_obstruction():
 a=exact_record_witness(0)['intrinsic_charge_profile']
 assert a['1/3']['ground_interval'][0]>0 and a['2/3']['ground_interval'][1]<0
 c=certificate();assert not c['exact_beta_threshold_obstruction']['n_to_p_e_antinu_ground_channel_energetically_allowed']
 assert 0<float(c['neutron_like_minus_proton_like_electron_units'])<1
