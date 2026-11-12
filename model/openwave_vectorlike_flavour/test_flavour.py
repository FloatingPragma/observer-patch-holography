import numpy as np
from .flavour import branches,ward_certificate,floor_certificate,certificate

def test_native_ward_and_complete_positive_floor():
 assert ward_certificate()['vectorlike_symmetric_cubic_anomaly_cancellations']==364
 assert len(floor_certificate()['rows'])==6

def test_full_singular_values_ordering_and_no_mixing_control():
 for m in (1e-4,.02,.3,2.):
  for y in (-2/3,1/3):
   lo,hi=branches(m,y);direct=np.linalg.svd([[m,0],[y,1]],compute_uv=False)
   assert np.allclose([float(hi),float(lo)],direct,rtol=1e-13)
   assert abs(float(lo*hi)-m)<1e-13
  assert branches(m,-2/3)[0]<branches(m,1/3)[0]
  assert abs(float(branches(m,0,g=0)[0])-min(m,1))<1e-13

def test_all_channels_outputs_and_same_electron_unit():
 r=certificate();assert len(r['quark_branches'])==6 and r['extra_branches_retained']==6
 unit=float(r['calibration']['source_energy_unit_GeV'])
 for x in r['quark_branches']:
  assert abs(sum(x['complete_channel_records']['rest_complete_original_leg_probabilities'])-1)<1e-12
  full=x['complete_channel_records']['full_CAR_addition_records']
  assert abs(sum(full['complete_unconditional_native_occupation_probabilities'])+full['creation_failure_probability']-1)<1e-11
  assert abs(float(x['conditional_light_GeV'])-unit*float(x['light_source_gap']))<1e-13
