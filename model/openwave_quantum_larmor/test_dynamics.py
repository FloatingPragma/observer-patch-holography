from .dynamics import branch_certificate,time_bound,records,prior

def test_exact_angular_sector_enclosure_resolves_finite_field_split():
 plus=branch_certificate('1/1000',1);minus=branch_certificate('1/1000',-1)
 assert plus['Jz_exact']=='1/2' and minus['Jz_exact']=='-1/2'
 assert plus['same_source_branch_energy_interval'][0]>minus['same_source_branch_energy_interval'][1]
 assert all(x['same_source_branch_energy_interval'][1]-x['same_source_branch_energy_interval'][0]<2e-12 for x in(plus,minus))

def test_bound_improves_at_fixed_precession_time():
 a=time_bound('1/1000',100000);b=time_bound('1/10000',1000000)
 assert b['full_source_reference_trace_distance_upper']<a['full_source_reference_trace_distance_upper']/9

def test_complete_field_and_failure_records_obey_uniform_bound():
 for row in records():
  assert abs(sum(x['probability']for x in row['complete_success_joint_Jx_photon_number_outputs'])+row['failure_probability']-1)<1e-11
  assert row['numerical_full_pure_state_trace_distance']<=row['bound']['full_source_reference_trace_distance_upper']
  assert len(row['full_joint_failure_state_coefficients'])==len(row['full_joint_success_state_coefficients'])==200
