import sympy as s
from .transfer import spinors,local_action,edge_moments,certificate
from model.defect_mass_join.source import car

def test_native_spinors_and_all_link_columns():
 states=spinors();matrices,supports=edge_moments()
 assert len(matrices)==3 and supports==[[3,15,18,18]]*6
 assert all(not local_action(car.gauge_current(a),state)for a in range(12)for state in states)
 wrong={(1<<10,0):s.Integer(1)}
 assert any(local_action(car.gauge_current(a),wrong)for a in range(12))

def test_boundary_control_and_full_probability_bound():
 r=certificate();assert r['boundary_control_third_moment']=='0'
 assert r['full_source_H3_receiver_amplitude']=='-3*sqrt(2)*I/64'
 lo,hi=r['original_receiver_probability_decimal_interval'];assert 0<lo<hi<1
 assert r['same_side_two_low_energy_band_or_ground_transferred']is False
