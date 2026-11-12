from .current import norm_controls,certificate

def test_global_current_control_is_fibre_bound_not_global_B_norm():
 d=norm_controls();assert d['probe_Hilbert_Schmidt_squared_exact']=='5'
 assert len(d['all_probe_and_field_matrix_controls'])==4

def test_infinite_field_current_and_long_time_complete_bound():
 d=certificate();lo,hi=d['rest_energy_normalized_infinite_g_interval'];assert 1.99827<lo<hi<1.99828 and hi-lo<1e-6
 assert d['global_differentiated_remainder_upper']<1e-9 and d['exact_Jz_target']=='1/2'
 assert d['infinite_ground_line_trace_distance_upper']<1e-9
 rows=d['complete_finite_time_Larmor_bounds'];assert rows[-1]['complete_source_arbitrary_reference_trace_distance_upper']<1.2e-5
 assert rows[-1]['complete_source_arbitrary_reference_trace_distance_upper']<rows[2]['complete_source_arbitrary_reference_trace_distance_upper']/90
