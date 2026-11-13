import sympy as s
from .flux import haar,class_moments,records,monodromy_controls,full_CAR_pulse,native_holes,complete_native_preparation

def test_exact_Haar_moments_and_no_singlet_postselection():
    assert haar(((0,0),(1,1),(2,2)),())==s.Rational(1,6)
    assert haar(((0,0),),((0,0),))==s.Rational(1,3)
    for r in range(3):
        m=class_moments(r)
        assert m['postselected_colour_singlet_probability_loop_branch_exact']=='5/27'
        assert len(m['all_nine_U_absolute_square_means_exact'])==9

def test_complete_spatial_flux_records_and_combined_turn():
    r=records(-1,1,2,s.pi/2)
    assert s.simplify(sum(s.sympify(x['probability_exact'])for x in r['all_spatial_flux_colour_pointer_records']))==1
    assert {x['quark_colour']for x in r['all_spatial_flux_colour_pointer_records']}=={0,1,2}
    assert all(x['combined_spatial_flux_and_centre_record_invariance']for x in monodromy_controls())

def test_literal_native_holes_and_full_CAR_transport():
    assert [r['native_charge_exact']for r in native_holes()]==['2/3','-1/3']
    p=full_CAR_pulse();assert len(p['all64_energy_values_diagnostic'])==64
    assert p['complete_Fock_quarter_turn_error']<2e-14

def test_preparation_includes_all_reference_failure_branches():
    p=complete_native_preparation()
    assert p['complete_active_Fock_dimension']==64
    assert p['exact_complete_Gram_sum']=='I_64'
    assert p['stated_filled_source_success_probability_exact']=='1'
    assert p['unprojected_pump_countercontrol_squared_norm_lower_exact']=='4/3'
