from fractions import Fraction as F
from .source import gauss_flux,electric,bound,A,B,UNIT

def test_full_integer_flux_dictionary_and_native_fractional_charge_hop():
    for q in(F(-2,3),F(-1,3),F(1,3),F(2,3),F(1)):
        before,r=gauss_flux([q,0,0,0],2)
        after,rr=gauss_flux([0,0,0,q],2)
        assert r==rr and after[-1]==before[-1]
        assert [after[i]-before[i]for i in range(3)]==[-3*q]*3

def test_actual_field_kernel_is_not_a_pairwise_inverse_distance_counterterm():
    assert electric([1,0,0])==0
    assert electric([0,0,1])==9*A
    assert electric([0,1,1])==F(45)*A/2
    # The retained reference fixes a physical boundary relation.
    assert electric([0,0,1],1)!=electric([0,0,1],0)

def test_complete_dynamic_ground_bounds_exclude_both_held_out_masses():
    helion=bound(3,2);alpha=bound(4,2)
    assert F(helion['full_dynamic_ground_GeV_interval'][1])<F('2.80839161112')
    assert F(alpha['full_dynamic_ground_GeV_interval'][0])>F('3.7273794118')
    assert helion['complete_bare_flux_zero_magnetic_dimension']==24
    assert alpha['complete_bare_flux_zero_magnetic_dimension']==96
    assert F(helion['uniform_charge_quartet_electric_mean_source_exact'])==F(3,4)
    assert F(alpha['uniform_charge_quartet_electric_mean_source_exact'])==F(39,32)
    assert all(row['all_integer_cycle_fluxes_retained']for row in(helion,alpha))
