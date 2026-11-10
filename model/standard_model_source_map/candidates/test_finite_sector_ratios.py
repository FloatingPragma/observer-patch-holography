"""Check correlated energy-zero bounds and their limited physical scope."""
import copy
from fractions import Fraction as F
from itertools import product

import pytest

from model.standard_model_source_map.candidates import finite_sector_ratios as ratios


def test_common_ground_correlation_is_retained_and_improves_independent_division():
    a,b,g=(F(8),F(10)),(F(5),F(7)),(F(1),F(3))
    lo,hi=ratios.ratio_enclosure(a,b,g)
    assert lo==F(7,6) and hi==F(7,2)
    assert hi < (a[1]-g[0])/(b[0]-g[1])
    for x,y,z in product((a[0],F(9),a[1]),(b[0],F(6),b[1]),(g[0],F(2),g[1])):
        assert lo <= (x-z)/(y-z) <= hi


def test_inverse_same_sector_and_bad_domains():
    a,b,g=(8,10),(5,7),(1,3)
    ab=ratios.ratio_enclosure(a,b,g);ba=ratios.ratio_enclosure(b,a,g)
    assert ab==(1/ba[1],1/ba[0])
    assert ratios.ratio_enclosure(a,a,g,True)==(1,1)
    with pytest.raises(ValueError):ratios.ratio_enclosure(a,b,(1,5))
    with pytest.raises(ValueError):ratios.ratio_enclosure(a,b,g,True)


def test_actual_source_intervals_ratios_and_nonpromotions():
    value=ratios.certificate();r=value['result']
    assert len(r['source_sector_rows'])==3 and len(r['ordered_ratio_rows'])==9
    assert all(row['total_native_Q6']==0 for row in r['source_sector_rows'])
    assert F(r['all_three_family_ratio_upper'])<F('21.16')
    assert F(r['odd_family_ratio_upper'])<F('15.31')
    for row in r['ordered_ratio_rows']:
        lo,hi=map(F,row['finite_excitation_ratio_interval']);assert 0<lo<=hi
        if row['numerator_id']==row['denominator_id']:assert lo==hi==1
    assert r['direct_common_scale_identification_with_charged_lepton_hierarchy_excluded_at_member_zero']
    for key in ('exclusion_applies_to_other_capacities_or_other_source_excitations',
                'exclusion_is_a_no_go_for_OPH','ratios_are_point_eigenvalue_calculations',
                'physical_dispersion_or_particle_species_identified',
                'finite_threshold_intervals_are_physical_rest_masses'):
        assert r[key] is False
    assert r['physical_mass_table_GeV']==[]
    assert not r['observational_comparison']['used_to_construct_source_or_bounds']


def test_changed_clock_rejected():
    a=ratios._receipt(ratios.SOURCE_INPUTS[3])['selected_member_bound']['result']
    b=copy.deepcopy(ratios._receipt(ratios.SOURCE_INPUTS[4])['selected_member_bound']['result'])
    b['clock_id']='different'
    with pytest.raises(ValueError,match='clock'):ratios.calculate(a,b)
