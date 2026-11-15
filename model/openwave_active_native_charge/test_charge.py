"""Controls for the active voltage, complete alphabet and physical dictionary."""
from fractions import Fraction as F
import pytest
from . import source as s


def test_native_complete_slots_and_charge_transport():
    for law in ('C','E','mixture'):
        rows=s.native.raw_local_branches(law)
        assert len(rows)==25
        for pair in s.PAIRS:
            for seam in ((0,5),(0,11),(5,11)):
                total=F(0)
                for r in rows:
                    y=s.target(pair,seam,r)
                    if y is not None:
                        total+=r.weight
                        for m in range(6):
                            z=s.lifted_column(pair,m,seam,r)
                            znext=s.lifted_column(pair,(m+1)%6,seam,r)
                            assert znext==(z[0],(z[1]+1)%6)
                assert total==1


def test_closed_histories_generate_six_sectors_and_voltage_control():
    odd=s.loop((0,11),[((0,5),(5,11)),((0,11),(0,5)),((5,11),(0,11))])
    third=s.loop((0,1),[((0,5),(1,5)),((5,11),(1,11)),((0,11),(0,1))])
    assert odd['voltage_mod6']==3 and third['voltage_mod6']==2
    assert odd['probability_exact']==third['probability_exact']=='1/8'
    assert s.voltage((0,1),(1,5),(0,5))==5
    assert s.voltage((0,1),(1,5),(0,5),extra=False)==3
    assert s.voltage((0,6),(5,6),(0,5))==0


def test_modular_half_is_not_additive_real_electric_charge():
    p=s.charge_boundaries()
    assert p['equal_class_split_solutions_mod6']['3']==[]
    assert p['equal_class_split_solutions_mod6']['2']==[1,4]
    assert p['representative_nonadditivity_control']['sum_assigned_charge_exact']=='2'
    assert p['representative_nonadditivity_control']['representative_of_sum_exact']=='0'
    assert not p['physical_charge_criterion_satisfied']
    assert not p['physical_quark_criterion_satisfied']


def test_memory_dephasing_and_incompatible_motion_are_rejected():
    p=s.localization_and_memory_error()
    assert sum(map(F,p['basis_memory_dephasing_character_distribution_exact']))==1
    assert p['source_reduced_trace_distance_exact']=='0'
    with pytest.raises(ValueError):s.voltage((0,1),(5,6),(0,5))
    with pytest.raises(ValueError):s.voltage((0,1),(1,5),(0,5),extra=1)
