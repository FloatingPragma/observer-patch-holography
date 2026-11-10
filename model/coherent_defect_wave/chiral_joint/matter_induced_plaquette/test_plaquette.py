"""Physical full-current loops, actual signs, complete tails and failures."""
import sympy as s
from .model import (BLOCKS,geometry_certificate,physical_checks,leading_histories,
                    loop_state,inner,vacuum_coordinate)
from .bounds import feshbach_bound,complete_records


def test_contractible_geometry_and_full_domain_current_identities():
    assert geometry_certificate()['Pauli_loop_product']=='-I2'
    for setting in BLOCKS:
        row=physical_checks(setting)
        assert row['full_119_power_4_loop_operator_Ward_equations']==96
        assert row['loop_local_Gauss_equations']==row['preparation_local_Gauss_equations']==96


def test_actual_four_and_five_factor_histories_include_colored_trace():
    for setting,d in (('singlet_ec',1),('colored_up',3)):
        row=leading_histories(setting)
        assert row['full_H_moments_0_through_5']==['0','0','0','0',str(-s.Rational(1,256*d)),str(-s.Rational(5,128*d))]
        assert len(row['nonzero_fourth_order_histories'])==1


def test_full_carrier_Feshbach_interval_is_nonzero():
    for setting in BLOCKS:
        row=feshbach_bound(setting);lower,upper=map(s.Rational,row['exact_loop_kernel_real_part_enclosure'])
        assert 0<lower<upper and row['P_plus_Q_equals_full_carrier_identity']


def test_complete_original_clock_records_have_positive_success_and_failure():
    row=complete_records();lo,hi=map(s.Rational,row['complete_four_outcome_probability_enclosure'])
    assert lo<=1<=hi
    assert len(row['settings'])==2
    for setting in row['settings']:
        for name in ('success','failure'):
            low,high=map(s.Rational,setting['original_setting_'+name+'_probability_interval'])
            assert 0<low<high<s.Rational(1,2)


def test_wrong_colored_normalization_and_orientation_are_distinguished():
    loop=loop_state('colored_up')
    assert inner({k:v/3 for k,v in loop.items()},{k:v/3 for k,v in loop.items()})==s.Rational(1,9)
    assert vacuum_coordinate(0,0,False)[0]!=vacuum_coordinate(0,0,True)[0]
