from fractions import Fraction as F
import math
from .source import certificate,sinc_interval,charge_form_factor

def test_physical_complete_records_and_exact_mass_current_transport():
    value=certificate();assert F(value['mass_ordering_GeV_exact'])>F('0.00051')
    for row in value['native_states']:
        records=row['complete_native_records']
        assert sum(F(x['success'])for x in records)==F(1,120)
        assert sum(F(x['failure'])for x in records)==F(119,120)
        assert all(F(x['failure'])>=0 for x in records)
        assert row['orientation_averaged_form_factors'][0]['exact_ground_form_factor_interval']==[row['charge_exact']]*2

def test_radius_ratio_eliminates_any_common_length():
    value=certificate()['shared_length_comparison']
    lower,upper=map(F,value['minus_neutron_over_proton_radius_squared_interval'])
    observed=tuple(map(F,value['observed_one_sigma_box_interval']))
    assert observed[1]<lower<upper and not value['physical_radius_target_used_to_fit_source_or_length']
    rows=certificate()['native_states']
    n,p=(F(r['mean_square_charge_radius_graph_exact'])for r in rows)
    for length in (F(1,100),F(1),F(100)):
        assert -n*length**2/(p*length**2)==-n/p

def test_taylor_remainder_and_zero_charge_retained():
    for x in map(F,('0','.1','1','3','-3')):
        lo,hi=sinc_interval(x,order=4)
        value=1 if x==0 else math.sin(float(x))/float(x)
        assert float(lo)-2e-16<=value<=float(hi)+2e-16
    assert charge_form_factor({F(1,3):F(1),F(2,3):-F(1)},F(0),F(1),F(1,100))==(0,0)
