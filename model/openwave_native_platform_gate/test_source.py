from .headless import run,independent_column
from fractions import Fraction as F

def test_complete_native_channel_and_altered_rule():
    x=run();assert len(x['all_endpoint_matrix_units'])==16
    assert x['maximum_load_conservation_error']==0 and x['altered_rule_negative_control_error']==1

def test_swapped_population_and_no_false_coherent_swap():
    assert independent_column(1,1)=={(1,1):F(1,2),(2,2):F(1,2)}
    assert independent_column(1,2)=={}
