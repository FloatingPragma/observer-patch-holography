from fractions import Fraction as F
import numpy as np
from . import averaging as a

def test_exact_full_current_orbit_and_detuning():
    v=a.exact_orbit();assert v['all_resonant_current_columns']==[['173/2','1'],['1','85']]
    assert v['exact_invariant_orbit_dimension']==2
    assert v['first_peak_receiver_probability']=='16/25'

def test_every_native_frequency_and_full_bound():
    v=a.frequency_certificate();assert v['native_quartic_words']==31464
    assert v['transformed_ordered_words']==55575 and v['frequency_labels']==299
    b=a.full_source_bound();assert F(b['uniform_full_source_unitary_error_upper'])<F(1,2500000)
    assert b['original_G_1_over_64_averaging_certified']is False
    assert F(b['all_other_source_outcomes_probability_upper'])<F(1,10**12)

def test_bogoliubov_orthogonality_and_full_energy():
    d=a.source()
    for i,terms in d['numeric_mapping'].items():
        assert abs(sum(z*z for _,_,z in terms)-1)<1e-14
    rows=a.small_block_records();assert abs(rows[-1]['receiver_probability']-.64)<3e-14
    assert all(abs(r['current_energy_expectation']-86.5)<1e-12 for r in rows)

def test_polynomial_reduction_retains_cancellations():
    p={((2,0,0,0,0,0),0):F(1),((0,2,0,0,0,0),0):F(1),((0,0,0,0,0,0),0):F(-1)}
    assert a.reduce_polynomial(p)=={}
    p[((0,0,0,0,0,0),0)]=F(-2)
    assert a.reduce_polynomial(p)=={((0,0,0,0,0,0),0):F(-1)}
