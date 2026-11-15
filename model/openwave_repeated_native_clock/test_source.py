from fractions import Fraction as F
import numpy as np
import pytest
from .algebra import run,exchange
from .controller import source,run as controller_run,propagation_constants,atan_bounds,phase_frequency_enclosure

def test_exact_two_cycle_complete_work_reset_and_depletion():
    for c,s in ((F(3,5),F(4,5)),(F(5,13),F(12,13))):
        r=run(1,3,c,s)
        assert sum(map(F,r['complete_64_occupation_probabilities_exact']))==1
        assert F(r['gap_switching_work_over_mass_exact'])==1-c**4/2
        assert r['source_excited_after_ideal_reset_exact']=='0'
        assert F(r['preparation_reference_final_occupation_exact'])==c**4/2

def test_complete_twirled_negative_control():
    rows=[run(a,b,twirl=True) for a,b in ((0,0),(1,2),(2,3))]
    assert all(r['complete_64_occupation_probabilities_exact']==rows[0]['complete_64_occupation_probabilities_exact'] for r in rows)
    assert all(r['second_prepared_source_A_expectation_exact']=='0' for r in rows)

def test_complete_static_matrix_domain_and_resources():
    H,pieces,exchanges,v=source(F(3477),3)
    assert H.shape==(32768,32768) and max(np.diff(H.indptr))<=4
    assert np.max(np.asarray(abs(H).sum(axis=1)))<40
    assert (H-H.T).nnz==0 and abs(np.vdot(v,v)-1)<1e-14
    assert set(pieces)=={'native_source','preparation','first_readout','battery','second_readout'}
    with pytest.raises(ValueError):source(F(4097))
    with pytest.raises(ValueError):source(F(1),4)

def test_rational_propagation_budget_and_phase_bounds():
    d,e0,e,r=propagation_constants()
    assert 0<d< F(1,10**10) and 0<e<F(1,10**7)
    import mpmath as mp
    mp.mp.dps=200
    for x in (F(1,10),F(4,5),F(3,2),F(-4,5)):
        lo,hi=atan_bounds(x);actual=mp.atan(mp.mpf(x.numerator)/x.denominator)
        assert mp.mpf(lo.numerator)/lo.denominator<=actual<=mp.mpf(hi.numerator)/hi.denominator
    assert phase_frequency_enclosure(0j,1j,1.,1.,F(1)) is None

def test_actual_repeated_controller_retains_failures_and_certified_bias():
    r=controller_run(F(1),2)
    assert len(r['complete_terminal_probabilities'])==512
    assert all(len(x)==64 for x in r['complete_terminal_probabilities'])
    assert r['reset_source_excitation_conditioned_upper']<F(1,1000)
    for flight in r['two_free_flight_certificates']:
        lo,hi=map(F,flight['certified_actual_phase_frequency_interval_exact'])
        assert float(lo)<=flight['phase_difference_frequency_diagnostic']<=float(hi)
        assert flight['visibility_lower']>0
        assert flight['continuous_phase_branch_without_winding_certified']
        assert flight['all_clock_outside_corridor_outcomes_retained']
    assert F(r['two_free_flight_certificates'][1]['certified_actual_phase_frequency_interval_exact'][1])<1
    assert abs(r['closed_total_energy_change_diagnostic'])<r['closed_total_energy_change_absolute_bound']
