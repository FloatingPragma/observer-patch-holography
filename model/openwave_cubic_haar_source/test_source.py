from fractions import Fraction as F
from .source import run,exact_checks
from .cube import spin_census,cg_entry,basis


def test_complete_gauss_census_and_independent_haar_checks():
    r=exact_checks()
    assert r['complete_dimension_through_energy_six']==32
    assert r['full_successor_commutators_checked']==86*15
    assert min(e for e,s in spin_census(26)if e>24)==26


def test_full_source_gaps_covariances_and_complete_instrument():
    r=run()
    for row in r['records']:
        assert F(row['full_source_gap_interval_exact'][0])>0
        assert row['complete_trial_output_dimension']>row['trial_dimension']
        assert len(row['all_fifteen_connected_covariances'])==15
        assert len(row['complete_joint_Wilson_instrument']['all_sixty_four_outcomes'])==64
    assert all(F(c['connected_covariance_interval_exact'][0])>0
               for c in r['records'][-1]['all_fifteen_connected_covariances']if not c['opposite'])


def test_native_coefficients_retain_spin_and_orientation():
    rows=basis();vacuum=rows[0][1]
    square=next(s for e,s,loops in rows if e==3)
    nonzero=[cg_entry(square,vacuum,f)for f in range(6)]
    assert sorted(abs(x)for x in nonzero)==[0,0,0,0,0,F(1,2)]
    assert all(cg_entry(square,vacuum,f)==cg_entry(vacuum,square,f)for f in range(6))
