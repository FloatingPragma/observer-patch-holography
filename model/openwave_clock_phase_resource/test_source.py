from fractions import Fraction as F
import numpy as np
import sympy as s
import pytest
from .source import (exact_source,complete_word_records,battery_basis_control,
                     relational_records,PROB)

def test_every_word_matches_dense_native_Kraus_and_work():
    N,A,ks=exact_source();dense=[np.asarray(k,dtype=complex)for k in ks]
    rows=complete_word_records(3)
    for case in rows:
        for row in case['all_battery_words']:
            v=np.eye(4,dtype=complex)[:,case['initial_native_occupation']]
            for event in row['battery_word']:v=dense[event]@v
            assert abs(np.vdot(v,v).real-float(F(row['probability_exact'])))<1e-14
            rho=np.outer(v,v.conj());h=np.asarray(N,dtype=complex)
            assert np.array_equal(h@rho,rho@h)
    with pytest.raises(ValueError):complete_word_records(8)

def test_energy_mixing_readout_is_the_phase_resource():
    r=battery_basis_control()
    assert sum(s.sympify(x['probability_exact'])for x in r['records'])==1
    assert all(not x['energy_invariant_branch']for x in r['records'][:2])
    assert r['records'][2]['probability_exact']=='0'
    assert F(r['readout_mean_work_over_mass_exact'])==PROB-F(1,2)

def test_relational_exchange_does_not_supply_free_mass_frequency():
    r=relational_records();assert len(r['records'])==4
    for row in r['records']:
        assert len(row['all_sixteen_native_occupation_probabilities_exact'])==16
        assert sum(map(s.sympify,row['all_sixteen_native_occupation_probabilities_exact']))==1
        assert row['free_total_energy_commutator_exact']=='0'
    assert r['records'][1]['relative_coherence_exact']!='0'
    assert r['free_native_frequency_of_relative_coherence_exact']=='0'
