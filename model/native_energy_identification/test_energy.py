"""Independent native-transfer and spectral-detector controls."""
from fractions import Fraction as F
from itertools import product
from math import comb, log

import numpy as np
import pytest

from . import energy as e, packet


@pytest.mark.parametrize('order',[1,2,5,16,32,64])
def test_complete_positive_binomial_law_including_endpoints(order):
    detector = e.plan(630,6,2,order)
    for x in (F(0),F(1,2),F(7,10),F(1)):
        actual = e.histogram(detector['power_coefficients'],tuple(x**j for j in range(order+1)))
        expected = tuple(comb(order,k)*x**k*(1-x)**(order-k) for k in range(order+1))
        assert actual == expected
        assert min(actual)>=0 and sum(actual)==1
    assert sum(detector['bin_amplifications'])==3**order-1


@pytest.mark.parametrize('order',[3,16])
def test_histogram_from_independently_enumerated_native_population(order):
    family, _, _, probe = packet.ground.local_case()
    source = e.exchange.sector.source
    states = source.sector_masks(12,2)
    index = {state:i for i,state in enumerate(states)}
    transition = np.zeros((len(states),len(states)))
    for mask in states:
        for out,weight in source.population_column(family,mask).items():
            transition[index[out],index[mask]] = float(weight)
    assert np.max(abs(transition-transition.T))<1e-14
    eigenvalues, vectors = np.linalg.eigh(transition)
    h = np.array([probe.get(tuple(i for i in range(12) if mask>>i&1),0) for mask in states])/2
    spectral_weights = abs(vectors.T@h)**2
    expected = np.array([np.sum(spectral_weights*comb(order,k)*eigenvalues**k*(1-eigenvalues)**(order-k))
                         for k in range(order+1)])
    exact = e.exact_case(family,probe,2,order)
    assert np.max(abs(np.array(list(map(float,exact['spectral_bin_weights'])))-expected))<2e-13
    assert sum(exact['spectral_bin_weights'])==1
    assert exact['spectral_bin_weights'][0]>0


@pytest.mark.parametrize('order',[4,16,32])
def test_energy_resolution_and_cdf_against_exact_conditional_laws(order):
    margin = F(1,4)
    detector = e.plan(630,6,2,order)
    budget = e.resolution_budget(630,6,2,order,margin)
    for numerator in range(10,21):
        x = F(numerator,20)
        law = [comb(order,k)*x**k*(1-x)**(order-k) for k in range(order+1)]
        bad = sum(law[k] for k in range(order+1) if abs(F(k,order)-x)>=margin)
        assert bad<=budget['bad_label_probability_upper']
        assert law[0]<=budget['zero_count_overflow_probability_upper']
        for k in range(1,order+1):
            if abs(F(k,order)-x)<margin:
                error = float(detector['energy_scale'])*abs(log(float(x)/float(F(k,order))))
                assert error<=float(budget['counted_energy_error_on_good_labels'])+1e-12
        for threshold in (F(1,2),F(3,4),F(1)):
            observed = sum(law[k] for k in range(order+1) if F(k,order)<=threshold)
            assert int(x<=threshold-margin)-budget['bad_label_probability_upper']<=observed
            assert observed<=int(x<=threshold+margin)+budget['bad_label_probability_upper']


@pytest.mark.parametrize('value',[F(1),F(15,16),F(1,2),F(1,16)])
def test_rational_log_enclosure(value):
    low,high = e.negative_log_interval(value,64)
    expected = -log(float(value))
    assert float(low)<=expected+1e-14
    assert float(high)>=expected-1e-14
    assert low<=high
    assert e.energy_labels(30,1,16)[0] is None


def test_no_clipping_of_inconsistent_noisy_moments():
    detector = e.plan(30,5,2,2)
    result = e.histogram(detector['power_coefficients'],(F(1),F(1,4),F(1,2)))
    assert result==(F(1),F(-1,2),F(1,2))
    assert sum(result)==1


def test_vector_error_for_every_extreme_moment_error():
    detector = e.plan(30,5,2,4)
    ideal = tuple(F(3,4)**j for j in range(5))
    expected = e.histogram(detector['power_coefficients'],ideal)
    error = F(1,1000)
    for signs in product((-1,1),repeat=4):
        perturbed = (F(1),)+tuple(ideal[j+1]+signs[j]*error for j in range(4))
        result = e.histogram(detector['power_coefficients'],perturbed)
        assert sum(result)==1
        assert sum(abs(a-b) for a,b in zip(result,expected,strict=True))/2<=detector['total_l1_amplification']*error/2


def test_sampling_budget_is_sufficient_and_uses_exact_m0():
    detector = e.plan(30,5,2,16)
    variance = F(2,33); error = F(1,100); failure = F(1,100)
    budget = e.estimation_budget(detector,variance,F(0),error,failure)
    threshold = 2*error/detector['total_l1_amplification']
    union_bound = F(16,budget['independent_repetitions_per_setting'])/(variance*threshold)**2
    assert union_bound<=failure
    assert budget['estimated_lag_settings']==16
    assert budget['total_independent_repetitions']==16*budget['independent_repetitions_per_setting']


def test_native_preparation_and_lag_program_share_one_actual_clock():
    family,first,second,probe = packet.ground.local_case()
    detector = e.plan(30,5,2,16)
    programs = packet.local_programs(first,second,detector['lag_attempts'])
    assert len(programs)==16
    assert programs[-1]['resources']['attempts_upper']==3871
    assert programs[-1]['resources']['AST_nodes']==3856
    assert len({row['program_digest'] for row in programs})==16


def test_capacity_independent_resolution_does_not_remove_attempt_growth():
    values = [e.plan(630*4**level,6,2,32) for level in (0,2,20,200)]
    assert all(value['energy_scale']<=96 for value in values)
    assert values[-1]['lag_attempts'][-1]>2**32
    assert all(value['power_coefficients']==values[0]['power_coefficients'] for value in values)


@pytest.mark.parametrize('args',[(True,5,2,16),(30,5,2,0),(19,5,2,2),(30,5,2,65),(30.0,5,2,16)])
def test_invalid_spectral_domain_rejected(args):
    with pytest.raises(ValueError):e.plan(*args)


@pytest.mark.parametrize('value',[F(0),F(-1),F(2),0.5,True])
def test_invalid_log_domain_rejected(value):
    with pytest.raises(ValueError):e.negative_log_interval(value)


def test_invalid_moment_table_rejected():
    coefficients = e.plan(30,5,2,2)['power_coefficients']
    for values in ((F(0),F(1),F(1)),(F(1),F(1)),(1,0.5,0.25)):
        with pytest.raises(ValueError):e.histogram(coefficients,values)
    bad = ((2,-2,1),(0,2,-2),(0,0,1))
    with pytest.raises(ValueError):e.histogram(bad,(F(1),F(1),F(1)))


def test_receipt_is_exact_fresh_native_replay():
    assert packet.verify()
