"""Independent native-path, tensor, spectral and retained-boundary controls."""
from collections import Counter
from fractions import Fraction as Q
import math

import numpy as np
import pytest

from . import certificate as c


def numeric_entropy(matrix):
    eigenvalues = np.linalg.eigvalsh(np.asarray(matrix,dtype=float))
    assert min(eigenvalues) > 0
    return -sum(x*math.log(x) for x in eigenvalues)


def numeric_log_pairing(matrix,reference):
    values,vectors = np.linalg.eigh(np.asarray(reference,dtype=float))
    logarithm = (vectors*np.log(values))@vectors.T
    return float(np.trace(np.asarray(matrix,dtype=float)@logarithm))


def test_all_native_two_attempt_paths_independently_reproduce_reference():
    # Call the endpoint function directly, then enumerate all 237600 paths.
    branches = [tuple(c.source.INDEX[c.source.native.native(state,edge,coin)[0]]
                      for edge in c.source.quantum_mixing.seams() for coin in (0,1))
                for state in c.STATES]
    actual = [Counter(last for middle in row for last in branches[middle]) for row in branches]
    C,Z,sigma,floor = c.source_data()
    assert all(Q(actual[i][j],Z) == sigma[i,j] for i in range(c.DIM) for j in range(c.DIM))
    assert sum(actual[i][i] for i in range(c.DIM)) == Z == 171720
    assert floor == Q(40,4293)


def test_tensor_partial_trace_independently_matches_completion():
    sigma = c.source_data()[2]
    expected = np.full(sigma.shape,Q(),dtype=object)
    for q in (0,1,2):
        local,exterior,lookup = c.layout(c.SEAM)[q]
        indices = [lookup[a,b] for a in local for b in exterior]
        tensor = sigma[np.ix_(indices,indices)].reshape(len(local),len(exterior),len(local),len(exterior))
        reduced = np.trace(tensor,axis1=1,axis2=3)
        assert np.array_equal(reduced,c.restrict(sigma,c.SEAM)[q])
        identity = np.array([[Q(int(i==j),len(exterior)) for j in range(len(exterior))]
                             for i in range(len(exterior))],dtype=object)
        expected[np.ix_(indices,indices)] = np.kron(reduced,identity)
    assert np.array_equal(expected,c.complete(sigma,c.SEAM))


def test_expectation_tower_on_arbitrary_correlated_matrix():
    matrix = np.array([[Q((i+1)*(j+1),1000000) for j in range(c.DIM)]
                       for i in range(c.DIM)],dtype=object)
    A = (c.SEAM[0],)
    B = c.SEAM
    assert np.array_equal(c.complete(c.complete(matrix,B),A),c.complete(matrix,A))
    assert np.array_equal(c.complete(c.complete(matrix,A),B),c.complete(matrix,A))
    assert c.trace(c.complete(matrix,A)) == c.trace(matrix)


def test_all_thirty_source_regions_have_exact_spectra_and_nonzero_reference_defect():
    sigma = c.source_data()[2]
    reference = c.restrict(sigma,c.SEAM)
    expected = sorted(v for B in reference.values() for v,P in c.spectrum(B))
    for region in c.source.quantum_mixing.seams():
        blocks = c.restrict(sigma,region)
        assert sorted(v for B in blocks.values() for v,P in c.spectrum(B)) == expected
        assert c.hs_squared(sigma-c.complete(sigma,region)) > 0
        charge = np.array([sum(state[p] for p in region) for state in c.STATES],dtype=object)
        assert c.hs_squared(sigma*charge-charge[:,None]*sigma) > 0


def test_central_shannon_is_in_bulk_and_completion_spectrum_is_independent():
    blocks = c.restrict(c.source_data()[2],c.SEAM)
    bulk,edge = c.bulk_entropy(blocks),c.edge_entropy(blocks,c.SEAM)
    assert bulk == c.bulk_sector_split(blocks)
    assert c.completion_entropy(blocks,c.SEAM) == c.combine((1,bulk),(1,edge))
    shannon = c.combine(*[(-c.trace(B),c.log_rational(c.trace(B))) for B in blocks.values()])
    assert c.expression_interval(shannon)[0] > 0
    assert bulk != c.combine((1,bulk),(-1,shannon))
    completed = c.complete(c.source_data()[2],c.SEAM)
    lo,hi = c.expression_interval(c.completion_entropy(blocks,c.SEAM))
    assert float(lo)-1e-12 <= numeric_entropy(completed) <= float(hi)+1e-12


def test_native_midpoint_filter_retains_all_return_and_nonreturn_weights():
    values = tuple(state[c.SEAM[0]]-state[c.SEAM[1]] for state in c.STATES)
    weights = c.complete_return_weights(values,Q(1,2))
    assert sum(map(Q,weights.values())) == 1
    assert Q(weights['accepted_return']) == Q(477,1320)
    assert weights['accepted_return'] == weights['rejected_return']


def test_occupation_fixed_variation_changes_entropy_without_static_cut_response():
    values,nullity = c.occupation_fixed_filter()
    tangent = c.filter_tangent(values)
    assert nullity == 54
    assert c.trace(tangent) == 0
    assert all(sum(state[p]*tangent[i,i] for i,state in enumerate(c.STATES)) == 0 for p in range(12))
    derivative = c.modular_pairing(c.restrict(tangent,c.SEAM),c.restrict(c.source_data()[2],c.SEAM))
    assert c.expression_interval(derivative)[0] > 0
    assert sum((a in c.SEAM)!=(b in c.SEAM) for a,b in c.source.quantum_mixing.seams()) == 8
    assert c.trace(c.restrict(tangent,c.SEAM)[2]) != 0


def test_noncommuting_rational_spectrum_response_has_positive_exact_remainder():
    reference = c.restrict(c.source_data()[2],c.SEAM)
    values = tuple(state[c.SEAM[0]]-state[c.SEAM[1]] for state in c.STATES)
    delta = c.restrict(c.filter_tangent(values),c.SEAM)
    theta = 3*reference[1][0,1]/(4*abs(delta[1][0,0]))
    assert 0 < theta < Q(1,2)
    assert c.hs_squared(reference[1]@delta[1]-delta[1]@reference[1]) > 0
    receipt = c.response_receipt(values,theta)
    assert receipt['bulk_entropy_derivative']['prime_log_coefficients'] == {}
    assert Q(receipt['relative_entropy']['outward_interval'][0]) > 0
    assert Q(receipt['bulk_entropy_change']['outward_interval'][1]) < 0


@pytest.mark.parametrize('theta',[Q(-1,2),Q(-1,10),Q(1,10),Q(1,2)])
def test_general_noncommuting_remainder_against_independent_eigendecomposition(theta):
    _,_,sigma,floor = c.source_data()
    # A different one-configuration filter, whose perturbed spectra need not be rational.
    values = tuple(Q(int(i==0))-sigma[0,0] for i in range(c.DIM))
    tangent = c.filter_tangent(values)
    reference,delta = c.restrict(sigma,c.SEAM),c.restrict(tangent,c.SEAM)
    rho = c.restrict(c.filtered_state(values,theta),c.SEAM)
    D = sum(-numeric_entropy(rho[q])-numeric_log_pairing(rho[q],reference[q]) for q in reference)
    entropy_delta = sum(numeric_entropy(rho[q])-numeric_entropy(reference[q]) for q in reference)
    derivative = sum(-numeric_log_pairing(delta[q],reference[q]) for q in reference)
    bound = theta**2*sum(c.hs_squared(B) for B in delta.values())/(2*(1-abs(theta))*floor)
    assert D > 0 and D <= float(bound)
    assert abs(entropy_delta-float(theta)*derivative+D) < 1e-12


def test_exact_log_intervals_bound_independent_scalar_values():
    for value in (Q(1,100),Q(9,20),Q(2),Q(113,7155),Q(1000)):
        lo,hi = c.log_interval(value)
        assert float(lo)-1e-14 <= math.log(float(value)) <= float(hi)+1e-14
        assert hi-lo < Q(1,10**35)


@pytest.mark.parametrize('region',[(),(1,0),(0,0),(12,),(True,),tuple(range(12))])
def test_invalid_region_rejected(region):
    with pytest.raises(ValueError,match='native port region'):
        c.layout(region)


def test_invalid_filters_and_nonpositive_state_parameters_rejected():
    with pytest.raises(ValueError,match='normalization'):
        c.filter_tangent((Q(1),)*c.DIM)
    with pytest.raises(ValueError,match='unit interval'):
        c.filter_tangent((Q(2),)*c.DIM)
    with pytest.raises(ValueError,match='exact rational'):
        c.filter_tangent((0.0,)*c.DIM)
    with pytest.raises(ValueError,match='filter parameter'):
        c.filtered_state((Q(),)*c.DIM,Q(1))


def test_receipt_retains_einstein_and_preparation_boundaries():
    receipt = c.strict_json.load(c.HERE/'results.json')
    assert receipt['scope']['count_derived_common_reference']
    assert receipt['scope']['quantitative_modular_first_law']
    assert not receipt['scope']['native_quantum_state_prepared']
    assert receipt['scope']['response_is_postselected_history_family']
    assert not receipt['scope']['causal_quantum_state_intervention_executed']
    assert not receipt['scope']['fixed_proper_volume_stationarity']
    assert not receipt['scope']['exterior_multiplicity_identified_with_area']
    assert not receipt['scope']['Einstein_premises_G1_to_G6']
