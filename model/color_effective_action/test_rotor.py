"""The distinct implemented rotor family has a uniform exact regulator limit."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import solve

from model.interacting_extension.chiral_lattice.rotor_response import build_flux_ring
from . import rotor


def compressed_inverse(q, cutoff, mu=8, sites=4):
    source = build_flux_ring(q, cutoff, L=sites)
    p = [source['labels'].index((x, 0)) for x in range(sites)]
    columns = np.eye(len(source['H']))[:, p]
    return solve(mu*np.eye(len(source['H']))+source['H'], columns, assume_a='pos')[p]


@pytest.mark.parametrize('q', (-4, -3, 1, 2, 6))
@pytest.mark.parametrize('sites', (3, 4, 5))
def test_every_source_column_and_Gauss_row_is_retained(q, sites):
    receipt = rotor.source_replay(q, 5, sites)
    assert receipt['complete_columns_checked'] == 11*sites
    assert receipt['Gauss_rows_checked'] == 11*sites*sites


@pytest.mark.parametrize('q', (-4, -3, 1, 2, 6))
def test_infinite_series_encloses_large_finite_source_inverse(q):
    matrix, series_error = rotor.local_resolvent(q, degree=12)
    actual = compressed_inverse(q, 24)
    comparison = rotor.cutoff_error(q, 24)
    error = series_error+F(comparison['local_compressed_resolvent_error'])
    assert np.linalg.norm(np.asarray(matrix, float)-actual, 2) < float(error)


def test_exit_and_return_degree_is_required_before_local_cutoff_changes():
    bound = rotor.cutoff_error(1, 1)
    assert bound['minimum_exit_insertions'] == 2
    assert bound['minimum_exit_and_return_insertions'] == 4
    a, _ = rotor.local_resolvent(1, degree=3)
    b, _ = rotor.local_resolvent(1, degree=3, cutoff=1)
    assert a == b
    a, _ = rotor.local_resolvent(1, degree=4)
    b, _ = rotor.local_resolvent(1, degree=4, cutoff=1)
    assert a != b


@pytest.mark.parametrize('q', (1, 6))
def test_local_action_enclosure_matches_the_complete_regulated_inverse(q):
    witness = rotor.action_witness(q)
    actual = np.linalg.inv(compressed_inverse(q, 48))
    cutoff_error = F(rotor.cutoff_error(q, 48)['local_compressed_resolvent_error'])
    floor = 1/(F(8)+F(3*q*q, 2)+3)
    transformed = cutoff_error/(floor*(floor-cutoff_error))
    for i, row in enumerate(witness['local_action_intervals']):
        for j, (lo, hi) in enumerate(row):
            assert float(F(lo)-transformed)-2e-12 <= actual[i, j] <= float(F(hi)+transformed)+2e-12
    assert witness['infinite_flux_limit_constructed_on_this_fixed_graph']
    assert not witness['native60_Higgs_source_replaced_or_identified']
    assert not witness['nonabelian_or_spatial_continuum_identification']


def test_bounds_are_independent_of_spatial_size_and_electric_cutoff_norm():
    assert rotor.cutoff_error(6, 48, sites=4) == rotor.cutoff_error(6, 48, sites=4000)
    for cutoff in (1, 4, 8, 12, 48):
        actual = rotor.cutoff_error(6, cutoff)
        count = (cutoff+6)//6
        assert actual['minimum_exit_insertions'] == count
        assert F(actual['local_compressed_resolvent_error']) == F(2, 5)*F(3, 8)**(2*count)


@pytest.mark.parametrize('kwargs', ({'mu':3}, {'initial_flux':5}, {'cutoff':-1}, {'q':1.0}))
def test_invalid_regulator_domains_fail(kwargs):
    args = {'q':1, 'cutoff':4, **kwargs}
    with pytest.raises((ValueError, TypeError)):
        rotor.cutoff_error(**args)
