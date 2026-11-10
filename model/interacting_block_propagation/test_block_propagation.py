"""Operator, finite-band, preparation and sharp-instability controls."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from .model import adjacency,bond_identity,load_parents,momentum_check


def result():
    from .model import calculate
    return calculate()


def test_exact_bond_and_full_complement():
    out=bond_identity()
    assert out['higher_sector_retained']
    assert out['total_excited_block_number_conserved']
    for branch in out['P_triplet_x_plus_P_triplet_y_plus_or_minus_exchange_positive'].values():
        assert all(F(e)>=0 for e in branch)


def test_frozen_native_parents_and_positive_line():
    parent,_=load_parents()
    line=parent['collective_spin_response']['same_source_interacting_records']['actual_interacting_energy_line']
    assert F(line['first_line_probability_conditioned_on_ground_preparation'][0])>0
    assert result()['native_block_operator_replayed']['actual_native_current_commutators']==36


def test_uniform_gap_and_isolated_whole_band():
    bounds=result()['source_bounds']
    d=F(bounds['block_gap_interval'][0])
    assert F(bounds['volume_uniform_gap_lower'])==3*d/4
    assert F(bounds['isolated_one_excitation_band_separation_lower'])>0
    assert F(bounds['all_other_excitation_energies_lower'])>F(bounds['one_excitation_band_upper'])


def test_every_finite_momentum_and_record():
    value=momentum_check(3)
    assert value['Fourier_eigenvector_max_error']<1e-12
    assert abs(value['position_probability_sum']-1)<1e-12
    assert min(value['complete_position_probabilities'])>=0


def test_large_hopping_really_destabilizes_vacuum():
    # The k=0 single excitation has energy Delta-zJ. The stability
    # requirement cannot be removed by increasing total volume.
    a=adjacency(3)
    constant=np.ones(27)/np.sqrt(27)
    h=np.eye(27)-a/5
    assert np.max(abs(h@constant+constant/5))<1e-14


def test_preparation_error_scales_with_volume():
    prep=result()['preparation_scaling']
    tolerance=F(prep['target_global_trace_distance'])
    r0=F(prep['initial_excited_odds_upper'])
    c=F(prep['geometric_odds_contraction'])
    for row in prep['rows']:
        n=row['binary_filter_depth'];volume=row['blocks']
        assert volume*r0*c**n<=tolerance**2
        if n:
            assert volume*r0*c**(n-1)>tolerance**2
        assert row['all_failure_branches_retained']


def test_spacing_limit_is_identified_correctly():
    out=result()['continuum_controls']
    assert out['limit_generator']=='sum_i k_i^2'
    assert not out['Lorentz_covariance_derived']
    assert out['maintaining_J_proportional_a^-2_violates_fixed_block_gap_control']


def test_receipt_exact_bounds_cannot_be_mutated():
    from copy import deepcopy
    from .certificate import replay_differences
    actual=result()
    broken=deepcopy(actual)
    broken['source_bounds']['volume_uniform_gap_lower']='1'
    assert replay_differences(broken,actual)==['/source_bounds/volume_uniform_gap_lower']


def test_repository_custody_rejects_bad_paths_and_hashes():
    import pytest
    from .certificate import validate_pins
    for pins in ({'../elsewhere':'0'},{'/tmp/elsewhere':'0'},
                 {'model/interacting_block_propagation/model.py':'0'}):
        with pytest.raises(ValueError):validate_pins(pins)
