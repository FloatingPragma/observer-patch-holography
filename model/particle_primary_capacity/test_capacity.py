from copy import deepcopy
from decimal import Decimal

import pytest
import sympy as s

from model.particle_common_contract import contract
from . import capacity


def test_native_top_wedge_witness_rejects_an_unjustified_charge_shift():
    result = capacity.singlet_rows()
    assert result['orthogonal_local_singlet_witnesses'] == 16
    shifted = list(capacity.native.native_currents())
    shifted[-1] += s.I*s.eye(15)
    with pytest.raises(ValueError, match='traceless'):
        capacity.singlet_rows(shifted)


def test_capacity_is_log_dimension_and_not_preparation_entropy_or_memory():
    row = capacity.member_capacity(contract.member(0))
    assert row['ambient_capacity_exact']['exact_log_terms'] == [
        {'base': 2, 'coefficient': 480}, {'base': 119, 'coefficient': 24}]
    assert Decimal('447.4096') < Decimal(row['ambient_capacity_exact']['approximate_nats_for_display']) < Decimal('447.4097')
    assert row['physical_capacity_lower']['exact_log_terms'] == [{'base': 2, 'coefficient': 32}]
    assert not row['exact_physical_Gauss_dimension_computed']
    assert not row['prepared_state_entropy_computed']
    assert not row['apparatus_and_retained_record_memory_included']


def test_higher_cutoff_is_a_bound_not_an_exact_physical_capacity():
    member = contract.member(1)
    row = capacity.member_capacity(member)
    assert row['ambient_capacity_exact'] is None
    assert row['Hamiltonian_factor_from_inverse_spacing'] == 2
    assert row['physical_coordinate_extent'] == '4'
    assert row['ambient_capacity_upper']['exact_log_terms'] == [
        {'base': 2, 'coefficient': 60*512}, {'base': 451, 'coefficient': 6*512}]
    bad = deepcopy(member); bad['source']['CAR_modes_per_cell'] = 90
    with pytest.raises(ValueError, match='differs'):
        capacity.member_capacity(bad)


def test_display_logarithm_cannot_accept_inexact_dimension_inputs():
    for terms in ([(2, True)], [(2, 1.5)], [(1, 3)], [(2, -1)]):
        with pytest.raises(ValueError):
            capacity.log_expression(terms)
