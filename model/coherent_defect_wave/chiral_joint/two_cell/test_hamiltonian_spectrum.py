"""Complete-physical-domain energy controls and bad-basis rejection."""
import numpy as np
import pytest
import sympy as s
from .hamiltonian import column,structural_checks,index,from_index,DIMENSION
from .spectrum import physical_block,numerical_block,exact_eigenbasis_error,sector_certificate


def test_uniform_source_contains_both_quantum_creation_terms():
    out=column((0,0,0,0))
    assert len(out)==8
    assert sum(key[2]!=0 for key in out)==4
    assert sum(key[3]!=0 for key in out)==4
    for j in (0,1,149,150,DIMENSION-1):assert index(from_index(j))==j
    with pytest.raises(ValueError):column((4096,0,0,0))


def test_full_elementary_source_and_primary_quantum_vertex_agree():
    result=structural_checks()
    assert result['hopping_CAR_link_columns_checked']==24576
    assert result['local_CAR_scalar_columns_checked']==320
    assert result['identical_to_primary_full_320_state_local_Higgs_operator']
    assert result['exact_local_quantum_Gauss_commutators']==12


def test_complete_Gauss_domain_ground_is_not_in_empty_vacuum_sector():
    result=sector_certificate(0)
    assert result['full_physical_dimension']==1053
    assert result['ground_is_certifiably_simple']
    assert result['ground_window_candidate_labels']==[{'L':2,'twoSz':0,'eigen_index':0}]
    low,high=map(s.Rational,result['global_ground_interval'])
    assert s.Rational(-456,100)<low<high<s.Rational(-455,100)
    assert s.Rational(result['ground_window_to_rest_gap_lower_bound'])>s.Rational(39,100)
    assert sum(len(b['energy_intervals']) for b in result['complete_conserved_blocks'])==1053


def test_both_charged_domains_retain_all_sectors_and_dark_global_ground():
    added=sector_certificate(-1);removed=sector_certificate(1)
    assert added['full_physical_dimension']==718
    assert removed['full_physical_dimension']==783
    assert {row['L'] for row in added['ground_window_candidate_labels']}=={1}
    assert removed['ground_is_certifiably_simple']
    assert removed['ground_window_candidate_labels']==[{'L':4,'twoSz':0,'eigen_index':0}]
    assert sum(len(b['energy_intervals']) for b in added['complete_conserved_blocks'])==718
    assert sum(len(b['energy_intervals']) for b in removed['complete_conserved_blocks'])==783


def test_complete_basis_certificate_rejects_missing_directions():
    _,_,g,k,_=physical_block(0,0,0)
    values,vectors=numerical_block(0,0,0)
    with pytest.raises(AssertionError):exact_eigenbasis_error(k,g,values,np.zeros_like(vectors))
    error,gram=exact_eigenbasis_error(k,g,values,vectors)
    assert error<s.Rational(1,10**7) and gram<1
