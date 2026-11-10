"""Full-domain certification and retained-apparatus regression controls."""
import pytest
import sympy as s
from .basis import metadata,physical_basis
from .spectrum import resolved_block,sector_certificate
from .symmetry import certificate_data
from .readout import charged_spectroscopy,probe_matrix


def test_exact_physical_dimensions_include_every_occupation_sector():
    assert metadata(0)['physical_dimension']==3424
    assert metadata(-1)['physical_dimension']==2456
    assert metadata(1)['physical_dimension']==2456
    b,g=physical_basis(0,2,0)
    assert b.cols==g.rows==502
    assert resolved_block(0,2,0,1)[2].rows==251
    assert resolved_block(0,2,0,-1)[2].rows==251
    with pytest.raises(ValueError):resolved_block(-1,1,1,1)


def test_unique_neutral_ground_uses_both_translation_sectors():
    data=sector_certificate(0)
    assert data['ground_certifiably_simple']
    assert data['ground_window']==[{'L':2,'twoSz':0,'translation_sign':1,'eigen_index':0}]
    assert sum(len(b['energy_intervals']) for b in data['blocks'])==3424
    low,high=map(s.Rational,data['global_ground_interval'])
    assert s.Rational(-434,100)<low<high<s.Rational(-433,100)
    assert s.Rational(data['ground_window_to_rest_gap_lower_bound'])>s.Rational(25,100)


def test_complete_conjugate_charged_sectors_have_the_same_ground_window():
    added=sector_certificate(-1);removed=sector_certificate(1)
    assert added['global_ground_interval']==removed['global_ground_interval']
    assert added['ground_window_rank']==removed['ground_window_rank']==2
    assert {r['L'] for r in added['ground_window']}=={1}
    assert {r['L'] for r in removed['ground_window']}=={3}
    for result in (added,removed):
        assert sum(len(b['energy_intervals']) for b in result['blocks'])==2456


def test_actual_conjugation_includes_quantum_links_and_boundary_charges():
    data=certificate_data()
    assert data['same_boundary_periodic_H_antiunitary_symmetry']
    assert data['actual_full_H_fixture_columns_checked']>1000
    assert data['actual_local_Gauss_checks']>=1000
    assert data['same_site_addition_removal_spectral_measures_match']
    with pytest.raises(ValueError):probe_matrix('naked_mass')


def test_same_ground_positive_charged_residue_and_complete_protocol():
    result=charged_spectroscopy()
    assert result['original_ground_preparation_success']=='1/3424'
    assert result['original_ground_preparation_failure']=='3423/3424'
    assert result['conditional_outcome_count_including_zero_phase_records']==2626
    assert result['original_attempt_outcome_count_including_preparation_failure']==2627
    lo,hi=map(s.Rational,result['complete_probability_enclosure'])
    assert lo<=1<=hi and hi-lo<s.Rational(1,10**5)
    first=[]
    for branch in result['charged_success_branches']:
        assert branch['exact_total_branch_probability']=='1/4'
        assert len(branch['complete_energy_windows'])==328
        record=branch['complete_energy_windows'][0]
        low,high=map(s.Rational,record['original_swap_probability_interval'])
        assert s.Rational(22,100)<low<=high<s.Rational(23,100)
        elo,ehi=map(s.Rational,record['energy_transfer_interval'])
        assert s.Rational(97,100)<elo<ehi<s.Rational(98,100)
        assert record['joint_source_reference_phase_splitting_factors']==['0','1/2','0','1/2']
        first.append(record['energy_transfer_interval'])
    assert first[0]==first[1]
    assert all(r['probability']=='1/4' and not r['energy_readout_performed'] for r in result['retained_failures'])
    assert not result['probabilities_renormalized']
    assert not result['source_only_particle_momentum_or_rest_mass_identified']
    assert not result['GeV_prediction']
