from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from ..interacting_vacuum import ground
from . import certificate as c,packet,comparison


def test_exact_square_root_and_decimal_enclosures():
    for x in [F(0),F(2),F(10**100+1,10**90),F(1,10**80)]:
        lo,hi=c.sqrt_bounds(x)
        assert lo*lo<=x<=hi*hi
    assert c.decimal_bounds(F(-1,3),F(1,3),3)==['-0.334','0.334']
    with pytest.raises(ValueError):c.sqrt_bounds(True)


def test_all_visible_lines_are_retained_and_positive(completed_spectrum):
    s=completed_spectrum['spectrum'];rows=s['rows']
    assert len(rows)==228 and s['signed_local_retarded_atoms']==456
    assert all(F(r['local_branch_weight_interval'][0])>0 for r in rows)
    assert min(F(r['local_branch_weight_interval'][0]) for r in rows)>F(4,10**21)
    assert min(r['weight_diagnostic'] for r in rows)<1e-20
    assert s['unresolved_local_weight']=='0'
    assert all(F(a['excitation_energy_interval'][1])<F(b['excitation_energy_interval'][0])
               for a,b in zip(rows,rows[1:]))


def test_complete_matrix_moments_enclose_source_identity(completed_spectrum):
    rows=completed_spectrum['spectrum']['rows'];a=ground.band_data()[0]
    eta=[(-1)**x.bit_count() for x in range(8)]
    for x in range(8):
        for y in range(8):
            bounds=[[F(0),F(0)] for _ in range(3)]
            for row in rows:
                blo,bhi=map(F,row['removal_spatial_residue_intervals'][x][y])
                elo,ehi=map(F,row['excitation_energy_interval'])
                for power,factor in enumerate([1+eta[x]*eta[y],eta[x]*eta[y]-1,1+eta[x]*eta[y]]):
                    vals=[factor*b*e**power for b in (blo,bhi) for e in (elo,ehi)]
                    bounds[power][0]+=min(vals);bounds[power][1]+=max(vals)
            expected=[F(int(x==y)),F(int(a[x,y]),16),F(67*int(x==y),256)]
            assert all(lo<=value<=hi for (lo,hi),value in zip(bounds,expected))


def test_visible_cyclic_rank_and_dark_states(completed_spectrum):
    s=completed_spectrum['spectrum']
    assert s['multiplicity_counts']=={'4':228}
    assert s['reachable_N7_dimension_all_linear_CAR_probes']==1824
    assert s['unobserved_N7_dimension']==9616
    assert s['all_N7_N9_eigenvalues_resolved'] is False


def test_ground_auxiliary_degeneracy_is_not_mislabelled():
    witness=packet.load_witness()
    r=c.eigensystem('ground',witness['ground_columns'],witness['ground_energies'],witness['ground_vectors'])
    assert r['vector_errors'][0] is not None
    assert sum(x is None for x in r['vector_errors'])==4


def test_corrupted_complete_basis_is_rejected():
    w=packet.load_witness();vectors=w['ground_vectors'].copy();vectors[:,0]=0
    with pytest.raises(ArithmeticError):c.eigensystem('ground',w['ground_columns'],w['ground_energies'],vectors)


def test_corrupted_energy_is_rejected():
    w=packet.load_witness();values=w['ground_energies'].copy();values[0]-=1<<52
    with pytest.raises(ArithmeticError):c.eigensystem('ground',w['ground_columns'],values,w['ground_vectors'])


def test_float_or_nonfinite_witness_is_rejected():
    w=packet.load_witness();vectors=w['ground_vectors'].astype(float);vectors[0,0]=float('nan')
    with pytest.raises(ValueError):c.eigensystem('ground',w['ground_columns'],w['ground_energies'],vectors)


def test_complete_finite_dictionary_is_excluded(completed_spectrum):
    r=completed_spectrum['comparison'];assert r['all_three_distinct_frequency_assignments_excluded'] is True
    assert F(r['comparisons'][0]['relative_error_lower'])>F(89,100)
    assert F(r['comparisons'][1]['relative_error_lower'])>F(99,100)
    assert r['physical_rest_masses_derived'] is False
    assert r['new_eligible_empirical_comparisons']==0


def test_minimax_boundary_and_admission():
    assert comparison.minimax_floor(F(2),F(4))==F(1,3)
    assert comparison.minimax_floor(F(4),F(2))==0
    with pytest.raises(ValueError):comparison.minimax_floor(True,F(2))


@pytest.mark.parametrize('key',['physical_particle_species_identified','native_observer_measurement_executed',
                                'finite_energies_identified_as_relativistic_rest_masses'])
def test_physical_promotion_is_rejected(completed_spectrum,key):
    value=deepcopy(completed_spectrum);value['spectrum'][key]=True
    with pytest.raises(ValueError):packet.validate(value)


def test_bool_zero_comparison_laundering_is_rejected(completed_spectrum):
    value=deepcopy(completed_spectrum);value['comparison']['new_eligible_empirical_comparisons']=False
    with pytest.raises(ValueError):packet.validate(value)


def test_integer_false_scope_laundering_is_rejected(completed_spectrum):
    value=deepcopy(completed_spectrum);value['spectrum']['native_observer_measurement_executed']=0
    with pytest.raises(ValueError):packet.validate(value)


def test_full_table_and_inherited_custody_are_explicit(completed_spectrum):
    text=packet.table(completed_spectrum)
    assert len(text.splitlines())==229
    assert completed_spectrum['inherited_inputs']['historical_custody_is_new_certificate_admission'] is False
    assert completed_spectrum['inherited_inputs']['old_all_number_Gram_execution_repeated_here'] is False
