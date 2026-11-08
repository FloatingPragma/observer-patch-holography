"""Actual source concentration, finite spectra and original readout controls."""
from fractions import Fraction as F
import numpy as np
import pytest
from . import uniform_mode as u, uniform_spectral as s


def test_actual_source_operator_identities_include_every_row():
    for n in s.SCALES:
        row=u.source_identity(n)
        assert row['all_source_rows']==256
        assert row['uniform_Ward_operator_identity_exact']
        assert row['double_commutator_first_moment_identity_exact']


def test_joint_uniform_concentration_has_a_zero_energy_centre_and_original_floor():
    previous=None
    for n in(2,4,16,256):
        row=u.concentration(n,F(1,20),1000)
        r=F(120,n)+F(1,n**6)
        assert row['scaled_inverse_weighted_second_moment_upper']==r
        assert row['scaled_third_inverse_weighted_energy_square_variance_upper']==r*r/4
        assert row['original_preparation_and_window_probability_lower']==F(1,160000)
        assert not row['nonzero_momentum_or_particle_shell_proved']
        assert not row['original_continuum_residue_lower_proved']
        if previous is not None:assert r<previous
        previous=r
    with pytest.raises(ValueError):u.concentration(2,F(3,4))


def test_source_locality_gives_the_finite_order_floor_without_allocating_a_ground():
    for n in(2,4):
        row=u.order_box(n)
        assert row['bulk_one_site_marginal']['half_trace_upper']<=F(1,40)
        assert row['excluded_boundary_fraction']<=F(1,64)
        assert row['finite_order_floor']>F(1,20)
        assert row['source_sites']==(2*row['source_radius']+1)**3
        assert not row['ground_matrix_allocated']
        assert not row['large_native_tape_executed']


@pytest.fixture(scope='module')
def result():return s.certificate()


def test_complete_actual_cube_spectra_have_positive_within_charge_weighted_width(result):
    previous=F()
    for row in result['finite_cases']:
        uniform=row['uniform_response']
        lo,hi=map(F,uniform['unscaled_weighted_energy_square_variance_interval'])
        assert F(1,100)<lo<=hi
        slo,shi=map(F,uniform['scaled_weighted_energy_square_variance_interval'])
        assert slo>previous;previous=slo
        assert uniform['charge_pair_exactly_equal_spectral_measures']
        assert F(uniform['each_charge_original_probability_interval'][0])>F(1,5)
        assert not row['reflection_characters_are_continuum_momentum']
        assert row['physical_mass'] is None


def test_all_original_charge_reflection_energy_outcomes_retained(result):
    for row in result['finite_cases']:
        channels=row['local_charge_reflection_energy_readout']
        assert {(x['source_charge'],x['even_cube_reflection_character'])for x in channels}=={(q,k)for q in(-1,1)for k in range(4)}
        records=[r for x in channels for r in x['records']]
        assert len(records)==row['complete_original_records']
        assert any(r['original_probability_interval']==['0.000000000000','0.000000000000']for r in records)
        lo=sum(F(r['original_probability_interval'][0])for r in records)
        hi=sum(F(r['original_probability_interval'][1])for r in records)
        assert lo<=1<=hi
        assert row['every_charge_energy_and_reflection_complement_retained']
        assert not row['native_spectral_pulse_tape_executed']


def test_spectral_proposal_mutations_are_rejected():
    with np.load(s.ROOT/s.WITNESS,allow_pickle=False)as data:
        e,v=data['e2'].copy(),data['v2'].copy()
    with pytest.raises(ValueError,match='complete'):s.verify(2,e[:-1],v)
    permuted=v.copy();permuted[:,[0,1]]=permuted[:,[1,0]]
    with pytest.raises(ValueError,match='isolated'):s.verify(2,e,permuted)
    v[0,0]+=1<<45
    with pytest.raises(ValueError,match='Gram'):s.verify(2,e,v)


def test_independent_dense_source_measure_is_inside_certified_intervals(result):
    h,scale=s.source.matrix(field=F(1,32));energy,v=np.linalg.eigh(h/scale)
    ground=v[:,0];signs=[-1 if x in s.source.cube()[1]else 1 for x in range(8)]
    values=np.array([sum(signs[x]*(1-2*((word>>x)&1))for x in range(8))/8 for word in range(256)])
    weights=(v.T@(values*ground))**2;gaps=energy[1:]-energy[0]
    moments={j:float(np.sum(weights[1:]*gaps**j))for j in(1,-1,-3)}
    row=result['finite_cases'][0]['uniform_response']
    for j,value in moments.items():
        lo,hi=map(float,row['unscaled_original_moments'][str(j)])
        assert lo<=value<=hi
    mean=moments[-1]/moments[-3]
    variance=moments[1]/moments[-3]-mean**2
    lo,hi=map(float,row['unscaled_weighted_energy_square_variance_interval'])
    assert lo<=variance<=hi


def test_fixed_box_gap_diverges_on_joint_clock_instead_of_becoming_continuum(result):
    g=F(result['same_cube_zero_field_certificate']['gap_interval'][0])
    assert g>F(4,5)
    for n in(2,4,16):assert u.fixed_box_gap(n,g)==n*n*g-F(8,n**3)
    assert not result['fixed_box_table_substituted_for_growing_volume_field']
