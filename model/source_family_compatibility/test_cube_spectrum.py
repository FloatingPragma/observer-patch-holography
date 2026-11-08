"""Independent controls for original source clock probabilities and errors."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy import sparse

from . import cube_spectrum as c


@pytest.fixture(scope='module')
def result():
    return c.certificate()


def test_machin_and_all_fourier_roots_have_the_declared_enclosure():
    lo,hi=c.pi_interval()
    decimal=F(314159265358979323846264338327950288,10**35)
    assert lo<decimal+F(1,10**35) and hi>decimal
    assert hi-lo<=F(1,2**78)
    roots=c.roots()
    for j,(real,imag)in enumerate(roots):
        target=np.exp(2j*np.pi*j/len(roots))
        assert abs(complex(real,imag)/2**c.ROOT_BITS-target)<2**-40
    assert roots[0]==(2**c.ROOT_BITS,0)
    assert roots[128]==(-2**c.ROOT_BITS,0)


def test_complete_finite_phase_instrument_and_original_weight():
    # A two-level algorithm control with original branch probability1/8.
    source={'clock':sparse.diags([16,32],dtype=np.int64,format='csr'),
            'row_bound':32,'vectors':np.tile([[2**23],[0]],(1,8)).astype(np.int64)}
    data={'norm_squared':2**46,'half_trace_distance_upper':F(0)}
    read=c.phase_records(data,source,8)
    expected=[]
    for j in range(8):
        amplitude=sum(np.exp(1j*t*(2*np.pi*j/8-.5))for t in range(8))/8
        expected.append(abs(amplitude)**2/8)
    for row in read['values']:
        assert all(abs(float(q)-p)<float(read['per_record_numerical_error'])for q,p in zip(row,expected,strict=True))
        assert abs(float(sum(row,F()))-1/8)<8*float(read['per_record_numerical_error'])
    assert int(np.argmax(expected))==1  # Catches reversal of the Fourier sign.
    # Orthogonality of the ancilla Fourier transform, on both input columns.
    effects=[]
    for j in range(8):
        diagonal=[sum(np.exp(1j*t*(2*np.pi*j/8-energy/2))for t in range(8))/8 for energy in(1,2)]
        op=np.diag(diagonal);effects.append(op.conj().T@op)
    assert np.linalg.norm(sum(effects)-np.eye(2))<1e-14


def test_dyadic_step_is_enclosed_and_overflow_is_rejected():
    h=sparse.csr_matrix(np.array([[0,16],[16,0]],dtype=np.int64))
    real=np.array([[2**46],[0]],dtype=np.int64);imag=np.zeros_like(real)
    r,i=c.phase_step(h,real,imag,16)
    exact=np.array([math.cos(.5),-1j*math.sin(.5)])
    actual=(r[:,0]+1j*i[:,0])/2**46
    assert np.linalg.norm(actual-exact)<float(c.step_error(2))
    with pytest.raises(OverflowError):
        c._action(h,np.full((2,1),np.iinfo(np.int64).max,dtype=np.int64),16)


def test_window_tail_controls_the_actual_wrapped_phase_distribution():
    n=256;distance=F(1,5);bound=float(c.window_tail(n,distance))
    assert c.window_tail(n,distance)==F(55,338)
    for energy in(.125,.301,.481,.735,1.1,5.625):
        phase=energy/2
        errors=(phase-2*np.pi*np.arange(n)/n+np.pi)%(2*np.pi)-np.pi
        numerator=np.sin(n*errors/2);denominator=n*np.sin(errors/2)
        probs=np.ones(n)
        np.divide(numerator,denominator,out=probs,where=np.abs(denominator)>1e-14)
        assert float(np.sum(probs[np.abs(errors)>=float(distance)/2]**2))<=bound
    assert c.window_tail(8,F(1,100))==1
    for bad in(0,3,False):
        with pytest.raises(ValueError):c.window_tail(bad,F(1,5))


def test_outward_rounding_includes_negative_and_positive_rationals():
    for value in(F(-7,13),F(7,13),F(1,10**6),F(0)):
        lo,hi=c.outward(value)
        assert lo<=value<=hi and hi-lo<=F(1,10**6)


def test_actual_full_source_reference_and_complete_record_count(result):
    assert result['source_full_Fock_dimension']==65536
    assert result['computed_charge_block_dimension']==3920
    assert result['c']==F(1,8) and result['U']==1
    assert result['phase_bins']==256 and result['clock_step']==F(1,2)
    assert result['ground_half_trace_error']==F(1,1920)
    assert result['phase_record_numerical_error_upper']==F(14,10**9)
    assert result['charge_and_spin_attachment']['exact_creator_maps']==16
    for row in result['profiles']:
        assert len(row['all_phase_record_probability_intervals'])==256
        ranges=[x['phase_records']for x in row['complete_phase_record_groups']]
        assert [j for lo,hi in ranges for j in range(lo,hi)]==list(range(256))
        assert all(0<=lo<=hi<=1 for lo,hi in row['all_phase_record_probability_intervals'])
        low=sum(lo for lo,_ in row['all_phase_record_probability_intervals'])
        high=sum(hi for _,hi in row['all_phase_record_probability_intervals'])
        assert low<=row['original_addition_weight'][1] and high>=row['original_addition_weight'][0]


def test_actual_profile_energy_means_are_separated_with_original_weights(result):
    first,second=result['profiles'][:2]
    assert first['original_addition_weight'][0]>F(583,1000)
    assert first['original_addition_weight'][1]<F(585,1000)
    assert second['original_addition_weight'][0]>F(415,1000)
    assert second['original_addition_weight'][1]<F(417,1000)
    assert first['conditional_excitation_mean'][1]<F(483,1000)
    assert second['conditional_excitation_mean'][0]>F(5228,10000)
    for row in result['profiles']:
        assert row['conditional_energy_variance'][0]>F(8,1000)
        assert row['conditional_energy_variance'][1]<F(15,1000)
        assert len(row['original_excitation_moment_intervals'])==3
        assert len(row['original_shifted_clock_moments_through_four'])==5


def test_true_source_spectral_windows_have_unconditional_probability_floors(result):
    for row in result['profiles']:
        window=row['true_spectral_windows'][0]
        assert window['excitation_energy_window'][0]==F(1,8)
        assert window['excitation_energy_window'][1]<F(936,1000)
        floor=F(466,1000)if row['profile_pi_multiple'][2]==0 else F(335,1000)
        assert window['original_weight_lower']>floor
        assert window['original_weight_lower']<row['original_addition_weight'][0]
        assert window['phase_tail_upper']==F(55,338)


def test_charge_spin_and_physical_readout_boundaries(result):
    assert result['removal_profile_correspondence']=='k maps to k xor (1,1,1)'
    assert result['spin_components_have_equal_ground_response']
    assert result['complete_addition_removal_instrument']
    assert result['physical_service_failures_retained_with_separate_error_budget']
    assert not result['phase_labels_are_exact_momentum']
    assert not result['phase_record_is_sharp_spectral_projector']
    assert not result['native_large_clock_tape_executed']
    assert result['physical_mass_assigned']is None
    assert not result['continuum_spectral_atom_inferred']


def test_certificate_returns_detached_nested_payloads(result):
    immutable=repr(result)
    changed=c.certificate();changed['profiles'][0]['original_addition_weight'][0]._numerator=0
    changed['profiles'][0]['true_spectral_windows'][0]['original_weight_lower']._numerator=0
    assert repr(c.certificate())==immutable and repr(result)==immutable


def test_output_envelopes_are_checked_and_cannot_supply_source_probabilities(result,monkeypatch):
    import copy
    corrupted=copy.deepcopy(c._ADMITTED_PROFILES)
    corrupted[0]['original_addition_weight']=('0','1/10')
    monkeypatch.setattr(c,'_ADMITTED_PROFILES',corrupted)
    with pytest.raises(ArithmeticError,match='fixed output envelope'):
        c._admit_profile(result['profiles'][0])


def test_low_precision_cannot_claim_the_fixed_readout_tolerance(monkeypatch):
    source={'clock':sparse.diags([16,32],dtype=np.int64,format='csr'),
            'row_bound':32,'vectors':np.tile([[2**23],[0]],(1,8)).astype(np.int64)}
    monkeypatch.setattr(c,'DEGREE',1)
    with pytest.raises(ArithmeticError):
        c.phase_records({'norm_squared':2**46,'half_trace_distance_upper':F(0)},source,8)
