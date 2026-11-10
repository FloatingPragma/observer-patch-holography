"""Cross-check the positive clock against literal source and spatial actions."""
from fractions import Fraction as F
import pytest
import sympy as s

from .model import (JointObserver,combine,dark_wave,local_action,
                    FockObserver,dark_pair,fock_local_action,lattice,
                    linear_combination,h_action,energy_moments,dark_plane_wave,
                    pair_energy_action,pair_hopping_laplacian,PairFockState,MODES,mode)
from ..species.spatial_fock import gamma as pair_gamma,positive_basis
from .certificate import certificate


def test_literal_dark_energy_and_nonadditive_pair():
    observer=JointObserver();state=dark_wave(observer,3)
    assert state.norm_squared()==1
    assert local_action(observer,state)==combine(state,state,2,0)
    fock=FockObserver();pair=dark_pair(fock)
    assert pair.norm_squared()==1
    assert fock.step(pair)==pair
    assert fock_local_action(fock,pair).norm_squared()==0


def test_positive_energy_form_and_complete_readout():
    source=lattice.packet();forward=lattice.onsite(source)
    minus=linear_combination((source,forward),(F(1,2),F(-1,2)))
    plus=linear_combination((source,forward),(F(1,2),F(1,2)))
    a,_=h_action(source,4)
    assert minus.norm_squared()+plus.norm_squared()==1
    assert lattice.overlap(source,a)==(2*minus.norm_squared(),0)
    assert minus.norm_squared()==F(59,60)


def test_energy_is_hermitian_and_whole_cell_translation_covariant():
    first=lattice.packet();second=lattice.onsite(first)
    h_first=h_action(first,4);h_second=h_action(second,4)
    for left,right in zip(h_first,h_second,strict=True):
        actual=lattice.overlap(first,right)
        reverse=lattice.overlap(second,left)
        assert actual==(reverse[0],-reverse[1])
    shift=(1,0,1)
    h_shifted=h_action(lattice.translate(second,4,shift),4)
    assert h_shifted==tuple(lattice.translate(state,4,shift) for state in h_second)


def test_full_spatial_action_has_exact_band_edges():
    for labels,energy in (((0,0,0),2),((1,1,1),2+3*s.pi/2),((2,2,2),2+3*s.pi)):
        state=dark_plane_wave(4,labels,4)
        moments=energy_moments(state,4)
        assert s.expand(moments['mean']-energy)==0
        assert moments['variance']==0


def test_source_energy_and_positive_time_response_are_not_mass_claims():
    result=certificate()
    record=result['finite_time_flux_measurement']
    assert 0<F(record['yes_probability_lower_bound'])<F(record['yes_probability_upper_bound'])<1
    assert result['dark_onsite_energy']=='2'
    assert not result['scope']['physical_rest_mass_identified']
    assert not result['clock_equals_discrete_seven_layer_tick']


def test_energy_domain_guards():
    observer=JointObserver()
    with pytest.raises(ValueError):dark_wave(observer,5)
    with pytest.raises(ValueError):dark_plane_wave(4,(0,0,4))
    with pytest.raises(ValueError):h_action(lattice.packet(),3)
    with pytest.raises(ValueError):energy_moments(linear_combination((lattice.packet(),),(2,)),4)


def test_two_cell_energy_preserves_both_constraints_and_pauli_signs():
    index=mode(5,3,0,0)
    state=positive_basis((1,3),(index,MODES+mode(6,0,0,0)))
    image=pair_energy_action(state)
    assert all(pair_gamma(image,cell)==image for cell in (0,1))
    # With both copies of one mode occupied, the hopping off-diagonal terms
    # vanish by Pauli blocking. Its dGamma(I-S) contribution is exactly2.
    same=positive_basis((1,3),(index,MODES+index))
    assert pair_hopping_laplacian(same)==PairFockState({k:2*v for k,v in same.amplitudes},same.norm_denominator)
    reference=PairFockState({((1,1),()):1})
    assert pair_energy_action(reference).norm_squared()==0
