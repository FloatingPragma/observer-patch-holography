"""Complete-carrier controls and independent finite symmetry consequences."""
from itertools import product
import numpy as np
import pytest
import sympy as s
from .symmetry import coordinate_certificate,word_certificate,transform,commutator
from .flatness import certificate
from model.openwave_neutral_candidate.source import translate_key,full_column


def test_actual_full_link_and_all_source_words():
    assert coordinate_certificate()['complete_coordinate_inversion_identities']==59
    assert coordinate_certificate()['complete_endpoint_current_inversion_identities']==24
    assert word_certificate()['all_axis_hopping_word_identities']==72


def test_required_factors_and_boundary_mutations():
    key=(1<<11,(0,)*24,(0,)*8)
    assert not commutator(key,0)
    assert commutator(key,0,parity=False)
    assert commutator(key,0,invert=False)
    assert commutator(key,0,action=lambda z:full_column(z,periodic=False))
    with pytest.raises(ValueError):transform(key,3)


def test_odd_parity_is_required_for_momentum_flip():
    for mask in (1<<11,(1<<11)|(1<<48)):
        key=(mask,tuple(j%119 for j in range(24)),tuple(j%5 for j in range(8)))
        ck,a=transform(key,1);tc,b=translate_key(ck,1,periodic=True)
        tk,c=translate_key(key,1,periodic=True);ct,d=transform(tk,1)
        assert tc==ct and a*b==(-1)**mask.bit_count()*c*d


def test_independent_symmetry_spectral_and_position_consequences():
    # A dense representation of the proven algebra, independent of the
    # native reduced assembly. Internal energies are deliberately unequal.
    positions=list(product((0,1),repeat=3));internal=np.array([[.3,.17j],[-.17j,.8]])
    h=np.kron(np.eye(8),internal);cs=[];ts=[]
    for axis in range(3):
        c=np.kron(np.diag([(-1)**p[axis] for p in positions]),np.eye(2))
        t=np.zeros((16,16))
        for y,p in enumerate(positions):
            z=list(p);z[axis]^=1;x=positions.index(tuple(z))
            t[2*x:2*x+2,2*y:2*y+2]=np.eye(2)
        assert np.allclose(c@t,-t@c) and np.allclose(c@h,h@c)
        cs.append(c);ts.append(t)
    spectra=[]
    for k in positions:
        e=np.eye(16)
        for i in range(3):e=e@(np.eye(16)+(-1)**k[i]*ts[i])/2
        vals,vec=np.linalg.eigh(e);basis=vec[:,vals>.5]
        spectra.append(np.linalg.eigvalsh(basis.conj().T@h@basis))
    assert all(np.allclose(v,spectra[0]) for v in spectra)
    assert spectra[0][1]-spectra[0][0]>.5
    # An admissible translation-invariant hopping perturbation that violates
    # C changes the eight momentum spectra, exposing the exact missing premise.
    mutant=h+.02*sum(ts)
    assert np.linalg.norm(cs[0]@mutant-mutant@cs[0])>.1
    assert np.linalg.eigvalsh(mutant)[0]<spectra[0][0]-.05


def test_full_rank_band_scope_and_exact_zero_remainder():
    result=certificate()
    assert result['reduced_algebra']['all_momentum_sector_ranks']==[128]*8
    assert result['exact_all_order_momentum_splitting']=='0'
    assert s.Rational(result['generic_cubic_remainder_upper_if_only_second_order_were_used'])>0
    assert result['all_order_flatness_remainder'].startswith('0 exactly')
    assert not result['stationary_neutral_energy_reinterpreted_as_a_moving_particle_mass']
