"""Full source identities, complete response reconstruction and hostile controls."""
import pytest
import sympy as s
from . import model,resolvent


def test_full_source_separate_spin_and_quark_response_cancellation():
    result=model.separate_spin()
    assert result['independent_quark_and_lepton_spin_H_commutators']==6
    assert result['opposite_sector_elementary_field_commutators']==90


def test_complete_charged_images_and_exact_target_annihilators():
    for family,dimensions,degree in (('ec',(17,7),9),('L',(33,12),17)):
        data=resolvent.response(family)
        assert tuple(b['complete_charged_irrep_dimension'] for b in data['blocks'])==dimensions
        assert data['q'].degree()==degree
        block=data['blocks'][0]
        with pytest.raises(ValueError):model.validate_annihilator(block['H'],data['q']+1)
    with pytest.raises(ValueError):model.transfer('ec',2)


def test_pair_resolvent_sign_and_two_spectral_sums():
    z=s.Symbol('z');q=s.Poly(z*(z-2),z,domain=s.QQ);n=s.Poly(-2,z,domain=s.QQ)
    # Signed weights (1,-1) at energies (0,2), offset 3.
    actual,proof=resolvent.pair_functional(q,n,(s.QQ(0),s.QQ(-2)),s.QQ(3))
    assert actual==s.QQ(8,105)
    assert actual==s.Rational(1,3)-s.Rational(2,5)+s.Rational(1,7)
    with pytest.raises(s.NotInvertible):resolvent.pair_functional(q,n,(s.QQ(0),s.QQ(-2)),s.QQ(0))
    with pytest.raises(ValueError):resolvent.pair_functional(q,n,(s.QQ(0),),s.QQ(3))


def test_exact_positive_exchange_and_both_native_lepton_channels():
    result=resolvent.exchange();lower,upper=map(s.Rational,result['exact_J_interval'])
    assert 0<lower<upper and upper-lower<=s.Rational(1,2**47)
    assert result['native_internal_multiplicities']=={'ec':1,'L':2}
    for family in ('ec','L'):
        value,proof=resolvent.pair_integral(family)
        q=proof['response']['denominator']
        check=(proof['reflected_denominator']*proof['inverse']).rem(q)
        assert check==s.Poly(1,q.gen,domain=q.domain)


def test_independent_complete_generalized_spectral_sum():
    """A floating check has a fixed tolerance and supplies no certified data."""
    import numpy as np
    from scipy.linalg import eigh
    from ..ground_composite_response.model import source
    sb=source(1)
    def array(matrix):return np.array(matrix.tolist(),dtype=float)
    energies,vectors=eigh(array(sb['gram']*sb['H']),array(sb['gram']))
    energy=energies[0];ground=vectors[:,0];value=0.0
    for family,multiplicity in (('ec',1),('L',2)):
        all_energies=[];all_weights=[]
        for spin in range(2):
            block=model.transfer(family,spin)
            levels,basis=eigh(array(block['gram']*block['H']),array(block['gram']))
            image=array(block['map'])@ground
            amplitudes=basis.T@array(block['gram'])@image
            all_energies.extend(levels);all_weights.extend(((-1)**spin)*amplitudes**2)
        levels=np.asarray(all_energies);weights=np.asarray(all_weights)
        denominators=1+levels[:,None]+levels[None,:]-2*energy
        assert denominators.min()>0
        value+=multiplicity*np.sum(weights[:,None]*weights[None,:]/denominators)/4
    lower,upper=map(s.Rational,resolvent.exchange()['exact_J_interval'])
    assert float(lower)-1e-10<value<float(upper)+1e-10
