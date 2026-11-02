from fractions import Fraction as F
import copy
import pytest
import sympy as s
from . import algebra, bounds, packet


def test_complete_continuity_and_torque():
    c=algebra.exact_checks()
    assert c['all_axes_continuity'] is True
    assert c['all_axes_staggered_torque'] is True


def test_positive_field_torque_cannot_be_omitted():
    sp=algebra._spins(2)
    field=-(sp[0][0]-sp[1][0])/4
    transverse=sp[0][1]+sp[1][1]
    actual=algebra._comm(field,transverse)
    assert not algebra._zero(actual)
    assert algebra._zero(actual+s.I*(sp[0][2]-sp[1][2])/4)


@pytest.mark.parametrize('index,spin,current,energy',[(0,F(1,6),F(1,3),2),(1,F(2,3),F(1,6),1)])
def test_exact_ring_spectral_measure(index,spin,current,energy):
    row=algebra.exact_checks()['four_site_ring_spectral_lines'][index]
    assert [int(r['energy']) for r in row['lines']]==[0,1,2,3]
    assert [(F(r['spin_weight']),F(r['current_weight'])) for r in row['lines'] if int(r['energy'])==energy]==[(spin,current)]
    assert sum(F(r['current_weight'])/int(r['energy']) for r in row['lines'] if int(r['energy']))==F(1,6)
    for r in row['lines']:
        assert F(r['energy'])**2*F(r['spin_weight'])==F(row['lattice_factor'])*F(r['current_weight'])


def test_open_bond_relaxation_cancels_positive_bare_twist():
    c=algebra.exact_checks()['open_bond_twist']
    assert F(c['bare_second_derivative'])==F(1,2)
    assert F(c['twice_inverse_current_moment'])==F(1,2)
    assert c['ground_energy_curvature']=='0'


def test_exact_cache_is_detached():
    x=algebra.exact_checks();x['four_site_ring_spectral_lines'][0]['lines'][0]['spin_weight']='1'
    assert algebra.exact_checks()['four_site_ring_spectral_lines'][0]['lines'][0]['spin_weight']=='0'


@pytest.mark.parametrize('m,k,r',[(F(1,1000),F(1,2),1),(F(1,2),F(1,8),F(1,4)),(F(1,10),F(1,4),8)])
def test_soft_current_fraction_has_full_response_denominator(m,k,r):
    c=bounds.soft_fraction(m,k,r)
    assert c['soft_first_moment_lower']==c['field_band_lower']*c['field_weight_lower']
    assert c['soft_current_fraction_lower']==c['soft_first_moment_lower']/c['total_first_moment_upper']
    assert c['soft_current_fraction_lower']>=c['universal_fraction_lower']>0
    assert c['open_source_band_lower']<c['field_band_lower']<c['field_band_upper']<c['open_source_band_upper']


@pytest.mark.parametrize('args',[(True,F(1,2),1),(F(1,2),F(1,16),1),(F(1,2),F(1,2),0),(0.1,F(1,2),1)])
def test_bad_fraction_inputs_fail(args):
    with pytest.raises(ValueError):bounds.soft_fraction(*args)


@pytest.mark.parametrize('degree,scale',[(2,2),(3,4),(4,16)])
def test_fractional_tail_bound(degree,scale):
    c=bounds.fractional_tail(F(1,4),degree,scale)
    assert c['alpha']==1-F(1,degree)
    assert c['cutoff']==scale**degree
    assert c['fractional_tail_upper']==F(1,4*scale)
    # Equality at the cutoff for a positive single-atom measure with first moment D.
    weight=F(1,4)/c['cutoff']
    assert weight*scale**(degree-1)==c['fractional_tail_upper']


def test_first_energy_endpoint_not_admitted():
    with pytest.raises(ValueError):bounds.fractional_tail(F(1,4),1,2)
    with pytest.raises(ValueError):bounds.fractional_tail(F(1,4),True,2)
    with pytest.raises(ValueError):bounds.fractional_tail(F(1,4),2,False)


def test_moments_and_continuity_allow_nonvanishing_energy_escape():
    for n in (3,16,1024):
        c=bounds.escaping_current(n)
        assert c['first_moment']==F(1,8)
        assert c['escaping_first_moment']==c['high_energy']*c['high_mass']==F(1,16)
        assert c['inverse_moment_upper']<F(1,12)
        assert c['low_continuous_mass']>=F(1,4)/12
        assert c['weak_limit_has_a_positive_energy_atom'] is False
        assert c['is_a_Heisenberg_ground_state_construction'] is False


def test_packet_source_scope_and_types():
    p=packet.build();packet.validate(p)
    for k in ['first_energy_uniform_integrability_proved','actual_ground_hard_current_tail_vanishes_proved',
              'sharp_mass_shell_proved','physical_species_or_particle_masses_identified']:
        bad=copy.deepcopy(p);bad['theorem'][k]=True
        with pytest.raises(ValueError):packet.validate(bad)
    bad=copy.deepcopy(p);bad['theorem']['sharp_mass_shell_proved']=0
    with pytest.raises(ValueError):packet.validate(bad)


def test_packet_rejects_spectral_weight_corruption():
    bad=packet.build()
    bad['exact_finite_checks']['four_site_ring_spectral_lines'][0]['lines'][2]['current_weight']='2/3'
    with pytest.raises(ValueError):packet.validate(bad)


def test_source_snapshot_race_is_rejected(monkeypatch):
    original=packet.bindings;calls=0
    def moving():
        nonlocal calls
        calls+=1
        result=original()
        if calls>1:result['model/__init__.py']='0'*64
        return result
    monkeypatch.setattr(packet,'bindings',moving)
    with pytest.raises(RuntimeError,match='source changed'):packet.build()
