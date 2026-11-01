from fractions import Fraction as F

import numpy as np
import pytest

from model.coarse_energy_ward.ward import Clifford, annihilator, commutator
from model.full_fock_field.field import FieldGeometry, FullFock
from model.thermodynamic_ground import local, volume


def test_complete_odd_and_even_local_algebra_stabilizes_on_actual_three_dimensional_source():
    result=local.source_certificate()
    assert result['complete_one_site_CAR_basis_checked']==16
    assert result['odd_probes_retained'] is True
    assert result['missing_incident_bond_control_nonzero'] is True
    assert result['quartic_majorana_terms']>0


def test_car_embedding_is_multiplicative_including_odd_fields():
    small,large=FieldGeometry((2,1,1)),FieldGeometry((2,2,1))
    a,b=annihilator(0),annihilator(3)
    assert (local.embed(a@b,small,large)-local.embed(a,small,large)@local.embed(b,small,large)).is_zero()
    assert (local.embed(a.adjoint(),small,large)-local.embed(a,small,large).adjoint()).is_zero()
    # Disjoint odd fields anticommute; silently commuting them is wrong.
    assert (a@b+b@a).is_zero()
    assert not commutator(a,b).is_zero()


def test_quartic_derivation_against_independent_full_fock_matrices():
    space=FullFock(4)
    a=[space.annihilator(j).toarray() for j in range(4)]
    n=[x.conj().T@x for x in a]
    d=n[0]@n[1]+n[2]@n[3]
    probe=a[0].conj().T@a[2]
    target=probe@(n[1]-n[3])
    assert np.array_equal(d@probe-probe@d,target)
    assert np.linalg.norm(target,2)==1
    # A vacuum expectation of the witness vanishes despite its nonzero
    # represented operator; functional faithfulness cannot be inferred.
    assert target[0,0]==0
    assert np.count_nonzero(target)>0


def test_ground_inequality_really_requires_full_number_minimum():
    geometry=FieldGeometry((2,1,1));space=FullFock(4)
    h=space.hamiltonian(geometry).toarray()
    # Exact binary-valued negative one-particle vector above the vacuum.
    v=np.zeros(16,complex);v[1]=1;v[8]=1j
    assert np.array_equal(h@v,-v/16)
    assert np.vdot(v,h@v)==-F(1,8)
    # N=0 is the canonical ground of its sector but fails the odd local
    # lowering condition; no such reference is used by the theorem.
    odd=space.annihilator(0).conj().T.toarray()+1j*space.annihilator(3).conj().T.toarray()
    q=odd.conj().T@(h@odd-odd@h)
    assert q[0,0]==-F(1,8)


def test_stationarity_budget_is_state_uniform_and_decreases_with_buffer():
    bounds=[local.stationarity_error(r,0,1,F(1,8)) for r in (16,32,64)]
    assert 0<bounds[2]<bounds[1]<bounds[0]<1
    assert local.stationarity_error(16,0,1,0)==0
    assert local.stationarity_error(16,0,1,F(-1,8))==bounds[0]
    assert local.stationarity_error(16,0,1,1,speed=0)==0


def test_averaging_boundary_fraction_retains_bad_translates():
    result=volume.averaged_dynamics_error(64,16,0,1,F(1,8))
    bad=1-F(48**3,64**3)
    assert result==2*local.stationarity_error(16,0,1,F(1,8))+2*bad
    assert result>2*bad
    assert volume.averaged_dynamics_error(4096,64,0,1,F(1,8))<F(1,10)
    assert volume.ground_inequality_defect(8,0,1)==F(7,2)*(1-F(7**3,8**3))
    assert volume.ground_inequality_defect(1,2,3)==F(21,2)


def test_translation_defect_matches_direct_symmetric_difference():
    from itertools import product
    for radius,shift in [(0,(1,0,0)),(2,(1,-2,0)),(2,(9,0,0)),(3,(0,0,0))]:
        box=set(product(range(-radius,radius+1),repeat=3))
        moved={tuple(a+b for a,b in zip(x,shift,strict=True)) for x in box}
        assert local.folner_translation_error(radius,shift)==F(len(box^moved),len(box))


def test_density_uses_unshifted_full_ground_and_not_assumed_half_filling():
    result=volume.density_certificate()
    assert result['number_density_lower']=='1/4'
    assert result['number_density_upper']=='1'
    assert result['particle_hole']['half_filling_concluded'] is False
    assert result['kinetic']['one_body_source_dimension']==54
    assert result['density_selected_uniquely'] is False


@pytest.mark.parametrize('call',[
    lambda:local.stationarity_error(2,2,1,F(1,8)),
    lambda:local.folner_translation_error(True,(0,0,0)),
    lambda:local.folner_translation_error(2,(0,.5,0)),
    lambda:local.folner_translation_error(2,(0,0,0),-1),
    lambda:volume.averaged_dynamics_error(2,2,0,1,F(1,8)),
    lambda:volume.ground_inequality_defect(1,0,2),
    lambda:volume.ground_inequality_defect(1,0,1,speed=.125),
    lambda:local.embed(Clifford({1<<100:1}),FieldGeometry((1,1,1)),FieldGeometry((2,1,1))),
])
def test_unadmitted_domains_fail_before_certificate(call):
    with pytest.raises(ValueError):call()
