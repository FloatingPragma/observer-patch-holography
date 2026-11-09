"""Independent source incidence and physical entropy comparison controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
import json

import numpy as np
import pytest

from . import model as m


def full_count(side,cells):
    selected=set(cells);inside=set();cross=set();plaquettes=[]
    for x in product(range(side),repeat=3):
        for i in range(3):
            y=tuple(t+(k==i) for k,t in enumerate(x))
            if y[i]>=side:continue
            edge=tuple(sorted((x,y)))
            if x in selected and y in selected:inside.add(edge)
            elif (x in selected)!=(y in selected):cross.add(edge)
        for i,j in combinations(range(3),2):
            if x[i]+1>=side or x[j]+1>=side:continue
            vertices=tuple(tuple(t+u*(k==i)+v*(k==j) for k,t in enumerate(x))
                for u,v in product((0,1),repeat=2))
            edges={tuple(sorted((u,v))) for u,v in combinations(vertices,2)
                   if sum(abs(a-b) for a,b in zip(u,v))==1}
            plaquettes.append((tuple(sorted(vertices)),edges))
    return cross,{p for p,edges in plaquettes if 0<len(edges&inside)<4}


def test_actual_mixed_supports_against_independent_full_box():
    for side in (2,3,4):
        allcells=tuple(product(range(side),repeat=3))
        regions=[((0,0,0),),((0,0,0),(1,0,0)),allcells,
                 tuple(x for x in allcells if sum(x)%2==0),
                 tuple(x for x in allcells if x[0]==0 or x[1]==0)]
        for cells in regions:
            value=m.interaction_cut(side,cells);b,p=full_count(side,cells)
            assert value['crossing_hopping_bonds']==len(b)
            assert set(value['mixed_plaquettes'])==p
            assert len(p)<=2*len(b)
        two=m.interaction_cut(side,((0,0,0),(1,0,0)))
        assert two['mixed_plaquette_count']==2
        assert two['crossing_hopping_bonds']==(4 if side==2 else 5)


def marginals(rho,dimensions,multiplicities):
    W,nr,no,_=m.information.singlet_embedding(dimensions,multiplicities)
    t=(W@rho@W.conj().T).reshape(nr,no,nr,no)
    return np.einsum('abcb->ac',t),np.einsum('abad->bd',t)


def test_comparison_preserves_physical_marginals_and_schur_entropy_identity():
    dims=[2,3];mult=[(2,2),(1,2)];size=6
    seed=np.arange(1,size*size+1).reshape(size,size)+1j*np.eye(size)
    rho=seed@seed.conj().T+np.eye(size);rho/=np.trace(rho)
    sigma,stats=m.decorrelate(rho,dims,mult)
    for a,b in zip(marginals(rho,dims,mult),marginals(sigma,dims,mult)):
        assert np.allclose(a,b,atol=2e-13)
    assert stats['center_coherence']>0
    assert stats['conditional_mutual_information']>0
    assert stats['entropy_gain']==pytest.approx(m.information.entropy(sigma)-m.information.entropy(rho))
    # The singlet remains entangled in its representation factors.
    stats2=m.information.information_control(sigma,dims,mult)
    assert abs(stats2['conditional_multiplicity_mutual_information'])<1e-12
    assert abs(stats2['center_coherence'])<1e-12
    assert stats2['representation_entropy']>0


def test_comparison_is_not_a_channel_on_unknown_states():
    r1=np.diag([.7,.1,.1,.1]);r2=np.diag([.1,.1,.1,.7])
    f=lambda r:m.decorrelate(r,[1],[(2,2)])[0]
    assert not np.allclose(f((r1+r2)/2),(f(r1)+f(r2))/2)


def test_noncommuting_boundary_control_uses_actual_polynomial_logarithm():
    H0=np.diag([0.,.2,.4,.6,.1,.3]);V=np.zeros((6,6),complex)
    V[0,3]=.2j;V[3,0]=-.2j;V[2,5]=V[5,2]=.15
    value=m.control(H0+V,H0,1.,4,2,[2,3],[(2,2),(1,2)])
    assert value['entropy_gain']>0
    assert value['relative_entropy']>=-1e-12
    assert abs(value['comparison_identity_residual'])<1e-12
    assert value['source_state_replaced_by_comparison'] is False


def test_polynomial_error_is_required_even_without_an_interaction():
    value=m.noninteracting_control()
    assert value['boundary_operator_norm']==0
    assert value['entropy_gain']>0
    assert F(value['exact_product_determinant'])==F(-58849,2802**2)
    assert value['entropy_gain']<=value['actual_modular_oscillation']
    assert abs(value['energy_change'])<1e-12


@pytest.mark.parametrize('level',(0,1,2,8))
def test_source_bound_uses_crossing_hopping_and_mixed_plaquettes_only(level):
    member=m.contract.member(level);bound=m.bind_member(member);m.verify_bound(member,bound)
    row=bound['result'];a=F(member['source']['spacing'])
    b=4 if level==0 else 5
    assert F(row['actual_crossing_interaction_norm_upper'])==(30*b+F(2,7))/a
    assert row['cut_interactions']['all_electric_terms_belong_to_one_regional_factor']
    assert not row['positive_finite_continuum_area_coefficient_proved']
    assert not row['comparison_replaces_original_source_preparation']
    assert row['physical_mass_table_GeV']==[]
    bad=deepcopy(bound);bad['result']['Einstein_dynamics_proved']=True
    bad['result_sha256']=m.contract.digest(bad['result'])
    with pytest.raises(ValueError):m.verify_bound(member,bad)
