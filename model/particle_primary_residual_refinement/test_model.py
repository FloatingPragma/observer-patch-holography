from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path

import pytest
import sympy as s

from model.particle_common_contract import contract
from model.particle_primary_variational_readout import variational
from model.particle_local_wilson_density import density
from model.particle_primary_residual_refinement import model as m


@lru_cache(maxsize=1)
def beta_basis_gamma():
    native=density.gamma_data();beta=native[3]
    positive=s.GramSchmidt((beta-s.eye(4)).nullspace(),True)
    negative=s.GramSchmidt((beta+s.eye(4)).nullspace(),True)
    change=s.Matrix.hstack(*(positive+negative))
    assert s.simplify(change.H*change)==s.eye(4)
    gamma=tuple(s.simplify(change.H*g*change) for g in native)
    assert gamma[3]==s.diag(1,1,-1,-1)
    return gamma


def hopping(sites,edges,signs=None):
    gamma=beta_basis_gamma();beta=gamma[3];h=s.zeros(4*sites)
    signs=signs or [1]*len(edges)
    for (x,y,axis),sign in zip(edges,signs):
        block=sign*(-beta+beta*gamma[axis])/2
        h[4*x:4*x+4,4*y:4*y+4]=block
        h[4*y:4*y+4,4*x:4*x+4]=block.H
    occupied=[4*x+j for x in range(sites) for j in range(4) if beta[j,j]==-1]
    assert len(occupied)==2*sites and beta==s.diag(*beta.diagonal())
    return h,occupied


def blocks(h,occupied):
    empty=[i for i in range(h.rows) if i not in occupied]
    A=h.extract(empty,empty);B=h.extract(empty,occupied);D=h.extract(occupied,occupied)
    C=A*B-B*D
    return s.trace(B.H*B),s.trace((B.H*B)**2),s.trace(C.H*C)


def car_action(h,state):
    """Literal full occupation action; no one-particle Hamiltonian substitution."""
    out={}
    terms=[(i,j,h[i,j]) for i in range(h.rows) for j in range(h.cols) if h[i,j]!=0]
    for mask,amplitude in state.items():
        for i,j,value in terms:
            if not mask&(1<<j):continue
            sign=(-1)**((mask&((1<<j)-1)).bit_count())
            reduced=mask^(1<<j)
            if reduced&(1<<i):continue
            sign*=(-1)**((reduced&((1<<i)-1)).bit_count())
            target=reduced|(1<<i)
            out[target]=out.get(target,0)+sign*value*amplitude
    return {key:s.expand(value) for key,value in out.items() if s.expand(value)!=0}


def norm2(state):return s.simplify(sum(s.conjugate(x)*x for x in state.values()))


def test_literal_two_site_CAR_has_all_pair_sectors():
    for axis in range(3):
        h,occ=hopping(2,[(0,1,axis)]);vac={sum(1<<i for i in occ):s.Integer(1)}
        one=car_action(h,vac);two=car_action(h,one)
        z,R,C=blocks(h,occ)
        assert z==1 and R==s.Rational(1,4) and C==0
        assert norm2(two)==3*z*z-2*R+C==s.Rational(5,2)
        sectors={k:{mask:value for mask,value in two.items()
                    if sum(bool(mask&(1<<i)) for i in range(h.rows) if i not in occ)==k}
                 for k in range(3)}
        assert [norm2(sectors[k]) for k in range(3)]==[1,0,s.Rational(3,2)]


def test_weighted_disjoint_native_CAR_mixed_identity():
    h0,occ0=hopping(2,[(0,1,1)])
    h=s.diag(h0,h0);weights=[s.Rational(1,3),s.Rational(2,5)]
    ht=s.diag(weights[0]*h0,weights[1]*h0)
    occ=occ0+[8+i for i in occ0];vac={sum(1<<i for i in occ):s.Integer(1)}
    actual=norm2(car_action(h,car_action(ht,vac)))
    zt=sum(t*t for t in weights);q=sum(weights)
    assert actual==2*zt+2*q*q-s.Rational(1,2)*zt
    assert h*ht==ht*h


@pytest.mark.parametrize('sites,edges',[
    (3,[(0,1,0),(1,2,0)]),
    (3,[(0,1,0),(1,2,1)]),
    (4,[(0,1,0),(0,2,1),(1,3,1),(2,3,0)]),
])
def test_exact_link_sign_average_checks_two_step_Haar_counts(sites,edges):
    # Each nontrivial elementary closed-loop monomial has an edge with odd
    # sign exponent. The exact finite sign average tests this cancellation.
    neighbors={x:[] for x in range(sites)}
    for x,y,axis in edges:
        neighbors[x].append((y,axis));neighbors[y].append((x,axis))
    degrees=[len(neighbors[x]) for x in range(sites)]
    S2=sum(d*d for d in degrees);B=len(edges)
    straight=sum(1 for x in neighbors for y,axis in neighbors[x]
                 for z,other in neighbors[y] if z!=x and axis==other)
    averages=[s.Integer(0)]*3
    for signs in product((-1,1),repeat=len(edges)):
        h,occ=hopping(sites,edges,signs)
        values=blocks(h,occ)
        averages=[a+v for a,v in zip(averages,values)]
    averages=[s.simplify(a/2**len(edges)) for a in averages]
    assert averages==[B,s.Rational(S2-B,4),s.Rational(S2-2*B+straight,4)]


@pytest.mark.parametrize('L',[2,3,4])
def test_open_cube_counts_from_literal_coordinates(L):
    vertices=list(product(range(L),repeat=3));degree={v:0 for v in vertices};straight=0;bonds=0
    for v in vertices:
        for axis in range(3):
            for sign in (-1,1):
                degree[v]+=int(0<=v[axis]+sign<L)
                straight+=int(0<=v[axis]+2*sign<L)
            bonds+=int(v[axis]+1<L)
    row=m.geometry(L)
    assert row['bonds']==bonds
    assert row['sum_squared_interior_degrees']==sum(d*d for d in degree.values())
    assert row['directed_straight_two_step_paths']==straight
    assert row['hopping_connected_cumulant_coefficient']==F(3*L*(L-1)*(5*L-4),2)


def test_fourth_character_integrals_have_independent_exact_singlet_count():
    row=m.fourth_character_control()
    assert row['integral_S_fourth']==714
    assert row['integral_S_squared_xr_xs']==row['independent_fourfold_tensor_singlet_matrix']
    assert row['integral_S_squared_xr_xs']==[
        [40,20,36,30,18],[20,36,34,24,16],[36,34,72,36,20],
        [30,24,36,44,16],[18,16,20,16,22]]


def test_unit_weights_recover_exact_full_parent_fourth_and_scaling():
    row=m.mixed_parent_moment(2,[1]*5,[1]*5)
    assert row['exact_hopping_mixed_fourth']==96930
    assert row['exact_magnetic_mixed_fourth']==F(6*714+30*300,210**4)
    other=m.mixed_parent_moment(2,[F(2,3)]*5,[F(2,3)]*5)
    assert other['exact_parent_mixed_fourth']==F(4,9)*row['exact_parent_mixed_fourth']
    zero=m.mixed_parent_moment(2,[0]*5,[0]*5)
    assert zero['exact_parent_mixed_fourth']==0


def test_actual_member_zero_residual_improves_without_changing_trial():
    member=contract.member(0);producer=variational.bind_member(member)
    row=m.refined_bound(member,producer)
    assert F(523,100)<F(row['source_residual_upper'])<F(524,100)
    assert F(row['certified_radius_bound_reduction_factor'])>19
    assert row['trial_rayleigh_quotient']==producer['result']['trial_rayleigh_quotient']
    assert row['link_degree_one_bound_uses_contractive_final_compression']
    assert not row['squared_residual_is_exact_on_selected_cutoff']
    assert not row['trial_subspace_invariant'] and not row['original_primary_probe_residue_computed']
    fake=deepcopy(producer);fake['result']['dimensionless_trial']['rational_hopping_coefficients'][0]='0'
    fake['result_sha256']=contract.digest(fake['result'])
    with pytest.raises(ValueError):m.refined_bound(member,fake)


def test_higher_word_cutoff_retains_complete_second_action():
    member=contract.member(1);row=m.refined_bound(member)
    assert row['squared_residual_is_exact_on_selected_cutoff']
    assert not row['link_degree_one_bound_uses_contractive_final_compression']
    assert not row['ground_isolation_or_particle_mass_proved']


def test_input_domains_and_sources():
    for value in (True,1,0,-1,2.0):
        with pytest.raises(ValueError):m.geometry(value)
    with pytest.raises(ValueError):m.mixed_parent_moment(2,[1]*4,[1]*5)
    with pytest.raises(ValueError):m.mixed_parent_moment(2,[1.0]*5,[1]*5)
    for path in m.SOURCE_INPUTS:assert Path(path).is_file(),path
