"""Literal source packet checks and adversarial spectral-scale controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
import sympy as s

from . import scale as v


@pytest.fixture(scope='module')
def member():
    return v.contract.member(1)


def test_all_word_casimir_floor_includes_characters_and_nonabelian_labels():
    # This is a conservative superset, without the other descent congruences.
    values=[]
    for p in range(4):
        for q in range(4):
            for j in range(5):
                for z in range(-12,13):
                    if p==q==j==0 and z%6:continue
                    c=v.casimir(p,q,j,z)
                    if (p,q,j,z)!=(0,0,0,0):values.append(c)
    assert min(values)==F(3,10)
    assert v.casimir(1,0,0,0)==F(2,3)
    assert v.casimir(0,0,1,0)==F(3,8)
    with pytest.raises(ValueError):v.casimir(0,0,0,1)
    with pytest.raises(ValueError):v.casimir(True,0,0,6)
    actual=v.vacuum_channels()
    assert min(F(r['trace_metric_Casimir']) for r in actual['all_native_blocks'])==F(3,10)
    assert actual['dimensionless_decoupled_physical_gap']==6


def _car_bond_in_negative_beta_frame(axis):
    beta=v.density.gamma_data()[3]
    columns=[];signs=[]
    for sign in (-1,1):
        for vector in (beta-sign*s.eye(4)).nullspace():
            columns.append(vector/s.sqrt((vector.H*vector)[0]));signs.append(sign)
    frame=s.Matrix.hstack(*columns)
    assert frame.H*frame==s.eye(4)
    gamma=v.density.gamma_data()[axis]
    hop=v.density.clean(frame.H*((-beta+beta*gamma)/2)*frame)
    one=s.zeros(4).row_join(hop).col_join(hop.H.row_join(s.zeros(4)))
    energies=[3*z for z in signs]*2
    masks=[sum(1<<i for i in chosen) for chosen in combinations(range(8),4)]
    index={mask:i for i,mask in enumerate(masks)}
    matrix=s.zeros(len(masks))
    def apply(mask,i,j):
        if not mask&(1<<j):return None
        sign=(-1)**((mask&((1<<j)-1)).bit_count())
        mask^=1<<j
        if mask&(1<<i):return None
        sign*=(-1)**((mask&((1<<i)-1)).bit_count())
        return mask|(1<<i),sign
    for col,mask in enumerate(masks):
        for (i,j),coefficient in one.todok().items():
            image=apply(mask,i,j)
            if image is not None:matrix[index[image[0]],col]+=image[1]*coefficient
    vacuum=sum(1<<i for i,e in enumerate(energies) if e<0)
    vector=matrix[:,index[vacuum]]
    diagonal=s.diag(*[sum(e for i,e in enumerate(energies) if mask&(1<<i)) for mask in masks])
    return matrix,diagonal,vector,index[vacuum]


@pytest.mark.parametrize('axis',[0,1,2])
def test_independent_full_number_sector_car_action_and_perturbative_energy(axis):
    hopping,onsite,vector,vacuum=_car_bond_in_negative_beta_frame(axis)
    assert hopping==hopping.H
    assert (vector.H*vector)[0]==1
    assert (onsite+12*s.eye(70))*vector==6*vector
    assert vector[vacuum]==0
    # Independent exact 70-state CAR matrix; numeric eigenvalues only check
    # the finite lambda expansion, not a primary spatial spectrum.
    h0=np.array(onsite,dtype=complex);pert=np.array(hopping,dtype=complex)
    step=1e-2
    shifted=np.linalg.eigvalsh(h0+step*pert)[0]+12
    assert abs(shifted/step**2+1/6)<1e-6


def test_native_packets_and_literal_open_box_counts(member):
    local=v.vacuum_channels();result=v.selected_expansion(member)
    assert {r['integer_Q6_charge']:F(r['bond_H0_excitation']) for r in local['all_native_blocks']}=={
        -4:F(48),-3:F(237,8),1:F(489,8),2:F(171,4),6:F(87,4)}
    assert F(local['bond_vacuum_variance'])==15
    assert F(local['plaquette_vacuum_variance'])==F(1,4410)
    assert F(local['bond_inverse_energy_weight'])<F(7,20)
    assert F(local['bond_inverse_square_weight'])<F(1,100)
    source=member['source'];side=source['side']
    bonds=sum(1 for x in range(side) for y in range(side) for z in range(side)
              for axis in range(3) if (x,y,z)[axis]+1<side)
    assert result['actual_interior_bonds']==bonds==source['interior_hopping_bonds']
    b=F(local['bond_inverse_energy_weight']);p=F(local['plaquette_inverse_energy_weight'])
    assert F(result['second_order_energy_coefficient'])==-b*bonds-p*source['interior_plaquettes']
    assert result['second_order_remainder_bound_at_primary_endpoint'] is None
    assert not result['anchor_state_is_selected_primary_preparation']
    assert not result['uniform_local_gap_stability_proved']


def test_volume_growth_prevents_global_dressing_inference():
    small=v.selected_expansion(v.contract.member(0))
    large=v.selected_expansion(v.contract.member(2))
    assert F(large['first_order_orthogonal_dressing_norm_squared'])>100*F(small['first_order_orthogonal_dressing_norm_squared'])
    assert F(large['strict_lambda_threshold_for_positive_global_Weyl_bound'])<F(small['strict_lambda_threshold_for_positive_global_Weyl_bound'])
    assert F(small['global_Weyl_gap_lower_at_primary_endpoint'])<0
    assert F(large['global_Weyl_gap_lower_at_primary_endpoint'])<0
    assert F(large['first_order_dressing_norm_squared_per_cell'])<F(3,100)


def test_selected_clock_window_and_conditional_gap_are_not_endpoint_proofs(member):
    a=F(member['source']['spacing'])
    result=v.fixed_band_scale(member)
    assert list(map(F,result['dimensionless_band']))==[a/4,4*a]
    assert F(result['dimensionless_selected_mean_ground_excess_upper'])==a*a/16
    control=v.conditional_gap_exclusion(member,'3')
    assert control['band_strictly_below_supplied_gap']
    assert F(control['conditional_original_band_weight_upper'])==a*a/48
    assert control['gap_is_external_unproved_premise_for_primary']
    assert not control['actual_primary_gap_or_collapse_proved']
    assert control['clock_id']==member['ids']['clock_id']
    for field in ('clock','state','reference'):
        altered=deepcopy(member);altered[field]['unregistered_change']=True
        with pytest.raises(ValueError):v.selected_expansion(altered)


def test_gap_endpoint_and_missing_gap_cannot_be_admitted():
    member=v.contract.member(0)
    endpoint=v.conditional_gap_exclusion(member,'4')
    assert not endpoint['band_strictly_below_supplied_gap']
    assert endpoint['conditional_original_band_weight_upper'] is None
    # H=diag(0,4), exact ground and exchange A: positive endpoint weight one.
    h=np.diag([0.,4.]);a=np.array([[0.,1.],[1.,0.]]);rho=np.diag([1.,0.])
    assert np.trace(rho@a@np.diag([0.,1.])@a)==1
    for invalid in (0,-1,True,0.5):
        with pytest.raises(ValueError):v.conditional_gap_exclusion(member,invalid)


def test_ground_degeneracy_and_original_excited_transition_weight():
    energies=(F(0),F(0),F(10),F(11))
    probabilities=(F(3,8),F(3,8),F(1,8),F(1,8))
    a=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])/2
    a[2,3]=a[3,2]=s.Rational(1,2)
    assert max(abs(float(x)) for x in a.eigenvals())<1
    original=sum(p*F(abs(a[j,i])**2) for i,p in enumerate(probabilities)
                 for j in range(4) if 0<energies[j]-energies[i]<10)
    excess=sum(p*e for p,e in zip(probabilities,energies))
    assert original==F(1,32)<=excess/10


def test_same_polynomial_shape_does_not_restore_postselected_weight():
    weights=[]
    for level in (0,2,5,8):
        result=v.polynomial_countercontrol(level)
        p=list(map(F,result['Q_squared_probabilities']))
        assert sum(p)==1 and p[0]>p[1]>p[2]>0
        assert F(result['mean_ground_excess'])<=F(result['cooling_tolerance'])
        original=F(result['original_positive_fixed_band_weight']);weights.append(original)
        assert 0<original<=F(result['conditional_gap_original_weight_upper'])
        assert F(result['postselected_positive_fraction'])>F(99,100)
        assert not result['is_selected_primary']
        assert not result['postselected_fraction_is_original_weight']
    assert all(left>right for left,right in zip(weights,weights[1:]))


def test_every_owned_dependency_exists():
    assert all((v.contract.ROOT/name).is_file() for name in v.SOURCE_INPUTS)
