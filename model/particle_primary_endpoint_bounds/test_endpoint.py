"""Full CAR and quantum-link controls for the actual endpoint inequality."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
import sympy as s

from . import endpoint as e


@pytest.fixture(scope='module')
def member():return e.contract.member(0)


def _number_sector(modes,particles):
    masks=[sum(1<<i for i in subset) for subset in combinations(range(modes),particles)]
    locate={mask:i for i,mask in enumerate(masks)}
    def monomial(i,j):
        out=np.zeros((len(masks),len(masks)),dtype=complex)
        for col,mask in enumerate(masks):
            if not mask&(1<<j):continue
            sign=(-1)**((mask&((1<<j)-1)).bit_count());reduced=mask^(1<<j)
            if reduced&(1<<i):continue
            sign*=(-1)**((reduced&((1<<i)-1)).bit_count())
            out[locate[reduced|(1<<i)],col]=sign
        return out
    return masks,monomial


@pytest.mark.parametrize('dimension',[1,2])
def test_independent_finite_unitary_design_and_full_occupation_haar_cap(dimension):
    d=dimension;masks,car=_number_sector(4*d,2*d)
    size=len(masks)
    a={(i,j):sum((car(flavor*d+i,2*d+flavor*d+j) for flavor in range(2)),np.zeros((size,size),complex))
       for i in range(d) for j in range(d)}
    if d==1:design=[np.ones((1,1),complex)]
    else:design=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    design=[phase*u for u in design for phase in (1,1j,-1,-1j)]
    average=np.zeros((size,size),complex)
    for u in design:
        hop=sum((u[i,j]*a[i,j] for i in range(d) for j in range(d)),np.zeros((size,size),complex))
        full=hop+hop.conj().T
        average+=full@full/len(design)
    def current(side,left,right):
        return sum((car(side*2*d+left*d+i,side*2*d+right*d+i) for i in range(d)),np.zeros((size,size),complex))
    ex={(a,b):current(0,a,b) for a in range(2) for b in range(2)}
    ey={(a,b):current(1,a,b) for a in range(2) for b in range(2)}
    nx=ex[0,0]+ex[1,1];ny=ey[0,0]+ey[1,1]
    formula=nx+ny-2*sum((ex[a,b]@ey[b,a] for a in range(2) for b in range(2)),np.zeros((size,size),complex))/d
    assert np.max(np.abs(average-formula))<1e-13
    assert abs(np.linalg.eigvalsh(average)[-1]-2*(d+1))<1e-12
    assert np.linalg.eigvalsh(average)[0]>-1e-12


def test_all_native_occupation_spin_maxima_are_inside_the_exact_cap():
    for d in (1,2,3,6):
        values=[e.haar_variance_max(d,x,y) for x in range(2*d+1) for y in range(2*d+1)]
        assert min(values)>=0 and max(values)==2*(d+1)
    with pytest.raises(ValueError):e.haar_variance_max(2,5,1)


def test_actual_eight_coordinate_polynomial_and_full_fock_minimum():
    result=e.actual_spin_polynomial()
    assert len(result['one_internal_column_checks'])==3
    for m in (F(1,2),F(1),F(3,2)):
        beta=e.density.gamma_data()[3]
        for gamma in e.density.gamma_data()[:3]:
            hop=(-beta+beta*gamma)/2
            h=(s.Rational(m)*beta).row_join(hop).col_join(hop.H.row_join(s.Rational(m)*beta))
            values=np.linalg.eigvalsh(np.array(h,dtype=complex))
            fock_floor=sum(value for value in values if value<0)+4*float(m)
            expected=-4*(np.sqrt(float(m*m+F(1,4)))-float(m))
            assert abs(fock_floor-expected)<1e-12


def test_quantum_link_compression_and_full_car_cluster_satisfy_schur_lower():
    # U(1) one-internal-column control; it is not the full primary source.
    masks,car=_number_sector(8,4);size=len(masks)
    beta=e.density.gamma_data()[3];gamma=e.density.gamma_data()[0]
    hop=np.array((-beta+beta*gamma)/2,dtype=complex)
    onsite=np.array(s.diag(beta,beta),dtype=complex)
    fermion=sum((onsite[i,j]*car(i,j) for i in range(8) for j in range(8)),np.zeros((size,size),complex))
    forward=sum((hop[i,j]*car(i,4+j) for i in range(4) for j in range(4)),np.zeros((size,size),complex))
    link=np.array([[0,0,0],[1,0,0],[0,1,0]],complex)
    electric=np.diag([63/4,0,63/4])
    h=np.kron(np.eye(3),fermion+4*np.eye(size))+np.kron(electric,np.eye(size))
    h+=np.kron(link,forward)+np.kron(link.conj().T,forward.conj().T)
    delta=4*(np.sqrt(1.25)-1)
    floor=-4/(63/4-delta)
    assert np.linalg.eigvalsh(h)[0]>=floor-1e-12
    # The coupling is present and produces a strictly negative endpoint level.
    assert np.linalg.eigvalsh(h)[0]<0
    p=np.kron(np.diag([0,1,0]),np.eye(size))
    q=np.eye(3*size)-p
    coupling=np.kron(link,forward)+np.kron(link.conj().T,forward.conj().T)
    assert np.max(np.abs(p@coupling@p))==0
    assert np.linalg.norm(p@coupling@q)>0


def test_actual_endpoint_bound_improves_source_floor_and_retains_complete_domain(member):
    value=e.bind_member(member);e.verify_bound(member,value)
    result=value['result'];lower=F(result['physical_ground_lower'])
    assert F(-741)<lower<F(-739)
    assert lower>F(result['inherited_full_source_lower'])
    assert F(result['physical_ground_upper_from_product_trial'])==-720
    assert result['full_physical_carrier'] and result['primary_endpoint_lambda']=='1'
    assert not result['selected_polynomial_state_changed']
    assert not result['source_reference_or_clock_changed']
    assert result['physical_mass_table_GeV']==[]
    changed=deepcopy(value);changed['result']['physical_ground_lower']='-723'
    with pytest.raises(ValueError):e.verify_bound(member,changed)
    changed_member=deepcopy(member);changed_member['clock']['fitted_scale']='1 GeV'
    with pytest.raises(ValueError):e.bind_member(changed_member)


def test_large_box_allocation_and_failed_schur_admission_are_explicit():
    member=e.contract.member(1)
    data=e.local_schur_data(member)
    assert data['maximum_interior_vertex_degree']==6
    assert data['onsite_coefficient_per_endpoint']==F(1,2)
    assert data['nonconstant_link_block_positive']
    failed=e.retained_penalty_bound(member,'1/2')
    assert not failed['local_bound']['nonconstant_link_block_positive']
    assert failed['physical_ground_lower'] is None
    assert failed['physical_second_eigenvalue_lower'] is None
    with pytest.raises(ValueError):e.local_schur_data(member,'1')


def test_retained_penalty_monotonicity_does_not_imply_global_gap(member):
    base=e.local_schur_data(member)
    for r in ('1/10','1/4'):
        value=e.local_schur_data(member,r)
        assert value['nonconstant_link_block_positive']
        assert value['per_bond_Schur_loss_upper']>base['per_bond_Schur_loss_upper']
        assert not e.retained_penalty_bound(member,r)['branch_isolation_proved']


def test_ground_floor_and_method_ceiling_use_opposite_rounding_directions(member):
    local=e.local_schur_data(member)
    m=local['onsite_coefficient_per_endpoint']
    lo=local['full_Fock_pointwise_lower_shift_lower']/60+m
    hi=local['full_Fock_pointwise_lower_shift_upper']/60+m
    assert lo*lo<=m*m+F(1,4)<=hi*hi
    assert local['per_bond_exact_Schur_loss_lower']<=local['per_bond_Schur_loss_upper']
    value=e.selected_bound(member)
    assert F(value['upper_ceiling_on_this_retained_penalty_method_second_eigenvalue_floors'])>=F(value['physical_ground_lower'])+6


def test_authenticated_variational_join_quantifies_failed_isolation(member):
    from model.particle_primary_variational_readout import variational
    upper=variational.bind_member(member)
    value=e.compare_with_variational(member,upper)
    lower,hi=map(F,value['full_physical_ground_energy_interval'])
    assert lower<hi<-723
    assert F(value['ceiling_minus_actual_variational_upper'])<0
    assert value['entire_retained_penalty_family_fails_with_this_upper']
    assert not value['actual_full_source_gap_ruled_out']
    changed=deepcopy(upper);changed['result']['ground_energy_upper']='-739'
    with pytest.raises(ValueError):e.compare_with_variational(member,changed)


def test_all_declared_source_paths_exist():
    assert all((e.contract.ROOT/name).is_file() for name in e.SOURCE_INPUTS)
