"""Independent CAR controls and authenticated actual-source bounds."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse

from . import model as o


def _bilinear(modes,i,j):
    rows=[];cols=[];values=[]
    for mask in range(1<<modes):
        if not mask&(1<<j):continue
        sign=(-1)**((mask&((1<<j)-1)).bit_count());reduced=mask^(1<<j)
        if reduced&(1<<i):continue
        sign*=(-1)**((reduced&((1<<i)-1)).bit_count())
        rows.append(reduced|(1<<i));cols.append(mask);values.append(sign)
    return sparse.csr_matrix((values,(rows,cols)),shape=(1<<modes,1<<modes),dtype=float)


def _site(d):
    modes=4*d;size=1<<modes;identity=sparse.eye(size,format='csr')
    currents={(a,b):sum((_bilinear(modes,a*d+i,b*d+i) for i in range(d)),sparse.csr_matrix((size,size)))
              for a in range(2) for b in range(2)}
    beta=sum((_bilinear(modes,j,2*d+j)+_bilinear(modes,2*d+j,j) for j in range(2*d)),sparse.csr_matrix((size,size)))
    defects=beta+2*d*identity
    vacuum={0:1.0}
    for j in range(2*d):
        changed={}
        for mask,value in vacuum.items():
            for mode,coefficient in ((j,1/np.sqrt(2)),(2*d+j,-1/np.sqrt(2))):
                if mask&(1<<mode):continue
                sign=(-1)**((mask&((1<<mode)-1)).bit_count())
                changed[mask|(1<<mode)]=changed.get(mask|(1<<mode),0)+sign*coefficient*value
        vacuum=changed
    vector=np.zeros(size)
    for mask,value in vacuum.items():vector[mask]=value
    total=sum((_bilinear(modes,j,j) for j in range(modes)),sparse.csr_matrix((size,size)))
    return identity,currents,defects,vector,total


def _two_sites(d):
    identity,e,D,omega,total=_site(d)
    number=e[0,0]+e[1,1]
    quartic=sparse.kron(number,identity)+sparse.kron(identity,number)
    quartic-=F(2,d).numerator/F(2,d).denominator*sum(
        (sparse.kron(e[a,b],e[b,a]) for a in range(2) for b in range(2)),
        sparse.csr_matrix((identity.shape[0]**2,identity.shape[0]**2)))
    defects=sparse.kron(D,identity)+sparse.kron(identity,D)
    return quartic.tocsr(),defects.tocsr(),np.kron(omega,omega),sparse.kron(total,identity).tocsr()


@pytest.mark.parametrize('d',[1,2])
def test_full_active_inactive_car_vacuum_variance_and_exact_four_defects(d):
    Q,D,omega,total_x=_two_sites(d)
    assert abs(omega@omega-1)<1e-13
    assert np.linalg.norm(D@omega)<1e-13
    image=Q@omega-d*omega
    assert abs(omega@(Q@omega)-d)<1e-13
    assert abs(image@image-1)<1e-13
    assert np.linalg.norm(D@image-4*image)<1e-12
    assert np.linalg.norm(Q@omega-d*omega)>.9
    assert sparse.linalg.norm(Q@total_x-total_x@Q)<1e-13
    # The quartic mixes beta defects; its expectation alone is no bound.
    assert sparse.linalg.norm(Q@D-D@Q)>1


@pytest.mark.parametrize('m,w',[(F(1),F(1,8)),(F(3,5),F(1,5)),(F(1,5),F(2))])
def test_full_fock_effective_spectrum_obeys_occupation_bound(m,w):
    Q,D,_,_=_two_sites(1)
    result=o.occupation_loss(1,m,w)
    eig=np.linalg.eigvalsh((float(m)*D-float(w)*Q).toarray())
    assert eig[0]>=-float(result['occupation_loss_upper'])-1e-12
    assert result['active_and_inactive_Dirac_rows_retained']
    assert eig[0]<-float(w)  # The omitted vacuum crossblock is measurable.


def test_quantum_link_and_all_eight_dirac_fock_modes_obey_improved_schur_bound():
    # A U(1) one-column control executes the coupled link and full Fock space.
    # It tests the operator inequality, without claiming a spatial source run.
    size=256;beta=np.array(o.endpoint.density.gamma_data()[3],dtype=complex)
    gamma=np.array(o.endpoint.density.gamma_data()[0],dtype=complex)
    spin=(-beta+beta@gamma)/2
    matter=sum((beta[i,j]*(_bilinear(8,i,j)+_bilinear(8,4+i,4+j))
                for i in range(4) for j in range(4)),sparse.csr_matrix((size,size),dtype=complex))
    forward=sum((spin[i,j]*_bilinear(8,i,4+j) for i in range(4) for j in range(4)),
                sparse.csr_matrix((size,size),dtype=complex))
    link=sparse.csr_matrix(np.array([[0,0,0],[1,0,0],[0,1,0]]))
    h=sparse.kron(sparse.eye(3),matter+4*sparse.eye(size))
    h+=sparse.kron(sparse.diags([63/4,0,63/4]),sparse.eye(size))
    h+=sparse.kron(link,forward)+sparse.kron(link.T,forward.conj().T)
    delta=4*(o.sqrt_upper(F(5,4))-1)
    w=1/(F(63,4)-delta)
    c=o.occupation_loss(1,F(1),w)['occupation_loss_upper']
    smallest=sparse.linalg.eigsh(h,k=1,which='SA',return_eigenvectors=False,tol=1e-12)[0]
    assert -float(c)-1e-11<=smallest<0
    assert c<4*w


@pytest.fixture(scope='module')
def member():return o.contract.member(0)


def test_actual_endpoint_floor_and_selected_defect_weight(member):
    receipt=o.bind_member(member);o.verify_bound(member,receipt)
    value=receipt['result']
    assert F('-727.499')<F(value['zero_retention_physical_ground_lower'])<F('-727.498')
    assert F(value['zero_retention_physical_ground_lower'])>F(value['inherited_endpoint_physical_ground_lower'])
    assert F(value['zero_retention_ground_enclosure_width'])<F('3.731')
    assert F('4.99')<F(value['defect_count_upper'])<5
    assert F(value['average_matched_q6_pair_weight_lower'])>F(11,16)
    assert F(value['average_matched_q6_pair_weight_simple_lower'])==F(11,16)
    assert value['matched_q6_pair_count']==16
    assert value['weight_is_uniform_address_average']
    assert not value['weight_is_conditioned_on_success']
    assert not value['each_individual_site_weight_claimed']
    assert not value['local_vacuum_is_used_as_selected_preparation']
    assert value['physical_mass_table_GeV']==[]


def test_pair_weight_operator_inequality_all_occupations():
    # Two beta+ and two beta- modes at a q6 singlet site, all 16 occupations.
    n=[_bilinear(4,j,j).toarray() for j in range(4)]
    identity=np.eye(16)
    weights=n[2]@(identity-n[0])+n[3]@(identity-n[1])
    D=(identity-n[2])+(identity-n[3])+n[0]+n[1]
    assert np.linalg.eigvalsh(weights-(2*identity-D))[0]>=0
    assert np.linalg.eigvalsh(D)[-1]==4


def test_independent_allocations_and_larger_member_failure_are_explicit(member):
    onsite=o.retained_bound(member,'2/5','0')
    assert F(onsite['local_bound']['retained_dimensionless_defect_coefficient'])==F(6,5)
    assert F(onsite['local_bound']['retained_H0_fraction_lower'])==0
    assert not onsite['endpoint_gap_or_particle_mass_proved']
    equal=o.retained_bound(member,'1/10','1/10')
    assert F(equal['physical_ordered_second_eigenvalue_lower'])-F(equal['physical_ground_lower'])==F(3,5)
    failed=o.retained_bound(o.contract.member(1),'2/5','0')
    assert not failed['local_bound']['nonconstant_link_block_positive']
    assert failed['physical_ground_lower'] is None
    for bad in ('1',-1,True,0.1):
        with pytest.raises(ValueError):o.cluster_data(member,bad)


def test_rounding_is_certified_by_exact_schur_inequalities(member):
    for r in ('0','1/5','2/5','1/2'):
        data=o.cluster_data(member,r)
        m=F(data['onsite_per_bond_endpoint'])
        for row in data['native_blocks']:
            c=row['occupation_loss_upper'];w=row['inverse_electric_margin'];d=row['dimension']
            assert (c-w*d)*row['vacuum_resolvent_denominator_lower']>=w*w
            assert c>=w*2*(d+1)-m
            assert row['balanced_complement_shifted_margin_lower']>0


def test_receipts_reject_source_clock_and_bound_mutations(member):
    receipt=o.bind_member(member)
    changed=deepcopy(receipt);changed['result']['defect_count_upper']='0'
    with pytest.raises(ValueError):o.verify_bound(member,changed)
    altered=deepcopy(member);altered['clock']['fitted_scale']='GeV'
    with pytest.raises(ValueError):o.bind_member(altered)
    with pytest.raises(ValueError):o.verify_bound(o.contract.member(1),receipt)


def test_all_evidence_paths_exist():
    assert all((o.contract.ROOT/name).is_file() for name in o.SOURCE_INPUTS)
