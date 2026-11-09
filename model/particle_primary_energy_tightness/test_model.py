"""Independent source columns, full-Fock ladders and retained-origin tests."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
from math import comb

import numpy as np
import pytest
from scipy import sparse

from . import model as m


def _bilinear(modes,i,j):
    rows=[];columns=[];values=[]
    for mask in range(1<<modes):
        if not mask&(1<<j):continue
        sign=(-1)**((mask&((1<<j)-1)).bit_count());state=mask^(1<<j)
        if state&(1<<i):continue
        sign*=(-1)**((state&((1<<i)-1)).bit_count())
        rows.append(state|(1<<i));columns.append(mask);values.append(sign)
    return sparse.csr_matrix((values,(rows,columns)),shape=(1<<modes,1<<modes),dtype=complex)


def _lift(one_body):
    modes=one_body.shape[0];result=sparse.csr_matrix((1<<modes,1<<modes),dtype=complex)
    for i,j in zip(*np.nonzero(one_body)):
        result+=one_body[i,j]*_bilinear(modes,int(i),int(j))
    return result


def test_actual_source_interior_gram_and_all_boundary_columns():
    value=m.source_column_control()
    assert len(value['all_boundary_patterns'])==27
    assert {row['degree'] for row in value['all_boundary_patterns']}=={3,4,5,6}
    assert sum(row['interior'] for row in value['all_boundary_patterns'])==1
    assert value['interior_native6_dimensionless_residual_norm_squared']==F(78,7)
    assert not value['link_coordinates_replaced_in_selected_source']
    assert not value['physical_anchor_is_selected_state']


def test_full_fock_interior_ladder_norm_is_larger_than_one_body_norm():
    # d=1: four disjoint raising pairs, all eight CAR modes retained.
    ladder=sum((_bilinear(8,2*j,2*j+1) for j in range(4)),sparse.csr_matrix((256,256),dtype=complex))
    residual=np.sqrt(1.5)*ladder
    square=(residual.conj().T@residual).toarray()
    assert abs(np.linalg.eigvalsh(square)[-1]-9)<1e-12
    # In the one-fermion sector the same squared norm is only3/2.
    rows=[i for i in range(256) if i.bit_count()==1]
    assert abs(np.linalg.eigvalsh(square[np.ix_(rows,rows)])[-1]-1.5)<1e-12


@pytest.mark.parametrize('native_dimension',[1,2,6])
def test_complete_collective_ladder_multiplicities_and_exact_maximum(native_dimension):
    pairs=4*native_dimension;dimension=0;maximum=F()
    for active in range(pairs+1):
        copies=comb(pairs,active)*2**(pairs-active)
        for q in range(active//2+1):
            twice_j=active-2*q
            mult=comb(active,q)-(comb(active,q-1) if q else 0)
            dimension+=copies*mult*(twice_j+1)
            j=F(twice_j,2)
            for twice_z in range(-twice_j,twice_j+1,2):
                z=F(twice_z,2)
                maximum=max(maximum,j*(j+1)-z*(z+1))
    assert dimension==4**pairs
    assert maximum==2*native_dimension*(2*native_dimension+1)
    normalized=3*maximum/(native_dimension*(native_dimension+1))
    assert normalized==F(6*(2*native_dimension+1),native_dimension+1)


def test_literal_car_and_quantum_link_anchor_residual_retains_word_degree_one():
    gamma=m.density.gamma_data();beta=np.array(gamma[3],dtype=complex)
    S=np.array((-gamma[3]+gamma[3]*gamma[0])/2,dtype=complex)
    pairs=[(np.array(p,dtype=complex).ravel(),np.array(q,dtype=complex).ravel())
           for p,q in m.response.beta_pairs()]
    B=sum((np.outer(p,q.conjugate()) for p,q in pairs),np.zeros((4,4),dtype=complex))
    central=np.zeros((8,8),dtype=complex);central[:4,:4]=B
    forward=np.zeros_like(central);forward[:4,4:]=S
    onsite=np.zeros_like(central);onsite[:4,:4]=3*beta;onsite[4:,4:]=3*beta
    b=_lift(central);hop=_lift(forward);h0=_lift(onsite)
    phi=np.zeros((8,4),dtype=complex)
    for j,(_,negative) in enumerate(pairs):
        phi[:4,j]=negative;phi[4:,j+2]=negative
    omega=np.zeros(256,dtype=complex)
    for occupied in combinations(range(8),4):
        mask=sum(1<<i for i in occupied);omega[mask]=np.linalg.det(phi[list(occupied),:])
    assert abs(np.vdot(omega,omega)-1)<1e-13
    shift=sparse.csr_matrix([[0,0,0],[1,0,0],[0,1,0]],dtype=complex)
    electric=sparse.diags([15.75,0,15.75],dtype=complex,format='csr')
    H=sparse.kron(sparse.eye(3),h0)+sparse.kron(electric,sparse.eye(256))
    H+=sparse.kron(shift,hop)+sparse.kron(shift.conj().T,hop.conj().T)
    probe=sparse.kron(sparse.eye(3),b)
    anchor=np.kron(np.array([0,1,0]),omega)
    output=probe@anchor;R=H@probe-probe@H-6*probe;residual=R@anchor
    assert abs(np.vdot(output,output)-2)<1e-12
    assert abs(np.vdot(residual,residual)-1)<1e-12
    assert abs(np.vdot(output,residual))<1e-13
    assert abs(np.vdot((H@probe-probe@H)@anchor,(H@probe-probe@H)@anchor)-73)<1e-11
    # The nonzero constant is independent of regulator spacing; H/a
    # multiplies the actual physical residual squared by1/a².
    for a in (1,F(1,2),F(1,8)):
        assert abs(np.vdot(residual/float(a),residual/float(a))-1/float(a*a))<1e-10


@pytest.mark.parametrize('side',[2,3,4,8])
def test_actual_site_sum_matches_the_boundary_and_interior_formula(side):
    total=F()
    for x in np.ndindex((side,side,side)):
        degree=sum(int(k>0)+int(k+1<side) for k in x)
        interior=all(0<k<side-1 for k in x)
        total+=F(78,7) if interior else F(24*degree,7)
    assert total/side**3==m.moment_coefficient(side,6)
    assert m.moment_coefficient(side,1)>0


def test_selected_sides_have_a_sharp_uniform_constant_without_a_continuum_inference():
    values=[m.moment_coefficient(2**(2*n+1),6) for n in range(30)]
    assert max(values)==values[1]==F(3141,224)
    assert values[0]==F(72,7)
    assert all(values[n+1]<values[n] for n in range(1,len(values)-1))
    assert values[-1]>F(78,7)


def test_countercontrol_satisfies_population_cooling_and_moment_bounds_but_escapes():
    for n in (0,2,5,9):
        a=F(1,2**n);E=6/a;Z=(1+4*E)**2+1
        excited=1/Z;ground=1-excited
        assert E*excited<a/16
        assert ground-excited>F(1,20)
        assert (E-6/a)**2*ground==0
        assert F(1,800)<E<F(46)/a
        assert ground>F(3,100)
        # All original weight lies at E=6/a. No fixed[1,10] response
        # survives this independent positive-polynomial sequence.
        if n>=2:assert E>10


def test_direct_origin_defect_estimate_is_not_a_probability_trimming_rule():
    fraction=F(15,16)**3;density=F(9658,1000)
    chi=(12-density/fraction)/42
    assert chi>F(1,200)
    # The old moving-band weight loses its proof when this volume is
    # simply subtracted; the direct commutator estimate remains positive.
    assert (F(3,100)-(1-fraction))/fraction<0


@pytest.fixture(scope='module')
def receipts():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,7)}


def test_actual_narrowed_source_bands_and_level_one_q6(receipts):
    for n in (0,1,2):
        result=receipts[n]['result'];window=result['native6_original_window']
        assert F(window['original_transfer_weight_lower'])>F(3,100)
        assert F(window['original_transfer_interval'][1])==46*2**n
        assert F(result['finite_reported_original_weight_lower'])>F(1,50)
        assert not result['fixed_energy_spectral_tightness_proved']
        assert not result['source_or_preparation_or_clock_changed']
        assert result['physical_mass_table_GeV']==[]
    q6=receipts[1]['result']['q6_scalar_original_window']
    assert F(q6['original_transfer_interval'][1])==48
    assert F(q6['original_transfer_weight_lower'])>F(1,20)


def test_actual_retained_origin_band_uses_the_same_member_and_positive_interior_weight(receipts):
    origin=receipts[7]['result']['retained_origin_response']
    assert F(origin['origin_geometry']['retained_fraction'])==F(15,16)**3
    assert F(origin['direct_selected_commutator_lower'])>F(7,1000)
    assert F(origin['original_window']['original_transfer_weight_lower'])>F(35,10000)
    assert origin['same_origin_register_for_complete_history']
    assert not origin['averaged_observable_substituted_for_origin_history']
    assert not origin['fixed_energy_uniform_weight_proved']
    assert receipts[1]['result']['retained_origin_response']['original_window'] is None


def test_physical_anchor_obstruction_is_not_a_selected_state_lower_bound(receipts):
    for n,receipt in receipts.items():
        result=receipt['result']
        assert F(result['physical_anchor_average_centered_residual_squared'])>0
        assert not result['physical_anchor_is_selected_state']
        assert not result['selected_state_weighted_residual_lower_claimed']
        assert not result['source_operator_norm_method_can_remove_inverse_spacing']


def test_actual_receipts_reject_source_or_result_mutation(receipts):
    member=m.contract.member(1);m.verify_bound(member,receipts[1])
    wrong=deepcopy(receipts[1]);wrong['result']['fixed_energy_spectral_tightness_proved']=True
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    wrong=deepcopy(member);wrong['state']['target_ground_excess']='0'
    with pytest.raises(ValueError):m.bind_member(wrong)
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(0),receipts[1])


def test_evidence_paths_exist():
    assert all((m.contract.ROOT/name).is_file() for name in m.SOURCE_INPUTS)


def test_source_number_cube_and_boundary_controls_are_actual():
    value=m.density_commutator_control()
    assert value['actual_Wilson_matrix_and_quantum_link_checked']
    assert value['source_rotation_control']['proper_spatial_rotations']==24
    assert value['physical_boundary_control']['interior_graph_cube_covariant']
    assert not value['ambient_dangling_link_permutation_used']
    assert sum(row['CAR_modes_per_cell'] for row in value['complete_native_blocks'])==60


@pytest.fixture(scope='module')
def density_link_control():
    # One actual four-Dirac-column block at both endpoints. Full Fock,
    # electric coordinates and omitted end-shift transitions are retained.
    g=m.density.gamma_data();beta=np.array(g[3],dtype=complex)
    S=np.array((-g[3]+g[3]*g[0])/2,dtype=complex)
    forward=np.zeros((8,8),complex);forward[:4,4:]=S
    onsite=np.zeros((8,8),complex);onsite[:4,:4]=3*beta;onsite[4:,4:]=3*beta
    hop=_lift(forward);h0=_lift(onsite)
    shift=sparse.csr_matrix([[0,0,0],[1,0,0],[0,1,0]],dtype=complex)
    J=sparse.kron(shift,hop)+sparse.kron(shift.conj().T,hop.conj().T)
    electric=sparse.diags([15.75,0,15.75],format='csr')
    H=J+sparse.kron(sparse.eye(3),h0)+sparse.kron(electric,sparse.eye(256))
    nx=np.array([(mask&15).bit_count() for mask in range(256)])
    ny=np.array([(mask>>4).bit_count() for mask in range(256)])
    X=sparse.diags(np.tile(nx-2,3),format='csr',dtype=complex)
    Y=sparse.diags(np.tile(ny-2,3),format='csr',dtype=complex)
    return H.tocsr(),J.tocsr(),X,Y


def test_full_car_quantum_link_double_commutator_and_gauss_density(density_link_control):
    H,J,X,Y=density_link_control
    def comm(A,B):return A@B-B@A
    error=comm(X,comm(H,X))+comm(Y,comm(H,Y))+2*J
    assert np.max(abs(error.data),initial=0)<1e-14
    assert np.max(abs(comm(X+Y,H).data),initial=0)<1e-14
    link_flux=sparse.kron(sparse.diags([-1,0,1],dtype=complex),sparse.eye(256))
    for gauss in (X-link_flux,Y+link_flux):
        assert np.max(abs(comm(H,gauss).data),initial=0)<1e-14
        assert np.max(abs(comm(X,gauss).data),initial=0)<1e-14
        assert np.max(abs(comm(Y,gauss).data),initial=0)<1e-14
    # Positive-β occupation is a different observable; retaining only
    # those modes does not inherit the all-mode number f-sum.
    beta=np.array(m.density.gamma_data()[3],dtype=complex)
    wrong=np.zeros((8,8),complex);wrong[:4,:4]=(np.eye(4)+beta)/2
    Z=sparse.kron(sparse.eye(3),_lift(wrong))
    incorrect=comm(Z,comm(H,Z))+J
    assert np.max(abs(incorrect.data),initial=0)>.1


def test_complete_density_instrument_original_fsum_and_retained_reference(density_link_control):
    H,J,X,Y=density_link_control;dimension=H.shape[0];normalization=4
    Kc=X/normalization;Ks=Y/normalization
    diagonal=1-Kc.diagonal().real**2-Ks.diagonal().real**2
    assert diagonal.min()>=0
    Kf=sparse.diags(np.sqrt(diagonal),format='csr')
    Q=64*sparse.eye(dimension)-H
    rho=(Q@Q).toarray();rho/=np.trace(rho)
    moment=0j
    for K in (Kc,Ks):
        moment+=np.trace(rho@(K@(H@K-K@H)).toarray())
    expected=-np.trace(rho@J.toarray())/normalization**2
    assert abs(moment-expected)<1e-13
    assert moment.real>0
    # Randomly selecting one completed setting has half this moment.
    setting_moment=sum(np.trace(rho@((K/np.sqrt(2))@(H@(K/np.sqrt(2))-
                       (K/np.sqrt(2))@H)).toarray()) for K in (Kc,Ks))
    assert abs(setting_moment-moment/2)<1e-13
    # The full quantum reference is checked on two nonorthogonal source
    # columns; no pointer or failure record is conditioned away.
    rng=np.random.default_rng(38)
    psi=rng.normal(size=(dimension,2))+1j*rng.normal(size=(dimension,2))
    psi/=np.linalg.norm(psi)
    output=sum(((K@psi).conj().T@(K@psi) for K in (Kc,Ks,Kf)),np.zeros((2,2),complex))
    assert np.max(abs(output-psi.conj().T@psi))<1e-13
    assert np.trace(output).real==pytest.approx(1)


def test_stationarity_and_complex_probe_factors_cannot_be_dropped():
    H=np.array([[0,1],[1,0]],complex);A=np.diag([0,1])
    psi=np.array([1,1j])/np.sqrt(2);rho=np.outer(psi,psi.conjugate())
    comm=H@A-A@H
    first=np.trace(rho@A@comm)
    half_double=np.trace(rho@(A@comm-comm@A))/2
    assert abs(first-half_double)>.4
    H=np.diag([0,2]);rho=np.diag([1,0]);raising=np.array([[0,0],[1,0]])
    first=np.trace(rho@raising.T@(H@raising-raising@H))
    adjoint=np.trace(rho@raising@(H@raising.T-raising.T@H))
    double=np.trace(rho@(raising.T@(H@raising-raising@H)-
                        (H@raising-raising@H)@raising.T))
    assert first==2 and adjoint==0 and double==first+adjoint
    assert first!=double/2


def test_selected_hopping_density_and_fsum_prove_only_their_normalized_tail(receipts):
    for n,receipt in receipts.items():
        row=receipt['result']['selected_density_sum_rule']
        assert F(row['uniform_negative_dimensionless_hopping_density_lower'])>F(3459,100000)
        assert F(row['selected_negative_dimensionless_hopping_density_lower'])>=F(
            row['uniform_negative_dimensionless_hopping_density_lower'])
        assert 0<F(row['original_signed_first_moment_lower'])<=F(
            row['original_signed_first_moment_upper'])<=F(row['simple_signed_first_moment_upper'])
        assert row['normalization']==30*(2**(2*n+1))**3
        assert row['complete_outcomes']==['cosine','sine','failure']
        assert row['random_setting_factor']==1
        assert row['every_fixed_positive_energy_tail_tends_to_zero']
        assert not row['phase_label_identified_with_conserved_momentum']
        assert not row['state_replaced_by_tile_trial']
    assert receipts[0]['result']['selected_density_sum_rule']['one_minus_cos_phase_enclosure']==['2','2']
    assert receipts[1]['result']['selected_density_sum_rule']['one_minus_cos_phase_enclosure']==['1','1']
    values=[]
    for n in range(20):
        a=F(1,2**n);V=2**(6*n+3);eps=F(1,2**(n+4))
        values.append(m.density_positive_tail(a/(6*V),eps,1))
    assert all(y<x for x,y in zip(values,values[1:]))
    assert values[-1]<F(1,1000000)


def test_original_negative_first_moment_is_controlled_without_postselection():
    # A Hermitian swap has both positive and negative transfers. The
    # complete diagonal/failure convention keeps their original weights.
    for E in (F(2),F(9)):
        excited=F(1,100);eps=E*excited
        signed=E*(1-2*excited);negative=E*excited
        assert negative==eps
        assert signed+negative==E*(1-excited)
        assert m.density_positive_tail(signed,eps,E)==1-excited
    for bad in (0,-1,True,.5):
        with pytest.raises(ValueError):m.density_positive_tail(1,0,bad)
