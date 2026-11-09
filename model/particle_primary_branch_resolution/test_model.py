"""Independent full-CAR packets, spectral controls and source authentication."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def _move(state,i,j):
    if not state&(1<<j) or (i!=j and state&(1<<i)):return None
    sign=1 if (state&((1<<j)-1)).bit_count()%2==0 else -1
    state^=1<<j
    if (state&((1<<i)-1)).bit_count()%2:sign=-sign
    return state|(1<<i),sign


def _apply(vector,terms):
    out={}
    for state,value in vector.items():
        for i,j,c in terms:
            hit=_move(state,i,j)
            if hit:out[hit[0]]=out.get(hit[0],0)+value*c*hit[1]
    return {key:value for key,value in out.items() if abs(value)>1e-14}


def test_literal_eight_cell_CAR_and_link_characters_give_the_twelve_packets():
    # All32 q6 CAR modes are retained. The dictionary only stores the actual
    # resulting trial vector; particle-hole sectors are not discarded.
    sites=list(product(range(2),repeat=3));position={x:i for i,x in enumerate(sites)}
    omega=sum(1<<(4*x+j) for x in range(8) for j in (2,3))
    phi=_apply({omega:1.},[(4*x+j,4*x+2+j,.25) for x in range(8) for j in range(2)])
    assert abs(sum(abs(v)**2 for v in phi.values())-1)<1e-14
    pair=m.response.beta_pairs();chart=s.Matrix.hstack(pair[0][0],pair[1][0],pair[0][1],pair[1][1])
    gamma=m.response.density.gamma_data();beta=gamma[3]
    totals={2:0.,4:0.};bond_count=0
    for x in sites:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);bond_count+=1
            spin=np.array(chart.H*((-beta+beta*gamma[axis])/2)*chart,dtype=complex)
            f=[(4*position[x]+i,4*position[y]+j,spin[i,j]) for i,j in product(range(4),repeat=2) if spin[i,j]!=0]
            b=[(j,i,np.conj(c)) for i,j,c in f]
            for terms in (f,b):
                # Different oriented link characters are orthogonal under
                # actual Haar integration. No phase is set to identity.
                for state,value in _apply(phi,terms).items():
                    D=sum(bool(state&(1<<(4*z+j))) for z in range(8) for j in (0,1))
                    D+=sum(not state&(1<<(4*z+j)) for z in range(8) for j in (2,3))
                    assert D in (2,4)
                    totals[D]+=abs(value)**2
    assert bond_count==12 and totals=={2:3.,4:10.5}
    local=m.odd_pair_packet_control()
    assert local['full_local_Fock_dimension']==256
    assert local['native_q6_bond_variance_on_global_pair']==F(9,8)


def test_native_specific_allocation_keeps_exact_full_number_blocks():
    allocation=m.q6_local_allocation();rows=allocation['native_blocks']
    assert [r['onsite_per_endpoint'] for r in rows]==[1,1,1,1,F(1,5)]
    assert allocation['nonconstant_link_margin_lower']>0
    q6=rows[-1]
    assert q6['full_Fock_dimension']==256 and q6['number_sector_count']==25
    assert q6['all_positive_pivots']
    assert m.refinement.q6_positive_shift(F(1,5),q6['inverse_electric_margin'],q6['occupation_loss_upper'])
    assert not m.refinement.q6_positive_shift(F(1,5),q6['inverse_electric_margin'],q6['dyadic_bracket'][0])
    assert allocation['other_native_occupations_and_inactive_rows_retained']


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_actual_selected_state_has_narrower_original_band_and_fixed_time_sign(bound):
    _,receipt=bound;r=receipt['result'];D=r['q6_selected_defects']
    assert F('1.839')<F(D['q6_defect_count_upper'])<F('1.840')
    assert F(D['matched_pair_original_success_lower'])>F(113,128)
    assert F(r['original_transfer_weight_lower'])>F('0.156')
    assert F(r['original_transfer_weight_simple_lower'])==F(3,20)
    assert r['original_transfer_band']==[4,8]
    assert F(r['real_original_correlation_upper'])<F('-0.389')
    assert F(r['original_real_contrast_lower'])>F('1.274')
    assert F(r['final_excitation']['original_final_excitation_weight_lower'])>F(7,80)
    assert not r['individual_coherent_A1u_window_is_this_averaged_window']
    assert not r['isolated_pole_or_experimental_threshold_resolved']
    assert r['physical_mass_table_GeV']==[]


def test_independent_numerical_diagonalization_of_actual_twelve_packet_compression(bound):
    r=bound[1]['result']['odd_endpoint_trial'];rows=r['hopping_packets']
    n=np.array([float(F(x['norm_squared'])) for x in rows]);d=np.array([float(F(x['relative_H0_energy'])) for x in rows])
    matrix=np.zeros((12,12));matrix[0,1:7]=np.sqrt(n);matrix[1:7,0]=np.sqrt(n)
    matrix[1:7,1:7]=np.diag(d)
    g=float(F(r['loop_packet_norm_squared']));matrix[0,7:]=np.sqrt(g);matrix[7:,0]=np.sqrt(g)
    matrix[7:,7:]=np.array([[float(F(v)) for v in row] for row in r['loop_packet_matrix']])
    values=np.linalg.eigvalsh(matrix)
    assert (values<0).sum()==1
    assert abs(values[0]-float(F(r['relative_Ritz_root_bracket'][0])))<1e-12
    trial=np.r_[1,-np.sqrt(n)*np.array([float(F(x)) for x in r['rational_hopping_coefficients']]),
                -np.sqrt(g)*np.array([float(F(x)) for x in r['rational_loop_coefficients']])]
    energy=trial@matrix@trial/(trial@trial)
    assert abs(energy-float(F(r['trial_energy_shift_relative_to_anchor_pair'])))<1e-12
    assert abs(float(F(r['odd_sector_absolute_ground_upper']))-(-714+energy))<1e-12
    assert F('5.913')<F(r['difference_of_two_trial_upper_bounds'])<F('5.914')
    assert not r['difference_of_trial_uppers_is_excitation_bound']
    assert F(r['dressed_trial_full_complement_residual_upper'])>100
    assert not r['trial_space_invariant'] and not r['isolated_branch_proved']


def test_bare_odd_trial_and_ground_trial_have_different_exact_variances(bound):
    r=bound[1]['result']['odd_endpoint_trial']
    assert F(r['anchor_pair_full_source_variance'])==F(363,2)+F(6,4410)
    ground_variance=F(180)+F(6,4410)
    assert F(r['anchor_pair_full_source_variance'])-ground_variance==F(3,2)
    assert r['cubic_irrep']=='A1' and r['spatial_inversion_parity']==-1
    assert not r['selected_state_substituted_by_trial']
    assert not r['odd_sector_bottom_is_strictly_above_full_ground_proved']


def test_complete_nonhermitian_pair_Ramsey_retains_failure_and_entangled_reference():
    B=m.channels.car_bilinear(2,0,1);I=s.eye(4)
    H=s.diag(0,3,-3,0);failure=I-B.H*B
    rows=m.apparatus.ramsey_instrument(H,(B,failure),'1/2',steps=64)
    assert len(rows)==4
    assert m.apparatus.clean(sum((r['K'].H*r['K'] for r in rows),s.zeros(4)))==I
    rho=s.diag(s.Rational(1,10),s.Rational(1,10),s.Rational(7,10),s.Rational(1,10))
    psi=s.zeros(16,1)
    for j in range(4):psi[5*j]=s.sqrt(rho[j,j])
    outputs=[s.kronecker_product(row['K'],I)*psi for row in rows]
    probabilities=[s.simplify((v.H*v)[0]) for v in outputs]
    assert s.simplify(sum(probabilities))==1
    assert sum(float(x) for x in probabilities[2:])>.29
    score=sum(row['original_score']*p for row,p in zip(rows,probabilities))
    U=m.apparatus.cayley(H,'1/2',64)
    direct=s.re(s.trace(rho*B.H*U.H*B*U))
    assert s.simplify(score-direct)==0
    assert float(score)<-.69
    # Score is unconditioned: dividing by the success mass changes its value.
    assert abs(float(score)/.7-float(score))>.2


def test_moment_and_negative_correlation_do_not_determine_one_line():
    mass=F(113,128)
    for energies,weights in (([6],[mass]),([F(9,2),F(15,2)],[mass/2,mass/2])):
        H=np.diag([0]+[float(x) for x in energies]);B=np.zeros_like(H)
        B[1:,0]=np.sqrt([float(x) for x in weights]);rho=np.diag([1.]+[0.]*len(energies))
        assert np.linalg.norm(B,2)<=1
        residual=H@B-B@H-6*B
        moment=np.trace(rho@residual.T@residual)
        assert moment<=float(F(102,35))
        Fhalf=sum(float(w)*np.exp(.5j*float(e)) for w,e in zip(weights,energies))
        assert Fhalf.real<-.375
        assert len(np.flatnonzero(B[:,0]))==len(energies)


def test_complete_clock_and_Ramsey_errors_preserve_public_original_bounds(bound):
    r=bound[1]['result'];p=r['complete_Ramsey_plan']
    assert F(p['source_time'])==F(1,2) and p['Cayley_steps']>0
    assert F(p['signed_original_score_error_upper'])<=F(1,512)
    assert F(r['finite_Ramsey_signed_original_score_upper'])<-F(191,512)
    assert F(r['finite_reported_original_weight_lower'])>=F(9,64)
    assert r['finite_reported_transfer_band']==['63/16','129/16']
    assert p['complete_Ramsey_outcomes']==34
    assert p['extra_compiler_or_preparation_error'] is None
    assert not p['source_energy_zero_changed'] and not p['complete_primary_histogram_executed']


def test_source_state_clock_and_numerical_promotions_are_rejected(bound):
    member,receipt=bound;m.verify_bound(member,receipt)
    bad=deepcopy(receipt);bad['result']['original_transfer_weight_lower']='1'
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    bad=deepcopy(receipt);bad['result']['odd_endpoint_trial']['isolated_branch_proved']=True
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    altered=deepcopy(member);altered['state']['target_ground_excess']='0'
    with pytest.raises(ValueError):m.bind_member(altered)
    altered=deepcopy(member);altered['clock']['seconds_per_tick']=1.0
    with pytest.raises(ValueError):m.bind_member(altered)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))


def test_declared_inputs_exist_and_proof_contains_no_control_bytes():
    assert all((m.contract.ROOT/path).is_file() for path in m.SOURCE_INPUTS)
    proof=(m.contract.ROOT/'model/particle_primary_branch_resolution/PROOF.md').read_bytes()
    assert all(x>=32 or x==10 for x in proof)


def test_exact_full_CAR_vacuum_loss_is_necessary_but_not_sufficient():
    # The balanced 36-state block retains active and inactive rows. Its
    # zero-defect vector tests every alleged scalar lower shift.
    _,_,D,Q=next(row for row in m.refinement.q6_integer_blocks() if row[0]==(2,2))
    D=s.Matrix(D);Q=s.Matrix(Q);kernel=D.nullspace()
    assert len(kernel)==1
    omega=kernel[0];norm=(omega.T*omega)[0]
    assert (omega.T*Q*omega)[0]==norm
    cross=(Q-s.eye(36))*omega
    assert (cross.T*cross)[0]==norm and D*cross==4*cross
    local=D-F(1,8)*Q
    assert (omega.T*local*omega)[0]==-norm/8
    # c=w meets the necessary vacuum expectation, but omits its coupling.
    assert not m.refinement.q6_positive_shift(1,F(1,8),F(1,8))
    assert np.linalg.eigvalsh(np.array(local,dtype=float))[0]<-.125


def test_actual_cube_flux_support_requires_four_nonconstant_edges():
    # Once the physical ghost and matter-singlet restrictions are imposed,
    # every nontrivial edge support has no degree-one vertex.
    sites=list(product(range(2),repeat=3));position={x:i for i,x in enumerate(sites)}
    edges=[]
    for x in sites:
        for axis in range(3):
            if x[axis]==0:
                y=list(x);y[axis]=1
                edges.append((position[x],position[tuple(y)]))
    smallest=12
    for mask in range(1,1<<len(edges)):
        degrees=[0]*8
        for j,(x,y) in enumerate(edges):
            if mask&(1<<j):degrees[x]+=1;degrees[y]+=1
        if 1 not in degrees:smallest=min(smallest,mask.bit_count())
    assert len(edges)==12 and smallest==4
    assert smallest*F(105,2)*m.scale.casimir(0,0,0,6)==63
    assert all(charge!=0 for _,_,charge in m.refinement.native_control())


def test_all_retention_case_ceilings_fail_against_authenticated_ground_upper(bound):
    member,receipt=bound;r=receipt['result']['retained_penalty_obstruction']
    ceilings=[F(row['comparison_floor_strict_upper']) for row in r['cases']]
    assert ceilings==[-F(18016725,24871),-F(198863986,274505)]
    assert max(ceilings)<-724<F(r['full_physical_ground_upper'])
    assert all(F(row['comparison_margin_strict_upper'])<-F(3,5) for row in r['cases'])
    assert r['all_admitted_retention_parameters_excluded_for_this_isolation_test']
    assert not r['ceiling_is_actual_excited_energy_upper']
    assert not r['source_has_no_gap_proved']
    assert r['sharper_ground_upper_could_change_this_comparison']
    # This is a comparison of a method ceiling with U, not proof of a gap:
    # a sufficiently sharper U would reverse the arithmetic comparison.
    assert -725<min(ceilings)
    bad=deepcopy(receipt)
    bad['result']['retained_penalty_obstruction']['source_has_no_gap_proved']=True
    with pytest.raises(ValueError):m.verify_bound(member,bad)
