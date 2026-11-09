"""Same-source pair concentration and a dressed odd variational sector.

All probabilities concern the selected polynomial state. The odd trial is
an independent variational vector and supplies no selected-state residue.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_refinement_response import model as refinement
from model.particle_primary_response_moments import model as response
from model.particle_primary_pair_channels import model as channels
from model.particle_primary_variational_readout import variational
from model.particle_primary_spectral_scale import scale
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_mass_readout import readout
from model.particle_filtered_probe import filter as apparatus
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_branch_resolution/__init__.py',
    'model/particle_primary_branch_resolution/__main__.py',
    'model/particle_primary_branch_resolution/model.py',
    'model/particle_primary_branch_resolution/test_model.py',
    'model/particle_primary_branch_resolution/PROOF.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_refinement_response/PROOF.md',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_filtered_probe/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
    'Lean/QFT/PeriodicMomentum.lean',
)


def _member_zero(member):
    contract.verify_member(member)
    if member['source']['side']!=2 or F(member['source']['spacing'])!=1:
        raise ValueError('the pair-window and odd trial certificate requires member zero')


@lru_cache(maxsize=1)
def q6_local_allocation():
    """Reuse full-occupation Schur and LDL bounds with degree-three allocation."""
    native=refinement.native_control()
    retained=(F(),)*4+(F(4,5),)
    onsite=tuple(1-r for r in retained)
    delta=4*sum((d*(sqrt_upper(m*m+F(1,4))-m)
                 for (d,_,_),m in zip(native,onsite)),F())
    margin=F(63,4)-delta
    if margin<=0:raise ArithmeticError('nonconstant link block lost positivity')
    rows=[]
    for (d,casimir,charge),m,r in zip(native,onsite,retained):
        w=1/(F(105,2)*casimir-delta)
        local=(refinement.exact_q6_loss(m,w) if d==1
               else occupation.occupation_loss(d,m,w))
        rows.append({'dimension':d,'Q6':charge,'onsite_per_endpoint':m,
                     'onsite_retained_fraction':r,'inverse_electric_margin':w,**local})
    return {'maximum_degree':3,'pointwise_fermion_shift_upper':delta,
        'nonconstant_link_margin_lower':margin,'native_blocks':rows,
        'per_bond_loss_upper':sum((x['occupation_loss_upper'] for x in rows),F()),
        'retained_q6_defect_coefficient':F(12,5),
        'other_native_occupations_and_inactive_rows_retained':True}


def q6_selected_defects(member):
    _member_zero(member)
    local=q6_local_allocation();upper=variational.bind_member(member)
    variational.verify_bound(member,upper)
    source=member['source'];V=source['side']**3
    energy=F(upper['result']['ground_energy_upper'])
    eps=F(member['state']['target_ground_excess'])
    loss=source['interior_hopping_bonds']*local['per_bond_loss_upper']
    loss+=F(source['interior_plaquettes'],7)
    defects=(energy+eps+90*V+loss)/F(12,5)
    if not 0<=defects<F(15,8):raise ArithmeticError('selected q6 defect bound failed')
    return {'allocation':local,'full_physical_ground_upper':energy,
        'ground_upper_producer_sha256':contract.digest(upper),
        'operator_floor_constant':-90*V-loss,
        'selected_ground_excess_upper':eps,'q6_defect_count_upper':defects,
        'simple_q6_defect_count_upper':F(15,8),
        'matched_pair_original_success_lower':1-defects/(2*V),
        'simple_original_success_lower':F(113,128)}


def _act(vector,terms):
    answer={}
    for state,value in vector.items():
        for i,j,coefficient in terms:
            move=refinement._bilinear(state,i,j)
            if move:
                answer[move[0]]=answer.get(move[0],0)+value*coefficient*move[1]
    return {key:s.simplify(value) for key,value in answer.items() if value!=0}


@lru_cache(maxsize=1)
def odd_pair_packet_control():
    """Exact CAR action on all occupations, applied to the local trial vector."""
    gamma=response.density.gamma_data();beta=gamma[3];pairs=response.beta_pairs()
    chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
    vacuum=sum(1<<j for j in (2,3,6,7))
    phi=_act({vacuum:s.S.One},[(off+j,off+2+j,s.S.One)
                               for off in (0,4) for j in range(2)])
    if sum(s.conjugate(v)*v for v in phi.values())!=4:
        raise ArithmeticError('the two-site scalar pair normalization changed')
    rows=[]
    for axis in range(3):
        hop=apparatus.clean(chart.H*((-beta+beta*gamma[axis])/2)*chart)
        forward=[(i,j+4,hop[i,j]) for i in range(4) for j in range(4) if hop[i,j]!=0]
        reverse=[(j+4,i,s.conjugate(hop[i,j])) for i in range(4) for j in range(4) if hop[i,j]!=0]
        parts={};vacuum_norm=s.S.Zero
        for terms in (forward,reverse):
            vacuum_norm+=sum(s.conjugate(v)*v for v in _act({vacuum:s.S.One},terms).values())
            for state,value in _act(phi,terms).items():
                D=sum(bool(state&(1<<j)) for j in (0,1,4,5))
                D+=sum(not state&(1<<j) for j in (2,3,6,7))
                parts[D]=parts.get(D,0)+s.conjugate(value)*value
        parts={key:s.simplify(value) for key,value in parts.items()}
        if vacuum_norm!=1 or parts!={4:s.Integer(2),2:s.Integer(4)}:
            raise ArithmeticError('actual odd pair electric packets changed')
        rows.append({'axis':axis,'vacuum_bond_variance':1,
                     'two_site_pair_norm_squared':4,
                     'one_pair_hopping_norm_squared':4,
                     'two_pair_hopping_norm_squared':2})
    return {'full_local_Fock_dimension':256,'axes':rows,
        'one_pair_per_bond_global_norm_squared':F(1,4),
        'two_pair_per_bond_global_norm_squared':F(7,8),
        'native_q6_bond_variance_on_global_pair':F(9,8)}


def odd_ritz(member):
    """Exact physical A1u endpoint compression with a full complement bound."""
    _member_zero(member);control=odd_pair_packet_control()
    source=member['source'];B=source['interior_hopping_bonds'];P=source['interior_plaquettes']
    native=scale.vacuum_channels()['all_native_blocks']
    rows=[{'Q6':row['integer_Q6_charge'],'norm_squared':F(B*row['dimension']),
           'relative_H0_energy':F(row['bond_H0_excitation'])}
          for row in native if row['integer_Q6_charge']!=6]
    rows += [{'Q6':6,'pair_number':1,'norm_squared':F(B,4),'relative_H0_energy':F(63,4)},
             {'Q6':6,'pair_number':2,'norm_squared':F(7*B,8),'relative_H0_energy':F(87,4)}]
    norm=[row['norm_squared'] for row in rows];diagonal=[row['relative_H0_energy'] for row in rows]
    g=F(P,22050);count=variational.character_control()['cubic_loop_singlet_matrix']
    correction=[[-F(v,420) for v in row] for row in count]
    loop=[[correction[i][j]+(F(native[i]['loop_H0_excitation']) if i==j else 0)
           for j in range(5)] for i in range(5)]
    if any(loop[i][i]<=sum(abs(loop[i][j]) for j in range(5) if j!=i) for i in range(5)):
        raise ArithmeticError('loop packet matrix lost positive dominance')
    def coefficients(d):
        h=[1/(x+d) for x in diagonal]
        l=variational.solve([[loop[i][j]+(d if i==j else 0) for j in range(5)] for i in range(5)],[1]*5)
        return h,l
    def schur(d):
        h,l=coefficients(d)
        return sum(x*t for x,t in zip(norm,h))+g*sum(l)
    lo,hi=F(),F(16)
    if schur(hi)>=hi:raise ArithmeticError('odd Schur upper bracket failed')
    while hi-lo>F(1,2**60):
        mid=(lo+hi)/2
        if schur(mid)>=mid:lo=mid
        else:hi=mid
    h,l=coefficients(lo);q=sum(x*t for x,t in zip(norm,h))+g*sum(l)
    square=sum(x*t*t for x,t in zip(norm,h))+g*sum(t*t for t in l)
    e=(-q-lo*square)/(1+square)
    if not -hi<=e<=-lo:raise ArithmeticError('odd rational Rayleigh quotient escaped its bracket')
    # This bound includes the entire omitted finite space. The weighted q6
    # projectors need not commute with V; no parent commutation is used.
    vnorm=30*B+F(P,7)
    projected_loop=[sum(correction[i][j]*l[j] for j in range(5)) for i in range(5)]
    outside=vnorm*vnorm*square-q*q-g*sum(x*x for x in projected_loop)
    if outside<0:raise ArithmeticError('complement norm bound is inconsistent')
    residual2=((q+e)**2+(lo+e)**2*square+outside)/(1+square)
    ground=variational.bind_member(member);variational.verify_bound(member,ground)
    lower=F(occupation.retained_bound(member)['physical_ground_lower'])
    absolute=-720+6+e
    return {'packet_dimension':12,'cubic_irrep':'A1','spatial_inversion_parity':-1,
        'trial_is_CAR_even':True,'local_CAR_packet_control':control,
        'hopping_packets':rows,'loop_packet_matrix':loop,'loop_packet_norm_squared':g,
        'relative_Ritz_root_bracket':[-hi,-lo],
        'rational_hopping_coefficients':h,'rational_loop_coefficients':l,
        'trial_norm_squared':1+square,'trial_energy_shift_relative_to_anchor_pair':e,
        'odd_sector_absolute_ground_upper':absolute,
        'full_physical_ground_lower':lower,
        'odd_sector_energy_above_full_ground_upper':absolute-lower,
        'odd_sector_bottom_is_strictly_above_full_ground_proved':False,
        'ground_Ritz_upper_for_comparison':F(ground['result']['ground_energy_upper']),
        'difference_of_two_trial_upper_bounds':absolute-F(ground['result']['ground_energy_upper']),
        'difference_of_trial_uppers_is_excitation_bound':False,
        'anchor_pair_full_source_variance':sum(norm)+5*g,
        'dressed_trial_full_complement_residual_squared_upper':residual2,
        'dressed_trial_full_complement_residual_upper':sqrt_upper(residual2),
        'trial_space_invariant':False,'selected_state_substituted_by_trial':False,
        'selected_residue_of_trial_proved':False,'isolated_branch_proved':False}


def ramsey_plan(member):
    """Finite rational Cayley implementation on the original energy origin."""
    _member_zero(member);source=member['preparation_plan']['source']
    norm=max(abs(F(source['spectral_lower'])),abs(F(source['spectral_upper'])))
    time=F(1,2);error=F(1,1024)
    needed=time**3*norm**3/(12*error)
    ceiling=-(-needed.numerator//needed.denominator)
    steps=isqrt(ceiling)
    if steps*steps<ceiling:steps+=1
    clock=apparatus.cayley_service_bound(time,time,norm,max(1,steps))
    if F(clock['unitary_operator_error'])>error:raise ArithmeticError('Ramsey evolution budget failed')
    return {'source_time':time,'source_norm_upper':norm,'Cayley_steps':max(1,steps),
        'clock_bound':clock,'complete_Ramsey_outcomes':34,
        'signed_original_score_error_upper':2*F(clock['unitary_operator_error']),
        'complete_all_reference_halftrace_error_upper':2*F(clock['unitary_operator_error']),
        'extra_compiler_or_preparation_error':None,'source_energy_zero_changed':False,
        'complete_primary_histogram_executed':False}


def retained_penalty_obstruction(member):
    """Exact ceiling on the scalar-loss Schur comparison, not on E1 itself.

    Native-specific onsite retentions and any admitted electric retention
    are allowed. Vacuum expectations bound every valid local scalar loss
    from below, even if its full occupation spectrum is solved exactly.
    """
    _member_zero(member)
    native=refinement.native_control()
    ground=variational.bind_member(member);variational.verify_bound(member,ground)
    upper=F(ground['result']['ground_energy_upper'])
    # Delta(0)>7 and Delta(1/2)>12 follow by positive rational squaring.
    if not F(67,60)**2<F(5,4) or not F(7,10)**2<F(1,2):
        raise ArithmeticError('pointwise-shift rational lower bounds failed')
    rows=[]
    for domain,delta,penalty in (('0 <= min(r_j) <= 1/2',7,3),
                                  ('1/2 <= min(r_j) < 1',12,6)):
        local=sum((F(d)/(F(105,2)*C-delta) for d,C,_ in native),F())
        loss=12*local+F(6,7)
        ceiling=-720-loss+penalty
        if not ceiling<-724<upper:
            raise ArithmeticError('retained-penalty exclusion lost its strict margin')
        rows.append({'retained_native_onsite_domain':domain,
            'pointwise_fermion_shift_strict_lower':delta,
            'unavoidable_per_bond_scalar_loss_strict_lower':local,
            'total_scalar_loss_strict_lower':loss,
            'attainable_scalar_odd_penalty_upper':penalty,
            'comparison_floor_strict_upper':ceiling,
            'comparison_margin_strict_upper':ceiling-upper})
    return {'allocation':'spatially uniform per-native onsite fractions; common electric fraction',
        'source_member_level':0,'source_endpoint_coupling':1,
        'local_loss_family':'m_j D_j - w_j Q_j >= -c_j I, w_j=(105(1-s)C_j/2-Delta)^(-1)',
        'vacuum_expectation_forces_c_j_at_least_w_j_d_j':True,
        'physical_one_matter_defect_excluded_by_total_Q6':True,
        'matter_zero_defect_nonconstant_flux_electric_floor':63,
        'common_retention_penalty_gap':'min(6*r,63*s)',
        'native_specific_retention_penalty_gap':'min(6*min(r_j),63*s)',
        'cases':rows,'full_physical_ground_upper':upper,
        'ground_upper_producer_sha256':contract.digest(ground),
        'all_admitted_retention_parameters_excluded_for_this_isolation_test':True,
        'exact_q6_local_diagonalization_evades_obstruction':False,
        'ceiling_is_actual_excited_energy_upper':False,
        'source_has_no_gap_proved':False,
        'sharper_ground_upper_could_change_this_comparison':True,
        'different_Schur_or_complement_estimates_excluded':False}


def selected_resolution(member):
    _member_zero(member);response.source_probe_control();channels.symmetry_summary()
    defects=q6_selected_defects(member);w=defects['matched_pair_original_success_lower']
    moment=F(102,35)  # sqrt(2)<99/70, verified by squaring.
    if F(99,70)**2<=2 or F(12,7)**2<moment:raise ArithmeticError('rational radical caps failed')
    band=w-moment/4
    if band<F(3,20):raise ArithmeticError('selected narrow-window floor failed')
    eps=F(member['state']['target_ground_excess'])
    final=readout.excitation_window((4,8),band,eps,1)
    # cos(omega/2)=-cos((omega-6)/2+3-pi), pi-3<1/7.
    correlation=-F(97,98)*w+moment/8+F(12,98)
    if correlation>=-F(3,8):raise ArithmeticError('fixed-time negative correlation failed')
    clock=spectroscopy.clock_plan(member,'1/32','1/512','1/512')
    clockloss=2*F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['two_complete_clock_error_upper'])
    if clockloss>F(1,128):raise ArithmeticError('two-clock budget failed')
    ramsey=ramsey_plan(member)
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_endpoint_coupling':'1','q6_selected_defects':defects,
        'probe':'uniform original average of sixteen matched q6 beta pairs',
        'full_Fock_centered_second_moment_upper':moment,'source_frequency_center':6,
        'original_transfer_band':[4,8],'original_transfer_weight_lower':band,
        'original_transfer_weight_simple_lower':F(3,20),
        'final_excitation':final,'final_excitation_simple_band':[4,9],
        'final_excitation_simple_weight_lower':F(7,80),
        'finite_two_clock_plan':clock,'finite_reported_transfer_band':[F(63,16),F(129,16)],
        'finite_reported_original_weight_lower':band-clockloss,
        'finite_reported_original_weight_simple_lower':F(9,64),
        'correlation_source_time':F(1,2),'real_original_correlation_upper':correlation,
        'real_original_correlation_simple_upper':-F(3,8),
        'original_real_contrast_lower':w-correlation,
        'complete_Ramsey_plan':ramsey,
        'finite_Ramsey_signed_original_score_upper':correlation+ramsey['signed_original_score_error_upper'],
        'finite_Ramsey_signed_original_score_simple_upper':-F(191,512),
        'odd_endpoint_trial':odd_ritz(member),
        'retained_penalty_obstruction':retained_penalty_obstruction(member),
        'source_state_reference_clock_unchanged':True,
        'all_failures_and_idle_reference_retained':True,
        'probabilities_conditioned_on_success':False,
        'individual_coherent_A1u_window_is_this_averaged_window':False,
        'isolated_pole_or_experimental_threshold_resolved':False,
        'source_momentum_or_rest_frame_established':False,
        'uniform_refinement_claimed':False,'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_resolution(member)
    return {'schema':'oph.primary_branch_resolution_bound.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('pair-resolution bound or actual member changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*refinement.SOURCE_INPUTS,
        *response.SOURCE_INPUTS,*variational.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS))
    if any(not (contract.ROOT/name).is_file() for name in paths):
        raise ValueError('a branch-resolution evidence path is absent')
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_branch_resolution_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'physical_mass_table_GeV':[],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('resolution inputs changed during replay')
    return result
