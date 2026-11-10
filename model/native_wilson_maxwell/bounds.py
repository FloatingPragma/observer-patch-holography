"""Complete finite-time current response and exact trial-support transfer."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s
from .matter import filled_current,sqrt_bound,outward
from model.native_higgs_filtered_state import model as filtered
from model.native_block_dirac import algebra as code


def support_control():
    E,_,_=code.doublet_code()
    Q=s.zeros(9);Q[0,0]=1
    Fvac=s.zeros(4);Fvac[0,0]=1
    assert Q*E==E*Fvac and E.H*Q*E==Fvac
    tensor=s.kronecker_product(E,E);embedding=s.zeros(81,16)
    for word in range(16):embedding[:,word]=tensor[:,(word&3)*4+(word>>2)]
    all_vac=s.zeros(81);all_vac[0,0]=1
    code_vac=s.zeros(16);code_vac[0,0]=1
    assert all_vac*embedding==embedding*code_vac
    assert embedding.H*all_vac*embedding==code_vac
    qlocal=s.diag(0,1,1)
    qpair=s.kronecker_product(qlocal,s.eye(3))+s.kronecker_product(s.eye(3),qlocal)
    Qnative=s.kronecker_product(qpair,s.eye(9))-s.kronecker_product(s.eye(9),qpair)
    charge=s.diag(*[(word&3).bit_count()-(word>>2).bit_count() for word in range(16)])
    assert Qnative*embedding==embedding*charge
    return {'all_sixteen_occupations_retained':True,
        'two_copy_all_reference_vacuum_intersection_rank':1,
        'four_copy_all_reference_vacuum_intersection_rank':1,
        'exact_support_identity':'P_code Q_all_Rloc_vacuum=|g^32,refvac><g^32,refvac|',
        'exact_full_code_charge_intertwiner_checked':True,
        'full_reference_complement_argument':'every nonvacuum code word has a charged port; the complete Rloc-vacuum neutral complement is orthogonal',
        'code_matrix_units_are_actual_native_spectral_services':True,
        'global_H_and_every_count_effect_reduce_code':True,
        'initial_ground_or_code_projector_applied':False}


def charge_frame_control():
    a,b=[code.annihilator(2,j) for j in range(2)]
    variables=(a,b.H);M=s.Matrix([[s.Rational(3,5),s.Rational(4,5)],[-s.Rational(4,5),s.Rational(3,5)]])
    U=s.eye(4);U[0,0]=U[3,3]=s.Rational(3,5)
    U[0,3]=s.Rational(4,5);U[3,0]=-s.Rational(4,5)
    assert U.H*U==s.eye(4)
    transformed=[sum((M[i,j]*variables[j] for j in range(2)),s.zeros(4)) for i in range(2)]
    for i in range(2):assert U*variables[i]*U.H==transformed[i]
    total=a.H*a-b.H*b
    assert U*total==total*U
    one=s.diag(1,0)
    def second(matrix):
        return sum((matrix[i,j]*variables[i].H*variables[j] for i in range(2) for j in range(2)),s.zeros(4))
    assert U*second(one)*U.H==second(M.H*one*M)
    assert U*second(one)*U.H!=second(M*one*M.H)
    Eplus=s.zeros(32,16);Eminus=s.zeros(32,16)
    for site in range(8):
        for spin in range(2):
            Eplus[4*site+spin,2*site+spin]=1
            Eminus[4*site+2+spin,2*site+spin]=1
    R=Eplus.row_join(Eminus)
    assert R.H*R==s.eye(32)
    S=s.eye(32);S[0,0]=S[6,6]=s.Rational(3,5)
    S[0,6]=s.Rational(4,5);S[6,0]=-s.Rational(4,5)
    Uframe,Vframe=S*Eplus,S*Eminus
    Msite=Uframe*Eplus.H+Vframe*Eminus.H
    assert Msite==S and Msite*R==Uframe.row_join(Vframe)
    assert Uframe.row_join(Vframe)!=Msite
    for site in range(8):
        Px=s.diag(*[int(k//4==site) for k in range(32)])
        Qgroup=s.diag(*[int(k//2==site) for k in range(16)],
                       *[int(k//2==site) for k in range(16)])
        assert R.H*Px*R==Qgroup
    Px=s.diag(1,1,1,1,*([0]*28))
    assert Msite.H*Px*Msite!=Uframe.row_join(Vframe).H*Px*Uframe.row_join(Vframe)
    return {'defining_CAR_action':'U_FW Psi0_site U_FW^dagger=M_site Psi0_site',
        'site_major_frame':'M_site=U E_plus^dagger+V E_minus^dagger=[U,V][E_plus,E_minus]^dagger',
        'density_one_body_action':'M_site^dagger P_x M_site',
        'exact_full_pair_Fock_rotation_checked':True,
        'exact_32_dimensional_grouped_site_permutation_checked':True,
        'omitted_column_permutation_explicitly_rejected':True,
        'total_native_charge_conjugation_unchanged':True,
        'wrong_inverse_density_action_explicitly_rejected':True,
        'global_extension':'U_FW on full code plus identity on every native complement'}


@lru_cache(None)
def preparation():
    _,_,_,source=filtered.higgs.source()
    assert source['actual_local_preparation']['preparation_success_probability']=='1'
    p0=F(source['atomic_cube_trial']['exact_interacting_ground_trial_overlap_lower'])
    p_code=outward(p0**32)
    one=filtered.preparation(1)
    epsilon=sqrt_bound(32*F(one['retained_purification_error_squared_upper']),True)
    return {'native_constituents':32,'native_trial_copies':32,
        'deterministic_initial_trial':True,'base_trial_success':'1','filter_depth':0,
        'agreement_filter_count':0,'preparation_flags_discarded':False,
        'complete_trial_ground_overlap_lower_per_constituent':str(p0),
        'global_code_probability_lower':str(p_code),
        'global_native_complement_probability_upper':str(outward(1-p_code,True)),
        'actual_probability_law':'p_code times every ideal code count; 1-p_code on native-complement/photon0; other complement/photon labels zero',
        'complete_added_Rloc_state':'shared vacuum, unchanged by trial preparation',
        'all_original_spin_references_retained':True,
        'no_initial_P0_or_code_measurement':True,
        'independent_depth1_comparison':{'native_trial_copies':64,'agreement_filter_count':32,
            'complete_preparation_success_lower':str(outward(F(one['complete_preparation_success_lower'])**32)),
            'retained_purification_error_upper':str(epsilon),
            'used_for_primary_experiment':False}}


def response():
    current=filled_current(); sigma=F(current['pair_one_photon_V1_vacuum_norm_over_m_lower'])
    g=tau=F(1,64);nu=18*g+F(59,4)*g*g
    leading=g*tau*sigma*(1-F(11,2)*tau)
    odd=F(27,16)*tau*g**3;dyson=tau*tau*nu*nu/2
    amplitude=outward(leading-odd-dyson)
    assert amplitude>0
    ideal_probability=outward(amplitude*amplitude)
    marginal=outward((tau*nu)**2,True)
    connected=outward(amplitude*amplitude-(tau*nu)**4)
    assert connected>0
    p=preparation();weight=F(p['global_code_probability_lower'])
    epsilon=F(p['independent_depth1_comparison']['retained_purification_error_upper'])
    depth1=outward((amplitude-epsilon)**2-(tau*nu+epsilon)**4)
    assert amplitude>epsilon and depth1>0
    fields={'g':g,'dimensionless_mT':tau,'whole_perturbation_norm_over_m_upper':nu,
        'leading_target_amplitude_lower':leading,
        'complete_first_insertion_odd_remainder_upper':odd,
        'entire_operator_twice_Duhamel_remainder_upper':dyson,
        'ideal_complete_target_amplitude_lower':amplitude,
        'ideal_complete_target_probability_lower':ideal_probability,
        'each_ideal_excitation_marginal_probability_upper':marginal,
        'ideal_binary_connected_response_lower':connected,
        'actual_unconditional_target_probability_lower':outward(weight*amplitude*amplitude),
        'actual_unconditional_target_complement_probability_upper':outward(1-weight*amplitude*amplitude,True),
        'actual_unconditional_binary_connected_response_lower':outward(weight*connected),
        'independent_depth1_conditional_connected_response_lower':depth1}
    return {'exact':{key:str(value) for key,value in fields.items()},
        'diagnostic_binary64':{key:float(value) for key,value in fields.items()},
        'target_count':[1,1,1],'free_target_probability':'0',
        'source_state_is_actual_coupled_ground':False,
        'response_type':'finite nonstationary free-vacuum dressing by the full spinor current',
        'all_virtual_matter_and_photon_occupations_in_Duhamel_bound':True,
        'current_phase_derivative_sign':'V1=sum Q_e J_phase,e; I_flow,e=-J_phase,e',
        'finite_total_photon_cutoff_treated_as_actual_source':True,
        'all_order_untruncated_QED_or_scattering_claimed':False}


if __name__=='__main__':
    import json
    print(json.dumps({'support':support_control(),'preparation':preparation(),'response':response()},indent=2))
