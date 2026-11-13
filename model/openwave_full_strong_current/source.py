"""Full autonomous five-body response, with a rigorous complete-carrier remainder."""
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict
import sympy as sy
from .native import *
from model.openwave_strange_shell_current.native import contraction,frame
from model.openwave_scattering_current_support.native import proton_fields,pion_field,product_column
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from model.openwave_capacity_stress.source import sqrt_interval,enclosure


def record(M):return [[str(v)for v in row]for row in M.tolist()]


def rest_floor(A,G):
    return min(F(A[i,i]/G[i,i])-sum((F(abs(A[i,j]))*sqrt_upper(F(1,G[i,i]*G[j,j]))for j in range(A.rows)if i!=j),F())for i in range(A.rows))


@lru_cache(None)
def source_identity():
    _,columns,A,G=triquark_frame();reverse={s:(i,v)for i,col in enumerate(columns)for s,v in col.items()};B=sy.zeros(100)
    for j,col in enumerate(columns):
        for shift in(-1,0,1):
            for s,v in action_vector(base_action,{translate(s,shift):v for s,v in col.items()}).items():
                if s in reverse:
                    i,c=reverse[s];B[i,j]+=sy.Rational(c*v)
    assert B==A
    _,mc,MA,MG=meson_frame()
    # Exact complete spin-zero charged-pion intertwiner, including spin-triplet complement.
    from model.openwave_cluster_bandwidth.meson import matrices
    shapes=meson_frame()[0];P=sy.zeros(48,3)
    for j,pos in enumerate(((0,0),(0,1),(1,0))):
        for sq,sa,sign in ((0,1,1),(1,0,-1)):
            P[shapes.index(((pos[0],sq,1),(pos[1],sa,0))),j]=sign
    M=matrices(F(1),F(0),repair=False)[0]
    assert P.T*MG*P==6*sy.eye(3) and MA*P==MG*P*M
    floor_tri=rest_floor(A,G);floor_meson=rest_floor(MA,MG);assert min(floor_tri,floor_meson)>0
    zeta=repair_matrix()[2]
    return {'triquark_floor':floor_tri,'meson_floor':floor_meson,'triquark_J_norm':zeta/floor_tri,
            'meson_J_norm':zeta/floor_meson,'triquark_blocks':correlated_data('triquark')[-1],
            'meson_blocks':correlated_data('meson')[-1]}


def support_norm(span):
    """||H P_span|| bound on the complete Nq4/Na1 carrier, not a finite cutoff."""
    p=coefficients();d=source_identity();N=5;C=F(4*N*N,3);S=F(N*(N+2),4)
    return (abs(p['mbar']+p['t'])*N+abs(p['delta'])*N/2+abs(p['A'])*N*N/9
        +abs(p['U']-p['kappa']/2)*C+p['kappa']*span*C
        +(abs(p['Js']-p['spin_gradient'])+6*abs(p['spin_gradient']))*S
        +(abs(p['Jw']-p['eta'])+2*abs(p['eta']))*S+2*p['t']*N
        +8*abs(p['pair_triplet'])+F(3,4)*abs(p['Jw']+p['eta'])*3*N*N
        +16*(d['triquark_J_norm']+d['meson_J_norm']))


@lru_cache(None)
def calculation():
    identity=source_identity();fields=proton_fields()
    Y=[native_vector(product_column(fields[0],pion_field(1,1))),native_vector(product_column(fields[5],pion_field(0,0)))]
    d=contraction();_,P=frame();U=d['matrix']*P
    vectors=[native_vector({mask:U[i,j]for i,mask in enumerate(d['raw_output_masks'])if U[i,j]})for j in range(7)]
    for vector in Y+vectors:assert not action_vector(colour_casimir,vector,(0,1))
    assert all(dot(y,v)==0 for y in Y for v in vectors)
    components={};counts={}
    for name,action in [('base',base_action),('triquark',lambda s:correlated_action(s,'triquark')),('meson',lambda s:correlated_action(s,'meson'))]:
        outputs=[action_vector(action,y)for y in Y]
        components[name]=sy.Matrix([[sy.Rational(dot(out,v))for v in vectors]for out in outputs])
        counts[name]=[len(out)for out in outputs]
    expected=sy.Matrix([[0,-54,18,0,36,0,0],[0,0,0,36,0,18,54]])*sy.Rational(coefficients()['t'])
    assert components['base']==expected
    full=sum(components.values(),sy.zeros(2,7));GY=sy.diag(*[dot(y,y)for y in Y]);G=P.T*P
    M=full.T*GY.inv()*full/6
    a2=F(M[1,1]/G[1,1]);a=sqrt_interval(a2)[0];assert a>0
    assert F(dot(vectors[1],vectors[1])/(6*G[1,1]))==6
    # Initial current output spans<=1; every full H successor spans<=3.
    B=support_norm(3)*support_norm(1)*sqrt_upper(F(6))
    tau=a/B;amplitude=tau*a-tau*tau*B/2;assert amplitude>0
    cab=d['Cabibbo_squared'];probability=cab*amplitude**2/1296
    return {'identity':identity,'components':components,'counts':counts,'full':full,'GY':GY,'G':G,'Gram':M,
            'a2':a2,'a':a,'B':B,'tau':tau,'amplitude':amplitude,'probability':probability,'cab':cab}


def run():
    d=calculation();return {'schema':'oph.full_autonomous_native_strong_current.v1',
        'verdict':'PASS_COMPLETE_NATIVE_N5_STRONG_RESPONSE_WITH_RIGOROUS_FINITE_TIME_LOWER_BOUND',
        'actual_source':{'carrier':'Infinite-chain full colour-singlet Nq=4, Na=1 light sector with every relative configuration and CAR/Pauli output.',
            'native_signed_antifundamental_join_checked':True,'complete_100_triquark_CAR_rest_matrix_identity':True,
            'complete_48_meson_CAR_rest_matrix_and_spin_zero_intertwiner':True,
            'common_correlated_rule':'Sum_x Delta C_x* (zeta M_rest^-1) Delta C_x on the complete triquark and meson frames; all other family and antitriquark annihilators vanish on this reducing sector.',
            'triquark_rest_floor_source':str(d['identity']['triquark_floor']),'meson_rest_floor_source':str(d['identity']['meson_floor']),
            'triquark_inverse_block_dimensions':list(d['identity']['triquark_blocks']),'meson_inverse_block_dimensions':list(d['identity']['meson_blocks']),
            'nuclear_three_occupied_triplets_term_exact':'0 at total occupation5','infinite_Gauss_string_electric_energy_retained':True,
            'finite_volume_boundary_or_representation_cutoff':False},
        'literal_current_projection':{'source_Gram':record(d['G']),'separated_target_Gram':record(d['GY']),
            'component_projection_matrices':{k:record(v)for k,v in d['components'].items()},
            'all_component_successor_counts_on_two_targets':d['counts'],'full_projection_matrix':record(d['full']),
            'full_projected_derivative_Gram_without_Cabibbo':record(d['Gram']),
            'direct_current_to_disjoint_targets_exact':'0',
            'actual_input':'Normalized second native Lambda coordinate localized at cells0/1, not the Lambda rest eigenstate.',
            'targets':'Native onsite proton/pion singlet fields at cells(0,1) and(1,0); all orthogonal N5 outcomes retained.',
            'second_coordinate_derivative_norm_squared_over_Cabibbo_squared':str(d['a2']),
            'second_coordinate_derivative_norm_squared_display':enclosure((d['a2'],d['a2'])),
            'full_autonomous_H5_used_for_derivative':True},
        'rigorous_finite_time':{'source_time_exact':str(d['tau']),
            'first_derivative_amplitude_lower_over_Cabibbo':str(d['a']),
            'H_squared_current_vector_norm_upper_over_Cabibbo':str(d['B']),
            'remainder_proof':'Spectral theorem: ||(exp(-i tau H)-I+i tau H)u|| <= tau^2 ||H^2 u||/2; ||H^2 u|| <= M3 M1 sqrt(6) |Cabibbo|.',
            'H_squared_vector_explicitly_enumerated':False,'complete_support_operator_bounds_used':True,
            'span_one_H_norm_upper':str(support_norm(1)),'span_three_H_norm_upper':str(support_norm(3)),
            'separated_amplitude_lower_over_Cabibbo':str(d['amplitude']),
            'probe_then_autonomous_strong_separated_probability_lower':str(d['probability']),
            'probe_then_autonomous_strong_separated_probability_lower_diagnostic':float(d['probability']),
            'nonzero_complete_autonomous_response_proved':True,'probability_is_lifetime_or_on_shell_rate':False},
        'complete_records':{'probe_success_exact':str(d['cab']/216),'probe_failure_exact':str(1-d['cab']/216),
            'terminal_outcomes':['two disjoint native proton-pion fields','every orthogonal N5 output','probe failure'],
            'orthogonal_success_probability':'probe success minus separated event probability; no branch omitted.',
            'strong_evolution_mean_energy_work_exact':'0 under its own full autonomous H5.',
            'probe_and_terminal_readout_work':'Tr[H(Phi(rho)-rho)] on each complete instrument; not asserted zero.',
            'idle_external_quantum_reference_retained':True},
        'new_empirical_calibration':False,'new_match_with_nature':False,
        'inherited_calibrations':'Electron GeV anchor, strong/axial inputs, common inverse-rest propagation and Cabibbo angle.',
        'outgoing_Moller_wave_operator_or_on_shell_weak_spectral_density_proved':False,
        'physical_3plus1_decay_or_lifetime_proved':False,
        'scope':'An actual full autonomous infinite-chain N5 response after the native component probe, for a specified localized Lambda coordinate. This closes the kinetic-only omission; it does not identify an asymptotic scattering amplitude or a physical decay rate.'}
