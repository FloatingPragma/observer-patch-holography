"""Original-weight link-dressed pair records on the unchanged primary state."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_refinement_response import model as refinement
from model.particle_primary_variational_readout import variational
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_mass_readout import readout
from model.particle_primary_channel_labels import channels
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_charged_string/__init__.py',
    'model/particle_primary_charged_string/__main__.py',
    'model/particle_primary_charged_string/model.py',
    'model/particle_primary_charged_string/test_model.py',
    'model/particle_primary_charged_string/PROOF.md',
    'model/particle_primary_charged_string/README.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_refinement_response/PROOF.md',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_endpoint_bounds/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
)


def _zero(member):
    contract.verify_member(member)
    if member['source']['level']!=0:
        raise ValueError('this charged-string occupation bound requires member zero')


@lru_cache(maxsize=1)
def retained_allocation():
    native=refinement.native_control();r=F(7,10);electric_retained=F(8,45)
    onsite=[1-(r if q==6 else 0) for d,C,q in native]
    shift=4*sum((d*(sqrt_upper(m*m+F(1,4))-m)
                   for (d,C,q),m in zip(native,onsite)),F())
    electric=105*(1-electric_retained)
    margin=F(3,20)*electric-shift
    if margin<=0:raise ArithmeticError('the full nonconstant-link block is not positive')
    rows=[]
    for (d,C,q),m in zip(native,onsite):
        inverse=1/(electric*C/2-shift)
        loss=(refinement.exact_q6_loss(m,inverse) if d==1
              else occupation.occupation_loss(d,m,inverse))
        rows.append({'dimension':d,'Q6':q,'Casimir':C,
                     'onsite_per_endpoint':m,'inverse_margin':inverse,**loss})
    loss=12*sum((F(x['occupation_loss_upper']) for x in rows),F())+F(6,7)
    k=48*r
    if k!=189*electric_retained:raise ArithmeticError('the union penalties do not match')
    return {'native_rows':rows,'retained_q6_onsite_fraction':r,
        'retained_electric_fraction':electric_retained,
        'pointwise_full_Fock_shift_upper':shift,
        'nonconstant_link_margin_lower':margin,
        'dimensionless_scalar_loss_upper':loss,'combined_bad_event_coefficient':k}


@lru_cache(maxsize=1)
def endpoint_control():
    """Exact endpoint charge and energy algebra on a finite quantum-link control."""
    c0=s.Matrix([[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]])
    c1=s.Matrix([[0,0,1,0],[0,0,0,-1],[0,0,0,0],[0,0,0,0]])
    nx=c0.H*c0;ny=c1.H*c1;pair=c0.H*c1
    P0=s.diag(0,1,0);shift=s.zeros(3);shift[2,1]=1
    E=s.diag(-1,0,1);B=s.kronecker_product(shift,pair)
    X=s.kronecker_product(s.eye(3),nx);Y=s.kronecker_product(s.eye(3),ny)
    electric=s.kronecker_product(E,s.eye(4))
    Gx=6*(X-electric);Gy=6*(Y+electric-s.eye(12))
    H0=3*(X-Y)+s.Rational(63,4)*electric**2
    def comm(a,b):return a*b-b*a
    if (comm(Gx,B)!=s.zeros(12) or comm(Gy,B)!=s.zeros(12)
            or comm(X,B)!=B or comm(Y,B)!=-B
            or comm(H0,B)!=s.Rational(87,4)*B
            or B.H*B!=s.kronecker_product(P0,ny*(s.eye(4)-nx))):
        raise ArithmeticError('charged endpoint or electric character identity failed')
    initial=s.eye(12)[:,6];target=B*initial
    if (target.H*target)[0]!=1 or any(G*v!=s.zeros(12,1)
            for G in (Gx,Gy) for v in (initial,target)):
        raise ArithmeticError('the endpoint control did not preserve physical Gauss support')
    return {'native_Q6':6,'weak_and_color_singlet':True,
        'endpoint_charge_increments':[6,-6],'net_charge_increment':0,
        'native_global_occupations_preserved':True,'CAR_parity':'even',
        'link_character_word_degree':1,'electric_character_cost':F(63,4),
        'onsite_pair_cost':6,'exact_reference_frequency':F(87,4),
        'control_dimension':12,'physical_source_histogram_executed':False}


def energy_bound(member):
    _zero(member);endpoint_control();channels.native_number_control()
    allocation=retained_allocation()
    upper=variational.bind_member(member);variational.verify_bound(member,upper)
    U=F(upper['result']['ground_energy_upper']);eps=F(member['state']['target_ground_excess'])
    bad=(U+eps+720+allocation['dimensionless_scalar_loss_upper'])/allocation['combined_bad_event_coefficient']
    weight=1-bad
    residual=sqrt_upper(F(40))+sqrt_upper(F(2))+F(4,7)
    if weight<=F(79,100) or residual**2>=70:
        raise ArithmeticError('the source charged-pair bound lost its margin')
    center=F(87,4);radius=F(12);band=(center-radius,center+radius)
    band_weight=weight-F(70)/radius**2
    if band_weight<=F(3,10):raise ArithmeticError('the positive charged-pair window failed')
    final=readout.excitation_window(band,F(3,10),eps,1)
    clocks=spectroscopy.clock_plan(member,'1/16','1/128','1/128')
    event_error=F(1,64)+F(clocks['two_complete_clock_error_upper'])
    recorded=F(3,10)-event_error
    if recorded<F(1,4):raise ArithmeticError('the complete charged energy record floor failed')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'allocation':allocation,'ground_upper_producer_sha256':upper['result_sha256'],
        'selected_mean_energy_upper':U+eps,
        'combined_bad_event':'D_Q6/16 + (5/18)*sum_link_Casimir',
        'combined_bad_event_mean_upper':bad,'original_success_weight_lower':weight,
        'selected_bond_hopping_residual_norm_upper':sqrt_upper(F(40)),
        'other_hopping_residual_norm_upper':sqrt_upper(F(2)),
        'magnetic_residual_norm_upper':F(4,7),
        'full_residual_norm_upper':residual,'centered_moment_upper':70,
        'exact_reference_frequency':center,'original_transfer_band':band,
        'original_transfer_weight_lower':band_weight,'original_transfer_simple_lower':F(3,10),
        'final_excitation':final,'complete_probe_outcomes':49,
        'successful_probe_count':48,'uniform_average_count':48,
        'complete_clock_plan':clocks,
        'reported_transfer_interval':[band[0]-F(1,8),band[1]+F(1,8)],
        'reported_transfer_original_weight_lower':recorded,
        'reported_transfer_simple_lower':F(1,4),
        'all_fine_energy_probe_records':49*int(clocks['clock_plan']['clock_size'])**2,
        'link_constant_projection_is_part_of_probe':True,
        'selected_preparation_replaced_by_constant_links':False,
        'original_weights_conditioned_on_link_success':False,
        'charge_record_intervention_has_its_own_disturbance':True,
        'energy_bound_transferred_to_charge_dephased_input':False,
        'source_reference_and_failures_retained':True,
        'physical_species_or_rest_mass_identified':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    value=energy_bound(member)
    return {'result':value,'result_sha256':contract.digest(value)}


def verify_bound(member,receipt):
    if receipt!=bind_member(member):raise ValueError('charged-string result or source changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_charged_string.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'endpoint_control':contract.serial(endpoint_control()),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('charged-string inputs changed during replay')
    return result
