"""Original selected-state odd response, exact CAR moments and refinement scale."""
from fractions import Fraction as F
from functools import lru_cache
import json

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_charge_conjugation import model as conjugation
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_refinement_response import model as refinement
from model.particle_primary_sector_threshold import model as sectors
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_mass_readout import readout
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
 'model/particle_primary_selected_odd_response/__init__.py',
 'model/particle_primary_selected_odd_response/__main__.py',
 'model/particle_primary_selected_odd_response/model.py',
 'model/particle_primary_selected_odd_response/test_model.py',
 'model/particle_primary_selected_odd_response/PROOF.md',
 'model/particle_primary_selected_odd_response/README.md',
 'model/particle_primary_charge_conjugation/PROOF.md',
 'model/particle_primary_occupation_bound/PROOF.md',
 'model/particle_primary_refinement_response/PROOF.md',
 'model/particle_primary_sector_threshold/PROOF.md',
 'model/particle_primary_joint_spectroscopy/PROOF.md',
 'model/particle_primary_current_tail/PROOF.md',
 'model/particle_common_contract/PROOF.md',
)


def creation_word(modes,count):
    """Exact full-CAR exterior creation, in the declared left-to-right order."""
    out=s.zeros(1<<count)
    for mask in range(1<<count):
        target=mask;sign=1
        for i in reversed(modes):
            if target&(1<<i):sign=0;break
            sign*=(-1)**((target&((1<<i)-1)).bit_count());target|=1<<i
        if sign:out[target,mask]=sign
    return out


@lru_cache(maxsize=1)
def exact_controls():
    """The six-CAR residual and the literal Wilson/Gauss identities."""
    R=s.zeros(64)
    for i in range(3):
        word=[0,1,2];word[i]=3+i
        R+=creation_word(word,6)
    G=R.H*R;I=s.eye(64)
    if G*(G-I)*(G-2*I)*(G-3*I)!=s.zeros(64):
        raise ArithmeticError('full-CAR residual polynomial failed')
    if G.trace()!=24 or (G-3*I).det()!=0:
        raise ArithmeticError('full-CAR residual spectrum changed')
    response=sectors.branch.response
    pairs=response.beta_pairs()
    plus=s.Matrix.hstack(pairs[0][0],pairs[1][0])
    g=response.density.gamma_data();beta=g[3]
    for gamma in g[:3]:
        hop=(-beta+beta*gamma)/2
        for h in (hop,hop.H):
            if response.apparatus.clean(plus.H*h.H*h*plus)!=s.eye(2)/2:
                raise ArithmeticError('positive-beta Wilson column Gram changed')
    C=conjugation.exact_controls()
    M=s.Matrix(C['actual_spin_matrix'])
    if response.apparatus.clean(beta*M*plus.conjugate()+M*plus.conjugate())!=s.zeros(4,2):
        raise ArithmeticError('positive occupation does not conjugate to negative holes')
    singlet=sectors.local_singlet((0,2,0,0,1))
    if singlet['norm_squared']!=1 or singlet['CAR_parity']!=1:
        raise ArithmeticError('native weak-pair charge-six monomial changed')
    density=response.density
    raw=json.loads((density.ROOT/density.SOURCE_INPUTS[-1]).read_text())
    masks=[sum(1<<j for j in raw['basis_subsets'][i]) for i in raw['even_indices']]
    ordered=[i for i in range(32) if i.bit_count() in (2,4)]
    permutation=[masks.index(mask) for mask in ordered]
    if permutation!=[0,1,4,2,5,7,10,3,6,8,11,9,12,13,14]:
        raise ArithmeticError('native Gauss and coordinate charts changed')
    z=(s.Integer(3)+4*s.I)/5
    u=s.Matrix([[s.Rational(3,5),4*s.I/5],[4*s.I/5,s.Rational(3,5)]])
    fundamental=s.diag(z**-2*s.diag(u,s.ones(1)),z**3*u)
    coordinate=density.native_representation(fundamental)
    native=coordinate.extract(permutation,permutation)
    blocks=dict(sectors.native_blocks());weak=blocks[-3];electron=blocks[6]
    if s.simplify(s.expand(native.extract(weak,weak).det()*native[electron[0],electron[0]]))!=1:
        raise ArithmeticError('actual nontrivial link does not preserve the odd singlet')
    return {'local_Gauss_kernel_dimension':1,'local_Gauss_generators_checked':singlet['Gauss_generators_checked'],
        'native_increment':[0,2,0,0,1],'native_Q6_order':[-4,-3,1,2,6],
        'local_creation_degree':3,'operator_CAR_parity':1,'operator_norm':1,
        'positive_beta_Wilson_column_Gram_per_neighbor':'I/2',
        'universal_residual_CAR_modes':6,'universal_residual_squared_norm':3,
        'universal_residual_squared_spectrum':[0,1,2,3],
        'universal_residual_Gram_trace':int(G.trace()),
        'C_maps_positive_number_to_negative_holes':True,
        'native_to_coordinate_basis_permutation':permutation,
        'nontrivial_coordinate_link_weak_determinant_times_singlet':1,
        'spin_choices_per_cell':4,'complete_rotation_multiplet_claimed':False,
        'pure_total_spin_assigned':False}


@lru_cache(maxsize=1)
def weak_e_allocation():
    """Full-source degree-six operator floor retaining only the needed rows."""
    retained=(F(),F(21,50),F(),F(),F(21,50))
    onsite=tuple((1-r)/2 for r in retained)
    delta=4*sum((d*(sqrt_upper(m*m+F(1,4))-m)
                 for (d,_,_),m in zip(refinement.native_control(),onsite)),F())
    if delta>=F(63,4):raise ArithmeticError('nonconstant link block lost positivity')
    rows=[]
    for (d,C,q),m,r in zip(refinement.native_control(),onsite,retained):
        w=1/(F(105,2)*C-delta)
        loss=refinement.exact_q6_loss(m,w) if d==1 else occupation.occupation_loss(d,m,w)
        rows.append({'Q6':q,'dimension':d,'retained_fraction':r,'onsite_per_endpoint':m,
                     'inverse_electric_margin':w,**loss})
    per_bond=sum((row['occupation_loss_upper'] for row in rows),F())
    seed=refinement.tile_seed()
    density=(seed['dimensionless_gain_per_cell']+3*per_bond+F(3,7)+F(1,128))/F(63,50)
    if not density<F(341,100):raise ArithmeticError('uniform selected weak-plus-e density failed')
    return {'maximum_vertex_degree':6,'retained_fractions':retained,
        'pointwise_full_Fock_shift_upper':delta,'nonconstant_link_margin_lower':F(63,4)-delta,
        'native_rows':rows,'per_bond_loss_upper':per_bond,'retained_defect_coefficient':F(63,50),
        'tiled_ground_upper':seed,'uniform_weak_plus_e_defect_density_upper':density,
        'uniform_simple_density_upper':F(341,100)}


def selected_response(member):
    contract.verify_member(member);control=exact_controls();allocation=weak_e_allocation()
    source=member['source'];n=source['level'];a=F(source['spacing']);V=source['side']**3
    eps=F(member['state']['target_ground_excess']);tile=refinement.tiled_upper(member)
    numerator=tile['dimensionless_ground_upper']+90*V
    numerator+=source['interior_hopping_bonds']*allocation['per_bond_loss_upper']
    numerator+=F(source['interior_plaquettes'],7)+a*eps
    defects=min(F(12*V),numerator/allocation['retained_defect_coefficient'])
    evidence={'weak_plus_e_retained_bound':defects}
    if n==0:
        value=occupation.bind_member(member);occupation.verify_bound(member,value)
        evidence['total_native_defect_bound']=F(value['result']['defect_count_upper'])
        evidence['occupation_receipt_sha256']=value['result_sha256']
        defects=min(defects,evidence['total_native_defect_bound'])
    success=max(F(),1-defects/(4*V))
    uniform_success=F(59,400)
    if success<uniform_success:raise ArithmeticError('selected original success floor failed')
    degree=3 if n==0 else 6
    moment=F(3*degree,2)/a**2;center=9/a
    lower=F(1,2)/a;upper=F(35,2)/a
    tail=moment/(center-lower)**2;weight=success-tail
    if weight<F(1,50):raise ArithmeticError('positive original cutoff-scale band failed')
    plan=spectroscopy.clock_plan(member,lower/4,'1/4096','1/4096')
    clock_loss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    reported=weight-clock_loss
    if reported<F(1,64):raise ArithmeticError('complete clock record floor failed')
    sharp=None
    if n==0:
        sharp_weight=success-moment/36
        if sharp_weight<F(23,32):raise ArithmeticError('member-zero odd band failed')
        sharp={'original_transfer_interval':[F(3),F(15)],
            'original_transfer_weight_lower':sharp_weight,'simple_original_weight_lower':F(23,32),
            'final_excitation_transport':readout.excitation_window((F(3),F(15)),sharp_weight,eps,F(1))}
    return contract.serial({'schema':'oph.primary_selected_odd_response_result.v1',
        'member_id':member['member_id'],**member['ids'],'source_controls':control,
        'weak_e_allocation':allocation,'selected_defect_evidence':evidence,
        'weak_plus_e_defect_count_upper':defects,'positive_beta_selected_count_upper':defects/2,
        'original_creation_success_lower':success,'uniform_simple_creation_success_lower':uniform_success,
        'complete_instrument':'uniform retained cell and two spin choices; charged odd memory exchange; success and failure',
        'probe_operators':'A_(x,alpha,b)=c+_(weak1,alpha)^dagger c+_(weak2,alpha)^dagger c+_(e,b)^dagger',
        'complete_address_spin_count':4*V,'complete_address_spin_success_failure_outcomes':8*V,
        'source_native_increment':[0,2,0,0,1],'compensating_memory_native_increment':[0,-2,0,0,-1],
        'total_source_memory_native_numbers_preserved':True,'total_source_memory_CAR_parity_preserved':True,
        'local_gauge_charge_zero':True,'odd_memory_prepared_as_apparatus_input':True,
        'source_marginal_replaced_by_trial':False,'all_failures_and_quantum_reference_retained':True,
        'energy_transfer_measure_uses_original_probability':True,'division_by_success_probability':False,
        'centered_second_moment_center':center,'centered_second_moment_upper':moment,
        'original_transfer_interval':[lower,upper],'original_transfer_weight_lower':weight,
        'uniform_simple_original_weight_lower':F(1,50),'complete_clock_plan':plan,
        'complete_clock_probability_error_upper':clock_loss,
        'finite_record_transfer_interval':[lower/2,upper+lower/2],
        'finite_record_original_weight_lower':reported,'uniform_simple_record_weight_lower':F(1,64),
        'final_excitation_lower_with_same_original_weight':lower,
        'member_zero_sharper_band':sharp,
        'charge_conjugate_hole_probe_has_identical_complete_energy_law':True,
        'charge_conjugate_source_increment':[0,-2,0,0,-1],
        'unscaled_odd_transfer_family_is_uniformly_energy_tight':False,
        'uniform_cutoff_scale_response_proved':True,
        'some_weight_escapes_every_fixed_upper_source_energy':True,
        'all_fixed_energy_weight_vanishes_proved':False,
        'critical_density_Weyl_survival_decided':False,
        'density_Weyls_preserve_all_native_numbers_and_CAR_parity':True,
        'fixed_physical_band_requires_dimensionless_difference_order_a':True,
        'source_reference_clock_changed':False,'exact_Gibbs_replacement':False,
        'isolated_rest_particle_or_Lorentz_spin_proved':False,'physical_mass_table_GeV':[]})


def bind_member(member):
    value=selected_response(member)
    return {'schema':'oph.primary_selected_odd_response.v1','result':value,'result_sha256':contract.digest(value)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('selected odd response or its actual source changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*occupation.SOURCE_INPUTS,*refinement.SOURCE_INPUTS,
        *conjugation.SOURCE_INPUTS,*sectors.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS,*readout.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_selected_odd_response_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],'source_sha256':pins,
        'fixed_finite_energy_survival_proved':False,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during replay')
    return value
