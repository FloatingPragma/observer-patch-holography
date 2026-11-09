"""Finite primary modular work and explicit geometric attachment domains."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
import json
from pathlib import Path

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_stress import stress
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_channel_labels import channels
from model.particle_local_wilson_density import density
from model.particle_spatial_ground_state import ground
from model.particle_filtered_probe import filter as apparatus

SOURCE_INPUTS = (
    'model/particle_primary_gravity_premises/__init__.py',
    'model/particle_primary_gravity_premises/__main__.py',
    'model/particle_primary_gravity_premises/model.py',
    'model/particle_primary_gravity_premises/test_model.py',
    'model/particle_primary_gravity_premises/PROOF.md',
    'model/particle_primary_gravity_premises/README.md',
    'model/particle_primary_gravity_premises/premises.json',
    'model/particle_primary_stress/PROOF.md',
    'model/particle_spatial_modular_response/PROOF.md',
    'model/particle_field_state_response/PROOF.md',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_capacity/PROOF.md',
    'model/particle_geometric_response/PROOF.md',
    'model/matched_einstein_limit/PROOF.md',
    'Lean/QFT/CoarseEnergyWard.lean',
    'Lean/Thermodynamics/EinsteinPremiseLink.lean',
    'Lean/Geometry/SourceOrderEinsteinComposition.lean',
    'Lean/Geometry/EventFrameSoldering.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean',
)


def polynomial_slope(beta, spectral_lower, spectral_upper, order, scale):
    """Exact slope bracket for -log(T_m(beta(U-H)/(2s))**(2s)/Z)."""
    inherited=ground.polynomial_bound(beta,spectral_lower,spectral_upper,order,scale)
    m=inherited['Taylor_order']
    if m<1:raise ValueError('positive Taylor order required for a strict modular slope')
    y=inherited['scaled_argument_upper'];b=inherited['beta']
    t=sum((y**k/F(factorial(k)) for k in range(m+1)),F())
    delta=y**m/(factorial(m)*t)
    if not 0<=delta<1:raise ArithmeticError('positive polynomial modular slope failed')
    return {'beta':b,'Taylor_order':m,'dyadic_scale':inherited['dyadic_scale'],
            'scaled_argument_upper':y,'relative_modular_slope_defect_upper':delta,
            'modular_slope_lower':b*(1-delta),'modular_slope_upper':b,
            'global_modular_affine_defect_range_upper':b*delta*inherited['spectral_width_upper'],
            'global_modular_clock_half_trace_error_per_unit_parameter':b*delta*inherited['spectral_width_upper']/2,
            'compared_source_time':'minus beta times the modular parameter',
            'global_modular_clock_comparison_retains_arbitrary_idle_reference':True,
            'regional_or_proper_clock_identification':False,
            'source_state_replaced_by_exact_Gibbs':False,
            'exact_derivative':'K_prime(E)=beta*T_(m-1)(y(E))/T_m(y(E))'}


def selected_modular_slope(member):
    contract.verify_member(member)
    data=member['state']
    return polynomial_slope(data['beta'],data['spectral_lower'],data['spectral_upper'],
                            data['Taylor_order'],data['dyadic_scale'])


def local_entropy_work_bound(member, theta='1/1000'):
    """Bound the actual local source intervention using the inherited energy flux."""
    contract.verify_member(member)
    angle=ground._signed(theta,'source intervention angle')
    cell=stress.cell_energy(member);flux=stress.regional_bounds(member,((0,0,0),))
    k=F(cell['norm_upper']);c=F(flux['outward_flux_norm_upper'])
    beta=F(member['state']['beta']);x=beta*c/k
    adaptive=1/(4*k*(1+x));adaptive_cap=beta*k*c*adaptive*adaptive
    if adaptive_cap>F(1,64):raise ArithmeticError('adaptive regional entropy budget failed')
    return {'generator':'h_origin','generator_owner_support':cell['owner_support'],
            'generator_norm_upper':k,'commutator_with_full_H_norm_upper':c,
            'flux_energy_index_region':[[0,0,0]],
            'source_angle':angle,'cyclic_source_work_upper':angle*angle*k*c,
            'regional_modular_slope_upper':beta*angle*angle*k*c,
            'adaptive_source_angle':adaptive,'adaptive_work_upper':adaptive*adaptive*k*c,
            'adaptive_regional_modular_slope_upper':adaptive_cap,
            'adaptive_regional_modular_slope_simple_upper':F(1,64),
            'regional_entropy_plus_relative_entropy_upper':'t times the regional modular slope cap',
            'unitary_domain':'ambient tensor-local unitary commuting with full physical Gauss projection',
            'all_boundary_center_sectors_retained':True,
            'selected_input_state_unchanged':True,
            'regional_modular_generator_equals_beta_local_energy':False,
            'adaptive_angle_has_uniform_nonzero_refinement_limit':False}


@lru_cache(maxsize=1)
def flip_control():
    """Actual native block and spin chart; full-Fock bond cap from its spectrum."""
    native=channels.native_number_control()
    if next(row['internal_dimension'] for row in native['blocks'] if row['Q6']==6)!=1:
        raise ArithmeticError('q6 native block changed')
    gamma=density.gamma_data();beta=gamma[3];flip=gamma[0]
    if flip.H*flip!=s.eye(4) or flip.H*beta*flip!=-beta:
        raise ArithmeticError('the declared spin flip fails beta reversal')
    spectra=[]
    for g in gamma[:3]:
        spin=apparatus.clean((-beta+beta*g)/2)
        bond=s.zeros(8);bond[:4,4:]=spin;bond[4:,:4]=spin.H
        values=bond.eigenvals()
        if values!={s.Integer(-1):2,s.Integer(0):4,s.Integer(1):2}:
            raise ArithmeticError('the actual dimension-one Wilson bond spectrum changed')
        spectra.append({'one_body_eigenvalue_multiplicities':{'-1':2,'0':4,'1':2},
                        'full_Fock_bond_norm':2,'conjugated_difference_norm_upper':4})
    return {'native_integer_Q6':6,'native_dimension':1,
            'spin_operation':'gamma_1 on q6 at every cell; identity on other native blocks and links',
            'beta_reversal_checked':True,'five_native_occupations_preserved':True,
            'source_Gauss_and_CAR_parity_preserved':True,'directional_bond_controls':spectra}


def selected_flip_bound(member):
    contract.verify_member(member)
    if member['source']['side']!=2:raise ValueError('nonzero selected flip witness is certified at member zero')
    source=member['source'];a=F(source['spacing']);volume=source['side']**3
    bonds=3*source['side']**2*(source['side']-1)
    inherited=occupation.bind_member(member);occupation.verify_bound(member,inherited)
    defect=F(inherited['result']['defect_count_upper'])
    lower=(6*(2*volume-defect)-4*bonds)/a;upper=F(12*volume+4*bonds)/a
    if lower<=18:raise ArithmeticError('the selected primary flip lacks the strict work witness')
    slope=selected_modular_slope(member)
    return {'source_operation':flip_control(),
            'occupation_producer_receipt_sha256':contract.digest(inherited),
            'selected_total_beta_defect_upper':defect,
            'changed_hopping_norm_upper':F(4*bonds)/a,
            'electric_and_magnetic_energy_unchanged':True,
            'cyclic_source_work_lower':lower,'cyclic_source_work_simple_lower':18,
            'cyclic_source_work_upper':upper,'work_observable_norm_upper':upper,
            'source_modular_relative_entropy_lower':slope['modular_slope_lower']*lower,
            'source_modular_relative_entropy_upper':slope['modular_slope_upper']*upper,
            'complete_positive_work_record_probability_lower':(1+lower/upper)/2,
            'complete_positive_work_record_probability_simple_lower':F(9,16),
            'work_record_labels':[-upper,upper],
            'work_records_are_trajectory_work_spectrum':False,
            'work_record_error_rule':'subtract complete-output half-trace error from event probability',
            'full_retained_purifier_relative_entropy':'positive_infinity',
            'source_relative_entropy_is_full_purifier_relative_entropy':False,
            'large_primary_work_histogram_executed':False}


def boundary_obstruction(member):
    contract.verify_member(member)
    side=member['source']['side'];a=F(member['source']['spacing'])
    return {'side':side,'spacing':a,'missing_face_sites_per_direction':side**2,
            'periodic_coframe_formula_is_directly_primary_admissible':False,
            'one_body_flat_action_difference':'beta_D sum_i(I-T_i_dagger T_i)/(2a)',
            'full_Fock_and_physical_norm_exact':F(45*side*side)/a,
            'negative_extremizer':'filled negative-beta complete M15 rows, constant links',
            'positive_extremizer':'filled positive-beta complete M15 rows, constant links',
            'extremizers_Gauss_admissible':True,'extremizer_span_invariant_under_H_claimed':False,
            'product_link_compression_changes_boundary_term':False,
            'boundary_completed_coframe_or_proper_geometry_constructed_here':False}


def open_shift_control(side=2):
    """Exact flat one-body control of the periodic-to-open boundary mismatch."""
    side=ground._i(side,'control side',2)
    if side>3:raise ValueError('exact matrix control uses side at most three')
    cells=tuple(product(range(side),repeat=3));index={x:i for i,x in enumerate(cells)}
    volume=len(cells);identity=s.eye(volume);gamma=density.gamma_data();beta=gamma[3]
    primary=s.kronecker_product(identity,3*beta);naive=s.zeros(4*volume)
    boundary=s.zeros(4*volume);missing=0
    for axis in range(3):
        shift=s.zeros(volume)
        for cell in cells:
            target=list(cell);target[axis]+=1;target=tuple(target)
            if target in index:shift[index[cell],index[target]]=1
        K=identity-shift;pi=identity-shift.H*shift;missing+=s.trace(pi)
        spin=apparatus.clean((-beta+beta*gamma[axis])/2)
        primary+=s.kronecker_product(shift,spin)+s.kronecker_product(shift.H,spin.H)
        naive+=s.kronecker_product((shift-shift.H)/2,beta*gamma[axis])
        naive+=s.kronecker_product(K.H*K,beta)/2
        boundary+=s.kronecker_product(pi,beta)/2
    if primary-naive!=boundary or missing!=3*side**2:
        raise ArithmeticError('open primary boundary defect identity failed')
    filled_negative_expectation=-15*missing
    return {'side':side,'one_body_columns_per_native_dimension':4*volume,
            'missing_directional_face_slots':int(missing),
            'open_vs_periodic_Gram_identity_checked':True,
            'native_filled_negative_beta_expectation':int(filled_negative_expectation),
            'full_native_physical_norm_at_unit_spacing':45*side**2}


def premise_map():
    value=json.loads(Path(__file__).with_name('premises.json').read_text())
    if value['canonical_target_node']!='geometric_gravity':raise ValueError('gravity map target changed')
    for row in value['theorem_inputs']:
        for evidence in row['evidence']:
            path=contract.ROOT/evidence['path']
            if not path.is_file() or evidence.get('symbol','') not in path.read_text():
                raise ValueError('missing geometric premise evidence: '+str(path))
    return value


def selected_bound(member):
    contract.verify_member(member)
    return contract.serial({'member_id':member['member_id'],**member['ids'],
            'modular_work_slope':selected_modular_slope(member),
            'local_regional_entropy_work':local_entropy_work_bound(member),
            'nonzero_global_source_flip':selected_flip_bound(member),
            'open_coframe_boundary_obstruction':boundary_obstruction(member),
            'complete_coin_and_reference_retained':True,
            'entropy_domain':'source global or owner-region marginal as explicitly specified',
            'coin_plus_region_entropy':'H2(t)+S(sigma_R)',
            'coin_plus_full_source_reference_entropy':'H2(t)',
            'proper_geometry_or_Einstein_equation_proved':False,
            'gravity_required_for_finite_mass_readout':False})


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_gravity_premises.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('gravity premise bound or primary source changed')


def certificate():
    mapping=premise_map()
    evidence=tuple(item['path'] for row in mapping['theorem_inputs'] for item in row['evidence'])
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*evidence,*occupation.SOURCE_INPUTS,*stress.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_gravity_premises_certificate.v1','theorem_premise_map':mapping,
            'selected_member_bound':bind_member(contract.member(0)),
            'open_coframe_controls':[open_shift_control(side) for side in (2,3)],
            'gravity_is_separate_from_finite_mass_readout':True,'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('gravity inputs changed during replay')
    return result
