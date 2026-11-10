"""Exact full-source graded responses and complete final-energy witnesses."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract,provenance
from model.particle_candidate_composites import model as identities
from model.particle_primary_composite_spin import model as spin
from model.particle_primary_composite_transport import model as transport
from model.particle_primary_odd_sum_rules import model as odd
from model.particle_primary_sector_threshold import model as sectors
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_spatial_modular_response import regional

SOURCE_INPUTS=tuple(dict.fromkeys((
 'model/particle_candidate_composite_spectrum/__init__.py',
 'model/particle_candidate_composite_spectrum/__main__.py',
 'model/particle_candidate_composite_spectrum/model.py',
 'model/particle_candidate_composite_spectrum/test_model.py',
 'model/particle_candidate_composite_spectrum/PROOF.md',
 'model/particle_candidate_composite_spectrum/README.md',
 *identities.SOURCE_INPUTS,*spin.SOURCE_INPUTS,*transport.SOURCE_INPUTS,
 *odd.SOURCE_INPUTS,*sectors.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS,*regional.SOURCE_INPUTS,
)))


def clean_state(vector):
    return {mask:v for mask,value in vector.items() if (v:=s.expand(value))!=0}


def add(out,vector,multiplier=1):
    for mask,value in vector.items():out[mask]=out.get(mask,0)+multiplier*value


def apply_word(vector,word):
    out={}
    for mask,value in vector.items():
        target,phase=spin.apply_operators(mask,word)
        if phase:out[target]=out.get(target,0)+phase*value
    return clean_state(out)


def apply_polynomial(vector,polynomial,site,adjoint=False):
    out={}
    for word,value in polynomial:
        if adjoint:word=tuple((i,not create) for i,create in reversed(word));value=s.conjugate(value)
        shifted=tuple((60*site+i,create) for i,create in word)
        add(out,apply_word(vector,shifted),value)
    return clean_state(out)


def norm_squared(vector):return s.simplify(sum(s.conjugate(c)*c for c in vector.values()))


@lru_cache(maxsize=None)
def irrep_polynomials(delta):
    """Literal cube projectors on the actual local gauge-singlet coefficient space."""
    delta=tuple(delta);local=spin.local_multiplet(delta);cube=spin.cube_restriction(delta)
    result=[]
    for label in cube['proper_binary_cube_irreps']:
        # Use the identity quaternion explicitly rather than the enumeration order.
        identity=next(q for q in spin.group.elements() if tuple(map(spin.channels._sym,q))==(1,0,0,0))
        dim=spin.channels._sym(spin.group.characters(identity)[label])
        P=spin.clean(dim*sum((s.conjugate(spin.channels._sym(spin.group.characters(q)[label]))*D
                             for q,D in zip(spin.group.elements(),cube['matrices'])),s.zeros(local['kernel_basis'].cols))/48)
        if P.H!=P or spin.clean(P*P)!=P:raise ArithmeticError('actual cube coefficient projector failed')
        W=s.Matrix.hstack(*s.GramSchmidt(P.columnspace(),orthonormal=True))
        if W.cols!=dim:raise ArithmeticError('the local irrep multiplicity changed')
        for j in range(W.cols):
            polynomial={}
            for k in range(W.rows):
                for word,value in spin.composite_polynomial(delta,k):
                    polynomial[word]=polynomial.get(word,0)+W[k,j]*value
            polynomial=tuple((word,s.expand(value)) for word,value in polynomial.items() if s.expand(value))
            if s.simplify(sum(s.conjugate(value)*value for _,value in polynomial))!=1:
                raise ArithmeticError('the irrep polynomial is not normalized on the physical anchor')
            result.append((label,j,polynomial))
    return tuple(result)


def derivative_coefficient(vector,polynomial,charge,left,right,conjugate):
    """One actual beta-preserving Wilson-link coefficient of [H,F_y]."""
    block=dict(sectors.native_blocks())[charge];out={}
    for word,value in polynomial:
        for position,(mode,create) in enumerate(word):
            row,native=divmod(mode,15)
            if native!=block[right] or (not create)!=conjugate:continue
            ops=[(60+i,c) for i,c in word]
            ops[position]=(15*row+block[left],create)
            add(out,apply_word(vector,tuple(ops)),-value/2)
    return clean_state(out)


def local_defects(mask,site):
    bits=(mask>>(60*site))&((1<<60)-1)
    return (bits^spin.OMEGA).bit_count()


@lru_cache(maxsize=None)
def graded_witness(delta,label,component):
    delta=tuple(delta);polynomial=next(p for ell,j,p in irrep_polynomials(delta) if (ell,j)==(label,component))
    t=len(polynomial);degree=sum(map(abs,delta));parity=degree%2
    omega=spin.OMEGA|(spin.OMEGA<<60)
    witness=apply_polynomial({omega:s.S.One},polynomial,0)
    if norm_squared(witness)!=1:raise ArithmeticError('the physical source witness lost norm one')
    annihilated=apply_polynomial(witness,polynomial,0,True)
    if annihilated!={omega:s.S.One}:raise ArithmeticError('the source polynomial adjoint did not recover its anchor')
    casimirs={q:F(C) for _,C,q in sectors.refinement.native_control()}
    rows=[];total=s.S.Zero;electric=s.S.Zero
    for (charge,block),increment in zip(sectors.native_blocks(),delta):
        if not increment:continue
        conjugate=increment<0;dimension=len(block)
        for left in range(dimension):
            for right in range(dimension):
                first=apply_polynomial(derivative_coefficient(witness,polynomial,charge,left,right,conjugate),polynomial,0,True)
                second=derivative_coefficient(annihilated,polynomial,charge,left,right,conjugate)
                add(first,second,-(-1)**parity);output=clean_state(first)
                for mask in output:
                    if (local_defects(mask,0),local_defects(mask,1))!=(1,degree-1):
                        raise ArithmeticError('graded output escaped the actual fixed-defect distribution')
                n=norm_squared(output)/dimension;total+=n
                electric+=n*s.Rational(105,2)*s.Rational(casimirs[charge].numerator,casimirs[charge].denominator)
                rows.append({'Q6':charge,'conjugate_link_coefficient':conjugate,
                    'left':left,'right':right,'native_dimension':dimension,
                    'matter_norm_squared':F(norm_squared(output)),
                    'Haar_norm_squared':F(n),'source_link_Casimir':casimirs[charge]})
    g=F(s.simplify(total));energy=F(-720+3*degree)+F(s.simplify(electric/total))
    if g<=0:raise ArithmeticError('the actual graded spectral witness vanished')
    return {'delta':delta,'probe_irrep':label,'component':component,
        'raw_polynomial_monomials':t,'raw_graded_witness_norm_squared':g,
        'contractive_probe_graded_witness_norm_squared':g/(t*t),
        'normalized_first_action_source_energy':energy,
        'one_link_electric_mean':energy+720-3*degree,
        'fixed_total_defect_count':degree,'output_defects_at_x_y':[1,degree-1],
        'source_Haar_coefficients':rows,'graded_commutator_sign':-(-1)**parity,
        'physical_link_words_admitted':True,
        'source_input_replaced_by_witness':False}


def spectral_window(lower,mean,ground_lower,ground_upper,plan):
    """Markov bound on the actual normalized first-action witness, before transfer to rho."""
    lower,mean,ground_lower,ground_upper=map(F,(lower,mean,ground_lower,ground_upper))
    eta=F(plan['energy_resolution']);loss=F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['one_complete_clock_error_upper'])
    if mean<=lower or lower-eta<=ground_upper or loss>=F(1,6):
        raise ValueError('source sector, first moment or complete clock bound is insufficient')
    upper=3*mean-2*lower
    return {'physical_sector_energy_lower':lower,'normalized_first_action_energy_mean':mean,
        'ideal_final_source_energy_interval':[lower,upper],
        'first_action_ideal_interval_weight_lower':F(2,3),
        'reported_final_source_energy_interval':[lower-eta,upper+eta],
        'reported_final_excitation_interval':[lower-eta-ground_upper,upper+eta-ground_lower],
        'one_complete_clock_tail_and_implementation_loss':loss,
        'first_action_reported_interval_weight_lower':F(2,3)-loss,
        'first_action_reported_interval_simple_lower':F(1,2),
        'energy_label_is_transfer':False,'initial_energy_conditioning_used':False}


def finite_time_plan(norm,density_exponent,gamma):
    R=F(norm);N=int(density_exponent);gamma=F(gamma)
    if R<1 or N!=density_exponent or N<0 or not 0<gamma<=1:
        raise ValueError('actual norm, exact density exponent and positive witness required')
    t=gamma/(64*R*R)
    near=gamma*t*t/16;far=9*R**4*t**4
    if far/near!=F(9,256)*gamma or 3*R*R*t>F(3,64):
        raise ArithmeticError('original spectral contrast constants changed')
    return {'source_norm_upper':R,'density_floor_negative_exponent':N,'common_witness_lower':gamma,
        'time':{'base':2,'negative_exponent':N,'rational_prefactor':t},
        'source_time_formula':'t=c*gamma/(64*R^2), c=2^(-N)',
        'same_source_clock_formula':'U_t=(I+i*t*H/2)^(-1)*(I-i*t*H/2)',
        'neighbor_original_joint_event_probability_lower':{'base':2,'negative_exponent':3*N,'rational_prefactor':near},
        'nonneighbor_original_joint_event_probability_upper':{'base':2,'negative_exponent':4*N,'rational_prefactor':far},
        'neighbor_minus_nonneighbor_original_joint_probability_lower':{'base':2,'negative_exponent':3*N,'rational_prefactor':near/2},
        'nonneighbor_to_neighbor_bound_ratio_upper':F(9,256)*gamma,
        'retained_uniform_setting_bit_joint_weight_multiplier':F(1,2),
        'all_forty_components_uniform_setting_joint_weight_multiplier':F(1,40),
        'uniform_in_capacity':False,'source_energy_zero_changed':False,
        'exponential_size_integer_or_full_source_matrix_allocated':False}


@lru_cache(maxsize=1)
def calculated_rows():
    rows=[]
    for family,delta in identities.FAMILIES:
        for sign in (1,-1):
            actual=tuple(sign*x for x in delta)
            for label,j,_ in irrep_polynomials(actual):
                rows.append({'family':family,'sign':sign,**graded_witness(actual,label,j)})
    if len(rows)!=40:raise ArithmeticError('the six signed local kernels do not have forty components')
    return tuple(rows)


def bind_member(member):
    contract.verify_member(member);sectors.branch._member_zero(member)
    lower=sectors.bind_member(member);sectors.verify_bound(member,lower)
    source=lower['result'];floor_payload=regional.common_contract_payload(member,1)
    floor=floor_payload['physical_density_floor_factored']
    if floor['base']!=2 or floor['operator']!='P_G':raise ArithmeticError('the original physical density floor changed')
    plan=spectroscopy.clock_plan(member,'1/64','1/1024','1/1024')
    floors={tuple(row['delta']):F(row['physical_sector_energy_lower'])
            for family in source['sector_results'] for row in family['signed_sectors']}
    rows=[]
    for raw in calculated_rows():
        row=dict(raw);delta=raw['delta']
        row['actual_final_native_occupation']=[x+y for x,y in zip(source['anchor_native_occupations'],delta)]
        row['complete_final_energy_event']=spectral_window(floors[delta],raw['normalized_first_action_source_energy'],
                source['global_ground_lower'],source['global_ground_upper'],plan)
        rows.append(row)
    gamma=min(F(1),*(r['contractive_probe_graded_witness_norm_squared'] for r in rows))
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'native_charge_order':[-4,-3,1,2,6],
        'sector_threshold_bound_sha256':contract.digest(lower),
        'physical_floor_payload_sha256':contract.digest(floor_payload),
        'physical_density_floor_factored':floor,'component_rows':rows,
        'complete_final_energy_clock':plan,
        'common_original_spatial_energy_contrast':finite_time_plan(plan['source_norm_upper'],floor['negative_exponent'],gamma),
        'complete_instrument':{'graded_success':'K=(A_x* U_t* A_y U_t - (-1)^p U_t* A_y U_t A_x*)/2',
            'failure':'sqrt(I-K* K)','intervention_outcomes':2,
            'followup_records':['all five final native occupations','every final energy-clock outcome'],
            'final_native_occupation_projectors_complete':True,
            'net_native_increment':[0]*5,'CAR_parity':0,'gauge_invariant':True,
            'original_input_density_replaced':False,'nonselective_output_density_claimed_unchanged':False,
            'all_failed_interventions_other_sectors_and_energy_records_retained':True,
            'quantum_reference_retained':True,'source_marginal_Loewner_floor_only':True,
            'scalar_event_probability_bound_is_not_joint_reference_domination':True},
        'nonzero_original_joint_spatial_final_energy_contrast_proved':True,
        'original_odd_neighboring_two_point_m1_value_or_sign_computed':False,
        'source_or_observer_conserved_momentum_identified':False,
        'physical_neutrino_or_neutron_identified':False,
        'isolated_line_or_relativistic_dispersion_proved':False,
        'named_physical_particles_identified':False,'physical_particle_map':[],
        'mass_predictions_GeV':[]})
    return {'schema':'oph.particle_candidate_composite_spectrum_member.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('composite spectral source, component or complete event changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.particle_candidate_composite_spectrum.v1',
        'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source inputs changed during composite replay')
    return value


def build():return certificate()
