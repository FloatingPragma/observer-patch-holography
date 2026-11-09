"""Faithful spatial restrictions, Gauss centers and actual native response.

The spatial family is the full Wilson product-cutoff field. Exact executed
controls use the208-state source's internal spinor/flux factors; those cuts
are explicitly not spatial balls. Regional matrix logarithms are enclosed
by rational polynomials, rather than assigned a geometric kernel.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as s

from model.particle_field_state_response import state_response as thermal
from model.particle_spatial_ground_state import ground
from model.particle_geometric_response import geometry
from model.native_gravity_identification import entropy as inherited_entropy
from model.particle_observer_outputs import field_route

ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_gauss_ground_state/gauss.py',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_spatial_ground_state/ground.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_field_state_response/state_response.py',
    'model/particle_field_state_response/PROOF.md',
    'model/particle_geometric_response/geometry.py',
    'model/particle_regional_geometry/PROOF.md',
    'model/prepared_regional_entropy/PROOF.md',
    'model/native_gravity_identification/entropy.py',
    'Lean/QFT/CPRestrictionNet.lean',
    'Lean/Thermodynamics/EinsteinPremiseLink.lean',
    'Lean/Thermodynamics/CapFirstLaw.lean',
)


def integer(value, name, minimum=0):
    return ground._i(value, name, minimum)


def rational(value, name, positive=False):
    return ground._q(value, name, positive)


def clean(matrix):
    return geometry.clean(matrix)


def fraction(value):
    value = s.Rational(value)
    return F(int(value.p), int(value.q))


def row_norm(matrix):
    return max(sum(abs(matrix[i,j]) for j in range(matrix.cols)) for i in range(matrix.rows))


def spatial_chart(mask, link_labels, cells, kept_cells):
    """Actual full60CAR/cell permutation and outgoing-link ownership.

    This computes basis labels/signs without allocating the Fock carrier.
    Word labels refer to any fixed orthonormal basis of the native V_ell.
    """
    cells=integer(cells,'cell count',1)
    if (type(mask) is not int or not 0 <= mask < (1 << (60*cells))
            or type(kept_cells) is not tuple or not kept_cells
            or tuple(sorted(set(kept_cells))) != kept_cells
            or any(type(i) is not int or not 0<=i<cells for i in kept_cells)
            or len(link_labels) != 3*cells
            or any(type(i) is not int or i<0 for i in link_labels)):
        raise ValueError('complete source labels and an ordered nonempty cell region required')
    kept=tuple(i for x in kept_cells for i in range(60*x,60*(x+1)))
    outside=tuple(i for x in range(cells) if x not in kept_cells
                  for i in range(60*x,60*(x+1)))
    left=sum(((mask>>i)&1)<<j for j,i in enumerate(kept))
    right=sum(((mask>>i)&1)<<j for j,i in enumerate(outside))
    inversions=sum(((mask>>i)&1)*((mask>>j)&1) for i in kept for j in outside if j<i)
    return {'inside_CAR_mask':left,'outside_CAR_mask':right,'CAR_sign':(-1)**inversions,
            'inside_mode_order':kept,'outside_mode_order':outside,
            'inside_link_words':tuple(link_labels[3*x+j] for x in kept_cells for j in range(3)),
            'outside_link_words':tuple(link_labels[3*x+j] for x in range(cells)
                                       if x not in kept_cells for j in range(3))}


def spatial_faithfulness(side=2, link_degree=1, region_cells=1, beta='1',
                         electric='105', magnetic='1/7', spacing='1'):
    """Spectral-width bound on the actual spatial Gibbs or positive polynomial state.

    Floors are factored dyadics to avoid allocating huge integers. They are
    finite-volume bounds, with no uniform low-temperature limit assertion.
    """
    source=ground.spatial_source_bounds(side=side,link_degree=link_degree,
        electric=electric,magnetic=magnetic,spacing=spacing,regulator='wilson')
    r=integer(region_cells,'region cell count',1)
    if r>source['cells']:raise ValueError('region exceeds the source box')
    b=rational(beta,'inverse temperature',True)
    width=source['spectral_width_upper']
    exponent=(2*b*width).__ceil__()
    regional_bits=(60+27*link_degree)*r
    return {'source':source,'beta':b,'region_cells':r,
            'regional_dimension_bits_upper':regional_bits,
            'ambient_marginal_minimum_factored':{'base':2,'negative_exponent':exponent+regional_bits},
            'ambient_full_minimum_factored':{'base':2,'negative_exponent':exponent+source['dimension_bits_upper']},
            'common_Boltzmann_floor_factored':{'base':2,'negative_exponent':exponent},
            'Gibbs_and_same_positive_Taylor_state_bound':True,
            'physical_Gauss_projection_preserves_ambient_full_rank':False,
            'floor_uniform_in_volume_or_zero_temperature_claimed':False}


def common_contract_payload(member,region_cells=1):
    """Authenticate the selected Gauss state and compute its regional support floor.

    Import is local so the contract can register this consumer without an
    import cycle. No warm native control or ambient full-rank bound enters.
    """
    from model.particle_common_contract import contract
    contract.verify_member(member)
    source=member['source'];state=member['state'];plan=member['preparation_plan']
    count=integer(region_cells,'region cell count',1)
    if count>source['side']**3:raise ValueError('region exceeds the declared box')
    beta=F(state['beta']);width=F(state['spectral_upper'])-F(state['spectral_lower'])
    exponent=(2*beta*width).__ceil__()+plan['source']['dimension_bits_upper']
    result={'region_cells':count,'beta':beta,'Taylor_order':state['Taylor_order'],
        'dyadic_scale':state['dyadic_scale'],'source_only_Gauss_support':True,
        'physical_density_floor_factored':{'base':2,'negative_exponent':exponent,
                                         'operator':'P_G'},
        'regional_density_floor_factored':{'base':2,'negative_exponent':exponent,
                                         'operator':'Tr_out(P_G)'},
        'modular_domain':'support of Tr_out(P_G), with physical boundary center retained',
        'regional_full_ambient_faithfulness_claimed':False,
        'boundary_multiplicities_computed_at_full_spatial_size':False,
        'large_regional_density_allocated':False,
        'warm_native_control_used_as_selected_state':False,
        'physical_geometric_kernel_or_area_identified':False}
    return serialize(result)


def common_contract_bound(member,region_cells=1):
    """Bind the registered one-cell consequence to the exact selected member."""
    from model.particle_common_contract import contract
    return contract.bind_result(member,
        'faithful regional support of the selected source-Gauss polynomial state',
        common_contract_payload(member,region_cells))


def boundary_blocks(sectors, beta_width='0'):
    """Observable-algebra and extended-state floors for supplied actual sectors.

    Each row is (label, inside multiplicity, outside multiplicity, irrep dim).
    Only the representation decomposition of the specified boundary action
    may supply these data. This function does not infer them from dimension.
    """
    if not sectors:raise ValueError('nonempty boundary decomposition required')
    checked=[];labels=set()
    for label,m,n,d in sectors:
        if label in labels:raise ValueError('duplicate boundary sector')
        labels.add(label)
        checked.append((label,integer(m,'inside multiplicity',1),
                        integer(n,'outside multiplicity'),integer(d,'irrep dimension',1)))
    rank=sum(m*n for _,m,n,_ in checked)
    if not rank:raise ValueError('the physical source support is empty')
    bw=rational(beta_width,'beta times spectral width')
    exponent=(2*bw).__ceil__()
    rows=[]
    for label,m,n,d in checked:
        rows.append({'label':label,'inside_multiplicity':m,'outside_multiplicity':n,
            'boundary_irrep_dimension':d,'supported':n>0,
            'partial_trace_physical_projector_eigenvalue':F(n,d),
            'extended_marginal_floor_factor':F(n,rank*d),
            'observable_block_floor_factor':F(n,rank),
            'common_floor_dyadic_negative_exponent':exponent,
            'central_edge_entropy_coefficient':'log('+str(d)+')'})
    return {'physical_rank':rank,'rows':rows,
            'regional_observable_algebra':'direct_sum_lambda B(M_R_lambda), with every boundary center weight',
            'entropy_relation':'S_extended=S_observable+sum_lambda p_lambda log d_lambda',
            'source_region_assumed_to_be_separately_Gauss_singlet':False,
            'edge_multiplicity_identified_with_geometric_area':False}


def native_boundary_pair():
    """Actual M15/conjugate-M15 singlets test the boundary multiplicities.

    This is a kinematic source-representation control, not the spatial
    Wilson prepared state. It imports the native charge blocks unchanged.
    """
    source=field_route.source_data.load();p=s.zeros(225);sectors=[]
    for block in source['blocks']:
        indices=block['indices'];d=block['dimension']
        sectors.append((str(block['q6Y']),1,1,d))
        for i in indices:
            for j in indices:p[15*i+i,15*j+j]=s.Rational(1,d)
    p=clean(p)
    if p*p!=p or s.trace(p)!=5:
        raise ArithmeticError('the native boundary singlet projector differs')
    marginal=partial_trace(p,15,15)
    expected=s.diag(*(s.Rational(1,next(b['dimension'] for b in source['blocks']
                                      if i in b['indices'])) for i in range(15)))
    if marginal!=expected:
        raise ArithmeticError('boundary multiplicities lost their irrep dimension')
    return {'block_data':boundary_blocks(sectors),'partial_projector_diagonal':list(marginal.diagonal()),
            'physical_singlet_rank':5,'full_pair_dimension':225,
            'observable_algebra':'C direct_sum C direct_sum C direct_sum C direct_sum C',
            'extended_regional_dimension':15,'spatial_Wilson_state_substituted':False}


def partial_trace(matrix, left, right, keep='left'):
    left=integer(left,'left dimension',1);right=integer(right,'right dimension',1)
    if matrix.shape != (left*right,left*right) or keep not in ('left','right'):
        raise ValueError('complete two-factor source matrix required')
    result=s.zeros(left if keep=='left' else right)
    for (i,j),v in matrix.todok().items():
        if keep=='left' and i%right==j%right:result[i//right,j//right]+=v
        if keep=='right' and i//right==j//right:result[i%right,j%right]+=v
    return clean(result)


def log_polynomial(sigma, order=12):
    """K=-log(sigma)=log(d)I+P, with rational ||P-P_n|| bound."""
    order=integer(order,'log order',1);sigma=clean(sigma)
    if sigma.rows!=sigma.cols or sigma!=sigma.H or s.trace(sigma)!=1:
        raise ValueError('a normalized Hermitian regional density required')
    d=sigma.rows;y=clean(s.eye(d)-d*sigma)
    q=fraction(row_norm(y))
    if q>=1:raise ValueError('this direct log series requires a certified contraction')
    power=s.eye(d);poly=s.zeros(d)
    for j in range(1,order+1):power=power*y;poly+=power/j
    tail=q**(order+1)/((order+1)*(1-q))
    majorant=sum(q**j/j for j in range(1,order+1))
    return {'matrix':clean(poly),'order':order,'contraction_norm_upper':q,
            'operator_error_upper':tail,'scalar_part':'log('+str(d)+')',
            'polynomial_norm_upper':majorant}


def fermion_swap(first=0,second=2):
    """Second-quantized mode permutation with its CAR signs, all16 masks."""
    if first==second or any(type(x)is not int or not 0<=x<4 for x in (first,second)):
        raise ValueError('two distinct native spinor modes required')
    permutation=list(range(4));permutation[first],permutation[second]=second,first
    result=s.zeros(16)
    for mask in range(16):
        occupied=[permutation[i] for i in range(4) if mask>>i&1]
        sign=(-1)**sum(occupied[i]>occupied[j] for i in range(len(occupied)) for j in range(i+1,len(occupied)))
        result[sum(1<<i for i in occupied),mask]=sign
    return clean(result)


@lru_cache(maxsize=2)
def native_state(physical=False):
    """The inherited warm rational preparation, optionally source-Gauss projected."""
    if type(physical)is not bool:raise ValueError('physical flag must be boolean')
    prep=thermal.prepared_density(order=12,beta='1/64');rho=prep['rho']
    masks=tuple(m for m in range(16) if not physical or m.bit_count()==2)
    indices=tuple(13*m+f for m in masks for f in range(13))
    kept=rho.extract(indices,indices);weight=s.trace(kept)
    kept=clean(kept/weight)
    raw_normalizer=s.factor(prep['norm_squared']*weight)
    h=thermal.exact_control()['H'].extract(indices,indices)
    if clean(kept*h-h*kept)!=s.zeros(len(indices)):
        raise ArithmeticError('the source preparation lost stationarity')
    return {'rho':kept,'H':h,'masks':masks,'indices':indices,
            'raw_Q_squared_normalizer':raw_normalizer,
            'minimum_source_weight':1/raw_normalizer,
            'source_Gauss_zero':physical,'dimension':len(indices),
            'original_source_dimension':208,'physical_dimension':78 if physical else None,
            'beta':F(1,64),'Taylor_order':12,
            'reference_retained':True,'postselection':False}


def modular_correction(poly, bare, error):
    """A trace-zero matrix-element witness against beta*bare plus a scalar."""
    residual=clean(poly-bare);d=residual.rows
    candidates=[]
    for i in range(d):
        for j in range(i+1,d):
            candidates.append((abs(residual[i,i]-residual[j,j]),'diagonal_difference',i,j,2))
            candidates.append((abs(residual[i,j]),'offdiagonal_entry',i,j,1))
    amount,kind,i,j,cost=max(candidates,key=lambda row:row[0])
    return {'witness_kind':kind,'indices':[i,j],'rational_witness_magnitude':amount,
            'log_approximation_error_upper':cost*error,
            'modular_equals_bare_local_energy_plus_scalar_excluded':amount>cost*error}


@lru_cache(maxsize=4)
def native_region(physical=False,keep='spinor'):
    if keep not in ('spinor','flux'):raise ValueError('select spinor or flux internal cut')
    state=native_state(physical);n=len(state['masks'])
    left=keep=='spinor';outside=13 if left else n
    sigma=partial_trace(state['rho'],n,13,'left' if left else 'right')
    log=log_polynomial(sigma)
    floor=outside/state['raw_Q_squared_normalizer']
    if floor<=0:raise ArithmeticError('the actual polynomial marginal has no positive floor')
    if left:
        u=fermion_swap().extract(state['masks'],state['masks'])
        global_u=s.kronecker_product(u,s.eye(13))
    else:
        u=s.eye(13);u[6,6]=u[12,12]=0;u[6,12]=u[12,6]=1
        global_u=s.kronecker_product(s.eye(n),u)
    u=clean(u);global_u=clean(global_u)
    if u.H*u!=s.eye(u.rows):raise ArithmeticError('complete regional intervention must be unitary')
    charges=tuple(6*(m.bit_count()-2) for m in state['masks'] for _ in range(13))
    parities=tuple((-1)**m.bit_count() for m in state['masks'] for _ in range(13))
    if any(charges[i]!=charges[j] or parities[i]!=parities[j]
           for (i,j),value in global_u.todok().items() if value):
        raise ArithmeticError('the regional apparatus violates the source gauge or parity constraint')
    after=clean(u*sigma*u.H);x=clean(after-sigma)
    full_x=clean(global_u*state['rho']*global_u.H-state['rho'])
    if partial_trace(full_x,n,13,'left' if left else 'right')!=x or s.trace(x)!=0:
        raise ArithmeticError('regional tangent does not reduce from the same source')
    hs=s.factor(s.trace(x*x));norm=fraction(row_norm(x));amplitude=F(1,8)
    # The entire convex path has the same floor, because both endpoints
    # are unitary conjugates. Reuse the resolvent bound with that floor.
    remainder=amplitude**2*fraction(hs)/(2*fraction(floor))
    generic_remainder=inherited_entropy.first_law_remainder(
        fraction(floor),fraction(hs),amplitude,norm)
    if remainder>generic_remainder:
        raise ArithmeticError('the convex-path floor weakened the inherited remainder')
    slope=s.factor(s.trace(x*log['matrix']));slope_error=2*log['operator_error_upper']
    B=log['polynomial_norm_upper']
    if B==0:raise ArithmeticError('the executed modular detector must have nonzero range')
    probabilities=[]
    for branch,weight,r in ((0,1-amplitude,sigma),(1,amplitude,after)):
        mean=s.factor(s.trace(r*log['matrix']))
        for sign in (-1,1):
            probability=s.factor(s.Rational(weight.numerator,weight.denominator)*(1+sign*mean/s.Rational(str(B)))/2)
            if probability<0:raise ArithmeticError('negative complete apparatus probability')
            probabilities.append({'preparation_coin':branch,'modular_pointer':sign,'original_probability':probability})
    if sum(r['original_probability'] for r in probabilities)!=1:
        raise ArithmeticError('the complete apparatus lost probability')
    bare=partial_trace(state['H'],n,13,'left' if left else 'right')/outside/64
    correction=modular_correction(log['matrix'],bare,log['operator_error_upper'])
    strain=geometry.rotor_strains()['isotropic'].extract(state['indices'],state['indices'])
    work=s.factor(s.trace(full_x*state['H']))
    strain_response=s.factor(s.trace(full_x*strain))
    return {'state':state,'sigma':sigma,'tangent':x,'unitary':u,'log':log,
            'internal_cut':keep,'spatial_region':False,'regional_dimension':sigma.rows,
            'complete_spinor_outcomes':16,'complete_flux_outcomes':13,
            'unsupported_spinor_outcomes':10 if physical else 0,
            'regional_minimum_weight':floor,'tangent_HS_squared':hs,
            'tangent_operator_norm_upper':norm,'amplitude':amplitude,
            'modular_slope_interval':[slope-slope_error,slope+slope_error],
            'finite_relative_entropy_upper':remainder,
            'inherited_general_tangent_relative_entropy_upper':generic_remainder,
            'entropy_change_interval':[amplitude*(slope-slope_error)-remainder,
                                       amplitude*(slope+slope_error)],
            'complete_pointer_records':probabilities,
            'apparatus_effect_sum_identity':True,'all_references_retained':True,
            'modular_correction':correction,'mean_energy_tangent':work,
            'declared_isotropic_strain_tangent':strain_response,
            'mechanical_pressure_tangent_at_declared_volume_one':-strain_response/3}


def native_summary(physical=False,keep='spinor'):
    r=native_region(physical,keep)
    omitted={'state','sigma','tangent','unitary','log'}
    return {**{k:v for k,v in r.items() if k not in omitted},
            'source_carrier_dimension':r['state']['dimension'],
            'source_Gauss_zero':physical,'beta':r['state']['beta'],
            'source_density_exactly_stationary':True,
            'varied_source_density_stationary_claimed':False,
            'joint_canonical_purification_stationary_claimed':False,
            'log_polynomial_order':r['log']['order'],
            'log_polynomial_operator_error_upper':r['log']['operator_error_upper'],
            'log_contraction_norm_upper':r['log']['contraction_norm_upper'],
            'same_H_and_clock':True,'original_source_dimension':208,
            'gauge_commuting_intervention':True,
            'entropy_is_regional_marginal_not_joint_pointer_entropy':True,
            'geometric_diamond_kernel_area_or_Einstein_premise_identified':False}


def serialize(value):
    if isinstance(value,dict):return {str(k):serialize(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serialize(v) for v in value]
    if isinstance(value,(F,s.Basic)):return str(value)
    return value


def certificate():
    sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_INPUTS}
    from model.particle_common_contract import contract
    return serialize({'schema':'oph.spatial_modular_response.v1','source_sha256':sources,
        'selected_common_state_consumer':common_contract_bound(contract.member(0)),
        'ambient_parameterized_theorem_example':spatial_faithfulness(),
        'native_boundary_pair_kinematic_control':native_boundary_pair(),
        'native_internal_controls':[native_summary(physical,keep)
            for physical in (False,True) for keep in ('spinor','flux')],
        'scope':{'spatial_restriction_applies_to_full_Wilson_product_cutoff':True,
            'actual_large_spatial_density_allocated':False,
            'physical_boundary_center_cannot_be_deleted':True,
            'finite_regional_entropy_first_law_on_actual_source':True,
            'physical_strain_identified_with_modular_diamond_energy':False,
            'fixed_proper_volume_area_variation_supplied':False,
            'Einstein_dynamics_derived':False}})
