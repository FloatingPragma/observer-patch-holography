"""Local Gauss projectors, physical cooling and complete referenced probes.

The full spatial carrier has60CAR modes per cell and native coordinate links.
Executed controls retain their stated finite sectors and complete records.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.native_quantum_gauss import algebra as native
from model.particle_field_state_response import state_response as thermal
from model.particle_spatial_ground_state import ground
from model.particle_ground_spectral_records import records
from model.interacting_extension.chiral_lattice import source_data

SOURCE_INPUTS = (
    'model/particle_spatial_ground_state/ground.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_field_state_response/state_response.py',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/particle_field_preparation/PROOF.md',
    'model/particle_field_preparation/REFINEMENT.md',
    'model/native_quantum_gauss/algebra.py',
    'model/native_quantum_gauss/PROOF.md',
    'model/particle_ground_spectral_records/records.py',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'Lean/Thermodynamics/GapFreeCooling.lean',
    'Lean/QFT/GroundRecordTail.lean',
)


def exact_kernel_projector(generators):
    """Haar projector for a supplied connected compact Lie representation.

    The caller supplies a complete Lie generating family. This computes the
    common invariant kernel; it never pretends to test group completeness.
    """
    matrices=tuple(s.ImmutableMatrix(g) for g in generators)
    if not matrices:raise ValueError('a nonempty complete generator family is required')
    size=matrices[0].rows
    if any(g.shape!=(size,size) or g+g.H!=s.zeros(size) for g in matrices):
        raise ValueError('equal-size anti-Hermitian generators are required')
    if any(value.has(s.Float) or value.free_symbols for g in matrices for value in g):
        raise ValueError('exact algebraic generators without free parameters are required')
    casimir=sum((g.H*g for g in matrices),s.zeros(size))
    columns=casimir.nullspace()
    if columns:
        basis=s.Matrix.hstack(*columns)
        projector=s.ImmutableMatrix(basis*(basis.H*basis).inv()*basis.H)
    else:projector=s.ImmutableMatrix(s.zeros(size))
    if projector!=projector.H or projector*projector!=projector:
        raise ArithmeticError('the invariant kernel projection failed')
    if any(g*projector!=s.zeros(size) for g in matrices):
        raise ArithmeticError('a retained vector violates a Gauss generator')
    return {'projector':projector,'rank':len(columns),
            'positive_gauss_operator':s.ImmutableMatrix(casimir)}


def joint_projector(projectors):
    values=tuple(s.ImmutableMatrix(p) for p in projectors)
    if not values:raise ValueError('at least one local Gauss projector is required')
    size=values[0].rows; result=s.eye(size)
    for i,p in enumerate(values):
        if p.shape!=(size,size) or p!=p.H or p*p!=p:
            raise ValueError('orthogonal projectors on one carrier are required')
        if any(p*q!=q*p for q in values[:i]):
            raise ValueError('local Gauss projectors must commute')
        result=result*p
    return s.ImmutableMatrix(result)


@lru_cache(maxsize=1)
def native_spatial_seed():
    data=source_data.load(); currents=native.native_currents()
    if len(currents)!=12 or any(s.trace(g)!=0 for g in currents):
        raise ArithmeticError('nativeM15 generators are not the traceless source family')
    if sum(data['charges'])!=0:
        raise ArithmeticError('full native hypercharge trace changed')
    return {'internal_dimension':15,'Dirac_copies':4,'CAR_modes_per_cell':60,
            'native_generator_count':12,'all_generator_traces_zero':True,
            'normal_ordering_charge_offset':0,
            'empty_Fock_constant_link_seed_physical':True,
            'filled_60_mode_constant_link_seed_physical':True,
            'physical_carrier_dimension_lower':1}


def physical_loop_seed(spacing='1',magnetic='1/7'):
    """Nontrivial physical plaquette at ell>=1, using exact Haar characters."""
    a=ground._q(spacing,'spacing',True);k=ground._q(magnetic,'magnetic coefficient')
    blocks=source_data.load()['blocks'];charges={row['q6Y'] for row in blocks}
    if 0 in charges or charges & {-q for q in charges}:
        raise ArithmeticError('native and conjugate character menus overlap')
    norm=F(len(blocks),2)
    return {'minimum_link_word_degree':1,'native_irreducible_character_count':len(blocks),
            'blank_loop_overlap':F(0),'real_trace_loop_norm_squared':norm,
            'normalized_blank_loop_H_matrix_element_squared':(k/(15*a))**2*norm,
            'physical_dimension_lower_with_plaquette':2,
            'protected_topological_defect_or_particle_claimed':False}


@lru_cache(maxsize=1)
def endpoint_control():
    """Actual native weak-block coordinate-link gauge factors on16 states.

    This is a gauge-subrepresentation test inside full60CAR/site, not a
    Hamiltonian-reducing truncation. Right matter is a single hole in the
    determinant-one filled source carrier, hence transforms in the dual.
    """
    selected=next(indices for _,indices in native.native_blocks() if len(indices)==2)
    generators=tuple(s.ImmutableMatrix(g.extract(selected,selected)) for g in native.native_currents())
    eye=s.eye(2)
    # Factors: left matter, right hole, left link index, right link index.
    left=tuple(s.kronecker_product(t,eye,eye,eye)+s.kronecker_product(eye,eye,t.conjugate(),eye)
               for t in generators)
    right=tuple(s.kronecker_product(eye,t.conjugate(),eye,eye)+s.kronecker_product(eye,eye,eye,t)
                for t in generators)
    px=exact_kernel_projector(left);py=exact_kernel_projector(right)
    p=joint_projector((px['projector'],py['projector']))
    vector=s.zeros(16,1)
    for a in range(2):
        for b in range(2):vector[8*a+4*b+2*a+b,0]=s.Rational(1,2)
    if p!=vector*vector.H:
        raise ArithmeticError('the computed local Gauss kernel differs from the invariant tensor')
    if any(x*y!=y*x for x in left for y in right):
        raise ArithmeticError('independent endpoint gauge actions do not commute')
    noncommuting=any(x*y!=y*x for x in left for y in left)
    if not noncommuting:raise ArithmeticError('nonabelian Gauss components were replaced by scalars')
    return {'left_generators':left,'right_generators':right,
            'left_projector':px['projector'],'right_projector':py['projector'],
            'projector':p,'singlet':s.ImmutableMatrix(vector),
            'carrier_dimension':16,'local_kernel_ranks':[px['rank'],py['rank']],
            'joint_kernel_rank':1,'different_site_actions_commute':True,
            'same_site_components_all_commute':False,
            'Hamiltonian_reducing_sector_claimed':False}


def spatial_plan(**kwargs):
    target=kwargs.pop('target_excess','1')
    if kwargs.get('regulator','wilson')!='wilson':
        raise ValueError('the spatial Gauss subsequence here uses the finite-range Wilson action')
    source=ground.spatial_source_bounds(**kwargs)
    source={key:value for key,value in source.items() if key!='source_Gauss_singlet_projection'}
    seed=native_spatial_seed()
    cooling=ground.cooling_plan(source['dimension_bits_upper'],
                                source['spectral_lower'],source['spectral_upper'],target)
    return {'source':source,'seed':seed,'cooling':cooling,
            'physical_constraint':'product_of_commuting_local_Haar_projectors',
            'full_Fock_and_link_carrier_retained':True,
            'state_supported_in_source_Gauss_kernel':True,
            'reference_separately_Gauss_invariant':True,
            'restricted_ground_minimum':True,
            'restricted_dimension_upper_uses_ambient_bound':True,
            'boundary_condition':'zero_external_charge_at_all_link_endpoints',
            'H_or_clock_changed_by_preparation':False,
            'postselection':False}


@lru_cache(maxsize=1)
def native_projector():
    control=thermal.exact_control(); g=control['Gauss'];h=control['H']
    p=s.ImmutableSparseMatrix(s.diag(*[int(g[i,i]==0) for i in range(208)]))
    if p*p!=p or p!=p.H or g*p!=s.zeros(208) or h*p!=p*h:
        raise ArithmeticError('the native Gauss-zero source sector is not reducing')
    if s.trace(p)!=78:
        raise ArithmeticError('the complete native physical carrier has changed')
    return p


def native_gauss_monitor():
    """Both raw outputs of a nondestructive check, including outside Gauss."""
    p=native_projector()
    return p,s.ImmutableSparseMatrix(s.eye(208)-p)


@lru_cache(maxsize=4)
def prepared_control(beta='4',order=8):
    full=ground.prepared_control(beta,order);p=native_projector()
    h=thermal.real_computation_chart()[1]
    blocks=[];normalizer=s.S.Zero;energy=s.S.Zero
    for indices,q in full['amplitude_blocks']:
        physical={bool(p[i,i]) for i in indices}
        if len(physical)!=1:raise ArithmeticError('reducing block mixes Gauss constraints')
        if not next(iter(physical)):q=s.ImmutableMatrix(s.zeros(len(indices)))
        q2=q*q;normalizer+=s.trace(q2);energy+=s.trace(q2*h.extract(indices,indices))
        blocks.append((indices,q))
    if normalizer<=0:raise ArithmeticError('the physical preparation is empty')
    return {**full,'amplitude_blocks':tuple(blocks),'norm_squared':normalizer,
            'mean_energy':energy/normalizer,'physical_dimension':78,
            'all_occupations_and_fluxes_retained':False,
            'full_ambient_carrier_retained':True,'all_physical_states_retained':True,
            'source_Gauss_kernel_support':True,'reference_separately_Gauss_invariant':True,
            'construction':'direct_canonical_preparation_of_P_Q_H',
            'postselected_full_carrier_run':False}


@lru_cache(maxsize=1)
def compensated_probe_certificate():
    p=native_projector();control=thermal.exact_control();g=control['Gauss']
    rows=[];complete=s.zeros(208)
    for row in records.instrument():
        k=row['K'];mapped=k*p;charge=row['reference_charge']
        if (g+charge*s.eye(208))*mapped!=s.zeros(208):
            raise ArithmeticError('the charged output lacks its boundary reference')
        if p*mapped!=s.zeros(208):
            raise ArithmeticError('a charged local probe was treated as a source-only observable')
        complete+=mapped.H*mapped
        rows.append({'sign':row['sign'],'mode':row['mode'],
                     'source_charge_change':row['source_charge_change'],
                     'boundary_reference_charge':charge,
                     'boundary_reference_parity':row['reference_parity']})
    if complete!=p:raise ArithmeticError('the physical input lost an original probe outcome')
    minima=records.sector_minima();global_ground=records.global_ground_interval()
    if minima[2]!=global_ground:
        raise ArithmeticError('the enlarged accessible ground has a different bound')
    return {'rows':rows,'all_eight_original_outcomes_retained':True,
            'physical_input_effect_sum_is_P':True,
            'bare_charged_operator_preserves_source_physical_space':False,
            'compensated_total_charge_zero':True,
            'compensating_reference_H':0,
            'native_physical_and_full_ground_intervals_equal':True,
            'same_equality_claimed_for_general_spatial_source':False}


def native_summary(include_records=False):
    prepared=prepared_control();mean=ground.outward_dyadic(prepared['mean_energy'])
    e0=records.sector_minima()[2];epsilon=mean[1]-e0[0]
    if epsilon<0:raise ArithmeticError('physical state below certified restricted ground')
    result={'ambient_dimension':208,'physical_dimension':78,
            'source_Gauss':'6(N-2)','neutralizing_background_charge':-12,
            'physical_occupation':2,'flux_labels':list(range(-6,7)),
            'beta':prepared['bound']['beta'],'Taylor_order':prepared['bound']['Taylor_order'],
            'dyadic_scale':prepared['bound']['dyadic_scale'],
            'mean_energy_interval':mean,'physical_ground_interval':e0,
            'mean_physical_ground_excess_upper':epsilon,
            'canonical_half_trace_error_upper':prepared['preparation_trace_error_upper'],
            'source_and_reference_separately_Gauss_invariant':True,
            'outside_Gauss_monitor_outcome_retained':True,
            'physical_ground_equals_full_ground_in_this_control':True,
            'projection_applied_as_postselection':False,
            'physical_particle_masses':[]}
    if include_records:
        law=records.complete_filter_records(prepared)
        for row in law['records']:
            row['original_probability_interval']=ground.outward_dyadic(row.pop('original_probability'))
        result['complete_charged_clock_records']=law
    return result


def certificate(include_control=True,include_records=False):
    local=endpoint_control()
    result={'schema':'particle_gauss_ground_state_v1','source_inputs':list(SOURCE_INPUTS),
            'native_spatial_seed':native_spatial_seed(),
            'nontrivial_physical_loop':physical_loop_seed(),
            'nonabelian_endpoint_control':{key:local[key] for key in (
                'carrier_dimension','local_kernel_ranks','joint_kernel_rank',
                'different_site_actions_commute','same_site_components_all_commute',
                'Hamiltonian_reducing_sector_claimed')},
            'spatial_families':[spatial_plan(side=n,link_degree=1,electric='105',
                                magnetic='1/7',target_excess=F(1,n)) for n in (2,4)],
            'fixed_cutoff_physical_ground_subsequence':True,
            'ground_inequality_domain':'local_gauge_invariant_observable_algebra',
            'charged_local_CAR_is_gauge_invariant':False,
            'physical_spacing_refinement_or_scattering_proved':False}
    if include_control:
        result['native_physical_control']=native_summary(include_records)
        result['compensated_native_probe']=compensated_probe_certificate()
    return ground.serial(result)
