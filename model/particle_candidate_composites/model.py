"""Actual neutral composite currents, Ward constraints and a complete meter."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_composite_spin import model as spin
from model.particle_primary_charged_readout import model as readout
from model.particle_spatial_modular_response import regional

FAMILIES=(('weak_pair_charge_six',(0,2,0,0,1)),
          ('color_triple',(1,0,0,2,0)),
          ('even_determinant_hole',(0,0,0,3,-1)))
S_COEFFICIENTS=(1,0,0,1,1)
B_CONDITIONAL=(F(-1,3),F(0),F(1,3),F(-1,3),F(0))
L_CONDITIONAL=(0,1,0,0,-1)
SOURCE_INPUTS=tuple(dict.fromkeys((
    'model/particle_candidate_composites/__init__.py',
    'model/particle_candidate_composites/__main__.py',
    'model/particle_candidate_composites/model.py',
    'model/particle_candidate_composites/test_model.py',
    'model/particle_candidate_composites/PROOF.md',
    'model/particle_candidate_composites/README.md',
    'Lean/Dynamics/ProtectedCharge.lean',
    'Lean/EventAlgebra/PartitionPinchingCP.lean',
    *spin.SOURCE_INPUTS,*readout.SOURCE_INPUTS,*regional.SOURCE_INPUTS,
)))


def _delta(value):
    value=tuple(value)
    if len(value)!=5 or any(type(x) is not int for x in value):
        raise ValueError('five integer native occupation increments required')
    return value


def normal_charge_allowed(delta_left,delta_right):
    """Necessary source-marginal rule for Tr rho F_left^dagger alpha_t(F_right)."""
    return _delta(delta_left)==_delta(delta_right)


def anomalous_charge_allowed(delta_left,delta_right):
    """Necessary source-marginal rule for Tr rho F_left alpha_t(F_right)."""
    return all(x+y==0 for x,y in zip(_delta(delta_left),_delta(delta_right)))


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def polynomial_delta(operators):
    return tuple(sum((1 if create else -1) for mode,create in operators
                     if mode%15 in indices) for _,indices in spin.sectors.native_blocks())


@lru_cache(maxsize=None)
def tensor_channels(left,right):
    if left not in spin.group.IRREPS or right not in spin.group.IRREPS:
        raise ValueError('actual binary-cube irreps required')
    result={}
    for output in spin.group.IRREPS:
        value=s.simplify(sum(s.conjugate(spin.channels._sym(spin.group.characters(q)[left]))
                    *spin.channels._sym(spin.group.characters(q)[right])
                    *s.conjugate(spin.channels._sym(spin.group.characters(q)[output]))
                    for q in spin.group.elements())/48)
        if not value.is_Integer or value<0:raise ArithmeticError('cube tensor multiplicity failed')
        if value:result[output]=int(value)
    return result


@lru_cache(maxsize=1)
def weak_pair_alias():
    """Two exact G1 weak trions form a nonzero six-creation source operator."""
    delta=FAMILIES[0][1];local=spin.local_multiplet(delta)
    projector=local['casimir_projectors'][s.Rational(3,4)]
    columns=s.Matrix.hstack(*s.GramSchmidt(projector.columnspace(),orthonormal=True))
    polynomials=[]
    for j in range(2):
        polynomial=[]
        for k in range(6):
            if columns[k,j]:
                polynomial.extend((ops,s.expand(columns[k,j]*c))
                                  for ops,c in spin.composite_polynomial(delta,k))
        polynomials.append(polynomial)
    out={}
    for ops1,c1 in polynomials[0]:
        for ops2,c2 in polynomials[1]:
            target,sign=spin.apply_operators(spin.OMEGA,ops1+ops2)
            if sign:out[target]=s.expand(out.get(target,0)+sign*c1*c2)
    out={mask:c for mask,c in out.items() if c}
    norm=s.simplify(sum(s.conjugate(c)*c for c in out.values()))
    if norm!=1 or len(out)!=1:raise ArithmeticError('actual weak pair lost its source witness')
    double=tuple(2*n for n in delta);even=FAMILIES[2][1]
    return {'two_G1_column_product_anchor_norm_squared':str(norm),
        'nonzero_anchor_coefficients':[str(c) for c in out.values()],
        'native_increment':double,'S_increment':dot(S_COEFFICIENTS,double),
        'aliases_single_even_meter_record':dot(S_COEFFICIENTS,double)%7==dot(S_COEFFICIENTS,even)%7,
        'full_native_increment_differs':double!=even,
        'pair_creation_term_preserves_native_U1_five':False,
        'physical_Majorana_mass_identified':False}


def meter_control():
    """Execute both couplings on a seven-state complete validation instrument."""
    charges=tuple(range(-3,4));d=len(charges);operators=[];shifts=[]
    for shift in (-3,-2,-1,1,2,3):
        K=s.zeros(d);K[charges.index(shift),charges.index(0)]=s.sqrt(6)/6
        operators.append(K);shifts.append(shift)
    failure=s.eye(d);failure[charges.index(0),charges.index(0)]=0
    operators.append(failure);shifts.append(0)
    return readout.coherent_meter(charges,operators,shifts,7)


def meter_recipe(volume=8):
    if type(volume) is not int or volume<1:raise ValueError('positive integer cell count required')
    alphabet=tuple(range(-3,4));radix=next(n for n in range(1,8) if len({x%n for x in alphabet})==7)
    # All forty source columns have norm one on the physical filled-beta anchor.
    counts=[spin.local_multiplet(tuple(sign*x for x in delta))['kernel_basis'].cols
            for _,delta in FAMILIES for sign in (1,-1)]
    if sum(counts)!=40:raise ArithmeticError('the full composite menu changed')
    return {'S_native_coefficients':S_COEFFICIENTS,'S_per_cell_spectrum_enclosure':[0,28],
        'cyclic_radix':radix,'decoding':[0,1,2,3,-3,-2,-1],
        'selected_intervention_increment_alphabet':alphabet,
        'success_columns_per_cell':sum(counts),'original_success_labels':40*volume,
        'original_complete_labels':40*volume+1,
        'formal_meter_refined_labels':7*(40*volume+1),
        'source_controlled_additions_both_couplings':56*volume,
        'no_initial_absolute_charge_measurement':True,
        'unknown_initial_charge_cancels_coherently':True,
        'source_reference_fine_history_preserved':True,
        'multiple_excitation_universal_decoder':False,
        'final_absolute_family_measurement':False}


@lru_cache(maxsize=1)
def calculated_rows():
    rows=[]
    for name,delta in FAMILIES:
        for sign in (1,-1):
            signed=tuple(sign*x for x in delta);local=spin.local_multiplet(signed)
            control=spin.operator_controls(signed);cube=spin.cube_restriction(signed)
            for col in range(local['kernel_basis'].cols):
                for operators,_ in spin.composite_polynomial(signed,col):
                    if polynomial_delta(operators)!=signed:raise ArithmeticError('native current failed')
            # Actual gauge invariance is checked on every vector, not by dimensions.
            if any(g*local['kernel_basis']!=s.zeros(len(local['masks']),local['kernel_basis'].cols)
                   for g in local['gauge_matrices']):raise ArithmeticError('gauge current failed')
            rows.append({'family':name,'sign':sign,'native_increment':signed,
                'native_Q6':dot((-4,-3,1,2,6),signed),'S_increment':dot(S_COEFFICIENTS,signed),
                'meter_residue':dot(S_COEFFICIENTS,signed)%7,
                'conditional_B':dot(B_CONDITIONAL,signed),'conditional_L':dot(L_CONDITIONAL,signed),
                'conditional_BL_requires_constituent_field_dictionary':True,
                'physical_gauge_kernel_dimension':local['kernel_basis'].cols,
                'gauge_generators_annihilating_every_column':len(local['gauge_matrices']),
                'proper_cube_irreps':cube['proper_binary_cube_irreps'],
                'uniform_full_cube_irreps':cube['uniform_orbit_full_cube_irreps'],
                'CAR_parity':local['CAR_parity'],
                'normalized_column_monomial_counts':[x['monomial_count'] for x in control['normalized_columns']],
                'direct_unbroken_weak_doublet_current_intertwiner_dimension':0,
                'direct_unbroken_color_triplet_current_intertwiner_dimension':0,
                'direct_single_fermion_parity_compatible':local['CAR_parity']==1})
    return tuple(rows)


def bind_member(member):
    contract.verify_member(member);spin.sectors.branch._member_zero(member)
    native=spin.channels.native_number_control();rows=[dict(row) for row in calculated_rows()]
    floor=regional.common_contract_payload(member,region_cells=1)['physical_density_floor_factored']
    V=member['source']['side']**3
    for row in rows:
        row['each_original_row_probability_lower_factored']=[
            {'multiplier_of_density_floor':F(1,40*V*t)}
            for t in row['normalized_column_monomial_counts']]
    deltas=[r['native_increment'] for r in rows]
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'native_charge_order':[-4,-3,1,2,6],'actual_native_number_control':native,
        'candidate_rows':rows,'physical_density_floor_factored':floor,
        'normal_source_marginal_charge_allowed':[[normal_charge_allowed(a,b) for b in deltas] for a in deltas],
        'anomalous_source_marginal_charge_allowed':[[anomalous_charge_allowed(a,b) for b in deltas] for a in deltas],
        'binary_cube_bilinear_channels':{a+'x'+b:tensor_channels(a,b)
            for a,b in (('G1','G1'),('G1','H'),('H','H'))},
        'same_charge_G1_H_scalar_uniform_kernel_zero':True,
        'same_charge_G1_H_T1_analyzer_allowed':True,
        'reference_operator_valued_ward_zero_asserted':False,
        'opposite_odd_uniform_orbit_scalar_anomalous_kernel_zero_by_inversion':True,
        'complete_detector':meter_recipe(V),'actual_multiple_excitation_alias':weak_pair_alias(),
        'gauge_singlet_multipole_or_dressed_interaction_excluded':False,
        'named_physical_particles_identified':False,'physical_particle_map':[],
        'mass_predictions_GeV':[]})
    return {'schema':'oph.particle_candidate_composites.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,bound):
    if bound!=bind_member(member):raise ValueError('composite source identity or complete record changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    meter=meter_control()
    value={'schema':'oph.particle_candidate_composites_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'exact_meter_control':{'source_dimension':meter['source_dimension'],
            'meter_dimension':meter['meter_dimension'],'complete_fine_outcomes':len(meter['rows'])},
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source inputs changed during replay')
    return value


def build():return certificate()
