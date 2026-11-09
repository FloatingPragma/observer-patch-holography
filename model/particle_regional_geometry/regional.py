"""Exact restricted states and modular obstructions on the actual Higgs source.

The regions are internal scalar/family subalgebras, not geometric balls.
The same faithful state, Hamiltonian and preparation are inherited from
thermal_clock; every occupation sector and reference is accounted for.
"""
from functools import lru_cache
import json

import sympy as s

from model.particle_observer_outputs import higgs_records as source
from model.particle_observer_outputs import thermal_clock as thermal
from model.particle_candidate_atlas.interacting_chiral import quantum_higgs_operators


SOURCE_INPUTS = (
    'model/particle_observer_outputs/higgs_records.py',
    'model/particle_observer_outputs/thermal_clock.py',
    'model/particle_observer_outputs/THERMAL_CLOCK.md',
    'model/particle_observer_outputs/readout_contract.py',
    'model/particle_candidate_atlas/interacting_chiral.py',
    'model/field_regions/README.md',
    'model/native_gravity_identification/PROOF.md',
    'model/matched_einstein_limit/PROOF.md',
    'Lean/QFT/CPRestrictionNet.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
)


def clean(matrix):
    return source._matrix(matrix)


@lru_cache(maxsize=3)
def family_split(family=0):
    if type(family) is not int or not 0 <= family < 3:
        raise ValueError('one of the three declared source families is required')
    selected = tuple(range(6*family,6*family+6))
    exterior = tuple(j for j in range(18) if j not in selected)
    rows = []
    for scalar, mask in thermal.thermal_fixture()['source']['basis']:
        local = sum(((mask>>bit)&1)<<j for j,bit in enumerate(selected))
        outside = sum(((mask>>bit)&1)<<j for j,bit in enumerate(exterior))
        inversions = sum((mask>>a)&(mask>>b)&1 for a in selected for b in exterior if b<a)
        rows.append((local,(scalar,outside),(-1)**inversions))
    if len(set((local,env) for local,env,_ in rows)) != 307:
        raise ArithmeticError('the signed CAR split is not an isometry')
    return tuple(rows)


def restrict_family(density, family=0, reference_dimension=1):
    """Trace the scalar and exterior CAR modes, retaining an arbitrary reference."""
    if type(reference_dimension) is not int or reference_dimension < 1:
        raise ValueError('a positive retained-reference dimension is required')
    density = s.ImmutableSparseMatrix(density)
    r = reference_dimension
    if density.shape != (307*r,307*r):
        raise ValueError('the complete source/reference density dimension is required')
    labels = family_split(family)
    entries = {}
    for (i,j),value in density.todok().items():
        si,ri = divmod(i,r); sj,rj = divmod(j,r)
        ai,ei,signi = labels[si]; aj,ej,signj = labels[sj]
        if ei == ej:
            key = ai*r+ri,aj*r+rj
            entries[key] = entries.get(key,0)+signi*signj*value
    return clean(s.ImmutableSparseMatrix(64*r,64*r,entries))


def inject_family(operator,family=0):
    """Dual compression; a homomorphism on the retained number-preserving algebra."""
    operator = s.ImmutableSparseMatrix(operator)
    if operator.shape != (64,64):
        raise ValueError('a full six-mode CAR operator is required')
    labels = family_split(family)
    entries = {}
    for i,(a,e,signi) in enumerate(labels):
        for j,(b,f,signj) in enumerate(labels):
            if e == f and operator[a,b] != 0:
                entries[i,j] = signi*signj*operator[a,b]
    return clean(s.ImmutableSparseMatrix(307,307,entries))


def trace_expectation(operator,family=0):
    """Trace-preserving conditional expectation onto the retained local algebra.

    Exterior multiplicities are133,24,2 in local number sectors0,1,2.
    This tracial expectation does not preserve the interacting Gibbs state.
    """
    local=restrict_family(operator,family)
    multiplicities=(133,24,2)
    entries={(a,b):value/multiplicities[a.bit_count()]
             for (a,b),value in local.todok().items()}
    return inject_family(s.ImmutableSparseMatrix(64,64,entries),family)


@lru_cache(maxsize=3)
def family_region(family=0):
    state = thermal.thermal_fixture()
    data = state['source']
    rho = restrict_family(state['rho'],family)
    y = s.Integer(family+1)
    alpha,beta = state['occupation_weights']['empty'],state['occupation_weights']['bright']
    w,g2 = state['per_state_weights']['dark'],data['g_squared']
    p0 = s.simplify(alpha+132*w+(beta-w)*(1-4*y*y/g2))
    p1,p2 = 24*w,2*w
    pb = s.simplify(p2+(beta-w)*2*y*y/g2)
    projectors = [s.ImmutableSparseMatrix(64,64,{(mask,mask):1 for mask in range(64)
                                                if mask.bit_count()==n}) for n in range(3)]
    pair = s.zeros(64)
    vectors = []
    for scalar in (1,2):
        entries = {}
        for j,(local,(h,out),sign) in enumerate(family_split(family)):
            if h==scalar and out==0 and data['bright_vector'][j,0]!=0:
                entries[local,0] = sign*data['bright_vector'][j,0]
        vector = s.ImmutableSparseMatrix(64,1,entries)
        if (vector.H*vector)[0] != 2*y*y:
            raise ArithmeticError('the literal family pair vector has changed')
        vectors.append(vector/s.sqrt(2*y*y))
        pair += vector*vector.H/(2*y*y)
    pair = clean(pair)
    P0,P1,P2 = projectors
    expected = clean(p0*P0+p1*P1+p2*P2+(pb-p2)*pair)
    if rho!=expected or s.trace(rho)!=1:
        raise ArithmeticError('the actual reduced source differs from the regional spectral formula')
    if any(p.is_positive is not True for p in (p0,p1,p2,pb)) or pb<=p2:
        raise ArithmeticError('regional positive support or pair splitting changed')
    dark = clean(P2-pair)
    spectra = ((p0,P0),(p1,P1),(pb,pair),(p2,dark))
    modular = clean(sum((-s.log(p)*P for p,P in spectra),s.zeros(64)))
    support = clean(P0+P1+P2)
    if s.trace(support)!=22 or s.trace(pair)!=2 or s.trace(dark)!=13:
        raise ArithmeticError('a regional occupation/degeneracy sector was omitted')
    # Left-spin Z commutes with the gauge action. It maps each weak
    # component of the spin-singlet pair to its orthogonal triplet partner.
    # Summing both weak components gives a gauge-invariant local probe.
    spin_z=s.diag(*[((mask>>0)&1)-((mask>>1)&1)+((mask>>2)&1)-((mask>>3)&1)
                   for mask in range(64)])
    dark_vectors=[clean(spin_z*b) for b in vectors]
    for b,dark_vector in zip(vectors,dark_vectors,strict=True):
        if (dark_vector.H*dark_vector)[0]!=1 or (b.H*dark_vector)[0]!=0 or clean(dark*dark_vector-dark_vector)!=s.zeros(64,1):
            raise ArithmeticError('the retained gauge-covariant dark probe differs')
    probe=clean(sum((b*d.H+d*b.H for b,d in zip(vectors,dark_vectors,strict=True)),s.zeros(64)))
    coefficient = s.log(p2/pb)
    oriented=clean(sum((b*d.H-d*b.H for b,d in zip(vectors,dark_vectors,strict=True)),s.zeros(64)))
    if clean(modular*probe-probe*modular-coefficient*oriented)!=s.zeros(64):
        raise ArithmeticError('the regional pair modular frequency differs')
    number = clean(P1+2*P2)
    global_number = inject_family(number,family)
    comm_vac = (data['H']*global_number-global_number*data['H'])*data['vacuum']
    norm_squared = s.simplify((comm_vac.H*comm_vac)[0])
    if norm_squared != 16*y*y:
        raise ArithmeticError('the actual source family-number dynamics changed')
    # This coherence connects central local-number sectors. Every bimodule
    # conditional expectation onto the local number algebra kills it.
    vector = clean(global_number*data['bright_vector']/(4*y))
    if (vector.H*vector)[0]!=1:
        raise ArithmeticError('the family bright probe is not normalized')
    cross = clean(data['vacuum']*vector.H+vector*data['vacuum'].H)
    cross_mean = s.simplify(s.trace(state['rho']*cross))
    if cross_mean != -s.Rational(131068,1077135)*y:
        raise ArithmeticError('the global faithful coherence witness changed')
    return {'family':family,'rho':rho,'support':support,'number':number,
            'occupation_projectors':projectors,'pair_projector':pair,'pair_vectors':vectors,
            'dark_projector':dark,'dark_vector':dark_vectors[0],'probe':probe,
            'modular_H':modular,'spectra':spectra,'probabilities_per_state':(p0,p1,pb,p2),
            'ranks':(1,6,2,13),'pair_modular_coefficient':coefficient,
            'global_number_commutator_vacuum_norm_squared':norm_squared,
            'cross_center_coherence':cross,'cross_center_mean':cross_mean}


def scalar_region():
    state = thermal.thermal_fixture()
    basis = state['source']['basis']
    entries = {}
    for (i,j),value in state['rho'].todok().items():
        hi,mi = basis[i]; hj,mj = basis[j]
        if mi==mj: entries[hi,hj]=entries.get((hi,hj),0)+value
    reduced = clean(s.ImmutableSparseMatrix(3,3,entries))
    empty = state['occupation_weights']['empty']
    each = s.simplify((1-empty)/2)
    if reduced!=s.diag(empty,each,each) or empty>=each:
        raise ArithmeticError('the scalar regional source law changed')
    return {'reduced_density':reduced,'per_state_probabilities':(empty,each,each),
            'algebra_preserving_the_retained_source_sector':'C direct_sum M2',
            'density_is_central_on_that_algebra':True,
            'regional_modular_flow_on_that_algebra':'identity',
            'formal_number_modular_coefficient':s.log(empty/each),
            'positive_multiple_of_bare_scalar_number_energy':False,
            'global_H_NH_commutator_vacuum_norm_squared':state['source']['g_squared']}


def entropy_response(family=0):
    """A source-admitted center-preserving pair/dark entropy variation."""
    r=family_region(family)
    state=thermal.thermal_fixture()
    delta=clean(r['pair_projector']/2-r['dark_projector']/13)
    lifted=clean(inject_family(delta,family)/2)
    if restrict_family(lifted,family)!=delta or s.trace(lifted)!=0:
        raise ArithmeticError('the same-source state variation has a different regional tangent')
    if s.trace(delta*delta)!=s.Rational(15,26):
        raise ArithmeticError('the complete regional tangent lost pair or dark degeneracy')
    slope=s.simplify(s.expand_log(s.trace(delta*r['modular_H']),force=True))
    if s.simplify(slope-r['pair_modular_coefficient'])!=0:
        raise ArithmeticError('the regional entropy slope differs from its actual pair coefficient')
    epsilon=state['per_state_weights']['excited']
    minimum=r['probabilities_per_state'][3]
    remainder=s.simplify(epsilon**2*s.Rational(15,26)/(2*(minimum-epsilon/2)))
    return {'regional_tangent':delta,'global_tangent':lifted,
            'source_global_minimum_density':epsilon,
            'declared_variation_size':epsilon,
            'perturbed_global_density_lower_bound':3*epsilon/4,
            'global_tangent_operator_norm':s.Rational(1,4),
            'regional_tangent_Hilbert_Schmidt_norm_squared':s.Rational(15,26),
            'regional_tangent_operator_norm':s.Rational(1,2),
            'regional_minimum_density':minimum,
            'entropy_first_derivative':slope,
            'finite_relative_entropy_remainder_upper':remainder,
            'local_number_center_probabilities_unchanged':True,
            'same_H_and_clock_reproducible_perturbed_preparation':True}


@lru_cache(maxsize=1)
def gauge_probe_certificate():
    region=family_region(0)
    operators=(region['rho'],region['pair_projector'],region['probe'])
    generators=quantum_higgs_operators()['generators']
    for full in generators:
        current=full[:64,:64]
        for operator in operators:
            if clean(current*operator-operator*current)!=s.zeros(64):
                raise ArithmeticError('the actual regional modular observable changes a source gauge current')
    return {'source_gauge_generators':len(generators),'operator_commutators_checked':len(generators)*len(operators),
            'regional_pair_probe_is_gauge_invariant':True,
            'probe_support':'both weak components of left-right pair; singlet/triplet spin interference',
            'spatial_rotation_invariance_of_apparatus_claimed':False}


def certificate():
    rows=[]
    for f in range(3):
        r=family_region(f)
        entropy=entropy_response(f)
        center=[s.simplify(s.trace(r['rho']*P)) for P in r['occupation_projectors']]
        rows.append({'family':f,'full_local_CAR_dimension':64,'positive_support_dimension':22,
            'outside_support_zero_outcomes':42,'ranks':r['ranks'],
            'eigenvalues_per_state':r['probabilities_per_state'],
            'number_center_probabilities':center,
            'pair_modular_coefficient':r['pair_modular_coefficient'],
            'global_number_commutator_vacuum_norm_squared':r['global_number_commutator_vacuum_norm_squared'],
            'cross_center_coherence_mean':r['cross_center_mean'],
            'entropy_pair_dark_response':{key:value for key,value in entropy.items()
                                         if not isinstance(value,s.MatrixBase)},
            'bare_onsite_or_number_polynomial_modular_identification':False,
            'global_state_preserving_bimodule_conditional_expectation_onto_region':False})
    return serial({'schema':'oph.particle_regional_geometry.v1',
        'same_source':'Y=diag(1,2,3), d=1, beta=log(2), full307 faithful Gibbs state',
        'same_preparation':'canonical conjugate-reference purification from thermal_clock',
        'conjugate_reference_retained_in_full_program':True,
        'scalar_region':scalar_region(),'family_regions':rows,
        'actual_regional_gauge_probe':gauge_probe_certificate(),
        'same_regional_CP_map_applies_with_arbitrary_references':True,
        'regional_homomorphism_domain':'number-preserving algebra C direct_sum M6 direct_sum M15 on positive support',
        'compression_is_homomorphism_on_all64_CAR_operators':False,
        'partial_trace_conditioning_or_weight_renormalization':False,
        'computed_positive_support_modular_H':True,
        'regional_modular_flow':'sum_ab (p_a/p_b)^(i*s) P_a A P_b on the22-dimensional positive support',
        'regional_tracial_conditional_expectation_constructed':True,
        'that_expectation_preserves_interacting_Gibbs_state':False,
        'spatial_ball_or_geometric_area_identified':False,
        'regional_diamond_kernel_premise_discharged':False,
        'bare_kernel_obstruction_scope':'specified internal regional algebras and the declared Gibbs source only',
        'stationary_repair_entropy_used_as_geometric_entropy':False,
        'source_inputs':list(SOURCE_INPUTS)})


def serial(value):
    if isinstance(value,s.MatrixBase):return [[str(x) for x in row] for row in value.tolist()]
    if isinstance(value,s.Basic):return str(value)
    if isinstance(value,dict):return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serial(v) for v in value]
    return value


if __name__=='__main__':
    print(json.dumps(certificate(),indent=2))
