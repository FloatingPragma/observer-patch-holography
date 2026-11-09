"""Coherent finite-aperture instruments on the selected physical source."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from model.particle_common_contract import contract
from model.particle_filtered_probe import filter as apparatus
from model.particle_primary_field_algebra import algebra
from model.particle_primary_spatial_algebra import spatial
from model.particle_physical_boundary_state import boundary

SOURCE_INPUTS = (
    'model/particle_primary_momentum/momentum.py',
    'model/particle_primary_momentum/PROOF.md',
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_field_algebra/algebra.py',
    'model/particle_primary_spatial_algebra/spatial.py',
    'model/particle_primary_spatial_algebra/PROOF.md',
    'model/particle_filtered_probe/filter.py',
    'model/particle_filtered_probe/PROOF.md',
    'model/particle_physical_boundary_state/boundary.py',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_joint_translation_readout/PROOF.md',
    'Lean/QFT/PeriodicMomentum.lean',
    'Lean/Variational/TranslationInvariantComposedInstance.lean',
)


def fourier_instrument(probes):
    """Execute a small exact coherent address/Fourier instrument, retaining failures."""
    probes = tuple(apparatus.exact_matrix(b) for b in probes)
    if not probes or len(probes)>16 or probes[0].rows>16:
        raise ValueError('exact control requires 1..16 probes of dimension at most16')
    m=len(probes); d=probes[0].rows
    if any(b.shape!=(d,d) for b in probes):
        raise ValueError('all probes must act on the same source')
    columns=[]
    for b in probes:
        u=apparatus.hermitian_dilation(b)
        columns.append((u[:d,:d],u[d:,:d]))
    phase=s.expand_complex(s.exp(-2*s.pi*s.I/m))
    rows=[]
    for branch in (0,1):
        for k in range(m):
            K=apparatus.clean(sum((phase**(k*j)*columns[j][branch]
                                 for j in range(m)),s.zeros(d))/m)
            rows.append({'branch':branch,'wavevector_index':k,'K':K})
    total=apparatus.clean(sum((r['K'].H*r['K'] for r in rows),s.zeros(d)))
    if total != s.eye(d):
        raise ArithmeticError('the spatial Fourier instrument lost a fine record')
    return {'rows':rows,'source_dimension':d,'addresses':m,
            'complete_effect':total,'outcomes':2*m,
            'source_translation_operator_used':False}


def aperture_bound(addresses=4, shifted_steps=1, spacing='1/4'):
    m=apparatus.integer(addresses,'address count',2)
    r=apparatus.integer(shifted_steps,'shifted address steps',0)
    h=apparatus.exact(spacing,'address spacing',True)
    return {'addresses':m,'shifted_steps':r,'address_spacing':h,
            'aperture_period_length':m*h,
            'operator_phase_covariance_error_upper':min(F(2),F(2*r,m)),
            'one_step_error_divided_by_spacing_upper':F(2,m)/h,
            'finite_aperture_has_exact_transfer_momentum':False}


def open_box_geometry(side=8):
    L=apparatus.integer(side,'open side',3)
    before=((L-2,0,0),(L-1,0,0))
    after=((L-1,0,0),(0,0,0))
    native_delta=tuple(y-x for x,y in zip(*after))
    # All source matter bonds have a positive unit coordinate displacement.
    allowed=native_delta in ((1,0,0),(0,1,0),(0,0,1))
    return {'side':L,'native_bond_before':before,'cyclic_image_bond':after,
            'cyclic_image_is_native_bond':allowed,
            'incoming_links_at_origin':0,
            'incoming_links_at_shifted_origin':1,
            'owned_outgoing_links_per_cell':3,
            'ghost_links_replaced_by_periodic_links':False,
            'cyclic_owner_relabelling_preserves_source_incidence':False}


def blank_loop_geometry(anchors, axes=(0,1)):
    """Check the unique-edge Haar orthogonality premise for actual plaquettes."""
    anchors=tuple(tuple(x) for x in anchors)
    if (not anchors or len(set(anchors))!=len(anchors)
            or any(len(x)!=3 or any(type(v) is not int for v in x) for x in anchors)
            or len(axes)!=2 or not 0<=axes[0]<axes[1]<3):
        raise ValueError('distinct integral plaquette anchors and an increasing plane required')
    i,j=axes
    def shift(x,axis):return tuple(v+int(k==axis) for k,v in enumerate(x))
    edges=[frozenset(((x,i),(shift(x,i),j),(shift(x,j),i),(x,j))) for x in anchors]
    witnesses=[]
    for k in range(len(edges)):
        for l in range(k):
            unique=sorted(edges[k]-edges[l])
            if not unique:raise ArithmeticError('distinct plaquettes lack an independent Haar edge')
            witnesses.append({'pair':(l,k),'unique_edge':unique[0]})
    return {'anchors':anchors,'edge_sets':[sorted(row) for row in edges],
            'pairwise_unique_edge_witnesses':witnesses,
            'normalized_blank_Gram_diagonal':F(1,90),
            'normalized_blank_Gram_offdiagonal':F(0),
            'uses_inherited_native_Haar_character_orthogonality':True}


@lru_cache(maxsize=1)
def exact_control():
    X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]])
    Z=s.diag(1,-1); rho=(s.eye(2)+Y/2)/2
    probes=(Z/2,X/2,-Z/2,-X/2)
    instrument=fourier_instrument(probes)
    rows=[]; d=2
    root=apparatus.positive_sqrt(rho)
    omega=s.Matrix([root[i,j] for j in range(d) for i in range(d)])
    joint=s.zeros(d*d)
    gram=s.Matrix([[s.trace(rho*b.H*c) for c in probes] for b in probes])
    for row in instrument['rows']:
        K=row['K']; vector=s.kronecker_product(s.eye(d),K)*omega
        joint+=vector*vector.H
        p=s.simplify((vector.H*vector)[0])
        assert p==s.simplify(s.trace(rho*K.H*K))
        rows.append({**row,'p':p})
    # Column vectorization puts the unchanged reference in the first factor.
    reference=s.Matrix([[sum(joint[j*d+i,k*d+i] for i in range(d))
                         for k in range(d)] for j in range(d)])
    if apparatus.clean(reference-rho.T)!=s.zeros(d):
        raise ArithmeticError('complete Fourier outputs changed the untouched reference')
    fifth=probes[0]
    shifted=fourier_instrument((*probes[1:],fifth))
    for k in range(4):
        phase=(-s.I)**k
        edge=s.conjugate(phase)*(fifth-probes[0])/4
        if apparatus.clean(shifted['rows'][k]['K']-s.conjugate(phase)*rows[k]['K']-edge)!=s.zeros(2):
            raise ArithmeticError('open aperture edge identity failed')
    return {'H':Y,'rho':rho,'probes':probes,'instrument':instrument,'rows':rows,
            'gram':gram,'reference_after':apparatus.clean(reference),
            'classical_address_mixture_blank_each_probability':F(1,16),
            'is_primary_matrix_execution':False}


def implementation_control():
    """A faithful noninvariant state can have spatial unitaries without a fixed vacuum."""
    X=s.Matrix([[0,1],[1,0]]); rho=s.diag(s.Rational(3,4),s.Rational(1,4))
    root=apparatus.positive_sqrt(rho); omega=s.Matrix([root[i,j] for j in range(2) for i in range(2)])
    U=(s.eye(2)+s.I*X)/s.sqrt(2)
    V=s.kronecker_product(s.conjugate(U),U)
    P=s.kronecker_product(s.eye(2),X)-s.kronecker_product(X.T,s.eye(2))
    variance=s.simplify((omega.H*P*P*omega)[0])
    if variance != 2-s.sqrt(3): raise ArithmeticError('implementation control changed')
    return {'rho':rho,'source_rotation':U,'standard_GNS_unitary':V,
            'prepared_vector':omega,'standard_GNS_generator':P,
            'prepared_vector_fixed':apparatus.clean(V*omega-omega)==s.zeros(4,1),
            'H':s.zeros(2),'mean_ground_excess':0,
            'momentum_variance':variance,'is_primary_source':False}


def common_payload(member, addresses=4, spacing='1/4', half_width='1/100000'):
    contract.verify_member(member)
    m=apparatus.integer(addresses,'addresses',2)
    if m>16: raise ValueError('the authenticated explicit address menu is capped at16')
    h=apparatus.exact(spacing,'spacing',True)
    if (h/F(member['source']['spacing'])).denominator != 1:
        raise ValueError('exactly equally spaced source anchors require an integer lattice step')
    T=apparatus.exact(half_width,'filter half width',True)
    coordinates=[(-m*h/2+j*h,F(0),F(0)) for j in range(m+1)]
    descriptors=[algebra.probe(member,x,T,kind='plaquette') for x in coordinates]
    anchors=[tuple(p['nearest_lower_owner_cell']) for p in descriptors]
    if not all(p['active'] for p in descriptors) or len(set(anchors))!=m+1:
        raise ValueError('all aperture and shifted endpoint probes must be distinct active source cells')
    parameters=apparatus.precision_parameters(member,half_width=str(T),ramsey_time='0')
    # Direct algebraic services, no optional compiler error is silently consumed.
    parameters['compiler_error']='0'
    base=member['preparation_plan']['source']
    error=apparatus.error_budget(base['spectral_lower'],base['spectral_upper'],**parameters)
    delta=F(error['filtered_probe_operator_error_upper'])
    column=apparatus._relative_sqrt_upper(delta*delta+2*delta)
    floor=boundary.common_payload(member)['physical_density_floor_factored']
    return contract.serial({
        'member_id':member['member_id'],**member['ids'],
        'aperture_probes':descriptors[:-1],'shifted_endpoint_probe':descriptors[-1],
        'filter_parameters':parameters,'filter_errors':error,
        'direct_complete_Kraus_formula':'K_(b,k)=M^-1 sum_j exp(-2pi i k j/M) C_(b,j)',
        'branch_operators':['B_j','sqrt(I-B_j^2)'],
        'fine_outcomes':2*m,'neutral_address_dimension':m,'neutral_probe_dimension':2,
        'source_state_or_idle_reference_replaced':False,
        'all_failures_and_original_reference_retained':True,
        'complete_error_to_ideal_filtered_aperture_upper':min(F(1),column),
        'separate_unfiltered_original_blank_weight_each_label_lower':{
            'coefficient':F(1,90*m),'base':2,
            'negative_exponent':floor['negative_exponent']},
        'unfiltered_blank_loop_geometry':blank_loop_geometry(anchors[:-1]),
        'unfiltered_positive_bound_transferred_to_fixed_time_filter':False,
        'aperture_edge':aperture_bound(m,1,h),
        'open_source_incidence':open_box_geometry(member['source']['side']),
        'source_translation_unitary_applied':False,
        'actual_source_momentum_PVM_proved':False,
        'wavevector_labels_are_conserved_particle_momenta':False,
        'selected_preparation_momentum_concentration_proved':False,
        'spatial_state_invariance_necessary_for_every_implementation':False,
        'same_GNS_spatial_implementation_proved':False,
        'full_primary_probabilities_computed':False,
        'source_operator_description_workspace_or_gate_cost_bounded':False})


def bind_member(member, **kwargs):
    result=common_payload(member,**kwargs)
    return {'schema':'oph.primary_spatial_fourier_instrument.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member, receipt, **kwargs):
    if receipt!=bind_member(member,**kwargs):
        raise ValueError('spatial Fourier receipt differs from source, aperture or scope')


def certificate():
    control=exact_control(); implementation=implementation_control()
    return contract.serial({'schema':'oph.primary_momentum_prerequisites.v1',
        'selected_aperture':bind_member(contract.member(2)),
        'exact_control':{'original_probabilities':[str(r['p']) for r in control['rows']],
            'fine_outcomes':8,'reference_unchanged':True,'is_primary':False},
        'implementation_without_invariant_preparation':{
            'mean_ground_excess':0,'momentum_variance':str(implementation['momentum_variance']),
            'prepared_vector_fixed':False,'is_primary':False},
        'conserved_primary_particle_momentum_or_mass':None})
