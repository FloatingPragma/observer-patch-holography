"""Full-Fock collective residuals on the actual quantum-link source.

The bounds sharpen finite energy intervals. Their inverse-spacing scale
does not establish tightness on a fixed physical energy window.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import isqrt

import sympy as s

from model.particle_common_contract import contract,provenance
from model.particle_local_wilson_density import density
from model.particle_primary_response_moments import model as response
from model.particle_primary_refinement_response import model as refinement
from model.particle_primary_pair_channels.model import passive_window
from model.particle_primary_mass_readout import readout
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_spatial_ground_state import ground
from model.particle_primary_channel_labels import channels


SOURCE_INPUTS=(
    'model/particle_primary_energy_tightness/__init__.py',
    'model/particle_primary_energy_tightness/__main__.py',
    'model/particle_primary_energy_tightness/model.py',
    'model/particle_primary_energy_tightness/test_model.py',
    'model/particle_primary_energy_tightness/PROOF.md',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_refinement_response/PROOF.md',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_primary_origin_momentum/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_primary_channel_labels/channels.py',
    'model/particle_comparison/electromagnetism.py',
    'Lean/QFT/GroundRecordTail.lean',
)


def _clean(value):return s.simplify(s.expand_complex(value))


@lru_cache(maxsize=1)
def source_column_control():
    """Execute the actual Dirac Gram matrices and physical anchor sum rule."""
    refinement.native_control()
    gamma=density.gamma_data();beta=gamma[3]
    hops=tuple((-beta+beta*g)/2 for g in gamma[:3])
    pairs=response.beta_pairs();B=sum((p*m.H for p,m in pairs),s.zeros(4))
    if beta*B-B*beta!=2*B:raise ArithmeticError('onsite pair frequency changed')
    rows=[]
    # Per axis: negative neighbour, positive neighbour, or both. This
    # includes every actual site boundary class on boxes of side>=2.
    for sides in product((-1,1,0),repeat=3):
        legs=[]
        for axis,side in enumerate(sides):
            if side in (-1,0):legs.append(hops[axis])
            if side in (1,0):legs.append(hops[axis].H)
        degree=len(legs)
        gram=sum((h.H*h for h in legs),s.zeros(4))
        for p,m in pairs:
            if _clean((p.H*gram*p)[0])!=F(degree,2) or _clean((m.H*gram*m)[0])!=F(degree,2):
                raise ArithmeticError('source column norms changed')
        interior=all(side==0 for side in sides)
        if interior and gram!=3*s.eye(4):raise ArithmeticError('interior columns are not orthogonal')
        # Slater anchor quadratic norm, evaluated one bond at a time.
        negative=(s.eye(4)-beta)/2;occupied=s.diag(negative,negative)
        probe=s.diag(B,s.zeros(4));total=s.Rational(0);cross=s.Rational(0)
        for leg in legs:
            hopping=s.zeros(4).row_join(leg.H).col_join(leg.row_join(s.zeros(4)))
            R=hopping*probe-probe*hopping
            total+=_clean(s.trace(occupied*R.H*(s.eye(8)-occupied)*R))
            cross+=_clean(s.trace(occupied*probe.H*(s.eye(8)-occupied)*R))
        if total!=degree or cross!=0:raise ArithmeticError('source Slater residual identity changed')
        rows.append({'axis_neighbour_pattern':list(sides),'degree':degree,
            'interior':interior,'one_column_squared_norm':str(F(degree,2)),
            'unnormalized_one_native_column_anchor_residual_squared':str(total),
            'onsite_residual_cross_term':'0'})
    return {'actual_Dirac_and_Wilson_matrices_checked':True,
        'all_boundary_patterns':rows,'interior_neighbour_Gram':'3 I_4',
        'collective_beta_pair_count':12,'interior_residual_raising_pair_count':24,
        'interior_collective_normalization_squared':42,
        'interior_native6_dimensionless_residual_norm_squared':F(78,7),
        'physical_anchor_is_selected_state':False,
        'link_coordinates_replaced_in_selected_source':False}


def moment_coefficient(side,native_dimension):
    L=ground._i(side,'side',2);d=ground._i(native_dimension,'native dimension',1)
    V=L**3;B=3*L*L*(L-1);I=(L-2)**3
    # Boundary: 4d degree/(d+1); interior exact ladder:6(2d+1)/(d+1).
    return F(8*d*B-6*(2*d-1)*I,(d+1)*V)


def uniform_primary_coefficient():
    # L0=2. For L>=8, x=1/L<=1/8 and the derivative
    # of78+252x-792x²+528x³ is >=252-198=54.
    return F(3141,224)


def _ceil_sqrt(value):
    value=F(value)
    if value<0:raise ValueError('nonnegative square required')
    n=isqrt(value.numerator//value.denominator)
    return n if n*n>=value else n+1


@lru_cache(maxsize=1)
def density_commutator_control():
    """Exact source-Wilson density identity with a retained quantum link."""
    channels.native_number_control();rotation=channels.rotation_control()
    graph=channels.boundary_graph_control(3)
    g=density.gamma_data();beta=g[3];S=(-beta+beta*g[0])/2
    forward=s.zeros(8);forward[:4,4:]=S
    shift=s.Matrix([[0,0,0],[1,0,0],[0,1,0]])
    J=s.kronecker_product(shift,forward)+s.kronecker_product(shift.H,forward.H)
    nx=s.diag(*([1]*4+[0]*4));ny=s.eye(8)-nx
    X=s.kronecker_product(s.eye(3),nx);Y=s.kronecker_product(s.eye(3),ny)
    H=J+s.kronecker_product(s.eye(3),s.diag(3*beta,3*beta))
    H+=s.kronecker_product(s.diag(s.Rational(63,4),0,s.Rational(63,4)),s.eye(8))
    def comm(A,B):return A*B-B*A
    if comm(X,comm(H,X))+comm(Y,comm(H,Y))!=-2*J:
        raise ArithmeticError('source number-density double commutator changed')
    if comm(X+Y,H)!=s.zeros(24):
        raise ArithmeticError('source total number is not conserved')
    return {'actual_Wilson_matrix_and_quantum_link_checked':True,
        'full_Fock_lift_uses_exact_CAR_bilinear_identity':True,
        'complete_native_blocks':channels.native_number_control()['blocks'],
        'source_rotation_control':rotation,'physical_boundary_control':graph,
        'ambient_dangling_link_permutation_used':False,
        'comparison_control_state_is_selected_primary_state':False}


def density_positive_tail(signed_first_moment_upper,ground_excess,threshold):
    """Original positive tail; success effects must sum to at most identity."""
    bound=ground._q(signed_first_moment_upper,'signed first moment upper')
    eps=ground._q(ground_excess,'ground excess')
    delta=ground._q(threshold,'positive threshold',positive=True)
    return min(F(1),(bound+eps)/delta)


def selected_density_sum_rule(member):
    """Same-state negative hopping density and complete real density readout."""
    contract.verify_member(member);control=density_commutator_control()
    source=member['source'];L=source['side'];V=L**3;a=F(source['spacing'])
    B=source['interior_hopping_bonds'];P=source['interior_plaquettes']
    eps=F(member['state']['target_ground_excess']);seed=refinement.tile_seed()
    gain=F(seed['dimensionless_gain_per_cell'])
    hopping=-gain-F(P,7*V)-a*eps/V
    uniform=-gain-F(3,7)-F(1,128)
    if hopping<uniform or uniform<=0:
        raise ArithmeticError('selected hopping density floor failed')
    level=source['level'];normalization=30*V
    if level==0:phase_lower=phase_upper=F(2)
    elif level==1:phase_lower=phase_upper=F(1)
    else:phase_lower=2*a*a;phase_upper=min(F(2),5*a*a)
    first_lower=phase_lower*V*hopping/(3*a*normalization**2)
    first_upper=phase_upper*10*B/(a*normalization**2)
    simple_upper=a/(6*V)
    if first_lower<=0 or first_lower>first_upper or first_upper>simple_upper:
        raise ArithmeticError('complete density first moment bounds failed')
    return {'source_control':control,'tile_seed_member_id':seed['member_id'],
        'tile_seed_result_sha256':seed['producer_result_sha256'],
        'selected_negative_dimensionless_hopping_density_lower':hopping,
        'uniform_negative_dimensionless_hopping_density_lower':uniform,
        'selected_negative_hopping_expectation_lower':V*hopping/a,
        'three_direction_hopping_expectations_equal':True,
        'number_operator':'N_x = sum of all 60 CAR occupations at x',
        'centered_number_operator':'n_x = N_x - 30 I',
        'phase_step':'pi * source_spacing',
        'one_minus_cos_phase_enclosure':[phase_lower,phase_upper],
        'normalization':normalization,'complete_outcomes':['cosine','sine','failure'],
        'random_setting_factor':1,'combined_success_effect_at_most_identity':True,
        'original_signed_first_moment_lower':first_lower,
        'original_signed_first_moment_upper':first_upper,
        'simple_signed_first_moment_upper':simple_upper,
        'original_negative_absolute_first_moment_upper':eps,
        'original_positive_first_moment_upper':first_upper+eps,
        'positive_weight_above_one_upper':density_positive_tail(first_upper,eps,1),
        'every_fixed_positive_energy_tail_tends_to_zero':True,
        'tail_conclusion_domain':'this volume-normalized total-density instrument',
        'source_reference_and_failures_retained':True,
        'instrument_preserves_input_state':False,
        'state_replaced_by_tile_trial':False,'full_histogram_executed':False,
        'phase_label_identified_with_conserved_momentum':False,
        'fixed_positive_band_noncollapse_proved':False}


def retained_origin_bound(member,inherited):
    """Same recorded origin for the source probe and both energy clocks."""
    from model.particle_primary_origin_momentum import model as origin
    geometry=origin.origin_geometry(member['source']['level'])
    V=member['source']['side']**3;count=geometry['origin_count']
    if geometry['source_side']!=member['source']['side'] or geometry['source_cells']!=V:
        raise ArithmeticError('retained origin and primary source differ')
    D=F(inherited['all_block_defects']['defect_count_upper'])
    chi=max(F(),(12-D/count)/42)
    a=F(member['source']['spacing']);beta=F(member['state']['beta'])
    interior=geometry['collar_cells']>=1
    coefficient=F(78,7) if interior else moment_coefficient(member['source']['side'],6)
    window=None
    if chi>0:
        lower=chi/(4*beta);upper=F(6+_ceil_sqrt(4*coefficient/chi))/a
        window=passive_window(chi,beta,6/a,coefficient/a**2,lower,upper)
        if window['original_transfer_weight_lower']<chi/2:
            raise ArithmeticError('retained-origin passive window budget failed')
    return {'origin_geometry':geometry,'same_origin_register_for_complete_history':True,
        'source_density_unchanged':True,'reference_idle_and_retained':True,
        'direct_selected_commutator_lower':chi,
        'centered_moment_coefficient':coefficient,'original_window':window,
        'one_success_and_one_failure_per_retained_origin':True,
        'averaged_observable_substituted_for_origin_history':False,
        'fixed_energy_uniform_weight_proved':False}


def selected_bound(member):
    contract.verify_member(member);source_column_control()
    inherited=refinement.bind_member(member);refinement.verify_bound(member,inherited)
    value=inherited['result'];source=member['source']
    L=source['side'];a=F(source['spacing']);beta=F(member['state']['beta'])
    eps=F(member['state']['target_ground_excess'])
    coefficient=moment_coefficient(L,6)
    if coefficient>uniform_primary_coefficient():raise ArithmeticError('uniform source moment failed')
    lower=1/(100*beta);upper=F(46)/a
    window=passive_window(F(value['collective_commutator_lower']),beta,6/a,coefficient/a**2,lower,upper)
    if window['original_transfer_weight_lower']<F(3,100):raise ArithmeticError('uniform narrowed band failed')
    final=readout.excitation_window((lower,upper),window['original_transfer_weight_lower'],eps,100*eps)
    q6_chi=F(value['q6_averaged_matched_pair_weight_lower']);q6=None
    if q6_chi>0:
        q6=passive_window(q6_chi,beta,6/a,moment_coefficient(L,1)/a**2,lower,F(24)/a)
    clock=spectroscopy.clock_plan(member,lower/4,'1/400','1/400')
    clock_loss=F(1,200)+F(clock['two_complete_clock_error_upper'])
    return contract.serial({'schema':'oph.primary_energy_tightness_result.v1',
        'member_id':member['member_id'],**member['ids'],
        'refinement_producer_result_sha256':inherited['result_sha256'],
        'source_column_control':source_column_control(),
        'dimensionless_centered_moment_upper':coefficient,
        'uniform_dimensionless_centered_moment_upper':uniform_primary_coefficient(),
        'native6_original_window':window,
        'original_weight_simple_lower':F(3,100),
        'final_excitation':final,'final_weight_simple_lower':F(1,50),
        'q6_scalar_original_window':q6,
        'finite_clock_plan':clock,
        'finite_reported_interval':[lower/2,upper+lower/2],
        'finite_reported_original_weight_lower':window['original_transfer_weight_lower']-clock_loss,
        'finite_reported_weight_simple_lower':F(1,50),
        'retained_origin_response':retained_origin_bound(member,value),
        'selected_density_sum_rule':selected_density_sum_rule(member),
        'physical_anchor_average_centered_residual_squared':F(2*source['interior_hopping_bonds'],7*L**3)/a**2,
        'physical_anchor_is_selected_state':False,
        'source_operator_norm_method_can_remove_inverse_spacing':False,
        'selected_state_weighted_residual_lower_claimed':False,
        'source_or_preparation_or_clock_changed':False,
        'actual_histogram_executed':False,'fixed_energy_spectral_tightness_proved':False,
        'isolated_branch_dispersion_or_particle_identification':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_energy_tightness_bound.v1','result':result,
        'result_sha256':contract.digest(result)}


def verify_bound(member,receipt):
    if receipt!=bind_member(member):raise ValueError('energy band or its source changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*refinement.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_energy_tightness_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2)],
        'retained_origin_member':bind_member(contract.member(7)),
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('energy inputs changed during replay')
    return result
