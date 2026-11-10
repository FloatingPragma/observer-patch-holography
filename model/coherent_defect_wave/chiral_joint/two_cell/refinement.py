"""Exact transfer failures and physical momentum controls for this source.

Finite comparison bounds use the complete residual, including off-image
components. No periodic ground or continuum spectrum is imported from the
open control.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import sympy as s
from ..model import clean
from .hamiltonian import column,local_yukawa,apply_state
from .basis import local_current_column
from .symmetries import (periodic_column,periodic_translation,periodic_current_column,
                         apply_column,norm2,subtract)


def blank_extension(state):
    return {(mask,link,0,hl,hr):value for (mask,link,hl,hr),value in state.items()}


def blank_compression(state):
    return {(mask,l01,hl,hr):value for (mask,l01,l10,hl,hr),value in state.items() if l10==0}


def full_transfer_residual(state,clock_scale=s.Integer(1),energy_shift=s.Integer(0)):
    target=apply_column(blank_extension(state),periodic_column)
    reference=apply_state(state)
    reference=clean({key:clock_scale*reference.get(key,0)+energy_shift*state.get(key,0)
                     for key in set(reference)|set(state)})
    return subtract(target,blank_extension(reference))


@lru_cache(maxsize=1)
def exact_transfer_control():
    # The actual local quantum pair/scalar state, not a naked charged mode.
    state={(mask,0,scalar,0):value/2 for (mask,scalar),value in local_yukawa(0,0).items()}
    assert norm2(state)==1
    initial=blank_extension(state)
    for site in range(2):
        for current in range(12):
            assert not apply_column(state,lambda key:local_current_column(key,site,current))
            assert not apply_column(initial,lambda key:periodic_current_column(key,site,current))
    residual=full_transfer_residual(state)
    outside={key:value for key,value in residual.items() if key[2]!=0}
    assert not blank_compression(residual)
    assert residual==outside and norm2(residual)==s.Rational(3,16)
    image=apply_column(initial,periodic_column)
    leak={key:value for key,value in image.items() if key[2]!=0}
    assert leak==residual
    # The omitted component survives every constant clock scale and offset.
    z,w=s.symbols('clock_scale energy_shift',real=True)
    arbitrary=full_transfer_residual(state,z,w)
    assert {key:value for key,value in arbitrary.items() if key[2]!=0}==outside
    time=s.Rational(1,10**6);bound=s.Integer(30)
    remainder=time**2*bound**2/2
    amplitude=time*s.sqrt(3)/4-remainder
    assert amplitude.is_positive is True
    lower=s.expand(amplitude**2)
    return {'open_dimension':614400,'periodic_dimension':3686400,
            'embedding':'append the blank second quantum link; retain both scalar memories',
            'input':'normalized actual left quantum-Yukawa pair/scalar state',
            'open_and_periodic_local_Gauss_checks':48,
            'compressed_residual_squared':'0','full_residual_squared':'3/16',
            'entire_residual_is_omitted_link_excitation':True,
            'off_image_residual_independent_of_constant_clock_scale_and_energy_shift':True,
            'full_domain_residual_norm_upper_bound':'15/2',
            'uniform_complete_readout_error_upper_bound':'min(1,15*abs(t)/2)',
            'periodic_full_H_norm_bound':'30','time':str(time),
            'first_order_unitary_remainder':str(remainder),
            'actual_new_link_yes_probability_lower_bound':str(lower),
            'actual_new_link_no_probability':'1-p_yes',
            'embedded_open_new_link_yes_probability':'0',
            'old_ground_or_spectrum_transferred_to_periodic':False}


def compensated_translation(key):
    # The last coordinate is the retained position of a fixed charge -1.
    *source,position=key
    if position not in (0,1):raise ValueError('compensator position must be0 or1')
    return {(*target,1-position):value for target,value in periodic_translation(tuple(source)).items()}


def compensated_current(key,site,current):
    *source,position=key
    boundary=[0,0];boundary[position]=-1
    return {(*target,position):value for target,value in
            periodic_current_column(tuple(source),site,current,tuple(boundary)).items()}


def compensated_hamiltonian(key):
    *source,position=key
    return {(*target,position):value for target,value in periodic_column(tuple(source)).items()}


def compensated_projector(state,label):
    if label not in range(4):raise ValueError('four complete momentum-phase outcomes')
    current=state;out=defaultdict(lambda:s.Integer(0))
    for power in range(4):
        coefficient=(-s.I)**(label*power)/4
        for key,value in current.items():out[key]+=coefficient*value
        current=apply_column(current,compensated_translation)
    return clean(out)


@lru_cache(maxsize=1)
def compensated_momentum_control():
    # Source e^c charge+1 and explicit colocated compensator charge-1.
    swap=s.Matrix([[0,1],[1,0]])
    qleft=s.diag(-1,0);qright=s.diag(0,-1)
    assert swap*qleft*swap.H==qright and swap*qright*swap.H==qleft
    state={(1<<6,0,0,0,0,1):s.Integer(1)}
    assert norm2(state)==1
    current=state;orbit=[]
    for _ in range(4):
        orbit.extend(current)
        for site in range(2):
            for generator in range(12):
                assert not apply_column(current,lambda key:compensated_current(key,site,generator))
        current=apply_column(current,compensated_translation)
    assert current==state
    twice=apply_column(apply_column(state,compensated_translation),compensated_translation)
    assert twice=={key:-value for key,value in state.items()}
    for key in orbit:
        assert apply_column(compensated_translation(key),compensated_hamiltonian)==apply_column(
            compensated_hamiltonian(key),compensated_translation)
    images=[compensated_projector(state,j) for j in range(4)]
    values=[s.simplify(norm2(image)) for image in images]
    assert values==[0,s.Rational(1,2),0,s.Rational(1,2)]
    total=defaultdict(lambda:s.Integer(0))
    for j,image in enumerate(images):
        for key,value in image.items():total[key]+=value
        for k in range(4):
            projected=compensated_projector(image,k)
            assert projected==(image if j==k else {})
        for site in range(2):
            for generator in range(12):
                assert not apply_column(image,lambda key:compensated_current(key,site,generator))
    assert clean(total)==state
    # In the principal branch, (2a/pi)P=Q_1-Q_3 on the odd source sector.
    scaled=subtract(images[1],images[3])
    squared=subtract(compensated_projector(scaled,1),compensated_projector(scaled,3))
    assert squared==state
    return {'full_source_plus_position_register_dimension':7372800,
            'compensator_position_states':2,'compensator_charge':'-1',
            'source_input':'one right spin-up e^c occupation, blank links and scalar memories',
            'reference_input':'fixed charge-1 at the same right cell; its position is retained',
            'joint_translation':'source antiperiodic shift and compensator-position swap',
            'joint_H':'actual periodic source H tensor identity on the position register',
            'all_24_Gauss_constraints_hold_for_input_and_every_phase_branch':True,
            'joint_translation_preserves_the_physical_Gauss_domain':True,
            'complete_phase_labels':['0','pi/2','pi','3pi/2'],
            'complete_probabilities':list(map(str,values)),
            'odd_sector_principal_momenta':['pi/(2a)','-pi/(2a)'],
            'odd_sector_P_squared':'pi^2/(4a^2) * identity',
            'source_and_reference_translated_together':True,
            'phase_observable':'joint source/reference translation',
            'source_only_particle_momentum_identified':False,
            'isolated_particle_pole_identified':False,
            'source_input_is_stationary_ground':False,
            'energy_or_relativistic_dispersion_measured_in_this_control':False}


def primitive_link_variance_control():
    # One actual e^c link channel acts on blank and its raised state.
    from ..model import active_terms
    assert active_terms()[3][1]==(0,0)
    r=s.symbols('r',nonnegative=True)
    u=s.Matrix([[0,0],[1,0]])
    combined=(u-r*s.eye(2)).H*(u-r*s.eye(2))+(u.H-r*s.eye(2)).H*(u.H-r*s.eye(2))
    difference=combined-r*r*s.eye(2)
    assert {s.factor(k):v for k,v in difference.eigenvals().items()}=={(r-1)**2:1,(r+1)**2:1}
    # On the other four actual link states U=Udagger=0, leaving r^2 I
    # after subtracting the claimed lower bound r^2.
    return {'actual_channel':'U_ec,ec=|raised><blank| in the six-state link',
            'both_residual_lower_bound':'||(U-u)psi||^2+||(Udagger-conj(u))psi||^2 >= |u|^2',
            'two_channel_eigenvalues_after_subtracting_bound':['(|u|-1)^2','(|u|+1)^2'],
            'four_other_link_state_eigenvalues_after_subtracting_bound':'|u|^2',
            'nonzero_c_number_transfer_can_have_vanishing_relative_residual':False,
            'scope':'fixed primitive link and product-background embedding; arbitrary dressed/interacting limits are not excluded'}


def odd_ring_momenta(side):
    if type(side) is not int or side<2 or side%2:
        raise ValueError('positive even number of cells at least2')
    # k/(pi/a), with antiperiodic field phases and the principal branch.
    return tuple(s.Rational(2*n+1,side) for n in range(-side//2,side//2))


def momentum_design(side):
    momenta=odd_ring_momenta(side);squares=tuple(sorted(set(k*k for k in momenta)))
    design=s.Matrix([[1,q] for q in squares])
    return {'side':side,'k_in_pi_over_a_units':list(map(str,momenta)),
            'distinct_squared_momenta':list(map(str,squares)),
            'mass_squared_speed_squared_design_rank':design.rank(),
            'independent_mass_shell_tests_after_two_parameter_fit':max(0,len(squares)-2),
            'executed_rest_momentum_present':False,
            'this_is_quantized_menu_not_computed_energy_spectrum':True}


def _rational(value):
    value=s.sympify(value)
    if not value.is_Rational or value.has(s.Float):
        raise ValueError('exact rational input required')
    return value


def shell_feasibility(samples):
    """Exact interval feasibility of E^2=M+C*p^2 with M,C nonnegative.

    samples contain (p_squared, E_squared_lower, E_squared_upper). An input
    is a specified energy-square envelope, not a source spectrum invented
    by this routine.
    """
    values=[tuple(map(_rational,row)) for row in samples]
    if not values or any(x<0 or lo<0 or hi<lo for x,lo,hi in values):
        raise ValueError('nonempty nonnegative exact squared-energy intervals required')
    constraints=[(-s.Integer(1),s.Integer(0),s.Integer(0)),
                 (s.Integer(0),-s.Integer(1),s.Integer(0))]
    for x,lo,hi in values:constraints.extend(((1,x,hi),(-1,-x,-lo)))
    vertices=set()
    for (a,b,c),(d,e,f) in combinations(constraints,2):
        determinant=a*e-b*d
        if determinant==0:continue
        mass=(c*e-b*f)/determinant;speed=(a*f-c*d)/determinant
        if all(A*mass+B*speed<=C for A,B,C in constraints):vertices.add((mass,speed))
    bounded=any(x>0 for x,_,_ in values)
    return {'feasible':bool(vertices),'bounded_parameter_region':bounded,
            'parameter_vertices':[[str(m),str(c)] for m,c in sorted(vertices)],
            'unique_parameters_on_this_menu':bounded and len(vertices)==1,
            'common_positive_speed_possible':bool(vertices) and
                (not bounded or any(c>0 for _,c in vertices)),
            'physical_particle_identification_implied':False}


def rational_error_budget(time,hbar,residual,initial=0,tail=0,effect=0,numerical=0):
    t,h,e,i,r,d,n=map(_rational,(time,hbar,residual,initial,tail,effect,numerical))
    if h<=0 or min(e,i,r,d,n)<0:raise ValueError('positive hbar and nonnegative errors required')
    return min(s.Integer(1),i+abs(t)*e/h+2*r+d+n)


def certificate_data():
    # Algebraic acceptance controls are labeled and never used as data from H.
    feasible=shell_feasibility([(1,5,5),(9,13,13),(25,29,29)])
    inconsistent=shell_feasibility([(1,5,5),(9,13,13),(25,30,30)])
    assert feasible['parameter_vertices']==[['4','1']] and not inconsistent['feasible']
    return {'source_transfer':exact_transfer_control(),
            'physical_charged_momentum':compensated_momentum_control(),
            'primitive_quantum_link_background':primitive_link_variance_control(),
            'antiperiodic_momentum_menus':[momentum_design(n) for n in (2,4,6)],
            'dispersion_validator_controls':{'algebraic_only_not_source_energies':True,
                'compatible_three_point_control':feasible,'incompatible_third_point_control':inconsistent},
            'same_source_complete_readout_bound':'min(1,initial+|t|*full_residual/hbar+2*tail+effect+numerical)',
            'finite_theorem_crosswalk':{
                'Lean/QFT/EffectiveQuantumComparison.lean':['interaction_picture_error','effect_error_of_interaction_residual','propagation_error_with_tail'],
                'Lean/QFT/PeriodicMomentum.lean':['polynomial_conserved','polynomials_commute'],
                'Lean/Geometry/MassShellKinematics.lean':['fourMomentum_time_sq_eq_mass_sq_add_spatial','fourMomentum_time_eq_mass_iff_spatial_zero'],
                'Lean/QFT/FiniteUnitaryScatteringNoGo.lean':['finite_unitary_ambient_powers_have_no_limit','identical_relative_evolution_tendsto']},
            'theorem_application_scope':{'finite_residual_and_complete_effect_bounds':'supplied for this finite comparison',
                'translation_polynomial_conservation':'supplied on the periodic joint source/reference domain',
                'mass_shell_physical_energy_momentum_assignment':'not supplied by the kinematic theorem',
                'asymptotic_scattering_reference_and_limit':'not supplied by the finite-unitary theorem'},
            'scope':{'actual_physical_joint_translation_apparatus_supplied':True,
                     'open_spectrum_transferred_to_periodic':False,'periodic_ground_energy_computed':False,
                     'relativistic_mass_shell_verified':False,'interacting_continuum_constructed':False,
                     'three_spatial_translation_directions_constructed_here':False}}
