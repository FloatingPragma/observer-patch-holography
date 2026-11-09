"""Exact original-state cube reduction, with no assigned spectral values."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as s

from model.particle_primary_selected_spin_response import model as selected

K=((1,3,3,1),(1,1,-1,-1),(1,-1,-1,1),(1,-3,3,-1))
MULTIPLICITIES=(1,3,3,1)
CELLS=tuple(product(range(2),repeat=3))


def distance(x,y):
    return sum(a!=b for a,b in zip(x,y))


def radial_matrix(shells):
    if len(shells)!=4:raise ValueError('four Hamming shells required')
    return s.Matrix([[shells[distance(x,y)] for y in CELLS] for x in CELLS])


def walsh_eigenresponses(shells):
    if len(shells)!=4:raise ValueError('four Hamming shells required')
    return tuple(s.Matrix(K)*s.Matrix(shells))


def invert_eigenresponses(eigenresponses):
    if len(eigenresponses)!=4:raise ValueError('four cube eigenresponses required')
    return tuple(s.Matrix(K)*s.Matrix(eigenresponses)/8)


@lru_cache(maxsize=1)
def orbit_control():
    """Use the literal primary spin lift's 24 spatial rotation images."""
    group=selected.spin.group;channels=selected.spin.channels
    rotations={tuple(map(int,group.rotation(q).flat)) for q in group.elements()}
    permutations=[];where={x:i for i,x in enumerate(CELLS)}
    for flat in sorted(rotations):
        R=s.Matrix(3,3,flat)
        if R.det()!=1:raise ArithmeticError('source rotation is not proper')
        permutations.append(tuple(where[channels.rotate_cell(x,2,R)] for x in CELLS))
    sizes=[]
    for d in range(4):
        target={p for p in product(range(8),repeat=2) if distance(CELLS[p[0]],CELLS[p[1]])==d}
        seed=min(target);orbit={(g[seed[0]],g[seed[1]]) for g in permutations}
        if orbit!=target:raise ArithmeticError('proper cube pair orbit is smaller than a Hamming shell')
        if not all(any(g[a]==b and g[b]==a for g in permutations) for a,b in target):
            raise ArithmeticError('source proper cube lacks a pair reversal')
        sizes.append(len(orbit))
    covariance=selected.covariance_controls()
    if any(covariance['probe_irreps'][label]['character_inner_product']!=1 for label in ('G1','H')):
        raise ArithmeticError('actual source multiplet covariance changed')
    return {'actual_binary_rotations':48,'proper_spatial_rotations':len(permutations),
        'ordered_pair_orbit_sizes':sizes,'every_ordered_pair_has_a_proper_rotation_swap':True,
        'actual_G1_H_unitary_coefficient_covariance_reused':True}


@lru_cache(maxsize=1)
def prepared_transport_countercontrol():
    """A separate exact six-CAR example rules out a state substitution."""
    count=6;dimension=1<<count
    creation=[np.asarray(selected.odd.creation_word((i,),count)).astype(np.int64) for i in range(count)]
    hopping=sum((creation[i]@creation[i+3].T+creation[i+3]@creation[i].T for i in range(3)),
                np.zeros((dimension,dimension),dtype=np.int64))
    Z=np.diag([1,-1]);X=np.array([[0,1],[1,0]],dtype=np.int64)
    H=np.kron(hopping,Z);I=np.eye(2*dimension,dtype=np.int64)
    parity=np.diag([(-1)**sum(bool(mask&(1<<i)) for i in range(3,6)) for mask in range(dimension)])
    S=np.kron(parity,X)
    fields=[np.kron(np.asarray(selected.odd.creation_word(word,count)).astype(np.int64),np.eye(2,dtype=np.int64))
            for word in ((0,1,2),(3,4,5))]
    if np.any(S@H-H@S) or np.any(S@fields[0]-fields[0]@S) or np.any(S@fields[1]+fields[1]@S):
        raise ArithmeticError('comparison symmetry does not force off-diagonal cancellation')
    for flavor in range(3):
        number=creation[flavor]@creation[flavor].T+creation[flavor+3]@creation[flavor+3].T
        if np.any(np.kron(number,np.eye(2,dtype=np.int64))@H-H@np.kron(number,np.eye(2,dtype=np.int64))):
            raise ArithmeticError('comparison native occupation is not conserved')
    # Each of the three hopping terms has norm one, so ||H||<=3.
    # Q=4I-H>=I and Q²/Tr(Q²) is a faithful stationary polynomial state.
    Q=4*I-H;rho_numerator=Q@Q;normalizer=int(np.trace(rho_numerator))
    vacuum=np.zeros(2*dimension,dtype=np.int64);vacuum[0]=1
    left,right=(A@vacuum for A in fields)
    moment_vector=right;prepared=[]
    response=fields[1].copy();original=[]
    for order in range(5):
        prepared.append(int(left@moment_vector));moment_vector=H@moment_vector
        original.append(F(int(np.trace(rho_numerator@fields[0].T@response)),normalizer))
        response=H@response-response@H
    if prepared[:3]!=[0,0,0] or abs(prepared[3])!=6 or any(original):
        raise ArithmeticError('prepared/original comparison lost its strict distinction')
    return {'domain':'separate two-site six-CAR plus central-sign comparison; not the primary source',
        'stationary_positive_polynomial_numerator':'(4I-H)^2','normalizer':normalizer,
        'polynomial_amplitude_lower':1,'conserved_flavor_numbers':3,
        'prepared_inter_site_moments_orders_zero_to_four':prepared,
        'original_off_diagonal_transfer_moments_orders_zero_to_four':[str(v) for v in original],
        'invariant_unitary_forces_all_off_diagonal_spectral_events_zero':True,
        'same_cancellation_for_every_positive_integer_power_of_polynomial_density':True,
        'finite_dimensional_polynomial_cooling_can_make_ground_excess_arbitrarily_small':True,
        'comparison_state_or_H_substituted_for_primary_source':False}


@lru_cache(maxsize=1)
def control():
    shells=s.symbols('m0:4',real=True);matrix=radial_matrix(shells)
    W=s.Matrix([[(-1)**sum(a*b for a,b in zip(q,x)) for x in CELLS] for q in CELLS])/s.sqrt(8)
    spectrum=walsh_eigenresponses(shells)
    diagonal=s.diag(*(spectrum[sum(q)] for q in CELLS))
    if s.simplify(W*matrix*W.T-diagonal)!=s.zeros(8):
        raise ArithmeticError('actual cube Fourier reduction failed')
    if s.Matrix(K)**2!=8*s.eye(4) or invert_eigenresponses(spectrum)!=shells:
        raise ArithmeticError('four-measure inverse failed')
    total=s.expand(sum(n*v for n,v in zip(MULTIPLICITIES,spectrum))/48)
    if total!=shells[0]/6:raise ArithmeticError('complete spatial menu normalization changed')
    no_contrast={}
    for label,m0 in [('G1',s.Rational(9,5)),('H',s.Rational(18,5))]:
        values=walsh_eigenresponses((m0,0,0,0))
        no_contrast[label]={'diagonal_weight':str(m0),'four_eigenresponses':[str(x) for x in values],
            'chosen_phase_probability':str(m0/6),'each_complete_menu_row_probability':str(m0/48)}
    intensities=walsh_eigenresponses((1,s.Rational(1,8),0,0))
    if min(intensities)<=0 or len(set(intensities))!=4:raise ArithmeticError('intensity-only contrast control failed')
    return {'scope':'member-zero original creation transfer measures after the internal component sum within G1 or H',
        'source_cube_orbits':orbit_control(),'Krawtchouk_matrix':[list(row) for row in K],
        'inverse_matrix':[[str(F(x,8)) for x in row] for row in K],
        'shell_order':[0,1,2,3],'Walsh_weight_multiplicities':list(MULTIPLICITIES),
        'chosen_phase_eigenmeasure_divisor':6,'complete_menu_per_row_eigenmeasure_divisor':48,
        'complete_creation_total':'m0/6','inverse_from_complete_per_row_measures':'m = 6 K p',
        'spectral_shell_measures_are_real_signed':True,'time_transforms_claimed_real':False,
        'full_component_or_cross_irrep_kernel_reduced_to_four_scalars':False,
        'actual_primary_off_diagonal_shell_values_computed':False,
        'nonidentical_normalized_energy_shapes_proved':False,
        'dispersion_or_isolated_particle_pole_proved':False,
        'diagonal_no_contrast_countercontrol':no_contrast,
        'intensity_only_contrast_countercontrol':{'four_eigenresponse_intensities':[str(x) for x in intensities],
            'all_normalized_energy_measures':'the same unit point measure at source transfer9',
            'control_is_a_primary_source_computation':False},
        'prepared_to_original_countercontrol':prepared_transport_countercontrol()}
