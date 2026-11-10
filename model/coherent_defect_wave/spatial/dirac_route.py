"""Declared Clifford routing alternative; distinct from the square-root SWAP.

The finite one-particle shift is executed in position space. Its mass angle
and four-component routing coin are construction inputs, not native masses.
"""
from functools import lru_cache
from itertools import product

import sympy as s

from model.common_source.overlap_regulator import gamma_matrices
from ..model import JointObserver
from ..gauss import gamma_amplitudes
from .family import expanded,weight,positive_local


@lru_cache(maxsize=1)
def clifford():
    raw,_ = gamma_matrices()
    def exact(value):
        real,imag = int(value.real),int(value.imag)
        if value != complex(real,imag):
            raise ValueError('consumed Clifford entries must be Gaussian integers')
        return s.Integer(real)+s.I*imag
    gamma = [s.Matrix(g).applyfunc(lambda v:exact(complex(v))) for g in raw]
    beta = gamma[3]
    alpha = tuple((s.I*beta*g).applyfunc(s.expand) for g in gamma[:3])
    generators = (*alpha,beta)
    for i,left in enumerate(generators):
        assert left == left.H
        for j,right in enumerate(generators):
            assert left*right+right*left == (2*s.eye(4) if i == j else s.zeros(4))
    return alpha,beta


def symbol(momentum,mass_angle=0):
    alpha,beta = clifford()
    out = s.eye(4)
    for k,a in zip(momentum,alpha,strict=True):
        out = (s.cos(k)*s.eye(4)-s.I*s.sin(k)*a)*out
    return ((s.cos(mass_angle)*s.eye(4)-s.I*s.sin(mass_angle)*beta)*out).applyfunc(s.expand)


def route_step(amplitudes,length,mass_angle=0,inverse=False):
    """Execute three nearest-neighbor spin-dependent shifts and one coin gate."""
    if type(length) is not int or length < 3:
        raise ValueError('periodic length at least three required')
    if type(inverse) is not bool:
        raise ValueError('inverse flag must be boolean')
    angle = s.sympify(mass_angle)
    if (isinstance(mass_angle,(float,complex)) or angle.has(s.Float)
            or angle.is_number is not True or angle.is_finite is not True
            or angle.is_real is not True):
        raise ValueError('exact finite real mass angle required')
    mass_angle = angle
    state = {}
    for (x,coin),value in amplitudes.items():
        if (len(x) != 3 or any(type(v) is not int or not 0 <= v < length for v in x)
                or type(coin) is not int or not 0 <= coin < 4):
            raise ValueError('a position in the torus and a four-component coin are required')
        number = s.sympify(value)
        if (isinstance(value,(float,complex)) or number.has(s.Float)
                or number.is_number is not True or number.is_finite is not True):
            raise ValueError('finite exact amplitudes required')
        if number != 0:
            state[(tuple(x),coin)] = number
    alpha,beta = clifford()
    def matrix_at_sites(source,matrix):
        out = {}
        for (x,column),value in source.items():
            for row in range(4):
                if matrix[row,column] != 0:
                    key = (x,row)
                    out[key] = expanded(out.get(key,0)+matrix[row,column]*value)
        return {key:value for key,value in out.items() if value != 0}
    def shift(source,axis,reverse):
        out = {}
        for sign in (1,-1):
            projector = (s.eye(4)+sign*alpha[axis])/2
            projected = matrix_at_sites(source,projector)
            for (x,row),value in projected.items():
                y = list(x)
                y[axis] = (y[axis]+sign*(-1 if reverse else 1))%length
                key = (tuple(y),row)
                out[key] = expanded(out.get(key,0)+value)
        return {key:value for key,value in out.items() if value != 0}
    coin = s.cos(mass_angle)*s.eye(4)-s.I*s.sin(mass_angle)*beta
    if inverse:
        state = matrix_at_sites(state,coin.H)
        for axis in (2,1,0):
            state = shift(state,axis,True)
    else:
        for axis in range(3):
            state = shift(state,axis,False)
        state = matrix_at_sites(state,coin)
    return state


def certificate():
    alpha,beta = clifford()
    p = s.symbols('p0:3',real=True)
    mass = s.symbols('m',real=True)
    h = mass*beta+sum((x*a for x,a in zip(p,alpha,strict=True)),s.zeros(4))
    assert (h*h-(mass**2+sum(x*x for x in p))*s.eye(4)).applyfunc(s.expand) == s.zeros(4)
    initial = {((0,0,0),0):s.Integer(1),((0,0,0),3):s.I}
    result = route_step(initial,4,s.pi/3)
    assert route_step(result,4,s.pi/3,inverse=True) == initial
    assert s.expand(sum(weight(v) for v in result.values())) == 2
    k = (s.pi/2,0,s.pi)
    phases = (s.I,1,-1)
    target = symbol(k,s.pi/3)
    for column in range(4):
        lifted = {(x,column):s.prod(phases[i]**x[i] for i in range(3))
                  for x in product(range(4),repeat=3)}
        evolved = route_step(lifted,4,s.pi/3)
        for x in product(range(4),repeat=3):
            phase = s.prod(phases[i]**x[i] for i in range(3))
            for row in range(4):
                assert s.expand(evolved.get((x,row),0)-phase*target[row,column]) == 0
    # Independent numerical controls test the proved product-formula bound.
    # These 4x4 calculations never diagonalize the 72,000-state native cell.
    import numpy as np
    from scipy.linalg import expm
    native_alpha = [np.array(a,complex) for a in alpha]
    native_beta = np.array(beta,complex)
    coefficients = (s.Rational(1,5),s.Rational(-1,4),s.Rational(1,3),s.Rational(1,3))
    components = [float(x)*a for x,a in zip(coefficients,(*native_alpha,native_beta),strict=True)]
    exact = expm(-1j*sum(components))
    constant = sum(abs(coefficients[i]*coefficients[j]) for i in range(4) for j in range(i+1,4))
    bounds = []
    errors = []
    for ticks in (50,100,200):
        one = np.eye(4,dtype=complex)
        for component in components:
            one = expm(-1j*component/ticks)@one
        error = float(np.linalg.norm(np.linalg.matrix_power(one,ticks)-exact,2))
        bound = constant/ticks
        assert error < float(bound)
        errors.append(error)
        bounds.append(str(bound))
    assert errors[2] < errors[1] < errors[0]
    # One fermion occupies either this cell or an untouched reference cell.
    # The fixed-number spatial qubit controls U without parity superposition.
    observer = JointObserver()
    native = positive_local(observer.basis_state(1,5,internal=3))
    moved = observer.step(native)
    assert moved.norm_squared() == 1
    assert gamma_amplitudes(moved.as_dict()) == moved.as_dict()
    overlap = expanded(sum(s.conjugate(value)*moved.as_dict().get(key,0)
                           for key,value in native.amplitudes)/native.norm_denominator)
    purity = expanded((1+weight(overlap))/2)
    assert purity < 1
    return {'schema':'oph.coherent_defect_wave.clifford_route.v1',
            'role':'declared alternative one-particle routing, not the square-root-SWAP spectrum',
            'clifford_source':'model/common_source/overlap_regulator.py:gamma_matrices',
            'coin_dimension':4,'exact_Clifford_pair_checks':16,
            'literal_position_space_Bloch_column_checks':4,
            'literal_inverse_and_norm_checked':True,
            'finite_shift_radius_per_tick':'one site per coordinate; three in Manhattan distance',
            'continuum_clock_convention':'spacing a and tick duration a, in chosen units c=1',
            'continuum_one_particle_generator':'sum_i alpha_i*p_i + beta*m',
            'exact_generator_square':'(|p|^2+m^2)*I4',
            'n_tick_operator_error_bound':'t*a*sum_(r<s)|h_r*h_s|; h=(px,py,pz,m), t=n*a',
            'compact_momentum_bound':'t*a*(3*P^2+3*abs(m)*P) when abs(p_i)<=P',
            'one_time_unit_test_error_bounds':bounds,
            'finite_native_coupling':{
                'gate':'C=(I-P_occ) tensor I_native + P_occ tensor U_native',
                'Gauss_constraint_preserved':True,
                'occupation_preserved':True,
                'native_preparation':'positive Gauss orbit of (g=1,p=5,j=3,r=0,address=0)',
                'native_overlap_U':str(overlap),
                'coherent_occupation_reference_purity':str(purity),
                'routing_reference':'one fermion coherently in the tested or reference cell; fixed total parity',
                'nonproduct_joint_gate':True,
                'weak_tick_or_relativistic_limit_of_fixed_U_control':False},
            'input_mass_selected_by_repair':False,
            'native_U_invariant_subspace_implements_this_routing':False,
            'interacting_relativistic_scaling_limit':False,
            'physical_spin_statistics_species_or_GeV_mass':False}
