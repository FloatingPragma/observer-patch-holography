"""Exact variational energy and current-coefficient bounds for this source."""
from fractions import Fraction as F
from itertools import product
import sympy as s

from model.spin_goldstone_sector import infrared
from model.centered_particle_spectrum.certificate import sqrt_bounds, decimal_bounds


def pair_certificate():
    x=s.Matrix([[0,1],[1,0]])/2
    y=s.Matrix([[0,-s.I],[s.I,0]])/2
    z=s.diag(1,-1)/2
    identity=s.eye(2)
    a=(s.sqrt(3)+s.sqrt(2))/s.sqrt(10)
    b=(s.sqrt(3)-s.sqrt(2))/s.sqrt(10)
    vector=s.Matrix([0,a,-b,0])
    h=sum((s.kronecker_product(axis,axis) for axis in (x,y,z)),s.zeros(4))
    norm=s.simplify((vector.H*vector)[0])
    energy=s.simplify((vector.H*h*vector)[0])
    left=s.simplify((vector.H*s.kronecker_product(z,identity)*vector)[0])
    right=s.simplify((vector.H*s.kronecker_product(identity,z)*vector)[0])
    if norm!=1 or energy!=-s.Rational(7,20) or left!=-right or left**2!=s.Rational(6,25):
        raise ArithmeticError('actual spin-one-half dimer expectation differs')
    angle=s.symbols('u', real=True)
    density=-s.Rational(3,4)-angle/4+s.Rational(5,8)*angle**2
    if s.expand(density+s.Rational(31,40)-s.Rational(5,8)*(angle-s.Rational(1,5))**2)!=0:
        raise ArithmeticError('dimer product energy minimizer differs')
    return {'pair_state':'[(sqrt(3)+sqrt(2))|up down>-(sqrt(3)-sqrt(2))|down up>]/sqrt(10)',
        'state_norm_squared':'1', 'twice_amplitude_product':'1/5',
        'intra_dimer_spin_dot_expectation':'-7/20',
        'inter_dimer_spin_dot_expectation':'-6/25',
        'staggered_single_site_spin_mean_squared':'6/25',
        'trial_spin_dot_energy_per_site':'-31/40',
        'trial_shifted_H_over_J_per_site':'-61/40'}


def torus_certificate(side=4):
    if type(side) is not int or side<4 or side%2:
        raise ValueError('even cubic torus side at least four required')
    vertices=tuple(product(range(side),repeat=3))
    edges=set()
    pairs=set()
    for v in vertices:
        if v[0]%2==0:
            pairs.add(tuple(sorted((v,((v[0]+1)%side,v[1],v[2])))))
        for axis in range(3):
            w=list(v);w[axis]=(w[axis]+1)%side
            edges.add(tuple(sorted((v,tuple(w)))))
    n=len(vertices)
    if len(edges)!=3*n or len(pairs)!=n//2 or not pairs<=edges:
        raise ArithmeticError('finite periodic source edge/dimer count differs')
    if any((sum(v)-sum(w))%2!=1 for v,w in edges):
        raise ArithmeticError('source torus has a non-bipartite edge')
    if set(sum((list(pair) for pair in pairs),[]))!=set(vertices):
        raise ArithmeticError('dimer covering misses a source site')
    intra, inter=F(-7,20), F(-6,25)
    density=(len(pairs)*intra+(len(edges)-len(pairs))*inter)/n
    if density!=F(-31,40):
        raise ArithmeticError('full periodic dimer trial energy differs')
    return {'side':side,'sites':n,'distinct_source_bonds':len(edges),
        'dimer_bonds':len(pairs),'inter_dimer_bonds':len(edges)-len(pairs),
        'periodic_wrap_bonds_retained':True,'spin_dot_energy_density':str(density)}


def certificate():
    pair=pair_certificate()
    ir=infrared.certificate(256)
    integral=F(ir['integral_upper'])
    e=F(31,40)
    low, high=sqrt_bounds(e/18)
    if not low*low<=e/18<=high*high:
        raise ArithmeticError('outward variational infrared square root required')
    bragg=e/9-high*integral
    bias=F(9)
    magnetization=2*bias*bragg/(1+bias*bias*bragg)
    if magnetization<=F(11,100):
        raise ArithmeticError('variational infrared input does not support the stated order')
    return {'exact_pair_source_certificate':pair,
        'periodic_source_count_controls':[torus_certificate(n) for n in (4,6)],
        'unshifted_ground_spin_energy_e_lower':str(e),
        'source_current_coefficient_kappa_lower':'1/120',
        'complete_inverse_current_measure_over_J_lower':'1/240',
        'source_current_coefficient_is_a_pole_speed':False,
        'unchanged_256_mesh_infrared_certificate':ir,
        'sqrt_e_over_18_outward_interval':[str(low),str(high)],
        'staggered_Bragg_weight_lower':str(bragg),
        'fixed_trial_bias':'9',
        'selected_magnetization_lower':str(magnetization),
        'selected_magnetization_decimal_lower':decimal_bounds(magnetization,magnetization,12)[0],
        'simple_strict_magnetization_lower':'11/100',
        'conservative_default_readout_order_lower':'1/10',
        'periodic_ground_projector_symmetry_changed':False,
        'empirical_ground_energy_or_magnetization_used':False}
