"""Exact CAR witnesses for a non-Gaussian multi-time ground field.

The finite matrices verify universal local identities.  The infinite-ground
claim uses the selected state and norm-differentiable original dynamics, as
proved in README.md; these matrices do not approximate that state.
"""
from collections import Counter
from fractions import Fraction as F

from ..field_symmetry import symmetry as s
from .algebra import anti


WORDS=(('b','ad','a','bd'),('b','ad','bd','a'),
       ('ad','b','a','bd'),('ad','b','bd','a'))


def _monomial(left,right):
    # Scalar two-point contractions commute; their operator arguments do not.
    return tuple(sorted((tuple(left),tuple(right))))


def wick_pairing_residual(pairing_signs=(1,-1,1)):
    """All Wick terms minus <{b,a*}><{a,b*}>, with anomalous terms retained."""
    value=Counter()
    for a,b,c,d in WORDS:
        for coefficient,left,right in ((pairing_signs[0],(a,b),(c,d)),
                                        (pairing_signs[1],(a,c),(b,d)),
                                        (pairing_signs[2],(a,d),(b,c))):
            value[_monomial(left,right)]+=coefficient
    for left in (('b','ad'),('ad','b')):
        for right in (('a','bd'),('bd','a')):
            value[_monomial(left,right)]-=1
    return {key:coefficient for key,coefficient in value.items() if coefficient}


def local_derivative_certificate():
    operators=[s._annihilator(2,j) for j in range(2)]
    numbers=[a.adjoint()@a for a in operators]
    contact=numbers[0]@numbers[1]
    for spin,a in enumerate(operators):
        derivative=s._scale(s.commutator(contact,a),1j)
        expected=s._scale(numbers[1-spin]@a,-1j)
        if not (derivative-expected).is_zero():
            raise ArithmeticError('contact Heisenberg derivative sign changed')
        if not (anti(derivative,a.adjoint())-
                s._scale(numbers[1-spin],-1j)).is_zero():
            raise ArithmeticError('centered anticommutator derivative changed')
        if not (numbers[spin]@numbers[spin]-numbers[spin]).is_zero():
            raise ArithmeticError('onsite occupation is not a projection')
    if wick_pairing_residual():
        raise ArithmeticError('full fermionic Wick cancellation failed')
    return {'onsite_Fock_dimension':4,'spins_checked':2,
            'contact_derivative':'delta(c_sigma)=-i U n_opposite c_sigma',
            'anticommutator_derivative':'{delta(c_sigma),c_sigma^dagger}=-i U n_opposite (plus scalar kinetic diagonal)',
            'centered_derivative':'D_prime(0)=-i U (n_opposite-r I)',
            'projection_variance':'omega((n_opposite-r I)^2)=r(1-r)',
            'four_point_words':[list(word) for word in WORDS],
            'cumulant_convention':'kappa(A,B,C,D)=omega(ABCD)-omega(AB)omega(CD)+omega(AC)omega(BD)-omega(AD)omega(BC)',
            'identity':'omega(D(t)^dagger D(t))=sum of the four displayed kappa words',
            'anomalous_pairings_retained':True,
            'formal_Wick_residual_terms':0}


def constants(speed=F(1,8),interaction=F(1)):
    if (type(speed) not in (int,F) or type(interaction) not in (int,F) or
            speed<0 or interaction<=0):
        raise ValueError('exact c>=0 and U>0 required')
    c,u=F(speed),F(interaction)
    m=9*c*c+12*c*u+u*u
    return {'kinetic_row_l1':3*c,'first_derivative_norm':3*c+u,
            'occupation_derivative_norm':6*c,'second_derivative_norm':m,
            'time_radius':u/(16*m),'variance_t_squared_lower':u*u/64}


def variance_lower(time,speed=F(1,8),interaction=F(1)):
    if type(time) not in (int,F):raise ValueError('exact nonzero time required')
    value=constants(speed,interaction);t=abs(F(time))
    if not 0<t<=value['time_radius']:
        raise ValueError('time outside the proved nonzero witness interval')
    return value['variance_t_squared_lower']*t*t


def quadratic_negative_control():
    """Exact paired/auxiliary example, distinct from the general proof."""
    a=[s._annihilator(4,j) for j in range(4)]
    term=a[0].adjoint()@a[2].adjoint()+a[0].adjoint()@a[3]
    h=term+term.adjoint()+a[2].adjoint()@a[2]
    current=a[0];identity=s._identity(16)
    for _ in range(5):
        value=anti(current,a[0].adjoint())
        coefficient=complex(value.matrix()[0,0])
        if not (value-s._scale(identity,coefficient)).is_zero():
            raise ArithmeticError('quadratic paired control lost scalar anticommutator')
        current=s._scale(s.commutator(h,current),1j)
    return {'Fock_dimension':16,'paired_and_auxiliary_modes_present':True,
            'derivative_orders_exactly_checked':[0,1,2,3,4],
            'scalar_anticommutators':True,
            'state_needed_for_scalar_operator_identity':False,
            'finite_example_is_general_proof':False}


def certificate():
    values=constants()
    radius=values['time_radius']
    return {'local_identity':local_derivative_certificate(),
            'quadratic_negative_control':quadratic_negative_control(),
            'source_parameters':{'c':'1/8','U':'1'},
            'uniform_bounds':{key:str(value) for key,value in values.items()},
            'bound_at_time_radius':str(variance_lower(radius)),
            'short_time_limit':'omega(D(t)^dagger D(t))/t^2 -> U^2 r(1-r) >=7 U^2/64',
            'nonzero_four_point_witness':'For every 0<|t|<=U/(16M), sum kappa=omega(D^dagger D)>=U^2 t^2/64>0; M=9c^2+12cU+U^2.',
            'state_domain':'the selected gauge-invariant hidden-SU2 thermodynamic ground, with 1/8<=r<=1/2',
            'non_Gaussian_multitime_ground_field_proved':True,
            'linear_CAR_quadratic_evolution_excluded':True,
            'quadratic_exclusion_includes_pairing_and_auxiliaries':True,
            'quadratic_exclusion_requires_Gaussian_reference':False,
            'nonlinear_observable_encodings_excluded':False,
            'local_state_faithfulness_assumed':False,
            'continuum_limit_or_scattering_proved':False,
            'physical_particle_species_identified':False}
