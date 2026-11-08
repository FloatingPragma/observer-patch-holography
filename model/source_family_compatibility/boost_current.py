"""Actual full-source energy-dipole current and its composite contribution."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from model.field_symmetry import symmetry
from model.full_fock_field.field import FieldGeometry
from .local_matter import annihilator


def exact(matrix,scale=1):
    return (s.Matrix(matrix.real.toarray().tolist())
            +s.I*s.Matrix(matrix.imag.toarray().tolist()))/scale


def simplify(matrix):return matrix.applyfunc(s.expand)


def lift(one,annihilators):
    result=s.zeros(annihilators[0].rows)
    for i in range(one.rows):
        for j in range(one.cols):
            if one[i,j]:result+=one[i,j]*annihilators[i].H*annihilators[j]
    return result


@lru_cache(maxsize=6)
def source_current(axis,contact=F(1)):
    if type(axis)is not int or axis not in range(3):raise ValueError('source axis0..2 required')
    if type(contact)not in (int,F)or contact<0:raise ValueError('exact nonnegative contact required')
    u=s.Rational(contact);c=s.Rational(1,8)
    # The source admits exactly represented binary rational coefficients.
    if F(float(contact))!=F(contact):raise ValueError('literal source coefficient must be binary rational')
    sides=tuple(3 if j==axis else 1 for j in range(3))
    source=symmetry.operators(FieldGeometry(sides,float(c),float(contact)))
    clock=exact(source.K,source.scale);dim=64;eye=s.eye(dim)
    a=tuple(annihilator(3,j)for j in range(6))
    n=tuple(x.H*x for x in a)
    b=tuple(2*(n[j^1]-eye/2)*a[j]for j in range(6))
    indices=[1<<j for j in range(6)]
    one=clock.extract(indices,indices)+u*s.eye(6)/2
    kinetic=lift(one,a)
    onsite=tuple(u*(n[2*x]*n[2*x+1]-(n[2*x]+n[2*x+1])/2)for x in range(3))
    if simplify(kinetic+sum(onsite,s.zeros(dim))-clock)!=s.zeros(dim):
        raise ArithmeticError('actual source onsite/kinetic split differs')
    X=s.diag(*(x for x in range(3)for _ in range(2)))
    dipole=lift((X*one+one*X)/2,a)+sum((x*onsite[x]for x in range(3)),s.zeros(dim))
    current=simplify(s.I*(clock*dipole-dipole*clock))
    free=simplify(lift(s.I*(one**2*X-X*one**2)/2,a))
    composite=s.zeros(dim)
    for x in range(2):
        y=x+1
        for sigma in range(2):
            for tau in range(2):
                i,j=2*x+sigma,2*y+tau
                coefficient=one[i,j]
                if coefficient:
                    forward=coefficient*(b[i].H*a[j]+a[i].H*b[j])
                    composite+=s.I*u*(y-x)*(forward-forward.H)/4
    composite=simplify(composite)
    if simplify(current-free-composite)!=s.zeros(dim):
        raise ArithmeticError('source dipole Ward identity has incorrect interaction coefficient')
    if current.H!=current or free.H!=free or composite.H!=composite:
        raise ArithmeticError('energy current must be Hermitian')
    # Number conservation and zero-contact limit concern every input column.
    number=sum(n,s.zeros(dim))
    if simplify(current*number-number*current)!=s.zeros(dim):
        raise ArithmeticError('source current changes occupation number')
    square_trace=s.simplify(s.trace(composite.H*composite))
    if bool(u) != bool(square_trace):raise ArithmeticError('contact deletion control lost sensitivity')
    # N=1 gives the exact centering contribution; it is not an interacting vacuum.
    expected_one=s.I*((one-u*s.eye(6)/2)**2*X-X*(one-u*s.eye(6)/2)**2)/2
    if simplify(current.extract(indices,indices)-expected_one)!=s.zeros(6):
        raise ArithmeticError('N=1 energy-reference control differs')
    return {'axis':axis,'geometry':list(sides),'full_Fock_dimension':dim,
            'source_c':str(c),'source_U':str(u),'source_coordinates':'integer site coordinates',
            'full_source_energy_dipole_identity':True,
            'composite_current_Hilbert_Schmidt_norm_squared':str(square_trace),
            'interaction_and_centering_remainder_nonzero':bool(square_trace),
            'all_occupation_sectors_checked':True,
            'number_conservation_checked':True,
            'one_particle_centering_control_checked':True,
            'energy_current_identified_with_geometric_momentum':False,
            'physical_boost_or_Lorentz_covariance_derived':False}


def alias_bound(spacing,momentum_radius,particle_number=1):
    if type(spacing)not in(int,F)or spacing<=0 or type(momentum_radius)not in(int,F)or momentum_radius<0:
        raise ValueError('exact positive spacing and nonnegative momentum radius required')
    if type(particle_number)is not int or particle_number<0:raise ValueError('finite particle number required')
    a,p=F(spacing),F(momentum_radius)
    return {'spacing':a,'momentum_radius':p,'particle_number':particle_number,
            'free_current_over_speed_squared_momentum_error':F(2,3)*particle_number*a*a*p**3,
            'same_window_source_translation_required':True,
            'interacting_history_momentum_error_bounded':False}


def certificate():
    return {'actual_Weyl_axes':[source_current(j)for j in range(3)],
            'zero_contact_controls':[source_current(j,F(0))for j in range(3)],
            'alias_example':{k:str(v)if isinstance(v,F)else v for k,v in alias_bound(F(1,16),2,3).items()},
            'source_composite_current_is_retained':True,
            'native_field_boost_or_Einstein_premise_completed':False}
