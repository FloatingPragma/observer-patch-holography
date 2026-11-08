"""Exact auxiliary-phase bounds for the joint-family second canonical pair."""
from fractions import Fraction as F
from math import isqrt
import sympy as s
from model.native_ordered_ground import source


def rational(value,positive=False):
    if type(value)not in(int,F)or value<0 or(positive and value==0):raise ValueError('nonnegative rational required')
    return F(value)


def inputs(n,n0,*norms):
    if type(n)is not int or n<2 or type(n0)is not int or n0<1:raise ValueError('family index at least two and positive support count required')
    return (F(1,n*n),*(rational(x)for x in norms))


def sqrt_upper(x,bits=40):
    x=rational(x)
    if type(bits)is not int or bits<0:raise ValueError('nonnegative precision required')
    denominator=1<<bits
    k=isqrt((x.numerator*denominator*denominator)//x.denominator)
    if k*k*x.denominator<x.numerator*denominator*denominator:k+=1
    return F(k,denominator)


def auxiliary(n,n0,f0,f1):
    p,f0,f1=inputs(n,n0,f0,f1)
    return {'p':p,'field_over_J':F(1,n**5),
        'operator_norm_upper':n0*f0/2,
        'mean':F(),
        'vacuum_norm_squared_upper':12*n0*f1*f1/F(n**5)+2*n0*f0*f0/F(n**6),
        'scaled_kick_kappa_squared':F(3,2)*n0*f1*f1/F(n**8)+n0*f0*f0/F(4*n**9),
        'Q_commutator_per_unit_test_amplitude':n0*f0*p*p/2,
        'P_commutator_per_unit_test_amplitude':n0*f0*p/2,
        'C_commutator_per_unit_test_amplitude':n0*f0*p*p/2,
        'A_commutator_per_unit_test_amplitude':n0*f0*p**3/2}


def longitudinal(n,n0,h0,h1,h2):
    p,h0,h1,h2=inputs(n,n0,h0,h1,h2)
    variance=12*n0*h1*h1/F(n);mean=3*n0*h2*p/8
    return {'p':p,'operator_norm_upper':n0*h0/(2*p),
        'variance_upper':variance,'absolute_mean_upper':mean,
        'vacuum_norm_squared_upper':variance+mean*mean,
        'scaled_kick_kappa_squared':F(3,2)*n0*h1*h1*p*p,
        'Q_commutator_per_unit_test_amplitude':n0*h0*p/2,
        'unbounded_C_strong_history_convergence_claimed':False}


def time_phase_squared(kappa_squared,absolute_time,phase_amplitude,history_energy_root):
    k2,t,a,e=(rational(x)for x in(kappa_squared,absolute_time,phase_amplitude,history_energy_root))
    return 4*t*a*sqrt_upper(k2)*e+10*t*a*a*k2


def C_P_rotation_remainder(n,n0,h0,g0,phase_C=F(1),phase_P=F(1)):
    p,h0,g0,a,b=inputs(n,n0,h0,g0,phase_C,phase_P)
    return F(n0,4)*(a*a*b*h0*h0*g0*p*p+a*b*b*h0*g0*g0*p)


def P_P_rotation_remainder(n,n0,g0,h0,phase_P2=F(1),phase_P3=F(1)):
    p,g0,h0,a,b=inputs(n,n0,g0,h0,phase_P2,phase_P3)
    amplitude=a*g0+b*h0
    if p*amplitude>1:raise ValueError('fourth-order onsite angle bound requires p*(a*g0+b*h0)<=1')
    return 2*n0*amplitude**4*p


def algebra_checks():
    i=s.I;half=s.Rational(1,2)
    sx=s.Matrix([[0,1],[1,0]])/2
    sy=s.Matrix([[0,-i],[i,0]])/2
    sz=s.diag(1,-1)/2
    def comm(a,b):return a*b-b*a
    assert comm(-sy,sz)==-i*sx
    assert comm(sz,-sy)==i*sx
    assert comm(sy,sz)==i*sx
    a,b,u=s.symbols('a b u',real=True)
    x=i*a*sy;y=i*b*sz;c=comm(x,y);d=comm(x+y,c)/2
    assert c==-i*a*b*sx
    assert d==-i*a*a*b*sz/2+i*a*b*b*sy/2
    def trunc(m):return m.applyfunc(lambda e:sum(s.expand(e).coeff(u,k)*u**k for k in range(4)))
    def exponential(m):return sum((m**k/s.factorial(k)for k in range(4)),s.zeros(2))
    f=s.eye(2)
    for generator in(x,y,-x,-y):f=trunc(f*exponential(u*generator))
    g=trunc(exponential(u*u*c)*exponential(u**3*d))
    assert s.simplify(f-g)==s.zeros(2)
    zthird=comm(x,c)/6+comm(y,c)/3
    split=trunc(trunc(trunc(exponential(u*x)*exponential(u*y))*exponential(-u*u*c/2))*exponential(u**3*zthird))
    assert s.simplify(trunc(exponential(u*(x+y)))-split)==s.zeros(2)
    q=s.symbols('q',nonnegative=True);cn=q*q/8;dn=q**3/16
    aa=[1,2*cn,2*cn+4*cn**2,12*cn**2+8*cn**3,12*cn**2+48*cn**3+16*cn**4]
    bb=[1,3*dn,6*dn+9*dn**2,6*dn+54*dn**2+27*dn**3,180*dn**2+324*dn**3+81*dn**4]
    fourth=s.Poly(s.expand(sum(s.binomial(4,k)*aa[k]*bb[4-k]for k in range(5))),q)
    assert all(power[0]>=4 and coefficient>=0 for power,coefficient in fourth.terms())
    total=sum(fourth.all_coeffs());assert total<4
    cz=q*q/16;dz=q**3/24
    az=[1,2*cz,2*cz+4*cz**2,12*cz**2+8*cz**3,12*cz**2+48*cz**3+16*cz**4]
    bz=[1,3*dz,6*dz+9*dz**2,6*dz+54*dz**2+27*dz**3,180*dz**2+324*dz**3+81*dz**4]
    split_fourth=s.Poly(s.expand(sum(s.factorial(4)/s.factorial(k)/s.factorial(l)/s.factorial(4-k-l)*(q/2)**k*az[l]*bz[4-k-l]for k in range(5)for l in range(5-k))),q)
    assert all(power[0]>=4 and coefficient>=0 for power,coefficient in split_fourth.terms())
    split_total=sum(split_fourth.all_coeffs());assert split_total<100
    # Exact physical-frame source double commutators on one bond.
    spins=[sx,sy,sz];site=lambda op,which:s.kronecker_product(op,s.eye(2))if which==0 else s.kronecker_product(s.eye(2),op)
    j,field,r,t=s.symbols('J B r t',real=True)
    exchange=j*(sum((site(v,0)*site(v,1)for v in spins),s.zeros(4))-s.eye(4)/4)
    selection=-field*(site(sx,0)-site(sx,1))
    native,scale=source.matrix(2,((0,1),),F(1),F(1,4))
    marshall=s.kronecker_product(s.eye(2),2*sz)
    assert marshall*s.Matrix(native.tolist())*marshall/scale==(exchange+selection).subs({j:1,field:s.Rational(1,4)})
    for axis in range(3):
        operator=r*site(spins[axis],0)+t*site(spins[axis],1)
        transverse=sum((site(spins[k],0)*site(spins[k],1)for k in range(3)if k!=axis),s.zeros(4))
        expected=-j*(r-t)**2*transverse
        if axis!=0:expected+=field*(r*r*site(sx,0)-t*t*site(sx,1))
        assert s.simplify(comm(operator,comm(exchange+selection,operator))-expected)==s.zeros(4)
    return {'same_sign_canonical_pair_commutators_checked':True,
        'mixed_P_commutator_minus_i_C_checked':True,
        'cubic_group_commutator_checked':True,
        'fourth_derivative_coefficient_sum':str(total),
        'fourth_derivative_bound_below_four':True,
        'rotated_source_momentum_cubic_splitting_checked':True,
        'rotated_momentum_fourth_derivative_sum':str(split_total),
        'rotated_momentum_fourth_derivative_bound_below_hundred':True,
        'all_three_source_double_commutators_checked':True,
        'actual_native_two_site_source_attached':True}


def certificate():
    rows=[]
    for n in(2,4,16,64,256):
        a=auxiliary(n,27,F(1),F(1));c=longitudinal(n,27,F(1),F(1),F(1))
        rows.append({'n':n,'auxiliary':a,'longitudinal':c,
            'auxiliary_time_phase_squared_upper':time_phase_squared(a['scaled_kick_kappa_squared'],F(1),F(1),F(2)),
            'longitudinal_time_phase_squared_upper':time_phase_squared(c['scaled_kick_kappa_squared'],F(1),F(1),F(2)),
            'C_P_group_commutator_remainder_upper':C_P_rotation_remainder(n,27,F(1),F(1)),
            'P_P_group_commutator_remainder_upper':P_P_rotation_remainder(n,27,F(1),F(1))})
    return {'family':'p=n^-2, B=J*n^-5 with actual H_B ground and clock',
        'proof_order':['auxiliary_transverse_phases','longitudinal_phases','mixed_momentum_Weyl_relation','regularity_and_internal_U1'],
        'source_algebra':algebra_checks(),'bounds':rows,
        'four_field_history_energy_bound_consumed_from_parent':True,
        'parent_vacuum_moment_bounds_used_only_for_regularity':True,
        'P3_regularity_assumed_before_CCR':False,
        'common_spatial_action_uses_even_source_translations':True,
        'joint_history_state_identified_with_original_ordered_limit':False,
        'physical_Lorentz_spin_or_mass_assigned':False}
