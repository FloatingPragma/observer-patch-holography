"""Exact rational bounds on the source weighted twist and its tails."""
from fractions import Fraction as F


def rational(x, name, *, positive=False):
    if type(x) not in (int, F) or x < 0 or (positive and x == 0):
        raise ValueError(name+' must be a nonnegative exact rational')
    return F(x)


def mesh(channel, spacing, radius, m0, m1, m2, m3=F(0)):
    """Bounds per unit of the actual ground correlation kappa.

    The test is C3 with support in the max-norm radius r box. mj bounds
    every derivative of order j (m0 bounds its value). The plateau of
    the source weight covers radius r+p. No kappa is fitted here.
    """
    if channel not in ('Q', 'P'):
        raise ValueError('one transverse Q or P kick required')
    p=rational(spacing,'spacing',positive=True)
    r=rational(radius,'radius',positive=True)
    m0,m1,m2,m3=[rational(x,'derivative bound')for x in (m0,m1,m2,m3)]
    if p>1:raise ValueError('require spacing at most one')
    v=(2*r+2*p)**3
    if channel=='Q':
        h=4*m0*m1+2*(r+p)*(m1*m1+m0*m2)
        quadrature=F(27,4)*p*p*h*v
        edge=F(3,2)*(r+p)*p*p*m0*m2*v
        cosine=2*(r+p)*v*p**4*m0**4
    else:
        h=4*m1*m2+2*(r+p)*(m2*m2+m1*m3)
        quadrature=F(27,16)*p*p*h*v
        edge=(r+p)*v*p*p*m1*m3/8
        cosine=(r+p)*v*p**4*m1**4/8
    return {'channel':channel,'spacing':p,'support_radius':r,
            'cutoff_plateau_radius_required':r+p,
            'midpoint_quadrature_error_per_kappa':quadrature,
            'edge_sampling_error_per_kappa':edge,
            'cosine_remainder_per_kappa':cosine,
            'total_error_per_kappa':quadrature+edge+cosine}


def schwartz_tail(channel, radius, power, amplitude, cutoff_gradient=F(2)):
    """Weighted continuum tail per kappa; pi is bounded above by22/7.

    Euclidean seminorms: |f| and (for P) |grad f| <= C(1+|X|)^(-N).
    The cutoff is one inside L, zero beyond2L and |grad zeta|<=b/L.
    """
    if channel not in ('Q','P'):raise ValueError('Q or P required')
    l=rational(radius,'radius',positive=True)
    c=rational(amplitude,'Schwartz seminorm')
    b=rational(cutoff_gradient,'cutoff derivative')
    if type(power)is not int or power<=2 or l<1:
        raise ValueError('integer power above two and radius at least one required')
    tail=l**(4-2*power)/F(2*power-4)
    value=(24 if channel=='Q'else 4*(1+b*b/l**2))*F(22,7)*c*c*tail
    return {'channel':channel,'radius':l,'power':power,'weighted_tail_per_kappa':value}


def mixed_kick(spacing, radius, f0, f1, f2, g0, g1, g2, g3):
    """Actual symmetric Q/P kick: signed polarization and a source error.

    The nonlinear rotation remainder uses the actual bond-covariance norm,
    independently of kappa. The two mesh errors remain per kappa.
    """
    q=mesh('Q',spacing,radius,f0,f1,f2)
    p=mesh('P',spacing,radius,g0,g1,g2,g3)
    spacing=q['spacing'];radius=q['support_radius']
    f0=rational(f0,'F amplitude');g0=rational(g0,'G amplitude')
    g1=rational(g1,'G derivative')
    t=spacing*g0+spacing**2*f0
    if t>1:raise ValueError('source small-angle bound requires t at most one')
    remainder=9*t*(radius+spacing)*(2*radius+2*spacing)**3*(g1*g1+4*f0*f0)
    return {'source_small_angle_upper':t,
            'nonlinear_rotation_absolute_error':remainder,
            'separate_mesh_errors_per_kappa':q['total_error_per_kappa']+p['total_error_per_kappa'],
            'mixed_quadratic_term_vanishes_by_actual_bond_covariance':True}


def charged_truncation(radius, power, amplitude, field_scale_squared=F(1)):
    """Square of a complete three-record preparation isometry error.

    Q3/Q2 have the same source variance. This is a state-specific bound
    on the actual vacuum preparation, never a universal operator norm.
    """
    l=rational(radius,'radius',positive=True)
    c=rational(amplitude,'Schwartz seminorm')
    scale=rational(field_scale_squared,'field amplitude squared')
    if type(power)is not int or power<=3 or l<1:
        raise ValueError('integer power above three and radius at least one required')
    l1=4*F(22,7)*c*l**(3-power)/F(power-3)
    l2_squared=4*F(22,7)*c*c*l**(3-2*power)/F(2*power-3)
    # Split Fourier momentum at1. 1/(4*pi²) <=1/36 since pi>3.
    variance=3*(l1*l1/36+l2_squared)
    return {'tail_L1_upper':l1,'tail_L2_squared_upper':l2_squared,
            'actual_configuration_variance_upper':variance,
            'complete_isometry_error_squared_upper':4*scale*variance,
            'complement_retained':True,'operator_norm_error_claimed':False}


def scalar_field_remainder(spacing, support_radius, first_derivative):
    p=rational(spacing,'spacing',positive=True)
    r=rational(support_radius,'radius',positive=True)
    m=rational(first_derivative,'first derivative')
    return {'Q_current_remainder_operator_norm_upper':m*(2*r+2*p)**3/4,
            'this_bound_vanishes_as_spacing_tends_to_zero':False,
            'source_four_spin_correlation_control_still_required':True}


def finite_readout(coefficients, scalar, reference_tv, created_tv):
    if type(scalar)not in (int,F)or any(type(a)not in (int,F)for a in coefficients):
        raise ValueError('exact signed source coefficients required')
    e0=rational(reference_tv,'reference TV error')
    e1=rational(created_tv,'created TV error')
    if e0>1 or e1>1:raise ValueError('TV error must not exceed one')
    coefficients=tuple(F(a)for a in coefficients if a!=0)
    m=len(coefficients)
    if not m:
        return {'nonconstant_terms':0,'assign_constant_to_every_record':True,
                'exact_increment':F(0),'increment_error_upper':F(0)}
    maximum=max(abs(F(scalar)+m*a*sign)for a in coefficients for sign in (-1,1))
    return {'nonconstant_terms':m,'assign_constant_to_every_record':False,
            'all_failure_records_retained_with_zero_score':True,
            'absolute_score_upper':maximum,'increment_error_upper':2*maximum*(e0+e1)}
