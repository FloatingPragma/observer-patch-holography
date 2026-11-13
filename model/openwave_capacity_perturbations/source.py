"""Exact sound-speed and finite perturbation controls for capacity histories."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import json
import sympy as s
from model.finite_area_response.response import Interval
from model.openwave_universal_energy_geometry.source import root_interval

ROOT=Path(__file__).resolve().parents[2]

def ca_polynomial(x,r):
    return (-2+(-3)*x+3*x*x+2*x*x*x+r*x*(1+x))/6

def causal_speed_interval(omega,r):
    x=root_interval(omega)
    return ca_polynomial(x,r)

def symbolic_identity():
    x,r=s.symbols('x r',positive=True)
    w=-(1+2*x)/3
    xp=x*(1-x*x)*(1+r+2*x)/2
    ca=s.factor(w-s.diff(w,x)*xp/(3*(1+w)))
    assert s.simplify(ca-ca_polynomial(x,r))==0
    assert s.simplify(s.diff(ca,x,2)-(1+2*x+r/3))==0
    return str(s.expand(ca))

def past_bound(parent):
    om=F(parent['inputs']['Omega_m0']);orr=F(parent['inputs']['Omega_r0']);x0=F(parent['inputs']['Omega_dark0'])
    one=next(row for row in parent['future_event_branch']if row['redshift']==1)
    one_hi=F(one['Omega_dark_interval'][1])
    assert x0<F(25,36) and one_hi<F(121,400)
    r1=2*orr/(om+2*orr);assert r1<F(1,1000)
    # c_a² increases with r and is convex in x. Therefore each rectangular
    # upper bound is the maximum at the two x endpoints and upper r.
    near=max(ca_polynomial(F(0),F(1,1000)),ca_polynomial(F(5,6),F(1,1000)))
    far=max(ca_polynomial(F(0),F(1)),ca_polynomial(F(11,20),F(1)))
    assert near<F(-1,5) and far<F(-1,4)
    return {'near_redshift_interval':'[0,1]','near_x_upper_exact':'5/6','near_r_upper_exact':'1/1000',
        'near_sound_speed_squared_upper_exact':str(near),'past_redshift_interval':'[1,infinity)',
        'past_x_upper_exact':'11/20','past_r_upper_exact':'1',
        'past_sound_speed_squared_upper_exact':str(far),'global_past_upper_exact':'-1/5',
        'used_parent_Omega_at_one_upper_exact':str(one_hi),
        'continuation_argument':'the positive conserved branch has dOmega/dlog(a)>0 and0<Omega<1; all earlier fractions are below the certified z=1 value'}

def periodic_sound_modes(n=8):
    L=s.Matrix(n,n,lambda i,j:2*int(i==j)-int(j==(i-1)%n)-int(j==(i+1)%n))
    eigenvalues=L.eigenvals();assert sum(eigenvalues.values())==n
    # The declared linear acoustic comparison has d²delta/dt²=-c_s² L delta.
    return [{'Laplacian_eigenvalue_exact':str(k),'multiplicity':int(v),
        'adiabatic_negative_speed_mode':'constant conserved mode'if k==0 else 'growing and decaying real exponents',
        'nonnegative_rest_speed_control':'constant or oscillatory; no real positive acoustic exponent'}
        for k,v in sorted(eigenvalues.items(),key=lambda kv:float(kv[0]))]

def canonical_completion():
    x,r=s.symbols('x r',positive=True)
    xp=x*(1-x*x)*(1+r+2*x)/2;rp=-r*(1-r)
    # In reduced Planck units, phi'=sqrt(2*x²*(1-x));
    # V=rho_x*(2+x)/3 and rho_x'=-2*(1-x)*rho_x.
    v1=-2*(1-x)*(2+x)+xp
    yprime_over_y=xp*(1/x-1/(2*(1-x)))
    numerator=s.diff(v1,x)*xp+s.diff(v1,r)*rp-2*(1-x)*v1-yprime_over_y*v1
    mass=s.Poly(s.cancel(numerator/(2*(1-x))),x,r)
    # Exact positive Bernstein basis on the entire unit square. Because
    # the basis sums to one, coefficient extrema bound the full polynomial.
    bernstein=[]
    for i in range(7):
        row=[]
        for j in range(3):
            coefficient=s.factor(sum(v*s.Rational(comb(i,k),comb(6,k))*s.Rational(comb(j,l),comb(2,l))
                for(k,l),v in mass.terms()if k<=i and l<=j))
            row.append(coefficient)
        bernstein.append(row)
    rebuilt=sum(bernstein[i][j]*s.binomial(6,i)*x**i*(1-x)**(6-i)*s.binomial(2,j)*r**j*(1-r)**(2-j)
                for i in range(7)for j in range(3))
    assert s.expand(rebuilt-mass.as_expr())==0
    assert min(v for row in bernstein for v in row)==2
    assert max(v for row in bernstein for v in row)==6
    # Positive scalar kinetic energy and these potential values reproduce
    # rho and p algebraically, hence conservation implies the scalar EOM.
    rho=s.symbols('rho',positive=True)
    kinetic=rho*(1-x)/3;potential=rho*(2+x)/3
    assert s.simplify(kinetic+potential-rho)==0
    assert s.simplify((kinetic-potential)/rho+(1+2*x)/3)==0
    return {'declared_action':'canonical minimally coupled scalar with the reconstructed potential V(phi)',
        'phi_derivative_per_log_a':'M_Pl*sqrt(2*Omega_x*(1-sqrt(Omega_x)))',
        'potential':'rho_x*(2+sqrt(Omega_x))/3',
        'kinetic_energy_density':'rho_x*(1-sqrt(Omega_x))/3',
        'field_origin':'phi(a=1)=0 is an irrelevant field-coordinate choice',
        'invertibility':'phi increases strictly with a on0<Omega_x<1, defining a single-valued parametric V(phi)',
        'Klein_Gordon_equation':'follows from exact rho/p reconstruction and separate scalar energy conservation where dot(phi)>0',
        'rest_frame_sound_speed_squared_exact':'1',
        'V_second_derivative_over_H_squared_polynomial':str(mass.as_expr()),
        'Bernstein_coefficients_exact':[[str(v)for v in row]for row in bernstein],
        'V_second_derivative_over_H_squared_bounds_exact':['2','6'],
        'positive_fixed_background_quadratic_gradient_and_potential':True,
        'radiation_dominated_potential_asymptotic':'V proportional to(phi-phi_early)^(-2)',
        'matter_dominated_small_Omega_potential_asymptotic':'V proportional to(phi-phi_matter_reference)^(-4) within the dust asymptotic regime',
        'ancestry':'constructed from the declared capacity history and inherited density amplitudes; not a source-derived new scalar particle',
        'source_scalar_transfer_proved':False,'full_gravity_perturbation_or_CMB_fit_proved':False}

def run():
    parent=json.loads((ROOT/'model/openwave_radiation_capacity/results.json').read_text())
    identity=symbolic_identity();bound=past_bound(parent);scalar=canonical_completion()
    om,orr=map(F,(parent['inputs']['Omega_m0'],parent['inputs']['Omega_r0']))
    rows=[]
    for r in parent['future_event_branch']:
        z=F(r['redshift']);ratio=orr*(1+z)/(om+orr*(1+z))
        speed=causal_speed_interval(Interval(*map(F,r['Omega_dark_interval'])),Interval.point(ratio))
        assert speed.hi<F(-1,5)
        rows.append({'redshift':r['redshift'],'adiabatic_sound_speed_squared_interval':speed.encoded(),
            'nonadiabatic_pressure_coefficient_for_zero_rest_speed_interval':(-speed).encoded()})
    acoustic=[]
    for k in (1,10,100):
        # Constant-coefficient diagnostic with declared damping2H. This is
        # not a fitted cosmological growth equation or a measured instability.
        rate=root_interval(Interval.point(1+F(k*k,5)))-1
        acoustic.append({'physical_wavenumber_over_H':k,'declared_damping_over_H':2,
            'growing_exponent_over_H_lower_interval':rate.encoded()})
    return {'schema':'oph.capacity_adiabatic_perturbation_boundary.v1','verdict':'PASS_CONSERVATION_AND_ADIABATIC_INSTABILITY_BOUNDARY',
      'inputs':parent['inputs'],'symbolic_adiabatic_sound_speed_squared':identity,
      'canonical_scalar_completion':scalar,
      'definition':'c_a^2=dot(p)/dot(rho)=w-(dw/dlog(a))/(3*(1+w)); rho separately conserved',
      'uniform_past_certificate':bound,'retained_background_samples':rows,
      'complete_finite_acoustic_control':{'periodic_cells':8,'all_modes':periodic_sound_modes(),
        'operator':'delta_tt=-c_s^2 L delta, L the complete periodic nearest-neighbour Laplacian',
        'zero_mode_retained':True,'positive_and_negative_growth_branches_retained':True},
      'frozen_damped_acoustic_diagnostics':acoustic,
      'capacity_controls':{'fixed_capacity':'rho and p are constant; dot(p)/dot(rho) is undefined, so no sound speed is inferred',
        'conserved_Hubble_area':'c_a^2=4r/[3(3+r)] lies in[0,1/3]; its distinct failed background comparison is retained',
        'comoving_area_constant_w':'c_a^2=-1/3 for the separately conserved barotropic control'},
      'required_completion':{'adiabatic_barotropic_attachment_rejected':True,
        'stable_rest_frame_closure':'delta_p=c_s^2 delta_rho with c_s^2>=0 needs nonadiabatic coefficient c_s^2-c_a^2>=1/5 throughout the past branch',
        'source_entropy_or_additional_field_required':'derive the stated nonadiabatic stress or another causal stable perturbation law from the same capacity source',
        'background_conservation_does_not_determine_rest_sound_speed':True},
      'scope':{'conditional_fluid_perturbation_result':True,'native_observer_perturbation_dynamics_derived':False,
        'every_OPH_capacity_model_excluded':False,'CMB_or_structure_likelihood_executed':False,
        'new_observational_match':False,'new_empirical_fit':False}}
