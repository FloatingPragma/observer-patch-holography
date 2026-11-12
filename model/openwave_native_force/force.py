"""Source derivatives, finite geometry energies and a universality obstruction."""
from fractions import Fraction as F
from math import isqrt
import numpy as np
from scipy.fft import dstn, idstn
import sympy as s

from model.native_equilibration import vacuum_repair as native
from model.finite_area_response import response as interval


def source_matrices():
    """Dimensionless M=H/sqrt(14) of the actual three-oscillator cutoff two."""
    r=s.symbols('r', positive=True)
    terms=tuple(s.Matrix(x.astype(int)) for x in native._terms(2))
    actual=sum((coefficient*term for coefficient,term in
                zip(native._coefficient_expressions(2,r),terms,strict=True)),s.zeros(8))/s.sqrt(14)
    actual=actual.applyfunc(s.simplify)
    d=s.diag(*(2*mask.bit_count()+3 for mask in range(8)))
    x=s.Matrix([[0,1],[1,0]]); eye=s.eye(2)
    c=s.kronecker_product(x,x,eye)+s.kronecker_product(eye,x,x)
    expected=(r/4+1/(7*r)+3*r*r/28)*d-c/(7*r)
    assert (actual-expected).applyfunc(s.simplify)==s.zeros(8)
    m=actual.subs(r,1)
    j=actual.diff(r).subs(r,1)
    assert m==d/2-c/7 and j==9*d/28+c/7
    return r,actual,m,j


def expectation(matrix, vector):
    vector=s.Matrix(vector)
    return s.cancel((vector.T*matrix*vector)[0]/(vector.T*vector)[0])


def positive_energy_nonuniversality():
    _,_,m,j=source_matrices()
    reference=s.eye(8)[:,0]
    states={'one_occupation':s.eye(8)[:,1],
            'coherent_plus':s.eye(8)[:,0]+s.eye(8)[:,3],
            'coherent_minus':s.eye(8)[:,0]-s.eye(8)[:,3],
            'small_negative_coherence':s.eye(8)[:,0]-s.eye(8)[:,3]/10}
    rows=[]
    for name,vector in states.items():
        energy=expectation(m,vector)-expectation(m,reference)
        charge=expectation(j,vector)-expectation(j,reference)
        assert energy>0
        rows.append({'state':name,'excess_energy_over_sqrt14':str(energy),
                     'geometry_charge_over_sqrt14':str(charge),
                     'charge_to_energy':str(charge/energy)})
    assert len({row['charge_to_energy'] for row in rows})==4
    assert F(rows[-1]['geometry_charge_over_sqrt14'])<0
    return rows


def square_root_upper(value):
    value=F(value)
    if value<0: raise ValueError('positive root bound required')
    scale=10**24
    n=isqrt(value.numerator*scale*scale//value.denominator)
    upper=F(n+1,scale)
    assert upper*upper>=value
    return upper


def ground_counterexample():
    """Certify opposite response charges above the actual source ground."""
    _,_,m,j=source_matrices()
    spectrum=interval.spectral_enclosures(m)
    values,vectors=np.linalg.eigh(np.asarray(m,float))
    v=vectors[:,0]
    derivative=np.asarray(j,float)
    u=derivative@v-(v@derivative@v)*v
    assert np.linalg.norm(u)>.01
    u/=np.linalg.norm(u)
    rational=lambda vector:s.Matrix([s.Rational(float(x)) for x in vector])
    rv=rational(v)
    energy=F(expectation(m,rv)); jv=F(expectation(j,rv))
    gap=spectrum[1].lo-spectrum[0].hi
    fidelity_error=(energy-spectrum[0].lo)/gap
    assert 0<=fidelity_error<F(1,10**18)
    j_error=2*F(89,28)*square_root_upper(fidelity_error)
    rows=[]
    for sign in (-1,1):
        candidate=rational(v+sign*.01*u)
        e=F(expectation(m,candidate)); q=F(expectation(j,candidate))-jv
        emin,emax=e-spectrum[0].hi,e-spectrum[0].lo
        qmin,qmax=q-j_error,q+j_error
        assert emin>0
        assert qmax<0 if sign<0 else qmin>0
        rows.append({'sign':sign,'binary_rational_state':[str(x) for x in candidate],
                     'excess_energy_interval':[str(emin),str(emax)],
                     'geometry_charge_interval':[str(qmin),str(qmax)]})
    return {'ground_energy_interval':[str(spectrum[0].lo),str(spectrum[0].hi)],
            'next_energy_lower':str(spectrum[1].lo),
            'ground_reference_infidelity_upper':str(fidelity_error),
            'ground_geometry_charge_error_upper':str(j_error),
            'opposite_charge_positive_energy_states':rows,
            'universality_for_all_positive_energy_preparations':False,
            'ground_eigensolver_used_as_uncertified_proof':False}


def laplacian(field):
    padded=np.pad(field,1)
    return (6*field-padded[:-2,1:-1,1:-1]-padded[2:,1:-1,1:-1]
            -padded[1:-1,:-2,1:-1]-padded[1:-1,2:,1:-1]
            -padded[1:-1,1:-1,:-2]-padded[1:-1,1:-1,2:])


def green_solve(charge, stiffness=1., mass_squared=0.):
    """Stationarity of declared positive geometry energy, not input Poisson data."""
    n=charge.shape[0]
    if charge.shape!=(n,n,n) or stiffness<=0 or mass_squared<0:
        raise ValueError('positive cubic elastic action required')
    axis=2-2*np.cos(np.pi*np.arange(1,n+1)/(n+1))
    symbol=stiffness*(axis[:,None,None]+axis[None,:,None]+axis[None,None,:])+mass_squared
    return -idstn(dstn(charge,type=1,norm='ortho')/symbol,type=1,norm='ortho')


def finite_energy_certificate(n, sites, charge=F(9,140000), bits=48,
                              geometry_window=F(1,50000)):
    """Exact dyadic witness and lower bound for a fully finite geometry carrier.

    Each geometry coordinate lies on a finite declared dyadic mesh.
    Residuals and energies are evaluated with exact integers. The lower
    bound follows from lambda_min(L)>=12/(n+1)^2. No dense state-space
    enumeration or claim of a native execution of the minimizer is made.
    """
    charge=F(charge); geometry_window=F(geometry_window); denominator=2**bits
    if type(n) is not int or n<3 or n%2!=1 or type(bits) is not int or not 1<=bits<=48:
        raise ValueError('odd cubic grid and supported dyadic precision required')
    if geometry_window<=0 or geometry_window>F(1,10):
        raise ValueError('geometry window must lie within the proved source domain')
    if 12*int(geometry_window*denominator+1)*charge.denominator+abs(charge.numerator)*denominator>=2**63:
        raise ValueError('exact integer residual would exceed its checked representation')
    source=np.zeros((n,n,n))
    for site in sites: source[site]+=float(charge)
    continuous=green_solve(source)
    integer=np.rint(continuous*denominator).astype(np.int64)
    assert F(int(np.max(np.abs(integer))),denominator)<geometry_window
    lq=laplacian(integer)
    qlq=sum(int(x)*int(y) for x,y in zip(integer.ravel(),lq.ravel(),strict=True))
    energy=F(qlq,2*denominator**2)+charge*sum(F(int(integer[site]),denominator) for site in sites)
    residual=charge.denominator*lq
    for site in sites: residual[site]+=charge.numerator*denominator
    maxres=max(abs(int(residual.min())),abs(int(residual.max())))
    norm2_upper=F(n**3*maxres**2,charge.denominator**2*denominator**2)
    spectral_lower=F(12,(n+1)**2)
    gap=norm2_upper/(2*spectral_lower)
    assert gap<F(1,10**15)
    return {'minimum_energy_interval':[str(energy-gap),str(energy)],
            'dyadic_bits':bits,'geometry_max_absolute':str(F(int(np.max(np.abs(integer))),denominator)),
            'finite_geometry_window':str(geometry_window),
            'residual_squared_norm_upper':str(norm2_upper),
            'relaxation_to_finite_minimum_error_upper':str(gap),
            'finite_grid_minimum_is_found_exactly':False,
            'finite_grid_minimum_is_enclosed':True}


def interval_subtract(a,b,c):
    a,b,c=[[F(x) for x in row] for row in (a,b,c)]
    return [str(a[0]-b[1]-c[1]),str(a[1]-b[0]-c[0])]


def force_run():
    _,_,_,j=source_matrices()
    # A 1/10000 mixture of |001><001| into the fixed |000> reference.
    native_charge=F(j[1,1]-j[0,0]); mixture=F(1,10000)
    charge=native_charge*mixture
    assert native_charge==F(9,14) and charge==F(9,140000)
    rows=[]
    for n in (31,63,95):
        center=(n//2,)*3
        density=np.zeros((n,n,n));density[center]=1
        green=-green_solve(density)
        for radius in (4,6,8):
            left=list(center);left[0]+=radius-1
            right=list(center);right[0]+=radius+1
            point=list(center);point[0]+=radius
            magnitude=charge**2*F(float((green[tuple(left)]-green[tuple(right)])/2))
            ratio=float(magnitude/charge**2)*4*np.pi*radius*radius
            assert magnitude>0 and .85<ratio<1.2
            rows.append({'grid':n,'separation':radius,
                         'attractive_radial_force_magnitude':float(magnitude),
                         'force_divided_by_j_squared_over_4pi_r_squared':ratio,
                         'cross_energy':float(-charge**2*F(float(green[tuple(point)])))})
    n=31; center=(15,15,15); other=(21,15,15)
    pair=finite_energy_certificate(n,(center,other),charge)
    one=finite_energy_certificate(n,(center,),charge)
    two=finite_energy_certificate(n,(other,),charge)
    cross=interval_subtract(pair['minimum_energy_interval'],one['minimum_energy_interval'],two['minimum_energy_interval'])
    assert F(cross[1])<0
    source=np.zeros((n,n,n));source[center]=float(charge)
    q=green_solve(source)
    q_stiff=green_solve(source,stiffness=2)
    assert np.max(np.abs(q_stiff-q/2))<1e-18
    screened=green_solve(source,mass_squared=1)
    assert abs(screened[other])<abs(q[other])/10
    # The chosen occupation contrast has an exact positive remainder.
    qsymbol=s.symbols('q',real=True)
    r,actual,_,_=source_matrices()
    contrast=actual[1,1]-actual[0,0]
    remainder=s.factor(contrast.subs(r,1+qsymbol)-1-s.Rational(9,14)*qsymbol)
    assert s.simplify(remainder-qsymbol**2*(7+3*qsymbol)/(14*(1+qsymbol)))==0
    qmax=F(pair['geometry_max_absolute'])
    window=F(1,50000)
    coefficient=(7-3*window)/(14*(1-window))
    nonlinear_energy_error=2*mixture*coefficient*qmax*qmax
    # Positive nonlinear remainders: pair upper and sum of single uppers
    # each shift the cross-energy interval by at most two local remainders.
    uniform_cross_error=2*mixture*coefficient*window*window
    full_cross=[str(F(cross[0])-uniform_cross_error),str(F(cross[1])+uniform_cross_error)]
    assert F(full_cross[1])<0
    neighbour_cross=[]
    for distance in (5,7):
        position=(15+distance,15,15)
        pair_d=finite_energy_certificate(n,(center,position),charge)
        second_d=finite_energy_certificate(n,(position,),charge)
        interval_d=interval_subtract(pair_d['minimum_energy_interval'],
                                    one['minimum_energy_interval'],second_d['minimum_energy_interval'])
        neighbour_cross.append([F(interval_d[0])-uniform_cross_error,
                                F(interval_d[1])+uniform_cross_error])
    full_force=[(neighbour_cross[1][0]-neighbour_cross[0][1])/2,
                (neighbour_cross[1][1]-neighbour_cross[0][0])/2]
    assert full_force[0]>0
    return {'native_dimensionless_geometry_charge':str(native_charge),
            'matter_mixture':str(mixture),'effective_charge':str(charge),
            'geometry_action':'(kappa/2) sum_edges(q_x-q_y)^2 + sum_x j_x q_x',
            'declared_stiffness':'1','declared_geometry_mass_squared':'0',
            'boundary':'zero Dirichlet geometry, unit cubic spacing',
            'force_rows':rows,
            'finite_carrier_energy_certificates':{'pair':pair,'source_one':one,'source_two':two},
            'finite_carrier_cross_energy_interval':cross,
            'nonlinear_native_matter_remainder_at_pair_witness_upper':str(nonlinear_energy_error),
            'exact_positive_native_occupation_remainder':str(remainder),
            'uniform_nonlinear_cross_energy_error_upper':str(uniform_cross_error),
            'full_nonlinear_native_cross_energy_interval':full_cross,
            'full_nonlinear_attractive_force_interval_at_six':[str(x) for x in full_force],
            'nonlinear_remainder_bounds_the_full_finite_geometry_minimum':True,
            'strength_changes_under_stiffness_change':True,
            'positive_geometry_mass_screens_long_range_response':True,
            'same_sign_native_charges_attract':True,
            'opposite_sign_native_charges_repel':True,
            'physical_Newton_constant_or_universal_mass_coupling_derived':False,
            'native_quantum_controller_executes_geometry_minimizer':False}
