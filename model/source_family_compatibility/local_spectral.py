"""Complete local Heisenberg-CAR spectral records with original probabilities."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, isqrt
import copy

import sympy as s

from .defect_budget import exact_source
from .local_matter import annihilator, boundary_upper
from .cube_spectrum import pi_interval


def _rational(value, *, positive=False):
    if type(value) not in (int,F) or value < 0 or (positive and value == 0):
        raise ValueError('nonnegative exact rational required')
    return F(value)


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('integer outside the admitted domain')
    return value


def _ceil(value):
    return -((-value.numerator)//value.denominator)


def leakage(count, step, distance):
    """Rational Fejer tail, using pi<22/7 and retaining original mass."""
    n=_integer(count,2); dt=_rational(step,positive=True)
    d=_rational(distance,positive=True)
    a=F(7*n,44)*dt*d
    return min(F(1),1/(2*(a-1))) if a>1 else F(1)


def local_filter_bound(count,step,radius,support_radius,support_sites,coupling,
                       marginal_error=0,service_error=0):
    n=_integer(count,2);dt=_rational(step,positive=True)
    c=_rational(coupling,positive=True)
    delta=_rational(marginal_error);service=_rational(service_error)
    if delta>1 or service>1:
        raise ValueError('half-trace errors must be at most one')
    time=(n-1)*dt
    boundary=boundary_upper(radius,support_radius,support_sites,time,c)
    # sqrt(2)<3/2 bounds the entire Stinespring difference, not each record.
    complete=min(F(1),F(3,2)*boundary+delta+service)
    return {'phase_count':n,'charge_phase_records':2*n,'time_step':dt,
        'source_duration':time,'cut_radius':radius,'support_radius':support_radius,
        'support_sites':support_sites,'boundary_operator_error':boundary,
        'marginal_halftrace_error':delta,'complete_service_error':service,
        'complete_output_halftrace_upper':complete,
        'error_has_phase_count_prefactor':False,
        'complete_records_and_untouched_reference_retained':True}


def plan_filter(count,step,support_radius,support_sites,coupling,tolerance):
    n=_integer(count,2);dt=_rational(step,positive=True)
    r=_integer(support_radius);_integer(support_sites,1)
    c=_rational(coupling,positive=True);eta=_rational(tolerance,positive=True)
    if eta>1 or support_sites>(2*r+1)**3:
        raise ValueError('valid support and accuracy at most one required')
    radius=r+max(1,_ceil(96*c*(n-1)*dt))
    result=local_filter_bound(n,dt,radius,r,support_sites,c)
    while result['complete_output_halftrace_upper']>eta:
        radius+=1
        result=local_filter_bound(n,dt,radius,r,support_sites,c)
    result.update({'requested_accuracy':eta,'cut_sites':(2*radius+1)**3,
        'cut_Fock_dimension_exponent':2*(2*radius+1)**3,
        'ground_marginal_prepared_by_planner':False,
        'large_native_tape_executed':False})
    return result


def window_from_records(count,step,first,last,distance,energy_cut,
                        probability_lower,original_mass_upper,uv_tail_upper,
                        complete_error=0):
    """A nonwrapped real window from a complete original record probability."""
    n=_integer(count,2);dt=_rational(step,positive=True)
    _integer(first);_integer(last)
    if first>last or last>n//2:
        raise ValueError('ordered phase centres in the first half circle required')
    d=_rational(distance,positive=True);cut=_rational(energy_cut,positive=True)
    p=_rational(probability_lower);mass=_rational(original_mass_upper)
    uv=_rational(uv_tail_upper);error=_rational(complete_error)
    if any(x>1 for x in (p,mass,uv,error)) or p>mass+error:
        raise ValueError('admissible original probability bounds required')
    pl,pu=pi_interval()
    if cut>pl/dt:
        raise ValueError('UV cut must be in the first half of the phase period')
    lo=max(F(0),2*pl*first/(n*dt)-d)
    hi=2*pu*last/(n*dt)+d
    tail=leakage(n,dt,d)
    return {'phase_records':[first,last+1],
        'true_positive_energy_window':[lo,hi],
        'UV_cut':cut,'original_UV_tail_upper':uv,
        'Fejer_leakage_upper':tail,'original_mass_upper':mass,
        'complete_output_error':error,
        'original_true_window_probability_lower':max(F(0),p-error-tail*mass-uv),
        'sharp_projection_implemented':False,'original_weights_retained':True}


def fourth_taper(length):
    """Exact fourth-convolution amplitudes, padded to a complete Fourier register."""
    m=_integer(length,1);count=4*m
    def choose3(x):return comb(x,3) if x>=3 else 0
    coefficients=tuple(sum((-1)**j*comb(4,j)*choose3(l-j*m+3)for j in range(5))
                       for l in range(count))
    square=sum(a*a for a in coefficients)
    if min(coefficients)<0 or sum(coefficients)!=m**4 or square*4<m**7:
        raise ArithmeticError('source time-register normalization failed')
    return coefficients,square


def fourth_taper_leakage(count,step,distance):
    n=_integer(count,4);dt=_rational(step,positive=True);d=_rational(distance,positive=True)
    if n%4:raise ValueError('fourfold padded Fourier register required')
    a=F(7*n,44)*dt*d
    return min(F(1),F(4096,7)/a**7) if a>=1 else F(1)


def default_inverse_window():
    """A complete default-source readout with a positive original true-window floor."""
    from .uv_tail import tail_budget
    p=F(25,100352);n=8192;step=F(1,16);d=F(1,4)
    uv=tail_budget(F(1,8),1,6,20)['original_total_charge_tail_upper']
    if uv!=F(1647837513,156250000000000) or uv>=p/16:
        raise ArithmeticError('independent source ultraviolet budget differs')
    coefficients,square=fourth_taper(n//4)
    tail=fourth_taper_leakage(n,step,d)
    if tail!=F(4096,7)*F(11,224)**7 or tail>=p/16:
        raise ArithmeticError('tapered clock leakage budget differs')
    pl,pu=pi_interval()
    # The last included centre is13pi/2. At E<=20 the excluded upper
    # centre is at least1/4 away; the inherited band starts above1/4.
    if F(13,2)*pl-20<d:
        raise ArithmeticError('known-band/group separation failed')
    plan=plan_filter(n,step,0,1,F(1,8),p/16)
    error=plan['complete_output_halftrace_upper']
    ideal_floor=(1-tail)*(p-uv)
    measured_floor=ideal_floor-error
    true_floor=measured_floor-error-tail/2-uv
    if true_floor<=3*p/4:
        raise ArithmeticError('original true-window floor failed')
    return {'source_c':F(1,8),'source_U':F(1),'source_U_over_c':F(8),
        'original_source_band_weight_lower':p,'original_total_UV_tail_above_20':uv,
        'source_UV_bound_is_conditional':False,'time_amplitude_convolutions':4,
        'rectangle_length':n//4,'phase_count':n,'charge_phase_records':2*n,
        'time_amplitude_sum':sum(coefficients),'time_amplitude_square_sum':square,
        'time_register_normalized_exactly':True,'phase_records':[0,1665],
        'phase_step':step,'excluded_centre_distance_lower':d,
        'taper_leakage_upper':tail,'each_ideal_group_probability_lower':ideal_floor,
        'each_local_measured_group_probability_lower':measured_floor,
        'true_energy_window':[F(0),F(13,2)*pu+F(1,4)],
        'each_original_true_window_probability_lower':true_floor,
        'published_each_original_true_window_probability_lower':3*p/4,
        'local_cut_plan':plan,'all_other_records_retained':True,
        'same_full_ground_marginal_required':True,
        'native_tapered_clock_executed':False,'physical_mass_assigned':None}


def spectral_score_bound(count,step,distance,energy_cut,original_mass_upper,
                         uv_tail_upper,lipschitz=1,score_norm=1,complete_error=0):
    """Bound a complete signed-centre score against the original energy measure."""
    n=_integer(count,2);dt=_rational(step,positive=True)
    d=_rational(distance,positive=True);cut=_rational(energy_cut,positive=True)
    mass=_rational(original_mass_upper);uv=_rational(uv_tail_upper)
    ell=_rational(lipschitz);norm=_rational(score_norm)
    eta=_rational(complete_error)
    if max(mass,uv,eta)>1:
        raise ValueError('original probability and output errors at most one required')
    pl,_=pi_interval()
    if cut+d>pl/dt:
        raise ValueError('energy cut and guard must fit inside the signed half-period')
    tail=leakage(n,dt,d)
    return {'original_spectral_score_error_upper':min(2*norm,ell*d*mass+2*norm*(tail*mass+uv+eta)),
        'Fejer_leakage_upper':tail,'original_UV_tail_upper':uv,
        'complete_output_error':eta,'conditional_renormalization_applied':False,
        'negative_signed_centres_clamped_to_zero':True}


def default_positive_event():
    # The inherited source inequality is for each charge separately.
    c,u,degree=F(1,8),F(1),6
    second=(u*u+degree*c*c)/4
    derivative=(u+degree*c)**2
    p=second**2/(32*derivative**2)
    n=4096;dt=F(1,16);distance=F(1,4)
    squared=(second/2,16*derivative**2/second)
    if p!=F(25,100352) or squared!=(F(35,256),F(2744,5)):
        raise ArithmeticError('inherited local charged moment band differs')
    if not F(1,4)**2<=squared[0] or not squared[1] <= F(47,2)**2:
        raise ArithmeticError('actual source band is outside the admitted interval')
    # First excluded upper centre >8pi>24; the other tail is beyond period.
    if not F(24)-F(47,2)>=distance:
        raise ArithmeticError('group separation failed')
    tail=leakage(n,dt,distance)
    if tail!=F(11,202):raise ArithmeticError('complete phase leakage differs')
    plan=plan_filter(n,dt,0,1,F(1,8),p/4)
    lower=(1-tail)*p-plan['complete_output_halftrace_upper']
    if lower<=p/2:raise ArithmeticError('positive local event floor failed')
    return {'source_c':F(1,8),'source_U':F(1),'source_U_over_c':F(8),
        'local_degree':6,'same_full_ground_required':True,
        'each_original_source_band_weight_lower':p,
        'source_energy_band_squared':list(squared),'phase_records':[0,1025],
        'excluded_centre_distance_lower':distance,'Fejer_leakage_upper':tail,
        'each_ideal_group_probability_lower':(1-tail)*p,
        'each_cut_group_probability_lower':lower,
        'published_each_cut_group_probability_lower':p/2,
        'all_other_charge_phase_records_retained':True,
        'outer_source_volume_independent':True,'local_cut_plan':plan,
        'cut_ground_substituted':False,'shrinking_unit_hopping_energy_floor':False}


def _simple(matrix):
    return matrix.applyfunc(s.simplify)


@lru_cache(maxsize=1)
def dimer_source():
    """Rebuild and verify the full source, ground and all spectral projectors."""
    h,_,_=exact_source(2,((0,1,1),),F(1,8),F(1))
    root=s.sqrt(17);e0=(4-root)/8
    v=(h-e0*s.eye(16)).nullspace()[0]
    rho=s.ImmutableMatrix(_simple(v*v.H/(v.H*v)[0]))
    energies=tuple(h.eigenvals())
    projections=[]
    for e in energies:
        p=s.eye(16)
        for other in energies:
            if other!=e:p=_simple(p*(h-other*s.eye(16))/(e-other))
        p=s.ImmutableMatrix(p)
        if _simple(p*p-p)!=s.zeros(16) or p!=p.H or _simple(h*p-e*p)!=s.zeros(16):
            raise ArithmeticError('full source projector identity failed')
        projections.append((s.simplify(e-e0),p))
    if sum((p for _,p in projections),s.zeros(16))!=s.eye(16):
        raise ArithmeticError('complete source spectrum missing')
    a=annihilator(2,0)
    if a.H*a+a*a.H!=s.eye(16):raise ArithmeticError('source CAR failed')
    physical_energies=(root/8-s.Rational(1,16),root/8+s.Rational(1,16))
    weights=(s.Rational(1,4)+1/(4*root),s.Rational(1,4)-1/(4*root))
    for op in (a.H,a):
        for e,p in projections:
            weight=s.simplify(s.trace(rho*op.H*p*op))
            expected=next((w for x,w in zip(physical_energies,weights,strict=True)
                           if s.simplify(e-x)==0),s.Integer(0))
            if s.simplify(weight-expected)!=0:
                raise ArithmeticError('actual original charged weight differs')
    return h,rho,a,tuple(projections)


def _round_dyadic(value,bits=60):
    scale=2**bits
    return F((value*scale+F(1,2)).numerator//(value*scale+F(1,2)).denominator,scale)


def _phase_approximation(midpoint,radius):
    """A rational complex phase with a certified 2^-56 complex error."""
    if abs(midpoint)>4 or radius<0:
        raise ValueError('phase midpoint/domain outside finite control')
    cosine=sum(((-1)**k*midpoint**(2*k)/factorial(2*k)for k in range(32)),F())
    sine=sum(((-1)**k*midpoint**(2*k+1)/factorial(2*k+1)for k in range(32)),F())
    remainder=max(abs(midpoint)**64/F(factorial(64)),
                  abs(midpoint)**65/F(factorial(65)))+radius
    re,im=_round_dyadic(cosine),_round_dyadic(sine)
    # Sum of coordinate errors bounds the complex modulus without a square root.
    if abs(re-cosine)+abs(im-sine)+2*remainder>F(1,2**56):
        raise ArithmeticError('rational phase admission failed')
    return re,im,F(1,2**56)


def _kernel_interval(energy_interval,j,count,step,time_amplitudes=None):
    lo,hi=energy_interval;pl,pu=pi_interval()
    signed=j if j<=count//2 else j-count
    midpoint=2*(pl+pu)/2*F(signed,count)-(lo+hi)/2*step
    radius=(pu-pl)*abs(F(signed,count))+(hi-lo)*step/2
    re,im,error=_phase_approximation(midpoint,radius)
    if count*error>F(1,2):raise ValueError('phase count exceeds propagation admission')
    coefficients=(1,)*count if time_amplitudes is None else tuple(time_amplitudes)
    if len(coefficients)!=count or any(type(x)is not int or x<0 for x in coefficients):
        raise ValueError('nonnegative integer time amplitudes required')
    square=sum(x*x for x in coefficients)
    if not square:raise ValueError('nonzero time register required')
    xr,xi=F(1),F(0);sr,si=F(0),F(0)
    for coefficient in coefficients:
        sr+=coefficient*xr;si+=coefficient*xi
        xr,xi=xr*re-xi*im,xr*im+xi*re
    estimate=(sr*sr+si*si)/(count*square)
    # Cauchy gives sum a_l/sqrt(n*S)<=1. Each power error is at most
    # 2l*error, so the complete weighted amplitude error is <=2(n-1)error.
    allowance=6*(count-1)*error
    return max(F(0),estimate-allowance),min(F(1),estimate+allowance)


def _outward(lo,hi,grid=10**12):
    return [F((lo*grid).numerator//(lo*grid).denominator,grid),
            F(_ceil(hi*grid),grid)]


def finite_records(*,tapered=False):
    dimer_source()
    scale=2**80;r=isqrt(17*scale*scale)
    low,high=F(r,scale),F(r+1,scale)
    if not low*low<17<high*high:raise ArithmeticError('root enclosure failed')
    energies=((low/8-F(1,16),high/8-F(1,16)),
              (low/8+F(1,16),high/8+F(1,16)))
    weights=((F(1,4)+1/(4*high),F(1,4)+1/(4*low)),
             (F(1,4)-1/(4*low),F(1,4)-1/(4*high)))
    count=32;step=F(1);rows=[]
    amplitudes,square=fourth_taper(8) if tapered else ((1,)*count,count)
    for j in range(count):
        kl,kh=(_kernel_interval(e,j,count,step,amplitudes)for e in energies)
        lo=weights[0][0]*kl[0]+weights[1][0]*kh[0]
        hi=weights[0][1]*kl[1]+weights[1][1]*kh[1]
        rows.append({'phase_record':j,'each_charge_original_probability':_outward(lo,hi)})
    if not sum(x['each_charge_original_probability'][0]for x in rows)<=F(1,2)<=sum(x['each_charge_original_probability'][1]for x in rows):
        raise ArithmeticError('complete original probability sum not enclosed')
    return {'source_c':F(1,8),'source_U':F(1),'source_U_over_c':F(8),
        'full_Fock_dimension':16,'all_source_energy_projectors':6,
        'local_mode':0,'same_ground_same_clock':True,'phase_count':count,
        'time_step':step,'source_duration':(count-1)*step,
        'time_amplitude_convolutions':4 if tapered else 1,
        'time_amplitude_square_sum':square,
        'charge_labels':['addition','removal'],'complete_record_count':2*count,
        'each_charge_original_total_weight':F(1,2),'complete_original_total_weight':F(1),
        'actual_positive_energy_intervals':[_outward(*x)for x in energies],
        'actual_original_weight_intervals':[_outward(*x)for x in weights],
        'all_phase_records':rows,'ground_energy_needed_by_instrument':False,
        'all_input_completeness':'Fourier Parseval and the full source CAR identity',
        'phase_clock_is_sharp_spectral_PVM':False,
        'excited_input_has_ground_spectral_interpretation':False,
        'native_transfer_or_large_clock_tape_executed':False,
        'physical_particle_or_mass_identified':False}


def _encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:_encode(v)for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [_encode(v)for v in value]
    return value


@lru_cache(maxsize=1)
def _certificate():
    return _encode({'actual_dimer_complete_clock':finite_records(),
        'actual_dimer_tapered_complete_clock':finite_records(tapered=True),
        'growing_default_source_local_positive_event':default_positive_event(),
        'default_complete_inverse_energy_window':default_inverse_window(),
        'local_clock_isometry_error':'sqrt(2)*B_R(T)',
        'phase_count_prefactor_in_local_error':False,
        'specified_ground_marginal_error_required':True,
        'infinite_ground_preparation_or_subsequence_modulus_derived':False,
        'shrinking_unit_hopping_spectral_weight_proved':False,
        'spatial_Walsh_label_identified_as_momentum':False,
        'physical_mass_prediction':None})


def certificate():
    # Encoded cached values contain only ordinary immutable scalars.
    return copy.deepcopy(_certificate())
