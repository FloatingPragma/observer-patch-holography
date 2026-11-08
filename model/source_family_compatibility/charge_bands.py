"""Original charged spectral weights of the centered full Hubbard source."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from .algebra import creation
from .defect_budget import exact_source


def rational(value, positive=False):
    if type(value) not in (int,F) or value<0 or (positive and value==0):
        raise ValueError('exact nonnegative rational required')
    return F(value)


def broad_band(coupling,contact,degree,degree_bound=None):
    c,u=rational(coupling,True),rational(contact)
    if type(degree) is not int or degree<1:
        raise ValueError('positive integer local degree required')
    maximum=degree if degree_bound is None else degree_bound
    if type(maximum) is not int or maximum<degree:
        raise ValueError('integer maximum degree at least local degree required')
    second=(u*u+degree*c*c)/4
    derivative=(u+maximum*c)**2
    return {'source_c':c,'source_U':u,'local_degree':degree,
        'maximum_degree':maximum,
        'each_original_preparation_weight':F(1,2),
        'each_original_second_moment':second/2,
        'double_commutator_norm_upper':derivative,
        'positive_energy_window_squared':[second/2,16*derivative**2/second],
        'each_original_band_probability_lower':second**2/(32*derivative**2),
        'centered_PH_and_hidden_spin_ground_required':True,
        'band_is_a_spectral_gap_or_particle_mass':False}


def _simple(matrix):
    return matrix.applyfunc(s.simplify)


@lru_cache(maxsize=None)
def dimer(coupling=F(1,8),contact=F(1)):
    c,u=rational(coupling,True),rational(contact,True)
    h,_,_=exact_source(2,((0,1,1),),c,u)
    cs,us=s.Rational(c),s.Rational(u)
    root=s.sqrt(us**2+4*cs**2); e0=(us-root)/2
    v=(h-e0*s.eye(16)).nullspace()[0]
    rho=_simple(v*v.H/(v.H*v)[0]); clock=h-e0*s.eye(16)
    creators=[creation(j) for j in range(4)]
    one=s.zeros(4)
    for a,b in ((0,2),(1,3)):one[a,b]=one[b,a]=cs/2
    moments=[]; power=s.eye(16)
    for k in range(5):
        moment=s.zeros(4)
        for i in range(4):
            for j in range(4):
                moment[i,j]=s.simplify(s.trace(rho*(creators[i].H*power*creators[j]
                    +(-1)**k*creators[j]*power*creators[i].H)))
        moments.append(moment);power=_simple(power*clock)
    identity=s.eye(4); alpha=us/2
    if moments[0]!=identity or moments[1]!=one or moments[2]!=one**2+alpha**2*identity:
        raise ArithmeticError('actual centered CAR moments differ')
    composite=_simple((moments[3]-one**3-2*alpha**2*one)/alpha**2)
    closure=_simple((moments[4]-one*moments[3]-moments[3]*one+one*moments[2]*one)/alpha**2
                    -alpha**2*identity-composite**2)
    if composite!=-3*one or closure!=s.zeros(4):
        raise ArithmeticError('complete source-generated two-response closure differs')
    if moments[3] != one**3-alpha**2*one:
        raise ArithmeticError('actual third spectral moment differs')
    second=(us**2+cs**2)/4
    for a in creators:
        derivative=_simple(h*a-a*h)
        if _simple(derivative.H*derivative+derivative*derivative.H)!=second*s.eye(16):
            raise ArithmeticError('full-source charged second-moment identity differs')
    expected_energies=[(root-cs)/2,(root+cs)/2]
    expected_weights=[s.Rational(1,4)+cs/(2*root),s.Rational(1,4)-cs/(2*root)]
    spectrum=list(h.eigenvals())
    projectors=[]
    for energy in spectrum:
        projector=s.eye(16)
        for other in spectrum:
            if energy!=other:projector=_simple(projector*(h-other*s.eye(16))/(energy-other))
        projectors.append((s.simplify(energy-e0),projector))
    if sum((p for _,p in projectors),s.zeros(16))!=s.eye(16):
        raise ArithmeticError('source energy resolution is incomplete')
    records=[]; total_effect=s.zeros(16)
    for mode,a in enumerate(creators):
        for label,op in (('addition',a),('removal',a.H)):
            for energy,projection in projectors:
                effect=_simple(op.H*projection*op/4);total_effect+=effect
                probability=s.simplify(s.trace(rho*effect))
                expected=next((w/4 for e,w in zip(expected_energies,expected_weights)
                               if s.simplify(e-energy)==0),s.Integer(0))
                if s.simplify(probability-expected)!=0:
                    raise ArithmeticError('complete source probability differs')
                records.append({'mode':mode,'charge':label,'energy':str(energy),
                                'original_probability':str(probability)})
    if _simple(total_effect)!=s.eye(16):raise ArithmeticError('incomplete charged instrument')
    return {'source_c':str(c),'source_U':str(u),'full_Fock_dimension':16,
        'centered_H_ground_energy':str(e0),
        'source_positive_charge_energies':[str(x) for x in expected_energies],
        'each_mode_each_charge_original_weights':[str(x) for x in expected_weights],
        'complete_records':records,'complete_record_count':len(records),
        'uniform_setting_probability':'1/4','zero_probability_records_retained':True,
        'original_probability_sum':'1',
        'source_retarded_moments_0_through_4':[[[str(x) for x in row] for row in m.tolist()] for m in moments],
        'charged_second_moment_sum_rule_on_whole_Fock_space':True,
        'actual_composite_clock':'B=-3h',
        'second_response_closes_on_this_dimer':True,
        'Dirac_composite_clock_B_equals_minus_h':False,
        'Dirac_third_moment_error':'-2*(U/2)^2*h',
        'spatial_continuum_or_measured_mass_assigned':False}


def certificate():
    def encode(value):
        if isinstance(value,F):return str(value)
        if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [encode(v) for v in value]
        return value
    return encode({'actual_dimer':dimer(),
        'cubic_readable_parameter_band':broad_band(F(1,8),1,6),
        'finite_and_infinite_ground_same_original_weight_bound':True,
        'charged_pole_mass_or_Lorentz_covariance_derived':False})
