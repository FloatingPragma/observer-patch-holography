"""Finite compensated matter hopping with transverse Peierls matrices."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import math
import sympy as s
import mpmath as mp

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
PARENT='model/affine_gauss_source/results.json'
from model.affine_gauss_source import model as parent
from model.affine_gauss_source.certificate import certificate as parent_certificate, validate_pins
ETA=F(1,24)


def rational(value):
    return s.Rational(value.numerator,value.denominator) if isinstance(value,F) else s.Rational(value)


def endpoint(interval,upper=False):
    sign,mantissa,exponent,_=interval._mpi_[int(upper)]
    return F((-1)**sign*mantissa)*F(2)**exponent


def ivalue(expression):
    expression=s.sympify(expression);iv=mp.iv
    if expression.is_Rational:return iv.mpf(int(expression.p))/int(expression.q)
    if expression.is_Add:return sum((ivalue(x) for x in expression.args),iv.mpf(0))
    if expression.is_Mul:return math.prod(ivalue(x) for x in expression.args)
    if expression.is_Pow:
        base,power=expression.args
        if power==s.Rational(1,2):return iv.sqrt(ivalue(base))
        if power==s.Rational(-1,2):return 1/iv.sqrt(ivalue(base))
        if power.is_Rational and power.q==2:return iv.sqrt(ivalue(base))**int(power.p)
        if power.is_Integer:return ivalue(base)**int(power)
    raise ValueError('unsupported exact algebraic interval expression: '+str(expression))


def rounded(interval,places=24):
    scale=10**places;low=endpoint(interval);high=endpoint(interval,True)
    return [str(F(low.numerator*scale//low.denominator,scale)),
            str(F(-(-high.numerator*scale//high.denominator),scale))]


def source():
    path=RER/PARENT
    receipt=json.loads(path.read_text())
    validate_pins(receipt['input_sha256'])
    if receipt!=parent_certificate():
        raise ValueError('fresh public static parent replay differs')
    return receipt,hashlib.sha256(path.read_bytes()).hexdigest()


def actual_edges():
    boundary=parent.geometry()[0]
    return [(0,j+1) for j in range(12)]+[(u+1,v+1) for u,v in boundary]


def universal_bounds():
    edges=actual_edges();degree=[sum(y in e for e in edges) for y in range(13)]
    assert degree==[12]+[6]*12
    gap=s.sqrt(4-s.sqrt(5))
    assert 2-12*rational(ETA)==s.Rational(3,2)>gap
    rows=[]
    charges,_=parent.charge_geometry()
    for y,d in enumerate(degree):
        residual_squared=ETA**2*d
        w=F(1,4);outside=residual_squared/w**2
        assert outside<=F(1,3)
        rows.append({'source_vertex':y,'degree':d,
            'charged_state_energy_mean':charges[y]['total_energy_above_actual_neutral_ground'],
            'source_static_electric_mean':charges[y]['longitudinal_electric_field'],
            'retained_reference_energy':'1',
            'squared_actual_H_energy_residual':str(residual_squared),
            'energy_bin_halfwidth':'1/4','conditional_energy_bin_outside_upper':str(outside),
            'unconditional_pump_and_energy_bin_hit_lower':str((1-outside)/4),
            'charge_pump_failure':'3/4'})
    t=F(1,4);sinc=1-F(10**2)*t*t/24
    first=ETA*t*sinc;remainder=24*ETA**2*t*t;amplitude=first-remainder
    assert (sinc,first,remainder,amplitude)==(F(71,96),F(71,9216),F(1,384),F(47,9216))
    tp=F(1,64)
    assert 1-F(10**2)*tp*tp/24>F(11,12)
    photon_first=ETA*tp/12;photon_remainder=24*ETA**2*tp*tp
    photon_amplitude=photon_first-photon_remainder
    assert photon_amplitude==F(13,294912)
    PT=parent.geometry()[7]
    assert set(PT.diagonal())=={s.Rational(165,217),s.Rational(151,217)}
    assert min(PT.diagonal())>s.Rational(2,3)
    return {'hopping_eta':str(ETA),'degree':degree,'hopping_form_lower':'-12 eta N',
        'physical_H_minus_vacuum_lower':'H_ph+E_C+(2-12 eta)N',
        'unique_actual_physical_ground_unchanged':True,
        'actual_physical_ground_gap':gap,
        'charged_sector_energy_lower':'3/2',
        'neutral_field_records_unchanged_exactly':True,
        'static_exact_charged_energy_line_law_transferred':False,
        'energy_records':rows,
        'neighbor_arrival':{'time':str(t),'charge_success_probability':'1/4',
            'first_Dyson_vector_norm_lower':str(first),'all_orders_remainder_upper':str(remainder),
            'actual_arrival_vector_norm_lower':str(amplitude),
            'unconditional_probability_lower':str(amplitude**2/4),
            'exact_complete_probability_law':['3/4','p=||P_x U_eta(t) psi_y||^2/4','1/4-p'],
            'arrival_probability_interval':[str(amplitude**2/4),'1/4'],
            'nonarrival_probability_interval':['0',str(F(1,4)-amplitude**2/4)],
            'complete_records':['charge_pump_failure','charge_success_arrival','charge_success_nonarrival']},
        'neighbor_arrival_and_transverse_occupation':{'time':str(tp),
            'Peierls_vacuum_transverse_excitation_probability_lower':'1/108',
            'first_Dyson_projected_norm_lower':str(photon_first),
            'all_orders_remainder_upper':str(photon_remainder),
            'actual_projected_norm_lower':str(photon_amplitude),
            'unconditional_probability_lower':str(photon_amplitude**2/4),
            'complete_records':['charge_pump_failure','arrival_vacuum_transverse','arrival_nonvacuum_transverse',
                'nonarrival_vacuum_transverse','nonarrival_nonvacuum_transverse'],
            'asymptotic_radiation_claim':False}}


def total_occupations(modes,cutoff):
    for number in range(cutoff+1):
        for word in itertools.combinations_with_replacement(range(modes),number):
            yield tuple(word.count(j) for j in range(modes))


@lru_cache(None)
def exact_peierls_data(edge):
    _,_,_,D,C,G,K,PT,powers,polynomials=parent.geometry()
    vector=s.eye(42)[:,edge]
    weights=[parent.projected_weight(vector,p,powers) for p in polynomials[1:]]
    frequencies=[s.sqrt(lam) for lam in parent.ROOTS[1:]]
    f_squared=[parent.simple(w/(2*omega)) for w,omega in zip(weights,frequencies)]
    variance=parent.simple(sum(f_squared))
    assert variance>0 and variance<s.Rational(1,2)
    assert parent.simple(sum(weights)-PT[edge,edge])==0
    return frequencies,f_squared,variance,PT[edge,edge]


def first_dyson_interval(edge,y,x,time):
    """Exact 4-level exponential plus finite frequency-shell multinomials.

    This computes the first integral only; the separate unitary Duhamel
    bound controls every higher hopping insertion of the actual H.
    """
    mp.iv.dps=80;iv=mp.iv
    frequencies,f_squared,variance,transverse=exact_peierls_data(edge)
    sigma=iv.sqrt(ivalue(variance));sqrt6=iv.sqrt(6)
    alpha=iv.sqrt(3-sqrt6);beta=iv.sqrt(3+sqrt6)
    ca=iv.cos(sigma*alpha);cb=iv.cos(sigma*beta)
    sa=iv.sin(sigma*alpha)/alpha;sb=iv.sin(sigma*beta)/beta
    odd3=(sb-sa)/(2*sqrt6);odd1=sa-(3-sqrt6)*odd3
    even2=(cb-ca)/(2*sqrt6);even0=ca-(3-sqrt6)*even2
    # Absolute amplitudes in occupations0..3; odd ones carry a common i.
    amplitudes=[even0+even2,odd1+3*odd3,iv.sqrt(2)*even2,sqrt6*odd3]
    norm=sum((a*a for a in amplitudes),iv.mpf(0))
    assert endpoint(norm)<=1<=endpoint(norm,True)
    probabilities=[ivalue(f/variance) for f in f_squared]
    assert endpoint(sum(probabilities,iv.mpf(0)))<=1<=endpoint(sum(probabilities,iv.mpf(0)),True)
    charges,_=parent.charge_geometry();delta0=charges[x]['static_electric_energy']-charges[y]['static_electric_energy']
    t=ivalue(rational(time));eta=ivalue(rational(ETA))
    first_squared=iv.mpf(0);photon_first_squared=iv.mpf(0);weight_sum=iv.mpf(0);count=0
    for n in total_occupations(8,3):
        number=sum(n)
        factor=math.factorial(number)//math.prod(math.factorial(value) for value in n)
        weight=amplitudes[number]**2*factor
        for p,value in zip(probabilities,n):weight*=p**value
        energy=s.sympify(delta0)+sum(value*w for value,w in zip(n,frequencies))
        delta=ivalue(energy)
        response=t*t if energy==0 else (2*iv.sin(delta*t/2)/delta)**2
        contribution=eta*eta*weight*response
        first_squared+=contribution
        if number:photon_first_squared+=contribution
        weight_sum+=weight;count+=1
    assert count==165 and endpoint(weight_sum)<=1<=endpoint(weight_sum,True)
    degree=sum(y in e for e in actual_edges())
    remainder=t*t*12*eta*eta*iv.sqrt(degree)/2
    actual=iv.sqrt(first_squared)-remainder
    photon_actual=iv.sqrt(photon_first_squared)-remainder
    assert endpoint(actual)>0 and endpoint(photon_actual)>0
    probability=actual*actual/4;photon_probability=photon_actual*photon_actual/4
    return {'actual_oriented_edge':edge,'source_vertex':y,'arrival_vertex':x,'time':str(time),
        'first_integral_frequency_configurations':count,
        'Peierls_vacuum_probability_norm_interval':rounded(norm),
        'Peierls_nonvacuum_probability_interval':rounded(1-amplitudes[0]**2),
        'first_arrival_norm_interval':rounded(iv.sqrt(first_squared)),
        'first_arrival_nonvacuum_norm_interval':rounded(iv.sqrt(photon_first_squared)),
        'all_orders_remainder_interval':rounded(remainder),
        'unconditional_arrival_probability_lower':str(endpoint(probability)),
        'unconditional_arrival_nonvacuum_probability_lower':str(endpoint(photon_probability)),
        'unconditional_arrival_probability_lower_decimal':str(mp.mpf(endpoint(probability).numerator)/endpoint(probability).denominator),
        'unconditional_arrival_nonvacuum_probability_lower_decimal':str(mp.mpf(endpoint(photon_probability).numerator)/endpoint(photon_probability).denominator),
        'actual_H_not_replaced_by_first_order_evolution':True}


def calculate():
    parent_receipt,sha=source();mp.mp.dps=50
    bounds=universal_bounds()
    edges=actual_edges();assert edges[0]==(0,1) and edges[12]==(1,2)
    controls=[first_dyson_interval(edge,1,x,time) for edge,x in ((0,0),(12,2)) for time in (F(1,4),F(1,64))]
    result={'schema':'oph.quantum_worker.finite_coulomb_mobile.v1',
        'parent_receipt_sha256':sha,
        'source':{'carrier':'same finite physical carrier as static parent; total photon occupation<=3',
            'Hamiltonian':'H_eta=H_static-eta sum_(u,v)[d_v^dag d_u exp(i A_e)+h.c.]',
            'eta':str(ETA),'d_y':'compensated annihilator including fixed-reference and longitudinal-register actions',
            'Peierls_unitary':'literal finite exp(i A_e); no canonical Weyl or BCH identities used',
            'reference_position':1,'reference_energy':'R','source_clock':'exp(-i t H_eta)',
            'actual_ground_preparation_success':'1','physical_ground_is_unique':True,
            'physical_CAR_modes':13,'physical_carrier_dimension':44695552,
            'one_charge_sector_dimension':13*5456,
            'retained_detector_pointer_dimension':12,
            'kinematic_source_reference_and_pointer_dimension':parent_receipt['architecture']['kinematic_source_reference_and_pointer_dimension'],
            'unused_pointer_outcomes_retained':True,
            'other_controller_storage_included':False,
            'all_source_Gauss_constraints_retained':True,
            'charge_current_continuity':'dot n=D^T J; D^T dot E=dot n on the physical source',
            'compensator_moves':False,'transverse_modes_free_when_N_zero':True,
            'physical_spin_or_measured_charge_assignment':False,
            'thermodynamic_Coulomb_phase_or_continuum_scattering':False,
            'source_is_PW119':False,
            'energy_arrival_and_joint_occupation_are_distinct_detector_settings':True},
        'uniform_exact_bounds':bounds,
        'complete_neutral_field_records_same_H':parent_receipt['complete_field_records'],
        'finite_Peierls_first_integral_controls':controls,
        'owned_source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}}
    return parent.serialize(result)
