"""Certified complete finite-resolution outputs from the frozen 78-line source."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

from tools import strict_json
from model.certified_reference import interval as iv
from . import kernels as k


ROOT=Path(__file__).resolve().parents[2]


def enclosure(lower,upper):
    return iv.Interval(iv.Interval.point(lower).lo,iv.Interval.point(upper).hi)


def load_source():
    result=strict_json.load(ROOT/'model/retarded_matter_response/results.json')
    for name,digest in result['source_bindings'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('frozen spectral source changed: '+name)
    section=result['resolved_local_spectrum'];rows=section['rows']
    if (result.get('all_pass') is not True or len(rows)!=78 or
            section['unresolved_local_probe_weight']!='0'):
        raise ValueError('complete certified local source spectrum required')
    channel=next(x for x in result['source_response']['channels']
                 if x['reference']=='ground' and x['channel']=='down')
    if [channel['moments'][str(n)][0][0] for n in range(3)]!=[
            ['1','0'],['1/2','0'],['131/256','0']]:
        raise ValueError('original source sum rules changed')
    values=[]
    for row in rows:
        lo,hi=map(F,row['gap_above_ground_interval'])
        low,high=map(F,row['normalized_weight_interval'])
        if not 0<lo<=hi<F(37,20) or not 0<low<=high<=1:
            raise ValueError('invalid source spectral enclosure')
        values.append((lo,hi,low,high))
    return values


@lru_cache(maxsize=2048)
def characteristic(time):
    time=k.exact(time,'sample time')
    real,imag=iv.Interval.point(0),iv.Interval.point(0)
    for lo,hi,low,high in load_source():
        sine,cosine=k.trig_interval(time*lo,time*hi)
        weight=enclosure(low,high)
        real+=weight*cosine;imag-=weight*sine
    return real,imag


def square(interval):
    lo,hi=interval.lo,interval.hi
    lower=0 if lo<=0<=hi else min(lo*lo,hi*hi)
    upper=max(lo*lo,hi*hi)
    return iv.Interval(lower//iv.SCALE,iv.ceiling(upper,iv.SCALE))


def kernel_at_energy(m,tau,index,energy):
    """Direct squared Fourier amplitude, independent of characteristic sums."""
    lo,hi=energy;s,c=k.trig_interval(tau*lo,tau*hi)
    phase_s,phase_c=k.phase(m,index)
    step_r=c*phase_c+s*phase_s;step_i=c*phase_s-s*phase_c
    power_r,power_i=iv.Interval.point(1),iv.Interval.point(0)
    total_r,total_i=iv.Interval.point(0),iv.Interval.point(0)
    for _ in range(m):
        total_r+=power_r;total_i+=power_i
        power_r,power_i=power_r*step_r-power_i*step_i,power_r*step_i+power_i*step_r
    probability=(square(total_r)+square(total_i))/(m*m)
    return iv.Interval(max(0,probability.lo),min(iv.SCALE,probability.hi))


def direct_distribution_and_leakage(m,tau,resolution):
    centers=k.energy_centers(m,tau);rows=load_source()
    probabilities=[iv.Interval.point(0) for _ in range(m)]
    leakage=iv.Interval.point(0);ambiguous=0
    for lo,hi,low,high in rows:
        weight=enclosure(low,high)
        for j,center in enumerate(centers):
            probability=weight*kernel_at_energy(m,tau,j,(lo,hi))
            probabilities[j]+=probability
            a,b=F(center.lo,iv.SCALE)-hi,F(center.hi,iv.SCALE)-lo
            minimum=0 if a<=0<=b else min(abs(a),abs(b));maximum=max(abs(a),abs(b))
            if minimum>=resolution:leakage+=probability
            elif maximum>=resolution:
                leakage+=iv.Interval(0,probability.hi);ambiguous+=1
    return probabilities,leakage,ambiguous


def outward_decimal(interval,digits=12):
    scale=10**digits
    low=interval.lo*scale//iv.SCALE;high=iv.ceiling(interval.hi*scale,iv.SCALE)
    def render(n):
        sign='-' if n<0 else '';n=abs(n)
        return f'{sign}{n//scale}.{n%scale:0{digits}d}'
    return [render(low),render(high)]


def readout(m,tau=F(1),resolution=F(1,2)):
    m=k.count(m);tau=k.exact(tau,'time step',True)
    if m>64:raise ValueError('retained finite example supports at most64 outcomes')
    probabilities=k.probabilities([characteristic(n*tau) for n in range(1,m)])
    direct,leakage,ambiguous=direct_distribution_and_leakage(m,tau,resolution)
    for first,second in zip(probabilities,direct,strict=True):
        if max(first.lo,second.lo)>min(first.hi,second.hi):
            raise ArithmeticError('independent positive-kernel and Fourier readouts disagree')
    centers=k.energy_centers(m,tau)
    bound=k.real_energy_leakage(m,tau,resolution,F(131,256),support_upper=F(37,20))
    if F(leakage.lo,iv.SCALE)>bound['total']:
        raise ArithmeticError('actual leakage violates uniform kernel bound')
    return {'outcome_count':m,'time_step':str(tau),'largest_sample_time':str((m-1)*tau),
            'resolution':str(resolution),'all_outcomes_retained':True,
            'source_lines_used':78,'direct_positive_and_Fourier_enclosures_agree':True,
            'joint_true_energy_readout_error_probability_interval':leakage.json(),
            'joint_error_probability_decimal':outward_decimal(leakage),
            'ambiguous_resolution_pairs':ambiguous,
            'uniform_error_bounds':{name:str(value) for name,value in bound.items()},
            'rows':[{'index':j,'energy_center_interval':center.json(),
                     'energy_center_decimal':outward_decimal(center),
                     'probability_interval':probability.json(),
                     'probability_decimal':outward_decimal(probability)}
                    for j,(center,probability) in enumerate(zip(centers,probabilities,strict=True))],
            'native_full_energy_POVM_executed':False,
            'finite_shot_records_generated':False,
            'infinite_ground_values_inferred':False,
            'energy_centers_are_relativistic_masses':False}


def native_first_harmonic():
    value=strict_json.load(ROOT/'model/ground_spectral_source/results.json')
    for name,digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('native readout source changed: '+name)
    settings=value['executed_native_readout']['settings']
    bounds=list(map(F,value['complete_native_error_certificate']['complete_native_state_and_probability_errors']))
    real,imag=characteristic(F(1))
    for index,part in enumerate((real,imag)):
        p=F.from_float(settings[index]['all_six_pointer_probabilities'][0])
        observed=2*p-1 if index==0 else 1-2*p
        if observed+2*bounds[index]<F(part.lo,iv.SCALE) or observed-2*bounds[index]>F(part.hi,iv.SCALE):
            raise ArithmeticError('native first harmonic inconsistent with certified source')
    return {'previous_complete_native_first_harmonic_agrees':True,
            'harmonic_time':'1','first_harmonic_source_interval':[real.json(),imag.json()],
            'all_higher_harmonics_natively_executed':False,
            'previous_native_error_is_not_zero':True}
