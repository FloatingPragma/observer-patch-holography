"""Source-blind two-stride interval reader; standard library and no I/O."""
import math

def cosine_band(record, eta, delta):
    if not 0 <= eta < 1 or not 0 <= delta <= .5:
        raise ValueError('invalid explicit certificate/error bound')
    for axis in ('x','y'):
        plus,minus=record[axis+'_plus'],record[axis+'_minus']
        if not all(math.isfinite(x) and 0 <= x <= 1 for x in (plus,minus)):
            raise ValueError('invalid binary probability')
        if abs(plus+minus-1)>1e-12:
            raise ValueError('incomplete binary outcome law')
    real=record['x_plus']-record['x_minus']
    lower=max(-1.,min(1.,(real-2*delta-eta)/(1-eta)))
    upper=max(-1.,min(1.,(real+2*delta+eta)/(1-eta)))
    return [math.acos(upper),math.acos(lower)]

def read(first, later, stride=8, eta=.01, probability_error=1e-5):
    if type(stride) is not int or stride<1:
        raise ValueError('stride must be a positive integer')
    initial=cosine_band(first,eta,probability_error)
    phase=cosine_band(later,eta,probability_error)
    branches=[]
    # Enumerate ALL cosine inverse branches on the entire unsigned chart [0,pi].
    for winding in range(-stride,stride+1):
        for sign in (-1,1):
            a,b=(phase if sign==1 else [-phase[1],-phase[0]])
            lower=max(0.,initial[0],(2*math.pi*winding+a)/stride)
            upper=min(math.pi,initial[1],(2*math.pi*winding+b)/stride)
            if lower<=upper:
                branches.append({'winding':winding,'sign':sign,'interval':[lower,upper]})
    merged=[]
    for lower,upper in sorted(x['interval'] for x in branches):
        if merged and lower<=merged[-1][1]+1e-14:
            merged[-1][1]=max(merged[-1][1],upper)
        else:merged.append([lower,upper])
    return {'first_interval':initial,'stride_phase_interval':phase,
            'all_compatible_branches':branches,'frequency_components':merged,
            'unique_component':len(merged)==1,'stride':stride,
            'eta_external_certificate':eta,'probability_error':probability_error,
            'finite_shot_confidence_claim':False}
