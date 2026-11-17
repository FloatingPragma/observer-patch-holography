"""Standalone readout excerpt of filed complete pair outcome laws.
External numerical source certificate eta=.01; no source density is supplied
here. This reproduces interval arithmetic, not native preparation or sampling.
"""
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

import json
SCALE=511*math.sqrt(5)/2
OLD_PAIR=1024.4508487023415
# Case tag, stride-1 X+/Y+, stride-8 X+/Y+; complementary outcomes follow.
DATA=[
    [{'M': 4, 'turns': 1, 'direction': 1, 'law': 'native'}, 0.39504593430543355, 0.9855666106041328, 0.416649503835124, 0.9898543054336295],
    [{'M': 4, 'turns': 1, 'direction': 1, 'law': 'swap'}, 0.3948276698140349, 0.9854960339103889, 0.41572761853962226, 0.9884956881567634],
    [{'M': 4, 'turns': 1, 'direction': -1, 'law': 'native'}, 0.39504593430584206, 0.9855666106042644, 0.41664950383685234, 0.9898543054361724],
    [{'M': 4, 'turns': 1, 'direction': -1, 'law': 'swap'}, 0.3948276698140349, 0.9854960339103891, 0.41572761853962215, 0.9884956881567633],
    [{'M': 4, 'turns': 2, 'direction': 1, 'law': 'native'}, 0.39504592564020335, 0.9855666078032951, 0.41664946718477225, 0.9898542514871652],
    [{'M': 4, 'turns': 2, 'direction': 1, 'law': 'swap'}, 0.39482766981403505, 0.9854960339103889, 0.4157276185396225, 0.9884956881567634],
    [{'M': 4, 'turns': 2, 'direction': -1, 'law': 'native'}, 0.39504592564121527, 0.9855666078036234, 0.4166494671890525, 0.9898542514934643],
    [{'M': 4, 'turns': 2, 'direction': -1, 'law': 'swap'}, 0.39482766981403483, 0.9854960339103891, 0.41572761853962215, 0.9884956881567633],
    [{'M': 8, 'turns': 1, 'direction': 1, 'law': 'native'}, 0.39504593589659, 0.9855666111184362, 0.4166495105650596, 0.9898543153395711],
    [{'M': 8, 'turns': 1, 'direction': 1, 'law': 'swap'}, 0.39482766981403494, 0.9854960339103889, 0.41572761853962215, 0.988495688156763],
    [{'M': 8, 'turns': 1, 'direction': -1, 'law': 'native'}, 0.3950459358965901, 0.9855666111184369, 0.4166495105650599, 0.9898543153395715],
    [{'M': 8, 'turns': 1, 'direction': -1, 'law': 'swap'}, 0.3948276698140349, 0.9854960339103889, 0.41572761853962187, 0.9884956881567628],
    [{'M': 8, 'turns': 2, 'direction': 1, 'law': 'native'}, 0.39504593572957997, 0.9855666110644539, 0.4166495098586763, 0.9898543142998294],
    [{'M': 8, 'turns': 2, 'direction': 1, 'law': 'swap'}, 0.3948276698140349, 0.9854960339103893, 0.4157276185396222, 0.9884956881567636],
    [{'M': 8, 'turns': 2, 'direction': -1, 'law': 'native'}, 0.39504593572958047, 0.9855666110644536, 0.41664950985867666, 0.9898543142998297],
    [{'M': 8, 'turns': 2, 'direction': -1, 'law': 'swap'}, 0.39482766981403483, 0.985496033910389, 0.4157276185396217, 0.9884956881567624],
 ]
def complete(x,y):
    return {'x_plus':x,'x_minus':1-x,'y_plus':y,'y_minus':1-y}
rows=[]
for case,x1,y1,x8,y8 in DATA:
    # 2e-5 envelopes measurement error <=1e-5 plus a reader allowance1e-5.
    result=read(complete(x1,y1),complete(x8,y8),probability_error=2e-5)
    assert result['unique_component']
    energy=[SCALE*x for x in result['frequency_components'][0]]
    assert energy[1]<OLD_PAIR
    rows.append({'case':case,'reader':result,'energy_interval_keV':energy})
print(json.dumps({'rows':rows,'source_certificate_verified_here':False,
    'eta_external_hypothesis':.01,'deterministic_probability_error':1e-5,
    'no_finite_shot_confidence_claim':True,'passed':True},indent=2))
