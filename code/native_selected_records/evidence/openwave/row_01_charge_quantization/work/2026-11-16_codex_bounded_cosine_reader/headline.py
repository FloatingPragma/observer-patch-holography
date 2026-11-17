"""Standalone stdlib-only bounded reader of filed finite outcome laws.
The fixed eta=0.01 certificate is an explicit hypothesis here; this code does
not reconstruct the source density that certifies it. No finite-shot claim.
"""
import math,json
ETA=.01
RECORDS=[[{'M': 4, 'turns': 1, 'direction': 1, 'law': 'native'}, 0.8122515617667198, 0.8883630411714384], [{'M': 4, 'turns': 1, 'direction': 1, 'law': 'swap'}, 0.8121546617228033, 0.8884270545457857], [{'M': 4, 'turns': 1, 'direction': -1, 'law': 'native'}, 0.8122515617669321, 0.8883630411712979], [{'M': 4, 'turns': 1, 'direction': -1, 'law': 'swap'}, 0.8121546617228034, 0.8884270545457857], [{'M': 4, 'turns': 2, 'direction': 1, 'law': 'native'}, 0.8122515579071918, 0.8883630437210907], [{'M': 4, 'turns': 2, 'direction': 1, 'law': 'swap'}, 0.8121546617228034, 0.8884270545457859], [{'M': 4, 'turns': 2, 'direction': -1, 'law': 'native'}, 0.8122515579077179, 0.888363043720743], [{'M': 4, 'turns': 2, 'direction': -1, 'law': 'swap'}, 0.8121546617228034, 0.8884270545457857], [{'M': 8, 'turns': 1, 'direction': 1, 'law': 'native'}, 0.8122515624754382, 0.8883630407032501], [{'M': 8, 'turns': 1, 'direction': 1, 'law': 'swap'}, 0.8121546617228033, 0.8884270545457857], [{'M': 8, 'turns': 1, 'direction': -1, 'law': 'native'}, 0.8122515624754383, 0.8883630407032501], [{'M': 8, 'turns': 1, 'direction': -1, 'law': 'swap'}, 0.8121546617228033, 0.8884270545457859], [{'M': 8, 'turns': 2, 'direction': 1, 'law': 'native'}, 0.8122515624010518, 0.8883630407523908], [{'M': 8, 'turns': 2, 'direction': 1, 'law': 'swap'}, 0.8121546617228034, 0.8884270545457857], [{'M': 8, 'turns': 2, 'direction': -1, 'law': 'native'}, 0.8122515624010518, 0.8883630407523908], [{'M': 8, 'turns': 2, 'direction': -1, 'law': 'swap'}, 0.8121546617228033, 0.8884270545457857]]

def read(x_plus,y_plus,delta=0):
    assert 0<=x_plus<=1 and 0<=y_plus<=1
    r=2*x_plus-1
    lower=max(-1,min(1,(r-2*delta-ETA)/(1-ETA)))
    upper=max(-1,min(1,(r+2*delta+ETA)/(1-ETA)))
    return [math.acos(upper),math.acos(lower)]

rows=[]
for case,x,y in RECORDS:
    interval=read(x,y)
    assert interval[1]-interval[0]<.06
    rows.append({'case':case,'binary_outcomes':{'x':[x,1-x],'y':[y,1-y]},
                 'frequency_interval':interval})
# The inequality covers arbitrary band frequencies, populations and other bands.
# These are mathematical sensitivity checks, not native experiment outcomes.
for omega in (0,.001,.3,1.,1.6,3.14,math.pi):
    for u in (0,.001,ETA):
        for other_cosine in (-1,0,1):
            r=(1-u)*math.cos(omega)+u*other_cosine
            lo,hi=read((1+r)/2,.5)
            assert lo-1e-12<=omega<=hi+1e-12
print(json.dumps({'hypothesis_eta':ETA,'chart':[0,math.pi],
    'rows':rows,'mathematical_boundary_checks':63,
    'source_certificate_verified_here':False,'finite_shot_claim':False,
    'passed':True},indent=2))
