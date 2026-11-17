"""Read oriented marker winding from retained current-decoded positions only."""
import gzip
import json
import math
from pathlib import Path

def xy(slot,m):
    triangle=((1.,0.),(-.5,math.sqrt(3)/2),(-.5,-math.sqrt(3)/2))
    v,k=divmod(slot,m)
    return tuple((1-k/m)*triangle[v][i]+k/m*triangle[(v+1)%3][i] for i in (0,1))

rows=[]
d=json.loads(gzip.decompress((Path(__file__).resolve().parent/'first_run/EVIDENCE.json.gz').read_bytes()))
for case in d['cases']:
    m=case['M'];old={'rot':2*m,'ref':2*m};angles={'rot':[],'ref':[]}
    for r in case['records']:
        if r['kind']!='marker':continue
        arm=r['arm'];positions={h['decoded'] for h in r['histories'] if h['outcome']=='hop'}
        assert len(positions)==1
        new=positions.pop();a=xy(old[arm],m);b=xy(new,m)
        angles[arm].append(math.degrees(math.atan2(a[0]*b[1]-a[1]*b[0],a[0]*b[0]+a[1]*b[1])))
        old[arm]=new
    rot=sum(angles['rot']);ref=sum(angles['ref']);maximum=max(abs(x) for v in angles.values() for x in v)
    assert abs(rot-case['direction']*360*case['turns'])<1e-9 and abs(ref)<1e-9 and maximum<90
    rows.append({'M':m,'turns':case['turns'],'direction':case['direction'],'law':case['law'],
                 'rotating_degrees':rot,'reference_degrees':ref,'maximum_local_degrees':maximum})
print(json.dumps({'reader':'current-decoded Cartesian frame; no supplied half-angle phase','cases':rows,'passed':True},indent=2))
