"""One private integer state continues through admission, calibration and return."""
import random


def step(q,edge,coin,law='native'):
    a,b=edge;d=q[a]-q[b]
    current=(d-(2*coin-1)*(abs(d)%2))//2 if law=='native' else d
    result=list(q);result[a]-=current;result[b]+=current
    return result,[a,b,coin,current]


def source(initial, edges, admission_seed, calibration_seed, probe_law='native', forced=None):
    if probe_law not in ('native','always_swap'):raise ValueError('probe law')
    state=list(initial);events=[];allowed={tuple(sorted(e)) for e in edges}
    admission_coins=random.Random(admission_seed);calibration_coins=random.Random(calibration_seed)
    counts={'native':0,'always_swap':0};forced=forced or {}

    def request(edge,stage):
        nonlocal state
        if tuple(sorted(edge)) not in allowed:raise ValueError('non-source edge')
        calibration=stage.startswith('calibration_')
        rng=calibration_coins if calibration else admission_coins
        coin=forced[stage] if stage in forced else rng.randrange(2)
        law=probe_law if stage=='calibration_probe' else 'native'
        state,record=step(state,edge,coin,law);counts[law]+=1
        events.append([stage,list(record)])
        return list(record)

    def audit():
        return {'final':list(state),'events':[[stage,list(record)] for stage,record in events],'literal_endpoint_calls':dict(counts)}
    return request,audit
