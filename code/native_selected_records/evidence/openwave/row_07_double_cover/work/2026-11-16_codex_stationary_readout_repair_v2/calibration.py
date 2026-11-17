"""Bounded current-confirmed routing for the declared source-column census.

The admitted unit starts and ends on target port1. Coin-basis settings are
declared coherent memory preparation, not inferred from fair random records.
"""
from fractions import Fraction as F
from native_source import native
from amplified_native import require

EDGES=((0,1),(1,2),(0,2))


def route(old,wanted,law):
    if old==wanted:return F(1),[]
    edge=tuple(sorted((old,wanted)));histories=[];live=[([],[],F(1))]
    for attempt in range(16):
        following=[]
        for coins,currents,probability in live:
            for coin in (0,1):
                q=tuple(int(x==old)for x in range(3));audit,event=native(q,edge,coin,law)
                j=event[3];decoded=edge[1]if j==1 else edge[0]if j==-1 else old
                require(audit[decoded]==1,'Calibration current decoder')
                cs,js,p=coins+[coin],currents+[j],probability/2
                if j or attempt==15:
                    histories.append({'coins':cs,'currents':js,'probability':str(p),
                        'decoded':decoded,'outcome':'arrived'if j else'timeout'})
                else:following.append((cs,js,p))
        live=following
    require(sum(F(r['probability'])for r in histories)==1,'Calibration route omitted outcomes')
    return sum(F(r['probability'])for r in histories if r['outcome']=='arrived'),histories


def census(law):
    position=1;reach=F(1);routes=[];columns=[]
    for wanted in range(3):
        success,histories=route(position,wanted,law)
        routes.append({'from':position,'to':wanted,'reach_before':str(reach),'histories':histories})
        reach*=success;position=wanted
        for edge in EDGES:
            for coin in (0,1):
                q=tuple(int(i==position)for i in range(3));_,event=native(q,edge,coin,law)
                a,b,c,j=event;decoded=list(q);decoded[a]-=j;decoded[b]+=j
                output=decoded.index(1);feedback=c^int(j!=0)
                _,reverse=native(tuple(decoded),edge,feedback,law)
                restored=list(decoded);restored[a]-=reverse[3];restored[b]+=reverse[3]
                require(tuple(restored)==q,'Calibration probe did not return admitted unit')
                columns.append({'input_position':position,'edge':list(edge),'coin_basis':coin,
                    'forward':event,'decoded_output':output,'reverse_probe':reverse,
                    'decoded_restored':restored.index(1)})
    success,histories=route(position,1,law)
    routes.append({'from':position,'to':1,'reach_before':str(reach),'histories':histories});reach*=success
    return {'law':law,'initial_position':1,'final_position':1,'routes':routes,'columns':columns,
        'success':str(reach),'fixed_column_forward_events':36,'maximum_route_forward_events':64,
        'meaning':'Exact finite column census with declared coin-basis inputs; no statistical norm estimation is claimed.',
        'weights':'Declared rational coherent selector amplitudes weight these current-derived endpoint columns.'}
