"""Opposing read-only reconstruction of the root row1 amendments."""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
import random
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
ROW=ROOT/'openwave-tuning/next/OBSTRUCTION_REPAIR/row1'
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROW))
import adaptive_phase as ap
import neutral_record as nr


def step(n,e,coin):
    a,b=nr.EDGES[e]
    total=n[a]+n[b]
    left=(total+coin)//2 if total%2 else total//2
    # For odd total native coin1 selects ceil; for even both select half.
    out=list(n);out[a]=left;out[b]=total-left
    return tuple(out),n[a]-left,n[a]-n[b]


def main():
    geom=json.loads((ROOT/'openwave-tuning/row_completion/ROW_01/GEOMETRY.json').read_text())
    rng=random.Random(112902)
    phase=[]
    # Analytic margin: .62+.46*pi < 2*pi follows even from pi>1.
    assert F(31,50)<F(77,50)
    for k in range(60):
        current=[rng.randrange(-512,513) for _ in range(30)]
        region=tuple(i for i in range(12) if rng.randrange(2))
        exact=sum(s*current[e] for v in region for e,s in geom['dual_port_boundaries'][v])
        for deformed in [False,True]:
            r=ap.read_phase(geom,current,region,deformed=deformed,sample_budget=200000)
            assert r['status']=='read'
            assert abs(r['phase_degree']-exact)<2e-9
            assert r['max_principal_increment']<math.pi/4
            phase.append(abs(r['phase_degree']-exact))
    assert ap.read_phase(geom,[10**80]*30,(0,),sample_budget=1000)['status']=='resource_stop'
    count=0
    for filling in [1,19]:
        for e in range(30):
            for sign in [-1,1]:
                loads=[filling]*12;a,b=nr.EDGES[e];loads[a]-=sign;loads[b]+=sign
                records=[dict(kind='prepare_pair',seam=e,current=sign)]
                order=list(range(30));rng.shuffle(order)
                schedules=[('C',[(rng.randrange(30),rng.randrange(2)) for _ in range(47)]),
                           ('E',[(i,rng.randrange(2)) for i in order])]
                for kind,seq in schedules:
                    for edge,c in seq:
                        loads,j,d=step(loads,edge,c)
                        rec=dict(kind=kind,seam=edge,coin=c,current=j)
                        if kind=='E':rec['contrast']=d
                        records.append(rec)
                tape=dict(records=records,sweep_order=order)
                for region in [[],list(range(12)),[e%12],list(range(0,12,2))]:
                    obs=nr.observe_source(filling,tape,region)
                    assert obs['charge']==sum(loads[v]-filling for v in region)
                    assert obs['native_attempts']==77
                    assert obs['preparation_interventions']==1
                    # Malformed received tapes must fail before pointer handling.
                    corrupt=deepcopy(tape);corrupt['records'][9]['current']+=2
                    try:nr.observe_source(filling,corrupt,region)
                    except ValueError:pass
                    else:raise AssertionError('mutated current admitted')
                    hidden=deepcopy(tape);hidden['final_loads']=[-999]*12
                    assert nr.observe_source(filling,hidden,region)==obs
                    count+=1
    examples=json.loads((ROW/'NEUTRAL_RECORD_EXAMPLES.json').read_text())
    effects=0
    for ex in examples[::7]:
        got=nr.admit(ex['filling'],ex['source_records'],ex['region'],ex['pointer'])
        assert got['admitted']
        hidden=deepcopy(ex['pointer']);hidden['final_loads']=None
        assert nr.admit(ex['filling'],ex['source_records'],ex['region'],hidden)==got
        for x,weight in enumerate(map(F,got['fine_record_effect'])):
            assert weight==(F(1,2**got['native_attempts']) if x in got['compatible_pointer_initials'] else 0)
        effects+=1
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
          [ROW/'adaptive_phase.py',ROW/'neutral_record.py',ROW/'NEUTRAL_RECORD_EXAMPLES.json',
           ROOT/'openwave-tuning/row_completion/ROW_01/GEOMETRY.json']}
    out=dict(schema='oph.obstruction_repair.row1.opposing.hadrons.v1',
        adaptive_increment_proof_pass=True,phase_reads=len(phase),phase_max_error=max(phase),
        independent_neutral_source_reads=count,retained_pointer_effect_cases=effects,
        sample_budget_stop_pass=True,malformed_source_rejected=True,hidden_source_and_pointer_ignored=True,
        component_review_pass=True,physical_charge_or_platform_validation=False,
        scope='The constructive decoder and known neutral-pair source protocol pass these independent checks. Faithful event records, scalar preparation filling, native integer alphabet and the spatial field dictionary remain declared. No source-independent physical charge quantization or coupling derivation follows.',
        pins=pins)
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='pins'},indent=2))


if __name__=='__main__':main()
