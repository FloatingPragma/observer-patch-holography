#!/usr/bin/env python3
"""Supplemental common-coin source coupling; frozen separately after review."""
from pathlib import Path
import argparse
import hashlib
import json
from native_endpoint import native

HERE=Path(__file__).resolve().parent


def execute():
    prereg=json.loads((HERE/'COUPLING_PREREG.json').read_text())
    for name,digest in prereg['files'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    seams=json.loads((HERE/'seams.json').read_text())
    records=[]; calls=0
    for a,b in seams:
        for sign in (-1,1):
            starts=[tuple(sign*int(k==p) for k in range(14)) for p in (a,b)]
            for coin in (0,1):
                results=[native(q,(a,b),coin) for q in starts]; calls+=2
                assert results[0][0]==results[1][0]
                assert results[0][1][1]!=results[1][1][1]
                assert starts[0]!=starts[1]
                for q,(out,(_,j)) in zip(starts,results):
                    observed=list(q); observed[a]-=j; observed[b]+=j
                    assert tuple(observed)==out
                    assert sum(out)==sign
                followups=[]
                for second_coin in (0,1):
                    future=[native(out,(a,b),second_coin) for out,_ in results]; calls+=2
                    assert future[0]==future[1]
                    followups.append({'coin':second_coin,'equal_full_native_output':True,
                                      'destination':list(future[0][0]),'current':future[0][1][1]})
                controls={}
                for law in ('hold','swap'):
                    outputs=[]
                    for q in starts:
                        out=list(q)
                        if law=='swap':out[a],out[b]=out[b],out[a]
                        outputs.append(out)
                    assert outputs[0]!=outputs[1]
                    controls[law]=outputs
                records.append({'edge':[a,b],'sign':sign,'coin':coin,
                                'initial_source_symbols':[list(q) for q in starts],
                                'native_common_endpoint':list(results[0][0]),
                                'distinct_first_currents':[x[1][1] for x in results],
                                'second_event':followups,'controls':controls})
    return {'cases':len(records),'literal_native_calls':calls,
            'conditional_on_each_coin':True,'future_source_processes_equal':True,
            'historical_current_records_retained_distinct':True,
            'same_group':True,'records':records}


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--record',action='store_true')
    args=parser.parse_args(); result=execute()
    if args.record:
        with (HERE/'COUPLING_RECEIPT.json').open('x')as h:
            h.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else: assert json.loads((HERE/'COUPLING_RECEIPT.json').read_text())==result
    print(json.dumps({k:v for k,v in result.items() if k!='records'},sort_keys=True))
