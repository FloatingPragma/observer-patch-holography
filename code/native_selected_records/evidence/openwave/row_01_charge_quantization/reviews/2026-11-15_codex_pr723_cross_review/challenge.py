#!/usr/bin/env python3
"""Adversarial current-reader challenge; independent floor/ceiling source."""
import ast
import hashlib
import json
from pathlib import Path
import random

HERE=Path(__file__).resolve().parent
PACKET=HERE.parents[1]/'work/2026-11-15_codex_pr723_unit_selection'
BOUND='35b67994df0c89d8478624815a95585293b99249da46b4dc90ba5fd27b72d5ae'

def endpoint(x,edge,coin,control=False):
    a,b=edge;y=list(x)
    if control:
        if coin:y[a],y[b]=x[b],x[a]
    else:
        s=x[a]+x[b];lo=s//2;hi=-((-s)//2)
        y[a],y[b]=(lo,hi)if coin==0 else(hi,lo)
    return y,x[a]-y[a]

def execute():
    assert hashlib.sha256((PACKET/'MANIFEST.sha256').read_bytes()).hexdigest()==BOUND
    edges=[tuple(e)for e in json.loads((PACKET/'seams.json').read_text())]
    seen={0};queue=[0];tree=[]
    for p in queue:
        for a,b in edges:
            q=b if a==p else a if b==p else None
            if q is not None and q not in seen:seen.add(q);queue.append(q);tree.append((p,q))
    text=(PACKET/'source.py').read_text();syntax=ast.parse(text)
    targets=[n for n in syntax.body if isinstance(n,ast.FunctionDef) and n.name in('decode_probe','decode_body_path')]
    assert {n.name: [a.arg for a in n.args.args]for n in targets}=={
        'decode_probe':['tape','cap'],'decode_body_path':['body_tape','terminal_relative']}
    namespace={'N':12,'TREE':tuple(tree)}
    exec(compile(ast.Module(body=targets,type_ignores=[]),str(PACKET/'source.py'),'exec'),namespace)
    decoder=namespace['decode_probe'];body_decoder=namespace['decode_body_path']
    calls=0
    def probe(initial,seed,cap=32,forced=None):
        nonlocal calls
        x=initial[:];rng=random.Random(seed);tape=[]
        for edge in tree:
            first=None
            for attempt in range(cap):
                coin=rng.randrange(2)if forced is None else forced
                x,j=endpoint(x,edge,coin);calls+=1;tape.append([*edge,coin,j])
                if first is None:first=coin
                elif first!=coin:break
            else:return tape
        return tape
    rng=random.Random(7230701);initials=[[rng.randint(-10000,10000)for _ in range(12)]for _ in range(128)]
    shift_equal=decoded=0
    for i,x in enumerate(initials):
        tape=probe(x,17031+i);other=probe([v+10**9 for v in x],17031+i)
        assert tape==other;shift_equal+=1
        value=decoder(tape);assert value==[v-x[0]for v in x];decoded+=1
    paths=0
    for initial in initials[:16]:
        for control in(False,True):
            x=initial[:];truth=[x[:]];tape=[]
            for _ in range(64):
                edge=rng.choice(edges);coin=rng.randrange(2);x,j=endpoint(x,edge,coin,control)
                calls+=1;tape.append([*edge,coin,j]);truth.append(x[:])
            terminal=probe(x,10031+paths);relative=decoder(terminal)
            recovered=body_decoder(tape,relative)
            assert [[v-y[0]for v in y]for y in truth]==[[v-y[0]for v in y]for y in recovered]
            if control:assert len({max(y)-min(y)for y in recovered})==1
            paths+=1
    timeout=probe([2]+[0]*11,1,cap=2,forced=0);assert decoder(timeout,2)is None
    zero=probe([0]*12,181);two=probe([2]*12,181)
    assert zero==two and decoder(zero)==decoder(two)==[0]*12
    return {'bound_manifest_sha256':BOUND,'signed_large_input_reconstructions':decoded,
            'exact_shift_invariant_tape_pairs':shift_equal,'native_and_control_backward_paths':paths,
            'independent_endpoint_calls':calls,'forced_timeout_preserved':True,
            'uniform_zero_two_operationally_indistinguishable':True,
            'decoder_has_no_total_or_preparation_parameter':True,'all_challenges_passed':True,
            'same_group':True,'G6_awarded':False}

if __name__=='__main__':print(json.dumps(execute(),sort_keys=True,indent=2))
