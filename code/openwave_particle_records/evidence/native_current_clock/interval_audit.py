#!/usr/bin/env python3
"""Independent directed-integer target audit; --record never overwrites."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
from model.native_effect_synthesis import intervals as iv,synthesis as syn
from model.native_observer_program import reflection
from model.common_source import slow_geometry

ZERO=iv.exact(0);ONE=iv.exact(1)
CZ=(ZERO,ZERO);CO=(ONE,ZERO)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def cadd(a,b):return iv.add(a[0],b[0]),iv.add(a[1],b[1])
def cmul(a,b):return iv.add(iv.mul(a[0],b[0]),iv.neg(iv.mul(a[1],b[1]))),iv.add(iv.mul(a[0],b[1]),iv.mul(a[1],b[0]))
def root(a):
    assert a[0]>=0
    return isqrt(a[0]*iv.Q),isqrt(a[1]*iv.Q)+1
def divide(a,b):
    assert b[0]>0
    values=[F(x*iv.Q,y) for x in a for y in b]
    lo,hi=min(values),max(values)
    return lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator)


def execute_interval(word,bits):
    size=1<<bits
    value=[[CO if i==j else CZ for j in range(size)] for i in range(size)]
    radius=F()
    for event in word:
        targets=event['targets'];d=1<<len(targets)
        gate=syn.primitive_target(event['kind'],event['angle'])
        output=[]
        for row in range(size):
            local=0;base=row
            for target in targets:
                local=2*local+((row>>target)&1);base&=~(1<<target)
            line=[CZ for _ in range(size)]
            for k in range(d):
                source=base
                for j,target in enumerate(targets):source|=((k>>(len(targets)-j-1))&1)<<target
                for col in range(size):line[col]=cadd(line[col],cmul(gate[local][k],value[source][col]))
            output.append(line)
        # Every exact binary-angle logical primitive is unitary. Therefore
        # its action preserves the previous Frobenius error radius. Recenter
        # the new interval product to an exact dyadic matrix and add only
        # the new enclosure/rounding radius, avoiding interval wrapping.
        center=[];fresh_squared=0
        for row in output:
            center_row=[]
            for entry in row:
                centered=[]
                for component in entry:
                    mid=(component[0]+component[1])//2
                    fresh_squared+=max(mid-component[0],component[1]-mid)**2
                    centered.append((mid,mid))
                center_row.append(tuple(centered))
            center.append(center_row)
        radius+=F(isqrt(fresh_squared)+1,iv.Q)
        value=center
    return value,radius


def error_upper(actual,wanted):
    total=0
    for ar,wr in zip(actual,wanted,strict=True):
        for a,w in zip(ar,wr,strict=True):
            for component in range(2):
                delta=iv.add(a[component],iv.neg(w[component]))
                total+=max(abs(x) for x in delta)**2
    return F(isqrt(total)+1,iv.Q)


def targets():
    sqrt5=root(iv.exact(5));raw=slow_geometry.raw_frame()
    frame=[[iv.add(iv.exact(a),iv.mul(iv.exact(b),sqrt5)) for a,b in row]for row in raw]
    spins=[]
    for a,b in ((0,5),(5,11),(0,11)):
        step=[iv.add(y,iv.neg(x))for x,y in zip(frame[a],frame[b])]
        squared=ZERO
        for x in step:squared=iv.add(squared,iv.mul(x,x))
        norm=root(squared);x,y,z=[divide(v,norm) for v in step]
        spins.append([[(z,ZERO),(x,iv.neg(y))],[(x,y),(iv.neg(z),ZERO)]])
    r3=root(iv.exact(3));one_root3=divide(ONE,r3)
    prep=[[(one_root3,ZERO)]if i<3 else[CZ]for i in range(4)]
    ref=[]
    for i in range(8):
        ref.append([(iv.exact(F(2,3)*int(i in(0,2,4)and j in(0,2,4))-int(i==j)),ZERO)for j in range(8)])
    multiplex=[[CO if i==j else CZ for j in range(8)]for i in range(8)]
    for e,s in enumerate(spins):
        for i in range(2):
            for j in range(2):multiplex[2*e+i][2*e+j]=s[i][j]
    return prep,ref,spins,multiplex


def run():
    words=json.loads((HERE/'SOURCE_WORDS.json').read_text())
    all_words=[words['preparation'],words['reflection'],words['spin_multiplex'],*words['spin_gates']]
    unique=sorted({(x['kind'],x['angle'])for word in all_words for x in word})
    for kind,angle in unique:assert syn.certify_primitive(kind,angle)
    prep,ref,spins,multiplex=targets()
    prepared,rp=execute_interval(words['preparation'],2)
    ep=error_upper([[row[0]]for row in prepared],prep)+rp
    reflected,rr=execute_interval(words['reflection'],3)
    er=error_upper(reflected,ref)+rr
    es=[];rs=[]
    for word,target in zip(words['spin_gates'],spins,strict=True):
        matrix,r=execute_interval(word,1);es.append(error_upper(matrix,target)+r);rs.append(r)
    mux,rm=execute_interval(words['spin_multiplex'],3)
    em=error_upper(mux,multiplex)+rm
    aggregate=ep+16*(er+em)
    pulses=len(words['preparation'])+16*(len(words['reflection'])+len(words['spin_multiplex']))
    bound=reflection.complete_bound(pulses,aggregate)
    old=F(1,10**12)
    checks={'all_full_native_primitives_certified':True,'aggregate_below_frozen_1e_minus9':aggregate<=F(1,10**9),
            'corrected_halfdiamond_below_frozen_2e_minus7':F(bound['complete_halfdiamond_upper'])<=F(2,10**7)}
    return {'schema':'oph.native_current_clock.interval_audit.v1','primitives':len(unique),
            'enclosure_algorithm':'exact dyadic centers and additive Frobenius radius under unitary pulses',
            'accumulated_enclosure_radii':{'preparation':str(rp),'reflection':str(rr),'spin_gates':list(map(str,rs)),'multiplex':str(rm)},
            'primitive_raw_gate':str(syn.RAW_GATE),'compiled_piece_pulses_through16':pulses,
            'geometric_error_upper_exact':{'preparation':str(ep),'reflection':str(er),'spin_gates':list(map(str,es)),
                                          'spin_multiplex':str(em),'aggregate16':str(aggregate)},
            'geometric_error_upper_decimal':{'preparation':float(ep),'reflection':float(er),'spin_gates':list(map(float,es)),
                                            'spin_multiplex':float(em),'aggregate16':float(aggregate)},
            'original_aggregate_allowance':str(old),'original_allowance_covers_certified_upper':old>=aggregate,
            'corrected_compiled_piece_bound':bound,'checks':checks,'all_checks_pass':all(checks.values()),
            'source_experiment_unchanged':True,'larger_current_oracle_M6_compiled':False,
            'physical_hardware_precision_authenticated':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args()
    bindings=json.loads((HERE/'INTERVAL_PREREG.json').read_text())['files']
    for name,h in bindings.items():assert sha(ROOT/name)==h,name
    for name,h in json.loads((HERE/'SOURCE_PINS.json').read_text()).items():assert sha(ROOT/name)==h,name
    out=HERE/'INTERVAL_AUDIT.json'
    if args.record and out.exists():raise FileExistsError('Interval evidence already exists')
    start=time.perf_counter();result=run()
    if args.record:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else:assert result==json.loads(out.read_text())
    print(json.dumps({'all_checks_pass':result['all_checks_pass'],'primitives':result['primitives'],
                      'geometric_error_upper_decimal':result['geometric_error_upper_decimal'],
                      'original_allowance_covers_certified_upper':result['original_allowance_covers_certified_upper'],
                      'corrected_halfdiamond':float(F(result['corrected_compiled_piece_bound']['complete_halfdiamond_upper'])),
                      'elapsed_seconds':time.perf_counter()-start},indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
