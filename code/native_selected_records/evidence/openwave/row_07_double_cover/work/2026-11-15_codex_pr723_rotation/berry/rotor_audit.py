#!/usr/bin/env python3
"""Replay frozen rotor tickets; check formula independently at high precision."""
from decimal import Decimal as D, localcontext
from pathlib import Path
import json
import math
import numpy as np
import run
import outcome_audit as audit

HERE=Path(__file__).resolve().parent
PI=D('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825')

def decimal_cos(x):
    with localcontext() as ctx:
        ctx.prec=85
        while x>PI:x-=2*PI
        while x<-PI:x+=2*PI
        term=D(1);total=term
        for n in range(1,300):
            term*=-x*x/D((2*n-1)*(2*n));total+=term
            if abs(term)<D('1e-80'):return +total
        raise AssertionError('cos series did not converge')

def check_tickets(row):
    with localcontext() as ctx:
        ctx.prec=80
        n=1 if row['coupling']=='polar' else 2
        theta=D(row['step'])*2*PI/96
        values=[D(64)*(1+decimal_cos(n*theta-2*PI*e/3)/2)/3 for e in range(3)]
        # Mathematical integers can be perturbed by truncated pi/Taylor error.
        for k,v in enumerate(values):
            nearest=v.to_integral_value()
            if abs(v-nearest)<D('1e-65'):values[k]=nearest
        floors=[int(v//1) for v in values]
        fractions=[v-f for v,f in zip(values,floors)]
        weights=list(row['native_ticket_weights'])
        added=[w-f for w,f in zip(weights,floors)]
        assert sum(weights)==64 and all(x in(0,1)for x in added)
        chosen=[i for i,a in enumerate(added)if a]
        unchosen=[i for i,a in enumerate(added)if not a]
        ties=[]
        for i in chosen:
            for j in unchosen:
                assert fractions[i]>=fractions[j]-D('1e-60')
                if abs(fractions[i]-fractions[j])<D('1e-60'):
                    ties.append([i,j])
        return ties

def execute(reexecute_source=False):
    receipt=json.loads((HERE/'ROTOR_RECEIPT.json').read_text())
    records=json.loads((HERE/'ROTOR_RECORDS.json').read_text())
    ties=[];sequences={};maps={};counts=0
    for coupling in('polar','apolar'):
        rows=[r for r in receipt['geometry_records']if r['coupling']==coupling]
        sequence=[]
        for row,table in zip(rows,records[coupling]):
            tie=check_tickets(row)
            if tie:ties.append({'coupling':coupling,'step':row['step'],'recorded_weights':row['native_ticket_weights'],'equally_valid_remainder_choices':tie})
            w=tuple(row['native_ticket_weights']);sequence.append(w)
            if reexecute_source:
                _,fresh,_=run.instrument(w,'native');assert fresh==table;counts+=len(fresh)
            maps[w]=audit.maps_from_records(table)
        sequences[coupling]=sequence
    original=(32,16,16)
    if original not in maps:raise AssertionError('source loop did not close')
    polar=sequences['polar'];apolar=sequences['apolar']
    two=polar+polar
    result=audit.path_compare([maps[w][0]for w in two],[maps[original][0]]*len(two))
    exact=audit.exact_compare(two,[original]*len(two))
    assert abs(exact['visibility']-result['relative_visibility'])<2e-10
    checks={'same_polar_program_twice_positive':result['loop_overlap_initial']>.99,
            'same_polar_program_twice_visibility_positive':result['relative_visibility']>.9}
    outcome={'polar96':audit.outcomes(polar,[original]*96,maps),
             'apolar96':audit.outcomes(apolar,[original]*96,maps),
             'polar192':audit.outcomes([w for w in polar for _ in range(2)],[original]*192,maps)}
    return {'portable_input':'recorded integer weights in ROTOR_RECEIPT.json; no libm weight regeneration',
            'decimal_precision':80,'cos_tolerance': '1e-60 (far above truncated-series/pi error)',
            'tie_cases':ties,'source_calls_reexecuted':counts,
            'two_turns':result,'two_turn_exact_check':exact,'checks':checks},outcome

if __name__=='__main__':
    result,outcomes=execute()
    for name,obj in [('ROTOR_AUDIT.json',result),('ROTOR_OUTCOMES.json',outcomes)]:
        with(HERE/name).open('x')as f:f.write(run.canonical(obj))
    print(run.canonical({'tie_cases':result['tie_cases'],'two_turns':result['two_turns'],'checks':result['checks']}),end='')
