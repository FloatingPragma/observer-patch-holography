#!/usr/bin/env python3
"""Declared polar-analyzer tickets, native source columns and exact sign audit."""
from pathlib import Path
import json
import math
import numpy as np
import run
import outcome_audit as audit
import independent_check as independent

HERE=Path(__file__).resolve().parent

def weights(theta,winding):
    real=[64*(1+.5*math.cos(winding*theta-2*math.pi*e/3))/3 for e in range(3)]
    values=[math.floor(v) for v in real]
    remain=64-sum(values)
    order=sorted(range(3),key=lambda e:(-(real[e]-values[e]),e))
    for e in order[:remain]:values[e]+=1
    assert sum(values)==64 and min(values)>0
    return tuple(values)

def execute():
    saved=json.loads((HERE/'RECEIPT.json').read_text())
    raw=json.loads((HERE/'RECORDS.json').read_text())
    loop=[tuple(w) for w in saved['weight_loop']]
    maps={w:audit.maps_from_records(t)[0] for w,t in zip(loop,raw['source']['native']['loop'])}
    refinement={}
    for dwell in (4,8):
        wa=[w for w in loop for _ in range(dwell)]
        wb=[w for w in loop+loop for _ in range(dwell//2)]
        assert sorted(wa)==sorted(wb)
        result=audit.path_compare([maps[w] for w in wa],[maps[w] for w in wb])
        exact=audit.exact_compare(wa,wb)
        assert abs(exact['visibility']-result['relative_visibility'])<2e-10
        refinement[str(dwell)]={'comparison':result,'exact_check':exact,'passes_original_visibility_threshold':result['relative_visibility']<-.9}
    geometry=[];records={};instruments={};products={}
    for coupling,winding in [('polar',1),('apolar',2)]:
        seq=[];ws=[];table=[];diags=[]
        for step in range(1,97):
            theta=step*2*math.pi/96;w=weights(theta,winding)
            p,r,d=run.instrument(w,'native');seq.append(p);ws.append(w);table.append(r);diags.append(d)
            geometry.append({'coupling':coupling,'step':step,'analyzer_theta':run.num(theta),
                             'analyzer_n':list(map(run.num,(math.cos(theta),math.sin(theta)))),
                             'native_ticket_weights':w})
        assert ws[-1]==(32,16,16)
        reference=maps[(32,16,16)];products[coupling]={}
        for dwell in (1,2):
            matrices=[p for p in seq for _ in range(dwell)]
            result=run.interference(matrices,reference)
            exact=audit.exact_compare([w for w in ws for _ in range(dwell)],[(32,16,16)]*len(matrices))
            assert abs(exact['visibility']-result['relative_interference_visibility'])<2e-10
            # Check each exact product sign without passing eigenspaces to source.
            aa=[independent.rational_matrix(w,'native')[0] for w in ws for _ in range(dwell)]
            rational_a,rational_b=independent.product(aa)
            result['exact_final_second_coordinate_sign']=(rational_b>0)-(rational_b<0)
            result['exact_visibility_checked']=True
            products[coupling][str(dwell)]=result
        records[coupling]=table;instruments[coupling]=diags
    checks={'polar_dwell2_negative':products['polar']['2']['normalized_source_overlap_with_initial']<-.99,
            'apolar_dwell2_positive':products['apolar']['2']['normalized_source_overlap_with_initial']>.99}
    return {'candidate':'A-native-triangle-optical-holonomy-v1 with declared polar apparatus extension',
            'extra_native_calls':2*96*384,'refinement':refinement,'rotor_products':products,
            'rotor_checks':checks,'geometry_records':geometry,'instruments':instruments,
            'native_rotor_implemented':False,'physical_rotation_dictionary_derived':False,
            'row_criterion_passed':False,'scope':'Externally declared polar versus apolar analyzer-to-scheduler map; native source response follows the tickets.'},records

if __name__=='__main__':
    result,records=execute()
    for name,obj in [('ROTOR_RECEIPT.json',result),('ROTOR_RECORDS.json',records)]:
        with(HERE/name).open('x')as f:f.write(run.canonical(obj))
    print(run.canonical({'refinement':{k:v['comparison'] for k,v in result['refinement'].items()},
                        'rotor_products':result['rotor_products'],'rotor_checks':result['rotor_checks']}),end='')
