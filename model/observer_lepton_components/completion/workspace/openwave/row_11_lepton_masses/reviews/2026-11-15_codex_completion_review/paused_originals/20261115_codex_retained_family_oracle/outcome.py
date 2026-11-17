#!/usr/bin/env python3
"""Read any prepared-input coordinate from the exact receipt, with no source calls."""
from pathlib import Path
import argparse,json
from fractions import Fraction as F
from headline import K,ZERO
HERE=Path(__file__).resolve().parent
def field(row):return K(tuple(map(F,row)))
def sparse(rows):return {int(i):field(v) for i,v in rows}

def read(r,source,count_x,count_y,obj,f_in,f_out,q_in,q_out):
    assert 0<=source<60 and 0<=count_x<60 and 0<=count_y<60
    assert 0<=obj<144 and f_in in range(3) and f_out in range(3)
    assert q_in in(0,1) and q_out in(0,1)
    fam=r['retained_family_query'];arm=r['whole_arm_native_query']
    prep=sparse(fam['prep_numerator']).get(obj,ZERO)
    initial=prep if f_in==f_out else ZERO
    if arm['input_states'][source]['tagged_port']==7:
        amplitude=initial/360 if q_out==1-q_in else ZERO
    else:
        count=r['retained_count_source']['retained_count_vectors']['xor']
        j=60*count_x+count_y;projected=ZERO
        for b in (0,1):
            h=F(count['q'+str(b)+'_numerator'][j])
            y=sparse(fam['query0_family_factors'][str(b)][f_in]).get(4*obj+f_out,ZERO)
            projected+=h*y
        if q_in==0:amplitude=(projected if q_out==0 else initial-projected)/360
        else:amplitude=(initial-projected if q_out==0 else -projected)/360
    return {'amplitude_coefficients':amplitude.data(),
        'original_probability_coefficients':(amplitude.conj()*amplitude).data(),
        'field_basis':['1','sqrt(2)','i','i sqrt(2)'],
        'physical_source_base_unchanged':source,'physical_family_output':f_out,
        'neutral_processor_output':3,'current_raw_output':0,'chi_work_output':0,
        'scope':'Prepared input only; exact receipt query, no source calls or new gate.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for n,d in (('source',0),('count-x',0),('count-y',0),('object',0),
                ('family-in',0),('family-out',0),('query-in',0),('query-out',0)):
        ap.add_argument('--'+n,type=int,default=d)
    a=ap.parse_args();r=json.loads((HERE/'RECEIPT.json').read_text())
    out=read(r,a.source,a.count_x,a.count_y,a.object,a.family_in,a.family_out,a.query_in,a.query_out)
    print(json.dumps(out,indent=2,sort_keys=True))
