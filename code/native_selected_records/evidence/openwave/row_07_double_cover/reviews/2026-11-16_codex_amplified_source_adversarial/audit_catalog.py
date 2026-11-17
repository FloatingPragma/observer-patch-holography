"""Independent complete source-catalog audit; no producer-module imports."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from circuits import (EDGES,TRI,array,near,need,native_column,reversal,hadamard,
                      weighted_source,contrast,amplified)
from exact_phase import certify_source,certify_five_selector

HERE=Path(__file__).resolve().parent


def characteristic(counts):
    c=sum(counts[0]);tr=sum(counts[i][i]for i in range(3))-c
    sq=sum(counts[i][j]*counts[j][i]for i in range(3)for j in range(3))-c*c
    return {'trace':str(tr),'square':str(sq),'discriminant':str(2*sq-tr*tr)}


def inspect(compiled,path):
    source=compiled['source_kind'];law=compiled['source_law'];records=compiled['source_records']
    need(law in ('native','swap'),'Undeclared endpoint law')
    signs={r['charge_sign']for r in records};need(len(signs)==1,'Mixed charge sectors')
    sign=next(iter(signs));need(sign in(-1,1),'Invalid source charge')
    phase=compiled['phase_record'];preboost=compiled['preboost']
    if source=='star':
        slot,m=compiled['slot'],compiled['M'];v,k=divmod(slot%(3*m),m)
        edges=[e for e in TRI if v in e]*(m-k)+[e for e in TRI if (v+1)%3 in e]*k
        mapping=(1,2,0)
        edges=[tuple(sorted(mapping[x]for x in e))for e in edges]
        need(compiled['coordinate_system']=='energy_target_order'and compiled['spin_marker_to_energy_vertex']==list(mapping),
             'Incorrect physical source-coordinate map')
        width=4*m;perm,counts,columns=reversal(edges,width,law,sign);h=hadamard(width)
        def primitive(x):
            blocks=np.einsum('ij,pjk->pik',h,x.reshape(3,width,-1),optimize=True).reshape(3*width,-1)
            blocks=blocks[perm].reshape(3,width,-1)
            return np.einsum('ij,pjk->pik',h,blocks,optimize=True).reshape(3*width,-1)
        if 'swap_fixed_native_calibration' in path:
            need(law=='swap','Fixed-native control does not use swap source')
            _,calibration_counts,_=reversal(edges,width,'native',sign)
        else:calibration_counts=counts
        selector=None
    else:
        need(source in('finite_weighted','contrast_lcu'),'Unknown source construction')
        selector_record=compiled['selector_phase']
        selector=certify_five_selector(selector_record)
        primitive,counts,columns,active,inactive=weighted_source(
            law,sign,selector_record['real'],selector_record['imag'])
        need(F(compiled['selector_active_weight'])==active and
             F(compiled['selector_inactive_weight'])==inactive,'Finite selector weights differ')
        if source=='contrast_lcu':primitive,counts=contrast(primitive,counts)
        width=len(primitive)//3;calibration_counts=counts
        near(primitive,array(compiled['U']),'Complete weighted/contrast dilation')
    need(compiled['width']==width,'Source workspace dimension differs')
    need([[F(x)for x in row]for row in compiled['counts']]==counts,'Current/choice source coefficients differ')
    near(array(compiled['P']),np.array(counts,dtype=float),'Source mean block')
    # Each ticket/coin/position has one explicit current record; no dropped or
    # duplicated source column is hidden by dict insertion or a shorter zip.
    expected={(p,t,c):(edge,j,out,cc)for p,t,edge,c,j,out,cc in columns}
    keys=[(r['input'],r['ticket'],r['forward'][2])for r in records]
    need(len(keys)==len(expected) and set(keys)==set(expected),'Incomplete/duplicate current catalog')
    for r,key in zip(records,keys):
        edge,j,out,cc=expected[key];p,t,c=key
        need(r['forward']==[edge[0],edge[1],c,j] and r['output']==out and
             r['feedback_coin']==cc,'Native forward column differs')
        back,j2,cc2=native_column(out,edge,cc,law,sign)
        need(back==p and j2==-j,'Literal reversal inverse identity')
        need(r['reverse']==[edge[0],edge[1],cc,j2] and r['after_two_forward_position']==p and
             r['after_inverse_position']==out,'V,V,V-adjoint microtrace differs')
        need(r['final_current_workspace']==0 and r['inverse_scratch_final']==0 and
             r['source_accesses']=={'forward':2,'inverse':1},'Current workspace/access accounting differs')
    need(phase['source_characteristic']==characteristic(calibration_counts),'Calibration characteristic differs')
    certified=certify_source(calibration_counts,phase,preboost)
    f,g,k=amplified(primitive,width,counts,F(phase['real']),F(phase['imag']),preboost)
    saved_f,saved_g=array(compiled['forward']),array(compiled['inverse'])
    need(saved_f.shape==f.shape and saved_g.shape==g.shape,'Missing full instrument flags')
    error=max(near(f,saved_f,'Complete forward instrument'),near(g,saved_g,'Complete adjoint instrument'),
              near(k,array(compiled['K']),'Actual adjoint-pair success block'))
    count=18 if preboost else 6
    need(compiled['native_reversals_per_filter']==count and
         compiled['source_forward_calls_per_filter']==3*count//2 and
         compiled['source_inverse_calls_per_filter']==3*count//2,'Source/adjoint query count differs')
    return {'path':path,'kind':source,'law':law,'sign':sign,'columns':len(records),
            'flags_each_orientation':width,'maximum_matrix_defect':error,
            'exact_phase_certificate':certified,'selector_certificate':selector}


def walk(value,path='$'):
    if isinstance(value,dict):
        if all(k in value for k in ('forward','inverse','phase_record','source_kind','counts')):
            yield path,value
        else:
            for key,item in value.items():yield from walk(item,path+'.'+str(key))
    elif isinstance(value,list):
        for i,item in enumerate(value):yield from walk(item,path+'['+str(i)+']')


def calculate(paths):
    rows=[];inputs={}
    for path in paths:
        raw=path.read_bytes();inputs[path.name]=hashlib.sha256(raw).hexdigest()
        data=json.loads(gzip.decompress(raw))if path.suffix=='.gz'else json.loads(raw)
        for name,compiled in walk(data,path.name):rows.append(inspect(compiled,name))
    need(rows,'No complete source catalogs found')
    return {'passed':True,'same_group':True,'external_G6':False,'no_producer_imports':True,
            'input_hashes':inputs,'compiled_catalogs':len(rows),'current_columns':sum(r['columns']for r in rows),
            'full_flags':sum(2*r['flags_each_orientation']for r in rows),
            'maximum_matrix_defect':max(r['maximum_matrix_defect']for r in rows),'rows':rows}


def replay_form(value):
    if isinstance(value,dict):
        output={}
        for key,item in value.items():
            if key=='maximum_matrix_defect':
                need(type(item)in(int,float)and 0<=item<2e-10,'Saved or fresh matrix defect failed')
                output[key]=0
            else:output[key]=replay_form(item)
        return output
    if isinstance(value,list):return [replay_form(item)for item in value]
    return value


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs',nargs='+',type=Path);parser.add_argument('--record-to',type=Path)
    parser.add_argument('--verify',type=Path,default=HERE/'CATALOG_AUDIT.json')
    args=parser.parse_args()
    if args.record_to:need(not args.record_to.exists(),'First-record destination exists')
    start=time.perf_counter();result=calculate(args.inputs)
    if args.record_to:
        with args.record_to.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:
        saved=json.loads(args.verify.read_text())
        # Verification results include arithmetic defects; compare these with a
        # fixed tiny tolerance, but no source coefficient or gate is omitted.
        need(replay_form(result)==replay_form(saved),'Saved complete catalog verdict differs')
    print(json.dumps({'passed':True,'catalogs':result['compiled_catalogs'],'columns':result['current_columns'],
                     'maximum_defect':result['maximum_matrix_defect'],'seconds':time.perf_counter()-start}))
