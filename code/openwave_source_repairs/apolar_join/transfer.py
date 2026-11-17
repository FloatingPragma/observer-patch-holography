#!/usr/bin/env python3
"""Exact complete-instrument bound, independently recomputed from source metadata."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,time
import numpy as np
HERE=Path(__file__).resolve().parent
def sha(data):return hashlib.sha256(data).hexdigest()
def read(name):return json.loads((HERE/name).read_text())
def run():
    reg=read('PREREG_TRANSFER.json')
    for n,h in reg['files'].items():assert sha((HERE/n).read_bytes())==h,n
    sharp=read('inputs/SHARP_FLOW_CERTIFICATE.json')
    firstreg=read('inputs/PREREG_FIRST_CELL.json')
    firstreceipt=read('inputs/first_cell_RECEIPT.json')
    assert sha((HERE/'inputs/PREREG_FIRST_CELL.json').read_bytes())==firstreceipt['prereg_sha256']
    assert sha((HERE/'inputs/source_columns.npz').read_bytes())==firstreg['flow_files']['source_columns.npz']
    src=np.load(HERE/'inputs/source_columns.npz')
    graph=read('inputs/graph.json');assert np.array_equal(src['seams'],np.array(graph['seams']))
    L=np.zeros((240,240),dtype=int)
    for i,j in graph['seams']:L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    assert np.array_equal(L,1260*np.eye(240,dtype=int)-src['counts'])
    catalogue_raw=gzip.decompress((HERE/'inputs/NATIVE_PRIMITIVES.json.gz').read_bytes())
    catalogue=json.loads(catalogue_raw);words=read('inputs/native_WORDS.json')
    flowreceipt=read('inputs/flow_RECEIPT.json')
    bindings={
      'attempts/stage2_002/RECEIPT.json':sha((HERE/'inputs/flow_RECEIPT.json').read_bytes()),
      'attempts/stage2_002/evidence/NATIVE_PRIMITIVES.json':sha(catalogue_raw),
      'attempts/stage2_002/evidence/native_M6_primitives.npz':sha((HERE/'inputs/native_M6_primitives.npz').read_bytes()),
      'attempts/stage2_002/evidence/native_WORDS.json':sha((HERE/'inputs/native_WORDS.json').read_bytes())}
    assert sharp['bindings']==bindings
    for key in ['NATIVE_PRIMITIVES.json','native_M6_primitives.npz','native_WORDS.json']:
        assert flowreceipt['artifacts']['evidence/'+key]==bindings['attempts/stage2_002/evidence/'+key]
    count=len(catalogue['keys']);assert count==196316==len(catalogue['exact_squared_error_bounds'])
    b=F(sharp['full_primitive_operator_error_upper'])
    largest=max(map(F,catalogue['exact_squared_error_bounds']));assert largest<=b*b
    row=[e for e in words['experiments'] if str(e['time'])=='4'];assert len(row)==1
    row=row[0];N=row['native_pulses'];assert N==3750720==sharp['native_pulses_per_U4']
    q=F(row['target_tree_phase_error_upper']);s=F(row['spectral_product_operator_upper'])
    assert q==F(sharp['target_tree_phase_error_upper']) and s==F(sharp['spectral_product_operator_upper'])
    x=N*b;e=(x+q)/(1-x);failure=N*b*b/(1-x)**2+e*e;gain=1/(1-x)**2+failure
    flow=(e*(2+e)+failure)/2+gain*s
    assert flow==F(sharp['complete_U4_halfdiamond_upper']) and gain==F(sharp['raw_map_gain_upper'])
    assert 512*flow==F(sharp['512_complete_U4_halfdiamond_upper'])
    source=read('inputs/first_cell_RECEIPT.json')['scientific']
    app=F(source['512epoch_apparatus_error_upper']);a=flow+app/512
    pump=512*a/(1-1024*a)
    field=read('FIELD_A_RESULT.json');assert field['all_F_A_gates']
    route=F(read('ROUTE_RECORDS.json')['route_failure_upper'])
    reader=F(field['terminal_8bit_error_upper']);native=pump+route+reader
    control_a=a+reader;control=512*control_a/(1-1024*control_a)+route+reader
    comparison=F(field['comparison_512_error_upper']);bad=F(source['512epoch_ideal_bad_weight_upper'])
    halfturn=2*(bad+native)
    control_observed=F(field['dephasing_coordinate_halfturn_TV'])
    control_lower=control_observed-2*(comparison+control)
    native_shape_lower=F(1,2)-(bad+native)
    native_Rx_lower=F(1,3)-2*(bad+native)
    assert native<F(field['reserved_native_error_allowance_NOT_YET_CERTIFIED'])
    assert halfturn<F(1,250)
    assert control_lower>F(1,100) and control_observed>10*(comparison+control)
    assert native_shape_lower>F(1,3) and native_Rx_lower>F(1,5)
    result={'candidate':'A-old240-compact-energy-apolar-join-v1','stage':'F-B certified native transfer',
      'native_implementation_transfer_complete':True,'original_raw_stage_F_executed':False,
      'same_group_review_not_G6':True,'fresh_native_flow_executions':0,
      'executed_primitive_catalogue_entries_consumed':count,'pulse_applications_per_flow':N,
      'flow_complete_halfdiamond_upper':str(flow),'apparatus512_upper':str(app),
      'gain_sensitive512pump_upper':str(pump),'route_failure_upper':str(route),
      'terminal8bit_upper':str(reader),'complete_native_error_upper':str(native),
      'complete_dephasing_control_error_upper':str(control),'certified_Q_comparison_error_upper':str(comparison),
      'ideal_bad_weight_upper':str(bad),'full_position_halfturn_TV_upper':str(halfturn),
      'active_dephasing_halfturn_TV_lower':str(control_lower),
      'isotropic_position_TV_lower':str(native_shape_lower),'Rx_position_TV_lower':str(native_Rx_lower),
      'same_density_source_sector':'P_D0/2; L=6 and Hlog=-1260log(209/210) on D0',
      'conditional_input_scope':'admitted actual source-unit checkpoint; unconditional over routing, pump and final read failures',
      'failure_budget_deleted':False,'historical_archive_invariance_claimed':False,
      'all_certified_stage_gates':True,'full_row_award_claimed':False,
      'bindings':{'field_result_sha256':sha((HERE/'FIELD_A_RESULT.json').read_bytes()),'sharp_certificate_sha256':sha((HERE/'inputs/SHARP_FLOW_CERTIFICATE.json').read_bytes()),'first_cell_receipt_sha256':sha((HERE/'inputs/first_cell_RECEIPT.json').read_bytes())}}
    return result
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args();started=time.perf_counter();result=run()
    if args.record:
        with (HERE/'TRANSFER_RESULT.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    else:assert read('TRANSFER_RESULT.json')==result
    print(json.dumps({'native_transfer_pass':True,'native_error':float(F(result['complete_native_error_upper'])),'halfturn_upper':float(F(result['full_position_halfturn_TV_upper'])),'control_lower':float(F(result['active_dephasing_halfturn_TV_lower'])),'seconds':time.perf_counter()-started}))
