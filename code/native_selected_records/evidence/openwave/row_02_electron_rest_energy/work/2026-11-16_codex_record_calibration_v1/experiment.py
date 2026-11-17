"""One continuing source per trial: unknown input, admission, observed calibration."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time
import traceback
from integration import graph,start
from compiler import catalog

HERE=Path(__file__).resolve().parent
PROGRESS=[]
SEEDS=(72371001,72371002)


def require(ok,message):
    if not ok:raise AssertionError(message)


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'))
def digest(value):return hashlib.sha256(canonical(value).encode()).hexdigest()


def run():
    begin=time.monotonic();PROGRESS.clear()
    inputs=json.loads(gzip.decompress((HERE/'ADMISSION_INPUTS.json.gz').read_bytes()))['rows']
    public=json.loads((HERE/'PUBLIC_GRAPHS.json').read_text());catalogs={};classes={};total_calls={'native':0,'always_swap':0}
    max_calls=0;ready=0;no_unit=0;calibration_calls=0;prefix_calls=0
    for input_index,old in enumerate(inputs):
        g=graph(old['B'],old['M'])
        require(json.loads(canonical(g))==public[f"{old['B']}:{old['M']}"],'public graph binding')
        for law in ('native','always_swap'):
            for seed in SEEDS:
                trial={'input_index_audit_only':input_index,'B':old['B'],'M':old['M'],
                       'admission_seed_audit_only':old['admission_seed_audit_only'],
                       'calibration_seed_audit_only':seed,'probe_law_audit_only':law,
                       'total_audit_only':old['total_audit_only'],'initial_audit_only':old['initial_audit_only']}
                PROGRESS.append(trial)
                request,audit,record=start(old['initial_audit_only'],g,old['admission_seed_audit_only'],seed,law)
                trial['observer']=record;trial['audit']=audit()
                entry=record['admission'];require(entry==old['expected_controller'],'actual admission replay differs')
                prefix=old['expected_audit']['events'];n=len(prefix)
                require(trial['audit']['events'][:n]==prefix,'actual admission native word differs')
                prefix_calls+=n
                for key,value in trial['audit']['literal_endpoint_calls'].items():total_calls[key]+=value
                if entry['outcome']=='no_unit':
                    require(record['outcome']=='not_admitted' and len(trial['audit']['events'])==n,'zero-input terminal')
                    no_unit+=1;continue
                require(entry['outcome']=='signed_target_unit' and record['outcome']=='ready','fair calibration incomplete')
                require(trial['audit']['final']==entry['observed_final'],'continuing particle or spectator not returned')
                observed=record['calibration'];blank=record['blank'];calls=observed['endpoint_calls']+blank['endpoint_calls']
                require(calls==len(trial['audit']['events'])-n and calls<=37568,'literal calibration budget')
                require(observed['return_position']==1 and len(observed['columns'])==18,'full column/return outcome')
                max_calls=max(max_calls,calls);calibration_calls+=calls;ready+=1
                signature=digest(observed['columns']);trial['measured_table_sha256']=signature
                if signature not in catalogs:catalogs[signature]=catalog(observed['columns'],blank['records'])
                key=law+':'+str(entry['sign'])
                if key in classes:require(classes[key]['table_sha256']==signature,'background/stream changes measured law')
                else:classes[key]={'table_sha256':signature,'columns':observed['columns'],'source_trial_index':len(PROGRESS)-1,
                                   'actual_blank_records':blank['records'],'observed_sign':entry['sign']}
        if (input_index+1)%100==0:print(json.dumps({'completed_inputs':input_index+1,'continuing_trials':len(PROGRESS)}),flush=True)
    require(len(PROGRESS)==5256 and ready==5184 and no_unit==72,'full registered input coverage')
    require(set(classes)=={'native:1','native:-1','always_swap:1','always_swap:-1'},'all observed source classes')
    for law in ('native','always_swap'):
        a=catalogs[classes[law+':1']['table_sha256']]['cooling5']['P']
        b=catalogs[classes[law+':-1']['table_sha256']]['cooling5']['P']
        require(a==b,'measured fair-coin mean sign dependence')
    native=catalogs[classes['native:1']['table_sha256']]['cooling5']['P']
    control=catalogs[classes['always_swap:1']['table_sha256']]['cooling5']['P']
    require(native!=control,'same observer failed to distinguish the probe law')
    # Inconclusive controls use the same continuing-source entry. Their biased
    # streams are declared detector controls, outside the fair-coin bound.
    old=next(r for r in inputs if r['total_audit_only']>0);g=graph(old['B'],old['M']);timeouts=[]
    for name,forced,expected in [('route',{'calibration_route':0},'route_timeout'),
                                 ('coupon',{'calibration_probe':0},'coupon_timeout'),
                                 ('blank',{'calibration_blank':0},'blank_timeout')]:
        request,audit,record=start(old['initial_audit_only'],g,old['admission_seed_audit_only'],SEEDS[0],'native',forced=forced)
        require(record['outcome']==expected,'timeout detector failed '+name)
        timeouts.append({'kind':name,'forced_coin_control':forced,'initial_audit_only':old['initial_audit_only'],
                         'B':old['B'],'M':old['M'],'observer':record,'audit':audit()})
    bounds={'cap':64,'core_calibration_failure_upper':str(F(595,2**64)),
            'with_blank_failure_upper':str(F(597,2**64)),'success_probability_lower':str(1-F(597,2**64)),
            'max_calibration_endpoint_calls':37568,'expected_endpoint_calls_upper':86,
            'conditional_expected_calls_upper':str(F(86)/(1-F(597,2**64)))}
    return {'candidate':'A-native-continuing-record-calibration-v1','passed':True,
            'trials':PROGRESS,'trial_count':len(PROGRESS),'ready_trials':ready,'no_unit_trials':no_unit,
            'representative_measured_tables':classes,'compiled_catalogs':catalogs,'timeout_controls':timeouts,
            'literal_trial_calls':total_calls,'literal_admission_prefix_calls':prefix_calls,
            'literal_calibration_calls':calibration_calls,'maximum_observed_calibration_calls':max_calls,
            'literal_timeout_control_calls':{law:sum(t['audit']['literal_endpoint_calls'][law] for t in timeouts) for law in ('native','always_swap')},
            'calibration_bounds':bounds,'seconds':time.monotonic()-begin,
            'quantum_coherence_demonstrated_by_calibration':False,'absolute_vacuum_derived_from_zero_current':False,
            'source_reset_after_admission':False,'source_reader_receives_hidden_input_labels':False,
            'external_G6':False,'official_awards_added':0}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record-to',type=Path);args=parser.parse_args()
    destination=args.record_to
    if destination is not None:destination.mkdir(parents=True,exist_ok=False)
    try:
        prereg=json.loads((HERE/'PREREG.json').read_text())
        for name,expected in prereg['files'].items():require(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,'frozen input '+name)
        result=run();require(result['passed'],'scientific verdict')
        if destination is None:
            old=json.loads(gzip.decompress((HERE/'records/EVIDENCE.json.gz').read_bytes()))
            left=dict(result);right=dict(old);left.pop('seconds');right.pop('seconds')
            require(left==right,'saved evidence differs')
            receipt=json.loads((HERE/'records/RECEIPT.json').read_text())
            require(receipt['passed'] and receipt['evidence_sha256']==hashlib.sha256((HERE/'records/EVIDENCE.json.gz').read_bytes()).hexdigest(),'saved receipt differs')
        else:
            blob=gzip.compress((canonical(result)+'\n').encode(),mtime=0)
            with (destination/'EVIDENCE.json.gz').open('xb') as stream:stream.write(blob)
            receipt={k:v for k,v in result.items() if k not in ('trials','compiled_catalogs','timeout_controls')}
            receipt['evidence_sha256']=hashlib.sha256(blob).hexdigest()
            with (destination/'RECEIPT.json').open('x') as stream:json.dump(receipt,stream,indent=2);stream.write('\n')
        print(json.dumps({k:v for k,v in result.items() if k not in ('trials','compiled_catalogs','timeout_controls','representative_measured_tables')},indent=2));return 0
    except Exception as error:
        failure={'passed':False,'error':repr(error),'traceback':traceback.format_exc(),'all_completed_or_current_trials':PROGRESS}
        if destination is not None:
            with (destination/'FAILURE.json').open('x') as stream:json.dump(failure,stream,indent=2);stream.write('\n')
        print(json.dumps({'passed':False,'error':repr(error)}));return 1


if __name__=='__main__':raise SystemExit(main())
