"""Portable primary: native source -> regional records -> native cycle apparatus.

Python standard library only. No official acceptance result is manufactured.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,itertools
from itertools import cycle
import native_charge_winding as spatial
import native_cycle_meter as meter
import joined_readout as joined
from native_kernel import native
HERE=Path(__file__).resolve().parent


def require(value,message):
    if not value:raise ValueError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def complete_law(reset,hop,hops):
    success=reset*hop**hops
    failure=1-success
    return {'success':str(success),'timeout':str(failure)}


def validate_law(law):
    require(F(law['success'])>=0 and F(law['timeout'])>=0,'negative outcome')
    require(F(law['success'])+F(law['timeout'])==1,'failure outcome omitted')


def execute(packet,retries=16):
    edges=tuple(tuple(e) for e in packet['seams'])
    spatial_result=spatial.evaluate(packet)
    require(all(spatial_result['checks'].values()),'spatial construction failed')
    hop=1-F(1,2**retries)
    reset=[meter.reset_success(v,retries) for v in range(12)]
    # Factorized apparatus law: every conditional label retains all failures.
    templates=[]
    for initial in range(12):
        for cls,winding in ((0,0),(1,1),(5,-1)):
            law=complete_law(reset[initial],hop,6*abs(winding));validate_law(law)
            replay=meter.execute_class(initial,cls,meter.favorable_coins(initial,(winding,),retries),retries)
            require(replay['success'] and replay['winding']==winding,'actual pointer readout')
            state=tuple(int(v==initial) for v in range(12))
            for event in replay['records']:
                state,(_,j)=native(state,tuple(event['seam']),event['coin'])
                require(j==event['current'],'apparatus record is not native')
            require(list(state)==replay['final_loads'],'apparatus final state')
            require(F(replay['fine_history_probability'])==F(1,2**len(replay['records'])),'apparatus history probability')
            templates.append(dict(initial_pointer_port=initial,counter_class=cls,decoded_spatial_winding=winding,decoded_relative_charge=-winding,unconditional_law=law,example_full_native_execution=replay))
    output=[];histories=events=noops=0;lost_current_detected=False;joined_checks=0;joined_examples=[]
    for setting,seed in enumerate(packet['spatial_result']['seed_records']):
        probability=F()
        for branch in seed['complete_prescribed_schedule_histories']:
            G=tuple(seed['initial_G']);state=G;counter=[0]*12;tape=[0]*30
            for event in branch['events']:
                e=event['seam_index'];a,b=edges[e]
                state,(_,j)=native(state,(a,b),event['coin']);tape[e]+=j
                counter[a]=(counter[a]+j)%6;counter[b]=(counter[b]-j)%6
                require(j==event['current'],'source current')
                require(counter==[(g-n)%6 for g,n in zip(G,state)],'counter/native equivariance')
                if j:
                    damaged=counter.copy();damaged[a]=(damaged[a]-j)%6;damaged[b]=(damaged[b]+j)%6
                    lost_current_detected |= damaged!=[(g-n)%6 for g,n in zip(G,state)]
                events+=1;noops+=int(j==0)
            require(tape==branch['final_tape'] and list(state)==branch['final_loads'],'source final state')
            schedule=[tuple(x['seam']) for x in seed['actual_history']]
            coins=[x['coin'] for x in branch['events']]
            require([tuple(x['seam']) for x in branch['events']]==schedule,'declared source schedule changed')
            actual_prefix=joined.source_prefix(G,schedule,coins)
            require(actual_prefix['counter']==counter and actual_prefix['source_final_loads_for_audit_only']==list(state),'joined prefix vs independent composition')
            p=F(branch['probability']);probability+=p
            # All twelve port contours. Arbitrary regions are obtained by the
            # same finite regional sum and the separately proven boundary law.
            readings=[]
            for port in range(12):
                cls=counter[port];require(cls in (0,1,5),'outside one-load charge sector')
                winding={0:0,1:1,5:-1}[cls]
                require(-winding==state[port]-G[port],'present relative charge')
                initial=(setting+port)%12
                actual=joined.execute_joined(G,schedule,coins,1<<port,initial,
                    meter.favorable_coins(initial,(winding,),retries),retries)
                require(actual['success'] and actual['relative_charge']==-winding and actual['source']['counter']==counter,'actual joined vs independent composition')
                require(F(actual['combined_history_probability'])==F(1,2**len(actual['combined_native_tape'])),'joined faircoin probability')
                require(actual['clock']['total_native_attempts']==len(actual['combined_native_tape']),'joined event clock')
                joined_checks+=1
                if len(joined_examples)<6:joined_examples.append(actual)
                readings.append(dict(region_ports=[port],regional_counter_class=cls,relative_charge=-winding,apparatus_template_indices=[3*v+(0 if cls==0 else 1 if cls==1 else 2) for v in range(12)]))
            output.append(dict(setting=setting,source_probability=str(p),source_coin_record=[e['coin'] for e in branch['events']],initial_G=list(G),present_n=list(state),counter=counter,contour_readings=readings))
            histories+=1
        require(probability==1,'source failure/history omitted')
    # Actual joined stopping branches; none reports a guessed charge.
    G=[1]+[0]*11;edge=(0,1);failure_traces=[]
    for name,initial,coins,reason,region in (
        ('discovery',0,joined.discovery_hold_coins(0,retries),'discovery_timeout',1),
        ('routing',0,cycle((0,)),'routing_timeout',1),
        ('measurement',1,cycle((1,)),'measurement_timeout',1),
        ('zero_class_discovery',0,joined.discovery_hold_coins(0,retries),'discovery_timeout',4095)):
        actual=joined.execute_joined(G,(edge,),(0,),region,initial,coins,retries)
        require(not actual['success'] and actual['reason']==reason and actual['relative_charge'] is None,'joined timeout erased or guessed')
        require(F(actual['combined_history_probability'])==F(1,2**len(actual['combined_native_tape'])),'joined failure fine probability')
        state=tuple(int(v==initial) for v in range(12))
        for event in actual['pointer']['records']:
            state,(_,j)=native(state,tuple(event['seam']),event['coin'])
            require(j==event['current'],'joined failure is not actual native')
        require(list(state)==actual['pointer']['final_loads'],'joined failure final state')
        failure_traces.append(dict(case=name,execution=actual))
    # Arbitrary separating regions on the complete primary state domain.
    region_checks=0
    for origin,occupied in itertools.product(range(12),repeat=2):
        q=[int(v==origin)-int(v==occupied) for v in range(12)]
        c=[v%6 for v in q]
        for region in range(4096):
            value=sum(q[v] for v in range(12) if region>>v&1)
            cls=sum(c[v] for v in range(12) if region>>v&1)%6
            require({0:0,1:1,5:-1}[cls]==value,'regional counter alias')
            region_checks+=1
    mutations={'forgotten_native_current_rejected':bool(lost_current_detected),'larger_charge_alias_rejected':{0:0,1:1,5:-1}[7%6]!=7}
    try:validate_law({'success':templates[1]['unconditional_law']['success'],'timeout':'0'})
    except ValueError:mutations['unconditional_failure_omission_rejected']=True
    else:mutations['unconditional_failure_omission_rejected']=False
    rejected=[meter.execute_class(0,c,iter(()),retries) for c in (2,3,4)]
    checks=dict(complete_native_source_histories=histories==456,actual_source_current_counter_readout=events>0 and noops>0,actual_native_pointer_all_unknown_inputs=len(templates)==36,external_joined_executor_matches_independent_composition=joined_checks==456*12,all_unconditional_apparatus_laws_complete=all(F(t['unconditional_law']['success'])+F(t['unconditional_law']['timeout'])==1 for t in templates),all_primary_regions=region_checks==589824,all_required_mutations_rejected=all(mutations.values()),actual_joined_timeout_stages_retained=len(failure_traces)==4,counter_complements_retained=all(not x['success'] and x['winding'] is None for x in rejected))
    return dict(schema='oph.row1.portable_primary.v1',checks=checks,spatial_checks=spatial_result['checks'],spatial_max_phase_error=spatial_result['maximum_phase_degree_error'],source_histories=histories,joined_executor_comparisons=joined_checks,joined_execution_examples=joined_examples,actual_joined_failure_traces=failure_traces,source_events=events,source_noops=noops,regional_counter_comparisons=region_checks,retries_per_requested_hop=retries,maximum_primary_meter_attempts=39*retries,worst_primary_meter_timeout=max(float(F(t['unconditional_law']['timeout'])) for t in templates),apparatus_templates=templates,source_to_apparatus_factorized_joint_law=output,retained_counter_complements=rejected,mutations=mutations,joint_probability_rule='P(source history,pointer outcome|setting,initial pointer)=source_probability * corresponding unconditional apparatus template probability; all source histories sum1 and each template success+timeout=1',scope='Actual literal integer native process; Q=n−retained G=−spatial winding. Native-only apparatus uses observed currents, finite controller and an additional one-load carrier. Bounded decoder applies to this primary one-load source. Spatial S1 representative is a declared diagnostic map. No coherent control service, C10 transfer, physical electric normalization or official acceptance claim.',full_criterion_met=False,official_rows_changed=[])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=HERE/'PRIMARY_RESULT.json');args=parser.parse_args()
    manifest=json.loads((HERE/'manifest.json').read_text())
    for name,pin in manifest['files_sha256'].items():require(sha(HERE/name)==pin,'public source drift: '+name)
    result=execute(json.loads((HERE/'source_packet.json').read_text()))
    require(all(result['checks'].values()),'primary failure')
    result['public_manifest_sha256']=sha(HERE/'manifest.json')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('checks','source_histories','source_events','source_noops','worst_primary_meter_timeout','mutations')},indent=2))

if __name__=='__main__':main()
