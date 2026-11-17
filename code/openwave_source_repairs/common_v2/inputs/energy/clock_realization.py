"""A common counted clock for both finite native processor programmes."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'inputs'))
from native_event import native
from finalize import finish
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    started=time.monotonic();reg=json.loads((HERE/'PREREG_CLOCK.json').read_text())
    for n,h in reg['files'].items():assert sha(HERE/n)==h,n
    flow=json.loads((HERE/'attempts/v2/RECEIPT.json').read_text())
    pulses={k:int(v['total_pulses'])for k,v in flow['words'].items()}
    # lcm(1260,1280)=80640. Input contains no source state or energy value.
    C=80640*((10*max(pulses.values())+80639)//80640)
    assert C%1260==0 and C%1280==0
    programmes=[]
    for duration in (F(1,10),F(4)):
        slots=C*duration;assert slots.denominator==1
        slots=int(slots);ticks=1260*duration;assert ticks.denominator==1
        for name,body in pulses.items():
            assert body%128==0 and slots%128==0 and slots>=body
            programmes.append({'object':name,'duration':str(duration),'total_clock_slots':slots,
              'inherited_count_ticks':int(ticks),'macrosteps':128,'slots_per_macrostep':slots//128,
              'processor_pulses_per_macrostep':body//128,'idle_slots_per_macrostep':(slots-body)//128,
              'represented_active_processor_pulses':body,'represented_idle_slots':slots-body,
              'complete_clock_record_law':{'alphabet':[[0,0],[1,0]],'length':slots,
                    'every_binary_history_probability':'2^(-'+str(slots)+')',
                    'summed_probability':'1','coin_dependent_scheduling':False},
              'on_processor_abort':'Retain first-failure flag, skip remaining source pulses, continue the same total zero-source clock slots.'})
    # Literal checks of the exact two native idle branches and longer
    # histories. The enormous deterministic repetitions above are not
    # claimed as these literal calls or sampled histories.
    records=[]
    for coin in (0,1):
        q,r=native((0,0),(0,1),coin)
        records.append({'input':[0,0],'coin':coin,'output':list(q),'record':list(r)})
        assert q==(0,0)and r==(coin,0)
    histories=[]
    for mask in range(16):
        q=(0,0);tape=[]
        for j in range(16):
            coin=(mask>>(j%4))&1;q,r=native(q,(0,1),coin);tape.append(list(r))
        assert q==(0,0)and all(j==0 for _,j in tape)
        histories.append({'seed_word':mask,'records':tape,'observed_tick_count':len(tape)})
    single={x['duration']:x for x in programmes if x['object']=='single'}
    pair={x['duration']:x for x in programmes if x['object']=='pair'}
    gates={'native_zero_source_branches_exact':all(x['output']==[0,0]and x['record']==[x['coin'],0]for x in records),
      'same_clock_count_for_both_objects':all(single[t]['total_clock_slots']==pair[t]['total_clock_slots']for t in single),
      'duration_ratio_exact40':single['4']['total_clock_slots']==40*single['1/10']['total_clock_slots'],
      'nonnegative_integral_padding':all(x['idle_slots_per_macrostep']>=0 for x in programmes),
      'source_tick_conversion_exact':all(x['total_clock_slots']==(C//1260)*x['inherited_count_ticks']for x in programmes),
      'finite_counter_capacity':all(x['total_clock_slots']<2**40 for x in programmes)}
    result={'candidate':reg['candidate'],'prereg_sha256':sha(HERE/'PREREG_CLOCK.json'),
      'clock_slots_per_time_unit':C,'clock_slots_per_inherited_count_tick':C//1260,
      'programmes':programmes,'literal_zero_source_calls':len(records)+sum(len(x['records'])for x in histories),
      'native_branch_records':records,'literal_histories':histories,
      'arbitrary_reference_idle_identity':'For each coin c the factorized branch is2^(-1/2) I_source,reference,spectator tensor|c,current0>; every finite repeated history is retained. Summing its branches is exactly the identity channel.',
      'supplied_architecture':['auxiliary zero source','fair native coin law','40-bit classical observer counter','slot scheduler and declared count-time scale'],
      'time_scope':'One observer count clock for both objects and both durations; not wall-clock time or total repair-plus-processor resource count. Counter and rate are supplied apparatus, not emergence credit.',
      'source_effect_error':'0 exactly under the native zero-source law and disjoint auxiliary scheduling',
      'source_intervals_and_words_unchanged':True,'gates':gates,'elapsed_seconds':time.monotonic()-started}
    out=HERE/'clock';out.mkdir(exist_ok=False)
    status=finish(out/'RECEIPT.json',result)
    print(json.dumps({'gates':gates,'clock_rate':C,'literal_calls':result['literal_zero_source_calls'],'seconds':result['elapsed_seconds']}))
    raise SystemExit(status)
if __name__=='__main__':main()
