"""Independent arithmetic and factorization audit of the native count clock."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,math
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    a=argparse.ArgumentParser();a.add_argument('source',type=Path);a.add_argument('--record',action='store_true');args=a.parse_args()
    p=args.source.resolve();r=read(p/'clock/RECEIPT.json');reg=read(p/'PREREG_CLOCK.json')
    assert r['all_gates_pass'] and all(r['gates'].values())
    assert sha(p/'PREREG_CLOCK.json')==r['prereg_sha256']
    for n,h in reg['files'].items():assert sha(p/n)==h,n
    flow=read(p/'attempts/v2/RECEIPT.json');body={k:v['total_pulses']for k,v in flow['words'].items()}
    multiple=math.lcm(1260,1280);rate=multiple*math.ceil(F(10*max(body.values()),multiple))
    assert rate==r['clock_slots_per_time_unit'] and rate//1260==r['clock_slots_per_inherited_count_tick']
    counts={}
    for row in r['programmes']:
        name=row['object'];t=F(row['duration']);slots=rate*t;assert slots.denominator==1
        N=int(slots);assert N<2**40 and N>=body[name]
        assert N==row['total_clock_slots'] and N==128*row['slots_per_macrostep']
        assert body[name]==128*row['processor_pulses_per_macrostep']==row['represented_active_processor_pulses']
        assert N-body[name]==128*row['idle_slots_per_macrostep']==row['represented_idle_slots']
        assert F(N,rate)==F(row['inherited_count_ticks'],1260)==t
        law=row['complete_clock_record_law'];assert law['alphabet']==[[0,0],[1,0]]
        assert law['length']==N and law['every_binary_history_probability']==f'2^(-{N})'
        assert law['summed_probability']=='1' and not law['coin_dependent_scheduling']
        counts[name,t]=N
    assert counts['single',F(4)]==counts['pair',F(4)]==40*counts['single',F(1,10)]==40*counts['pair',F(1,10)]
    for row in r['native_branch_records']:
        coin=row['coin'];total=sum(row['input']);left=total//2+coin*(total%2);right=total-left
        assert row['input']==[0,0] and row['output']==[left,right]==[0,0] and row['record']==[coin,0]
    for row in r['literal_histories']:
        assert len(row['records'])==row['observed_tick_count']==16
        assert row['records']==[[(row['seed_word']>>(j%4))&1,0]for j in range(16)]
    assert r['literal_zero_source_calls']==2+16*16==258
    # Each independent native clock branch has Kraus factor 1/sqrt(2) I.
    # Multiplication gives 2^(-N/2) I, with exactly2^N orthogonal records.
    # This symbolic integer-exponent identity avoids materializing huge tapes.
    result={'clock_receipt_sha256':sha(p/'clock/RECEIPT.json'),'rate':rate,
        'slots_per_inherited_tick':rate//1260,'duration_ratio':40,'programmes_checked':4,
        'literal_zero_source_records_checked':258,'arbitrary_reference_factorization_checked':True,
        'all_huge_histories_materialized':False,'literal_new_native_calls':0,'producer_imports':False,
        'clock_scale_and_counter_supplied':True}
    if args.record:
        with(HERE/'CLOCK_RESULT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    else:assert read(HERE/'CLOCK_RESULT.json')==result
    print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
