"""48 new native initialization calls, then exact archived-word counter extension."""
import gzip,hashlib,importlib.util,json,sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
SOURCE='43bec8c495e9740e87947bbf3f0d0988c6fc084c4f1794ee5c6733241056a457'
def dump(name,value):(ROOT/name).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def run():
    pin=ROOT/'inputs/native.py'
    assert hashlib.sha256(pin.read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    assert hashlib.sha256((ROOT/'inputs/compatible_records.jsonl.gz').read_bytes()).hexdigest()==SOURCE
    spec=importlib.util.spec_from_file_location('character_native_pin',pin);mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    init=[];counts=[0]*8;calls=0
    for start in (0,1):
      for word in range(8):
        loads=tuple(int(i==start) for i in range(12));js=[]
        for i in range(3):
            coin=(word>>i)&1;loads,(c,j)=mod.native(loads,(0,1),coin);calls+=1
            assert c==coin;js.append(j)
        counter=sum(abs(j)<<i for i,j in enumerate(js));counts[counter]+=1
        init.append({'start':start,'word':word,'currents':js,'counter':counter,'end':loads.index(1)})
    assert counts==[2]*8 and calls==48
    dump('INITIALIZATION.json',{'rows':init,'counts':counts,'native_calls':calls,'new_preparations':16})
    stats={name:{'accepted':0,'rejected':0,'phase_sum':0,'phase_negative':0,'transitions':[[0]*8 for _ in range(8)]} for name in ('native','hold_mutation')}
    rows=0
    with gzip.open(ROOT/'inputs/compatible_records.jsonl.gz','rt') as source,(ROOT/'COUNTERS.bin.gz').open('wb') as sink:
      with gzip.GzipFile(filename='',mode='wb',fileobj=sink,mtime=0) as tape:
        for line in source:
            mode,start,word,currents,values,accepted=json.loads(line)
            assert values==[(1 if abs(currents[4*k+3]) else -1) if all(abs(x)==1 for x in currents[4*k:4*k+3]) else 0 for k in range(3)]
            assert accepted==(sum(values)==0)
            for initial in range(8):
                final=(initial+values[0])%8 if accepted else initial
                phase=1 if (final-initial)%2==0 else -1
                tape.write(bytes((initial,final,(2 if accepted else 0)+(1 if phase==1 else 0))))
                s=stats[mode];rows+=1
                if accepted:
                    s['accepted']+=1;s['phase_sum']+=phase;s['phase_negative']+=phase<0;s['transitions'][initial][final]+=1
                else:s['rejected']+=1
    native=stats['native'];control=stats['hold_mutation'];lam=F(native['phase_sum'],native['accepted'])
    matrix=[[F(x,5656) for x in row] for row in native['transitions']]
    target=[[F(99,101) if i==j else F(1,101) if (i-j)%8 in (1,7) else F(0) for j in range(8)] for i in range(8)]
    assert matrix==target
    checks={'uniform_native_counter':counts==[2]*8,'exact_connected_kernel':matrix==target,
       'positive_character':lam==F(679,707),'hold_character_one':control['phase_sum']==control['accepted']==65536,
       'all_rejected_words_retained':native['rejected']==20288 and control['rejected']==0,
       'native_negative_phases_retained':native['phase_negative']==896,'record_rows':rows==131072,
       'changed_law_gap':1-lam==F(28,707)}
    result={'candidate':'A-native-harmonic-character-v1','all_pass':all(checks.values()),'checks':checks,
       'new_native_calls':48,'new_preparations':16,'new_full_source_word_calls':0,
       'reused_native_calls':98304,'reused_mutation_events':98304,'derived_counter_rows':rows,
       'source_archive_sha256':SOURCE,'counter_archive_sha256':hashlib.sha256((ROOT/'COUNTERS.bin.gz').read_bytes()).hexdigest(),
       'cases':stats,'conditional_kernel':[[str(x) for x in row] for row in matrix],
       'k4_character':str(lam),'all_mode_positive_floor':'97/101','actual_joint_counter_tape':False,
       'physical_Maxwell_execution':False,'G6_pass':False}
    dump('RECEIPT.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['cases','conditional_kernel']},sort_keys=True))
    return result['all_pass']
if __name__=='__main__':raise SystemExit(0 if run() else 1)
