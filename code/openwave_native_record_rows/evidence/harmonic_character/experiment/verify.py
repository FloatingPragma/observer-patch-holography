"""Independent source-record and counter reader; imports no producer."""
import gzip,hashlib,importlib.util,json,sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
def load(name):return json.loads((ROOT/name).read_text())
def run():
    spec=importlib.util.spec_from_file_location('character_reader_native',ROOT/'inputs/native.py');mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    assert hashlib.sha256((ROOT/'inputs/native.py').read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    init=load('INITIALIZATION.json');idx=0;counts=[0]*8;calls=0
    for start in range(2):
      for word in range(8):
        location=start;currents=[]
        for bit in range(3):
            coin=word//(2**bit)%2
            out,record=mod.native(tuple(int(k==location) for k in range(12)),(0,1),coin);calls+=1
            assert record==(coin,(1-location)-coin)
            location=out.index(1);currents.append(record[1])
        counter=sum(abs(current)*2**bit for bit,current in enumerate(currents));counts[counter]+=1
        assert init['rows'][idx]=={'start':start,'word':word,'currents':currents,'counter':counter,'end':location};idx+=1
    assert calls==48 and counts==init['counts']==[2]*8
    matrices={m:[[0]*8 for _ in range(8)] for m in ('native','hold_mutation')};reject={m:0 for m in matrices};phase={m:0 for m in matrices};negative={m:0 for m in matrices};rows=events=0
    with gzip.open(ROOT/'inputs/compatible_records.jsonl.gz','rt') as source,gzip.open(ROOT/'COUNTERS.bin.gz','rb') as counters:
      for mode in ('native','hold_mutation'):
       for initial in range(2):
        for word in range(4096):
            left=1-initial;js=[]
            for bit in range(12):
                coin=word//(2**bit)%2;j=left-coin if mode=='native' else 0
                left-=j;js.append(j);events+=1
            xs=[abs(js[4*k]*js[4*k+1]*js[4*k+2])*(2*abs(js[4*k+3])-1) for k in range(3)];accepted=sum(xs)==0
            assert json.loads(next(source))==[mode,initial,word,js,xs,accepted]
            for c in range(8):
                final=(c+xs[0])%8 if accepted else c;sign=(-1)**((final-c)%2)
                record=counters.read(3);assert record==bytes((c,final,2*accepted+int(sign>0)));rows+=1
                if accepted:matrices[mode][c][final]+=1;phase[mode]+=sign;negative[mode]+=sign<0
                else:reject[mode]+=1
      assert source.read()=='' and counters.read(1)==b''
    receipt=load('RECEIPT.json')
    for mode in matrices:
        s=receipt['cases'][mode];assert s['transitions']==matrices[mode] and s['phase_sum']==phase[mode] and s['phase_negative']==negative[mode] and s['rejected']==reject[mode]
        assert sum(map(sum,matrices[mode]))==s['accepted']
    P=[[F(x,5656) for x in row] for row in matrices['native']]
    for i,row in enumerate(P):
        assert sum(row)==1 and row[i]==F(99,101)
        assert sum(row[j] for j in range(8) if j!=i)==F(2,101)
        for j in range(8):assert P[i][j]==P[j][i]
        # Exact alternating-character eigenvector check, no trigonometry.
        assert sum(row[j]*(-1)**j for j in range(8))==F(97,101)*(-1)**i
    assert F(phase['native'],sum(map(sum,matrices['native'])))==F(679,707)
    assert phase['hold_mutation']==65536 and reject['native']==20288 and negative['native']==896
    assert rows==131072 and events==196608
    result={'all_pass':True,'new_native_verification_calls':calls,'archived_source_and_control_events_checked':events,
      'counter_rows_checked':rows,'imports_producer':False,'exact_positive_floor':'97/101',
      'all_eigenvalues_positive_by_symmetric_diagonal_dominance':True,'same_group':True,'G6_pass':False,
      'receipt_sha256':hashlib.sha256((ROOT/'RECEIPT.json').read_bytes()).hexdigest()}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
