"""Independent floor/ceiling source continuation; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,random,time
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();reg=json.loads((P/'PREREG.json').read_text())
 for n,h in reg['sha256'].items():assert sha(P/n)==h
 inputs=json.loads((P/'HANDOFFS.json').read_text());ev=json.loads(gzip.decompress((P/'RECORDS.json.gz').read_bytes()));r=json.loads((P/'RECEIPT.json').read_text())
 assert sha(P/'RECORDS.json.gz')==r['evidence_sha256']
 original=P.parent/'2026-11-15_codex_signed_pair_loading'
 if original.exists():
  assert sha(original/'MANIFEST.sha256')==inputs['manifest_sha256']
  assert sha(original/'RECORDS.json.gz')==inputs['evidence_sha256']
  oe=json.loads(gzip.decompress((original/'RECORDS.json.gz').read_bytes()))
  original_rows={x['index']:x for x in oe['trials']if x['law']=='native'}
  for h in inputs['handoffs']:assert h['actual_source']==original_rows[h['index']]['source_audit_only']==h['observed_source']
 byindex={h['index']:h for h in inputs['handoffs']};counts={'native':0,'hold':0};success=0
 graph=json.loads((P/'graph.json').read_text())['seams'];edges={tuple(e)for e in graph}
 for trial in ev['trials']+[ev['forced_native_timeout']]:
  forced='index'not in trial
  if forced:x=[0]*240;x[0]=1;x[14]=-1;rng=None
  else:x=byindex[trial['index']]['observed_source'][:];rng=random.Random(72351000+trial['index'])
  law=trial['law'];assert len(trial['planned_moves'])<=8
  assert trial['planned_moves'][-2:]==[[0,1,1],[1,0,1]]
  for move,planned in zip(trial['moves'],trial['planned_moves']):
   source,target,sign=planned;assert(source,target,sign)==(move['source'],move['target'],move['sign'])
   assert x[source]==sign and x[target]==0 and source//12==target//12
   assert 1<=len(move['records'])<=32
   for k,(a,b,c,j)in enumerate(move['records']):
    assert(a,b)in edges and {a,b}=={source,target}
    if not forced:assert c==rng.randrange(2)
    before=x[a]
    if law=='native':
     s=x[a]+x[b];x[a]=s//2+c*(s%2);x[b]=s-x[a];assert j==before-x[a]
    else:assert j==0
    counts[law]+=1
    assert sum(abs(v)for v in x)==2 and sum(x)==0 and all(v==0 for v in x[24:])
    if k<len(move['records'])-1:assert j==0
   assert move['arrived']==bool(j)
   if move['arrived']:assert x[source]==0 and x[target]==sign
   else:assert len(move['records'])==32
  complete=len(trial['moves'])==len(trial['planned_moves'])and all(m['arrived']for m in trial['moves'])
  assert complete==trial['success'] and x==trial['observed_source']==trial['source_audit_only']
  if complete:
   expected=[0]*240;expected[0]=1;expected[14]=-1
   assert x==expected and trial['pair_code']==14 and trial['positive_frame']==[0,1];success+=1
  else:assert trial['terminal']=='timeout'and trial['pair_code']is None
 assert success==12 and len(ev['trials'])==24 and not ev['forced_native_timeout']['success']
 assert counts['native']==r['literal_native_calls']==150 and counts['hold']==r['hold_body_calls']==384
 assert F(r['per_trial_timeout_upper'])==1-(1-F(1,2**32))**8
 for target,lo in [(0,0),(14,12)]:
  distance={target:0};queue=[target]
  for a in queue:
   for u,v in edges:
    b=v if u==a else u if v==a else None
    if b is not None and lo<=b<lo+12 and b not in distance:distance[b]=distance[a]+1;queue.append(b)
  assert len(distance)==12 and max(distance.values())<=3
 assert r['all_gates_pass']and all(r['gates'].values())
 print(json.dumps({'independent_check':True,'producer_imports':False,'native_events':counts['native'],
  'hold_events':counts['hold'],'actual_routed_pairs':success,'seconds':time.monotonic()-t},indent=2))
if __name__=='__main__':main()
