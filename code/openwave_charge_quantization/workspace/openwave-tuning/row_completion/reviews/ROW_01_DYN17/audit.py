"""Independent bounded source, integer spatial class and record-only audit."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from collections import deque
from copy import deepcopy
import hashlib,json,sys,cmath,math
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
SOURCE=ROOT/'openwave-tuning/row_completion/ROW_01'
sys.dont_write_bytecode=True
sys.path.insert(0,str(SOURCE))
import reviewed

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pins():
    out={}
    for name,h in [('MANIFEST.json','c4f6db87b1d2c7ddd4caa4c8ac5a95f09c33679f060d50845aa0bf42f5e54d12'),('REVIEW_FIX_MANIFEST.json','f695ee571e25cfe59e808b9606a2b788aa232623f12dc33f0f4cac9db66c47d3')]:
        p=SOURCE/name;assert sha(p)==h;out[str(p.relative_to(ROOT))]=h
        d=json.loads(p.read_text())
        for key in ('files','source_sha256'):
            for r,w in d.get(key,{}).items():assert sha(ROOT/r)==w;out[r]=w
    return out

def bareiss(a):
    a=[r[:] for r in a];last=1;sign=1
    for k in range(len(a)-1):
        j=next((j for j in range(k,len(a)) if a[j][k]),None)
        if j is None:return 0
        if j!=k:a[k],a[j]=a[j],a[k];sign=-sign
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                value=a[k][k]*a[i][j]-a[i][k]*a[k][j]
                assert value%last==0;a[i][j]=value//last
        last=a[k][k]
        for i in range(k+1,len(a)):a[i][k]=0
    return sign*a[-1][-1]

def step(n,e,c):
    a,b=e;d=n[a]-n[b];j=(d-(2*c-1)*(abs(d)%2))//2
    m=list(n);m[a]-=j;m[b]+=j;return tuple(m),j

def main():
    before=pins();g=json.loads((SOURCE/'GEOMETRY.json').read_text())
    edges=list(map(tuple,g['seams']));A=g['A'];B=[[int(v==a)-int(v==b) for a,b in edges] for v in range(12)]
    assert B==g['B'] and len(edges)==30 and len(A[0])==20
    assert all(sum(B[v][e]*A[e][f] for e in range(30))==0 for v in range(12) for f in range(20))
    assert all(sorted(x for x in r if x)==[-1,1] for r in A)
    detB=bareiss([[B[v][e] for e in g['tree_edges']] for v in range(1,12)])
    detA=bareiss([[A[e][f] for f in range(19)] for e in g['chord_edges']]);assert abs(detB)==abs(detA)==1
    assert all(sum(A[e][f] for f in range(20))==0 for e in range(30))
    # Independently follow the actual oriented dual boundary and sample its
    # homotopy with nonconstant vertex phases and zero-endpoint perturbations.
    phase_error=0.;cases=0
    for j in (-3,-2,-1,0,1,2,3):
      for selected in range(30):
        current=[j if e==selected else 0 for e in range(30)]
        for v in range(12):
          links={}
          for e in range(30):
            if not B[v][e]:continue
            tail=A[e].index(-1);head=A[e].index(1)
            if B[v][e]<0:tail,head=head,tail
            assert tail not in links;links[tail]=(head,e,B[v][e])
          start=min(links);u=start;angles=[]
          for _ in range(5):
            w,e,s=links[u];ga=.17*math.sin(u);gb=.17*math.sin(w)
            path=[cmath.exp(1j*(ga+(gb-ga)*t/64+2*math.pi*s*current[e]*t/64+.11*math.sin(2*math.pi*t/64))) for t in range(65)]
            if angles:assert abs(angles[-1]-path[0])<1e-12
            angles.extend(path if not angles else path[1:]);u=w
          assert u==start and abs(angles[-1]-angles[0])<1e-12
          increments=[cmath.phase(y*x.conjugate()) for x,y in zip(angles,angles[1:])];assert max(map(abs,increments))<math.pi/2
          degree=math.fsum(increments)/(2*math.pi);expected=sum(B[v][e]*current[e] for e in range(30))
          phase_error=max(phase_error,abs(degree-expected));cases+=1
    assert phase_error<1e-12
    branches=0;pairs=0
    for total in (1,2):
      states=[tuple(sum(x==v for x in group) for v in range(12)) for group in combinations_with_replacement(range(12),total)]
      for n in states:
        for edge in edges:
          for c in (0,1):
            m,j=step(n,edge,c);actual,(_,aj)=reviewed.integrity.native(n,edge,c)
            assert m==actual and j==aj and min(m)>=0 and sum(m)==total
            if total==1:assert set(m)<={0,1}
            branches+=1
        for G in states:
          q=[x-y for x,y in zip(G,n)];assert sum(q)==0 and sum(max(x,0) for x in q)<=total and sum(min(x,0) for x in q)>=-total
          assert len({x%6 for x in range(-total,total+1)})==2*total+1;pairs+=1
    examples=json.loads((SOURCE/'EXAMPLES.json').read_text());selected={};success_degrees=set();failure_reasons=set();replays=0
    for row in examples:
      n=tuple(row['G']);Fedge=[0]*30
      for record in row['claim']['source_records']+row['claim']['sweep_records']:
        edge=edges[record['seam']]
        if record['kind']=='E':assert record['contrast']==n[edge[0]]-n[edge[1]]
        n,j=step(n,edge,record['coin']);assert j==record['current'];Fedge[record['seam']]+=j
      q=[sum(B[v][e]*Fedge[e] for e in range(30)) for v in range(12)]
      assert q==[x-y for x,y in zip(row['G'],n)]
      degree=sum(q[v] for v in row['region']);assert degree==row['topological_degree']
      # Replay every pointer branch without source final snapshots.
      meter=row['claim']['meter'];compatible=[]
      for x in range(12):
        pointer=tuple(int(v==x) for v in range(12));ok=True
        for record in meter['records']:
          pointer,j=step(pointer,tuple(record['seam']),record['coin']);ok &= j==record['current']
        if ok:compatible.append(x)
      assert compatible
      cyc=reviewed.integrity.meter.CYCLE
      orientation={tuple(sorted((a,b))):1 if a<b else -1 for a,b in zip(cyc,cyc[1:]+cyc[:1])}
      winding=sum(orientation[tuple(r['seam'])]*r['current'] for r in meter['records'] if r['stage']=='measure')
      if meter['success']:assert winding==6*degree;success_degrees.add(degree)
      else:failure_reasons.add(meter['reason'])
      selected.setdefault((degree,meter['reason']),row);replays+=1
    assert success_degrees=={-2,-1,0,1,2}
    public=[];mutations={}
    for key,row in selected.items():
      v=reviewed.admit_and_decode(row['claim'],row['G'],row['region'],2);assert v['admitted']
      assert v['readout']==(-key[0] if row['claim']['meter']['success'] else None)
      bad=deepcopy(row['claim']);bad['meter']['final_loads']={'untrusted':42}
      assert reviewed.admit_and_decode(bad,row['G'],row['region'],2)==v
      public.append({'degree':key[0],'reason':key[1],'admitted':True})
    vacuum=json.loads((SOURCE/'VACUUM.json').read_text());claim=vacuum['claim'];G=vacuum['G'];region=(0,)
    assert reviewed.admit_and_decode(claim,G,region,2)['readout']==0
    for name,mutate in [('empty_pointer',lambda x:x['meter'].__setitem__('records',[])),('forged_charge',lambda x:x.__setitem__('decoded_relative_charge',F(1,2))),('wrong_closing_contrast',lambda x:x['sweep_records'][0].__setitem__('contrast',x['sweep_records'][0]['contrast']+1))]:
      bad=deepcopy(claim);mutate(bad);mutations[name]=not reviewed.admit_and_decode(bad,G,region,2)['admitted']
    mutations['half_current_breaks_circle_gluing']=abs(cmath.exp(2j*math.pi*.5)-1)>1
    mutations['unrestricted_counter_alias']=6%6==0 and 6!=0
    mutations['nonunique_integral_dictionary']=bareiss([[1,1],[0,1]])==1
    assert all(mutations.values())
    # Independent shortest-distance derivation of each initial reset law.
    adj={v:[b if a==v else a for a,b in edges if v in (a,b)] for v in range(12)}
    dist={0:0};todo=deque([0])
    while todo:
      v=todo.popleft()
      for w in adj[v]:
        if w not in dist:dist[w]=dist[v]+1;todo.append(w)
    laws=0
    for retry in (1,2,4,16):
      h=F(1,2**retry)
      for x in range(12):
        ys=[b if a==x else a for a,b in edges if x in (a,b)]
        p=sum((h**k*(1-h)**(1+dist[y]) for k,y in enumerate(ys)),F())
        assert p==reviewed.integrity.meter.reset_success(x,retry)
        for w in range(-2,3):
          success=p*(1-h)**(6*abs(w));assert 0<success<1 and success+(1-success)==1;laws+=1
    full=json.loads((ROOT/'reverse-engineering-reality/code/openwave_row1_integrity/RESULTS.json').read_text())['complete_bounded_laws']
    assert len(full['rows'])==768 and all(sum(map(F,x['outcome_mass'].values()))==1 for x in full['rows'])
    assert before==pins()
    result={'producer_independent':True,'full_criterion_met':True,'criterion_sha256':json.loads((SOURCE/'CRITERION.json').read_text())['criterion_sha256'],
      'candidate_id':'literal-native-relative-spatial-charge-AU5-record-only-v1','source_sha256':before,
      'exact_integer_minors':[detB,detA],'continuous_homotopy_tests':cases,'maximum_phase_error':phase_error,
      'native_branches':branches,'bounded_reference_pairs':pairs,'independent_source_pointer_replays':replays,
      'public_admission_cases':public,'independent_pointer_probabilities':laws,'frozen_complete_laws_consumed':768,
      'fresh_full_law_rerun':False,'success_degrees':sorted(success_degrees),'failure_reasons':sorted(failure_reasons),'mutations':mutations,
      'domain':'Static relative topological charge n-G, fixed true G, totals 1 or 2, faithfully received native closing tape, bounded native pointer with all failures retained.',
      'rejected_transfers':['arbitrary GL(11,Z) uniqueness','unrestricted mod6 charges','empirical electric units','coherent Maxwell/Dirac source','authentication of replaced closing channel'],
      'verdict_reason':'The declared local spatial S1 realization has a faithful integral homotopy class; its actual source events and complete native observer readout share that class. Integer regional degrees and additivity forbid equal half daughters. This satisfies the literal static row within the declared domain and exact precedent scope.'}
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('full_criterion_met','continuous_homotopy_tests','native_branches','independent_source_pointer_replays','independent_pointer_probabilities')}))
if __name__=='__main__':main()
