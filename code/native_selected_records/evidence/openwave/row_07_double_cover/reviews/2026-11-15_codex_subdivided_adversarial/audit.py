"""Independent same-group scientific audit; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
import json,math,sys,hashlib
import numpy as np
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[2]/'work/2026-11-15_codex_subdivided_rotation'
receipt=sys.argv[2] if len(sys.argv)>2 else 'CONFIRM_RECEIPT.json'
stationary=receipt=='FIXED_REFERENCE_RECEIPT.json'
r=json.loads((BASE/receipt).read_text())
EDGES=((0,1),(1,2),(0,2));I=np.eye(3,dtype=int)
CHECKS=0

def need(ok,label):
 global CHECKS
 CHECKS+=1
 if not ok:raise AssertionError(label)

def endpoint(q,e,c,law):
 a,b=e;y=list(q);d=q[a]-q[b]
 j=(d-(2*c-1)*(abs(d)%2))//2 if law=='native' else d
 y[a]-=j;y[b]+=j
 return y,j

def fine(edges,law,blank=False):
 maps=[]
 for e in edges:
  for c in (0,1):
   block=np.zeros((3,3,3),dtype=int)
   for p in range(3):
    out,j=endpoint([int(x==p) for x in range(3)],e,c,law) if not blank else ([int(x==p) for x in range(3)],0)
    block[j+1,out.index(1),p]=1
   maps.extend(block)
 return maps

def eraser(maps):
 count=len(maps)//3;prefix=np.zeros((3,3),dtype=int);out=[]
 for k,x in enumerate(maps):
  if k:out.append((prefix-k*x,count*k*(k+1)))
  prefix=prefix+x
 return [(prefix,count*len(maps))]+out

def contrast(edges,law):
 return [(b+s*a,4*d) for (a,d),(b,db) in zip(eraser(fine(edges,law)),eraser(fine(edges,law,True))) for s in (1,-1)]

def address(slot,m):
 v,k=divmod(slot,m);w=(v+1)%3
 return [e for e in EDGES if v in e]*(m-k)+[e for e in EDGES if w in e]*k

def mv(a,v):return tuple(sum(F(a[i,j])*v[j] for j in range(3)) for i in range(3))
def norm(v):return sum(x*x for x in v)
def inner(x,y):return sum(a*b for a,b in zip(x,y))

def xy(slot,m):
 verts=np.array([[1.,0.],[-.5,math.sqrt(3)/2],[-.5,-math.sqrt(3)/2]])
 v,k=divmod(slot%(3*m),m);return (1-k/m)*verts[v]+(k/m)*verts[(v+1)%3]

def delta(a,b):return math.atan2(a[0]*b[1]-a[1]*b[0],a@b)

# All-input preparation, including untouched-reference extension.
prep_bounds={};catalog={}
for law in ('native','swap'):
 warm=fine(EDGES,law);accepted,d=contrast((EDGES[0],),law)[1]
 effect=sum((x.T@accepted.T@accepted@x for x in warm),np.zeros((3,3),int))
 E=[[F(int(effect[i,j]),6*d) for j in range(3)] for i in range(3)]
 expected=[[F(5,144),F(-1,72),0],[F(-1,72),F(5,144),0],[0,0,F(1,72)]] if law=='native' else [[F(int(i==j),9) for j in range(3)] for i in range(3)]
 need(E==expected,'all-input preparation effect')
 prep_bounds[law]={'effect':[[str(x) for x in row] for row in E],'min_eigenvalue':str(F(1,72) if law=='native' else F(1,9))}
for key,c in r['instrument_catalog'].items():
 m,law=key.split(':');m=int(m);tables=[]
 for slot,encoded in enumerate(c['target']):
  ops=eraser(fine(address(slot,m),law))
  for (matrix,den),actual in zip(ops,encoded):
   need(matrix.tolist()==[[int(v) for v in row] for row in actual['numerator']] and den==actual['squared_denominator'],'source-derived catalog')
  need(len(ops)==len(encoded),'full catalog flags')
  total=sum(fine(address(slot,m),law));P=np.array([[F(int(x),4*m) for x in row] for row in total],dtype=object)
  tables.append((P,ops))
 plane=contrast(EDGES,law)
 for (matrix,den),actual in zip(plane,c['plane_filter']):
  need(matrix.tolist()==[[int(v) for v in row] for row in actual['numerator']] and den==actual['squared_denominator'],'plane filter source attachment')
 for ops in [*(x[1] for x in tables),plane]:
  gram=[[sum(F(int((a.T@a)[i,j]),d) for a,d in ops) for j in range(3)] for i in range(3)]
  need(gram==[[F(int(i==j)) for j in range(3)] for i in range(3)],'complete instrument')
 catalog[key]=(tables,plane)
need(len(r['cases'])==(4 if stationary else 16) and not r['gate_failures'],'complete confirmation grid')
expected_cases={(m,t,d,l) for m in((8,) if stationary else (4,8)) for t in(1,2) for d in((1,) if stationary else (-1,1)) for l in('native','swap')}
need({(c['M'],c['turns'],c['direction'],c['law']) for c in r['cases']}==expected_cases,'case menu')
results=[];flags=0;histories=0;frame_rows=0
basis=np.array([[1/math.sqrt(2),1/math.sqrt(6)],[-1/math.sqrt(2),1/math.sqrt(6)],[0,-2/math.sqrt(6)]])
for case in r['cases']:
 m,t,direction,law=case['M'],case['turns'],case['direction'],case['law'];tables,plane=catalog[f'{m}:{law}']
 positions={'rotating':2*m,'reference':2*m};angles={x:[] for x in positions};anisotropy={x:[] for x in positions}
 for arm in positions:
  T=basis.T@np.array(tables[positions[arm]][0],float)@basis;anisotropy[arm].append(np.array([T[0,0]-T[1,1],2*T[0,1]]))
 a=b=(F(1),F(-1),F(0));reach=F(5,144) if law=='native' else F(1,9);rejected=1-reach
 for row in case['records']:
  need(F(row['reach'])==reach,'reach before stage')
  step=row['step']
  if row['kind']=='marker':
   arm=row['arm'];old=positions[arm];edge=row['address'];success=F(0);destinations=set()
   need((edge[1]-edge[0]) in (1,3*m-1),'actual adjacent marker seam')
   for h in row['histories']:
    p=old;seen=False
    need(len(h['coins'])==len(h['currents']) and 1<=len(h['coins'])<=8,'capped marker tape')
    for j,(coin,current) in enumerate(zip(h['coins'],h['currents'])):
     need(coin in(0,1) and not seen,'native coin/stopping boundary')
     q=[int(x==p) for x in range(3*m)];out,actual=endpoint(q,edge,coin,law)
     need(actual==current,'marker current');p=out.index(1);seen=bool(current)
    need(p==h['decoded_position'],'record-only marker decoder')
    mass=F(1,2**len(h['coins']));need(mass==F(h['probability']),'marker word probability')
    if h['outcome']=='hop':need(seen,'actual moved marker');success+=mass;destinations.add(p)
    else:need(not seen and len(h['coins'])==8,'actual timeout')
    histories+=1
   need(len(destinations)==1,'unambiguous accepted marker location');new=destinations.pop()
   offset=min(step%m,m-step%m);request=(2*m+direction*(step if arm=='rotating' else offset))%(3*m)
   need(new==request,'same controller decoded route')
   ang=math.degrees(delta(xy(old,m),xy(new,m)));angles[arm].append(ang)
   need(abs(ang-row['angle_degrees'])<1e-11 and np.max(abs(xy(new,m)-row['decoded_coordinate']))<1e-12,'frame from actual positions')
   need(abs(ang)<90,'no coarse spatial alias');positions[arm]=new;frame_rows+=1
   rejected+=reach*(1-success);need(F(row['timeout'])==reach*(1-success),'all timeout mass');reach*=success
   continue
  need([[str(x) for x in a],[str(x) for x in b]]==row['audit_input_vectors'],'source continuation state')
  normalizer=max(abs(x) for x in a+b);aa=np.array([float(x/normalizer) for x in a]);bb=np.array([float(x/normalizer) for x in b]);denorm=aa@aa+bb@bb
  if row['kind']=='target':
   need(row['rotating_slot']==positions['rotating'] and row['reference_slot']==positions['reference'],'source acts at recorded frame')
   P,left=tables[positions['rotating']];Q,right=tables[positions['reference']]
   predicted=[((g@aa)@(g@aa)+(h@bb)@(h@bb))/(den*denorm) for (g,den),(h,den2) in zip(left,right)]
   if row['event']==0:
    for arm,mat in [('rotating',P),('reference',Q)]:
     T=basis.T@np.array(mat,float)@basis;anisotropy[arm].append(np.array([T[0,0]-T[1,1],2*T[0,1]]))
    # Full 3x3 spectra match: same trace and determinant on both arms.
    if not stationary:need(abs(np.linalg.det(np.array(P,float))-np.linalg.det(np.array(Q,float)))<1e-13,'matched instantaneous spectrum')
   outa,outb=mv(P,a),mv(Q,b);accept=(norm(outa)+norm(outb))/(3*(norm(a)+norm(b)));a,b=outa,outb;accepted=0
  else:
   predicted=[((g@aa)@(g@aa)+(g@bb)@(g@bb))/(den*denorm) for g,den in plane]
   accept=F(1,48 if law=='native' else 12);accepted=1
  actual=[float(F(x)) for x in row['flag_probabilities']]
  need(len(actual)==len(predicted) and np.max(abs(np.array(actual)-predicted))<1e-12,'every flag population')
  need(F(row['flag_probabilities'][accepted])==accept,'exact selected mass')
  need(F(row['first_failure'])==reach*(1-accept),'all first failures');rejected+=reach*(1-accept);reach*=accept;flags+=len(actual)
 visibility=2*inner(a,b)/(norm(a)+norm(b));need(visibility==F(case['exact_visibility']),'independent signed interference')
 need(reach==F(case['success']) and reach+rejected==1,'complete total outcome')
 wind={arm:sum(delta(x,y) for x,y in zip(v,v[1:]))/(2*math.pi) for arm,v in anisotropy.items()}
 need(abs(sum(angles['rotating'])-direction*360*t)<1e-10 and abs(sum(angles['reference']))<1e-10,'actual frame winding')
 need(abs(wind['rotating'])==t or abs(abs(wind['rotating'])-t)<1e-10,'source anisotropy winding')
 need(abs(wind['reference'])<1e-10,'reference anisotropy unwound')
 results.append({'M':m,'turns':t,'direction':direction,'law':law,'visibility':float(visibility),
  'frame_turn_degrees':sum(angles['rotating']),'max_frame_increment':max(map(abs,angles['rotating'])),
  'anisotropy_winding':wind,'maximum_anisotropy_increment_degrees':max(abs(math.degrees(delta(x,y))) for x,y in zip(anisotropy['rotating'],anisotropy['rotating'][1:])),
  'log10_success':math.log10(reach.numerator)-math.log10(reach.denominator)})
print(json.dumps({'receipt':receipt,'verdict':'PASS','checks':CHECKS,'preparation_all_input':prep_bounds,'all_outcome_flags':flags,'marker_histories':histories,'actual_frame_steps':frame_rows,'cases':results,'same_group_only':True,'producer_imports':False},indent=2))
