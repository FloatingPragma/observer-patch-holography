"""Conditional driven relative-count state from atomic record preparation.

A uniform signed record event followed by complete native repair yields a
unique stationary law on relative heights. The preparation and readout
projection are additional hypotheses, not a source-selected physical vacuum.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as sp
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SPEC=HERE/'native_record_renewal_spec.json'
OUTPUT=HERE.parent/'runs/flavor/native_record_renewal.json'
PINS=(
 'code/a5_closure/manifests/echosahedral_federation_reference.json',
 'code/a5_closure/manifests/record_counting_mechanism_reference.json',
 'code/a5_closure/manifests/a3_scheduler_kernel_reference.json',
 'code/a5_closure/source_repair_generator_certificate.py',
 'code/a5_closure/record_counting_mechanism_certificate.py',
 'code/particles/flavor/native_record_renewal_spec.json',
 'code/particles/flavor/native_record_renewal.py',
)
sys.path.insert(0,str(ROOT/'code/a5_closure'))
import source_repair_generator_certificate as source

def build_payload():
 spec=json.loads(SPEC.read_text())
 record=json.loads((ROOT/PINS[1]).read_text())['mechanism']
 scheduler=json.loads((ROOT/PINS[2]).read_text())['a3_selection']
 assert record['repair_move']=='conservative unit transfer across a seam with difference magnitude at least two'
 assert record['event_grammar']=='append-only atomic signed events with checkpoints'
 assert scheduler['reference']=='uniform_on_moves' and scheduler['selected_kernel_probability_per_seam']=='1/30'
 carrier_path=ROOT/'code/a5_closure/manifests/echosahedral_federation_reference.json'
 carrier=json.loads(carrier_path.read_text())['carrier'];labels=carrier['ports']
 edges=sorted(tuple(sorted((labels.index(i),labels.index(j)))) for i,j in carrier['edges'])
 adj=[sorted(j if i==v else i for i,j in edges if v in (i,j)) for v in range(12)]
 d=np.full((12,12),13,dtype=np.int64);np.fill_diagonal(d,0)
 for i,j in edges:d[i,j]=d[j,i]=1
 for k in range(12):d=np.minimum(d,d[:,k,None]+d[None,k,:])
 faces=[tuple(labels.index(v) for v in f) for f in carrier['oriented_faces']]
 positive={f[i:]+f[:i] for f in faces for i in range(3)}
 group=source.proper_autos(12,edges,faces)
 assert len(group)==60
 # Anchor one height and traverse graph constraints, rather than enumerate a box.
 order=sorted(range(12),key=lambda i:(d[0,i],i));assigned={0:0};height_states=set()
 def visit(k):
  if k==12:
   m=min(assigned.values());height_states.add(tuple(assigned[i]-m for i in range(12)));return
  i=order[k];known=[assigned[j] for j in adj[i] if j in assigned]
  for value in range(max([-int(d[0,i])]+[x-1 for x in known]),min([int(d[0,i])]+[x+1 for x in known])+1):
   assigned[i]=value;visit(k+1)
  assigned.pop(i,None)
 visit(1);states=sorted(height_states);assert len(states)==6077
 normalize=lambda z:tuple(x-min(z) for x in z)
 def orbit(z):return {tuple(z[p[i]] for i in range(12)) for p in group}
 index={};reps=[];sizes=[]
 for z in states:
  if z in index:continue
  members=orbit(z);a=len(reps);reps.append(z);sizes.append(len(members))
  for y in members:index[y]=a
 assert len(reps)==136

 def response(z):
  """Full endpoint law; downhill defect routing proved from native transfer."""
  @lru_cache(None)
  def endpoints(v,sign):
   legal=[j for j in adj[v] if sign*(z[v]-z[j])==1]
   if not legal:return {v:F(1)}
   result=Counter()
   for j in legal:
    for end,weight in endpoints(j,sign).items():result[end]+=weight/len(legal)
   assert sum(result.values())==1
   return dict(result)
  row=Counter();max_cost=0;expected_cost=F(0)
  for sign in [-1,1]:
   for start in range(12):
    for end,weight in endpoints(start,sign).items():
     y=list(z);y[end]+=sign;y=normalize(y)
     assert y in index and all(abs(y[i]-y[j])<=1 for i,j in edges)
     row[y]+=weight/24
     moves=sign*(z[start]-z[end]);assert 0<=moves<=3
     max_cost=max(max_cost,moves);expected_cost+=weight*moves/24
  assert sum(row.values())==1
  return row,max_cost,expected_cost

 rows=[];max_moves=0;cost=[];full_representative_rows=[]
 for a,z in enumerate(reps):
  full,m,mean=response(z);row=Counter()
  for y,p in full.items():row[index[y]]+=p
  rows.append(row);full_representative_rows.append(full);max_moves=max(max_moves,m);cost.append(mean)

 # Exact invariant distribution in the proper-rotation quotient.
 A=sp.zeros(136);rhs=sp.zeros(136,1)
 for a,row in enumerate(rows):
  for b,p in row.items():A[b,a]+=sp.Rational(p.numerator,p.denominator)
  A[a,a]-=1
 A[135,:]=sp.ones(1,136);rhs[135]=1
 solution=A.inv(method='DM')*rhs
 assert all(x>0 for x in solution) and sum(solution)==1
 assert all(sum(solution[a]*sp.Rational(rows[a].get(b,0)) for a in range(136))==solution[b] for b in range(136))

 # Independent finite direct-step connectivity witness and period.
 for z in states:
  current=list(z)
  while max(current):
   i=current.index(max(current));before=tuple(current);current[i]-=1
   assert all(abs(current[a]-current[b])<=1 for a,b in edges)
   assert normalize(current) in index
   assert all(before[i]>=before[j] for j in adj[i])
   assert all(current[i]<=current[j] for j in adj[i])
 assert all((sum(reps[b])-sum(reps[a]))%2==1 for a,row in enumerate(rows) for b in row)

 # Exact spectral weights using integer30Q matrices (not mass eigenvalues).
 sqrt5=sp.sqrt(5);phi=(1+sqrt5)/2
 vertices=[sp.Matrix(v) for v in [(0,1,phi),(0,1,-phi),(0,-1,phi),(0,-1,-phi),(1,phi,0),(1,-phi,0),(-1,phi,0),(-1,-phi,0),(phi,0,1),(phi,0,-1),(-phi,0,1),(-phi,0,-1)]]
 assert all(sp.simplify(v.dot(v)-(phi+2))==0 for v in vertices)
 assert [(i,j) for i in range(12) for j in range(i+1,12)
         if sp.simplify((vertices[i]-vertices[j]).dot(vertices[i]-vertices[j]))==4]==edges
 assert all((vertices[i]==-vertices[j])==(d[i,j]==3) for i in range(12) for j in range(12))
 basis=[(v*v.T/(phi+2)-sp.eye(3)/3).applyfunc(sp.simplify) for v in vertices]
 field=sp.QQ.algebraic_field(sqrt5);zero=field.zero
 basisf=[[field.from_sympy(x) for x in q] for q in basis]
 spectra=Counter();residues=Counter();q2mean=sp.S.Zero;q4radial=sp.S.Zero;oddmean=sp.S.Zero;cov=[sp.S.Zero]*4
 orbitrows=[]
 for a,z in enumerate(reps):
  q=[sum((z[i]*basisf[i][j] for i in range(12)),zero) for j in range(9)]
  square=[sum((q[3*i+k]*q[3*k+j] for k in range(3)),zero) for i in range(3) for j in range(3)]
  t2=field.to_sympy(square[0]+square[4]+square[8]);t3=sp.expand(field.to_sympy(sum((square[3*i+k]*q[3*k+i] for i in range(3) for k in range(3)),zero)))
  key=';'.join(map(str,(t2,t3.coeff(sqrt5,0),t3.coeff(sqrt5))))
  p=solution[a];spectra[key]+=p;residues[str(sum(z)%12)]+=p;q2mean+=p*t2;q4radial+=p*t2*t2
  centered=[sp.Rational(12*x-sum(z),12) for x in z]
  for distance in range(4):
   pairs_at_distance=[(i,j) for i in range(12) for j in range(12) if d[i,j]==distance]
   cov[distance]+=p*sum(centered[i]*centered[j] for i,j in pairs_at_distance)/len(pairs_at_distance)
  pairs=[(i,j) for i in range(12) for j in range(i+1,12) if d[i,j]==3]
  oddnorm=sp.Rational(sum((z[i]-z[j])**2 for i,j in pairs),2);oddmean+=p*oddnorm
  orbitrows.append({'state':list(z),'orbit_size':sizes[a],'stationary_weight':str(p),'per_state_weight':str(p/sizes[a]),'spectral_key':key,'mean_accepted_moves':str(cost[a])})
 assert set(residues.values())=={sp.Rational(1,12)}
 # Charge conjugation pairs have equal exact stationary weights.
 assert all(solution[a]==solution[index[normalize(tuple(-x for x in z))]] for a,z in enumerate(reps))
 mean_cost=sum(solution[a]*sp.Rational(cost[a]) for a in range(136))
 assert mean_cost==sp.Rational(11,24)
 q2mean=sp.factor(q2mean);oddmean=sp.factor(oddmean)
 gaussian_defect=sp.factor(q4radial-sp.Rational(7,5)*q2mean*q2mean)
 assert gaussian_defect!=0
 current=None
 for a,z in enumerate(reps):
  for y,p in sorted(full_representative_rows[a].items()):
   b=index[y];reverse=response(y)[0].get(z,F(0))
   forward_flow=solution[a]*sp.Rational(p)/sizes[a]
   reverse_flow=solution[b]*sp.Rational(reverse)/sizes[b]
   if forward_flow!=reverse_flow:
    current={'from':list(z),'to':list(y),'forward_probability':str(p),
             'reverse_probability':str(reverse),'forward_stationary_flow':str(forward_flow),
             'reverse_stationary_flow':str(reverse_flow),'net_current':str(forward_flow-reverse_flow)}
    break
  if current:break
 assert current is not None
 transitions=[{str(b):str(p) for b,p in sorted(row.items())} for row in rows]
 b3=sp.expand(cov[0]+sqrt5*(cov[1]-cov[2])-cov[3]);b3p=sp.expand(cov[0]-sqrt5*(cov[1]-cov[2])-cov[3]);b5=sp.factor(cov[0]-cov[1]-cov[2]+cov[3])
 assert sp.simplify(q2mean-8*b5)==0 and sp.simplify(oddmean-3*(b3+b3p))==0
 out={'schema':'oph.native_record_renewal.v1','specification':spec,'scope':spec['scope'],'source_pins':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},'complete_height_sha256':hashlib.sha256(b''.join(bytes(z) for z in states)).hexdigest(),'stationary_port_covariance':{'by_distance':list(map(str,cov)),'P3':str(b3),'P3prime':str(b3p),'P5':str(b5)},'mean_trace_Q3':'0','centered_count_odd_moments':'zero' ,'carrier_sha256':hashlib.sha256(carrier_path.read_bytes()).hexdigest(),
 'states':6077,'A5_orbits':136,'quotient_nonzero_entries':sum(map(len,rows)),
 'full_kernel_nonzero_entries':sum(sizes[a]*len(row) for a,row in enumerate(full_representative_rows)),
 'irreducible':True,'period':2,'max_accepted_moves_per_preparation':max_moves,
 'detailed_balance':False,'stationary_current_witness':current,
 'fourth_moment':{'mean_squared_trace_Q2':str(q4radial),'isotropic_Gaussian_radial_defect':str(gaussian_defect)},
 'stationary_residue_weights':{k:str(p) for k,p in sorted(residues.items())},
 'stationary_Q_spectral_weights':{k:str(p) for k,p in sorted(spectra.items())},
 'stationary_mean_Q':'0','stationary_mean_trace_Q2':str(q2mean),'stationary_quadrupole_covariance_scalar':str(q2mean/5),
 'stationary_odd_port_square_mean':str(oddmean),'stationary_mean_accepted_moves_per_preparation':str(mean_cost),
 'numeric':{'mean_trace_Q2':float(q2mean),'Q_covariance_scalar':float(q2mean/5),'odd_port_square_mean':float(oddmean),'mean_accepted_moves':float(mean_cost),'mean_squared_trace_Q2':float(q4radial),'isotropic_Gaussian_radial_defect':float(gaussian_defect)},
 'orbits':orbitrows,'orbit_kernel':transitions}
 return out

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
 rendered=json.dumps(build_payload(),indent=2,sort_keys=True)+'\n'
 if args.check:
  if OUTPUT.read_text()!=rendered:raise SystemExit('Native record renewal differs; regenerate')
  print('Native record renewal verified')
 else:
  OUTPUT.write_text(rendered);print(OUTPUT)
if __name__=='__main__':main()
