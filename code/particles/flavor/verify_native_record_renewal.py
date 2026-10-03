"""Independent full-state replay of the conditional native record renewal.

No producer or source Python module is imported. A ternary-box census and
NetworkX incidence symmetries replace the producer's anchored traversal;
actual native transfers replace its downhill-defect routing. The supplied
stationary weights are certified on every state, not accepted from a solver.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import networkx as nx
import numpy as np
import sympy as sp
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SPEC=HERE/'native_record_renewal_spec.json'
OUTPUT=HERE.parent/'runs/flavor/native_record_renewal.json'
REVIEWED_SPEC_SHA256='ac910d0b2f52000e4611bc36a10b68abe03a894bb310fcc18ba7765a4f519f50'
PINS=(
 'code/a5_closure/manifests/echosahedral_federation_reference.json',
 'code/a5_closure/manifests/record_counting_mechanism_reference.json',
 'code/a5_closure/manifests/a3_scheduler_kernel_reference.json',
 'code/a5_closure/source_repair_generator_certificate.py',
 'code/a5_closure/record_counting_mechanism_certificate.py',
 'code/particles/flavor/native_record_renewal_spec.json',
 'code/particles/flavor/native_record_renewal.py',
)
def require(condition,message):
 if not condition:raise ValueError(message)
def strict_load(path):
 def pairs(items):
  result={}
  for k,v in items:
   require(k not in result,'duplicate JSON key');result[k]=v
  return result
 def nonfinite(value):raise ValueError('nonfinite JSON value')
 return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=nonfinite)
def serialized(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def reviewed_spec():
 require(hashlib.sha256(SPEC.read_bytes()).hexdigest()==REVIEWED_SPEC_SHA256,'unreviewed renewal specification')
 return strict_load(SPEC)
def rational(value):
 require(type(value) is str,'rational fields must be exact strings')
 try:q=F(value)
 except (ValueError,ZeroDivisionError) as e:raise ValueError('invalid rational') from e
 require(str(q)==value,'noncanonical rational');return q

def normal(z):return tuple(x-min(z) for x in z)

@lru_cache(None)
def independent_data():
 spec=reviewed_spec();carrier=strict_load(ROOT/spec['carrier'])['carrier'];labels=carrier['ports']
 record=strict_load(ROOT/PINS[1])['mechanism'];scheduler=strict_load(ROOT/PINS[2])['a3_selection']
 require(record['repair_move']=='conservative unit transfer across a seam with difference magnitude at least two'
         and record['event_grammar']=='append-only atomic signed events with checkpoints','source integer grammar')
 require(scheduler['reference']=='uniform_on_moves' and scheduler['selected_kernel_probability_per_seam']=='1/30','source seam reference')
 edges=sorted(tuple(sorted((labels.index(i),labels.index(j)))) for i,j in carrier['edges'])
 graph=nx.Graph();graph.add_nodes_from(range(12));graph.add_edges_from(edges)
 distance=dict(nx.all_pairs_shortest_path_length(graph))
 require(all(Counter(distance[i].values())=={0:1,1:5,2:5,3:1} for i in range(12)),'icosahedral distances')
 faces=[tuple(labels.index(i) for i in f) for f in carrier['oriented_faces']]
 oriented={f[i:]+f[:i] for f in faces for i in range(3)}
 group=[tuple(p[i] for i in range(12)) for p in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter()
        if all(tuple(p[i] for i in f) in oriented for f in faces)]
 require(len(group)==60,'proper rotation group')
 raw=np.array(list(product(range(3),repeat=12)),dtype=np.int8);keep=raw.min(axis=1)==0
 for i,j in edges:keep &= np.abs(raw[:,i]-raw[:,j])<=1
 cones=[tuple(3-distance[i][j] for j in range(12)) for i in range(12)]
 states=sorted([tuple(map(int,z)) for z in raw[keep]]+cones);require(len(states)==6077,'global count classes')
 state_index={z:i for i,z in enumerate(states)};index={};reps=[];sizes=[]
 for z in states:
  if z in index:continue
  members={tuple(z[p[i]] for i in range(12)) for p in group};a=len(reps);reps.append(z);sizes.append(len(members))
  for y in members:index[y]=a
 require(len(reps)==136,'global A5 orbits')

 @lru_cache(None)
 def settle(z):
  legal=[(i,j) if z[i]>z[j] else (j,i) for i,j in edges if abs(z[i]-z[j])>=2]
  if not legal:return ((normal(z),F(1)),),0,F(0)
  result=Counter();longest=0;mean=F(1)
  for i,j in legal:
   require(z[i]-z[j]==2,'one-event repairs have difference exactly2')
   y=list(z);y[i]-=1;y[j]+=1;y=tuple(y)
   require(sum(y)==sum(z) and sum(t*t for t in z)-sum(t*t for t in y)==2,'strict native unit descent')
   terminal,depth,cost=settle(y);longest=max(longest,depth+1);mean+=cost/len(legal)
   for end,p in terminal:result[end]+=p/len(legal)
  return tuple(sorted(result.items())),longest,mean

 full_rows=[];cost=[];max_moves=0;orbit_rows=[None]*136
 for z in states:
  row=Counter();mean=F(0)
  for sign in (-1,1):
   for port in range(12):
    kicked=list(z);kicked[port]+=sign
    terminal,depth,steps=settle(tuple(kicked));max_moves=max(max_moves,depth);mean+=steps/24
    for y,p in terminal:row[state_index[y]]+=p/24
  require(sum(row.values())==1,'full kernel row normalization')
  require(all((sum(states[j])-sum(z))%2==1 for j in row),'alternating total parity')
  full_rows.append(dict(row));cost.append(mean)
  aggregated=Counter()
  for j,p in row.items():aggregated[index[states[j]]]+=p
  a=index[z]
  if orbit_rows[a] is None:orbit_rows[a]=dict(aggregated)
  else:require(orbit_rows[a]==dict(aggregated),'all-state strong lumpability')
 require(max_moves==3,'three-step sharp bound')
 # Source-independent constructive irreducibility certificate. These direct
 # event paths are bidirectional and have probability at least1/24 per edge.
 for z in states:
  current=list(z)
  while max(current):
   i=current.index(max(current));before=normal(current);current[i]-=1;after=normal(current)
   require(all(abs(current[a]-current[b])<=1 for a,b in edges),'settled maximum removal')
   require(full_rows[state_index[before]].get(state_index[after],0)>=F(1,24),'forward direct path')
   require(full_rows[state_index[after]].get(state_index[before],0)>=F(1,24),'reverse direct path')
 require(any(full_rows[j].get(i,0)>0 for i,row in enumerate(full_rows) for j in row),'two-cycle')

 # Independent integer spectral evaluation,30Q=A+sqrt(5)B.
 sqrt5=sp.sqrt(5);phi=(1+sqrt5)/2
 vv=[sp.Matrix(v) for v in [(0,1,phi),(0,1,-phi),(0,-1,phi),(0,-1,-phi),(1,phi,0),(1,-phi,0),(-1,phi,0),(-1,-phi,0),(phi,0,1),(phi,0,-1),(-phi,0,1),(-phi,0,-1)]]
 require(all(sp.simplify(v.dot(v)-(phi+2))==0 for v in vv),'unit frame normalization')
 geometric_edges=[]
 for i in range(12):
  for j in range(i+1,12):
   square=sp.expand(sum((vv[i][k]-vv[j][k])**2 for k in range(3)))
   if sp.simplify(square-4)==0:geometric_edges.append((i,j))
   require((vv[i]==-vv[j])==(distance[i][j]==3),'frame antipodes match graph')
 require(geometric_edges==edges,'vertex labels match pinned nearest-neighbor edges')
 aa=[];bb=[]
 for v in vv:
  q=(30*(v*v.T/(phi+2)-sp.eye(3)/3)).applyfunc(sp.simplify);ar=[];br=[]
  for x in q:
   x=sp.expand(x);b=x.coeff(sqrt5);a=sp.simplify(x-b*sqrt5)
   require(a.is_Integer and b.is_Integer,'integer quadratic-field coefficients');ar.append(int(a));br.append(int(b))
  aa.append(ar);bb.append(br)
 heights=np.array(states,dtype=np.int64);A=(heights@np.array(aa)).reshape((-1,3,3));B=(heights@np.array(bb)).reshape((-1,3,3))
 require(max(abs(A).max(),abs(B).max())<1000,'exact integer arithmetic bound')
 tr2=np.einsum('nij,nji->n',A,A)+5*np.einsum('nij,nji->n',B,B)
 require(np.all(np.einsum('nij,nji->n',A,B)==0),'rational quadratic invariant')
 tr3a=np.einsum('nij,njk,nki->n',A,A,A)+15*np.einsum('nij,njk,nki->n',A,B,B)
 tr3b=3*np.einsum('nij,njk,nki->n',A,A,B)+5*np.einsum('nij,njk,nki->n',B,B,B)
 invariants=[(F(int(tr2[i]),900),F(int(tr3a[i]),27000),F(int(tr3b[i]),27000)) for i in range(len(states))]
 require(len(set(invariants))==31,'individual spectral menu')
 return {'states':states,'state_index':state_index,'orbits':reps,'sizes':sizes,'index':index,'kernel':full_rows,
         'orbit_kernel':orbit_rows,'cost':cost,'max_moves':max_moves,'invariants':invariants,'distance':distance,
         'digest':hashlib.sha256(b''.join(bytes(z) for z in states)).hexdigest()}

@lru_cache(None)
def independently_verified_payload(weights):
 d=independent_data();states=d['states'];reps=d['orbits'];sizes=d['sizes'];index=d['index'];state_index=d['state_index'];rows=d['orbit_kernel']
 require(len(weights)==136 and all(x>0 for x in weights) and sum(weights)==1,'positive normalized orbit weights')
 require(all(sum(weights[a]*rows[a].get(b,0) for a in range(136))==weights[b] for b in range(136)),'exact quotient stationarity')
 per_state=[weights[index[z]]/sizes[index[z]] for z in states]
 incoming=[F(0)]*len(states)
 for i,row in enumerate(d['kernel']):
  for j,p in row.items():incoming[j]+=per_state[i]*p
 require(incoming==per_state,'exact stationarity on all6077 states')
 require(all(weights[a]==weights[index[normal(tuple(-x for x in z))]] for a,z in enumerate(reps)),'sign symmetry')
 spectra=Counter();residues=Counter();t2mean=F(0);t4radial=F(0);oddmean=F(0);costmean=F(0);cov=[F(0)]*4;orbitrows=[]
 pairs=[(i,j) for i in range(12) for j in range(i+1,12) if d['distance'][i][j]==3]
 for a,z in enumerate(reps):
  i=state_index[z];p=weights[a];t2,t3a,t3b=d['invariants'][i];key=';'.join(map(str,(t2,t3a,t3b)))
  spectra[key]+=p;residues[str(sum(z)%12)]+=p;t2mean+=p*t2;t4radial+=p*t2*t2;costmean+=p*d['cost'][i]
  oddmean+=p*F(sum((z[j]-z[k])**2 for j,k in pairs),2)
  centered=[F(12*x-sum(z),12) for x in z]
  for distance in range(4):
   pairs_at_distance=[(j,k) for j in range(12) for k in range(12) if d['distance'][j][k]==distance]
   cov[distance]+=p*sum(centered[j]*centered[k] for j,k in pairs_at_distance)/len(pairs_at_distance)
  orbitrows.append({'state':list(z),'orbit_size':sizes[a],'stationary_weight':str(p),'per_state_weight':str(p/sizes[a]),'spectral_key':key,'mean_accepted_moves':str(d['cost'][i])})
 require(set(residues.values())=={F(1,12)} and costmean==F(11,24),'residue and stationarity balance')
 gaussian_defect=t4radial-F(7,5)*t2mean*t2mean;require(gaussian_defect!=0,'non-Gaussian radial fourth moment')
 sqrt5=sp.sqrt(5);b3=sp.expand(sp.Rational(cov[0]-cov[3])+sqrt5*sp.Rational(cov[1]-cov[2]));b3p=sp.expand(sp.Rational(cov[0]-cov[3])-sqrt5*sp.Rational(cov[1]-cov[2]));b5=cov[0]-cov[1]-cov[2]+cov[3]
 require(t2mean==8*b5 and oddmean==3*(2*(cov[0]-cov[3])),'port/quadrupole covariance normalization')
 current=None
 for a,z in enumerate(reps):
  i=state_index[z]
  for j,p in sorted(d['kernel'][i].items(),key=lambda x:states[x[0]]):
   reverse=d['kernel'][j].get(i,F(0));forward_flow=per_state[i]*p;reverse_flow=per_state[j]*reverse
   if forward_flow!=reverse_flow:
    current={'from':list(z),'to':list(states[j]),'forward_probability':str(p),'reverse_probability':str(reverse),
             'forward_stationary_flow':str(forward_flow),'reverse_stationary_flow':str(reverse_flow),'net_current':str(forward_flow-reverse_flow)};break
  if current:break
 require(current is not None,'nonzero stationary current')
 spec=reviewed_spec()
 return {'schema':'oph.native_record_renewal.v1','specification':spec,'scope':spec['scope'],
  'source_pins':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},'complete_height_sha256':d['digest'],
  'carrier_sha256':hashlib.sha256((ROOT/spec['carrier']).read_bytes()).hexdigest(),
  'states':6077,'A5_orbits':136,'quotient_nonzero_entries':sum(map(len,rows)),'full_kernel_nonzero_entries':sum(map(len,d['kernel'])),
  'irreducible':True,'period':2,'max_accepted_moves_per_preparation':d['max_moves'],
  'detailed_balance':False,'stationary_current_witness':current,
  'fourth_moment':{'mean_squared_trace_Q2':str(t4radial),'isotropic_Gaussian_radial_defect':str(gaussian_defect)},
  'stationary_residue_weights':{k:str(p) for k,p in sorted(residues.items())},'stationary_Q_spectral_weights':{k:str(p) for k,p in sorted(spectra.items())},
  'stationary_mean_Q':'0','mean_trace_Q3':'0','centered_count_odd_moments':'zero',
  'stationary_mean_trace_Q2':str(t2mean),'stationary_quadrupole_covariance_scalar':str(t2mean/5),
  'stationary_odd_port_square_mean':str(oddmean),'stationary_mean_accepted_moves_per_preparation':str(costmean),
  'stationary_port_covariance':{'by_distance':list(map(str,cov)),'P3':str(b3),'P3prime':str(b3p),'P5':str(b5)},
  'numeric':{'mean_trace_Q2':float(t2mean),'Q_covariance_scalar':float(t2mean/5),'odd_port_square_mean':float(oddmean),'mean_accepted_moves':float(costmean),'mean_squared_trace_Q2':float(t4radial),'isotropic_Gaussian_radial_defect':float(gaussian_defect)},
  'orbits':orbitrows,'orbit_kernel':[{str(b):str(p) for b,p in sorted(row.items())} for row in rows]}

def verify_payload(payload):
 spec=reviewed_spec()
 require(payload.get('specification')==spec and payload.get('scope')==spec['scope'],'renewal scope or specification drift')
 require(payload.get('source_pins')=={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},'live source pin drift')
 try:weights=tuple(rational(row['stationary_weight']) for row in payload['orbits'])
 except (KeyError,TypeError) as e:raise ValueError('orbit weights missing') from e
 expected=independently_verified_payload(weights)
 require(serialized(payload)==serialized(expected),'native record renewal receipt mismatch')
 return {'verified':True,'independent_states':6077,'independent_A5_orbits':136,
         'unique_stationary_relative_state':True,'period':2,'detailed_balance':False,
         'mean_accepted_repairs':'11/24','source_selected_physical_vacuum':False}

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('receipt',nargs='?',type=Path,default=OUTPUT);args=parser.parse_args()
 print(json.dumps(verify_payload(strict_load(args.receipt)),indent=2))
if __name__=='__main__':main()
