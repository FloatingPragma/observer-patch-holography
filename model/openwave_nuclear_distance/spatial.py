"""Full native three-quark singlets on two-, three- and four-cell graphs."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,combinations_with_replacement,permutations,product
from math import comb,sqrt
from collections import Counter
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_native_composites.spectrum import one_body,native
from model.openwave_axial_state.source import parameters,INPUTS

@lru_cache(None)
def carrier(L):
 orbitals=list(product(range(L),range(2),range(2)));idx={v:i for i,v in enumerate(orbitals)};labels=list(combinations_with_replacement(range(4*L),3));lookup={v:i for i,v in enumerate(labels)};counts=np.array([[lab.count(i)for i in range(4*L)]for lab in labels],dtype=int)
 def E(i,j):
  rr=[];cc=[];vv=[]
  for k,row in enumerate(counts):
   if not row[j]:continue
   out=row.copy();out[j]-=1;factor=sqrt(row[j]*(out[i]+1));out[i]+=1;lab=tuple(a for a,c in enumerate(out)for _ in range(c));rr.append(lookup[lab]);cc.append(k);vv.append(factor)
  return sp.csr_matrix((vv,(rr,cc)),shape=(len(labels),len(labels)))
 def angular(cells,spin):
  R=sp.csr_matrix((len(labels),len(labels)));z=np.zeros(len(labels))
  for x,other in product(cells,range(2)):
   a=(x,0,other)if spin else(x,other,0);b=(x,1,other)if spin else(x,other,1);R+=E(idx[a],idx[b]);z+=counts[:,idx[a]]-counts[:,idx[b]]
  return sp.diags(z*(z+2)/4)+R.T@R
 n=np.array([counts[:,[i for i,v in enumerate(orbitals)if v[0]==x]].sum(axis=1)for x in range(L)]).T
 # Explicit isometry to the full native occupation carrier.
 full=tuple(sum(1<<i for i in c)for c in combinations(range(12*L),3));findex={v:i for i,v in enumerate(full)};rr=[];cc=[];vv=[]
 for j,lab in enumerate(labels):
  assignments=sorted(set(permutations(lab)));den=sqrt(len(assignments))
  for assignment in assignments:
   occupied=[3*a+c for c,a in enumerate(assignment)];mask=sum(1<<i for i in occupied);sign=(-1)**sum(occupied[i]>occupied[k]for i in range(3)for k in range(i+1,3));rr.append(findex[mask]);cc.append(j);vv.append(sign/den)
 B=sp.csr_matrix((vv,(rr,cc)),shape=(len(full),len(labels)));assert np.max(abs((B.T@B-sp.eye(len(labels))).data),initial=0)<1e-14
 # Verify every native global colour current annihilates this explicit basis.
 qg=native()[0]
 colour_error=0.
 for a in range(8):
  A=np.zeros((12*L,12*L),complex)
  for x,s,f,c,d in product(range(L),range(2),range(2),range(3),range(3)):A[3*idx[x,s,f]+c,3*idx[x,s,f]+d]=qg[f][a][c,d]
  D=one_body(A,full)@B;colour_error=max(colour_error,float(np.max(abs(D.data),initial=0)))
 assert colour_error<2e-14
 # Check one nontrivial offsite and one spin generator against the actual CAR action.
 errs=[]
 for i,j in((idx[0,0,0],idx[L-1,1,1]),(idx[0,0,0],idx[0,1,0])):
  A=np.zeros((12*L,12*L))
  for c in range(3):A[3*i+c,3*j+c]=1
  D=one_body(A,full)@B-B@E(i,j);errs.append(float(np.max(abs(D.data),initial=0)))
 assert max(errs)<2e-14
 return dict(L=L,orbitals=orbitals,index=idx,labels=labels,counts=counts,n=n,E=E,angular=angular,full=full,B=B,colour_error=colour_error,generator_errors=errs)


def graph(L,edges,kappa):
 d=carrier(L);p=parameters();N=len(d['labels']);n=d['n'];idx=d['index'];E=d['E'];angular=d['angular'];counts=d['counts'];K=sp.csr_matrix((N,N));edge_colour=np.zeros(N)
 H=sp.diags(3*float(p['mbar'])+float(p['delta'])/2*(3-2*counts[:,[i for i,o in enumerate(d['orbitals'])if o[2]==0]].sum(axis=1))+float(p['A'])/9*np.sum(n*n,axis=1)+float(p['U'])*2/3*np.sum(n*(3-n),axis=1))
 for x in range(L):H+=float(p['Js'])*angular((x,),True)+float(p['Jw'])*angular((x,),False)
 for x,y in edges:
  for s,f in product(range(2),repeat=2):
   i,j=idx[x,s,f],idx[y,s,f];K+=E(i,i)+E(j,j)-E(i,j)-E(j,i)
  H+=float(p['eta'])*angular((x,y),False);nn=n[:,x]+n[:,y];edge_colour+=2/3*nn*(3-nn)
 H+=float(p['t'])*K+sp.diags(kappa*edge_colour)
 all_levels=eigh(H.toarray(),eigvals_only=True)
 # Complete fixed-charge/spin-projection block, then the exact S=I=1/2 spectral subspace.
 nu=counts[:,[i for i,o in enumerate(d['orbitals'])if o[2]==0]].sum(axis=1);ns=counts[:,[i for i,o in enumerate(d['orbitals'])if o[1]==0]].sum(axis=1);ids=np.flatnonzero((nu==1)&(ns==2));S=angular(tuple(range(L)),True).toarray()[np.ix_(ids,ids)];I=angular(tuple(range(L)),False).toarray()[np.ix_(ids,ids)];eye=np.eye(len(ids));F=(S-.75*eye)@(S-.75*eye)+(I-.75*eye)@(I-.75*eye);ef,V=eigh(F);P=V[:,ef<1e-9];h=H.toarray()[np.ix_(ids,ids)];e,v=eigh(P.T@h@P);state=np.zeros(N);state[ids]=P@v[:,0];prob=state*state;assert abs(prob.sum()-1)<2e-13
 Ip=sp.csr_matrix((N,N));Jp=Ip.copy()
 for x,ss in product(range(L),range(2)):
  term=E(idx[x,ss,0],idx[x,ss,1]);Ip+=term;Jp+=(1 if ss==0 else -1)*term
 native_gV=float(np.linalg.norm(Ip@state));native_gA=float(np.vdot(Ip@state,Jp@state).real)
 records=Counter()
 edge_set={tuple(sorted(e))for e in edges}
 for i,row in enumerate(n):
  occ=np.flatnonzero(row)
  key='colocated'if len(occ)==1 else'quark_diquark_adjacent'if len(occ)==2 and tuple(occ)in edge_set else'quark_diquark_nonadjacent'if len(occ)==2 else'three_occupied_cells'
  records[key]+=float(prob[i])
 # Explicit complete success/failure instrument on the uniform full native N3 carrier.
 raw=np.asarray(d['B']@state).ravel();fullN=len(raw);success=(abs(raw)**2/fullN);failure=(1-abs(raw)**2)/fullN;assert abs(success.sum()+failure.sum()-1)<2e-13
 return {'vertices':L,'edges':list(map(list,edges)),'edge_colour_coefficient':kappa,'full_N3_dimension':len(d['full']),'complete_colour_singlet_dimension':N,'all_colour_singlet_levels':all_levels.tolist(),'selected_spin_isospin_half_dimension':len(e),'selected_all_levels':e.tolist(),'selected_ground_GeV':float(e[0])*float(INPUTS['electron_unit_GeV']),'native_gV':native_gV,'native_gA':native_gA,'distance_class_records':dict(records),'native_current_isometry_errors':d['generator_errors'],'global_colour_annihilation_error':d['colour_error'],
  'original_record_instrument':{'preparation_success':1/fullN,'preparation_failure':1-1/fullN,'all_success_native_occupation_probabilities':success.tolist(),'all_failure_native_occupation_probabilities':failure.tolist(),'mean_work_on_uniform_input':0,'discarded_outcomes':0,'reference':'idle opposite-charge reference, unchanged by charge-preserving complete instruments'},
  'source_parameters_fixed':{k:str(v)for k,v in p.items()},'arithmetic':'Exact native-current isometry formula, floating generator checks and complete diagonalization diagnostics.'}

@lru_cache(None)
def spatial_certificate():
 rows=[]
 for L,edges in((2,((0,1),)),(3,((0,1),(1,2),(2,0))),(4,((0,1),(1,2),(2,3),(3,0)))):
  for k in(0,1):rows.append(graph(L,edges,k))
 p=parameters();return {'same_local_law_graphs':rows,'exact_colour_identity':'For native globally singlet3Q: C_X²=(2/3)n_X(3-n_X), for any cell subset X. This follows the antisymmetric colour epsilon and T_i dot T_j=-2/3.',
  'four_cycle_relative_potential':{'quark_diquark_adjacent_edge_cost':'8*kappa/3','quark_diquark_opposite_edge_cost':'16*kappa/3','adjacent_relative_attraction':'8*kappa/3','two_cell_edge_cost_for_all_singlets':'0','interpretation':'Actual finite-range diagonal colour-current attraction on the same four-cycle source; no disconnected-copy convention is used.'},
  'failed_local_colour_complement_shortcut':{'three_quark_coarse_gap':str(F(8,3)*p['U']+p['A']/3-(p['A']+F(15,4)*(p['Js']+p['Jw']+2*p['eta'])+6*p['t'])),'reason':'The co-located trial/complement floor used by a large-U local-band proof is not separated on this split-site source. Its local-band and nuclear graph enclosures do not transfer.'},
  'scope':'Complete fixed-graph native three-quark singlet spectra and original occupation records. Geometry changes the rest bands at fixed parameters. No continuum spatial scale, asymptotic cluster binding or full six-quark graph spectrum is inferred.'}
