"""Explicit record-only integer lift and dual-cell readout on frozen C10 source."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4];C10=ROOT/'openwave-tuning/next/C10'
sys.path[:0]=[str(C10),str(ROOT/'openwave-tuning/next/C11')]
from joint_source import JointNativeSource
from spectral_readout import GroupedEnergyInstrument
from boundary_charge import decode
EDGES=((0,3),(1,3),(1,4),(2,4),(2,5),(0,5))
G=np.array([1,1,1,0,0,0],dtype=np.int64)
D=np.zeros((6,6),dtype=np.int64)
for e,(a,b) in enumerate(EDGES):D[a,e]=1;D[b,e]=-1
TREE=tuple(range(5));ROOT_VERTEX=0

def lift_record(records,boundary):
 """Controller inputs are only exact Fmod6/B and declared G/Pauli domain."""
 if type(boundary) is not int or not -3<=boundary<=3:raise ValueError('integer boundary B in[-3,3] required')
 if len(records)!=6 or any(type(x) is not int or not 0<=x<6 for x in records):raise ValueError('six canonical Z6 half-records required')
 f=np.asarray(records,dtype=np.int64);number=3-boundary
 raw=G-D@f;raw[0]-=boundary
 # Singleton interval decoder; totalN and Pauli bounds are explicit premises.
 occupations=np.array([decode(int(x%6),number,1,6,6) for x in raw],dtype=np.int64)
 if sum(occupations)!=number:raise ValueError('record/background violates complete Gauss occupation sector')
 desired=G-occupations;desired[0]-=boundary
 delta=desired-D@f
 if np.any(delta%6) or sum(delta)!=0:raise ValueError('nonneutral integer-lift correction')
 residual=delta//6;correction=np.zeros(6,dtype=np.int64)
 # Tree0..4 is path0-3-1-4-2-5. Leaf elimination, no n oracle.
 order=(0,3,1,4,2,5);parents={3:(0,0),1:(3,1),4:(1,2),2:(4,3),5:(2,4)}
 for vertex in reversed(order[1:]):
  parent,edge=parents[vertex];value=int(residual[vertex]);correction[edge]=value*D[vertex,edge]
  residual[parent]+=value;residual[vertex]=0
 if residual[0]:raise ValueError('tree lift residual')
 lifted=f+6*correction;degree=D@lifted
 if not np.array_equal(degree,desired):raise ArithmeticError('integer lift failed')
 return {'lifted_records':lifted,'occupations':occupations,'degrees':degree,'tree_correction':correction}

def grouped_poststates(rho,labels):
 result=[]
 for label in sorted(set(labels)):
  ids=np.array([i for i,x in enumerate(labels) if x==label]);post=np.zeros_like(rho)
  post[np.ix_(ids,ids)]=rho[np.ix_(ids,ids)]
  result.append((label,ids,post))
 return result

def execute():
 model=JointNativeSource();source=model.source;N=model.n;boundaries=3-model.numbers
 lifted=[];occup=[];degrees=[]
 for (mask,f),boundary in zip(model.states,boundaries):
  data=lift_record(list(f),int(boundary));lifted.append(data['lifted_records']);occup.append(data['occupations']);degrees.append(data['degrees'])
 lifted=np.asarray(lifted);occup=np.asarray(occup);degrees=np.asarray(degrees)
 actual=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in model.states]);assert np.array_equal(occup,actual)
 # Dualsphere:2vertices,6meridianedges,6bigonfaces. No discarded multiedges.
 s=np.array([1,-1,1,-1,1,-1]);dual_boundary=np.vstack((-s,s));face_boundary=D.T
 assert np.array_equal(dual_boundary@face_boundary,np.zeros((2,6),int))
 assert np.linalg.matrix_rank(dual_boundary)==1 and np.linalg.matrix_rank(face_boundary)==5
 # Each edgecurve exp(2pi i Fhat*t) has integral degree Fhat and commonvertexvalue1.
 # Faceboundary degrees are matrixcomposition, never programmed from expectedn.
 assert np.array_equal(lifted@face_boundary, degrees)
 # Exact full-H/native support update covariance, including finite modular wrap.
 hcount=kcount=wraps=0
 for row,col in zip(*np.nonzero(model.hamiltonian)):
  if row==col:continue
  assert np.array_equal(degrees[row]-degrees[col],-(occup[row]-occup[col]))
  assert boundaries[row]==boundaries[col];hcount+=1
  wraps+=int(np.any(abs(np.array(model.states[row][1])-np.array(model.states[col][1]))>1))
 native_operator_count=0
 for law in ('C','E','MIX'):
  for ops in source.operators[law]:
   for meta,k in ops:
    native_operator_count+=1
    for row,col in zip(*np.nonzero(k)):
     assert boundaries[row]==boundaries[col]
     assert np.array_equal(degrees[row]-degrees[col],-(occup[row]-occup[col]));kcount+=1
 # Root reservoir is crucial rejecting control: spatial compensateddegreesunchanged.
 charged_transitions=0
 for meta,k in source.reservoir_instrument():
  for row,col in zip(*np.nonzero(k)):
   assert np.array_equal(degrees[row],degrees[col])
   if boundaries[row]!=boundaries[col]:charged_transitions+=1
 assert charged_transitions>0 and wraps>0
 # Explicit boundary-register geometry attachment: edge boundary6 -> root0,
 # flux-B. Its dualedge is a selfloop; no new dynamical register is introduced.
 Daug=np.zeros((7,7),dtype=np.int64);Daug[:6,:6]=D;Daug[6,6]=1;Daug[0,6]=-1
 Faug=np.column_stack((lifted,-boundaries));waug=Faug@Daug.T
 target=np.column_stack((G[None,:]-occup,-boundaries))
 assert np.array_equal(waug,target)
 dual_aug=np.column_stack((dual_boundary,np.zeros(2,dtype=np.int64)))
 assert not np.any(dual_aug@Daug.T) and np.linalg.matrix_rank(Daug)==6
 # Reservoir shiftsglobalcharge via existingB edge, now visible atboundaryface.
 for meta,k in source.reservoir_instrument():
  for row,col in zip(*np.nonzero(k)):
   assert waug[row,6]-waug[col,6]==model.numbers[row]-model.numbers[col]
 # All64 regional projectors match actualoccupation, with Bcontrol explicit.
 comparisons=0
 for bits in range(64):
  region=[v for v in range(6) if bits>>v&1];root=int(0 in region)
  decoded=G[region].sum()-degrees[:,region].sum(axis=1)-boundaries*root
  assert np.array_equal(decoded,occup[:,region].sum(axis=1));comparisons+=N
 # Actual complete charged-preparation source output, no selected species.
 archive=np.load(C10/'joint_source_states.npz');assert np.array_equal(archive['hamiltonian'],model.hamiltonian)
 rho=sum(archive[f'global_allN_reservoir_{i}_{j}'] for i in (0,1) for j in (0,1));assert abs(np.trace(rho)-1)<1e-12
 region=(0,1,3);winding=degrees[:,region].sum(axis=1);q=occup[:,region].sum(axis=1)
 labels=[(int(b),int(w)) for b,w in zip(boundaries,winding)];arrays={'lifted_records':lifted,'decoded_occupations':occup,'face_degrees':degrees,'boundary_B':boundaries,'augmented_lifted_records':Faug,'augmented_face_degrees':waug,'input_state':rho};rows=[]
 # Complete declaredjoint apparatus includeszero outcomes for all B/q pairs.
 for boundary in range(-3,4):
  for count in range(4):
   degree=int(sum(G[list(region)])-count-boundary);ids=np.flatnonzero((boundaries==boundary)&(winding==degree));post=np.zeros_like(rho,dtype=complex);post[np.ix_(ids,ids)]=rho[np.ix_(ids,ids)]
   key=f'joint_B{boundary}_Q{count}';arrays[key]=post
   rows.append({'B':boundary,'spatial_degree':degree,'regional_occupation':count,'probability':float(np.trace(post).real),'indices':ids.tolist(),'poststate_key':key})
 augmented_readout=[]
 for mask in range(64):
  pattern=np.array([(mask>>v)&1 for v in range(6)]);ids=np.flatnonzero(np.all(occup==pattern,axis=1));coefficient=rho[np.ix_(ids,ids)];key=f'augmented_whole_coefficient_{mask}';arrays[key]=coefficient
  signature=np.r_[G-pattern,mask.bit_count()-3]
  assert np.all(waug[ids]==signature)
  augmented_readout.append({'face_degrees':signature.tolist(),'probability':float(np.trace(coefficient).real),'indices':ids.tolist(),'coefficient_key':key,'reconstruction':'Embed(indices,coefficient); completecomplexwithin-Ffiber matrix retained'})
 assert abs(sum(x['probability'] for x in augmented_readout)-1)<1e-12
 coarse=[]
 for count in range(4):
  ids=np.flatnonzero(q==count);post=np.zeros_like(rho,dtype=complex);post[np.ix_(ids,ids)]=rho[np.ix_(ids,ids)];key=f'coarse_Q{count}';arrays[key]=post
  coarse.append({'regional_occupation':count,'probability':float(np.trace(post).real),'indices':ids.tolist(),'poststate_key':key})
 assert abs(sum(r['probability'] for r in rows)-1)<1e-12 and abs(sum(r['probability'] for r in coarse)-1)<1e-12
 selected=max(rows,key=lambda r:r['probability']);post=arrays[selected['poststate_key']]
 after,ledger=GroupedEnergyInstrument(model).continue_branch(post);arrays['continued_joint_poststate']=after
 assert abs(sum(r['joint_probability'] for r in ledger)-selected['probability'])<1e-12
 # FineFforgetting destroyswithinclasscoherence, unlike declared groupedPVM.
 chosen=next(ids for _,ids,_ in grouped_poststates(rho,labels) if len(ids)>1);a,b=map(int,chosen[:2]);control=np.zeros(N,complex);control[a]=1/np.sqrt(2);control[b]=1j/np.sqrt(2)
 arrays['fine_forget_control_state']=np.outer(control,control.conj())
 fine_distance=.5*float(np.linalg.norm(np.outer(control,control.conj())-np.diag(abs(control)**2),ord='nuc'));assert abs(fine_distance-.5)<1e-12
 # SameoccupationdifferentB control showsjoint(B,w) isnot coarseQ alone.
 pair=next((int(i),int(j)) for i in range(N) for j in range(i+1,N) if q[i]==q[j] and boundaries[i]!=boundaries[j] and (model.numbers[i]-model.numbers[j])%2==0)
 c=np.zeros(N,complex);c[pair[0]]=1/np.sqrt(2);c[pair[1]]=1j/np.sqrt(2);cross=np.outer(c,c.conj());joint=sum(x[2] for x in grouped_poststates(cross,labels));coarse_control=sum(x[2] for x in grouped_poststates(cross,list(map(int,q))));arrays['crossB_coarse_control_poststate']=coarse_control;arrays['crossB_joint_control_poststate']=joint
 b_distance=float(np.linalg.norm(joint-coarse_control,ord='nuc')/2);assert abs(b_distance-.5)<1e-12
 paths=[Path(__file__),C10/'joint_source.py',C10/'charged_reservoir.py',C10/'modular_field_adapter.py',C10/'spectral_readout.py',ROOT/'openwave-tuning/next/C11/boundary_charge.py',C10/'joint_source_states.npz']
 np.savez_compressed(HERE/'quantum_spatial_states.npz',**arrays)
 return {'schema':'oph.quantum.spatial_record_bridge.v1','checks':{'all384_record_only_decoders_equal_actualPaulioccupation':True,'exact_integer_tree_lift_and_neutrality':True,'dual2vertex6edge6bigon_boundary_squared_zero':True,'allH_and_allnative_support_updates_covariant_includingwrap':True,'all64_regional_readout_projectors_equal_actualoccupation_withB':True,'complete_joint_and_coarse_Luders_poststates_retained':True,'full_sameH_native_continuation_complete':True,'fine_record_and_extraB_dephasing_controls_detected':True,'charged_reservoir_changes_globalQ_without_spatialdegree':True,'explicit_existingB_boundary_attachment_exposes_globalQ':True},'parameters':{'source':'frozen C10 allnumber384/Z6 interactingCAR','fixed_G':G.tolist(),'boundary_root':0,'primal_halfedges':EDGES,'tree_edges':TREE,'region':region,'modulus':6},'topology':{'dual_boundary1':dual_boundary.tolist(),'dual_boundary2':face_boundary.tolist(),'cell_counts':[2,6,6],'Euler_characteristic':2,'face_degree_definition':'w=D Fhat, actual lifted edgecurve degrees; w=G-n-Bdelta_root','global_charge_definition':'Q=N-3=-B; sum spatialface degrees0. Not universalQ=f(spatialw).','apparatus_attachment':'Oriented sphere cellulation dual to actual6cycle,6bigonfaces/multiedges retained. Fixedtree integer representative chosen explicitly; no newsourceedge/dynamics.'},'boundary_attachment':{'primal_edge':[6,0],'actual_register_flux':'-B, existing integerboundaryregister; no newdynamicalfield','dual_boundary1':dual_aug.tolist(),'dual_boundary2':Daug.T.tolist(),'dual_cell_counts':[2,7,7],'source_face_degrees':'G-n','boundary_face_degree':'-B=N-3=globalQ','apparatus_status':'Explicit appended boundarypuncture geometry for existingroot/reservoirservice; not intrinsic old6cyclegeometry or uniqueemergentphysicalspace.'},'source_checks':{'H_support_transitions':hcount,'native_support_transitions':kcount,'native_operator_count':native_operator_count,'wrapping_H_transitions':wraps,'all_region_basis_comparisons':comparisons,'charged_reservoir_transitions_with_unchanged_spatial_degree':charged_transitions},'complete_augmented_whole_degree_readout':augmented_readout,'complete_joint_B_winding_readout':rows,'complete_coarse_regional_occupation_readout':coarse,'selected_joint_continuation':{'branch':selected,'complete_outcomes':ledger,'poststate_key':'continued_joint_poststate','same_declared_coherent_time':.37,'commit_probability':.29,'law':'MIX'},'mutations':{'fine_record_forgetting_trace_distance':fine_distance,'fine_control_basis_indices':[a,b],'jointB_vs_coarseQ_crossnumber_distance':b_distance,'crossnumber_control_basis_indices':pair,'noB_input_and_mod6_arbitrary_lift_not_admitted':True},'scope':'Record-only topology/readout adapter on same quantumsource, not an intrinsic electromagnetic field derivation. FixedGauss, Pauli interval, integerboundaryB, tree lift and dualcell apparatus explicit. GroupedPVMs are added coherent apparatus choices, not claimed nativephysicalcompiler. Joint(B,w) preserveswithinclasscoherences but dephasesB; coarseQ may be computed coherently without exposingB. Spatialcharge isboundarycompensated; reservoir globalchargecounterexample retained.','state_archive':'quantum_spatial_states.npz','state_archive_sha256':hashlib.sha256((HERE/'quantum_spatial_states.npz').read_bytes()).hexdigest(),'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'official_rows_changed':[]}
if __name__=='__main__':
 r=execute();(HERE/'QUANTUM_SPATIAL.json').write_text(json.dumps(r,indent=2)+'\n');print(r['checks']);print(r['source_checks'])
