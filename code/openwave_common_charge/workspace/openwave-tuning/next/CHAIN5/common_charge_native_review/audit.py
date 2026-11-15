"""Independent complete native review and a separately authored terminal annex.
No common_charge_topology execution, readout or topology routine is imported.
The existing primitive compiler and this reviewer's frozen raw-M6 executor are
reused as declared dependencies; all producer numbers are recomputed here.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import defaultdict,Counter
import hashlib,importlib.util,json,math,random,sys
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
WORK=HERE.parents[3]
SHARED=WORK/'openwave-tuning/next/CHAIN5/dirac_kg_opposing_review'
P='openwave-tuning/next/CHAIN5/common_charge_topology/'
FIELD='openwave-tuning/next/CHAIN3/native_field_connection/'

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,d):
 (HERE/name).write_text(json.dumps(d,indent=2,default=str)+'\n');print(name,'PASS',flush=True)
def frozen(name):
 p=HERE/'inputs'/name
 return p if p.exists() else SHARED/'inputs'/name

def custody():
 d=read(HERE/'INPUTS.json')
 assert sha(SHARED/'MANIFEST.json')==d['shared_review_manifest']
 assert sha(SHARED/'audit.py')==d['shared_executor_sha256']
 for group in ('copied_files','shared_frozen_files'):
  for name,digest in d[group].items():assert sha(frozen(name))==digest,name
 return {k:len(d[k]) for k in ('copied_files','shared_frozen_files')}

spec=importlib.util.spec_from_file_location('earlier_independent_native_executor',SHARED/'audit.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)

class Executor(prior.IndependentNative):
 """Fresh complete per-column first-export replay, independent of producer."""
 def full(self,tape,bits):
  dim=1<<bits;M=np.eye(dim,dtype=complex);failed=np.zeros((dim,dim),complex);rows=[];bound=F()
  for j,event in enumerate(tape):
   signature=json.dumps(event,sort_keys=True)
   if signature not in self.certified:self.compiler.native.verify_pulse(event);self.certified.add(signature)
   U=self.response(tuple(event['controls_exact']));ix=self.groups(bits,tuple(event['targets']));old=M[ix];d=ix.shape[1]
   lost=np.einsum('ij,bjk->bik',U[d:,:d],old).reshape(-1,dim);effect=lost.conj().T@lost
   failed+=effect;M[ix]=np.einsum('ij,bjk->bik',U[:d,:d],old)
   rows.append({'pulse':j,'targets':event['targets'],'basis_column_probabilities':np.diag(effect).real.tolist(),'trace':float(np.trace(effect).real)})
   bound+=F(event['full_response_operator_error_exact']);self.pulses+=1
  gram=float(np.linalg.norm(M.conj().T@M+failed-np.eye(dim),2));assert gram<1e-10
  return M,failed,bound,gram,rows

def tensor(effects):
 a=np.ones((1,1),complex)
 for E in effects:a=np.kron(E,a)
 return a

def reduced(out,cellbits=3):
 dim=out.shape[0]//(1<<cellbits);rho=np.zeros((dim,dim),complex)
 for cell in range(1<<cellbits):
  block=out[cell::(1<<cellbits)];rho+=block@block.conj().T
 return rho

def pointer(native,tape):
 # The generic earlier review reconstructs endpoint source C/E effects from
 # the actual native event formula, not from producer readout JSON.
 effects,summary=prior.pointer_effects(native,tape)
 _,failed,_,_,_=native.run(tape,3,np.eye(8));failure=failed[np.ix_([0,4],[0,4])]
 for law in effects:
  assert np.linalg.norm(sum((E for _,E in effects[law]),np.zeros((2,2),complex))+failure-np.eye(2),2)<1e-12
 return effects,failure,summary

def complete_law(rho,effects,pfailed,laws,full=False):
 choices=[effects[law] for law in laws];coarse=np.zeros(1<<len(laws));fine=[]
 for rows in product(*choices):
  code=sum(meta['bit']<<j for j,(meta,E) in enumerate(rows));prob=float(np.trace(tensor([E for _,E in rows])@rho).real)
  assert prob>-2e-13;coarse[code]+=prob
  if full:fine.append({'records':[meta for meta,E in rows],'probability':prob})
 failures=[]
 for j in range(len(laws)):
  E=tensor([sum((e for _,e in choices[k]),np.zeros((2,2),complex)) if k<j else pfailed if k==j else np.eye(2) for k in range(len(laws))])
  failures.append(float(np.trace(E@rho).real))
 gram=abs(float(coarse.sum()+sum(failures)-np.trace(rho).real));assert gram<1e-10
 return {'coarse':coarse.tolist(),'fine':fine,'first_pointer_exports':failures,'completeness_residual':gram}

def geometry(c):
 # Independently rebuild incidence from actual graph seams; use source face
 # list solely to specify the native complex, never bridge.exact_geometry.
 graph,producerB,A,periods,winding,labels=c.g.source.field.cell_complex.build(2)
 edges=list(graph.seams);B=s.zeros(96,288)
 for j,(a,b) in enumerate(edges):B[a,j]=1;B[b,j]=-1
 assert B==producerB and B*A==s.zeros(96,208)
 data=read(frozen(FIELD+'COORDINATES.json'))
 Ps=s.Matrix.hstack(*(s.Matrix([s.Rational(x) for x in p['particular']]) for p in data['charge_sectors']))
 H=s.Matrix.hstack(*(s.Matrix([s.Rational(x) for x in row]) for row in data['harmonic_orthogonal_rational_vectors']))
 Q=s.Matrix(c.g.data()['Q']);K=A*A.T
 charges=s.zeros(96,8)
 for i in range(8):charges[0,i]+=1;charges[12*i,i]-=1
 assert B*Ps==charges and B*H==s.zeros(96,3) and A.T*H==s.zeros(208,3)
 assert B*K==s.zeros(96,288)
 # Each exact source route is reconstructed from its literal path.
 routes=[];perms={axis:[] for axis in range(3)}
 for row in data['currents']:
  J=s.zeros(288,1)
  for a,b in zip(row['path'],row['path'][1:]):
   j=edges.index(tuple(sorted((a,b))));J[j]+=1 if a<b else -1
  assert list(J)==row['integer_current']
  old,new=row['cell'],row['destination'];D=Ps[:,new]-Ps[:,old]
  assert B*J==charges[:,new]-charges[:,old] and B*(J-D)==s.zeros(96,1)
  assert A.T*D==s.zeros(208,1) and H.T*D==s.zeros(3,1)
  assert (D.T*(J-D))[0]==0
  assert new==old^(1<<(2-row['axis']))
  routes.append({'cell':old,'axis':row['axis'],'sign':row['sign'],'destination':new,'native_path':row['path'],'exact_Gauss_and_orthogonal_current_join':True})
 regions=read(frozen(P+'TOPOLOGY.json'))['regions'];fresh={'alternating_roots':[0,24,48,72],'nonroot_ball':[37],'mixed_roots_ports':[12,25,48,71,84]}
 cuts=[]
 for name,vs in {**{k:v['vertices'] for k,v in regions.items()},**fresh}.items():
  indicator=s.Matrix([int(v in vs) for v in range(96)]);cut=B.T*indicator;values=[int(0 in vs)-int(12*c in vs) for c in range(8)]
  assert list((cut.T*Ps))==values and cut.T*A==s.zeros(1,208) and cut.T*H==s.zeros(1,3)
  if name in regions:assert values==regions[name]['flux_eigenvalues_by_cell']
  cuts.append({'name':name,'vertices':vs,'flux_eigenvalues':values,'boundary_edges':[[j,int(x)] for j,x in enumerate(cut) if x]})
 # Two-disjoint daughters are orthogonal projectors, not fractions of one
 # eigenvalue. Equal expectations are explicitly an allowed control.
 half=(Ps[:,1]+Ps[:,2])/2;halfdiv=B*half;assert halfdiv[0]==1 and halfdiv[12]==halfdiv[24]==-s.Rational(1,2)
 BQ=B*Q;qv,qj=next((v,j) for v in range(96) for j in range(190) if BQ[v,j]);amp=1/(2*BQ[qv,qj]);assert (BQ[qv,qj]*amp)==s.Rational(1,2)
 topo=read(frozen(P+'TOPOLOGY.json'));sphere=topo['sphere_complex'];D=s.Matrix(sphere['boundary']);ones=s.ones(len(sphere['oriented_faces']),1)
 assert D*ones==s.zeros(len(sphere['edges']),1) and D.rank()==11
 assert math.gcd(*sphere['primitive_fundamental_cycle'])==1
 # Fresh exact field sample, all96 primitive surface periods in each sector.
 transverse=A[:,3]/19+A[:,111]/23-H[:,0]/29+H[:,2]/31;surfacechecks=0
 for cell in range(8):
  E=Ps[:,cell]+transverse
  for v in range(96):
   flux=sum(B[v,j]*E[j] for j in range(288));assert flux==charges[v,cell];surfacechecks+=1
 globalhalf=H*(H.T*H).inv()*H.T*winding[:,0]/2
 assert B*globalhalf==s.zeros(96,1) and A.T*globalhalf==s.zeros(208,1) and periods*globalhalf==s.Matrix([s.Rational(1,2),0,0])
 return {'routes':routes,'fresh_and_archived_cuts':cuts,'fresh_primitive_surface_checks':surfacechecks,'sphere_boundary_rank':11,
  'exact_polar_frame_transverse_reason':'B K=0; positive polar normalization of K Q Lambda^-1 preserves its rank190 range im(A). Rounded Q is not the exact field.',
  'half_daughter_control':{'background':1,'daughter_fluxes':['-1/2','-1/2'],'transition_endpoints':[-1,-1],'total_relative_charge':1},
  'rounded_Q_control':{'vertex':qv,'mode':qj,'BQ_entry':str(BQ[qv,qj]),'amplitude':str(amp)},
  'global_harmonic_control_periods':['1/2','0','0'],'one_macro_destination_by_initial_cell':[i^7 for i in range(8)],
  'no_retained_current_history':True,'full_source_g':1},regions

def reset_checks(c):
 # Source controller is reused as definition. Every emitted source update
 # and C/E label is independently replayed using the integer current formula.
 from model.native_clock import configuration as config
 states=tuple(config.native.states());seams=tuple(config.quantum.seams());target=(0,0,1,1,0,0,0,0,0,0,0,0)
 checks=0
 for seam in seams:
  for law in ('C','E'):
   groups=defaultdict(dict);coverage=Counter()
   for state in states:
    for coin in (0,1):
     a,b=seam;d=state[a]-state[b];j=(d-(2*coin-1)*(abs(d)%2))//2
     out=list(state);out[a]-=j;out[b]+=j;out=tuple(out);label=j if law=='C' else d
     groups[coin,label][state]=out;coverage[state]+=1;checks+=1
   assert set(coverage.values())=={2} and len(coverage)==78
   assert all(len(g)==len(set(g.values())) for g in groups.values())
 successes=[];failures=[];totalrecords=0
 for law in ('C','E'):
  for k,state in enumerate(states):
   rng=random.Random(17011+991*k+(law=='E'));coins=(rng.randrange(2) for _ in range(6000))
   execution=config.execute(state,coins,full=True,mode='reset',target=target,law=law)
   current=state
   for seam,coin,label in execution.records:
    a,b=seam;d=current[a]-current[b];j=(d-(2*coin-1)*(abs(d)%2))//2
    assert label==(j if law=='C' else d)
    out=list(current);out[a]-=j;out[b]+=j;current=tuple(out)
   assert current==execution.output and execution.outcome.success and current==target
   support=config.word_support(execution,full=True);assert support=={state:target}
   successes.append({'law':law,'input':state,'attempts':len(execution.records),'records':execution.records,'support_columns':1})
   totalrecords+=len(execution.records)
  execution=config.execute(states[0],iter([0]*6000),full=True,mode='reset',target=target,law=law)
  assert not execution.outcome.success
  failures.append({'law':law,'input':states[0],'reason':execution.outcome.reason,'records':execution.records,'output':execution.output})
 receipt=read(c.native.PROVIDER/'COMMON_WALK_RECEIPT.json')['preparation_reset']['certificate']
 assert receipt['input_columns']==78 and receipt['failure_upper']=='35/2147483648'
 probabilities=config.success_probabilities(full=True,mode='reset',target=target)
 assert len(probabilities)==78 and max(1-x for x in probabilities)<=F(receipt['failure_upper'])
 save('FRESH_RESET_RECORDS.json',{'successful_tapes':successes,'timeout_controls':failures,'exact_success_probabilities_all78':list(map(str,probabilities))})
 return {'independent_native_source_column_coin_checks':checks,'fresh_reset_successes':len(successes),'fresh_reset_record_count':totalrecords,
  'explicit_timeout_controls':len(failures),'base_uniform_failure_exact':receipt['failure_upper'],'exact_all78_success_probability_recomputed':True,
  'arbitrary_reference_scope':'Conclusive native tape has one input-column support. Discarding its private resolved input/tape, success is blank times subnormalized original reference; full first-timeout instrument is retained.'}

def replay(c,native,regions):
 programs=read(frozen(P+'NATIVE_PROGRAMS.json'));archfail=read(frozen(P+'NATIVE_FIRST_FAILURES.json'));records=read(frozen(P+'NATIVE_RECORDS.json'))
 effects,pfailed,psummary=pointer(native,programs['native_READY_pointer']);replayed={};summaries=[];maxfaildiff=0.
 for name,tape in programs.items():
  bits=3 if name=='native_READY_pointer' else 6 if '+' in name else 5
  U,failed,bound,gram,stages=native.full(tape,bits);arch=archfail[name]
  expected=np.array(arch['failure_effect_real'])+1j*np.array(arch['failure_effect_imag'])
  assert np.linalg.norm(failed-expected,2)<1e-12 and len(stages)==len(arch['first_export_records'])
  for row,old in zip(stages,arch['first_export_records']):
   error=max(abs(np.array(row['basis_column_probabilities'])-old['probability_on_every_input_basis_column']));maxfaildiff=max(maxfaildiff,float(error));assert error<1e-13
  replayed[name]=(U,failed,bound,gram);summaries.append({'program':name,'pulses':len(tape),'all_input_columns':1<<bits,'complete_Gram_residual':gram,'error_bound_exact':str(bound)})
 for law,rows in records['fine_pointer_effects'].items():
  for old,(meta,E) in zip(rows,effects[law]):
   assert old['record']==meta and np.linalg.norm(E-np.array(old['effect_real'])-1j*np.array(old['effect_imag']),2)<1e-12
 plus_tape=c.events(c.words.remap(c.native.primitive_templates()['hadamard'],(0,)))
 plus,_,plusbound,_,_=native.run(plus_tape,1,np.eye(2));plus=plus[:,0]
 prefixes={x['macros']:x for x in read(frozen(FIELD+'RESULT.json'))['source_contracts']}
 improved=read(frozen('openwave-tuning/next/CHAIN4/native_field_opposing/BOUNDS.json'))
 newprefix={x['macros']:F(x['new_complete_record_TV_exact']) for x in improved['rows']};budgetrows=[];finecount=0;arrayerror=0.
 for row in records['experiments']:
  names=row['regions'];r=len(names);n=row['macros'];clock=3+r;dim=1<<(4+r);name='+'.join(names)+'_'+row['clock_setting'];U,failed,bound,gram=replayed[name]
  ticks=c.events([('Z',(0,),c.native.angle(0,0,F(-1,128)))]*n)
  T,_,tb,_,_=native.run(ticks,1,np.eye(2));cv=T@plus
  initial=np.zeros((dim,8),complex)
  for cell in range(8):initial[cell,cell]=cv[0]/math.sqrt(8);initial[cell+(1<<clock),cell]=cv[1]/math.sqrt(8)
  rho=reduced(U@initial)
  for laws,arch in row['laws'].items():
   out=complete_law(rho,effects,pfailed,laws,full=True);assert len(out['fine'])==len(arch['fine_records'])
   for fresh,old in zip(out['fine'],arch['fine_records']):
    assert fresh['records']==old['records'];diff=abs(fresh['probability']-old['probability']);arrayerror=max(arrayerror,diff);assert diff<1e-11;finecount+=1
   assert max(abs(np.array(out['coarse'])-arch['all_joint_binary_probabilities']))<1e-11
   assert max(abs(np.array(out['first_pointer_exports'])-arch['first_pointer_export_probability_by_pointer']))<1e-12
  reset=F(3*r+2)*F(35,2**31)**2;meter=bound+(r+1)*F(psummary['word_bound_exact'])+reset
  assert meter==F(row['additional_meter_bound_exact'])
  old=F(prefixes[n]['complete_source_record_TV_bound_exact']);assert old+meter==F(row['complete_common_source_record_TV_bound_exact'])
  # Explicitly price the displayed 100-bit diagnostic source-clock surrogate.
  # It is not silently identical to the 237/258-bit inherited source clock.
  surrogate=plusbound+tb;safe=newprefix[n]+meter+surrogate;assert safe<old+meter
  budgetrows.append({'regions':names,'setting':row['clock_setting'],'macros':n,'meter_bound_exact':str(meter),'diagnostic_clock_surrogate_extra_exact':str(surrogate),
    'improved_plus_surrogate_complete_exact':str(safe),'published_complete_exact':str(old+meter),'surrogate_absorbed_by_independently_reviewed_Gaussian_slack':True})
 # Fresh rejected source controls, including a non-half fractional mutant.
 controls={};controltapes={}
 ref=regions['reference_singleton']['flux_eigenvalues_by_cell']
 wrong=c.events(annex_meter(c,[[x-1 for x in ref]],'X'));wrongU,_,wrongbound,wronggram,_=native.full(wrong,5)
 separation=float(np.linalg.norm(wrongU-replayed['reference_singleton_X'][0],2));assert separation>1.9
 controls['wrong_reference']={'all_columns':32,'operator_separation':separation,'word_error_exact':str(wrongbound),'Gram_residual':wronggram};controltapes['wrong_reference']=wrong
 particle=regions['particle_singleton']['flux_eigenvalues_by_cell'];third=c.events(annex_meter(c,[[F(2,3)*x for x in particle]],'X'))
 thirdU,_,thirdbound,thirdgram,_=native.full(third,5);thirdprob=float(sum(abs(thirdU[x,1])**2 for x in range(32) if x&8));assert abs(thirdprob-.75)<1e-11
 controls['one_third_charge_2pi_mutant']={'all_columns':32,'wrong_bit_probability':thirdprob,'ideal_wrong_bit_probability':'3/4','word_error_exact':str(thirdbound),'Gram_residual':thirdgram};controltapes['one_third_charge_2pi_mutant']=third
 integer=float(sum(abs(replayed['integer_2pi_loop'][0][x,1])**2 for x in range(32) if x&8));half=float(sum(abs(replayed['half_charge_2pi_mutant'][0][x,1])**2 for x in range(32) if x&8))
 assert integer<1e-20 and half>1-1e-12
 controls['archived_integer_and_half']={'integer_wrong_bit_probability':integer,'half_wrong_bit_probability':half}
 save('FRESH_CONTROL_PROGRAMS.json',controltapes)
 return {'fresh_and_replayed_rejection_controls':controls,'program_replays':summaries,'archived_experiments':len(records['experiments']),'fine_probabilities_recomputed':finecount,
  'fine_probability_max_difference':arrayerror,'per_pulse_all_column_failure_max_difference':maxfaildiff,'pointer':psummary,
  'complete_bound_rows':budgetrows,'full_field_native_prefix_executed':False},effects,pfailed,psummary,programs

# Named new producer annex starts here. This is not an independent review of
# itself. It must receive a separate review before any acceptance uses it.
def annex_meter(c,values,clocksetting):
 r=len(values);clock=3+r;word=[]
 for j,vals in enumerate(values):
  probe=3+j;word+=list(c.words.remap(c.native.primitive_templates()['hadamard'],(probe,)))
  # Independent full padded 16-state diagonal Fourier transform. Exp(i*pi
  # Phi) only on probe1; both dirty probe columns are included in the word.
  diag=[F(vals[x&7])*((x>>3)&1) for x in range(16)]
  for mask in range(16):
   coeff=sum((v*(-1)**((x&mask).bit_count()) for x,v in enumerate(diag)),F())/16
   label=''.join('Z' if mask>>b&1 else 'I' for b in range(4))
   word+=list(c.words.remap(c.native.pauli(label,c.native.angle(coeff)),(0,1,2,probe)))
  word+=list(c.words.remap(c.native.primitive_templates()['hadamard'],(probe,)))
 if clocksetting=='Y':word.append(('Z',(clock,),c.native.angle(F(1,4))))
 word+=list(c.words.remap(c.native.primitive_templates()['hadamard'],(clock,)))
 return word

def ideal_terminal_operator(native,bits,values,setting,n,prepname):
 """Independent direct target, including all dirty-probe/clock columns."""
 dim=1<<bits;M=np.eye(dim,dtype=complex);clock=bits-1
 X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
 H=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
 def one(G,b):
  ix=native.groups(bits,(b,));M[ix]=np.einsum('ij,bjk->bik',G,M[ix])
 if prepname=='uniform_plane':
  for b in range(3):one(H,b)
 elif prepname=='basis_cell2':one(X,1)
 else:
  one(math.cos(math.pi/7)*np.eye(2)+1j*math.sin(math.pi/7)*Y,0);one(H,2)
 one(H,clock)
 for macro in range(n):
  for b in range(3):one(X,b)
  one(np.diag(np.exp([-1j/128,1j/128])),clock)
 for j,vals in enumerate(values):
  one(H,3+j)
  # Source Fourier observable on each literal occupied-cell code.
  for code in range(dim):
   if code>>(3+j)&1:M[code]*=np.exp(1j*math.pi*float(vals[code&7]))
  one(H,3+j)
 if setting=='Y':one(np.diag(np.exp([1j*math.pi/4,-1j*math.pi/4])),clock)
 one(H,clock)
 return M

def annex(c,native,regions,effects,pfailed,psummary):
 programs={};failures={};rows=[];maxlaw=0.;maxoperator=0.;worst=F();primitive=c.native.primitive_templates()
 selections=[('fresh_alternating_and_daughter',[[0,24,48,72],[12,36,60,84]]),
             ('original_daughters',[regions['daughter_A']['vertices'],regions['daughter_B']['vertices']]),
             ('single_particle',[[12]])]
 preparations=[('uniform_plane',list(c.words.remap(primitive['hadamard'],(0,)))+list(c.words.remap(primitive['hadamard'],(1,)))+list(c.words.remap(primitive['hadamard'],(2,))),[F(1,8)]*8),
   ('basis_cell2',list(c.words.remap(primitive['flip'],(1,))),[F(int(i==2)) for i in range(8)]),
   ('nonuniform_entanglement_sensitive', [('Y',(0,),c.native.angle(F(1,7)))]+list(c.words.remap(primitive['hadamard'],(2,))),None)]
 for name,vertices in selections:
  values=[[int(0 in vs)-int(12*cell in vs) for cell in range(8)] for vs in vertices];r=len(values);clock=3+r;bits=4+r;dim=1<<bits
  for setting in ('X','Y'):
   meter=annex_meter(c,values,setting)
   for n in (1,2):
    prefix=[]
    for macro in range(n):
     for b in (0,1,2):prefix+=list(c.words.remap(primitive['flip'],(b,)))
     prefix.append(('Z',(clock,),c.native.angle(0,0,F(-1,128))))
    for prepname,prep,weights in preparations:
     key=f'{name}_{setting}_n{n}_{prepname}'
     prepword=prep+list(c.words.remap(primitive['hadamard'],(clock,)))
     tape=c.events(prepword+prefix+meter);U,failed,bound,gram,stages=native.full(tape,bits)
     ideal=ideal_terminal_operator(native,bits,values,setting,n,prepname);operator_error=float(np.linalg.norm(U-ideal,2));assert operator_error<1e-10;maxoperator=max(maxoperator,operator_error)
     programs[key]=tape;failures[key]={'first_export_records':stages,'failure_effect_real':failed.real.tolist(),'failure_effect_imag':failed.imag.tolist(),'complete_Gram_residual':gram}
     rho=reduced(U[:,[0]]);source_failure=float(failed[0,0].real)
     # Independent exact ideal preparation probabilities.
     if weights is None:
      weights=[0.]*8
      for b0,b2 in product((0,1),repeat=2):weights[b0+4*b2]=(math.sin(math.pi/7)**2 if b0 else math.cos(math.pi/7)**2)/2
     expected=np.zeros(1<<(r+1));mean=(math.cos if setting=='X' else math.sin)(n/64)
     for initial,w in enumerate(weights):
      cell=initial^(7 if n%2 else 0);code=sum(abs(v[cell])<<j for j,v in enumerate(values))
      expected[code]+=float(w)*(1+mean)/2;expected[code+(1<<r)]+=float(w)*(1-mean)/2
     lawsout={}
     for laws in product(('C','E'),repeat=r+1):
      result=complete_law(rho,effects,pfailed,laws,full=True);diff=max(abs(np.array(result['coarse'])-expected));maxlaw=max(maxlaw,float(diff));assert diff<3e-10
      assert abs(sum(result['coarse'])+sum(result['first_pointer_exports'])+source_failure-1)<1e-10
      lawsout[''.join(laws)]=result
     registers=6+3*r;timeout=F(registers)*F(35,2**31)**2
     total=bound+(r+1)*F(psummary['word_bound_exact'])+timeout;worst=max(worst,total)
     rows.append({'name':key,'regions':vertices,'flux_values_by_final_cell':values,'macros':n,'physical_time_exact':str(F(n,8)),
       'clock_phase_exact':str(F(n,64)),'preparation':prepname,'ideal_initial_cell_weights':list(map(float,weights)),
       'complete_prep_prefix_meter_pulses':len(tape),'prefix_pulses':7*n,'all_column_ideal_operator_error':operator_error,'all_input_columns':dim,'complete_Gram_residual':gram,
       'first_prep_prefix_meter_export_probability':source_failure,'laws':lawsout,'reset_registers':registers,'reset_blocks_per_register':2,
       'reset_timeout_union_bound_exact':str(timeout),'all_native_pulse_error_bound_exact':str(bound+(r+1)*F(psummary['word_bound_exact'])),
       'complete_terminal_record_TV_bound_exact':str(total),'complete_terminal_record_TV_bound':float(total),
       'reset_success_weights_not_renormalized':True,'original_huge_field_native_tape_executed':False})
 save('TERMINAL_ANNEX_PROGRAMS.json',programs);save('TERMINAL_ANNEX_FAILURES.json',failures);save('TERMINAL_ANNEX_RECORDS.json',rows)
 # A full-rank cell+field+reference state tests the exact discarded-system
 # channel statement with generic noncommuting controlled internal rotations.
 rng=np.random.default_rng(73917);internal=4;reference=3
 state=rng.normal(size=(8,internal,reference))+1j*rng.normal(size=(8,internal,reference));state/=np.linalg.norm(state)
 original=state.copy()
 for axis in range(3):
  out=np.zeros_like(state)
  for cell in range(8):
   A=rng.normal(size=(internal,internal))+1j*rng.normal(size=(internal,internal));V,_=np.linalg.qr(A)
   out[cell^(1<<(2-axis))]=V@state[cell]
  state=out
 error=0.
 for cell in range(8):
  # Trace internal system, retaining quantum reference, not merely scalars.
  actual=state[cell].T@state[cell].conj();target=original[cell^7].T@original[cell^7].conj();error=max(error,float(np.linalg.norm(actual-target,2)))
 assert error<1e-13
 result={'schema':'oph.terminal_charge_heisenberg_annex.v1','producer_role':'New annex authored by this reviewer; separate root audit required.',
  'complete_native_terminal_experiments':len(rows),'native_pulse_count_range':[min(len(p) for p in programs.values()),max(len(p) for p in programs.values())],
  'all_input_columns_including_dirty_probes_and_clock':True,'max_all_column_ideal_operator_error':maxoperator,'max_coarse_record_error':maxlaw,'max_complete_terminal_TV_exact':str(worst),'max_complete_terminal_TV':float(worst),
  'numerical_roundoff_scope':'Binary64 diagnostic discrepancies are not the exact-rational pulse error bounds. Exact certificates, not the floating residual, certify the stated TV bound.',
  'generic_entangled_reference_channel_diagnostic_error':error,'original_source_g':1,'source_manifest':'719036fe46769e68202731f765045dc3d75ebc8fe6e901e85319ce435776779f',
  'exact_ideal_terminal_effect_identity':'U_g1^dagger P_c U_g1 = P_(c xor7) on L2; after n, P_(c xor(7*(n mod2))). Source clock is unchanged h/8 per macro.',
  'observable_domain':'Terminal regional charge/electric flux and source clock; discard cell, spin, coexact and harmonic fields; arbitrary retained quantum reference. No intermediate measurement.',
  'not_claimed':['Equality of source unitaries','Equality of full native failure tapes','Any L>2 reduction','Any interacting dispersion/mass pole','Execution of huge original field prefix'],
  'complete_reset_semantics':'All first reset timeouts and fine source words are retained. Displayed pulse/pointer components precede subnormalized reset-success weights; budget uses conditional uniform first-failure union bound.',
  'ideal_reset_comparison':'Private reset fine records are marginalized for ideal blank-state comparison, and retained in reproducibility records.'}
 save('TERMINAL_ANNEX_RESULT.json',result);return result

def main():
 pins=custody();c=prior.load_compiler();native=Executor(c)
 geo,regions=geometry(c);save('EXACT_GEOMETRY_REVIEW.json',geo)
 resets=reset_checks(c);save('RESET_REVIEW.json',resets)
 replayed,effects,pfailed,psummary,programs=replay(c,native,regions);save('NATIVE_REPLAY_REVIEW.json',replayed)
 extra=annex(c,native,regions,effects,pfailed,psummary)
 assert custody()==pins
 save('RESULT.json',{'all_checks_pass':True,'custody':pins,'producer_manifest':'b45e64648b2379cdc0629a8789f74ca9a798a211e58e078d5c9f573e1513eb4f',
  'raw_native_pulses_executed':native.pulses,'distinct_exact_pulse_certificates_checked':len(native.certified),
  'signed_actual_current_routes':len(geo['routes']),'fresh_local_surface_checks':geo['fresh_primitive_surface_checks'],
  'archived_native_experiments':replayed['archived_experiments'],'archived_fine_probabilities':replayed['fine_probabilities_recomputed'],
  'new_terminal_annex_experiments':extra['complete_native_terminal_experiments'],'new_terminal_annex_max_TV':extra['max_complete_terminal_TV'],
  'new_terminal_annex_requires_independent_review':True,'huge_original_full_field_native_prefix_executed':False,
  'upstream_disclosure':'Reviewer authored Coulomb consumer27, Gaussian bound/task32 review, fixedQ budget and earlier KG opposing review; did not author field25 or common_charge_topology b45. Reused primitive compiler/native executor are source/execution dependencies, not independent derivations of those components.',
  'official_row_award_by_this_code':False})
if __name__=='__main__':main()
