#!/usr/bin/env python3
"""Exact conditional matter identification from the existing native response.

The native M6 generators are imported unchanged. The documented spin lift,
trace redistribution, auxiliary exterior construction and proper-code contract
are declared inputs. No field-name table, empirical mass or coupling is used
by the eigencharge/Casimir decoder. Internal representation is not physical
Weyl statistics, propagation, a selected generation count or a mass pole.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib,importlib.util,itertools,json,sys
from collections import Counter
from pathlib import Path
sys.dont_write_bytecode=True
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PUBLIC=ROOT
I=s.I

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
 sys.modules[name]=m;spec.loader.exec_module(m);return m

def clean(a):return a.applyfunc(s.expand)
def zero(a):return all(s.expand(v)==0 for v in a)
def unit(n,i,j):
 a=s.zeros(n);a[i,j]=1;return a

def spin_lift(d):
 """Existing matter_transport: A plus spin(hat w), weak trace -Tr(A)/2."""
 a=d[:3,:3];b=d[3:,3:]
 assert zero(b+b.T) and s.trace(b)==0
 w=(b[2,1],b[0,2],b[1,0])
 pauli=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1))
 weak=sum((-I*w[k]*pauli[k]/2 for k in range(3)),s.zeros(2))-s.trace(a)*s.eye(2)/2
 return s.diag(a,weak)

def fock_data():
 basis=sorted(itertools.chain.from_iterable(itertools.combinations(range(5),k) for k in range(6)),key=lambda x:(len(x),x))
 index={x:i for i,x in enumerate(basis)}
 creation=[]
 for k in range(5):
  a=s.zeros(32)
  for j,x in enumerate(basis):
   if k not in x:a[index[tuple(sorted((*x,k)))],j]=(-1)**sum(v<k for v in x)
  creation.append(s.SparseMatrix(a))
 for a in range(5):
  for b in range(5):
   assert creation[a]*creation[b]+creation[b]*creation[a]==s.zeros(32)
   assert creation[a].H*creation[b]+creation[b]*creation[a].H==s.eye(32)*int(a==b)
 products={(i,j):creation[i]*creation[j].H for i in range(5) for j in range(5)}
 def dg(a):
  out=s.SparseMatrix(32,32,{})
  for (i,j),v in a.todok().items():out+=v*products[i,j]
  return clean(out)
 return basis,dg

def restrict(a,indices):return a.extract(indices,indices)
def exact_values(a):
 roots=s.polys.polytools.ground_roots(a.charpoly().as_expr())
 assert sum(roots.values())==a.rows
 return sorted(roots,key=lambda x:float(x))

def projectors(a):
 vals=exact_values(a);out={}
 for v in vals:
  p=s.eye(a.rows)
  for w in vals:
   if w!=v:p=p*(a-w*s.eye(a.rows))/(v-w)
  p=clean(p);assert p*p==p and p.H==p
  out[v]=p
 assert sum(out.values(),s.zeros(a.rows))==s.eye(a.rows)
 return out

def decode(operators):
 """Only operators in an arbitrary fixed basis, no Fock subsets/field names."""
 q,c3,c2,h1,h2,t3=(operators[k] for k in ('q','color_C2','weak_C2','H1','H2','T3'))
 weights1,weights2=projectors(h1),projectors(h2)
 rows=[]
 for charge,p in projectors(q).items():
  dim=s.trace(p);assert dim.is_Integer and dim>0
  color=s.simplify(s.trace(p*c3)/dim);weak=s.simplify(s.trace(p*c2)/dim)
  assert p*c3==color*p and p*c2==weak*p
  weights=[]
  for a,pa in weights1.items():
   for b,pb in weights2.items():
    mult=s.trace(p*pa*pb)
    assert mult.is_Integer and mult>=0
    if mult:weights.append({'H1':str(a),'H2':str(b),'multiplicity':int(mult)})
  highest=[w for w in weights if s.Rational(w['H1'])>=0 and s.Rational(w['H2'])>=0]
  assert len(highest)==1
  weak_weights=[]
  for v,pv in projectors(t3).items():
   mult=s.trace(p*pv)
   if mult:weak_weights.append({'T3':str(v),'multiplicity':int(mult),'electric_Q':str(v+charge/6)})
  rows.append({'q6Y':str(charge),'Y':str(charge/6),'dimension':int(dim),
   'color_C2':str(color),'weak_C2':str(weak),
   'color_highest_weight':[int(highest[0]['H1']),int(highest[0]['H2'])],
   'color_Cartan_weight_multiplicities':weights,'weak_and_electric_weights':weak_weights,
   'complete_uniform_state_probability':str(dim/q.rows)})
 return rows

def encoded_controls(source):
 """Exact inversion onto actual first-two-coordinate proper code."""
 gs=source.source_generators()
 gram=s.Matrix(12,12,lambda i,j:s.expand(-s.trace(gs[i]*gs[j])))
 inverse=gram.inv().applyfunc(s.simplify)
 code=s.eye(6)[:,:2]
 pauli=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1))
 rows=[]
 for name,sigma in zip(('x','y','z'),pauli):
  target=s.zeros(6);target[:2,:2]=-I*sigma/2
  rhs=s.Matrix([-s.trace(d*target) for d in gs])
  h=(inverse*rhs).applyfunc(s.simplify)
  actual=clean(sum((c*d for c,d in zip(h,gs)),s.zeros(6)))
  assert actual==target and code.H*actual*code==-I*sigma/2
  rows.append({'axis':name,'port_coefficients':[str(v) for v in h]})
 # Existing complete proper-code service retains complement failure.
 ks=[s.Matrix(k).applyfunc(lambda x:s.nsimplify(x)) for k in source.code_transfer_kraus(2)]
 assert sum((k.H*k for k in ks),s.zeros(6))==s.eye(6)
 assert ks[0]*code==s.eye(3)[:,:2]
 assert all(k*code==s.zeros(3,2) for k in ks[1:])
 weak_native=s.zeros(6);weak_native[3:6,3:6]=s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
 assert code.H*weak_native*code==s.zeros(2)
 return {'exact_source_control_rows':rows,'complete_code_transfer_kraus_count':len(ks),
  'success_effect_rank':2,'failure_effect_rank':4,'completeness_exact':True,
  'code_services_are_declared_not_unconstrained_entropy_outputs':True,
  'chosen_code_is_in_first_u3_block_while_native_weak_ideal_is_in_second_block':True,
  'native_weak_ideal_restricted_to_chosen_code_is_zero':True,
  'programmable_spinor_controls_are_not_an_inherited_same_current_intertwiner':True,
  'new_full_15_state_native_pulse_experiment_executed':False}

def sparse(a):
 return {'shape':list(a.shape),'entries':[[i,j,str(v)] for (i,j),v in a.todok().items()]}

def main():
 source=load('sm_identification_native_source',PUBLIC/'code/m1_source_realization/model.py')
 ds=source.source_generators()
 gauge=load('sm_identification_probe_decoder',HERE/'gauge_particle_identification.py')
 observed=gauge.decoded_generators()
 assert all(zero(d-s.trace(d)*s.eye(6)/6-o) for d,o in zip(ds,observed,strict=True))
 assert clean(sum(ds,s.zeros(6)))==s.diag(I,I,I,0,0,0)
 lifted=[spin_lift(d) for d in ds]
 assert s.Matrix.hstack(*(d.reshape(25,1) for d in lifted)).rank()==12
 for a,b in itertools.combinations(range(12),2):
  assert zero(spin_lift(clean(ds[a]*ds[b]-ds[b]*ds[a]))-(lifted[a]*lifted[b]-lifted[b]*lifted[a]))
 # Same map on the central source direction; integer normalization is derived
 # by clearing denominators and making the weak charge positive (conjugation).
 center=clean(sum(lifted,s.zeros(5))/I)
 assert center==s.diag(1,1,1,-s.Rational(3,2),-s.Rational(3,2))
 denominators=[v.q for v in center.diagonal()]
 scale=s.ilcm(*denominators)
 ints=[int(v*scale) for v in center.diagonal()]
 gcd=s.igcd(*ints);qone=-scale*center/gcd
 assert qone==s.diag(-2,-2,-2,3,3)
 basis,dg=fock_data()
 even=[i for i,x in enumerate(basis) if len(x)%2==0 and x]
 odd=[i for i,x in enumerate(basis) if len(x)%2==1 and len(x)<5]
 assert len(even)==len(odd)==15
 # Invariants are constructed from the lifted Lie-algebra blocks rather than
 # supplied five-field tables. Color C2 is normalized with Tr(t_a t_b)=delta/2.
 e={ (i,j):dg(unit(5,i,j)) for i in range(3) for j in range(3)}
 nc=sum((e[i,i] for i in range(3)),s.zeros(32))
 c3=sum((e[i,j]*e[j,i]/2 for i in range(3) for j in range(3)),s.zeros(32))-nc*nc/6
 sigmas=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1))
 ts=[]
 for sigma in sigmas:
  a=s.zeros(5);a[3:,3:]=sigma/2;ts.append(dg(a))
 c2=sum((t*t for t in ts),s.zeros(32))
 q=dg(qone);h1=e[0,0]-e[1,1];h2=e[1,1]-e[2,2]
 full={'q':q,'color_C2':c3,'weak_C2':c2,'H1':h1,'H2':h2,'T3':ts[2]}
 # Invariance and Lie functor on a spanning canonical 12-generator set.
 one=[]
 for i,j in itertools.combinations(range(3),2):one.extend((unit(5,i,j)-unit(5,j,i),I*(unit(5,i,j)+unit(5,j,i))))
 one.extend((I*(unit(5,0,0)-unit(5,1,1)),I*(unit(5,1,1)-unit(5,2,2))))
 for sigma in sigmas:
  a=s.zeros(5);a[3:,3:]=-I*sigma/2;one.append(a)
 one.append(I*qone)
 fullrep=[dg(a) for a in one]
 for a,b in itertools.combinations(range(12),2):assert zero(fullrep[a]*fullrep[b]-fullrep[b]*fullrep[a]-dg(one[a]*one[b]-one[b]*one[a]))
 p=s.diag(*[int(i in even) for i in range(32)])
 assert all(p*x==x*p for x in fullrep)
 evenops={k:restrict(v,even) for k,v in full.items()}
 oddops={k:restrict(v,odd) for k,v in full.items()}
 # Independent fixed rational change of basis mixes unlike gauge sectors.
 z=s.eye(15)
 for a,b in ((0,5),(1,12),(2,9),(4,10),(6,14)):
  r=s.eye(15);r[a,a]=r[b,b]=s.Rational(3,5);r[a,b]=s.Rational(4,5);r[b,a]=-s.Rational(4,5);z=r*z
 assert z.T*z==s.eye(15)
 decoded=decode({k:z*v*z.T for k,v in evenops.items()})
 assert decoded==decode(evenops)
 conjugate=decode(oddops)
 assert {r['q6Y'] for r in decoded}.isdisjoint({r['q6Y'] for r in conjugate})
 # Complete joint Cartan/charge record law for maximally mixed selected code.
 records=[]
 for i in even:
  records.append({'q6Y':str(q[i,i]),'H1':str(h1[i,i]),'H2':str(h2[i,i]),
   'T3':str(ts[2][i,i]),'probability':'1/15'})
 assert sum((s.Rational(r['probability']) for r in records),s.Rational(0))==1
 # Omitted neutral invariant lines are visible, not dropped as failed shots.
 invariants=[i for i in range(32) if all(x[:,i]==s.zeros(32,1) for x in fullrep)]
 assert [basis[i] for i in invariants]==[(),(0,1,2,3,4)]
 a,b=s.symbols('a b');general=dg(s.diag(a,a,a,b,b))
 gravity=s.expand(s.trace(restrict(general,even)))
 assert gravity==24*a+16*b
 native_z=s.zeros(6);native_z[3:6,3:6]=s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
 assert s.exp(2*s.pi*native_z)==s.eye(6)
 lift_z=spin_lift(native_z)
 assert s.exp(2*s.pi*lift_z)==s.diag(1,1,1,-1,-1)
 artifact=HERE/'results/representation_identification_operators.json'
 artifact.write_text(json.dumps({'basis_subsets':basis,'even_indices':even,'odd_indices':odd,
  'native_pulse_generators':[sparse(d) for d in ds],
  'spin_lifted_one_body_generators':[sparse(d) for d in lifted],
  'even_identification_operators':{k:sparse(v) for k,v in evenops.items()},
  'basis_change':sparse(z)},indent=2)+'\n', encoding='utf-8', newline='\n')
 paths=['code/source_selection_model/response.json','code/source_selection_model/response.py',
 'code/source_selection_model/DERIVATION.md','code/m1_source_realization/model.py',
 'docs/research/COHERENT_SOURCE_CLOCKS.md','code/a5_closure/super_tannakian_matter_lift_certificate.py',
 'code/a5_closure/manifests/super_tannakian_matter_reference.json',
 'Lean/Screen/ExteriorSelection.lean','Lean/Screen/MatterGrammarIndexBridge.lean',
 'Lean/Screen/A5FamilyBand.lean','Lean/Screen/FermionSectorAssembly.lean',
 'paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex',
 'paper/tex_fragments/STANDARD_MODEL_DERIVATION_SCOPE.tex','code/b15_matter_freeze/matter_class_freeze_v1.md']
 sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
 bindings={str((PUBLIC/p).relative_to(ROOT)):sha(PUBLIC/p) for p in paths}
 bindings[str((HERE/'gauge_particle_identification.py').relative_to(ROOT))]=sha(HERE/'gauge_particle_identification.py')
 result={'status':'EXACT_CONDITIONAL_SOURCE_SPIN_LIFT_AND_MATTER_IDENTIFICATION_PASS',
 'producer_sha256':sha(__file__),'source_sha256':bindings,
 'artifacts':{artifact.name:sha(artifact)},'empirical_mass_or_coupling_inputs':[],
 'native_source_to_lift':{'source_directions':12,'source_size':6,'lift_size':5,
  'all_generators_match_independent_matrix_probe_decoder_modulo_scalar':True,
  'all_66_bracket_pairs_exact':True,'source_center_on_lift':[str(x) for x in center.diagonal()],
  'lift_is_faithful_at_Lie_algebra_level':True,
  'primitive_integer_charge_from_center':[str(x) for x in qone.diagonal()],
  'central_trace_redistribution_declared':True,'spin_half_lift_declared':True,
  'native_2pi_weak_rotation':'I6','lifted_2pi_weak_rotation':'diag(1,1,1,-1,-1)',
  'lift_is_not_a_subrepresentation_of_original_native_M6':True},
 'auxiliary_exterior':{'CAR_checks':50,'dimension':32,'even_minus_vacuum_rank':15,
  'odd_minus_top_rank':15,'canonical_dgamma_bracket_checks':66,
  'full_invariant_lines':[list(basis[i]) for i in invariants],
  'even_auxiliary_CAR_parity':1,'odd_auxiliary_CAR_parity':-1,
  'physical_all_left_Weyl_fermion_grading_inferred_from_auxiliary_parity':False,
  'conjugate_selected_modules_have_zero_current_intertwiner_space':True,
  'general_even_package_gravitational_trace':str(gravity),
  'anomaly_balance_and_primitive_integrality_separate_contracts':True},
 'decoded_even_internal_multiplets':decoded,'decoded_conjugate_internal_multiplets':conjugate,
 'complete_cartan_charge_record_law':records,
 'charge_readout_contract':{'integer_charge':'primitive integral central character in the declared lifted representation',
  'Y_convention':'q/6; this fixes the normalization convention, not the gauge coupling',
  'electric_Q_convention':'T3+Y; physical unbroken-generator/scalar identification remains a separate premise',
  'uniform_selected_code_state_is_declared_for_the_ideal_record_law':True,
  'source_chooses_central_compactification_or_global_form':False},
 'decoder_uses_supplied_field_name_table':False,'fixed_exact_basis_change_identification_invariant':True,
 'proper_code_implementation':encoded_controls(source),
 'countercontrols':{'retain_even_vacuum_adds_gauge_neutral_invariant_state':True,
  'neutral_invariant_not_removed_by_charge_traces':True,
  'charge_conjugate_representation_also_realized':True,
  'native_weak_2pi_identity_becomes_nontrivial_on_lift':True},
 'scope':{'source_selected_entire_matter_module':False,'physical_Spin_or_fermion_statistics_selected':False,
  'propagating_particle_poles_computed':False,'generation_count_selected':False,
  'laboratory_species_identified':False,'new_native_detector_execution_claimed':False,
  'positive_result':'Five internal gauge multiplets recovered from exact source-linked conditional representation operators without a field-name table.'}}
 out=HERE/'results/representation_identification.json';out.write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
 print(json.dumps({'status':result['status'],'decoded_even_internal_multiplets':decoded},indent=2))

if __name__=='__main__':main()
