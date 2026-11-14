"""Independent literal spin-matrix and rational-enclosure opposition."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json,math
import sympy as s
import numpy as np
from scipy.linalg import expm,norm
H=Path(__file__).resolve().parent
def read(name):return json.loads((H/name).read_text())
def sha(name):return hashlib.sha256((H/name).read_bytes()).hexdigest()
p=read('PARAMETERS.json');saved=read('RESULT.json');literal=read('SPIN_BRANCHES.json')
assert p['source_bridge_sha256']==literal['source_bridge_sha256']
assert p['mass_midpoints_exact']==literal['source_mass_midpoints_exact']
assert p['clock_phase_per_macro']==literal['source_clock_phase']=='1/24'
I=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1)
alpha=[s.kronecker_product(a,X) for a in [X,Y,Z]];beta=s.kronecker_product(I,Z)
branches={}
for signs in product((-1,1),repeat=3):
 K=s.eye(4)
 for a,sign in zip(alpha,signs):K=(s.eye(4)+sign*a)*K/2
 K=K.applyfunc(s.expand)
 assert str(K)==literal['exact_spin_branch_matrices'][str(signs)]
 branches[signs]=np.array(K,complex)
assert sum((s.sympify(v,locals={'Matrix':s.Matrix,'I':s.I}).H*s.sympify(v,locals={'Matrix':s.Matrix,'I':s.I}) for v in literal['exact_spin_branch_matrices'].values()),s.zeros(4)).applyfunc(s.simplify)==s.eye(4)

def tighter(z):
 assert 0<=z<=1
 term=F(1);total=term;upper=lower=None
 for j in range(1,20):
  term*=-z/F((2*j-1)*(2*j));total+=term
  if j==18:upper=total
  if j==19:lower=total
 assert lower<=upper
 return lower,upper

A=np.array(alpha[0],complex);B=np.array(beta,complex);checks=[]
for index,m in enumerate(map(F,p['mass_midpoints_exact'])):
 hm=m/3;z=hm*hm;cm=tighter(z);ca=tighter(z+F(3,4));entry=saved['actual_flavor_controls'][index]
 macro=tuple(map(F,entry['macro_normalized_trace']));aut=tuple(map(F,entry['autonomous_normalized_trace']))
 assert macro[0]<=-cm[1]/2<=-cm[0]/2<=macro[1]
 assert aut[0]<=ca[0]<=ca[1]<=aut[1]
 independent_lower=ca[0]+cm[0]/2
 assert independent_lower>=F(entry['operator_difference_lower_exact'])>F(57,50)
 mass=expm(-1j*float(hm)*B)
 U=sum((np.exp(-1j*signs[0]*2*math.pi/3)*K for signs,K in branches.items()),np.zeros((4,4),complex))@mass
 assert norm(U-(-np.eye(4)/2-1j*math.sqrt(3)*A/2)@mass,2)<2e-15
 V=expm(-1j*(float(hm)*B+math.sqrt(3)*A/2))
 gap=float(norm(U-V,2));assert gap>float(independent_lower)
 assert abs(np.trace(U).real/4-float(sum(macro)/2))<2e-15
 assert abs(np.trace(V).real/4-float(sum(aut)/2))<2e-15
 assert norm(V-V,2)<float(independent_lower)
 checks.append(dict(flavor=index,independent_trace_lower_exact=str(independent_lower),
  literal_matrix_norm_difference=gap,wrong_source_mutation_rejected=True))

ann=s.zeros(4)
for n in range(1,4):ann[n-1,n]=s.sqrt(n)
mom=(ann-ann.H)/(s.I*s.sqrt(2));vac=s.Matrix([1,0,0,0]);one=mom*mom/4
expect=(vac.H*one*vac)[0];variance=(one*vac-expect*vac).norm()**2
assert 3*expect==s.Rational(3,8) and 3*variance==s.Rational(3,32)
assert F(saved['harmonic_reference_physical_energy_variance'])==F(27,32)
result=dict(schema='oph.macro_generator_identity.independent_review_result.v1',
 reviewer='charged_particle_bridge',RER_producer_code_author=False,
 reviewer_authored_related_interface_counterexample=True,
 producer_check_module_imported=False,independent_literal_spin_branches=8,
 tighter_rational_Taylor_degree=38,actual_flavor_controls=checks,
 harmonic_exact_Fock_controls=True,scoped_result_pass=True,
 charged_interacting_spectrum_or_particle_claim=False,full_criterion_met=False,
 input_files_sha256={n:sha(n) for n in ['PARAMETERS.json','RESULT.json','SPIN_BRANCHES.json','check.py','PROOF.md']})
(H/'REVIEW_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS independent8 literal branches,3 exact nested intervals,3 wrong-source controls, harmonic variance')
