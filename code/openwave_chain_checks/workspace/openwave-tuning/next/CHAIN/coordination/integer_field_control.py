"""Exact integer-versus-real dual-field minimization control on the carrier."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
source=ROOT/'reverse-engineering-reality/code/openwave_row1_native/source_packet.json'
packet=json.loads(source.read_text());edges=packet['seams'];size=12
B=[[int(v==a)-int(v==b) for a,b in edges] for v in range(size)]
L=[[sum(B[a][e]*B[b][e] for e in range(len(edges))) for b in range(size)] for a in range(size)]
q=[row[0] for row in B]
# Fix the additive potential at the last vertex; exact rational elimination.
a=[[Q(L[i][j]) for j in range(size-1)]+[Q(q[i])] for i in range(size-1)]
for k in range(size-1):
 pivot=next(i for i in range(k,size-1) if a[i][k]);a[k],a[pivot]=a[pivot],a[k]
 z=a[k][k];a[k]=[x/z for x in a[k]]
 for i in range(size-1):
  if i!=k:
   z=a[i][k];a[i]=[x-z*y for x,y in zip(a[i],a[k])]
u=[row[-1] for row in a]+[Q(0)]
f=[u[a]-u[b] for a,b in edges]
assert [sum(B[v][e]*f[e] for e in range(len(edges))) for v in range(size)]==q
energy=sum(x*x for x in f)/2
assert energy==Q(11,60) and any(x.denominator!=1 for x in f)
# Every admissible integer field is nonzero, hence its integer squared norm >=1.
# Unit flow on the chosen adjacent seam is admissible and attains equality.
integer_energy=Q(1,2)
result={'schema':'oph.integer_field_minimization_control.v1','source_sha256':{str(source.relative_to(ROOT)):hashlib.sha256(source.read_bytes()).hexdigest(),str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},'unit_dipole_seam':edges[0],'charge':q,'real_minimizer':[str(x) for x in f],'real_minimum_energy':str(energy),'integer_minimum_energy':str(integer_energy),'checks':{'exact_real_gauss_solution':True,'real_minimizer_not_integer':True,'integer_lower_bound_attained_by_unit_flow':True,'integer_equals_real_minimizer_claim_rejected':integer_energy>energy},'scope':'For integer face loads x, F-boundary(x) remains integer and cannot equal this real Coulomb minimizer in one configuration. Real-valued face variables, an ensemble/time-averaged field or a proved refinement are different explicit domains. Topological degree is unaffected by this distinction; no force-law conclusion.', 'official_rows_changed':[]}
(HERE/'INTEGER_FIELD_CONTROL.json').write_text(json.dumps(result,indent=2)+'\n');print(result['real_minimum_energy'],result['integer_minimum_energy'],result['checks'])
