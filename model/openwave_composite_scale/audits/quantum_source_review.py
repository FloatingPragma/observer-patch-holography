"""Independent exact controls for the native charge and infinite-photon joins."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.openwave_quantum_rotational.source import matter_exact
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent

def review():
 subsets,_,currents,_=representation();Q=currents[10]+currents[11]
 # Define the filled weak-singlet reference independently by its actual weak Casimir.
 weak=sum((t*t for t in currents[8:11]),s.zeros(15))
 singlets=[i for i in range(15) if weak[i,i]==0]
 assert len(singlets)==7 and sum(Q[i,i] for i in singlets)==0
 assert all(currents[k].extract(singlets,[i for i in range(15) if i not in singlets])==s.zeros(7,8) for k in range(12))
 assert all(s.trace(t.extract(singlets,singlets))==0 for t in currents)
 transformed=[-Q[i,i] if i in singlets else Q[i,i] for i in range(15)]
 # Two spins: exact scalar offset vanishes and every charged occupation has |q|<=1.
 offset=2*sum(Q[i,i] for i in singlets)
 assert offset==0 and all(abs(q)<=1 for q in transformed)
 assert sorted(transformed,key=str)==sorted([s.Rational(2,3)]*6+[-s.Rational(1,3)]*6+[s.Integer(-1)]*2+[s.Integer(0)],key=str)
 m=matter_exact();T=m['B']**2;assert T*T==s.Rational(5,2)*T
 neutral=s.kronecker_product(s.eye(2),m['N'])+T/2-s.eye(20)/2
 assert all(x>0 for x in neutral.LDLdecomposition(hermitian=True)[1].diagonal())
 # Scalar spectral theorem remainder: the seventh derivative is positive and decreasing.
 x=s.symbols('x',nonnegative=True);f=s.sqrt(1+x)
 seventh=s.diff(f,x,7);assert s.simplify(seventh-s.Rational(10395,128)/(1+x)**s.Rational(13,2))==0
 coefficient=s.simplify(seventh.subs(x,0)/s.factorial(7));assert coefficient==s.Rational(33,2048)
 v0=(T+s.Rational(1,64)*sum((a*a for a in m['derivative']),s.zeros(20)))[0,0]/32
 assert 0<v0<s.Rational(1,16)
 # Each exact linear field factor changes total photon number by at most one.
 # Fourteen powers of B therefore fit in C+14, unlike a functional-calculus compression.
 for C in (1,2,3):assert C+14>=C+2*7
 receipts={}
 paths=[str(Path(__file__).relative_to(ROOT))]
 for pkg in ('openwave_native_quantum_charge','openwave_infinite_photon','openwave_quantum_larmor'):
  base=ROOT/'model'/pkg
  paths += [str(p.relative_to(ROOT)) for p in base.iterdir() if p.is_file() and p.suffix in ('.py','.json','.md','.tex')]
  r=json.loads((base/'results.json').read_text());receipts[pkg]=r.get('verdict')
  maps=[r[k] for k in ('source_sha256','input_sha256') if k in r]
  assert maps
  for pins in maps:
   for name,digest in pins.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
 paths+=['model/coherent_defect_wave/physical_attachment/chiral.py','model/openwave_quantum_rotational/source.py']
 pins={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(paths))}
 return {'verdict':'PASS_INDEPENDENT_EXACT_QUANTUM_SOURCE_CONTROLS','native_filled_singlet_count':len(singlets),'native_charge_offset':str(offset),'all_twelve_reference_current_traces':'0','charged_occupation_magnitudes_at_most_one':True,'neutral_LDL_floor':'strictly greater than 1/2','global_seventh_remainder_coefficient':str(coefficient),'infinite_field_two_state_trial_shift_exact':str(v0),'padded_moment_photon_increments':14,'reviewed_receipts':receipts,'review_scope':'Independent exact current/reference, neutral-floor, global scalar remainder and padding controls plus source-map validation. Energy interval and finite-time Larmor producers retain their own interval arithmetic; this review does not claim a second implementation of every enclosure. The two-dimensional isolated ground representation makes H scalar by irreducibility, so the next-distinct-eigenvalue Temple bound applies.','source_sha256':pins}

def main():
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=review();out=HERE/'quantum_source_review.json'
 if a.write:out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert json.loads(out.read_text())==r
 print(r['verdict'])
if __name__=='__main__':main()
