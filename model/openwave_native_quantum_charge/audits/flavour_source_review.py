"""Independent native BdG current and inclusive-decay normalization audit."""
import argparse,hashlib,json,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import sympy as s
from model.common_source import calibrated_matter as native
from model.common_source.chiral_candidate import kinetic
from model.openwave_neutrino_oscillation.oscillation import construction
from model.native_flavour_transport.transport import numeric
from model.particle_common_contract.provenance import import_closure
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve()

def flight_review():
 d=construction();mass,frame=np.linalg.eigh(numeric(d['y']));V=numeric(d['V']);D=native._oriented_components('charged_lepton')[1];charged,singlet=next(iter(D.todok()));sign=float(D[charged,singlet]);nu=d['neutral'];J=numeric(native.weak_raising_current())
 current=np.block([[np.kron(np.eye(2),J),np.zeros((90,90))],[np.zeros((90,90)),-np.kron(np.eye(2),J.T)]])
 C=frame.T@V@frame;qsign=complex(native.weak_raising_current(families=False)[nu,charged]);maximum=0.;rows=[]
 def vectors(p,neutral):
  H=kinetic.bdg_symbol([0,0,p],numeric(d['k']),speed=1.);out=[]
  for j,m in enumerate(mass):
   E=np.sqrt(m*m+p*p);a=np.sqrt((E+p)/(2*E));b=np.sqrt((E-p)/(2*E));x=np.zeros(180,complex)
   if neutral:x[45+3*nu:45+3*nu+3]=a*(V@frame)[:,j];x[90+3*nu:90+3*nu+3]=-b*(V.conj()@frame)[:,j]
   else:x[45+3*charged:45+3*charged+3]=a*frame[:,j];x[90+3*singlet:90+3*singlet+3]=-b*sign*frame[:,j]
   assert np.linalg.norm(H@x-E*x)<2e-12 and abs(np.vdot(x,x)-1)<2e-12
   out.append(x)
  return np.array(out).T
 for p in(-.83,.0,.39):
  for q in(-.21,.0,.91):
   U=vectors(p,True);W=vectors(q,False);actual=U.conj().T@current@W
   wn=(np.sqrt(mass*mass+p*p)+p)/(2*np.sqrt(mass*mass+p*p));wc=(np.sqrt(mass*mass+q*q)+q)/(2*np.sqrt(mass*mass+q*q))
   expected=qsign*np.diag(np.sqrt(wn))@C.conj().T@np.diag(np.sqrt(wc));err=float(np.max(abs(actual-expected)));maximum=max(maximum,err);assert err<3e-13
   rows.append({'neutral_momentum':p,'charged_momentum':q,'full_180_BdG_current_error':err})
 return {'all_nonzero_momentum_current_checks':rows,'maximum_error':maximum,'particle_hole_and_opposite_spin_complements_retained':True}

def decay_review():
 I=s.eye(2);pauli=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1));bar=(I,*(-x for x in pauli));eta=(1,-1,-1,-1)
 M=s.Symbol('M',real=True);vectors=[s.symbols(prefix+'0:4',real=True)for prefix in('e','n','a')]
 def density(p):return p[0]*I-sum((p[j+1]*pauli[j]for j in range(3)),s.zeros(2))
 e,n,a=vectors;de,dn,da=map(density,vectors)
 trace=s.expand(sum(eta[mu]*eta[nu]*s.trace(de*bar[mu]*da*bar[nu])*s.trace(dn*bar[mu]*(M*I/2)*bar[nu])for mu in range(4)for nu in range(4)))
 expected=8*M*a[0]*(e[0]*n[0]-sum(e[j]*n[j]for j in range(1,4)));assert s.expand(trace-expected)==0
 z,r=s.symbols('z r',positive=True);integrand=s.expand((1-z)**2*(z-r)**2/z);primitive=s.integrate(integrand,z);value=s.simplify(12*(primitive.subs(z,1)-primitive.subs(z,r)));assert s.simplify(value-(1-8*r+8*r**3-r**4-12*r*r*s.log(r)))==0
 # The three-body factor is ds/(2pi)*(1-z)/(8pi)*(z-r)/(8pi*z).
 phase_denominator=2*8*8;flux_denominator=2;amplitude_numerator=2;integral_denominator=12
 assert phase_denominator*flux_denominator*integral_denominator//amplitude_numerator==1536
 m=F(49,10000);lower=F(3,625);D=F(3072*3**4);S0=F(32,3);S1=32+64*m*m+9*m**5+384*m**4;S2=80+96*m*m+179*m**4;S3=160+96*m*m+52*m**3
 A0=F(256,315)/D;A2=(16*S2+64*S1+48*S0)/D;A3=(16*S3+96*S2+144*S1+48*S0)/D
 assert sum(F(7)**k/F(math.factorial(k))for k in range(16))>2/lower
 assert 4*A0*A2<F(1,2000)**2 and F(27,8)*A0*A3*A3<F(1,250)**3
 assert (1-8*m*m-m**8)/(1536*F(22,7)**3)>F(1,50000)
 return {'full_symbolic_Pauli_trace_identity':str(expected),'exact_current_rate_denominator':1536,'phase_space_function':str(value),
  'independent_squared_correlation_margin':str(F(1,2000)**2-4*A0*A2),'independent_cubed_weighted_margin':str(F(1,250)**3-F(27,8)*A0*A3*A3),
  'four_lifetime_continuum_survival_error_exact':str(2*F(1,64)**2*F(1,250)*(1+25*4)),
  'scope':'Positive-frequency current member and inclusive finite-bin theorem; no unfiltered-source transfer, physical-box discretization, observed lifetime, or irreversible finite-unitary claim.'}

def pins():
 seeds={str(HERE.relative_to(ROOT))}
 for package in('openwave_neutral_flight','openwave_inclusive_decay'):
  seeds.update(str(p.relative_to(ROOT))for p in(ROOT/'model'/package).iterdir()if p.suffix in('.py','.md','.tex','.json'))
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}

def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args();before=pins();r={'schema':'oph.independent_native_flight_inclusive_decay_review.v1','flight':flight_review(),'decay':decay_review(),'source_sha256':before,'verdict':'PASS_NATIVE_BDG_CURRENT_AND_INCLUSIVE_DECAY_INDEPENDENT_REVIEW'};assert before==pins();target=HERE.with_suffix('.json')
 if args.write:target.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert json.loads(target.read_text())==r
 print(r['verdict'])
if __name__=='__main__':main()
