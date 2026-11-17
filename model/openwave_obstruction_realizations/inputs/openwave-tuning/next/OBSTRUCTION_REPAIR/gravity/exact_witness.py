"""Exact nonzero q^2 witness inside the existing Dirac source family."""
import json
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
X=s.Matrix([[0,1],[1,0]])
Y=s.Matrix([[0,-s.I],[s.I,0]])
Z=s.diag(1,-1)
I=s.eye(4)
alpha=[s.kronecker_product(p,X) for p in (X,Y,Z)]
beta=s.kronecker_product(s.eye(2),Z)


def clean(a): return a.applyfunc(s.expand)
def mul(a,b): return clean(a*b)


f=[(I-s.I*a)/s.sqrt(2) for a in alpha]
coin=s.sqrt(399)*I/20-s.I*beta/20
U=mul(mul(mul(f[2],f[1]),f[0]),coin)
dU=mul(mul(mul(f[2],f[1]),(-I-s.I*alpha[0])/s.sqrt(2)),coin)
H=clean(s.I*(U-U.adjoint())/2)
D=clean(s.I*(dU-dU.adjoint())/2)
K=clean(mul(H,D)+mul(D,H))
H2=mul(H,H);H3=mul(H2,H)
tr=lambda a:s.simplify(s.trace(a))
values=[tr(mul(K,K)),tr(mul(mul(mul(H,K),H),K)),
        tr(mul(mul(mul(H3,K),H),K)),tr(mul(mul(mul(H3,K),H3),K))]
expected=[s.Rational(1,400),s.Rational(599,64000000),
          s.Rational(1037401,102400000000),s.Rational(446561,40960000000)]
assert values==expected
assert clean(mul(H2,H2)-s.Rational(1199,800)*H2+s.Rational(319201,640000)*I)==s.zeros(4)
# sign(H)=[(lambda_+ + lambda_- + sqrt(lambda_+lambda_-))H-H^3]
#         /[sqrt(lambda_+lambda_-)*sqrt(lambda_++lambda_-+2sqrt(lambda_+lambda_-))].
d=s.sqrt(319201)/800
A=s.Rational(1199,800)+d
invariant=s.radsimp(s.simplify(values[0]-(A*A*values[1]-2*A*values[2]+values[3])/(d*d*(s.Rational(1199,800)+2*d))))
target=s.Rational(80199,32240000)+s.Rational(79401,10291040240000)*s.sqrt(319201)
assert s.simplify(invariant-target)==0
out=dict(source='Existing ordered Dirac source, k=(pi/4,pi/4,pi/4), sin(mu)=1/20',
         mu=str(s.asin(s.Rational(1,20))),trace_values=list(map(str,values)),
         positive_interband_invariant=str(invariant),invariant_decimal=str(s.N(invariant,30)),
         local_chi_q2_lower_bound=str(invariant/64),
         conclusion='Strictly positive q^2 susceptibility on an open momentum neighborhood; not an Einstein coefficient',
         all_checks_pass=True)
(HERE/'EXACT_WITNESS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
