"""Exact native-current two-photon coefficient and finite-time lower bound."""
from fractions import Fraction as F
import sympy as s
from .radiation import apply

def certificate():
    sx=s.Matrix([[0,1],[1,0]]);sy=s.Matrix([[0,-s.I],[s.I,0]]);tz=s.diag(1,-1)
    frames=[]
    for pi in (s.Integer(0),s.Rational(1,2),-s.Rational(1,2)):
        cols=[]
        for sign in (-1,1):
            for spin in (0,1):
                h=1 if spin==0 else -1;E=sign*s.sqrt(1+pi*pi)
                small=s.Matrix([1,E+h*pi]);small/=s.sqrt((small.T*small)[0])
                v=s.zeros(4,1);v[spin]=small[0];v[spin+2]=small[1];cols.append(v)
        frames.append(s.Matrix.hstack(*cols))
    matrices={}
    for q in (1,2):
        for pol,sigma in enumerate((sx,sy)):
            J=s.zeros(12);alpha=-s.kronecker_product(tz,sigma)
            for p in range(3):
                out=(p-q)%3;J[4*out:4*out+4,4*p:4*p+4]=(frames[out].H*alpha*frames[p]).applyfunc(s.simplify)
            matrices[q,pol]=J
    ground=sum(1<<i for i in (0,1,4,5,8,9));pair={}
    for spin in(0,1):
        out,sign=apply(ground,spin+2,spin);pair[out]=sign/s.sqrt(2)
    def act(J,state):
        result={}
        for i in range(12):
            for j in range(12):
                if J[i,j]==0:continue
                for mask,amp in state.items():
                    out=apply(mask,i,j)
                    if out is not None:result[out[0]]=result.get(out[0],0)+out[1]*J[i,j]*amp
        return {k:s.simplify(v)for k,v in result.items()if v!=0}
    def overlap_after(A,B):
        v=act(B,pair);w=act(A.H,{ground:s.Integer(1)})
        return s.simplify(sum(s.conjugate(w.get(k,0))*v[k]for k in v))
    rows=[]
    for x in range(2):
        for y in range(2):
            A=matrices[1,x];B=matrices[2,y]
            amp=s.simplify(overlap_after(A,B)+overlap_after(B,A))
            expected=0 if x==y else (1 if x==0 else-1)*2*s.I*s.sqrt(s.Rational(2,5))
            assert s.simplify(amp-expected)==0,(x,y,amp,expected)
            rows.append({'polarization_plus':x,'polarization_minus':y,'H_squared_amplitude_divided_by_eta_squared':str(amp)})
    # With cap2, ||H0||<21, ||eta V||<=48 eta sqrt(2)<4 at eta1/20.
    # B25 is an outward rational full-sector bound. Unitary Taylor
    # integral remainder for order2 is ||R3||<=B^3 t^3/6.
    eta=F(1,20);t=F(1,10**8);B=F(25)
    # 2sqrt(2/5)>6/5, so each of the two amplitudes has leading
    # magnitude >=3 eta^2 t^2/5. Bound both via one uniform norm error.
    leading=F(3,5)*eta**2*t**2;tail=B**3*t**3/6
    assert leading>tail
    prob=2*(leading-tail)**2
    return {'exact_opposite_momentum_two_photon_coefficients':rows,
        'full_cap2_norm_bound':'25 native masses','eta':'1/20','time_in_inverse_native_mass':str(t),
        'each_leading_amplitude_lower':str(leading),'Taylor_amplitude_error_upper':str(tail),
        'original_two_photon_radiation_probability_lower':str(prob),
        'original_two_photon_radiation_probability_lower_float':float(prob),
        'uses_finite_boson_cap_not_canonical_infinite_CCR':True}
if __name__=='__main__':
 import json
 print(json.dumps(certificate(),indent=2))
