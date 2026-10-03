"""Independent one-loop complex top pole via Feynman-parameter quadrature."""
import cmath,json,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

def A(x,q):return x*(math.log(x/q**2)-1) if x else 0.
def B(x,y,s,q):
    coeff=[s,x-y-s,y]
    roots=sorted(r.real for r in np.roots(coeff) if abs(r.imag)<1e-10 and 0<r.real<1)
    points=[0,*roots,1]
    val=0j
    for lo,hi in zip(points,points[1:]):
        def f(z):return s*z*z+(x-y-s)*z+y
        re=-quad(lambda z:math.log(abs(f(z))/q**2),lo,hi,epsabs=3e-12,epsrel=3e-12)[0]
        im=math.pi*(hi-lo) if f((hi+lo)/2)<0 else 0
        val+=complex(re,im)
    return val

def evaluate(w):
 q,v,gy,g,g3,yt,lam=(w[k] for k in ('Q','v','gY','g2','g3','yt','lambda'))
 t=yt*yt*v*v/2;h=2*lam*v*v;W=g*g*v*v/4;Z=(g*g+gy*gy)*v*v/4
 at,ah,aw,az=[A(x,q) for x in (t,h,W,Z)]
 b0w,bht,btz=B(0,W,t,q),B(h,t,t,q),B(t,Z,t,q)
 bfv0w=(W-t)*b0w+aw+t+((W*t-t*t)*b0w-t*aw)/(2*W)
 bfvtz=(Z-t)*btz+az-at+t
 dq=g3*g3*t*(32/3-8*math.log(t/q**2))
 de=8/9*(g*g*gy*gy/(g*g+gy*gy))*(t-3*at)+yt*yt/2*((h-4*t)*bht+ah-2*at)+g*g/2*bfv0w
 de+=((g*g+gy*gy)/4-2*gy*gy/3+8*gy**4/(9*(g*g+gy*gy)))*bfvtz
 de+=4*gy*gy/3*(2*gy*gy/(3*(g*g+gy*gy))-.5)*t*(3*btz-2)
 pole=cmath.sqrt(t+(dq+de)/(16*math.pi**2))
 d1=(3*g**4+2*g*g*gy*gy+gy**4)*v*v/8+3*lam*ah-6*yt*yt*at+3*g*g*aw/2+3*(g*g+gy*gy)*az/4
 return {'mass_GeV':pole.real,'width_GeV':-2*pole.imag,'m2_minimum_1loop':-lam*v*v-d1/(16*math.pi**2)}

def verify(receipt):
    """Check independent one-loop top complex pole and potential minimum.

    This is an analytic-formula numerical check, not certification of SMDR's
    two-loop library, whose separately pinned full replay remains required.
    """
    rows=receipt['forward']['rows']
    if len(rows)!=9:raise ValueError('incomplete scale/order control menu')
    for mult in (1,2,4):
        selected=[r for r in rows if r['scale_multiplier']==mult and r['matching_order']==1]
        if len(selected)!=1:raise ValueError('missing or duplicated one-loop row')
        row=selected[0]; predicted=evaluate(row['working'])
        for key in ('mass_GeV','width_GeV'):
            actual=row['top_method0'][key]
            if not math.isfinite(actual) or not math.isclose(predicted[key],actual,rel_tol=2e-10,abs_tol=2e-8):
                raise ValueError(f'independent one-loop top {key} mismatch')
        if mult==1:
            for order,value in [('0',-row['working']['lambda']*row['working']['v']**2),('1',predicted['m2_minimum_1loop'])]:
                if not math.isclose(value,receipt['forward']['m2_minimum_controls'][order],rel_tol=2e-10,abs_tol=2e-7):
                    raise ValueError('independent tadpole minimum mismatch')
    return True


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt',nargs='?',type=Path,default=Path(__file__).resolve().parents[1]/'runs/calibration/source_ew_vev_matching.json')
    args=parser.parse_args();verify(json.loads(args.receipt.read_text()))
    print('Independent one-loop top pole and tadpole: PASS')

if __name__=='__main__':main()
