"""Detector checks only: valid one-unit optical alternatives, not discoveries."""
import json
import math
import numpy as np
import core

rows=[]
for k in (0,2,6):
    c=core.clock(k);A,R,W=c['A'],c['R'],c['W'];eps=c['epsilon'];states=[]
    for vector,mu in ((np.array([0.,1.,-1.])/math.sqrt(2),.4),
                      (np.array([-2.,1.,1.])/math.sqrt(6),.6)):
        lam=1-mu*eps**2;af=A@vector;j=(R@af-lam*af)/math.sqrt(1-lam*lam)
        states.append((af-1j*j)/math.sqrt(2))
    f,g=states;omega=2*math.asin(eps/math.sqrt(5))/eps
    for name,rho,expect_contains,expect_sign in (
        ('positive_f',np.outer(f,f.conj()),True,True),
        ('negative_f',np.outer(f.conj(),f),True,False),
        ('higher_g',np.outer(g,g.conj()),False,True)):
        raw=core.ramsey(rho,W,2**k);read=core.cosine_reader(raw)
        contains=read['lower']<=omega<=read['upper'];sign=read['imaginary_first_moment']>.5
        assert contains==expect_contains and sign==expect_sign
        off=1-float(np.vdot(f,rho@f).real+np.vdot(f.conj(),rho@f.conj()).real)
        rows.append({'k':k,'plant':name,'interval_contains_compact':contains,'positive_sign':sign,
                     'off_band_weight':off,'fixed_certificate_pass':off<=.002+1e-10})
assert all(r['fixed_certificate_pass']==(r['plant']!='higher_g') for r in rows)
print(json.dumps({'type':'sensitivity controls on representable source sector','rows':rows,'passed':True},indent=2))
