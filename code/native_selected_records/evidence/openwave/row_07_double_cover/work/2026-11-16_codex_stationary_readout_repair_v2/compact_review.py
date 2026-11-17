"""Standalone, stdlib-only verification of the frozen v5 k=0 compression.

No files, network, subprocesses or existing package imports. This checks
current-derived means -> amplifier polynomial -> actual finite counter/decode
mixture -> stationary interference. Full source words/complements, measured
calibration custody and rigorous outward bounds remain separate packet checks.
Coherent inverse/phase/choice services, unit entry and spatial embedding are
declared. No arbitrary-axis spin or laboratory noise measurement is claimed.
"""
import cmath
import json
import math
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
from fractions import Fraction as F


def check(value,message):
    if not value:raise RuntimeError(message)


def eye():return [[int(i==j)for j in range(3)]for i in range(3)]
def add(a,b):return [[a[i][j]+b[i][j]for j in range(3)]for i in range(3)]
def scale(a,s):return [[s*x for x in row]for row in a]
def mul(a,b):return [[sum(a[i][k]*b[k][j]for k in range(3))for j in range(3)]for i in range(3)]
def adj(a):return [[complex(a[j][i]).conjugate()for j in range(3)]for i in range(3)]
def trace(a):return sum(a[i][i]for i in range(3))
def mv(a,v):return [sum(a[i][j]*v[j]for j in range(3))for i in range(3)]
def outer(v):return [[x*complex(y).conjugate()for y in v]for x in v]
def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)


def mean(edges,law,sign=1):
    p=[[F(0)for _ in range(3)]for _ in range(3)]
    for x in range(3):
        for a,b in edges:
            for coin in(0,1):
                q=[sign*int(j==x)for j in range(3)];difference=q[a]-q[b]
                current=(difference-(2*coin-1)*(abs(difference)%2))//2 if law=='native'else difference
                q[a]-=current;q[b]+=current
                p[q.index(sign)][x]+=F(1,2*len(edges))
    return p


def rational_phase(t):
    n=int((t*(1<<40)).to_integral_value(rounding=ROUND_HALF_EVEN));v=F(n,1<<40)
    return (1-v*v)/(1+v*v),2*v/(1+v*v)


def match(p):
    constant=sum(p[0]);t=trace(p)-constant
    d=2*(trace(mul(p,p))-constant**2)-t*t;a=abs(t)
    # Exact boundary/preboost classification precedes the decimal proposal.
    boundary=a<=1 and d==(1-a)**2
    boost=a<1 and d<(1-a)**2
    with localcontext()as ctx:
        ctx.prec=90;s=(dec(a)+dec(d).sqrt())/2
        z=abs(4*s**3-3*s)if boost else s
        phase=(F(-1),F(0))if boundary else rational_phase(1/(4*z*z-1).sqrt())
    return phase,boost


def amplifier(p):
    (re,im),boost=match(p);alpha=2*(1-re)
    matrix=[[float(x)for x in row]for row in p]
    b=add(scale(mul(mul(matrix,matrix),matrix),4),scale(matrix,-3))if boost else matrix
    square=mul(b,b);remainder=add(eye(),scale(square,-float(alpha)))
    k=add(eye(),scale(mul(add(eye(),scale(square,-1)),mul(remainder,remainder)),-1))
    return k,alpha,boost


def finite_weighted(law):
    with localcontext()as ctx:
        ctx.prec=90;re,im=rational_phase((Decimal(2)/3).sqrt())
    zr,zi=re*re-im*im,2*re*im
    active=5*((5*zr+6*re-3)**2+(5*zi+6*im)**2)/512
    base=mean([(0,1),(0,1),(0,2),(0,2),(1,2)],law)
    return add(scale(base,active),scale(eye(),1-active)),active


def preparation(law):
    p,active=finite_weighted(law)
    c=scale(add(eye(),scale(p,-1)),F(1,2))
    kc,_,_=amplifier(c);kp,_,_=amplifier(p)
    v=mv(kc,[0,1,0])
    for _ in range(24):v=mv(kp,v)
    probability=sum(abs(x)**2 for x in v)
    fweight=abs(v[1]-v[2])**2/(2*probability)
    return v,probability,fweight,active


def counter(v):
    # A†W^tA=T_t(P) follows from the source-current reversal and its choice
    # reflection. D_q is the actual successful choice decoder of QPE_q.
    p=mean([(0,1),(0,1),(0,2),(0,2),(1,2)],'native')
    p=[[float(x)for x in row]for row in p];q=256
    polynomials=[eye(),p]
    for _ in range(2,q):polynomials.append(add(scale(mul(p,polynomials[-1]),2),scale(polynomials[-2],-1)))
    density=scale(eye(),0j);metric=scale(eye(),0.);complete=scale(eye(),0.)
    for bin_index in range(q):
        step=cmath.exp(-2j*math.pi*bin_index/q);phase=1+0j
        decoded=scale(eye(),0j);effect=scale(eye(),q)
        for t,poly in enumerate(polynomials):
            decoded=add(decoded,scale(poly,phase/q))
            if t:effect=add(effect,scale(poly,2*(q-t)*phase.real))
            phase*=step
        effect=scale(effect,1/q**2);complete=add(complete,effect)
        if 0<bin_index<q//2:
            density=add(density,outer(mv(decoded,v)));metric=add(metric,effect)
    check(max(abs(complete[i][j]-int(i==j))for i in range(3)for j in range(3))<1e-9,'Counter completeness')
    prefix=trace(density).real
    return scale(density,1/prefix),metric,prefix


def frame(slot,law):
    vertex,step=divmod(slot%96,32);other=(vertex+1)%3
    edges=((0,1),(1,2),(0,2));mapping=(1,2,0)
    tickets=[e for e in edges if vertex in e]*(32-step)+[e for e in edges if other in e]*step
    physical=[tuple(sorted((mapping[a],mapping[b])))for a,b in tickets]
    p=mean(physical,law);k,alpha,boost=amplifier(p)
    if law=='native':return k,None
    u=(trace(mul(p,p))-1)/2
    if boost:u=u*(4*u-3)**2
    return k,1-(1-u)*(1-alpha*u)**2


def reader(binary):
    check(all(x>=0 for x in binary),'Negative path probability')
    return (binary[0]-binary[1])/sum(binary)


def run():
    v,prepared,fidelity,active=preparation('native')
    _,_,changed_fidelity,_=preparation('swap')
    check(fidelity>.999999 and changed_fidelity<.7,'Source mode-selection control')
    rho,metric,prefix=counter(v);rows=[]
    for law in('native','swap'):
        catalog=[frame(slot,law)for slot in range(96)]
        for stride,h in((1,2),(2,4)):
            for turns in(1,2):
                for direction in(-1,1):
                    a=eye();b=eye();sa=F(1);sb=F(1)
                    for step in range(1,96*turns+1):
                        if step%stride:continue
                        slot=(64+direction*step)%96
                        for _ in range(h):
                            a=mul(catalog[slot][0],a);b=mul(catalog[64][0],b)
                            if law=='swap':sa*=catalog[slot][1];sb*=catalog[64][1]
                    if law=='native':
                        aa=trace(mul(metric,mul(mul(a,rho),adj(a)))).real
                        bb=trace(mul(metric,mul(mul(b,rho),adj(b)))).real
                        cross=trace(mul(metric,mul(mul(a,rho),adj(b)))).real
                        path=[(aa+bb+2*cross)/4,(aa+bb-2*cross)/4]
                        probability=sum(path);binary=[x/probability for x in path]
                        exact_minus_zero=None
                    else:
                        plus=(sa+sb)**2/(2*(sa*sa+sb*sb));minus=(sa-sb)**2/(2*(sa*sa+sb*sb))
                        binary=[float(plus),float(minus)];exact_minus_zero=minus==0
                        probability=float((sa*sa+sb*sb)/2)*trace(mul(metric,rho)).real
                    visibility=reader(binary);selected=1 if turns==1 else 0
                    cost=prefix*probability*binary[selected]*(1-2**-16)**(96*turns)*float(1-F(3,8)**8)**2
                    if law=='native':
                        check(visibility<-.98 if turns==1 else visibility>(.96 if stride==1 else .95),'Sign/return')
                        check(cost>.03,'Complete conditional acceptance')
                    else:check(visibility>.999999,'Symmetry-preserving control')
                    rows.append({'law':law,'stride':stride,'turns':turns,'direction':direction,
                        'visibility':visibility,'selected_clock_given_entry_and_calibration':cost,
                        'minus_probability':binary[1],'minus_exact_zero':exact_minus_zero})
    check(mean([(0,1),(0,2),(1,2)],'native',-1)==mean([(0,1),(0,2),(1,2)],'native',1),'Signed mean conjugacy')
    return {'scope':'Independent source/word compression at k=0; full flags and outward certificate are separate.',
        'actual_preparation_success':prepared,'native_selected_weight':fidelity,'swap_selected_weight':changed_fidelity,
        'finite_selector_inactive_weight':str(1-active),'rows':rows,'gates_pass':True}


if __name__=='__main__':print(json.dumps(run(),indent=2,allow_nan=False))
