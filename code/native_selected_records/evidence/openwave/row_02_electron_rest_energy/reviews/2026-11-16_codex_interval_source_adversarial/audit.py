"""Independent exact source binding and high-precision finite224 certificate audit.

Imports no producer module. Point containment checks supplement, rather than
replace, the reviewed outward-arithmetic and exact invariant-sector proofs.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
from mpmath import mp

mp.dps=130
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'work/2026-11-16_codex_interval_source_certificate'
COUNTS={'independent_formula_evaluations':0,'checks':0,'point_enclosures':0}
TRI=((0,1),(1,2),(0,2))
TICKETS=((0,1),(0,1),(0,2),(0,2),(1,2))


def require(ok, message):
    COUNTS['checks']+=1
    if not ok:raise AssertionError(message)


def event(position, edge, coin, law='native'):
    COUNTS['independent_formula_evaluations']+=1
    q=[int(i==position) for i in range(3)];a,b=edge
    if law=='native':
        total=q[a]+q[b];q[a]=(total+int(coin==1))//2;q[b]=total-q[a]
    else:q[a],q[b]=q[b],q[a]
    return q.index(1)


def mean(edges,law='native'):
    p=[[F(0) for _ in range(3)] for _ in range(3)]
    for e in edges:
        for c in (0,1):
            for i in range(3):p[event(i,e,c,law)][i]+=F(1,2*len(edges))
    return p


def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(a,v):return [dot(row,v) for row in a]
def exact_matrix(a):return mp.matrix([[mp.mpf(x.numerator)/x.denominator for x in row] for row in a])


def exact_source():
    p=mean(TICKETS);vectors=((0,1,-1),(-2,1,1));mu=[]
    require(all(sum(row)==1 for row in p),'row stochastic')
    require(all(sum(p[i][j] for i in range(3))==1 for j in range(3)),'column stochastic')
    for v in vectors:
        y=mv(p,v);i=next(i for i,x in enumerate(v) if x);value=y[i]/v[i]
        require(y==[value*x for x in v],'exact native optical eigenvector')
        mu.append(value)
    require([10*(1-x) for x in mu]==[4,6],'filed spectral constants follow from native mean')
    require(dot(vectors[0],vectors[1])==0 and all(sum(v)==0 for v in vectors),'orthogonal optical sector')
    permutation=[]
    for position in range(3):
        for choice in range(11):
            if choice==10:new,other=position,10
            else:
                edge=TICKETS[choice//2];coin=choice%2;new=event(position,edge,coin)
                current=int(position==edge[0])-int(new==edge[0])
                other=choice^int(current!=0)
            permutation.append(11*new+other)
    require(sorted(permutation)==list(range(33)),'native R permutation')
    require(all(permutation[permutation[i]]==i for i in range(33)),'native R involution')
    blocks=[]
    for k in range(7):
        eps2=F(1,4**k);transfer=[[F(0) for _ in range(3)] for _ in range(3)]
        weights=[eps2/10]*10+[1-eps2]
        require(sum(weights)==1,'native choice isometry normalization')
        for i in range(3):
            for choice in range(11):
                j,other=divmod(permutation[11*i+choice],11)
                require(weights[other]==weights[choice],'native reversal preserves choice amplitude')
                transfer[j][i]+=weights[choice]
        require(transfer==[[eps2*p[i][j]+(1-eps2)*int(i==j) for j in range(3)] for i in range(3)],'exact A*RA')
        for v,value in zip(vectors,mu):
            frequency_cos=1-eps2*(1-value)
            require(mv(transfer,v)==[frequency_cos*x for x in v],'all-k source eigenblock')
        blocks.append({'k':k,'A_star_R_A':[[str(x) for x in row] for row in transfer]})
    # Every actual rotation mean preserves the full optical plane; it does not
    # require replacement by the compact ray between instruments.
    slots=0
    for m in (4,8):
        for law in ('native','swap'):
            for slot in range(3*m):
                v,r=divmod(slot,m);w=(v+1)%3
                edges=[e for e in TRI if v in e]*(m-r)+[e for e in TRI if w in e]*r
                q=mean(edges,law)
                require(mv(q,(1,1,1))==[1,1,1],'slot fixes uniform mode')
                require(all(sum(q[i][j] for i in range(3))==1 for j in range(3)),'slot preserves optical plane')
                slots+=1
    return p,mu,{'native_P':[[str(x) for x in row] for row in p],
                'optical_eigenvalues':[str(x) for x in mu],'native_R':permutation,
                'all_k_transfer_blocks':blocks,'exact_slot_checks':slots,
                'walk_block_proof':'a=A u,b=(Ra-pa)/sqrt(1-p²); R²=I and A*RA u=pu give A*b=0, Ra=pa+tb,Rb=ta-pb; W=R(2AA*-I) is [[p,-t],[t,p]], and decoder A* maps a to u,b to0.'}


E=mp.matrix([[0,-2/mp.sqrt(6)],[1/mp.sqrt(2),1/mp.sqrt(6)],[-1/mp.sqrt(2),1/mp.sqrt(6)]])
B=mp.zeros(4)
for i in range(2):
    B[i,2*i]=B[i,2*i+1]=1/mp.sqrt(2)
    B[2+i,2*i]=-1j/mp.sqrt(2);B[2+i,2*i+1]=1j/mp.sqrt(2)


def lift(x):
    result=mp.zeros(4)
    for i in range(2):
        for j in range(2):result[i,j]=x[i,j]
    return result


def finite_counter(k,mu):
    q=256*2**k;z=[]
    for value in mu:
        # Exact rational conversion, never binary-float coefficients.
        p=1-(1-mp.mpf(value.numerator)/value.denominator)/(4**k)
        z.extend([p+1j*mp.sqrt(1-p*p),p-1j*mp.sqrt(1-p*p)])
    numerator=[(1-x**q)/q for x in z];root=mp.exp(-2j*mp.pi/q);phase=mp.mpc(1);gram=mp.zeros(4)
    for r in range(1,q//2):
        phase*=root;f=[a/(1-b*phase) for a,b in zip(numerator,z)]
        for i in range(4):
            for j in range(4):gram[i,j]+=f[i]*f[j].conjugate()
    return gram,z


def cp(gram,rho):
    x=B.H*rho*B
    for i in range(4):
        for j in range(4):x[i,j]*=gram[i,j]
    return B*x*B.H


def tr(x):return sum(x[i,i] for i in range(x.rows))


@lru_cache(None)
def operator(m,law,slot):
    v,r=divmod(slot%(3*m),m);w=(v+1)%3
    edges=[e for e in TRI if v in e]*(m-r)+[e for e in TRI if w in e]*r
    p=exact_matrix(mean(edges,law));g=mp.matrix([[E[1,j] for j in range(2)],[E[2,j] for j in range(2)],[E[0,j] for j in range(2)]])
    return (g.H*p*g)**4/9


def dyadic(value):
    sign,mantissa,exponent=value
    return (-1 if sign else 1)*mp.mpf(int(mantissa))*mp.power(2,exponent)


def contains(record,point,label):
    lo=dyadic(record['lo_dyadic']);hi=dyadic(record['hi_dyadic'])
    require(lo<=hi,'ordered exact endpoints '+label)
    require(lo<=point<=hi,'independent point outside certificate '+label)
    COUNTS['point_enclosures']+=1


def compatible(first,later,eta=mp.mpf('.01'),delta=mp.mpf('.00002'),stride=8):
    def principal(p):
        r=2*p-1;a=max(-1,(r-2*delta-eta)/(1-eta));b=min(1,(r+2*delta+eta)/(1-eta))
        require(a<=b,'independent principal domain');return mp.acos(b),mp.acos(a)
    lo,hi=principal(first);a,b=principal(later);components=[]
    for j in range(-stride,stride+1):
        for sign in (1,-1):
            x,y=sorted(((2*mp.pi*j+sign*a)/stride,(2*mp.pi*j+sign*b)/stride))
            x=max(x,lo);y=min(y,hi)
            if x<=y:components.append((x,y))
    merged=[]
    for x,y in sorted(components):
        if merged and x<=merged[-1][1]:merged[-1]=(merged[-1][0],max(y,merged[-1][1]))
        else:merged.append((x,y))
    return merged


def point_case(saved,gram,z,position):
    k,m,turns,direction=(saved[n] for n in ('k','M','turns','direction'))
    reference,law=saved['reference'],saved['law'];frames=3*m*turns;ar=mp.eye(2);af=mp.eye(2)
    for step in range(1,frames+1):
        rot=(2*m+direction*step)%(3*m);r=step%m
        ref=(2*m+direction*min(r,m-r))%(3*m) if reference=='matched' else 2*m
        ar=operator(m,law,rot)*ar;af=operator(m,law,ref)*af
    arms=[ar,af];blocks=[[cp(gram,lift(arms[a]*position*arms[b].H/2)) for b in range(2)] for a in range(2)]
    total=tr(blocks[0][0]+blocks[1][1]).real;visibility=2*tr(blocks[0][1]).real/total
    sign=-1 if law=='native' and turns==1 else 1
    require(sign*visibility>0,'native/control sign')
    tau=(blocks[0][0]+blocks[1][1]+sign*(blocks[0][1]+blocks[1][0]))/2;success=tr(tau).real
    require(success>0,'conditioning probability');spectral=B.H*tau*B/success;weights=[spectral[i,i].real for i in range(4)]
    require(all(x>=0 for x in weights),'positive exact-mode weights')
    single=weights[2]+weights[3];pair=1-weights[0]**2-weights[1]**2
    require(single<mp.mpf('.002') and pair<mp.mpf('.01'),'source impurity hypotheses')
    plane=mp.mpf(1)/(48 if law=='native' else 12)
    cost=success*plane**frames*(1-mp.mpf(1)/256)**(frames*(2 if reference=='matched' else 1))
    cost*=(1-(mp.mpf(3)/8)**8)**11*(1-(mp.mpf(1)/4)**8)**frames
    results={'visibility':visibility,'single_off_band':single,'pair_off_band':pair,
             'selected_clock_probability':cost,'independent_pair_probability':cost**2,
             'log10_selected_clock_probability':mp.log10(cost)}
    for name,value in results.items():contains(saved[name],value,name)
    records={}
    for stride in (1,8):
        zpair=sum(w*x**(stride*2**k) for w,x in zip(weights,z))**2
        values={'x_plus':(1+zpair.real)/2,'x_minus':(1-zpair.real)/2,'y_plus':(1+zpair.imag)/2,'y_minus':(1-zpair.imag)/2}
        for name,value in values.items():contains(saved['pair_records'][str(stride)][name],value,'binary '+name)
        records[stride]=values
    intervals=compatible(records[1]['x_plus'],records[8]['x_plus'])
    require(len(intervals)==1,'all-branch uniqueness')
    for endpoint in intervals[0]:contains(saved['pair_energy_interval_keV'],endpoint*511*mp.sqrt(5)/2,'energy branch endpoint')
    return {'k':k,'M':m,'law':law,'reference':reference,'turns':turns,'direction':direction,
            'single_off_band':mp.nstr(single,30),'pair_off_band':mp.nstr(pair,30)}


def run(ks):
    start=time.monotonic();p,mu,binding=exact_source();matrix=exact_matrix(p)
    certificate=json.loads((SOURCE/'CERTIFICATE.json').read_text());rows=[]
    for k in ks:
        gram,z=finite_counter(k,mu);v=(mp.eye(3)-matrix)*mp.matrix([0,1,0])/(2*mp.sqrt(3))
        v=(matrix/mp.sqrt(3))**8*v;rho=cp(gram,lift(E.H*(v*v.H)*E));position=rho[:2,:2]
        for case in certificate['cases']:
            if case['k']==k:rows.append(point_case(case,gram,z,position))
        print(json.dumps({'k':k,'checked_cases':len(rows)}),flush=True)
    return {'passed':True,'same_group':True,'external_G6':False,'producer_imports':False,
            'precision_digits':mp.dps,'source_certificate_sha256':hashlib.sha256((SOURCE/'CERTIFICATE.json').read_bytes()).hexdigest(),
            'exact_binding':binding,'cases':rows,'counts':COUNTS,'seconds':time.monotonic()-start,
            'scope':'Exact source identities plus independent high-precision point containment; outward theorem reviewed separately. Conditional on filed classical prefixes.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--k',type=int,nargs='+',default=list(range(7)))
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--record-to',type=Path);mode.add_argument('--verify',type=Path,default=HERE/'AUDIT.json')
    args=parser.parse_args()
    require(all(0<=k<=6 for k in args.k),'audit k domain')
    require(args.record_to is None or not args.record_to.exists(),'refusing overwrite')
    result=run(args.k)
    if args.record_to:
        with args.record_to.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:
        expected=json.loads(args.verify.read_text());a=dict(result);b=dict(expected);a.pop('seconds');b.pop('seconds')
        if a!=b:raise AssertionError('saved audit differs')
    print(json.dumps({k:v for k,v in result.items() if k not in ('exact_binding','cases')}))
