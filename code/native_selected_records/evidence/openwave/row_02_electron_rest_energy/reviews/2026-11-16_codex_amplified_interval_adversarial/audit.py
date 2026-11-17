"""Independent amplified-source point audit; no producer imports."""
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
SOURCE=HERE.parents[1]/'work/2026-11-16_codex_amplified_source_interval_v1'
COUNTS={'checks':0,'point_enclosures':0,'independent_formula_evaluations':0}
TRI=((0,1),(1,2),(0,2))
TICKETS=((0,1),(0,1),(0,2),(0,2),(1,2))
def require(ok,message):
    COUNTS['checks']+=1
    if not ok:raise ValueError(message)
def rational(x):
    q=F(x);return mp.mpf(q.numerator)/q.denominator
def exact_matrix(rows):return mp.matrix([[rational(x) for x in row] for row in rows])
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
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


def literal_table(sign,law):
    rows={}
    for p in range(3):
        for edge in TRI:
            for coin in (0,1):
                q=[sign*int(i==p) for i in range(3)];a,b=edge
                total=q[a]+q[b]
                newa=q[b] if law=='swap' else (total+coin)//2
                current=q[a]-newa;q[a]=newa;q[b]=total-newa
                require(q.count(sign)==1 and q.count(0)==2,'literal signed unit')
                rows[p,edge,coin]=(q.index(sign),current)
    return rows


def mean(table,tickets):
    p=[[F(0) for _ in range(3)] for _ in range(3)]
    for edge in tickets:
        for c in (0,1):
            for i in range(3):p[table[i,edge,c][0]][i]+=F(1,2*len(tickets))
    return p


def source_binding(data):
    native_means=[];all_R=[];phase_count=0;frame_count=0
    for sk,source in data['sources'].items():
        sign=int(sk);native=literal_table(sign,'native');control=literal_table(sign,'swap')
        supplied={(r['input_position'],tuple(r['edge']),r['coin']):(r['decoded_output'],r['current']) for r in source['native_columns']}
        require(len(source['native_columns'])==18 and supplied==native,'actual measured native columns')
        p=mean(native,TICKETS);native_means.append(p);R=[]
        for pos in range(3):
            for h in range(11):
                if h==10:R.append(11*pos+h)
                else:
                    y,j=native[pos,TICKETS[h//2],h%2];R.append(11*y+(h^int(j!=0)))
        require(sorted(R)==list(range(33)) and all(R[R[i]]==i for i in range(33)),'signed reversal permutation')
        all_R.append(R)
        for k in range(7):
            weights=[F(1,10*4**k)]*10+[1-F(1,4**k)];P=[[F(0) for _ in range(3)] for _ in range(3)]
            for pos in range(3):
                for h in range(11):
                    y,hh=divmod(R[11*pos+h],11);require(weights[h]==weights[hh],'source choice amplitude preservation');P[y][pos]+=weights[h]
            for v,lam in [((0,1,-1),4),((-2,1,1),6)]:
                require([sum(P[i][j]*v[j] for j in range(3)) for i in range(3)]==[(1-F(lam,10*4**k))*x for x in v],'exact source native clock eigenvalues')
        prep=source['preparation'];selector=prep['cool']['selector_phase'];z=F(selector['real']);s=F(selector['imag'])
        require(z*z+s*s==1,'selector rational unit phase')
        a=5*((5*(z*z-s*s)+6*z-3)**2+(10*z*s+6*s)**2)/512
        inactive=3*((5*(z*z-s*s)-2*z+5)**2+(10*z*s-2*s)**2)/512
        require(a+inactive==1 and 0<a<=1,'complete finite five-choice weights')
        expected=[[a*p[i][j]+inactive*int(i==j) for j in range(3)] for i in range(3)]
        require([[F(x) for x in r] for r in prep['cool']['counts']]==expected,'actual finite cooling mean')
        require([[F(x) for x in r] for r in prep['contrast']['counts']]==[[(F(i==j)-expected[i][j])/2 for j in range(3)] for i in range(3)],'actual contrast mean')
        records=list(prep.values())
        for law,catalog in source['catalogs'].items():
            require(len(catalog)==96 and [r['slot'] for r in catalog]==list(range(96)),'frame catalog coverage')
            table=native if law=='native' else control
            for slot,record in enumerate(catalog):
                vertex,t=divmod(slot,32);other=(vertex+1)%3;mapping=(1,2,0)
                spin=[e for e in TRI if vertex in e]*(32-t)+[e for e in TRI if other in e]*t
                tickets=[tuple(sorted((mapping[a],mapping[b]))) for a,b in spin]
                expected=mean(table,tickets);require([[F(x) for x in r] for r in record['counts']]==expected,'frame mean from native/swap source')
                require(record['coordinate_system']=='energy_target_order','actual common coordinate chart')
                require(all(sum(row)==1 for row in expected) and all(sum(expected[i][j] for i in range(3))==1 for j in range(3)),'entire optical plane invariant')
                frame_count+=1
            records+=catalog
        for record in records:
            re,im=F(record['phase_record']['real']),F(record['phase_record']['imag'])
            require(re*re+im*im==1,'actual matching phase unitarity')
            # |(2z−1)t+(z−1)^2 t^3|² =1−(1−t²)(1−alpha*t²)².
            alpha=2*(1-re);b_re=(re-1)**2-im*im;b_im=2*(re-1)*im
            require((2*re-1)**2+(2*im)**2==1+2*alpha,'forward polynomial linear norm')
            require(2*((2*re-1)*b_re+2*im*b_im)==-(2*alpha+alpha*alpha),'forward polynomial cross norm')
            require(b_re*b_re+b_im*b_im==alpha*alpha,'forward polynomial cubic norm')
            phase_count+=1
    require(native_means[0]==native_means[1],'signed source means agree')
    gauge=[11*p+(h^1 if h<10 else h) for p in range(3) for h in range(11)]
    require(all_R[0]==[gauge[all_R[1][gauge[i]]] for i in range(33)],'negative arc is actual coin-gauge conjugate')
    return {'signed_native_columns':36,'signed_arc_columns':66,'exact_frame_means':frame_count,'exact_phase_polynomials':phase_count,
            'source_clock_lambdas':[4,6],'native_mean':[[str(v) for v in r] for r in native_means[0]],
            'admission_attempt_bound':38360,'calibration_attempt_bound':37568,'entry_attempt_bound_each':75928,
            'entry_failure_upper_each':str(F(1219,2**64))}


@lru_cache(None)
def compiled(record_json):
    record=json.loads(record_json);P=exact_matrix(record['counts']);T=E.H*P*E
    if record['preboost']:T=4*T**3-3*T
    z=mp.mpc(rational(record['phase_record']['real']),rational(record['phase_record']['imag']))
    amplitude=(2*z-1)*T+(z-1)**2*T**3
    return amplitude.H*amplitude


def prepared(source,gram):
    p=source['preparation'];v=compiled(canon(p['cool']))**24*compiled(canon(p['contrast']))*E.H*mp.matrix([0,1,0])
    rho=v*v.H;prob=tr(rho).real;projected=cp(gram,lift(rho))
    return projected[:2,:2],prob,[rho[i,i].real/prob for i in range(2)]


@lru_cache(None)
def arm_pair(catalog_json,h,turns,direction,stride):
    catalog=json.loads(catalog_json);moving=mp.eye(2);fixed=mp.eye(2)
    for step in range(1,96*turns+1):
        if step%stride==0:
            moving=compiled(canon(catalog[(64+direction*step)%96]))**h*moving
            fixed=compiled(canon(catalog[64]))**h*fixed
    return moving,fixed


def point_case(saved,source,position,gram,z):
    k,law,h,turns,direction,stride=(saved[key] for key in ('k','law','h','turns','direction','filter_stride'))
    a,b=arm_pair(canon(source['catalogs'][law]),h,turns,direction,stride)
    aa,bb=cp(gram,lift(a*position*a.H)),cp(gram,lift(b*position*b.H));cross=cp(gram,lift(a*position*b.H))
    analyzer=tr(aa+bb).real/2;visibility=tr(cross).real/analyzer
    parity=-1 if turns%2 else 1
    selected=(a+parity*b)/2;tau=cp(gram,lift(selected*position*selected.H));prob=tr(tau).real
    require(prob>0,'actual fixed-port probability retained')
    spec=B.H*tau*B/prob;weights=[spec[i,i].real for i in range(4)]
    require(all(w>0 for w in weights) and abs(sum(weights)-1)<mp.mpf('1e-115'),'four signed weights')
    cost=prob*(1-mp.mpf(1)/2**16)**(96*turns)*(1-(mp.mpf(3)/8)**8)**2
    unconditional=cost*(1-mp.mpf(1219)/2**64);joint=[[x*y for y in weights] for x in weights]
    values={'visibility':visibility,'selected_clock_probability_given_classical_entry':cost,
        'selected_clock_probability_unconditional_lower':unconditional,'selected_clock_probability_unconditional_upper':cost,
        'independent_pair_probability_unconditional_lower':unconditional**2,'independent_pair_probability_unconditional_upper':cost**2,
        'analyzer_selected_port_probability':prob/analyzer,'single_off_band':weights[2]+weights[3],
        'pair_off_band':1-joint[0][0]-joint[1][1]}
    for name,value in values.items():contains(saved[name],value,name)
    require(saved['selected_path']==('minus' if parity==-1 else 'plus'),'same fixed parity on all laws')
    for i in range(4):
        contains(saved['signed_mode_weights'][i],weights[i],'signed weight')
        for j in range(4):contains(saved['joint_signed_mode_weights'][i][j],joint[i][j],'joint signed weight')
    records={}
    for step in (1,8):
        moment=sum(w*x**(step*2**k) for w,x in zip(weights,z))**2
        rec={'x_plus':(1+moment.real)/2,'x_minus':(1-moment.real)/2,'y_plus':(1+moment.imag)/2,'y_minus':(1-moment.imag)/2}
        for name,value in rec.items():contains(saved['pair_records'][str(step)][name],value,'actual pair binary')
        records[step]=rec
    if values['pair_off_band']<mp.mpf('.01'):
        intervals=compatible(records[1]['x_plus'],records[8]['x_plus']);require(len(intervals)==1,'independent all-branch reader uniqueness')
        for endpoint in intervals[0]:
            contains(saved['pair_frequency_interval'],endpoint,'reader endpoint')
            contains(saved['pair_energy_interval_keV'],endpoint*511*mp.sqrt(5)/2,'anchored energy endpoint')
    return {key:saved[key] for key in ('k','law','h','turns','direction','filter_stride','charge_sign')}|{'selected_cost':mp.nstr(cost,35),'visibility':mp.nstr(visibility,35),'four_signed_weights':[mp.nstr(w,35) for w in weights]}


def run(path,ks):
    started=time.monotonic();data=json.loads((SOURCE/'SOURCE_INPUTS.json').read_text());source=source_binding(data);certificate=json.loads(path.read_text())
    rows=[];preps=[]
    for k in ks:
        gram,z=finite_counter(k,[F(3,5),F(2,5)])
        for sk,program in data['sources'].items():
            position,prob,weights=prepared(program,gram)
            saved=next(r for r in certificate['preparations'] if (r['k'],r['charge_sign'])==(k,int(sk)))
            contains(saved['success'],prob,'preparation success')
            contains(saved['counter_decode_success_given_preparation'],tr(position).real/prob,'counter decode success')
            for i in range(2):contains(saved['optical_weights'][i],weights[i],'preparation optical weight')
            preps.append({'k':k,'sign':int(sk),'success':mp.nstr(prob,35)})
            for case in certificate['cases']:
                if (case['k'],case['charge_sign'])==(k,int(sk)):rows.append(point_case(case,program,position,gram,z))
        print(json.dumps({'k':k,'checked_cases':len(rows)}),flush=True)
    require(len(rows)==sum(48 if k in (0,6) else 32 for k in ks),'independent exact case multiplicity')
    return {'passed':True,'producer_imports':False,'same_group':True,'external_G6':False,'precision_digits':mp.dps,
            'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'SOURCE_INPUTS_sha256':hashlib.sha256((SOURCE/'SOURCE_INPUTS.json').read_bytes()).hexdigest(),
            'source_binding':source,'preparations':preps,'cases':rows,'counts':dict(COUNTS),'seconds':time.monotonic()-started,
            'scope':'Independent point containment and exact native/operator binding; outward and conditioning proof reviewed separately. Declared coherent services remain hypotheses.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,default=SOURCE/'CERTIFICATE.json');p.add_argument('--k',nargs='+',type=int,default=list(range(7)));p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None and args.record_to.exists():raise FileExistsError(args.record_to)
    result=run(args.certificate,args.k)
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'AUDIT.json').read_text());a=dict(result);b=dict(saved);a.pop('seconds');b.pop('seconds');require(a==b,'saved audit result differs')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','preparations','source_binding')},indent=2))


