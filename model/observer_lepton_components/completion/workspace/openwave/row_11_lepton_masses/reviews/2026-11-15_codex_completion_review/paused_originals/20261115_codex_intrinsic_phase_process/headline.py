#!/usr/bin/env python3
"""Exact positive native count/reference process; no physical mass claim."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from decimal import Decimal, localcontext
import hashlib, json, sys, argparse, time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
from model.observer import native
PRIOR=ROOT/'openwave/row_11_lepton_masses/work/20261114_codex_closed_clock/RECEIPT.json'
PIN='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
KS=(1,2,3,6); CODES=((0,0),(0,1),(1,0),(1,1))
STATES=tuple(product(range(4),range(8)))
CALLS=Counter()

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x): return json.dumps(x,sort_keys=True,indent=2)+'\n'
def fs(x): return str(x)
def matdata(a): return [[str(x) for x in row] for row in a]
def multiply(a,b):
    n=len(a); out=[[F(0) for _ in range(n)] for _ in range(n)]
    for i,row in enumerate(a):
        for k,x in enumerate(row):
            if x:
                for j,y in enumerate(b[k]):
                    if y:out[i][j]+=x*y
    return out
def apply(a,v):return [sum(x*y for x,y in zip(row,v)) for row in a]
def rowapply(v,a):return [sum(v[i]*a[i][j] for i in range(len(v))) for j in range(len(v))]
def dot(v,w):return sum(x*y for x,y in zip(v,w))
def parity(mask,bits):return -1 if (mask&bits).bit_count()%2 else 1
def jcode(c):
    r,o=CODES[c];return r*(1-2*o)
def fdec(q):
    with localcontext() as c:
        c.prec=45
        return str(Decimal(q.numerator)/Decimal(q.denominator))
def loginterval(x,n=100):
    lo=sum((x**a/F(a) for a in range(1,n+1)),F(0))
    return lo,lo+x**(n+1)/F(n+1)/(1-x)
def hold(loads,seam,coin):return loads,(coin,0)

def run():
    CALLS.clear();checks={}
    def check(name,cond):
        checks[name]=bool(cond)
        if not cond:raise AssertionError(name)
    sources=json.loads((HERE/'SOURCES.json').read_text())
    check('frozen_sources',all(sha(ROOT/p)==h for p,h in sources.items()))
    check('native_pin',sha(Path(native.__file__))==PIN)
    for line in (HERE/'PREREG.txt').read_text().splitlines():
        if '  ' in line:
            hv,p=line.split('  ',1);assert sha(HERE/p)==hv,p
    old=json.loads(PRIOR.read_text());addresses=tuple(tuple(r['seam']) for r in old['count_source_address_records'][0])
    check('sixty_pinned_addresses',len(addresses)==60 and len(set(addresses))==60 and all(addresses[i^1]==addresses[i][::-1] for i in range(60)))
    double=tuple(2 if i==2 else 0 for i in range(12));zero=(0,)*12
    roles=('phase1_X','phase1_Y','phase2_X','phase2_Y','refresh_before','refresh_after')
    records={};labels={};holdrecords={}
    for role in roles:
        rows=[];coarse=[];hc=[]
        for a,seam in enumerate(addresses):
            pair=[]
            for coin in (0,1):
                out,rec=native.native(double,seam,coin);CALLS['double_native']+=1
                j=rec[1];pair.append((out,j))
                r=abs(j);o=int(j<0) if j else a%2
                assert j in (-1,0,1) and r*(1-2*o)==j
                rows.append(dict(address=a,coin=coin,current=j,record_code=[r,o],final_load=out,weight='1/120'))
                if role.startswith('phase'):
                    hout,hrec=hold(double,seam,coin);CALLS['changed_hold']+=1
                    hc.append(dict(address=a,coin=coin,current=hrec[1],final_load=hout,weight='1/120'))
            assert pair[0]==pair[1]
            coarse.append(dict(address=a,current=pair[0][1],orientation=a%2))
        records[role]=rows;labels[role]=coarse
        if hc:holdrecords[role]=hc
    check('all720_literal_double_columns',CALLS['double_native']==720)
    pair_tables={};count_ps=[]
    for slot in (1,2):
        rows=[]
        for x in labels[f'phase{slot}_X']:
            for y in labels[f'phase{slot}_Y']:
                a=int(x['current']==1);b=int(y['current']!=0)
                rows.append([x['address'],y['address'],a,b,a^b,a*b])
        pp=F(sum(r[4] for r in rows),len(rows));pa=F(sum(r[5] for r in rows),len(rows))
        count_ps.append(pp)
        check(f'phase{slot}_all3600_count_records',len(rows)==3600 and pp==F(2,9) and pa==F(1,72))
        pair_tables[str(slot)]=rows
    p=count_ps[0]
    pi=[]
    for role in ('refresh_before','refresh_after'):
        weights=Counter(tuple(row['record_code']) for row in records[role])
        one=[F(weights[c],120) for c in CODES]
        check('actual_'+role+'_law',one==[F(5,12),F(5,12),F(1,12),F(1,12)])
        pi.append(one)
    pi=pi[0]
    coin_rows=[];coin_words=[]
    for k in KS:
        for role,width in [('reference',k+1),('phase1',k+1),('phase2',k+1),('before',3*k),('after',3*k)]:
            for bit in range(width):
                for coin in (0,1):
                    out,rec=native.native(zero,(0,1),coin);CALLS['zero_coin_native']+=1
                    assert out==zero and rec==(coin,0)
                    coin_rows.append(dict(k=k,role=role,bit=bit,coin=coin,current=0,weight='1/2'))
            coin_words.append(dict(k=k,role=role,width=width,all_zero_weight=str(F(1,2**width)),
                                   every_word_weight=str(F(1,2**width)),word_count=2**width,
                                   decoder='Each original bit word is retained; flag1 exactly on the all-zero word.'))
    check('all240_literal_zero_coin_columns',CALLS['zero_coin_native']==240)
    rho=[pi[c]/8 for c,b in STATES]
    output={}
    for k in KS:
        e=F(1,2**k);delta=e**3
        t4=[[F(i==j)*(1-delta)+delta*pi[j] for j in range(4)] for i in range(4)]
        T=[[t4[c][d]*int(b==z) for d,z in STATES] for c,b in STATES]
        def dkernel(pp):
            D=[]
            for c,b in STATES:
                j=jcode(c);flip=(e/2,e*pp/2 if j<=0 else F(0),e*pp/2 if j==-1 else F(0))
                row=[]
                for d,z in STATES:
                    value=F(c==d)
                    for bit,prob in enumerate(flip):value*=prob if ((b^z)>>bit)&1 else 1-prob
                    row.append(value)
                D.append(row)
            return D
        D=dkernel(p)
        P=multiply(multiply(T,D),T)
        # Independent direct grouped source law: both refresh outcomes and all
        # three independent counter flags. No matrix multiplication on this side.
        direct=[]
        for c,b in STATES:
            row=[]
            for d,z in STATES:
                ans=F(0)
                for mid in range(4):
                    j=jcode(mid);flip=(e/2,e*p/2 if j<=0 else F(0),e*p/2 if j==-1 else F(0))
                    q=t4[c][mid]*t4[mid][d]
                    for bit,prob in enumerate(flip):q*=prob if ((b^z)>>bit)&1 else 1-prob
                    ans+=q
                row.append(ans)
            direct.append(row)
        check(f'k{k}_full32_original_kernel',P==direct and all(sum(r)==1 for r in P) and all(x>0 for r in P for x in r))
        check(f'k{k}_stationary_detailed_balance',rowapply(rho,P)==rho and all(rho[i]*P[i][j]==rho[j]*P[j][i] for i in range(32) for j in range(32)))
        chars={}
        for a in range(8):
            diag=[]
            for c in range(4):
                j=jcode(c);power=int(bool(a&2) and j<=0)+int(bool(a&4) and j==-1)
                diag.append((1-e)**int(bool(a&1))*(1-e*p)**power)
            expected=[[sum(t4[i][m]*diag[m]*t4[m][j] for m in range(4)) for j in range(4)] for i in range(4)]
            measured=[]
            for c in range(4):
                measured.append([sum(P[STATES.index((c,0))][STATES.index((d,b))]*parity(a,b) for b in range(8)) for d in range(4)])
            assert measured==expected
            chars[str(a)]={'T_D_T':matdata(expected),'middle_diagonal':list(map(str,diag))}
        dmin=(1-e)*(1-e*p)**2;m=(1-delta)**2*dmin
        check(f'k{k}_positive_complete_character_blocks',len(chars)==8 and dmin>0 and m>0)
        clock=[F(parity(1,b)) for c,b in STATES]
        check(f'k{k}_same_process_reference_eigenvector',apply(P,clock)==[(1-e)*v for v in clock])
        corr=[];v=clock
        for n in range(9):
            cv=dot([w*x for w,x in zip(rho,clock)],v);assert cv==(1-e)**n
            corr.append(str(cv));v=apply(P,v)
        trace=sum(P[i][i] for i in range(32));return1=trace/32
        moments=[]
        for i in range(32):
            for j in range(32):
                event=P[i][j]*P[j][i]/32
                reconstructed=event/(P[i][j]*return1)
                assert reconstructed==P[j][i]/trace
                moments.append(dict(i=i,j=j,word_coefficient=str(P[i][j]),closed_probability=str(event),moment=str(reconstructed)))
        check(f'k{k}_all1024_history_matrix_units',len(moments)==1024 and 0<return1<1)
        bounds={}
        glo,ghi=loginterval(e);hlo,hhi=loginterval(e*p)
        rlo,rhi=hlo/ghi,hhi/glo
        difference=p*(1-p)*e/(2*(1-e)*(1-p*e))
        total=2*delta/(m*e)+2*difference
        assert 0<rlo<=rhi<p and p-rlo<=difference
        bounds=dict(reference_gap_interval=[str(glo),str(ghi)],phase_gap_interval=[str(hlo),str(hhi)],
                    ratio_interval=[str(rlo),str(rhi)],ratio_decimal=[fdec(rlo),fdec(rhi)],
                    positive_transfer_lower_bound=str(m),ratio_error_bound=str(difference),
                    joint_generator_to_limit_norm_bound=str(total),joint_bound_decimal=fdec(total),
                    bound_status='Analytical proof; logarithm enclosures use100 positive terms plus explicit tail.')
        controls={}
        for name,pp in [('count_HOLD',F(0)),('count_AND',F(1,72)),('count_one_plant',F(1))]:
            dd=dkernel(pp);changed=multiply(multiply(T,dd),T)
            assert all(sum(row)==1 for row in changed)
            if name=='count_HOLD':assert all(changed[i][j]==0 for i,(_,b) in enumerate(STATES) for j,(_,z) in enumerate(STATES) if (b^z)&6)
            # A frozen conditional current-memory field isolates the actual
            # phase-update stage, with the reference setting unchanged.
            inputi=STATES.index((CODES.index((1,1)),0))
            one=sum(dd[inputi][i] for i,(_,b) in enumerate(STATES) if b&2)
            assert one==e*pp/2
            controls[name]={'full_kernel':matdata(changed),'phase1_flip_given_Jminus':str(one),'all_outcomes':True}
        check(f'k{k}_frozen_native_hold_AND_and_plant',F(controls['count_HOLD']['phase1_flip_given_Jminus'])==0 < e*p/2 < F(controls['count_one_plant']['phase1_flip_given_Jminus']))
        # First reference-mark stop. The mark ignores its last parity bit.
        # Before the mark b0 stays fixed; at the mark it is uniformly redrawn.
        smaller=tuple(product(range(4),range(4)))
        A=[]
        for c,b in smaller:
            i=STATES.index((c,2*b))
            A.append([sum(P[i][STATES.index((d,2*z+b0))] for b0 in (0,1)) for d,z in smaller])
        mu=[F(i==0) for i in range(16)];durations=[]
        for n in range(1,9):
            mu=rowapply(mu,A);weight=e*(1-e)**(n-1)
            final=[weight*mu[smaller.index((c,b//2))]/2 for c,b in STATES]
            assert sum(final)==weight
            durations.append({'duration':n,'probability':str(weight),'full32_terminal_distribution':list(map(str,final))})
        tail=(1-e)**8;timeout=[tail*mu[smaller.index((c,b//2))] if b%2==0 else F(0) for c,b in STATES]
        check(f'k{k}_complete_first_reference_mark',sum(F(r['probability']) for r in durations)+sum(timeout)==1)
        lower=delta**2*e**3*p**2/F(1152)
        assert min(x for row in P for x in row)>=lower
        output[str(k)]={'epsilon':str(e),'delta':str(delta),'stationary_law':list(map(str,rho)),
          'refresh4':matdata(t4),'counter32':matdata(D),'kernel32':matdata(P),'characters':chars,
          'reference_correlations':corr,'return_L1':str(return1),'nonreturn_L1':str(1-return1),
          'history_matrix_units':moments,'smallest_entry_lower_bound':str(lower),'bounds':bounds,
          'controls':controls,'first_mark':{'durations':durations,'timeout_weight':str(tail),'timeout32':list(map(str,timeout)),
             'unbounded_duration_law':'epsilon*(1-epsilon)^(n-1), n>=1; every finite n retained',
             'phase_transfer':'epsilon*A*[I-(1-epsilon)*A]^-1; A is the same16-state phase/current transfer.'}}
    return {'schema':'oph.row11.native_count_reference.v1','candidate':'A-native-count-reference-clock','class':'N+',
      'checks':checks,'all_pass':all(checks.values()),'source_bindings':sources,'literal_call_counts':dict(CALLS),
      'actual_native_calls':CALLS['double_native']+CALLS['zero_coin_native'],'source_roles':records,
      'changed_hold_source_roles':holdrecords,'count_pair_tables':pair_tables,'zero_coin_columns':coin_rows,
      'coin_word_contracts':coin_words,'original_count_probability':str(p),'refresh_law':list(map(str,pi)),
      'states':[{'record_code':CODES[c],'parity_bits_integer':b} for c,b in STATES],'refinements':output,
      'source_compression':'720 fresh literal double columns and240 fresh zero-coin columns; complete independent source words are product cylinders. Every original coin/address/current outcome has a disclosed weight and deterministic memory map. Product words are not counted as additional native calls.',
      'physical_row_established':False,'family_coherent_kick_executed':False,'mass_values_used':False,
      'interpretation':'The positive transfer logarithm is reconstructed from native histories. Its conditional phase is an analytic continuation with a shared reference gap. Choosing a one-reference-time family action and realizing its analyzer on this source remain separate obligations.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',type=Path);ap.add_argument('--check',type=Path);args=ap.parse_args()
    start=time.perf_counter();r=run();blob=canon(r)
    if args.write:
        with args.write.open('x') as f:f.write(blob)
    if args.check and args.check.read_text()!=blob:raise SystemExit('byte replay mismatch')
    print(canon({'all_pass':r['all_pass'],'checks':len(r['checks']),'actual_native_calls':r['actual_native_calls'],
      'changed_law_calls':r['literal_call_counts']['changed_hold'],'states':32,
      'native_count':r['original_count_probability'],
      'ratio_intervals':{k:v['bounds']['ratio_decimal'] for k,v in r['refinements'].items()},
      'generator_error_bounds':{k:v['bounds']['joint_bound_decimal'] for k,v in r['refinements'].items()},
      'seconds':time.perf_counter()-start,'physical_row_established':False}))
