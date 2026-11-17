"""Contact-controlled native boundary filter; complete flags and energy bookkeeping."""
from __future__ import annotations
import argparse, gzip, hashlib, importlib.util, itertools, json, math, random, sys, time
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FILTER=HERE.parent/'gravity_native_coherent_filter_research'
STATIC=HERE.parent/'compact_static_adjudication'
REST=HERE.parent/'localized_rest_join'
sys.dont_write_bytecode=True

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
iv=module('boundary_filter_intervals',FILTER/'probe.py')
source=module('boundary_contact_source',STATIC/'encoded_source.py')
native=source.native.native
B=6; ETA=F(1,2); FLAG_BITS=64; EXPONENTS=(24,28,32); MOMENT_ORDER=16

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canonical(obj):return json.dumps(obj,sort_keys=True,indent=2)+'\n'
def check_pins():
    for name in ('SOURCE_PINS.json','PREREG.json','PREREG_AMENDMENT.json'):
        for rel,digest in json.loads((HERE/name).read_text())['files'].items():
            assert sha(ROOT/rel)==digest,rel

def native_boundary_moments(exponent,alpha):
    """All finite paths, occupation plus independent native-coin completion."""
    n=1<<exponent
    scale=iv.inv(iv.mul(iv.rat(18*(1<<(exponent//2))),iv.sqrt_iv(iv.rat(13))))
    noise=iv.inv(iv.sqrt_iv(iv.rat(4*iv.W*n)))
    q=[]
    for t,x in itertools.product(range(3),repeat=2):
        cos=iv.rat(1) if x==0 else iv.rat(-1,2)
        q.append(iv.mul(scale,cos) if t==0 else iv.ZERO)
    D=[]
    for value in q:
        coeff=[iv.ZERO]*(iv.DEG+1)
        for sign in (-1,1):
            value2=iv.add(value,iv.mul(iv.rat(sign),noise))
            term=iv.exp_iv(iv.mul(iv.rat(alpha),value2))
            coeff[0]=iv.add(coeff[0],iv.mul(iv.rat(1,2),term))
            for k in range(1,iv.DEG+1):
                term=iv.mul(term,value2);term=(term[0]//k,iv.ceildiv(term[1],k))
                coeff[k]=iv.add(coeff[k],iv.mul(iv.rat(1,2),term))
        D.append(coeff)
    A=[[[iv.ZERO for _ in range(9)] for _ in range(9)] for _ in range(iv.DEG+1)]
    for i,(t,x) in enumerate(itertools.product(range(3),repeat=2)):
        nb={((t+1)%3)*3+x,((t-1)%3)*3+x,t*3+(x+1)%3,t*3+(x-1)%3}
        for j in range(9):
            p=iv.rat(2*iv.W-4,2*iv.W) if i==j else (iv.rat(1,2*iv.W) if j in nb else iv.ZERO)
            for k in range(iv.DEG+1):A[k][i][j]=iv.mul(p,D[j][k])
    for _ in range(exponent):A=iv.matrix_product(A,A)
    raw=[]
    for k in range(iv.DEG+1):
        total=iv.ZERO
        for row in A[k]:
            for value in row:total=iv.add(total,value)
        raw.append(iv.mul(total,iv.rat(math.factorial(k),9)))
    norm=[iv.div(v,raw[0]) for v in raw]
    return {'Z':raw[0],'mean':norm[1],'second':norm[2],'fourth':norm[4],
            'variance':iv.sub(norm[2],iv.sq(norm[1]))}

def finite_filters():
    rows=[]
    for exponent in EXPONENTS:
        main=native_boundary_moments(exponent,2)
        tail=native_boundary_moments(exponent,14)
        delta=iv.mul(iv.exp_iv(iv.rat(-6*B)),iv.div(tail['Z'],main['Z']))
        base=iv.mul(iv.rat(1,2),iv.mul(iv.exp_iv(iv.rat(-B)),main['Z']))
        p=(iv.mul(base,iv.sub(iv.ONE,(0,delta[1])))[0],base[1])
        coin_tv=iv.div(iv.rat(1,1<<(FLAG_BITS-1)),p)
        qmax=iv.add(iv.div(iv.rat(1<<(exponent//2)),iv.mul(iv.rat(18),iv.sqrt_iv(iv.rat(13)))),
                      iv.sqrt_iv(iv.rat((1<<exponent),4*iv.W)))
        mean_clip=iv.add(iv.mul(delta,iv.sqrt_iv(tail['second'])),
                         iv.div(iv.mul(delta,iv.sqrt_iv(main['second'])),iv.sub(iv.ONE,delta)))
        second_clip=iv.add(iv.mul(delta,iv.sqrt_iv(tail['fourth'])),
                           iv.div(iv.mul(delta,main['second']),iv.sub(iv.ONE,delta)))
        me=iv.add(mean_clip,iv.mul(iv.rat(2),iv.mul(coin_tv,qmax)))
        se=iv.add(second_clip,iv.mul(iv.rat(2),iv.mul(coin_tv,iv.sq(qmax))))
        mag=max(abs(main['mean'][0]),abs(main['mean'][1]))
        ve=iv.add(se,iv.add(iv.mul((2*mag,2*mag),me),iv.sq(me)))
        checks={'mean_limit':iv.inside(main['mean'],F(4,9),F(1,200)),
                'variance_limit':iv.inside(main['variance'],F(2,9),F(1,200)),
                'tail_bound':delta[1]<iv.S//10000,
                'mean_filter_error':me[1]<iv.S//10000,
                'variance_filter_error':ve[1]<iv.S//10000,
                'coin_error':coin_tv[1]<iv.S//1000000}
        rows.append({'n':1<<exponent,'unclipped':{k:iv.show(v) for k,v in main.items()},
                     'tilt14':{k:iv.show(v) for k,v in tail.items()},
                     'clip_tv_bound':iv.show(delta),'original_acceptance_probability':iv.show(p),
                     'original_reject_probability':[str(1-F(iv.show(p)[1])),str(1-F(iv.show(p)[0]))],
                     'mean_filter_error_bound':iv.show(me),'variance_filter_error_bound':iv.show(ve),
                     'coin_tv_bound':iv.show(coin_tv),
                     'expected_native_calls_per_accept_upper':iv.show(iv.div(iv.rat(2*(1<<exponent)+FLAG_BITS+1),iv.sub(p,iv.rat(1,1<<(FLAG_BITS-1))))),
                     'checks':checks})
    reverse=native_boundary_moments(32,-2)
    off=native_boundary_moments(32,0)
    assert iv.inside(reverse['mean'],F(-4,9),F(1,200))
    assert iv.inside(off['mean'],F(0),F(1,10**12))
    assert all(rows[-1]['checks'].values())
    return {'rows':rows,'reverse':{k:iv.show(v) for k,v in reverse.items()},
            'disabled':{k:iv.show(v) for k,v in off.items()},
            'mean_target':'4/9','variance_target':'2/9','finite_source_exactly_gaussian':False}

def contact_energy(graph):
    seams=[tuple(x) for x in graph['seams']];nb=[[]for _ in range(240)]
    for a,b in seams:nb[a].append(b);nb[b].append(a)
    seed=tuple(seams[0]);v={seed:1};moments=[1];sizes=[1]
    for _ in range(MOMENT_ORDER):
        out=defaultdict(int)
        for (p,r),value in v.items():
            out[p,r]+=(len(nb[p])+len(nb[r]))*value
            for q in nb[p]:
                if q!=r:out[q,r]-=value
            for q in nb[r]:
                if q!=p:out[p,q]-=value
        v={k:x for k,x in out.items() if x};moments.append(v.get(seed,0));sizes.append(len(v))
    lower=sum((F(moments[n],n*1260**(n-1)) for n in range(1,MOMENT_ORDER+1)),F(0))
    upper=lower+F(1260)*F(24,1260)**(MOMENT_ORDER+1)/((MOMENT_ORDER+1)*(1-F(24,1260)))
    law=[];m=MOMENT_ORDER
    for k in range(m+1):
        value=F(math.comb(m,k),24**m)*sum((F(math.comb(m-k,r)*(-1)**r*24**(m-k-r)*moments[k+r]) for r in range(m-k+1)),F(0))
        assert value>=0;law.append(value)
    assert sum(law)==1
    return {'preparation':list(seed),'pair_coordinate':256*seed[0]+seed[1],
            'spectral_measure':'<a,b|1_B(-1260 log(I-L_D/1260))|a,b>',
            'LD_moments':moments,'krylov_support_sizes':sizes,
            'mean_energy_interval_exact':[str(lower),str(upper)],
            'mean_energy_decimal_diagnostic':float(lower),
            'LD_spectrum_bound':[0,24],'energy_upper':'-1260 log(1-24/1260)',
            'bernstein_energy_pointer_probabilities':list(map(str,law)),
            'pointer_is_full_spectral_measure':False,'pointer_is_declared_functional_calculus_service':True,
            'post_contact_energy_distribution':{'zero':1},
            'uses_compact_pair_E12':False}

def contact_and_field_tapes(graph):
    sites,edges=iv.graph();contact=[];calls=0
    for e,seam in enumerate(graph['seams']):
        a,b=seam
        for sign in (-1,1):
            label=256*(a if sign==1 else b)+(b if sign==1 else a)
            for coin in (0,1):
                out,j,flag=source.apply_pair_column(None,seam,coin,label);calls+=1
                assert (out,j,flag)==(0,sign,'body')
                contact.append([e,sign,coin,label,out,j])
    trials=[]
    for i in range(20):
        rng=random.Random(17000+i);pos=rng.randrange(iv.M);start=pos;body=[];coins=[];qcount=F(0);noise=0
        for step in range(64):
            t,x,*_=sites[pos]
            if t==0:qcount+=F(1) if x==0 else F(-1,2)
            e=rng.randrange(iv.W);c=rng.randrange(2);a,b=edges[e]
            initial=tuple(int(k==pos)for k in range(iv.M));out,(_,j)=native(initial,(a,b),c);calls+=1
            if j:pos=b if pos==a else a
            assert out[pos]==1;body.append([e,c,j,pos])
            edge=next((u,v) for u,v in edges if pos not in (u,v));c=rng.randrange(2)
            _,(_,j)=native(out,edge,c);calls+=1;assert j==0
            noise+=2*c-1;coins.append([list(edge),c,j])
        q=iv.add(iv.div(iv.rat(qcount.numerator,qcount.denominator),iv.mul(iv.rat(144),iv.sqrt_iv(iv.rat(13)))),
                 iv.div(iv.rat(noise),iv.sqrt_iv(iv.rat(4*iv.W*64))))
        flag_edge=next((a,b)for a,b in edges if pos not in(a,b));draw=0;flags=[]
        for k in range(FLAG_BITS):
            c=0 if i==18 else(1 if i==19 else rng.randrange(2))
            initial=tuple(int(z==pos)for z in range(iv.M));out,(_,j)=native(initial,flag_edge,c);calls+=1
            assert out==initial and j==0;draw=2*draw+c;flags.append([c,j])
        outcomes=[];trigger_records=[]
        ca,cb=graph['seams'][0]
        for trigger in (-1,0,1):
            if trigger:
                label=256*(ca if trigger==1 else cb)+(cb if trigger==1 else ca)
                output,current,flag=source.apply_pair_column(None,(ca,cb),i%2,label);calls+=1
                assert output==0 and current==trigger and flag=='body'
                trigger_records.append({'law':'native','input_pair_label':label,'seam':[ca,cb],
                                        'coin':i%2,'output_pair_label':output,'current':current})
            else:
                label=256*ca+cb;current=0
                trigger_records.append({'law':'changed hold','input_pair_label':label,'seam':[ca,cb],
                                        'coin':i%2,'output_pair_label':label,'current':0})
            assert current==trigger
            exponent=iv.mul(iv.rat(2*trigger),q)
            clipped=(min(exponent[0],B*iv.S),min(exponent[1],B*iv.S))
            prob=iv.mul(iv.rat(1,2),iv.exp_iv(iv.sub(clipped,iv.rat(B))))
            threshold=prob[0]*(1<<FLAG_BITS)//iv.S
            assert prob[1]-prob[0]<iv.S//(1<<FLAG_BITS)
            outcomes.append({'contact_current':trigger,'probability_interval':iv.show(prob),
                             'threshold':threshold,'accept':draw<threshold})
        trials.append({'seed':17000+i,'initial':start,'body':body,'completion_coins':coins,
                       'q_count':str(qcount),'noise_sum':noise,'q_interval':iv.show(q),'flag_edge':list(flag_edge),
                       'flags':flags,'draw':draw,'outcomes':outcomes,'trigger_records':trigger_records,
                       'flag_mode':'forced-zero'if i==18 else('forced-one'if i==19 else'seeded'),
                       'forced_flag_branch_probability':'2^-64'if i>=18 else None})
    assert all(x['accept']for x in trials[18]['outcomes'])
    assert all(not x['accept']for x in trials[19]['outcomes'])
    return {'contact_rows':contact,'trials':trials,'literal_native_calls':calls,
            'contact_branch_count':len(contact),'fresh_trial_contact_calls':40,'changed_hold_trigger_records':20,'seeded_accepts_by_trigger':{str(j):sum(t['outcomes'][j+1]['accept']for t in trials[:18])for j in (-1,0,1)},
            'same_field_tapes_for_trigger_controls':True}

def gaussian_instrument():
    mp.mp.dps=80;s2=mp.mpf(2)/9;Om=mp.sqrt(21)/2;E=mp.acosh(mp.mpf(5)/2)
    bb=Om**2*s2-1/(4*s2);cut=mp.mpf(B)/2
    def rho(q):return mp.exp(-q*q/(2*s2))/mp.sqrt(2*mp.pi*s2)
    def p(q):return mp.exp(min(2*q,B)-B)/2
    def vals(q,accepted):
        prob=p(q);f=mp.sqrt(prob if accepted else 1-prob)
        fp=(f if accepted else -prob/f)if q<cut else mp.mpf(0)
        return f,fp
    reports=[]
    for accepted in (True,False):
        def integ(kind):
            def integrand(q):
                f,fp=vals(q,accepted)
                value={'prob':f*f,'q':q*f*f,'q2':q*q*f*f,
                       'p2':(fp-q*f/(2*s2))**2+bb*f*f}[kind]
                return rho(q)*value
            return mp.quad(integrand,[-mp.inf,cut,mp.inf])
        P,Q,Q2,P2=[integ(k)for k in ('prob','q','q2','p2')]
        reports.append({'outcome':'accept'if accepted else'reject','probability':str(P),
                        'unnormalized_q':str(Q),'unnormalized_q2':str(Q2),'unnormalized_p2':str(P2),
                        'mean_q':str(Q/P),'variance_q':str(Q2/P-(Q/P)**2),
                        'mean_p':'0','symmetrized_qp':'0',
                        'unnormalized_energy':str(E/(2*Om)*(P2+Om**2*Q2))})
    totalP=sum(mp.mpf(x['probability'])for x in reports)
    totalQ=sum(mp.mpf(x['unnormalized_q'])for x in reports)
    totalH=sum(mp.mpf(x['unnormalized_energy'])for x in reports)
    work=E/(2*Om)*mp.quad(lambda q:rho(q)*p(q)/(1-p(q)),[-mp.inf,cut])
    assert abs(totalP-1)<mp.mpf('1e-60') and abs(totalQ)<mp.mpf('1e-60')
    assert abs(totalH-E*Om*s2-work)<mp.mpf('1e-60')
    assert work>0
    return {'energy_E':str(E),'canonical_Omega':str(Om),'thermal_beta':3,
            'q_variance':'2/9','branches':reports,'baseline_energy':str(E*Om*s2),
            'apparatus_average_work':str(work),'total_nonselective_energy':str(totalH),
            'eta_one_rejection_energy':'diverges logarithmically at finite cutoff',
            'eta_half_complete_energy_finite':True,'gaussian_quadratures_are_high_precision_diagnostics':True,
            'native_finite_filter_is_causal_quantum_service':False,
            'joint_matter_field_energy_conservation_established':False,
            'physical_neutrino_or_PMNS_claim':False}

def run(record):
    check_pins();start=time.perf_counter()
    graph=json.loads((REST/'evidence/graph.json').read_text())
    energy=contact_energy(graph);print('Exact contact energy moments and complete coarse pointer passed.',flush=True)
    tapes=contact_and_field_tapes(graph);print('Native contact/body/completion/flag tapes passed.',flush=True)
    finite=finite_filters();print('Finite source interval moments and filter bounds passed.',flush=True)
    quantum=gaussian_instrument();print('Complete Gaussian boundary instrument and work identity passed.',flush=True)
    scientific={'candidate':'A-contact-boundary-field-filter-v1','class':'N+ record protocol with conditional OS instrument',
                'contact_energy':energy,'finite_filter':finite,'gaussian_instrument':quantum,
                'literal_native_calls':tapes['literal_native_calls'],'contact_branches':tapes['contact_branch_count'],
                'seeded_accepts_by_trigger':tapes['seeded_accepts_by_trigger'],
                'eta':'1/2','clip_B':B,'flag_bits':FLAG_BITS,'full_rows_solved':[],
                'nature_values_used':[],'all_registered_component_checks_pass':True}
    files={'RECEIPT.json':canonical(scientific).encode(),
           'NATIVE_TAPES.json.gz':gzip.compress(canonical(tapes).encode(),mtime=0)}
    if record:
        for name,data in files.items():
            if (HERE/name).exists():raise FileExistsError(name)
            (HERE/name).write_bytes(data)
    else:
        assert json.loads((HERE/'RECEIPT.json').read_text())==scientific
        assert json.loads(gzip.decompress((HERE/'NATIVE_TAPES.json.gz').read_bytes()))==tapes
    check_pins();print(canonical({'all_component_checks_pass':True,'seconds':time.perf_counter()-start,
                                'literal_native_calls':tapes['literal_native_calls'],'full_rows_solved':[]}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');run(ap.parse_args().record)
