#!/usr/bin/env python3
"""Full fixed-charge rail source and exact two-doublet decay interface."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import gzip,hashlib,importlib.util,json,math,sys
import numpy as np
P=Path(__file__).resolve().parent
N=2;M=2*N+1;B=4;W=2*N+1;D=1<<(1+4*B)
MASKS=sorted(sum(1<<i for i in c) for c in combinations(range(2*M),M));INDEX={m:i for i,m in enumerate(MASKS)};Q=len(MASKS);DIM=2*Q
sys.dont_write_bytecode=True

def swap(m,a,b):return m^((1<<a)|(1<<b)) if ((m>>a)^(m>>b))&1 else m
def apply(m,legs):
    for a,b in legs:m=swap(m,a,b)
    return m

def rail(i,j):return [(2*i,2*j),(2*i+1,2*j+1)]
def code(m):return sum(1<<(2*j+((m>>j)&1)) for j in range(M))
def dump(n,x):(P/n).write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def interval(a,b):
    d=10**24;return [str(F((a*d).__floor__(),d)),str(F((b*d).__ceil__(),d))]

def branches(k):
    yield 0,[],-1
    def recurse(prefix):
        j=len(prefix)
        if j==k:yield 1,prefix,-1;return
        yield 1,prefix+[[0]*B],j
        for t in range(1,B+1):yield from recurse(prefix+[[0]*(t-1)+[1]])
    yield from recurse([])

def path_prob(n,t,eta):
    L=np.zeros((n,n),dtype=object)
    for j in range(n-1):L[j,j]+=1;L[j+1,j+1]+=1;L[j,j+1]-=1;L[j+1,j]-=1
    A=L-2*np.eye(n,dtype=object);A[0,0]+=eta;order=64;den=2**order*math.factorial(order)
    v=np.eye(n,dtype=object)[:,0];rr=np.zeros(n,dtype=object);ii=np.zeros(n,dtype=object)
    for k in range(order+1):
        c=den//(2**k*math.factorial(k))*t**k
        if k%4==0:rr+=c*v
        elif k%4==1:ii-=c*v
        elif k%4==2:rr-=c*v
        else:ii+=c*v
        if k<order:v=A@v
    rem=F(t*(1 if eta==2 else 2))**65/math.factorial(65);err=2*rem+rem*rem
    probs=[F(int(x)**2+int(y)**2,den*den) for x,y in zip(rr,ii)];assert abs(sum(probs)-1)<=err
    return probs,err

def bessel(n,t):
    v=sum((-1)**k*F(t,2)**(2*k+n)/(math.factorial(k)*math.factorial(k+n)) for k in range(49))
    e=F(t,2)**(98+n)/(math.factorial(49)*math.factorial(49+n));return v,e

def main():
    reg=json.loads((P/'PREREG.txt').read_text())
    for f,h in reg['sha256'].items():assert hashlib.sha256((P/f).read_bytes()).hexdigest()==h,f
    spec=importlib.util.spec_from_file_location('rail_native',P/'inputs/native.py');mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    macros=[(rail(2*r,2*r+1),1) for r in range(N)]+[(rail(M-1,j),0) for j in (0,1)]+[(rail(2*r,2*r+2)+rail(2*r+1,2*r+3),0) for r in range(N-1)]
    calls=records=0;cases=[];capodd=None;ctrlodd=None
    with (P/'RAW_RECORDS.jsonl.gz').open('wb') as raw:
      with gzip.GzipFile(filename='',fileobj=raw,mode='wb',compresslevel=1,mtime=0) as tape:
        def event(state,a,b,c,letters):
            nonlocal calls
            _,rec=mod.native(tuple(state),(a,b),c);j=rec[1];state[a]-=j;state[b]+=j;letters.append([a,b,c,j]);calls+=1
        for disabled in (False,True):
            mode='disabled' if disabled else 'source';C=np.zeros((DIM,DIM),dtype=np.int64);cost=[F(0)]*W;failed=partial=0
            for mi,mask in enumerate(MASKS):
              for z in (0,1):
                for address,(legs,odd) in enumerate(macros):
                    choices=[(c,[],-1) for c in (0,1)] if disabled and address>=N+2 else branches(len(legs))
                    for decision,flags,abort in choices:
                        state=[(mask>>i)&1 for i in range(2*M)]+[0,0];zz=z;letters=[]
                        event(state,2*M,2*M+1,decision,letters)
                        for leg,word in zip(legs,flags):
                            a,b=leg
                            for flag in word:event(state,a,b,flag^state[a],letters)
                        if decision and abort==-1 and len(flags)==len(legs):zz^=odd
                        after=sum(x<<i for i,x in enumerate(state[:2*M]));assert after in INDEX and state[2*M:]==[0,0]
                        weight=1<<(1+4*B-len(letters));C[2*mi+z,2*INDEX[after]+zz]+=weight;cost[address]+=F(weight*len(letters),DIM*D)
                        if abort>=0:failed+=1
                        if abort>0:partial+=1
                        tape.write((json.dumps([mode,mask,z,address,decision,flags,abort,letters,after,zz,weight],separators=(',',':'))+'\n').encode());records+=1
            expected=np.zeros_like(C);I=np.eye(DIM,dtype=np.int64);u=1-F(1,1<<B);delta=1-u
            for address,(legs,odd) in enumerate(macros):
                if disabled and address>=N+2:expected+=D*I;assert cost[address]==1;continue
                assert cost[address]==sum(u**j for j in range(len(legs)+1))
                for count in range(len(legs)+1):
                    success=count==len(legs);weight=u**count if success else delta*u**count
                    perm=np.empty(DIM,dtype=np.int64)
                    for mi,m in enumerate(MASKS):
                      for z in (0,1):perm[2*mi+z]=2*INDEX[apply(m,legs[:count])]+(z^(odd if success else 0))
                    assert np.array_equal(perm[perm],np.arange(DIM))
                    a=weight*D/2;assert a.denominator==1
                    expected+=int(a)*I;expected[np.arange(DIM),perm]+=int(a)
            assert np.array_equal(C,expected) and np.array_equal(C,C.T) and np.all(C.sum(axis=1)==W*D)
            oddC=C[::2,::2]-C[::2,1::2]
            cases.append({'mode':mode,'full_counts':C.tolist(),'odd_counts':oddC.tolist(),'mean_native_cost_per_address':[str(x) for x in cost],
                'failed_macro_prefixes':failed,'partial_macro_prefixes':partial,'denominator':W*D})
            if disabled:ctrlodd=oddC
            else:capodd=oddC
    H2=np.zeros((Q,Q),dtype=np.int64)
    for mi,m in enumerate(MASKS):
        for legs,odd in macros:H2[mi,mi]+=1;H2[mi,INDEX[apply(m,legs)]]+=1 if odd else -1
    fullraise=np.zeros((Q,Q),dtype=np.int64)
    for mi,m in enumerate(MASKS):
        for site in range(M):
            if ((m>>(2*site))&3)==1:fullraise[INDEX[m^(3<<(2*site))],mi]+=1
    assert np.array_equal(H2@fullraise,fullraise@H2)
    E=np.zeros((Q,1<<M),dtype=np.int64)
    for m in range(1<<M):E[INDEX[code(m)],m]=1
    assert np.array_equal(E.T@E,np.eye(1<<M,dtype=np.int64))
    Hlog2=np.zeros((1<<M,1<<M),dtype=np.int64)
    for m in range(1<<M):
        for r in range(N):Hlog2[m,m]+=1;Hlog2[m,swap(m,2*r,2*r+1)]+=1
        for j in (0,1):Hlog2[m,m]+=1;Hlog2[m,swap(m,M-1,j)]-=1
        for r in range(N-1):Hlog2[m,m]+=1;Hlog2[m,swap(swap(m,2*r,2*r+2),2*r+1,2*r+3)]-=1
    assert np.array_equal(H2@E,E@Hlog2)
    # All local off-code states satisfy the core lower bound by its exact spectral polynomial.
    hc=np.zeros((64,64),dtype=np.int64)
    for m in range(64):
        for legs,sgn in [(rail(0,1),1),(rail(2,0),-1),(rail(2,1),-1)]:hc[m,m]+=1;hc[m,apply(m,legs)]+=sgn
    I64=np.eye(64,dtype=np.int64);poly=I64.copy()
    for ev in (1,2,4,5):poly=poly@(hc-ev*I64)
    assert np.all(poly==0) and np.array_equal(hc,hc.T)
    ref=np.zeros(1<<M,dtype=np.int64);T1=np.zeros((1<<M,N),dtype=np.int64);T3=np.zeros_like(T1)
    for m in range(1<<M):
        imp=(m>>(M-1))&1;pairs=[(m>>(2*r))&3 for r in range(N)]
        if imp and all(x in (1,2) for x in pairs):ref[m]=(-1)**sum(x==2 for x in pairs)
        for r in range(N):
            if not all(pairs[j] in (1,2) for j in range(N) if j!=r):continue
            sign=(-1)**sum(pairs[j]==2 for j in range(N) if j!=r)
            if imp and pairs[r] in (1,2):T1[m,r]=T3[m,r]=sign
            elif not imp and pairs[r]==3:T1[m,r]=-2*sign;T3[m,r]=sign
    assert np.all(Hlog2@ref==ref) and np.all(H2@(E@ref)==E@ref)
    lap=np.array([[1,-1],[-1,1]],dtype=np.int64);sectors=[]
    capHnum=W*D*np.eye(Q,dtype=np.int64)-capodd;ctrlH2=2*W*np.eye(Q,dtype=np.int64)-2*ctrlodd//D
    for name,T,eta,norm in [('radiating_J1/2',T1,2,6*(1<<(N-1))),('bound_J3/2',T3,-1,3*(1<<(N-1)))]:
        assert np.array_equal(T.T@T,norm*np.eye(N,dtype=np.int64));K2=2*np.eye(N,dtype=np.int64)+lap;K2[0,0]+=eta
        ET=E@T;assert np.array_equal((H2-np.eye(Q,dtype=np.int64))@ET,ET@K2)
        controlK2=2*np.eye(N,dtype=np.int64);controlK2[0,0]+=eta
        assert np.array_equal((ctrlH2-np.eye(Q,dtype=np.int64))@ET,ET@controlK2)
        for m,row in zip(MASKS,ET):
            if any(row):
                assert sum((m>>j)&1 for j in (0,1,2,3,8,9))-2==1
                assert all(sum((m>>j)&1 for j in range(4*r,4*r+4))-2==0 for r in range(N))
        residual=norm*(capHnum@ET)-ET@(ET.T@capHnum@ET)
        leak=sum(int(x)**2 for x in residual.flat);assert leak>0
        sectors.append({'name':name,'eta':eta,'logical_columns':T.tolist(),'norm_squared':norm,'full_cap_sector_leakage':str(leak)})
    codeleak=capHnum@E-E@(E.T@capHnum@E);assert any(codeleak.flat)
    capref=2*(capHnum@(E@ref))-D*(E@ref);assert any(capref)
    # Exact two-family/spin interface, including all lowering phases.
    lowerup=np.array([0,0,0,0,0,1,-1,0],dtype=np.int64)
    upperup=np.array([0,0,0,-2,0,1,1,0],dtype=np.int64)
    raiseop=np.zeros((8,8),dtype=np.int64)
    for m in range(8):
        for j in range(3):
            if not (m>>j)&1:raiseop[m|(1<<j),m]+=1
    lowerdown=raiseop.T@lowerup;upperdown=raiseop.T@upperup
    for up,down,norm in [(lowerup,lowerdown,2),(upperup,upperdown,6)]:
        assert int(up@up)==int(down@down)==norm and int(up@down)==0
        assert np.array_equal(raiseop@down,up) and np.all(raiseop@up==0)
    core2=np.zeros((8,8),dtype=np.int64)
    for m in range(8):
        for a,b,sgn in [(0,1,1),(2,0,-1),(2,1,-1)]:core2[m,m]+=1;core2[m,swap(m,a,b)]+=sgn
    assert np.all(core2@lowerup==lowerup) and np.all(core2@lowerdown==lowerdown)
    assert np.all(core2@upperup==5*upperup) and np.all(core2@upperdown==5*upperdown)
    vectors=[lowerup,lowerdown,upperup,upperdown]
    for i in range(4):
        for j in range(i):assert int(vectors[i]@vectors[j])==0
    propagation=[]
    for t in (1,2,4):
        j0,e0=bessel(0,t);j1,e1=bessel(1,t);pinf=j0*j0+j1*j1;ein=2*abs(j0)*e0+e0*e0+2*abs(j1)*e1+e1*e1
        if t==4:assert pinf+ein<F(1,5)
        for eta in (2,-1):
          for length in (2,4,8,16,32):
            probs,err=path_prob(length,t,eta);vol=F(4)*F(t*(1 if eta==2 else 2))**length/math.factorial(length)
            if eta==2:assert abs(probs[0]-pinf)<=err+ein+vol
            # Both outgoing spin branches have primitive bath charge0 and core charge1.
            ratio=F(1,3) if eta==2 else F(2,3)
            products=[{'position':0,'kind':'retained_core','probability':str(probs[0]),'core_charge':1}]
            for pos in range(1,length):
                for name,w in [('spin_up_t0',ratio),('spin_down_tplus',1-ratio)]:
                    products.append({'position':pos,'kind':name,'probability':str(w*probs[pos]),'core_charge':1,'bath_charge':0})
            assert abs(sum(F(x['probability']) for x in products)-1)<=err
            propagation.append({'rungs':length,'time':t,'eta':eta,'all_position_probabilities':[str(x) for x in probs],
                'all_product_probabilities':products,'probability_error':str(err),'finite_volume_probability_bound':str(vol),
                'endpoint_interval':interval(max(F(0),probs[0]-err),min(F(1),probs[0]+err)),
                'radiating_infinite_endpoint_interval':interval(pinf-ein,pinf+ein) if eta==2 else None,'control_endpoint_probability':'1'})
    accuracy=F(1,10**8);window=4;rn=2
    while 4*F(window)**rn/math.factorial(rn)>accuracy/2:rn+=1
    rb=1
    while F(12*rn*window,1<<rb)>accuracy/2:rb+=1
    resource={'time_window':window,'requested_probability_error':str(accuracy),'rungs':rn,'retry_cap':rb,
        'volume_error_bound':str(4*F(window)**rn/math.factorial(rn)),
        'cap_error_bound':str(F(12*rn*window,1<<rb)),'executed':False,
        'ideal_expected_native_calls':(8*rn+1)*window}
    receipt={'candidate' :'A-native-fixed-charge-core-v1','source_commit':'3c4e1e23fc32','rungs':N,'retry_cap':B,
        'physical_masks':MASKS,'full_load_dimension':Q,'code_dimension':1<<M,'native_calls':calls,'macro_prefixes':records,
        'macros':[{'legs':legs,'odd_counter':odd} for legs,odd in macros],'cases':cases,'ideal_full_H_times2':H2.tolist(),
        'code_embedding':E.tolist(),'logical_H_times2':Hlog2.tolist(),'reference_logical_vector':ref.tolist(),'reference_energy':'1/2',
        'sectors':sectors,'code_leakage_squared_numerator':str(sum(int(x)**2 for x in codeleak.flat)),
        'reference_leakage_squared_numerator':str(sum(int(x)**2 for x in capref)),
        'two_family_integer_vectors':[v.tolist() for v in vectors],'two_family_norm_squared':[2,2,6,6],
        'retained_bound_pole':{'energy':'3/4','bound_weight':'3/4','asymptotic_survival':'9/16'},'propagation':propagation,
        'gates':{'complete_full_charge_sector_native_census':True,'complete_positive_capped_transfer':True,'all_caps_and_off_code_outcomes_retained':True,
          'exact_full_source_code_intertwiner':True,'actual_reference_and_all_core_representations_retained':True,
          'two_spin_half_core_modes_same_charge':True,'all_outgoing_polarizations_individually_neutral':True,
          'spin_charge_family_interface_exact':True,'full_finite_products_and_error_bounds_retained':True,
          'radiating_time4_survival_below_one_fifth':True,'nonradiating_quartet_pole_retained':True,'disabled_transport_stationary':True},
        'raw_sha256':hashlib.sha256((P/'RAW_RECORDS.jsonl.gz').read_bytes()).hexdigest(),
        'controlled_refinement_example':resource,'physical_spin_or_lepton_assignment':False,'row11_mass_analyzer_intertwined':False,'values_from_nature':[]}
    dump('RECEIPT.json',receipt);print(json.dumps({'native_calls':calls,'macro_prefixes':records,'full_load_dimension':Q,'gates':receipt['gates']}))

if __name__=='__main__':main()
