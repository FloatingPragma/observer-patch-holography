#!/usr/bin/env python3
"""Independent complete finite prefix and fixed-charge interface audit."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import gzip,hashlib,json,math
import numpy as np
P=Path(__file__).resolve().parent

def swap(m,a,b):return m^((1<<a)|(1<<b)) if ((m>>a)^(m>>b))&1 else m

def main():
    r=json.loads((P/'RECEIPT.json').read_text());N=2;M=5;B=4;D=1<<17;W=5
    masks=sorted(sum(1<<i for i in c) for c in combinations(range(10),5));index={m:i for i,m in enumerate(masks)};Q=252;dim=504
    macros=[([(0,2),(1,3)],1),([(4,6),(5,7)],1),([(8,0),(9,1)],0),([(8,2),(9,3)],0),([(0,4),(1,5),(2,6),(3,7)],0)]
    assert masks==r['physical_masks'] and r['code_dimension']==32
    calls=records=0;Cs={};seen={};weights={};costs={};failed={};partial={}
    with gzip.open(P/'RAW_RECORDS.jsonl.gz','rt') as tape:
      for line in tape:
        row=json.loads(line);mode,m,z,address,decision,flags,abort,letters,after,zz,weight=row
        assert mode in ('source','disabled') and m in index and z in (0,1) and address in range(W) and decision in (0,1)
        legs,odd=macros[address];k=0 if mode=='disabled' and address==4 else len(legs)
        if decision==0:assert flags==[] and abort==-1
        elif k==0:assert flags==[] and abort==-1
        elif abort==-1:
            assert len(flags)==k
            assert all(1<=len(f)<=B and f==[0]*(len(f)-1)+[1] for f in flags)
        else:
            assert abort in range(k) and len(flags)==abort+1 and flags[-1]==[0]*B
            assert all(1<=len(f)<=B and f==[0]*(len(f)-1)+[1] for f in flags[:-1])
        key=(mode,m,z,address);signature=(decision,tuple(tuple(f) for f in flags),abort)
        if key not in seen:seen[key]=set();weights[key]=0
        assert signature not in seen[key];seen[key].add(signature);weights[key]+=weight
        state=[(m>>i)&1 for i in range(10)]+[0,0];expected=[]
        def event(a,b,c):
            nonlocal calls
            j=state[a]-c if state[a]!=state[b] else 0;state[a]-=j;state[b]+=j;expected.append([a,b,c,j]);calls+=1
        event(10,11,decision)
        for (a,b),word in zip(legs,flags):
            for flag in word:event(a,b,flag^state[a])
        out=sum(x<<i for i,x in enumerate(state[:10]));ez=z^(odd if decision and k>0 and abort==-1 else 0)
        assert expected==letters and out==after and ez==zz and state[10:]==[0,0]
        assert weight==1<<(17-len(letters))
        if mode not in Cs:Cs[mode]=np.zeros((dim,dim),dtype=np.int64);costs[mode]=[F(0)]*W;failed[mode]=partial[mode]=0
        Cs[mode][2*index[m]+z,2*index[out]+ez]+=weight;costs[mode][address]+=F(weight*len(letters),dim*D)
        failed[mode]+=int(abort>=0);partial[mode]+=int(abort>0);records+=1
    assert len(seen)==2*Q*2*W
    for (mode,m,z,address),entries in seen.items():
        k=0 if mode=='disabled' and address==4 else len(macros[address][0])
        assert len(entries)==1+sum(B**j for j in range(k+1)) and weights[(mode,m,z,address)]==D
    for case in r['cases']:
        C=Cs[case['mode']];assert C.tolist()==case['full_counts'] and np.array_equal(C,C.T) and np.all(C.sum(axis=1)==W*D)
        assert (C[::2,::2]-C[::2,1::2]).tolist()==case['odd_counts']
        assert [str(x) for x in costs[case['mode']]]==case['mean_native_cost_per_address']
        assert failed[case['mode']]==case['failed_macro_prefixes'] and partial[case['mode']]==case['partial_macro_prefixes']
    assert calls==r['native_calls'] and records==r['macro_prefixes']
    H2=np.zeros((Q,Q),dtype=np.int64)
    for i,m in enumerate(masks):
        for legs,odd in macros:
            out=m
            for a,b in legs:out=swap(out,a,b)
            H2[i,i]+=1;H2[i,index[out]]+=1 if odd else -1
    assert H2.tolist()==r['ideal_full_H_times2']
    E=np.zeros((Q,32),dtype=np.int64)
    for m in range(32):E[index[sum(1<<(2*j+((m>>j)&1)) for j in range(5))],m]=1
    assert E.tolist()==r['code_embedding'] and np.array_equal(E.T@E,np.eye(32,dtype=np.int64))
    Hlog=np.array(r['logical_H_times2'],dtype=np.int64);assert np.array_equal(H2@E,E@Hlog)
    ref=np.array(r['reference_logical_vector'],dtype=np.int64);assert np.all(H2@(E@ref)==E@ref)
    for sector in r['sectors']:
        T=np.array(sector['logical_columns'],dtype=np.int64);norm=sector['norm_squared'];eta=sector['eta']
        assert np.array_equal(T.T@T,norm*np.eye(2,dtype=np.int64))
        K=np.array([[3+eta,-1],[-1,3]],dtype=np.int64);assert np.array_equal((H2-np.eye(Q,dtype=np.int64))@(E@T),(E@T)@K)
        for m,row in zip(masks,E@T):
            if any(row):
                assert sum((m>>i)&1 for i in (0,1,2,3,8,9))-2==1
                assert sum((m>>i)&1 for i in (4,5,6,7))-2==0
    vectors=[np.array(x,dtype=np.int64) for x in r['two_family_integer_vectors']]
    raising=np.zeros((8,8),dtype=np.int64);hc=np.zeros((8,8),dtype=np.int64)
    for m in range(8):
        for j in range(3):
            if not (m>>j)&1:raising[m|(1<<j),m]+=1
        for a,b,sign in [(0,1,1),(2,0,-1),(2,1,-1)]:hc[m,m]+=1;hc[m,swap(m,a,b)]+=sign
    for i,(norm,energy) in enumerate(zip([2,2,6,6],[1,1,5,5])):
        assert int(vectors[i]@vectors[i])==norm and np.all(hc@vectors[i]==energy*vectors[i])
        for j in range(i):assert int(vectors[i]@vectors[j])==0
    assert np.array_equal(raising@vectors[1],vectors[0]) and np.array_equal(raising@vectors[3],vectors[2])
    for p in r['propagation']:
        n,t,eta=p['rungs'],p['time'],p['eta'];A=np.zeros((n,n),dtype=object)
        for i in range(n-1):A[i,i]+=1;A[i+1,i+1]+=1;A[i,i+1]-=1;A[i+1,i]-=1
        A-=2*np.eye(n,dtype=object);A[0,0]+=eta;den=2**64*math.factorial(64)
        v=np.eye(n,dtype=object)[:,0];rr=np.zeros(n,dtype=object);ii=np.zeros(n,dtype=object)
        for k in range(65):
            coefficient=den//(2**k*math.factorial(k))*t**k
            if k%4==0:rr+=coefficient*v
            elif k%4==1:ii-=coefficient*v
            elif k%4==2:rr-=coefficient*v
            else:ii+=coefficient*v
            if k<64:v=A@v
        probs=[F(int(x)**2+int(y)**2,den*den) for x,y in zip(rr,ii)]
        assert [str(x) for x in probs]==p['all_position_probabilities']
        rem=F(t*(1 if eta==2 else 2))**65/math.factorial(65);err=2*rem+rem*rem
        assert F(p['probability_error'])==err and abs(sum(probs)-1)<=err
        assert F(p['finite_volume_probability_bound'])==4*F(t*(1 if eta==2 else 2))**n/math.factorial(n)
        ratio=F(1,3) if eta==2 else F(2,3)
        products=[{'position':0,'kind':'retained_core','probability':str(probs[0]),'core_charge':1}]
        for j in range(1,n):
            for name,w in [('spin_up_t0',ratio),('spin_down_tplus',1-ratio)]:products.append({'position':j,'kind':name,'probability':str(w*probs[j]),'core_charge':1,'bath_charge':0})
        assert products==p['all_product_probabilities']
        if eta==2:
            vals=[];errors=[]
            for order in (0,1):
                vals.append(sum((-1)**k*F(t,2)**(2*k+order)/(math.factorial(k)*math.factorial(k+order)) for k in range(49)))
                errors.append(F(t,2)**(98+order)/(math.factorial(49)*math.factorial(49+order)))
            exact=sum(x*x for x in vals);error=sum(2*abs(x)*e+e*e for x,e in zip(vals,errors));lo,hi=map(F,p['radiating_infinite_endpoint_interval'])
            assert lo<=exact-error<=exact+error<=hi
            assert abs(probs[0]-exact)<=err+error+F(p['finite_volume_probability_bound'])
    q=F(1,2);assert 1-q/2==F(3,4) and 2-(q+1/q)/2==F(3,4) and F(3,4)/(1-q*q)==1
    out={'pass':True,'same_group':True,'G6_pass':False,'all_native_calls_accounted':calls,'all_macro_prefixes_verified':records,
         'all_full_charge_source_states_retained':Q,'all_off_code_states_retained':220,'full_source_matrices_reconstructed':True,
         'all_finite_product_probabilities_recomputed':True,'all_outgoing_polarizations_individually_neutral':True,
         'two_family_charge_and_spin_interface_verified':True,'physical_spin_or_lepton_assignment':False,
         'receipt_sha256':hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest()}
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
