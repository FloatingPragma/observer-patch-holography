"""Bounded source-local flow verification; no producer-module imports."""
from fractions import Fraction as F
from functools import lru_cache
import argparse,hashlib,json,math,time,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
U=F(1,2**53)
sys.set_int_max_str_digits(0)

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sqrt_upper(x):
    scale=10**32
    return F(math.isqrt(x.numerator*scale*scale//x.denominator)+1,scale)
@lru_cache(None)
def trig(x,cosine=False):
    x=F(x);assert abs(x)<=4
    term=F(1) if cosine else x;total=term
    # Alternating tail after48 terms decreases strictly for |x|<=4.
    for k in range(1,48):
        term=-term*x*x/F((2*k-1)*(2*k) if cosine else (2*k)*(2*k+1));total+=term
    next_term=abs(term*x*x/F(95*96 if cosine else 96*97))
    return total,next_term
def trig_error(x):
    c,ce=trig(F(x),True);s,se=trig(F(x),False)
    return abs(F(math.cos(x))-c)+ce+abs(F(math.sin(x))-s)+se
def phase(x):
    c,_=trig(F(x),True);s,_=trig(F(x),False)
    return complex(float(c),float(s))
def raw_pulse_bound(raw,kind,angle):
    co,ce=trig(F(angle),True);si,se=trig(F(angle),False)
    center=[[(F(int(i==j)),F())for j in range(6)]for i in range(6)]
    errors=[[(F(),F())for j in range(6)]for i in range(6)]
    def put(i,j,re,im,er=F(),ei=F()):center[i][j]=(re,im);errors[i][j]=(er,ei)
    if kind in ('I','J'):
        for i in range(2 if kind=='I'else 3):put(i,i,co,si,ce,se)
    elif kind=='Z':put(0,0,co,si,ce,se);put(1,1,co,-si,ce,se)
    elif kind=='X':
        put(0,0,co,F(),ce);put(1,1,co,F(),ce)
        put(0,1,F(),si,ei=se);put(1,0,F(),si,ei=se)
    elif kind=='Y':
        put(0,0,co,F(),ce);put(1,1,co,F(),ce)
        put(0,1,si,F(),se);put(1,0,-si,F(),se)
    elif kind=='H':
        scale=2**240;lo=F(math.isqrt(scale*scale//2),scale);hi=lo+F(1,scale)
        root=(lo+hi)/2;re=(hi-lo)/2;v=si*root;ve=abs(si)*re+abs(root)*se+se*re
        put(0,0,co,v,ce,ve);put(1,1,co,-v,ce,ve)
        put(0,1,F(),v,ei=ve);put(1,0,F(),v,ei=ve)
    else:raise AssertionError(kind)
    sq=sum(((abs(F(raw[i,j].real)-center[i][j][0])+errors[i][j][0])**2+
            (abs(F(raw[i,j].imag)-center[i][j][1])+errors[i][j][1])**2
            for i in range(6)for j in range(6)),F())
    return sqrt_upper(sq)

def local_identity():
    states=[(p,n)for p in range(4)for n in range(4)if p!=n]
    twice=np.zeros((12,12),int);vac=[]
    for j,(p,n)in enumerate(states):
        for coin in (0,1):
            q=[0]*4;q[p]=1;q[n]=-1;total=q[0]+q[1]
            old=q[0];q[0]=total//2+coin*(total%2);q[1]=total-q[0]
            if not any(q):vac.append((p,n,coin,old-q[0]))
            else:twice[states.index((q.index(1),q.index(-1))),j]+=1
    H=2*np.eye(12,dtype=int)-twice
    Le=np.zeros((4,4),int);Le[0,0]=Le[1,1]=1;Le[0,1]=Le[1,0]=-1
    full=np.kron(Le,np.eye(4,dtype=int))+np.kron(np.eye(4,dtype=int),Le)
    inds=[4*a+b for a,b in states]
    assert np.array_equal(H,full[np.ix_(inds,inds)])
    assert np.array_equal(H,H.T) and np.array_equal(H@H,2*H)
    reflection=np.eye(12,dtype=int)-H
    assert np.array_equal(reflection@reflection,np.eye(12)) and len(vac)==4
    wedges=[]
    for a in range(4):
        for b in range(a+1,4):
            w=np.zeros(16,int);w[4*a+b]=1;w[4*b+a]=-1;wedges.append(w)
    A=np.column_stack(wedges);out=full@A
    assert not np.any(out[[0,5,10,15]]) and np.array_equal(out.reshape(4,4,6),-out.reshape(4,4,6).transpose(1,0,2))
    assert np.array_equal(H@A[inds],out[inds])
    return states,twice,H

def matching_once(state,edges,h,pair):
    n=state.shape[0];partner=np.arange(n);moving=np.zeros(n,bool)
    for a,b in edges:partner[a]=b;partner[b]=a;moving[a]=moving[b]=True
    z=phase(-2*h);d=np.ones(n,complex);o=np.zeros(n,complex)
    d[moving]=(1+z)/2;o[moving]=(1-z)/2
    if not pair:return d[:,None]*state+o[:,None]*state[partner]
    # Independent four-term tensor contraction, followed by native contact block.
    out=(d[:,None,None]*d[None,:,None]*state+
         o[:,None,None]*d[None,:,None]*state[partner,:,:]+
         d[:,None,None]*o[None,:,None]*state[:,partner,:]+
         o[:,None,None]*o[None,:,None]*state[partner,:,:][:,partner,:])
    for a,b in edges:
        for i,j in ((a,a),(a,b),(b,a),(b,b)):out[i,j]=z*state[i,j]
    return out
def one_step(initial,program,pair):
    out=initial.copy()
    for edges,h in program:out=matching_once(out,edges,h,pair)
    return out
def decode(short,long,epsilon):
    radius=2*math.sqrt(2)*epsilon/(1-epsilon);angle=radius/(1-radius)
    e=-math.atan2(short[1],short[0])/.1;coarse=[e-10*angle,e+10*angle]
    raw=-math.atan2(long[1],long[0]);lo=math.ceil((4*coarse[0]-raw-angle)/(2*math.pi));hi=math.floor((4*coarse[1]-raw+angle)/(2*math.pi))
    assert lo==hi
    value=(raw+2*math.pi*lo)/4
    return [value-angle/4,value+angle/4],lo

def main():
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path)
    parser.add_argument('--receipt',default='attempts/v2/RECEIPT.json');parser.add_argument('--observer')
    parser.add_argument('--record',action='store_true');args=parser.parse_args()
    started=time.perf_counter();source=args.source.resolve();r=read(source/args.receipt)
    prereg=read(source/'PREREG.json');assert sha(source/'PREREG.json')==r['prereg_sha256']
    for n,h in prereg['files'].items():assert sha(source/n)==h,n
    graph=[tuple(e)for e in read(source/'inputs/graph.json')['seams']]
    assert len(graph)==630 and len(set(map(frozenset,graph)))==630
    cert=read(source/'inputs/certificate.json');basis=np.asarray(cert['integer_modes'],dtype=int)
    L=np.zeros((240,240),int)
    for a,b in graph:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    assert max(L.diagonal())==6
    energies=[F(int(v@(L@v)),int(v@v))for v in basis]
    assert all(np.array_equal(L@v,int(e)*v)for v,e in zip(basis,energies))
    assert [str(e)for e in energies]==r['source_census']['single_generators']
    assert not np.any(basis[:2,12:]) and not np.any(basis[2:,:12])
    for i,v in enumerate(basis[:2]):
        for j,w in enumerate(basis[2:]):
            psi=np.outer(v,w);expected=np.outer(L@v,w)+np.outer(v,L@w)
            assert not np.any(np.diag(expected))
            assert np.array_equal(expected,int(energies[i]+energies[j+2])*psi)
            assert r['source_census']['pair_generators'][2*i+j]['generator']==str(energies[i]+energies[j+2])
    states,twice,H=local_identity()
    assert np.array_equal(H,r['local_charged_census']['generator'])
    assert np.array_equal(twice,r['local_charged_census']['twice_transfer'])
    frames=[]
    for first,second in [(basis[0],basis[1]),(basis[2],basis[3])]:
        assert int(first@first)==int(second@second)==4 and int(first@second)==2
        frames.append(np.column_stack((first/2,(2*second-first)/math.sqrt(12))))
    pair=np.stack([np.outer(frames[0][:,a],frames[1][:,b])for a in range(2)for b in range(2)],axis=-1)
    archive=source/Path(args.receipt).parent
    matrices=np.load(archive/'NATIVE_M6.npz');primitive=[]
    for i,key in enumerate(r['matrix_catalog_keys']):
        kind,angle=key.split('/');bound=raw_pulse_bound(matrices[str(i)],kind,float.fromhex(angle))
        assert bound<=F(r['native_catalog']['raw_frobenius_bounds'][key])+F(1,10**28)
        primitive.append(bound)
    pulse=F(r['native_catalog']['complete_per_pulse_allowance'])
    assert 2*max(primitive)+F(1,10**15)<=pulse+F(1,10**28)
    n=r['source_programs']['4']['repetitions'];independent_certs={};duration_errors={}
    # Exact frame division comparison, including the product-frame gain term.
    scale=10**90;root=math.isqrt(12*scale*scale);rl,rh=F(root,scale),F(root+1,scale)
    delta=max(abs(F(1/math.sqrt(12))-1/rl),abs(F(1/math.sqrt(12))-1/rh))
    assert delta<=U
    frame=100*delta+16*U;actual_frame=10*delta+6*U
    assert 2*actual_frame<=frame
    for t,desc in r['source_programs'].items():
        groups=[[tuple(e)for e in g]for g in desc['edge_groups']]
        flat=[e for g in groups for e in g];assert sorted(flat)==sorted(graph)
        assert all(len({v for e in g for v in e})==2*len(g)for g in groups)
        program=[(groups[x['group']],float.fromhex(x['dt_hex']))for x in desc['one_step']]
        assert len(program)==10*len(groups)
        sums={i:sum((F(float.fromhex(x['dt_hex']))for x in desc['one_step']if x['group']==i),F())for i in range(len(groups))}
        assert len(set(sums.values()))==1
        duration_errors[t]=str(abs(n*sums[0]-F(t)))
        independent_certs[t]={}
        for name,initial,is_pair,E in [('single',frames[0],False,energies[0]),('pair',pair,True,energies[0]+energies[2])]:
            actual=one_step(initial,program,is_pair)
            idealphase=phase(-float(E)*float(t)/n)
            residual=float(np.linalg.norm(actual-idealphase*initial))
            directed=r['certificates'][t][name]['directed_rounding']
            tau=max([trig_error(-2*h)for _,h in program]+[trig_error(-float(E)*float(t)/n)])
            assert tau<=F(directed['directed_trig_L1_error_upper'])+F(1,10**70)
            norm=F(2)if is_pair else F(3,2);depth=len(program)
            local=256*U+16*tau
            assert depth*n*local<F(1,2)
            independent_guard=norm*depth*local/(1-depth*local)+2*frame+F(directed['clock_rounding_upper'])+16*U
            residual_rounding=F(4*initial.size,2**53-4*initial.size)*F(str(residual))+F(1,10**30)
            independent_upper=n*(F(str(residual))+independent_guard+residual_rounding)
            claimed=F(str(r['certificates'][t][name]['selected_subspace_operator_error_upper']))
            assert independent_upper<=claimed,(t,name,float(independent_upper),float(claimed))
            independent_certs[t][name]={'one_step_residual':residual,'independent_upper':str(independent_upper),'producer_upper':str(claimed)}
    source_outputs=np.load(archive/'SOURCE_OUTPUTS.npz');decoded={};budget_results={}
    for name,initial in [('single',frames[0]),('pair',pair)]:
        previous=r['complete_budgets'][name];words=r['words'][name]
        x=words['total_pulses']*pulse;err=x/(1-x)
        full=(err*(2+err)+words['total_pulses']*pulse*pulse/(1-x)**2+err*err)/2
        assert full==F(previous['source_word'])
        comparison=max(F(str(r['certificates'][t][name]['selected_subspace_operator_error_upper']))for t in r['certificates'])
        assert comparison==F(previous['selected_source_flow'])
        stored=max(2*F(r['certificates'][t][name]['directed_rounding']['stored_full_propagation_arithmetic_upper'])+F(r['certificates'][t][name]['directed_rounding']['final_inner_product_readout_upper'])for t in r['certificates'])
        assert stored==F(previous['stored_propagation_and_readout'])
        reset=F(35,2**31);pointer=F(200000000001,20000000000000000000000)
        apparatus=words['initial_and_terminal_scratch_bits']*(reset+2*pointer)+2*pointer+4*pulse
        assert apparatus==F(previous['reset_read'])
        inherited=F(cert['common_full_preparation_and_apparatus_allowance'])
        assert inherited==F(previous['inherited'])
        total=inherited+full+comparison+stored+apparatus;total=total/(1-2*total)
        assert total==F(previous['complete'])<F(1,1000)
        records={}
        for t in ('0.1','4'):
            final=source_outputs[name+'_'+t.replace('.','p')]
            z=np.sum(initial.conj()*final)/initial.shape[-1]
            norm=float(np.sum(abs(final)**2)/initial.shape[-1])
            record=r['complete_unfiltered_reference_record_centres'][t][name]
            assert abs(z-complex(*record['fringe_centre']))<1e-12
            assert abs(norm-record['final_norm_centre'])<1e-12
            records[t]=[z.real,z.imag]
        interval,winding=decode(records['0.1'],records['4'],float(total));decoded[name]=interval
        assert np.max(abs(np.asarray(interval)-r['decoder_intervals'][name]['energy_interval']))<1e-12
        assert winding==r['decoder_intervals'][name]['winding_from_short']
        budget_results[name]=str(total)
    ratio=[decoded['pair'][0]/decoded['single'][1],decoded['pair'][1]/decoded['single'][0]]
    assert np.max(abs(np.asarray(ratio)-r['pair_to_single_interval']))<1e-12 and ratio[0]<2<ratio[1]
    assert 511*ratio[1]<1024.45085
    observer_checked=False
    if args.observer:
        observed=read(source/args.observer);assert observed['flow_receipt_sha256']==sha(source/args.receipt)
        assert observed['all_gates_pass']
        observer_checked=True
    result=dict(producer_receipt_sha256=sha(source/args.receipt),producer_prereg_sha256=sha(source/'PREREG.json'),
        source_single_values=[str(e)for e in energies],local_charged_columns=24,full_local_dirichlet_identity=True,
        exterior_square_identity=True,exact_source_norm_bounds=[12,24],native_matrices_checked=len(primitive),
        directed_primitive_bound_max=str(max(primitive)),one_step_certificates=independent_certs,
        program_duration_rounding=duration_errors,complete_budgets=budget_results,
        independently_decoded_pair_keV=[511*x for x in ratio],observer_receipt_bound=observer_checked,
        full_128step_reexecution=False,full_72273600_circuit_reexecution=False,
        producer_imports=False,literal_new_native_calls=0,seconds=time.perf_counter()-started)
    if args.record:
        with(HERE/'RESULT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        old=read(HERE/'RESULT.json');a=dict(old);b=dict(result);a.pop('seconds');b.pop('seconds');assert a==b
    print(json.dumps({k:v for k,v in result.items()if k not in ['one_step_certificates','complete_budgets','directed_primitive_bound_max']},sort_keys=True))
if __name__=='__main__':main()
