"""Independent full source-workspace and signed-clock interface audits."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from circuits import array,need,near

HERE=Path(__file__).resolve().parent
LOADS=list(itertools.product(range(-2,4),repeat=3));LI={q:i for i,q in enumerate(LOADS)}
def load(p):return json.loads(gzip.decompress(p.read_bytes()))if p.suffix=='.gz'else json.loads(p.read_text())
def pack(q,c,j):return (LI[tuple(q)]*2+c)*11+j%11
def unpack(i):
    x,j=divmod(int(i),11);q,c=divmod(x,2);return LOADS[q],c,j if j<=5 else j-11

def workspace(folder):
    saved=load(folder/'WORKSPACES.json.gz');seen=set();columns=0;complements=0;states=0
    for e in saved['extensions']:
        edge=tuple(e['edge']);law=e['law'];key=(edge,law);need(key not in seen,'Duplicate extension');seen.add(key)
        v=np.full(4752,-1,dtype=int);used=set();native=[]
        for q in LOADS:
            for c in(0,1):
                a,b=edge;d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2 if law=='native'else d
                out=list(q);out[a]-=j;out[b]+=j;i=pack(q,c,0);o=pack(out,c,j)
                need(o not in used,'Native image collision');used.add(o);v[i]=o
                native.append({'input':i,'output':o,'loads':list(q),'event':[a,b,c,j]})
        need(native==e['native_image'],'Literal full signed source image')
        buckets={}
        for o in range(4752):
            if o not in used:
                q,c,_=unpack(o);buckets.setdefault((sum(q),c),[]).append(o)
        it={s:iter(x)for s,x in buckets.items()};completion=[]
        for i in range(4752):
            if v[i]<0:
                q,c,_=unpack(i);v[i]=next(it[sum(q),c]);completion.append([i,int(v[i])])
        need(completion==e['orthogonal_complement_permutation'],'Total-load/coin conserving complement')
        need(np.array_equal(np.sort(v),np.arange(4752)),'Full V permutation');inverse=np.argsort(v)
        # Independent index-state simulation, including nonzero initial current
        # registers. The native inverse is the inverse of the declared full V.
        result=[]
        for i in range(4752*11):
            source,scratch=divmod(i,11);q,c,current=unpack(v[source]);c^=int(current!=0)
            second=v[pack(q,c,scratch)];q2,c2,scratch2=unpack(second)
            current=(current+scratch2)%11;restored=inverse[pack(q2,c2,scratch2)]
            q3,c3,scratch3=unpack(restored);result.append(pack(q3,c3,current)*11+scratch3%11)
        need(sorted(result)==list(range(52272)),'Full two-register unitary word')
        digest=hashlib.sha256(np.array(result,dtype='<i8').tobytes()).hexdigest()
        need(digest==e['full_workspace_permutation_sha256'],'Full V,V,V-adjoint word differs')
        columns+=len(native);complements+=len(completion);states+=len(result)
    need(seen=={(e,l)for e in((0,1),(1,2),(0,2))for l in('native','swap')},'Workspace law/edge coverage')
    expected={(law,edge,bg,sign,pos,coin) for law in('native','swap')for edge in((0,1),(1,2),(0,2))for bg in range(-1,3)for sign in(-1,1)for pos in range(3)for coin in(0,1)}
    actual=[(r['law'],tuple(r['edge']),r['background'],r['sign'],r['input'].index(r['background']+r['sign']),r['coin'])for r in saved['microtraces']]
    need(len(actual)==len(expected)==saved['count'] and set(actual)==expected,'Signed background trace coverage')
    for r in saved['microtraces']:
        q=r['input'];a,b=r['edge'];c=r['coin'];d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2 if r['law']=='native'else d
        out=q.copy();out[a]-=j;out[b]+=j;cc=c^int(j!=0)
        need(r['after_first']==[out,c,j]and r['after_second']==[q,cc,-j]and r['after_inverse']==[out,cc,0],'Actual signed background microtrace')
        need(r['first_current_after_add']==0 and r['forward_accesses']==2 and r['inverse_accesses']==1,'Current erasure/query count')
    return {'extensions':len(seen),'literal_columns':columns,'complement_columns':complements,'complete_two_register_states':states,'background_microtraces':len(actual)}

def bridge(folder):
    flip=np.arange(33)
    for p in range(3):
        for c in range(10):flip[11*p+c]=11*p+(c^1)
    f=np.eye(33)[flip];maximum=0.;clockcount=0;densitycount=0;paircount=0
    def check(x,y,label):
        nonlocal maximum
        maximum=max(maximum,near(x,y,label,2e-9))
    catalog=load(folder/'CATALOG.json.gz');prep=load(folder/'PREPARATION.json.gz')
    sourcecount=0
    for law in('native','swap','swap_fixed_native_calibration'):
        for pos,neg in zip(catalog['1'][law],catalog['-1'][law]):
            need(pos['counts']==neg['counts']and pos['phase_record']==neg['phase_record'],'Sign-dependent measured norm/phase')
            check(array(pos['K']),array(neg['K']),'Sign-independent accepted source polynomial');sourcecount+=1
    for law in('native','swap'):
        check(array(prep['1'][law]['final_factor']),array(prep['-1'][law]['final_factor']),'Actual preparation sign consistency')
    for k in range(7):
        pos=load(folder/f'INPUT_s1_k{k}.json.gz');neg=load(folder/f'INPUT_s-1_k{k}.json.gz')
        for name in('R','W'):check(array(neg['clock'][name]),f@array(pos['clock'][name])@f,'Signed clock conjugacy '+name)
        check(f@array(pos['clock']['A']),array(pos['clock']['A']),'Choice embedding invariant under coin flip')
        check(array(neg['clock']['A']),array(pos['clock']['A']),'Actual two signed embeddings')
        fp=array(pos['source']['output_optical_factor']);fn=array(neg['source']['output_optical_factor'])
        check(fp@fp.conj().T,fn@fn.conj().T,'Actual decoded mixed density signs');clockcount+=1
        for stride,turns,direction in itertools.product((1,2),(1,2),(-1,1)):
            for law in(('native','swap','swap_fixed_native_calibration')if k in(0,6)else('native','swap')):
                suffix=f'_k{k}_r{stride}_t{turns}_d{direction}_{law}.json.gz'
                p=load(folder/('CASE_s1'+suffix));n=load(folder/('CASE_s-1'+suffix));paircount+=1
                check(p['returned']['path_probabilities_plus_minus'],n['returned']['path_probabilities_plus_minus'],'Both signed path laws')
                for pp,nn in zip(p['returned']['conditional_clocks'],n['returned']['conditional_clocks']):
                    check(array(nn['arc_density']),f@array(pp['arc_density'])@f,'Full returned arc-state signed bridge');densitycount+=1
    return {'coin_flip_permutation':flip.tolist(),'source_catalog_pairs':sourcecount,'clock_refinements':clockcount,'corresponding_case_pairs':paircount,'returned_density_pairs':densitycount,'maximum_matrix_defect':maximum,'meaning':'Exact literal signed conjugacy and actual finite density covariance support all four independently admitted charge-sign products; no marginal-only composition inference.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('--record-to',type=Path);p.add_argument('--verify',type=Path,default=HERE/'WORKSPACE_BRIDGE.json');args=p.parse_args()
    if args.record_to:need(not args.record_to.exists(),'Exclusive first-record path exists')
    start=time.perf_counter();result={'passed':True,'same_group':True,'external_G6':False,'index_sha256':hashlib.sha256((args.folder/'INDEX.json').read_bytes()).hexdigest(),'workspace':workspace(args.folder),'signed_bridge':bridge(args.folder)}
    if args.record_to:
        with args.record_to.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    else:
        saved=load(args.verify)
        old=saved['signed_bridge']['maximum_matrix_defect']
        need(type(old)in(int,float) and 0<=old<2e-9,'Saved bridge residual exceeds its fixed bound')
        saved['signed_bridge']['maximum_matrix_defect']=result['signed_bridge']['maximum_matrix_defect']
        need(saved==result,'Saved bridge evidence mismatch')
    print(json.dumps({'passed':True,'workspace':result['workspace'],'signed_bridge':{k:v for k,v in result['signed_bridge'].items()if k!='coin_flip_permutation'},'seconds':time.perf_counter()-start}))
