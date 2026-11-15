"""Independent occupation/tree-coordinate review; no producer edits."""
from pathlib import Path
import sys,hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.many_load_native import ManyLoadNative

def require(ok,message):
    if not ok:raise ValueError(message)

def compositions(total,slots):
    if slots==1:yield (total,);return
    for x in range(total+1):
        for tail in compositions(total-x,slots-1):yield (x,)+tail

def run():
    c=Construction(3,((0,1),(0,2)),.23);model=ManyLoadNative(c,(3,0,0,0,0),7)
    occ=list(compositions(3,5));require(len(occ)==35==model.n,'tree occupation sector incomplete')
    # Subdivision tree rooted0: leafports1,2; transit3,4; incidence+atport.
    def records(u):
        f1=-u[1];f3=-u[2]
        f0=u[3]-f1;f2=u[4]-f3
        return tuple(x%7 for x in (f0,f1,f2,f3))
    order=[model.index[(u,records(u))] for u in occ]
    d=np.zeros((35,35));lookup={u:i for i,u in enumerate(occ)}
    for i,u in enumerate(occ):
        d[i,i]=.23*(sum(u[:3])-sum(u[3:]))
        for port,transit,sign in [(0,3,1),(1,3,-1),(0,4,1),(2,4,-1)]:
            for a,b in [(port,transit),(transit,port)]:
                if u[a]:
                    v=list(u);v[a]-=1;v[b]+=1
                    d[lookup[tuple(v)],i]=sign*np.sqrt(u[a]*(u[b]+1))
    herr=float(np.max(abs(d-model.hamiltonian[np.ix_(order,order)])))
    kerr=0.;groups_checked=0
    for law in ('C','E','MIX'):
        for seam,(a,b) in enumerate(c.seams):
            expected={};laws=('C','E') if law=='MIX' else (law,)
            for named in laws:
                for coin in (0,1):
                    for i,u in enumerate(occ):
                        difference=u[a]-u[b]
                        # Independent balanced integer output: coin0 floor(a+b)/2,
                        # coin1 ceil(a+b)/2; transfer follows conservation.
                        target_a=(u[a]+u[b])//2 if coin==0 else (u[a]+u[b]+1)//2
                        transfer=u[a]-target_a
                        v=list(u);v[a]=target_a;v[b]=u[a]+u[b]-target_a
                        key=(named,coin,transfer if named=='C' else difference)
                        expected.setdefault(key,np.zeros((35,35)))[lookup[tuple(v)],i]=1/np.sqrt(2*len(laws))
            actual={(r['law'],r['coin'],r['record_value']):k[np.ix_(order,order)] for r,k in model.operators[law][seam]}
            require(set(actual)==set(expected),'different complete C/E alphabet')
            for key in actual:kerr=max(kerr,float(np.max(abs(actual[key]-expected[key]))));groups_checked+=1
    require(herr<1e-12 and kerr<1e-12,'occupation/tree source mismatch')
    # Whole-position regional charge, independent exact integer lift on the tree.
    cuterr=0
    for u in occ:
        state=(u,records(u))
        for mask in range(32):
            region={i for i in range(5) if mask>>i&1};q=sum(u[i] for i in region);g=3*int(0 in region)
            residue=model.closed_cut_record(state,region)
            cuterr=max(cuterr,abs(model.decode_bounded(residue)-(g-q)))
    require(cuterr==0,'bounded modular cut lift differs')
    return {'schema':'oph.many_load.independent_review.v1','outcome':'PASS_WITH_DECLARED_SCOPE','derivation':'Independent all35 occupation compositions; unique tree half-records; direct second-quantized signed hopping; balanced endpoint arithmetic for every C/E/MIX Kraus group; all32 regional cuts.',
            'Hamiltonian_max_error':herr,'every_Kraus_max_error':kerr,'Kraus_groups_checked':groups_checked,'cut_integer_max_error':cuterr,
            'scope':'Nonnegative Bose fixedtotal3, modulus7, tree; no signedsea or manyload uniformfiber theorem. Cyclic one-load reduction is producer-tested and not independently rerun here. No physical species/statistics inference.',
            'full_criterion_met':False,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*[ROOT/'candidate/oph_openwave'/name for name in ['many_load_native.py','construction.py','uniform_native.py','modular_native.py']]]}}
if __name__=='__main__':
    result=run();Path(__file__).with_name('MANY_LOAD_NATIVE_INDEPENDENT_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
