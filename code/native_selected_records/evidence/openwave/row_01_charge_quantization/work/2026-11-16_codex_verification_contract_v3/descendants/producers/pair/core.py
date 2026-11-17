"""Common integer source/controller; finite records, complete instruments, no package imports."""

# V3 verification contract: these checks are unconditional, also under python -O.
def _require(condition, message="Scientific check failed"):
    if not condition:
        raise AssertionError(message)

import argparse
from collections import deque
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import time
import numpy as np

CELL = ((0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),
        (2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),
        (4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),
        (7,8),(7,10),(8,9),(10,11))
TRI = ((0,1),(1,2),(0,2))
CALLS = {'native':0, 'swap_hold':0, 'swap':0}


def native(q, edge, coin, law='native'):
    CALLS[law] += 1
    a,b = edge
    d = q[a]-q[b]
    j = ((d-(2*coin-1)*(abs(d)%2))//2 if law=='native'
         else d if law=='swap' else coin*d)
    out = list(q); out[a] -= j; out[b] += j
    return tuple(out), [a,b,coin,j]


def graph(bound, m):
    n=12*bound; p,q=n,n+1
    core=[(12*c+a,12*c+b) for c in range(bound) for a,b in CELL]
    core += [(12*c+v,12*(c+1)+v) for c in range(bound-1) for v in range(12)]
    target=[0,p,q]; copy=[n+2,n+3,n+4]; blank=[n+5,n+6]
    rings=[list(range(n+7,n+7+3*m)),list(range(n+7+3*m,n+7+6*m))]
    edges=list(core)
    for vertices in (target,copy):
        edges += [(vertices[a],vertices[b]) for a,b in TRI]
    edges += list(zip(target,copy))+[tuple(blank)]
    for ring in rings:
        edges += list(zip(ring,ring[1:]+ring[:1]))
    bridges=[(0,blank[0]),(0,rings[0][0]),(0,rings[1][0])]
    edges += bridges
    return {'N':n+7+6*m,'core_N':n,'core':core,'edges':edges,'target':target,
            'copy':copy,'blank':blank,'rings':rings,'idle_bridges':bridges}


def tree(n, edges):
    seen={0}; queue=[0]; result=[]
    for a in queue:
        for u,v in edges:
            b=v if u==a else u if v==a else None
            if b is not None and b not in seen:
                seen.add(b);queue.append(b);result.append((a,b))
    _require(len(seen)==n, 'Scientific check at original line 56')
    return result


class Source:
    """Controller receives only event records; state is inaccessible by its interface."""
    def __init__(self, q): self.__q=tuple(q)
    def event(self, edge, coin, law='native'):
        self.__q, record=native(self.__q,edge,coin,law)
        return record


def observe(n, edges, body, probe, cap):
    change=[0]*n; relative=[None]*n; relative[0]=0; cursor=0
    for a,b in tree(n,edges):
        correction=change[a]-change[b]; first=None
        for _ in range(cap):
            if cursor==len(probe): return {'outcome':'probe_timeout'}
            aa,bb,c,j=probe[cursor];cursor+=1
            _require((aa,bb)==(a,b) and c in (0,1) and type(j) is int, 'Scientific check at original line 75')
            change[a]-=j;change[b]+=j
            if first is None: first=(c,j)
            elif c!=first[0]:
                relative[b]=relative[a]-(2*first[1]+j-correction);break
        else: return {'outcome':'probe_timeout'}
    _require(cursor==len(probe), 'Scientific check at original line 81')
    initial=relative.copy()
    for a,b,c,j in reversed(body): initial[a]+=j;initial[b]-=j
    offset=-initial[12];initial=[x+offset for x in initial]
    live=initial.copy();first=None
    for t,record in enumerate(body+[None]):
        if first is None and (set(live)<={0,1} or set(live)<={-1,0}): first=t
        if record:
            a,b,c,j=record;live[a]-=j;live[b]+=j
    before=live.copy()
    for a,b,c,j in probe:live[a]-=j;live[b]+=j
    return {'outcome':'decoded','initial':initial,'before_probe':before,
            'after_probe':live,'first_unit_event':first,'total':sum(initial)}


def shortest(start,goal,edges):
    queue=[start];prev={start:None}
    for a in queue:
        for u,v in edges:
            b=v if a==u else u if a==v else None
            if b is not None and b not in prev:prev[b]=a;queue.append(b)
    path=[goal]
    while path[-1]!=start:path.append(prev[path[-1]])
    return path[::-1]


def charge_route(g, law, total, seed, horizon=8192, cap=32):
    n=g['core_N'];initial=[2,total-2]+[0]*(g['N']-2)
    for ring in g['rings']:initial[ring[2*(len(ring)//3)]]=1
    source=Source(initial);rng=random.Random(seed)
    body=[source.event(g['core'][rng.randrange(len(g['core']))],rng.randrange(2),law)
          for _ in range(horizon)]
    probe=[]
    for edge in tree(n,g['core']):
        first=None
        for _ in range(cap):
            record=source.event(edge,rng.randrange(2));probe.append(record)
            if first is None:first=record[2]
            elif first!=record[2]:break
        else:break
    reading=observe(n,g['core'],body,probe,cap)
    result={'law':law,'total_input_audit':total,'seed':seed,'body':body,'probe':probe,
            'reading':reading,'initial_audit_only':initial,'route':[]}
    if reading['outcome']!='decoded':return result
    _require(reading['initial']==initial[:n], 'Scientific check at original line 125')
    q=reading['after_probe']
    if q.count(1)!=1 or any(x not in (0,1) for x in q):
        result['outcome']='ineligible_target_sector';return result
    here=q.index(1);goal=g['target'][1]
    path=shortest(here,goal,g['core']+[(0,goal)])
    for there in path[1:]:
        edge=tuple(sorted((here,there)))
        for _ in range(cap):
            record=source.event(edge,rng.randrange(2));result['route'].append(record)
            a,b,c,j=record
            if j:
                decoded=b if j==1 else a
                _require(decoded==there, 'Scientific check at original line 138');here=decoded;break
        else:
            result['outcome']='route_timeout';return result
    result['outcome']='target_basis_unit';result['decoded_target']=here
    return result


def fine(edges,law,vacuum=False):
    branches=[];columns=[]
    for e,edge in enumerate(edges):
        for c in (0,1):
            matrices=np.zeros((3,3,3))
            for p in range(3):
                q=tuple(int(x==p) for x in range(3))
                out,rec=native((0,0,0) if vacuum else q,edge,c,law)
                j=rec[3];y=p if vacuum else out.index(1)
                matrices[j+1,y,p]=1
                columns.append([e,c,p,j,y])
            branches.extend(matrices)
    return np.asarray(branches)/math.sqrt(2*len(edges)),columns


def helmert(branches):
    count=len(branches);out=[sum(branches)/math.sqrt(count)];prefix=branches[0].copy()
    for k in range(1,count):
        out.append((prefix-k*branches[k])/math.sqrt(k*(k+1)));prefix+=branches[k]
    return np.asarray(out)


def instrument(edges,law): return helmert(fine(edges,law)[0])


def contrast(edges,law):
    target=instrument(edges,law);blank=helmert(fine(edges,law,True)[0])
    return np.asarray([(b+sign*a)/2 for a,b in zip(target,blank) for sign in (1,-1)])


def routing(slot,m):
    v,k=divmod(slot%(3*m),m);w=(v+1)%3
    return [e for e in TRI if v in e]*(m-k)+[e for e in TRI if w in e]*k


def encode(a):
    a=np.asarray(a)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()} if np.iscomplexobj(a) else a.tolist()


def completeness(ops):
    return float(np.max(abs(np.einsum('kba,kbc->ac',ops.conj(),ops)-np.eye(ops.shape[2]))))


def channel(rho,ops,accept):
    outs=np.einsum('kxa,ab,kyb->kxy',ops,rho,ops.conj())
    probabilities=np.trace(outs,axis1=1,axis2=2).real
    selected=sum(outs[k] for k in accept);p=float(np.trace(selected).real)
    return selected/p,p,probabilities


def clock(k):
    edges=((0,1),(0,1),(0,2),(0,2),(1,2));R=np.eye(33);columns=[]
    for p in range(3):
        for e,edge in enumerate(edges):
            for c in (0,1):
                q=tuple(int(x==p) for x in range(3));out,rec=native(q,edge,c)
                j=rec[3];y=out.index(1);cc=c^int(j!=0)
                i=11*p+2*e+c;o=11*y+2*e+cc;R[:,i]=0;R[o,i]=1
                columns.append([p,e,c,j,y,cc])
    eps=2.**(-k);a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps*eps)]
    A=np.kron(np.eye(3),a[:,None]);P=A.T@R@A
    C=np.kron(np.eye(3),2*np.outer(a,a)-np.eye(11));W=R@C
    v=a.copy();v[0]-=1;H=np.eye(11)-2*np.outer(v,v)/(v@v)
    decoders=np.asarray([np.kron(np.eye(3),H[j:j+1]) for j in range(11)])
    Q=256*2**k;waves=[];live=A.astype(complex)
    for _ in range(Q):waves.append(live.copy());live=W@live
    ops=np.fft.fft(np.asarray(waves),axis=0)/Q
    return {'R':R,'A':A,'P':P,'W':W,'a':a,'decode':decoders,'qpe':ops,'Q':Q,
            'epsilon':eps,'columns':columns,'edges':edges}


def prepare(cl):
    edges=cl['edges'];ops=contrast(edges,'native')
    rho=np.diag([0.,1.,0.]);records=[];reach=0.
    rho,p,probs=channel(rho,ops,[1]);reach=math.log10(p)
    records.append({'catalog':'cooling_contrast','accepted':1,'probabilities':probs.tolist(),
                    'input_density':encode(np.diag([0.,1.,0.])),'log10_reach_before':0.})
    continuing=instrument(edges,'native')
    for n in range(8):
        before=rho.copy();rho,p,probs=channel(rho,continuing,[0])
        records.append({'catalog':'cooling','accepted':0,'probabilities':probs.tolist(),
                        'input_density':encode(before),'log10_reach_before':reach})
        reach+=math.log10(p)
    position=rho.copy();before=rho.copy()
    rho,p,probs=channel(rho,cl['qpe'],range(1,cl['Q']//2))
    records.append({'catalog':'qpe','accepted':list(range(1,cl['Q']//2)),
                    'probabilities':probs.tolist(),'input_density':encode(before),
                    'log10_reach_before':reach})
    reach+=math.log10(p)
    arc=rho.copy();before=rho.copy();rho,p,probs=channel(rho,cl['decode'],[0])
    records.append({'catalog':'decode','accepted':0,'probabilities':probs.tolist(),
                    'input_density':encode(before),'log10_reach_before':reach})
    reach+=math.log10(p)
    return rho,arc,position,reach,records,ops,continuing


def hop(old,desired,ring,graph_size,cap=8):
    edge=tuple(sorted((old,desired)));physical=tuple(ring[x] for x in edge)
    live=[([],[],1.)];history=[]
    for attempt in range(cap):
        nxt=[]
        for coins,currents,p in live:
            for c in (0,1):
                q=[0]*graph_size;q[ring[old]]=1
                out,rec=native(q,physical,c)
                j=rec[3];co=coins+[c];cu=currents+[j];prob=p/2
                if j or attempt==cap-1:
                    decoded=edge[1] if j==1 else edge[0] if j==-1 else old
                    history.append({'coins':co,'currents':cu,'probability':prob,
                                    'decoded':decoded,'outcome':'hop' if j else 'timeout'})
                else:nxt.append((co,cu,prob))
        live=nxt
    success=sum(r['probability'] for r in history if r['outcome']=='hop')
    return success,{'edge':list(edge),'physical_edge':list(physical),'histories':history}


def apply_pair(rho,rot,ref,index):
    # All maps have one COMMON eraser flag; neither arm is measured separately.
    probs=np.zeros(len(rot));outs=[]
    for a,b in zip(rot,ref):
        op=np.zeros((6,6),complex);op[:3,:3]=a;op[3:,3:]=b
        out=op@rho@op.conj().T;outs.append(out);probs[len(outs)-1]=np.trace(out).real
    p=float(probs[index]);return outs[index]/p,p,probs


def ramsey(rho,W,stride):
    unit=np.linalg.matrix_power(W,stride);power=np.eye(len(W),dtype=complex);records=[]
    for n in range(17):
        z=np.trace(rho@power)
        records.append({'time':n,'x_plus':float((1+z.real)/2),'y_plus':float((1+z.imag)/2)})
        power=unit@power
    return records


def read_poles(records):
    z=np.array([2*r['x_plus']-1+1j*(2*r['y_plus']-1) for r in records])
    mat=np.array([[z[n+j] for j in range(4)] for n in range(4)])
    coefficients=np.linalg.solve(mat,-z[4:8]);poles=np.roots(np.r_[1.,coefficients[::-1]])
    order=np.argsort(np.angle(poles));poles=poles[order]
    weights=np.linalg.solve(np.array([[p**n for p in poles] for n in range(4)]),z[:4])
    predicted=np.array([sum(w*p**n for p,w in zip(poles,weights)) for n in range(17)])
    positive=[(float(np.angle(p)),float(w.real)) for p,w in zip(poles,weights)
              if np.angle(p)>0 and w.real>.01]
    return {'frequency':min(positive)[0],'all_poles':encode(poles),'all_weights':encode(weights),
            'heldout_error':float(max(abs(predicted[8:]-z[8:]))),'fit_condition':float(np.linalg.cond(mat))}


def spin_run(rho,m,turns,direction,law,catalog,g):
    rho=np.kron(np.ones((2,2))/2,rho);reach=0.;records=[];rot=ref=2*m
    for step in range(1,3*m*turns+1):
        k=step%m;wanted=((2*m+direction*step)%(3*m),(2*m+direction*min(k,m-k))%(3*m))
        for arm,old,new in zip(('rot','ref'),(rot,ref),wanted):
            p,record=hop(old,new,g['rings'][0 if arm=='rot' else 1],g['N'])
            records.append({'kind':'marker','step':step,'arm':arm,'log10_reach_before':reach,**record})
            reach+=math.log10(p)
        rot,ref=wanted
        rops=catalog[f'{m}:{law}:target:{rot}'];fops=catalog[f'{m}:{law}:target:{ref}']
        for event in range(4):
            before=rho.copy();rho,p,probs=apply_pair(rho,rops,fops,0)
            records.append({'kind':'target','step':step,'event':event,'rot':rot,'ref':ref,
                'input_density':encode(before),'probabilities':probs.tolist(),
                'log10_reach_before':reach,'accepted':0})
            reach+=math.log10(p)
        plane=catalog[f'{law}:plane'];before=rho.copy();rho,p,probs=apply_pair(rho,plane,plane,1)
        records.append({'kind':'plane','step':step,'input_density':encode(before),
                        'probabilities':probs.tolist(),'log10_reach_before':reach,'accepted':1})
        reach+=math.log10(p)
    return rho,reach,records


def reencode(rho,cl):
    blocks=rho.reshape(2,3,2,3);ops=cl['qpe'];pos=ops[1:cl['Q']//2]
    output=np.einsum('kxa,iajb,kyb->ixjy',pos,blocks,pos.conj(),optimize=True).reshape(66,66)
    gram=np.einsum('kxa,kxb->kab',ops.conj(),ops)
    marginal=blocks[0,:,0,:]+blocks[1,:,1,:]
    probs=np.einsum('kab,ba->k',gram,marginal).real
    p=float(np.trace(output).real)
    return output/p,p,probs


def bridge_audit(cl):
    f=np.array([0.,1.,-1.])/math.sqrt(2);g=np.array([-2.,1.,1.])/math.sqrt(6)
    basis=np.stack([f,g],axis=1);A=cl['A'];R=cl['R'];P=cl['P']
    lambdas=np.diag(basis.T@P@basis)
    J=(R@A@basis-A@basis@np.diag(lambdas))@np.diag(1/np.sqrt(1-lambdas*lambdas))
    plus=(A@basis-1j*J)/math.sqrt(2)
    # This Bell state checks off-diagonal reference blocks, not only two rays.
    bell=np.concatenate([plus[:,0],plus[:,1]])/math.sqrt(2)
    op=np.kron(np.eye(2),cl['decode'][0]);decoded=op@bell
    expected=np.concatenate([f,g])/2
    return {'positive_lift':encode(plus),'decode_identity_error':float(np.max(abs(cl['decode'][0]@plus-basis/math.sqrt(2)))),
        'entangled_reference_error':float(np.max(abs(decoded-expected))),
        'isometry_error':float(np.max(abs(plus.conj().T@plus-np.eye(2)))),
        'all_decoder_probability':float(sum(np.linalg.norm(np.kron(np.eye(2),d)@bell)**2 for d in cl['decode']))}


def work(cl):
    edges=((0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(0,3),(1,4),(2,5))
    primitives=[];adj={i:[] for i in range(12)};columns=[]
    for e in edges:
        for edge in (e,e[::-1]):
            perm=[]
            for p in range(6):
                for c in (0,1):
                    out,rec=native(tuple(int(x==p) for x in range(6)),edge,c)
                    j=rec[3];y=out.index(1);cc=c^int(j!=0);perm.append(2*y+cc)
                    columns.append([list(edge),p,c,j,y,cc])
            moved=[i for i in range(12) if perm[i]!=i];a,b=moved;k=len(primitives)
            adj[a].append((b,k));adj[b].append((a,k));primitives.append((edge,perm))
    word=[]
    for x,y in ((0,3),(1,4),(2,5)):
        queue=[2*x];prev={2*x:None}
        for p in queue:
            for q,k in adj[p]:
                if q not in prev:prev[q]=(p,k);queue.append(q)
        path=[];p=2*y
        while p!=2*x:p,k=prev[p];path.append(k)
        path.reverse();word+=path+path[-2::-1]
    mapping=[];traces=[]
    for p0 in range(6):
        p,c=p0,0;trace=[]
        for k in word:
            edge=primitives[k][0];out,rec=native(tuple(int(x==p) for x in range(6)),edge,c)
            j=rec[3];y=out.index(1);cc=c^int(j!=0);trace.append([k,p,c,j,y,cc]);p,c=y,cc
        mapping.append([p,c]);traces.append(trace)
    T=np.zeros((66,66))
    for p,(y,c) in enumerate(mapping):
        for h in range(11):T[11*y+h,11*p+h]=1
    common=np.kron(np.eye(2),cl['W'])
    return T,{'word':word,'columns':columns,'traces':traces,'final_mapping':mapping,
              'word_length':len(word),'clock_commutator':float(np.max(abs(T@common-common@T)))}


def embedding(g,k):
    edge_set={tuple(sorted(e)) for e in g['edges']};rows=[]
    mappings={'energy':g['target'],'spin':[g['target'][1],g['target'][2],0]}
    for kind,mapping in mappings.items():
        for edge in TRI:
            physical=tuple(mapping[x] for x in edge)
            _require(tuple(sorted(physical)) in edge_set, 'Scientific check at original line 385')
            for law in ('native','swap'):
                for p in range(3):
                    for c in (0,1):
                        local,rec=native(tuple(int(x==p) for x in range(3)),edge,c,law)
                        q=[0]*g['N'];q[mapping[p]]=1
                        for ring in g['rings']:q[ring[0]]=1
                        out,global_rec=native(q,physical,c,law)
                        expected=q.copy();expected[mapping[p]]=0;expected[mapping[local.index(1)]]=1
                        _require(tuple(expected)==out and rec[3]==global_rec[3], 'Scientific check at original line 394')
                        rows.append([kind,law,list(physical),p,c,rec[3],local.index(1)])
    for edge in g['core']+g['idle_bridges']:
        _require(tuple(sorted(edge)) in edge_set, 'Scientific check at original line 397')
    tree(g['N'],g['edges'])
    tickets=[]
    for number in range(4**k):
        bits=[(number>>j)&1 for j in range(2*k)];records=[]
        for c in bits:
            q=[0]*g['N'];out,rec=native(q,tuple(g['blank']),c);records.append(rec)
            _require(not any(out) and rec[3]==0, 'Scientific check at original line 404')
        tickets.append({'bits':bits,'records':records,'active':not any(bits),'probability':1/(4**k)})
    return {'columns':rows,'blank_selector_tickets':tickets,'connected':True}


def run(k):
    for law in CALLS: CALLS[law]=0
    start=time.perf_counter();cl=clock(k);g4=graph(2,4);g8=graph(2,8)
    charge=[]
    for m,g in ((4,g4),(8,g8)):
        for total in (0,1,2):
            for law in ('native','swap_hold'):
                record=charge_route(g,law,total,72315000+total)
                record['M']=m;charge.append(record)
    position,arc,cooled,prep_log,prep_records,cool_contrast,cooling=prepare(cl)
    # energy order(center,p,q) to spin order(p,q,center), same physical excitation.
    perm=np.eye(3)[[1,2,0]];spin_input=perm@position@perm.T
    catalog={'cooling_contrast':cool_contrast,'cooling':cooling,'qpe':cl['qpe'],'decode':cl['decode']}
    for law in ('native','swap'):
        catalog[f'{law}:plane']=contrast(TRI,law)
        for m in (4,8):
            for slot in range(3*m):catalog[f'{m}:{law}:target:{slot}']=instrument(routing(slot,m),law)
    bridge=bridge_audit(cl);T,transport=work(cl);cases=[];embeddings=[embedding(g4,k),embedding(g8,k)]
    for m in (4,8):
        for turns in (1,2):
            for direction in (1,-1):
                for law in ('native','swap'):
                    spin,logp,records=spin_run(spin_input,m,turns,direction,law,catalog,g4 if m==4 else g8)
                    change=np.kron(np.eye(2),perm.T);energy=change@spin@change.T
                    final,p,qprobs=reencode(energy,cl);logp+=math.log10(p)
                    fblocks=final.reshape(2,33,2,33)
                    raw_path=[float(np.trace((fblocks[0,:,0,:]+fblocks[1,:,1,:]+sign*(fblocks[0,:,1,:]+fblocks[1,:,0,:]))/2).real) for sign in (1,-1)]
                    visibility=2*raw_path[0]-1
                    selected=0 if turns==2 or law=='swap' else 1;sign=1 if selected==0 else -1
                    arc_out=(fblocks[0,:,0,:]+fblocks[1,:,1,:]+sign*(fblocks[0,:,1,:]+fblocks[1,:,0,:]))/(2*raw_path[selected])
                    raw=ramsey(arc_out,cl['W'],2**k)
                    try: reader=read_poles(raw)
                    except (np.linalg.LinAlgError,ValueError) as error:
                        reader={'failed':type(error).__name__,'frequency':None,'heldout_error':None}
                    plus=np.asarray(bridge['positive_lift']['real'])+1j*np.asarray(bridge['positive_lift']['imag'])
                    fidelity=float(np.vdot(plus[:,0],arc_out@plus[:,0]).real)
                    before=np.zeros((66,66),complex);before[:33,:33]=arc_out;after=T@before@T.T
                    workerror=float(np.max(abs(after[33:,33:]-arc_out)))
                    cases.append({'M':m,'turns':turns,'direction':direction,'law':law,
                        'records':records,'spin_final_density':encode(spin),'reencode_input':encode(energy),
                        'all_reencode_bin_probabilities':qprobs.tolist(),'final_density':encode(final),
                        'path_probabilities':raw_path,'visibility':visibility,'clock_path_flag':selected,
                        'clock_density':encode(arc_out),'ramsey_records':raw,'reader':reader,
                        'positive_fidelity':fidelity,'transport_actual_state_error':workerror,
                        'stable_reader':cosine_reader(raw),
                        'audited_off_compact_band_weight':float(1-fidelity-np.vdot(plus[:,0].conj(),arc_out@plus[:,0].conj()).real),
                        'log10_conditional_success':logp,'log10_total_with_preparation':prep_log+logp})
    expected=2*math.asin(cl['epsilon']/math.sqrt(5))/cl['epsilon']
    gates={'charge_selection':all(r['reading']['first_unit_event'] is not None for r in charge if r['law']=='native'),
        'charge_control_retains_multiload':all(r['reading']['first_unit_event'] is None for r in charge if r['law']=='swap_hold'),
        'record_driven_same_target':all(r.get('outcome')=='target_basis_unit' for r in charge if r['law']=='native' and r['total_input_audit']==1),
        'ineligible_inputs_retained':all(r.get('outcome')=='ineligible_target_sector' for r in charge if r['total_input_audit']!=1),
        'bridge_reference':bridge['entangled_reference_error']<1e-12,
        'complete_instruments':max(completeness(v) for v in catalog.values())<1e-10,
        'native_sign':all((c['visibility']<-.9 if c['turns']==1 else c['visibility']>.9) for c in cases if c['law']=='native'),
        'swap_control_positive':all(abs(c['visibility']-1)<1e-10 for c in cases if c['law']=='swap'),
        'same_object_frequency':all(c['stable_reader']['lower']<=expected<=c['stable_reader']['upper'] and c['stable_reader']['width']<.01 and c['stable_reader']['imaginary_first_moment']>.5 for c in cases),
        'fixed_off_band_bound':all(-1e-10<=c['audited_off_compact_band_weight']<=.002+1e-10 for c in cases),
        'heldout_real_moments':all(max(abs(2*r['x_plus']-1-math.cos(r['time']*expected)) for r in c['ramsey_records'][2:])<=.004+1e-10 for c in cases),
        'positive_mode_retained':min(c['positive_fidelity'] for c in cases)>.98,
        'native_copy_transport':transport['clock_commutator']<1e-12 and max(c['transport_actual_state_error'] for c in cases)<1e-12}
    output={'candidate':'A-native-common-refinement-family-v1','graphs':[g4,g8],'charge_records':charge,
        'clock':{'k':k,'epsilon':cl['epsilon'],'Q':cl['Q'],'columns':cl['columns'],'R':encode(cl['R']),'W':encode(cl['W']),'A':encode(cl['A'])},
        'bridge':bridge,'embeddings':embeddings,'preparation_records':prep_records,'cooled_position_density':encode(cooled),
        'positive_counter_density':encode(arc),'decoded_position_density':encode(position),
        'catalog':{k:encode(v) for k,v in catalog.items()},'transport':transport,'cases':cases,
        'gates':gates,'passed':all(gates.values()),'literal_native_calls':dict(CALLS),
        'seconds':time.perf_counter()-start,'expected_frequency':expected,'same_group':True,'physical_award':False}
    return output


def cosine_reader(records, eta=.002):
    """Record-only stable interval; eta is fixed before this family execution."""
    real=2*records[1]['x_plus']-1
    low=max(-1.,min(1.,(real-eta)/(1-eta)))
    high=max(-1.,min(1.,(real+eta)/(1-eta)))
    a,b=math.acos(high),math.acos(low)
    return {'lower':a,'upper':b,'width':b-a,'fixed_eta':eta,
            'imaginary_first_moment':2*records[1]['y_plus']-1}
