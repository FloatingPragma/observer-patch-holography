"""Source-local continuation. No selected mode or eigenvalue is an input."""
from fractions import Fraction
import hashlib,json,math
import numpy as np

def colour(edges):
    groups=[];used=[]
    for a,b in edges:
        for k,s in enumerate(used):
            if a not in s and b not in s:break
        else:k=len(groups);groups.append([]);used.append(set())
        groups[k].append((a,b));used[k].update((a,b))
    return [np.array(x,dtype=int) for x in groups]

def schedule(groups,time,steps):
    # Fixed composition coefficients, independent of source spectrum.
    p=float(1/(4-4**(1/3)));coefficients=[p,p,1-4*p,p,p]
    return [(g,float(c*time/(2*steps))) for c in coefficients
            for g in list(groups)+list(reversed(groups))]

def phase(angle):return complex(math.cos(angle),math.sin(angle))

def matching(state,edges,dt,pair):
    """Exact ideal source block, evaluated in compressed matching form.

    A source edge contributes 2(I-T_e^charged). Contact columns of T_e
    vanish; their coherent word is the diagonal phase exp(-2i dt).
    Every other charged block is the native averaging block.
    """
    z=phase(-2*dt);d=(1+z)/2;o=(1-z)/2;a,b=edges.T
    if not pair:
        x=state[a].copy();y=state[b].copy()
        state[a]=d*x+o*y;state[b]=o*x+d*y
        return
    aa=state[a,a].copy();ab=state[a,b].copy()
    ba=state[b,a].copy();bb=state[b,b].copy()
    x=state[a].copy();y=state[b].copy()
    state[a]=d*x+o*y;state[b]=o*x+d*y
    x=state[:,a].copy();y=state[:,b].copy()
    state[:,a]=d*x+o*y;state[:,b]=o*x+d*y
    # Native killed contact block, not an unrestricted tensor walk.
    state[a,b]=z*ab;state[b,a]=z*ba
    # Collisions are an initial terminal input type; this is a unitary
    # diagonal extension on that unreachable chart, not a physical vacuum.
    state[a,a]=z*aa;state[b,b]=z*bb

def step(state,program,pair):
    out=np.array(state,dtype=complex,copy=True)
    for g,h in program:matching(out,g,h,pair)
    return out

def propagate(state,program,steps,pair,progress=None):
    out=np.array(state,dtype=complex,copy=True)
    for n in range(steps):
        for g,h in program:matching(out,g,h,pair)
        if progress and (n+1)%16==0:progress(n+1)
    return out

def description(edges,time,steps):
    groups=colour(edges);program=schedule(groups,time,steps)
    obj={'edge_groups':[x.tolist() for x in groups],
         'one_step':[{'group':next(i for i,x in enumerate(groups) if x is g),
                      'dt_hex':h.hex()}for g,h in program],
         'repetitions':steps,'time':str(time),
         'phase_formula':'-2*dt for native local generator 2(I-T_e^charged)'}
    text=json.dumps(obj,sort_keys=True,separators=(',',':'))
    return program,obj,hashlib.sha256(text.encode()).hexdigest()

def endpoint(q,a,b,c):
    out=list(q);s=out[a]+out[b];old=out[a]
    out[a]=s//2+c*(s%2);out[b]=s-out[a]
    return out,old-out[a]

def pair_endpoint(p,n,a,b,c):
    q={p:1,n:-1};x=q.get(a,0);s=x+q.get(b,0)
    q[a]=s//2+c*(s%2);q[b]=s-q[a];j=x-q[a]
    q={k:v for k,v in q.items() if v}
    if not q:return None,j
    return (next(k for k,v in q.items()if v==1),
            next(k for k,v in q.items()if v==-1)),j

def local_pair_census():
    # Four labels exhaust all endpoint occupancy types up to relabelling.
    states=[(p,n)for p in range(4)for n in range(4)if p!=n]
    T=np.zeros((12,12),dtype=int);vac=[];records=[]
    for j,(p,n)in enumerate(states):
        for c in (0,1):
            q,current=pair_endpoint(p,n,0,1,c)
            records.append([p,n,c,current,q])
            if q is None:vac.append([p,n,c,current])
            else:T[states.index(q),j]+=1
    L=2*np.eye(12,dtype=int)-T
    assert np.array_equal(L,L.T) and np.array_equal(L@L,2*L)
    return {'states':states,'twice_transfer':T.tolist(),
            'generator':L.tolist(),'vacuum_records':vac,'fine_records':records}

def mcx_count(k):return 2*k-3

def swap_counts(a,b,extra_controls=1):
    path=[a]
    for j in range(8):
        if(a^b)>>j&1:path.append(path[-1]^(1<<j))
    segments=list(zip(path[:-1],path[1:]));ccx=flip=0
    for x,y in segments+segments[-2::-1]:
        target=(x^y).bit_length()-1
        zeros=sum(not(x>>j&1)for j in range(8)if j!=target)
        ccx+=mcx_count(7+extra_controls);flip+=2*zeros
    return ccx,flip

def route_audit(edges):
    """Execute every clean Boolean column of each generated source route.

    These are exact circuit permutations, separately from native M6
    pulse catalogue execution; their complete native comparison is charged.
    """
    checked=0;coupled=0;tables={}
    for extra in (1,2):
        bits=8+extra;initial=np.arange(1<<bits,dtype=np.int64)
        for a,b in edges:
            out=initial.copy();path=[a]
            for j in range(8):
                if(a^b)>>j&1:path.append(path[-1]^(1<<j))
            pairs=list(zip(path[:-1],path[1:]))
            for x,y in pairs+pairs[-2::-1]:
                target=(x^y).bit_length()-1
                controls=[j for j in range(8)if j!=target]+list(range(8,bits))
                zeros=[j for j in range(8)if j!=target and not(x>>j&1)]
                for j in zeros:out^=1<<j
                scratch=list(range(bits,bits+len(controls)-2))
                chain=[(controls[0],controls[1],scratch[0])]
                for j in range(2,len(controls)-1):chain.append((scratch[j-2],controls[j],scratch[j-1]))
                for u,v,w in chain+[(scratch[-1],controls[-1],target)]+list(reversed(chain)):
                    out^=(((out>>u)&1)&((out>>v)&1))<<w
                for j in reversed(zeros):out^=1<<j
            expected=initial.copy();active=(initial>>8)==(1<<extra)-1
            expected[active&((initial&255)==a)]+=b-a
            expected[active&((initial&255)==b)]+=a-b
            assert np.array_equal(out,expected),(a,b,extra)
            checked+=len(initial)
            if extra==2:tables[(a,b)]=out&255
    grid=np.arange(240*240);p0,n0=np.divmod(grid,240);valid=p0!=n0
    p0=np.tile(p0[valid],2);n0=np.tile(n0[valid],2)
    parity=np.repeat(np.arange(2),240*239)
    for a,b in edges:
        ap=(p0==a)|(p0==b);an=(n0==a)|(n0==b)
        # Execute the composed membership / controlled-route / contact-sign
        # / uncompute circuit on every charged pair and parity column.
        p=tables[(a,b)][p0+256*parity+512*(~an)]
        n=tables[(a,b)][n0+256*parity+512*(~ap)]
        sign=1-2*(parity*ap*an)
        scratch_p=ap^((p==a)|(p==b));scratch_n=an^((n==a)|(n==b))
        wantp=p0.copy();wantn=n0.copy();active=parity.astype(bool)
        choose=active&ap&~an;wantp[choose]=a+b-p0[choose]
        choose=active&an&~ap;wantn[choose]=a+b-n0[choose]
        wantsign=np.where(active&ap&an,-1,1)
        assert np.array_equal(p,wantp)and np.array_equal(n,wantn)and np.array_equal(sign,wantsign)
        assert not np.any(scratch_p)and not np.any(scratch_n)
        coupled+=len(p)
    return {'exact_clean_route_truth_columns':checked,'source_edges':len(edges),
            'complete_coupled_pair_reflection_columns':coupled,
            'native_raw_route_replays':0,'scratch_restored':True,
            'type_scope':'Admitted classically tagged charged chart p!=n; no vacuum/charged coherent type superposition is asserted.'}

def word_counts(edges,steps,pair,ccx_pulses):
    # Each symmetric step uses every edge ten times. Two parity-controlled
    # source reflections surround a two-pulse reference/ancilla phase.
    per_reflection=0
    for a,b in edges:
        ccx,flip=swap_counts(a,b,2 if pair else 1)
        if pair:
            # Both source addresses; controls require other-membership=0.
            ccx*=2;flip=2*flip+4
            # Compute/uncompute two membership flags, each two exact
            # eight-bit equalities, using clean scratch.
            ccx+=8*mcx_count(8)
            flip+=8*((8-int(a).bit_count())+(8-int(b).bit_count()))
            # Contact sign: CCZ on parity and both membership flags.
            ccx+=1;extra=2
        else:extra=0
        per_reflection+=ccx*ccx_pulses+2*flip+extra
    local=2*per_reflection+6*len(edges)
    return {'per_step_pulses':10*local,'total_pulses':10*local*steps,
            'source_local_edge_words':10*len(edges)*steps,
            'initial_and_terminal_scratch_bits':24 if pair else 16,
            'compressed_routing_semantics':'Exact Boolean clean-input routes, with native routine comparison and all first-failure prefix bounds.'}
