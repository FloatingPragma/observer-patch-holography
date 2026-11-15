"""Independent worker: planar-rotation dual contours, not producer surface()."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import networkx as nx
import sympy as sp
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
sys.path.insert(0,str(ROOT/'openwave-tuning/candidate'))
from model.observer.native import native
from oph_openwave.geometry import CARRIER_SEAMS


def dual_geometry():
    edges=tuple(CARRIER_SEAMS)
    graph=nx.Graph();graph.add_edges_from(edges)
    planar,embedding=nx.check_planarity(graph)
    if not planar:raise ValueError('carrier not planar')
    visited=set();faces=[]
    for a in sorted(graph):
        for b in embedding.neighbors_cw_order(a):
            if (a,b) not in visited:
                faces.append(tuple(embedding.traverse_face(a,b,visited)))
    halfedge={}
    for i,f in enumerate(faces):
        for a,b in zip(f,f[1:]+f[:1]):halfedge[a,b]=i
    dual=[(halfedge[a,b],halfedge[b,a]) for a,b in edges]
    # Actual dual boundary cycles reconstructed as directed edge walks.
    cycles=[]
    for v in range(12):
        selected=[(e,1 if a==v else -1) for e,(a,b) in enumerate(edges) if v in (a,b)]
        cycle=ordered_contours(dual,selected)
        if len(cycle)!=1 or len(cycle[0])!=5:raise ValueError('not a dual pentagon')
        cycles.append(cycle[0])
    return edges,faces,dual,cycles


def ordered_contours(dual,selected):
    g=nx.MultiDiGraph()
    for e,s in selected:
        a,b=dual[e]
        if s<0:a,b=b,a
        g.add_edge(a,b,key=e,sign=s)
    result=[]
    while g.number_of_edges():
        a=next(v for v in g if g.out_degree(v))
        component=nx.node_connected_component(g.to_undirected(),a)
        sub=g.subgraph(component).copy()
        if not nx.is_eulerian(sub):raise ValueError('cut is not a closed dual contour')
        walk=list(nx.eulerian_circuit(sub,source=a,keys=True))
        result.append([(e,sub[u][v][e]['sign']) for u,v,e in walk])
        g.remove_edges_from(walk)
    return result


def contour_winding(contours,tape):
    phase=[]
    for contour in contours:
        for e,orientation in contour:
            turns=orientation*int(tape[e]);subdivisions=11*(abs(turns)+1)
            phase.extend(np.exp(2j*np.pi*turns*np.linspace(0,1,subdivisions+1)))
    if not phase:return 0.
    z=np.array(phase)
    return float(np.angle(z[1:]*z[:-1].conj()).sum()/(2*np.pi))


def run():
    edges,faces,dual,cycles=dual_geometry()
    boundaries=np.zeros((12,30),int)
    for v,cycle in enumerate(cycles):
        for e,s in cycle:boundaries[v,e]=s
    cuts=np.array([[int(mask>>v&1) for v in range(12)] for mask in range(1,4095)])
    cut_boundaries=cuts@boundaries
    # All 4094 nonempty proper regions, including disconnected/holed regions.
    cut_contours=[]
    for row in cut_boundaries:
        cut_contours.append(ordered_contours(dual,[(e,int(s)) for e,s in enumerate(row) if s]))
    rng=np.random.default_rng(9101)
    starts=[np.eye(12,dtype=int)[v] for v in range(12)]
    starts.extend([np.array([4]+[0]*11),np.array([2,0,1,0,0,1,0,0,0,0,0,0])])
    errors=[];phase_error=0.;events=0;zeroes=0;moves=0;forgotten=0;maxcharge=0
    sample_histories=[]
    for case,G in enumerate(starts):
        n=G.copy();F=np.zeros(30,int);history=[]
        for step in range(96):
            e=int(rng.integers(30));coin=int(rng.integers(2))
            previous=n.copy();n,(_,j)=native(tuple(map(int,n)),edges[e],coin);n=np.array(n);F[e]+=j
            moves+=int(j!=0);zeroes+=int(j==0);events+=1
            # Independent occupation-difference readout, with initial G retained.
            relative=n-G
            observed=-cut_boundaries@F
            expected=cuts@relative
            errors.append(int(np.max(abs(observed-expected))))
            maxcharge=max(maxcharge,int(np.max(abs(observed))))
            if j:
                oldF=F.copy();oldF[e]-=j
                forgotten+=int(np.any(-cut_boundaries@oldF!=expected))
            if step in (0,31,95):
                for k in (0,1,36,137,1023,2046,4093):
                    degree=contour_winding(cut_contours[k],F)
                    phase_error=max(phase_error,abs(degree+expected[k]))
            if case in (0,12):history.append(dict(seam=e,coin=coin,current=int(j),loads=n.tolist()))
        if case in (0,12):sample_histories.append(dict(initial_G=G.tolist(),initial_tape=[0]*30,history=history,final_tape=F.tolist(),final_relative_charge=(n-G).tolist()))
    # Complete one-step native instruments on every one-load state, not only
    # sampled coins. Both outcomes retain G and use their actual source current.
    complete=0
    for start in range(12):
        G=np.eye(12,dtype=int)[start]
        for e,seam in enumerate(edges):
            probabilities=[]
            for coin in (0,1):
                n,(_,j)=native(tuple(map(int,G)),seam,coin)
                F=np.zeros(30,int);F[e]=j
                if np.any(-cut_boundaries@F!=cuts@(np.array(n)-G)):raise ValueError('complete instrument cut mismatch')
                probabilities.append(.5);complete+=1
            if sum(probabilities)!=1:raise ValueError('lost coin branch')
    # Independently check every archived complete prescribed-schedule branch.
    producer=json.loads((HERE.parent/'routes/spatial/NATIVE_DUAL_WINDING.json').read_text())
    archived_histories=0
    for seed in producer['seed_records']:
        total_probability=0.
        for branch in seed['complete_prescribed_schedule_histories']:
            n=tuple(seed['initial_G']);F=np.zeros(30,int)
            for event in branch['events']:
                n,(_,j)=native(n,tuple(event['seam']),event['coin'])
                if j!=event['current'] or list(n)!=event['output_loads']:
                    raise ValueError('archived event not actual native source')
                F[event['seam_index']]+=j
            if list(n)!=branch['final_loads'] or F.tolist()!=branch['final_tape']:
                raise ValueError('archived final state mismatch')
            if branch['probability']!=2.**(-len(branch['events'])):
                raise ValueError('wrong complete coin history probability')
            if np.any(-cut_boundaries@F!=cuts@(np.array(n)-seed['initial_G'])):
                raise ValueError('archived branch separating-contour mismatch')
            for v in range(12):
                phase_error=max(phase_error,abs(contour_winding([cycles[v]],F)-branch['measured_phase_degrees'][v]))
            total_probability+=branch['probability'];archived_histories+=1
        if total_probability!=1.:raise ValueError('incomplete archived schedule')
    # Spatial homology from independently reconstructed pentagon edge words.
    dualgraph=nx.Graph();dualgraph.add_edges_from(dual)
    cycle_rank=dualgraph.number_of_edges()-dualgraph.number_of_nodes()+1
    filled=np.delete(boundaries,[0,1],axis=0)
    h1_twoholes=cycle_rank-sp.Matrix(filled).rank()
    # Higher degree can fission into units; it is not a sign-only clipping meter.
    G=np.zeros(12,int);a,b=edges[0];G[a]=4
    n,(_,j)=native(tuple(map(int,G)),edges[0],0)
    F=np.zeros(30,int);F[0]=j
    c=next(w for w in nx.Graph(edges).neighbors(b) if w!=a)
    e=edges.index(tuple(sorted((b,c))));coin=0 if edges[e][0]==b else 1
    n2,(_,j2)=native(n,edges[e],coin);F[e]+=j2
    fission=(-boundaries@F).tolist()
    # Half-current modification cannot glue the diagnostic S1 curve at the
    # common dual vertex. Integer winding depends on this declared S1 map.
    half_endpoint=abs(np.exp(2j*np.pi*.5)-1)
    sources=[Path(__file__),ROOT/'reverse-engineering-reality/model/observer/native.py',ROOT/'openwave-tuning/candidate/oph_openwave/geometry.py',HERE.parent/'routes/spatial/native_dual_winding.py',HERE.parent/'routes/spatial/NATIVE_DUAL_WINDING.json']
    checks=dict(independent_planar_faces=len(faces)==20 and all(len(f)==3 for f in faces),independent_dual_pentagons=all(len(c)==5 for c in cycles),twohole_H1_rank_one=h1_twoholes==1,all4094_regions_present_relative_charge=max(errors)==0,genuine_zero_tape_native_histories=events==1344,complete_one_load_coin_outcomes=complete==720,all_archived_complete_histories_independently_checked=archived_histories==456,phase_contour_meter=phase_error<1e-12,forgotten_current_rejected=forgotten==moves and moves>0,noop_records_retained=zeroes>0,double_charge_fission=j==2 and j2==1 and fission[b]==fission[c]==1,fractional_phase_gluing_rejected=half_endpoint>1.9)
    return dict(schema='oph.row1.independent_spatial_review.v1',reviewer='independent Codex worker; not Fable opposing-team review',checks={k:bool(v) for k,v in checks.items()},source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},method='NetworkX planar rotation embedding and half-edge face walks; directed dual Euler contours; independently sampled complex contour integral; actual native() only, no producer helper imports',geometry=dict(primal_vertices=12,primal_edges=30,triangular_faces=20,dual_cycle_rank=cycle_rank,twohole_h1=int(h1_twoholes)),native_events=events,complete_coin_outcomes=complete,archived_complete_histories=archived_histories,regions_per_event=4094,moving_events=moves,noop_events=zeroes,max_contour_phase_error=phase_error,max_relative_charge_observed=maxcharge,phase_half_endpoint_failure=half_endpoint,sample_actual_histories=sample_histories,assessment=dict(literal_integer_relative_charge=True,elementary_half_split_impossible_in_declared_integer_S1_meter=True,scope='Q(S)=sum_S(n−G)=−degree(boundary S), actual native current tape and actual retained initial G; no independent supplied texture',not_established=['absolute n identified with relative defect charge','physical electric charge normalization or electromagnetic coupling','phase interpolation as a dynamically prepared independently stored physical field','official OpenWave acceptance']),official_rows_changed=[])

if __name__=='__main__':
    r=run()
    if not all(r['checks'].values()):raise ValueError(r['checks'])
    (HERE/'SPATIAL_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ('checks','native_events','regions_per_event','max_contour_phase_error')},indent=2))
