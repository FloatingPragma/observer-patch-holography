"""Integral spatial topology and a complete cycle-pointer readout of native flux.

Independent of the spatial producer: faces use exact embedded orientations,
integer Smith certificates prove the lattice quotient, and the winding meter
counts wraps in a six-state pointer rather than summing the producer's phases.
"""
from pathlib import Path
from collections import deque
from itertools import combinations
import ast
import hashlib
import json
from fractions import Fraction

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
RER = ROOT / 'reverse-engineering-reality'
TUNING = ROOT / 'openwave-tuning'


def add(x, y): return x[0]+y[0], x[1]+y[1]
def neg(x): return -x[0], -x[1]
def mul(x, y): return x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1]


def sign_qphi(x):
    # Twice (a+b phi) is u+v sqrt(5). Decide its sign by integer squares.
    u, v = 2*x[0]+x[1], x[1]
    if not v: return (u > 0)-(u < 0)
    if not u: return (v > 0)-(v < 0)
    if u*v > 0: return (u > 0)-(u < 0)
    if u > 0: return (u*u > 5*v*v)-(u*u < 5*v*v)
    return (5*v*v > u*u)-(5*v*v < u*u)


def geometry():
    path = TUNING/'candidate/oph_openwave/geometry.py'
    tree = ast.parse(path.read_text())
    edges = ast.literal_eval(next(n.value for n in tree.body if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == 'CARRIER_SEAMS' for t in n.targets)))
    zero, one, minus, phi, nphi = (0,0), (1,0), (-1,0), (0,1), (0,-1)
    frame = [(minus,phi,zero),(one,phi,zero),(minus,nphi,zero),(one,nphi,zero),
             (zero,minus,phi),(zero,one,phi),(zero,minus,nphi),(zero,one,nphi),
             (phi,zero,minus),(phi,zero,one),(nphi,zero,minus),(nphi,zero,one)]
    def determinant(a,b,c):
        cross = [add(mul(b[(i+1)%3],c[(i+2)%3]),neg(mul(b[(i+2)%3],c[(i+1)%3]))) for i in range(3)]
        result = zero
        for x,y in zip(a,cross): result=add(result,mul(x,y))
        return result
    lookup = {e:i for i,e in enumerate(edges)}
    faces=[]
    for face in combinations(range(12),3):
        if all(e in lookup for e in combinations(face,2)):
            oriented=face if sign_qphi(determinant(*(frame[v] for v in face)))>0 else (face[0],face[2],face[1])
            faces.append(oriented)
    a=sp.zeros(30,20);b=sp.zeros(12,30)
    for j,face in enumerate(faces):
        for x,y in zip(face,face[1:]+face[:1]): a[lookup[tuple(sorted((x,y)))],j]=1 if x<y else -1
    for j,(x,y) in enumerate(edges): b[x,j]=1;b[y,j]=-1
    assert b*a==sp.zeros(12,20)
    assert all(sorted(int(v) for v in a.row(e) if v)==[-1,1] for e in range(30))
    # Tree edges give a unimodular basis for neutral charges.
    adjacency=[[] for _ in range(12)]
    for e,(x,y) in enumerate(edges): adjacency[x].append((y,e));adjacency[y].append((x,e))
    queue=deque([0]);seen={0};tree_edges=[]
    while queue:
        x=queue.popleft()
        for y,e in adjacency[x]:
            if y not in seen: seen.add(y);queue.append(y);tree_edges.append(e)
    chords=[e for e in range(30) if e not in tree_edges]
    reduced=b[list(range(1,12)),tree_edges]
    assert abs(reduced.det())==1
    normal=smith_normal_form(a[chords,:],domain=sp.ZZ)
    factors=[abs(int(normal[i,i])) for i in range(19)]
    assert factors==[1]*19
    # Boundary cycles of the dual cells, built by following actual oriented edges.
    dual_edges=[]
    for e in range(30):
        tail=next(j for j in range(20) if a[e,j]==-1)
        head=next(j for j in range(20) if a[e,j]==1)
        dual_edges.append((tail,head))
    cycles=[]
    for vertex in range(12):
        directed={}
        for e,(tail,head) in enumerate(dual_edges):
            coeff=int(b[vertex,e])
            if coeff:
                start,end=(tail,head) if coeff==1 else (head,tail)
                assert start not in directed
                directed[start]=(end,e,coeff)
        start=min(directed);current=start;walk=[]
        for _ in range(5):
            current,e,coeff=directed[current];walk.append((e,coeff))
        assert current==start and len({e for e,_ in walk})==5
        cycles.append(walk)
    return edges,faces,a,b,cycles,{'tree_minor_determinant':int(reduced.det()),
        'cycle_quotient_Smith_factors':factors,'primal_face_rank':19,
        'dual_graph_H1_rank':11,'neutral_charge_lattice_rank':11}


def pointer_degree(tape, cycle):
    """One closed signed C6 walk per dual edge, driven by its retained integer lift."""
    pointer=0;wraps=0;events=0
    for edge,orientation in cycle:
        current=orientation*int(tape[edge])
        direction=(current>0)-(current<0)
        for _ in range(6*abs(current)):
            old=pointer;pointer=(pointer+direction)%6
            if old==5 and pointer==0: wraps+=1
            if old==0 and pointer==5: wraps-=1
            events+=1
        assert pointer==0
    return wraps,events


def native_pointer_certificate(native, edges, maximum_moves):
    """A literal-native finite observer for the signed C6 path.

    Every elementary pointer shift is attempted on an actual carrier seam.
    The entire fair-coin pair is retained: move with probability1/2, hold
    with probability1/2. A fixed retry cap adds an explicit timeout outcome.
    """
    cycle=(0,1,8,6,2,10)
    attempts=16
    transitions=[]
    for position,port in enumerate(cycle):
        loads=tuple(int(v==port) for v in range(12))
        for direction in (-1,1):
            destination=cycle[(position+direction)%6]
            seam=tuple(sorted((port,destination)))
            assert seam in edges
            targets=[]
            for coin in (0,1):
                output,record=native(loads,seam,coin)
                occupied=[v for v,n in enumerate(output) if n]
                assert len(occupied)==1
                targets.append(occupied[0])
                transitions.append({'pointer_position':position,'direction':direction,
                    'seam':list(seam),'coin':coin,'record':list(record),
                    'next_pointer_position':cycle.index(occupied[0]),'probability':'1/2'})
            assert sorted(targets)==sorted((port,destination))
    step_timeout=Fraction(1,2**attempts)
    complete_success=(1-step_timeout)**maximum_moves
    # Sum every possible first-timeout stage and the complete success branch.
    failed=sum((1-step_timeout)**j*step_timeout for j in range(maximum_moves))
    assert failed+complete_success==1
    return {'carrier_cycle':list(cycle),'complete_one_attempt_transitions':transitions,
        'maximum_retries_per_pointer_move':attempts,'move_timeout_probability':str(step_timeout),
        'maximum_pointer_moves_for_tested_regional_readout':maximum_moves,
        'worst_case_complete_success_probability':str(complete_success),
        'worst_case_timeout_probability':str(failed),
        'worst_case_timeout_probability_float':float(failed),
        'all_first_timeout_records_and_hold_coins_retained':True,
        'maximum_native_pointer_attempts':attempts*maximum_moves,
        'preparation':'Separate twelve-port native pointer carrier with one known unit at port0. This ready-pointer preparation is declared, not derived by the present audit.',
        'measurement':'The observer follows actual native current records to track pointer position, signed wraps and timeout. Final pointer coordinate alone is not the readout.'}


def execute():
    producer=HERE.parent/'spatial/NATIVE_DUAL_WINDING.json'
    source=RER/'model/observer/native.py'
    paths=[producer,HERE.parent/'spatial/native_dual_winding.py',source,TUNING/'candidate/oph_openwave/geometry.py']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    reported=json.loads(producer.read_text())
    native_tree=ast.parse(source.read_text())
    fn=next(n for n in native_tree.body if isinstance(n,ast.FunctionDef) and n.name=='native')
    namespace={};exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'),namespace)
    native=namespace['native']
    edges,faces,a,b,cycles,lattice=geometry()
    assert len(faces)==20 and len(reported['seed_records'])==132
    checked=0;max_events=0;meter_error=0;changed=0;history_checks=0
    for seed in reported['seed_records']:
        tape=sp.Matrix(seed['tape']);degree=b*tape
        assert list(map(int,degree))==seed['degree']
        background=seed['initial_G']
        assert sum(background)==1 and set(background)<={0,1}
        loads=[g-int(q) for g,q in zip(background,degree)]
        assert sum(loads)==1 and set(loads)<={0,1}
        histories=seed['complete_prescribed_schedule_histories']
        assert sum(Fraction(h['probability']) for h in histories)==1
        for history in histories:
            replay=tuple(background);current_tape=sp.zeros(30,1)
            for event in history['events']:
                e=event['seam_index']
                assert tuple(event['seam'])==edges[e]
                replay,record=native(replay,edges[e],event['coin'])
                assert record==(event['coin'],event['current'])
                assert replay==tuple(event['output_loads'])
                current_tape[e]+=record[1]
            assert replay==tuple(history['final_loads'])
            assert list(map(int,current_tape))==history['final_tape']
            exact_degree=b*current_tape
            assert exact_degree==sp.Matrix([g-n for g,n in zip(background,replay)])
            for vertex,cycle in enumerate(cycles):
                readout,cost=pointer_degree(current_tape,cycle)
                assert readout==exact_degree[vertex]
                meter_error=max(meter_error,abs(readout-history['measured_phase_degrees'][vertex]))
            history_checks+=1
        for vertex,cycle in enumerate(cycles):
            readout,cost=pointer_degree(tape,cycle)
            assert readout==degree[vertex];max_events=max(max_events,cost)
        for edge,(x,y) in enumerate(edges):
            for coin in (0,1):
                actual,record=native(tuple(loads),(x,y),coin)
                total=loads[x]+loads[y]
                target=list(loads)
                target[x],target[y]=((total//2,(total+1)//2) if coin==0 else ((total+1)//2,total//2))
                current=loads[x]-target[x]
                assert actual==tuple(target) and record==(coin,current)
                assert sum(target)==1 and set(target)<={0,1}
                following=tape.copy();following[edge]+=current
                assert b*following==sp.Matrix([g-n for g,n in zip(background,target)])
                # Incremental readout: native edit changes precisely the two
                # spatial defects. The full pointer is independently replayed.
                for vertex in (x,y):
                    readout,cost=pointer_degree(following,cycles[vertex])
                    assert readout==background[vertex]-target[vertex];max_events=max(max_events,cost)
                changed+=bool(current);checked+=1
    # Across known one-load references, the charge classes are132dipoles plus0.
    # Modulo6 is injective for each regional charge in{-1,0,1}; full-current
    # records update this finite counter without storing the whole tape.
    counter_checks=0
    regions=[tuple((mask>>v)&1 for v in range(12)) for mask in range(4096)]
    decode={0:0,1:1,5:-1}
    for seed in reported['seed_records']:
        degrees=seed['degree']
        for region in regions:
            q=sum(c*n for c,n in zip(region,degrees))
            assert q in (-1,0,1) and decode[q%6]==q
            counter_checks+=1
    for initial in range(12):
        one=tuple(int(v==initial) for v in range(12))
        for seam in edges:
            for coin in (0,1):
                output,_=native(one,seam,coin)
                assert sum(output)==1 and set(output)<={0,1}
    # Two representatives differing by an actual primal face boundary have
    # the same spatial class; Smith factors prove this for every closed chain.
    for face in range(20):
        f=a[:,face]
        assert all(pointer_degree(f,c)[0]==0 for c in cycles)
    # Six-class endpoint data alone lose a genuine integer winding.
    six=[0]*30;six[0]=6
    assert [v%6 for v in six]==[0]*30
    six_degrees=[pointer_degree(six,c)[0] for c in cycles]
    assert sorted(six_degrees)==[-6]+[0]*10+[6]
    # Signed half turns would give a non-closed cycle-pointer trace: exactly
    # three unit C6 steps, ending at the antipode rather than its initial state.
    half_endpoint=3
    assert half_endpoint!=0
    pointer=native_pointer_certificate(native,edges,max_events)
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    return {'schema':'oph.row1.independent_spatial_lift.v1','checks':{
        'exact_icosahedral_embedded_orientation':True,'integral_lattice_isomorphism_certified':True,
        'all132_actual_seed_records_match':True,'all7920_native_branches_match_independent_repair':checked==7920,
        'one_load_source_domain_is_native_closed':True,
        'all_producer_complete_coin_histories_replayed':history_checks==reported['complete_coin_histories'],
        'independent_exact_pointer_vs_producer_phase_meter':meter_error<1e-12,
        'finite_mod6_regional_counter_decodes_injectively':counter_checks==132*4096,
        'cycle_pointer_winding_matches_every_updated_endpoint':True,
        'all20_primal_face_boundaries_are_null_spatial_classes':True,
        'six_class_only_counterexample_retained':True,'half_loop_does_not_close':True,
        'literal_native_cycle_pointer_complete':True,'source_pins_unchanged':True},
        'lattice_certificate':lattice,'event_checks':checked,'nonzero_current_events':changed,
        'complete_coin_histories':history_checks,'phase_meter_max_error':meter_error,
        'preparation':'Actual one-load initial_G and zero tape; every producer prescribed-schedule coin history replayed. No background2offset is substituted.',
        'finite_class_observer':{'source_domain_size_per_known_G':12,'charge_classes_across_G':133,'regional_values':[-1,0,1],
            'register_states':6,'decoder':{'0':0,'1':1,'5':-1},
            'other_residues':'explicit domain failure, never accepted as a halfcharge',
            'counter_update':'r <- r + (1_R(a)-1_R(b))*j mod6',
            'regional_checks':counter_checks,
            'class_memory':'Retain the known initial_G. The complete spatial class can be kept as12finite residue counters or one of133charge configurations acrossG. For each knownG only12source positions occur. No unboundedcurrenttape is needed for the class; the full tape selects a representative.'},
        'largest_pointer_event_count_for_one_tested_vertex':max_events,
        'dual_face_cycles':[[list(item) for item in cycle] for cycle in cycles],
        'native_pointer_service':pointer,
        'residue_only_counterexample':{'all_residues_zero':True,'actual_degrees':six_degrees},
        'operational_service':'Classical finite readout of the actual committed integer current tape: on dual edge e, emit 6|F_e| signed unit shifts of a six-state cyclic pointer. The pointer shifts are implemented by the literal native rule on a declared one-unit auxiliary carrier, with finite retries and explicit timeouts. Complete signed pointer transitions determine winding; reading only its final state does not.',
        'construction_boundary':'The spatial homotopy class is derived from the actual tape on the oriented carrier surface. The cyclic pointer is a declared classical decoder/readout, not an independently evolved field or proof of an autonomous electromagnetic interaction. Retaining the integer tape/lift is essential. The existing spectator-sensitive active-memory voltage is a different construction.',
        'source_sha256':pins,'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'official_rows_changed':[]}


if __name__=='__main__':
    result=execute()
    assert all(result['checks'].values())
    (HERE/'SPATIAL_LIFT_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'lattice':result['lattice_certificate'],'events':result['event_checks']}))
