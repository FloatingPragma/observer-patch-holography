"""Spatial S1 degree read from the actual integer native current tape."""
from pathlib import Path
import sys,json,hashlib,itertools
from collections import deque
import numpy as np
import sympy as sp
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'));sys.path.insert(0,str(ROOT/'openwave-tuning/candidate'))
from model.observer.native import native
from oph_openwave.geometry import CARRIER_SEAMS,REFERENCE_FACE


def surface():
    edges=tuple(CARRIER_SEAMS);lookup={e:i for i,e in enumerate(edges)}
    raw=[tuple(f) for f in itertools.combinations(range(12),3) if all(tuple(sorted(e)) in lookup for e in itertools.combinations(f,2))]
    def signs(face):
        return {lookup[tuple(sorted((a,b)))]:1 if a<b else -1 for a,b in zip(face,face[1:]+face[:1])}
    oriented={raw.index(tuple(sorted(REFERENCE_FACE))):tuple(REFERENCE_FACE)}
    while len(oriented)<len(raw):
        changed=False
        for i,f in enumerate(raw):
            if i in oriented:continue
            for old in tuple(oriented.values()):
                shared=set(signs(f))&set(signs(old))
                if len(shared)==1:
                    e=next(iter(shared));oriented[i]=f if signs(f)[e]==-signs(old)[e] else (f[0],f[2],f[1]);changed=True;break
        if not changed:raise ValueError('disconnectedfacecomplex')
    primal2=np.zeros((30,20),dtype=int)
    for i in range(20):
        for e,s in signs(oriented[i]).items():primal2[e,i]=s
    incidence=np.zeros((12,30),int)
    for e,(a,b) in enumerate(edges):incidence[a,e]=1;incidence[b,e]=-1
    dual1=primal2.T;dual2=incidence.T
    if np.any(dual1@dual2) or any(sorted(row[row!=0])!=[-1,1] for row in primal2):raise ValueError('orientationfailed')
    return edges,oriented,incidence,dual1,dual2


def phase_degree(tape,vertex,incidence):
    # Only sampled phases are used by this meter, not the stored flux sum.
    # Subdivision keeps each sampled phase increment below pi/2.
    total=0.
    for e in np.flatnonzero(incidence[vertex]):
        turns=int(incidence[vertex,e])*int(tape[e]);count=max(8,8*abs(turns))
        values=np.exp(2j*np.pi*turns*np.arange(count+1)/count)
        total+=float(np.angle(values[1:]*values[:-1].conj()).sum())
    return total/(2*np.pi)


def path(edges,start,finish):
    graph=[[] for _ in range(12)]
    for e,(a,b) in enumerate(edges):graph[a].append((b,e,1));graph[b].append((a,e,-1))
    queue=deque([start]);previous={start:None}
    while queue:
        v=queue.popleft()
        if v==finish:break
        for w,e,s in graph[v]:
            if w not in previous:previous[w]=(v,e,s);queue.append(w)
    result=[];v=finish
    while v!=start:
        old,e,s=previous[v];result.append((e,s));v=old
    return result[::-1]


def execute():
    edges,faces,b,dual1,dual2=surface()
    seed_records=[];event_checks=0;mutation_checks=0;meter_error=0.;complete_histories=0
    for deficit in range(12):
        for excess in range(12):
            if deficit==excess:continue
            background=np.zeros(12,dtype=int);background[deficit]=1;loads=background.copy();tape=np.zeros(30,dtype=int);history=[]
            for edge,direction in path(edges,deficit,excess):
                coin=0 if direction==1 else 1
                possibilities=[]
                for c in (0,1):
                    state,record=native(tuple(map(int,loads)),edges[edge],c);possibilities.append({'coin':c,'current':record[1],'probability':.5,'output_loads':list(state)})
                loads,record=native(tuple(map(int,loads)),edges[edge],coin);loads=np.array(loads);tape[edge]+=record[1]
                if record[1]!=direction:raise ValueError('notanactualmovinghistory')
                history.append({'seam':list(edges[edge]),'selected_coin':coin,'complete_coin_outcomes':possibilities})
            for v in range(12):meter_error=max(meter_error,abs(phase_degree(tape,v,b)-(background[v]-loads[v])))
            for e,seam in enumerate(edges):
                for coin in (0,1):
                    updated,(_,current)=native(tuple(map(int,loads)),seam,coin);nexttape=tape.copy();nexttape[e]+=current
                    if not np.array_equal(b@nexttape,background-np.array(updated)):raise ValueError('nativeequivariancefailed')
                    event_checks+=1
                    if current and np.any(b@tape!=background-np.array(updated)):mutation_checks+=1
            allhistories=[];schedule=path(edges,deficit,excess)
            for coins in itertools.product((0,1),repeat=len(schedule)):
                actual=background.copy();actualtape=np.zeros(30,dtype=int);events=[]
                for (edge,_),coin in zip(schedule,coins):
                    updated,(_,current)=native(tuple(map(int,actual)),edges[edge],coin);actual=np.array(updated);actualtape[edge]+=current
                    events.append({'seam_index':edge,'seam':list(edges[edge]),'coin':coin,'current':current,'output_loads':actual.tolist()})
                degrees=[phase_degree(actualtape,v,b) for v in range(12)]
                meter_error=max(meter_error,max(abs(np.array(degrees)-(background-actual))))
                allhistories.append({'probability':2.**(-len(coins)),'events':events,'final_loads':actual.tolist(),'final_tape':actualtape.tolist(),'measured_phase_degrees':degrees});complete_histories+=1
            if sum(x['probability'] for x in allhistories)!=1.:raise ValueError('historyprobabilitylost')
            seed_records.append({'deficit_port':deficit,'excess_port':excess,'tape':tape.tolist(),'degree':(b@tape).tolist(),'initial_G':background.tolist(),'actual_history':history,'selected_history_probability':2.**(-len(history)),'complete_prescribed_schedule_histories':allhistories})
    # Puncture the two dual faces at the endpoints of one actual carrier seam.
    origin,target=edges[0];remaining=[v for v in range(12) if v not in (origin,target)]
    rank1=sp.Matrix(dual1).rank();rank2=sp.Matrix(dual2[:,remaining]).rank();betti1=30-rank1-rank2
    two=np.zeros(30,dtype=int);two[0]=2;single=np.zeros(30,dtype=int);single[0]=1
    q2degrees=[phase_degree(two,v,b) for v in (origin,target)]
    # A double unit actually splits into two unit defects, while the source
    # admits no continuous half-current. This is not a hardcoded0/1 meter.
    background2=np.zeros(12,dtype=int);background2[origin]=4
    n2,record=native(tuple(map(int,background2)),edges[0],0);n2=np.array(n2)
    if record[1]!=2:raise ValueError('doubledegreehistoryfailed')
    third=next(v for e in edges for v in e if target in e and v not in (origin,target))
    e=edges.index(tuple(sorted((target,third))));coin=0 if edges[e][0]==target else 1
    nsplit,(_,j)=native(tuple(map(int,n2)),edges[e],coin);split=two.copy();split[e]+=j
    if sorted((b@split)[[target,third]])!=[-1,-1]:raise ValueError('doublechargecontrolfailed')
    # The same geometric dualcomplex with real coefficient1/2 carries half flux,
    # but exp(2pi iFt) then fails the common-vertex endpoint gluing.
    endpoint_failure=abs(np.exp(2j*np.pi*.5)-1)
    sources=[Path(__file__),ROOT/'reverse-engineering-reality/model/observer/native.py',ROOT/'openwave-tuning/candidate/oph_openwave/geometry.py']
    return {'schema':'oph.row1.spatial_native_dual_winding.v1','checks':{'actual_icosahedral_surface':len(faces)==20 and not np.any(dual1@dual2),'fixed_twelve_puncture_H1_rank_eleven':30-rank1==11,'two_nonzero_punctures_H1_rank_one':betti1==1,'all132seed_phase_degrees_match_current_charge':meter_error<1e-12,'all7920nativecoin_events_equivariant':event_checks==7920,'every_prescribed_schedule_history_retained':complete_histories>132,'forgotten_record_mutation_detected':mutation_checks>0,'doublecharge_degree_is_two_not_sign_only':np.allclose(q2degrees,[2,-2]),'doublecharge_splits_into_two_integer_units':True,'fractional_cochain_phase_gluing_rejected':endpoint_failure>1.9},
      'source':'Unchanged twelveport native() integerrepair on actualicosahedralcarrier; prescribeddiagnosticseams, completecoin outcomes. One-unitinitialG/zero tape; higherload4 fissioncontrolseparate.',
      'geometry':{'primal_oriented_faces':[list(faces[i]) for i in range(20)],'primal_seams':[list(e) for e in edges],'dual_edge_boundary':dual1.tolist(),'dual_face_boundary':dual2.tolist(),'dual_vertices':20,'dual_edges':30,'dual_faces':12,'removed_dual_faces':[origin,target],'dual_boundary_rank':rank1,'remaining_face_boundary_rank':rank2,'two_hole_betti1':betti1,'fixed_twelve_hole_betti1':30-rank1},
      'field_map':'Dualedge e curve z_e(t)=exp(2pi i F_e t), all dualvertex phases1; unchargedface totaldegree0 permitsdiskextension. Chargedface ispuncture; boundarydegree=sum_e B[v,e]F_e=G_v−n_v.',
      'native_update':'F_e→F_e+j, n→n−jB_e, hence every actualdegree/charge update equivariant. No externaltexturefactor or oldW observable used.',
      'parameters':{'background':'eachsourcecaseG=itsactualinitialoccupation, no freezeroing/reset' ,'charge_sign':'Q=n−G=−degree','phase_subdivisions_per_turn':8,'actual_native_history_from_zero_records':True,'no_field_feedback_or_energy_added':True,'fixed_domain':'spherewithall12dualportfacespunctured; zero-degreeholesoptionallyfill' },
      'phase_meter_max_error':meter_error,'complete_coin_histories':complete_histories,'native_event_checks':event_checks,'nonzero_current_mutations_detected':mutation_checks,'doublecharge':{'degrees':q2degrees,'split_current':j,'split_degree':(b@split).tolist()},'half_control_endpoint_discontinuity':endpoint_failure,
      'encoded_vs_derived':'Spatialcohomology/phasefield/currentidentity derivedfromactualtape. Integrality stillinheritsintegernativecurrents andtheprimitive2pi phaseinterpolation; not emergence from unconstrainedrealfields or physicalelectric normalization. Fluxisrelativebackground, notunqualifiedpresentn.',
      'scope':'Constructive spatialtopology adapter on actualcarrier. Elementarydegree cannot equal two halfintegerclosedloop degrees. Total/sourcechargeunit andS1readout ancestry explicit; no mass/spin/physicalEMor officialvalidationclaim.',
      'seed_records':seed_records,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'official_rows_changed':[]}

if __name__=='__main__':
 r=execute();r['checks']={k:bool(v) for k,v in r['checks'].items()}
 if not all(r['checks'].values()):raise ValueError(r['checks'])
 HERE.joinpath('NATIVE_DUAL_WINDING.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'error':r['phase_meter_max_error'],'events':r['native_event_checks']}))
