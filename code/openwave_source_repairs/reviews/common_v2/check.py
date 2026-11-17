"""Independent finite common-source interface audit; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,math,sys,time
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gain(x):return x/(1-2*x)
def ceil15(x):return F(-((-x.numerator*10**15)//x.denominator),10**15)
def pump(Q,seed,dephase=False):
    escape=np.eye(12)-Q.conj().T@Q
    M=np.kron(Q,Q.conj())+np.outer(np.eye(12).reshape(144)/12,escape.T.reshape(144))
    if dephase:M*=np.eye(12).reshape(144)[:,None]
    initial=np.zeros(144,complex);initial[13*seed]=1
    return (np.linalg.matrix_power(M,512)@initial).reshape(12,12)
def flag_state(columns,rho,E):
    A=columns[:,:,:12].reshape(512,12)
    block=(A@rho@A.conj().T).reshape(2,256,2,256)
    return np.einsum('ab,biaj->ij',E,block)
def addresses(rho,effects):
    law=[]
    for word in range(256):
        E=np.ones((1,1),complex)
        for bit in range(7,-1,-1):E=np.kron(E,effects[(word>>bit)&1])
        law.append(float(np.sum(E.T*rho).real))
    return np.asarray(law)
def exact_R(d):
    a,b,c=int(d[:,0]@d[:,0]),int(d[:,0]@d[:,1]),int(d[:,1]@d[:,1]);det=a*c-b*b
    return [[F(c*int(d[i,0])*int(d[j,0])-b*(int(d[i,0])*int(d[j,1])+int(d[i,1])*int(d[j,0]))+a*int(d[i,1])*int(d[j,1]),det)for j in range(12)]for i in range(12)]
def exact_chart():
    n,p,z,f,g=(-1,0),(1,0),(0,0),(0,1),(0,-1)
    body=[(n,f,z),(p,f,z),(n,g,z),(p,g,z),(z,n,f),(z,p,f),(z,n,g),(z,p,g),(f,z,n),(f,z,p),(g,z,n),(g,z,p)]
    points=[]
    for cell in range(20):
        for x,y,z in body:
            if cell==1:y,z=tuple(-v for v in z),y
            points.append((x,(y[0]+10*cell,y[1]),z))
    lookup={p:i for i,p in enumerate(points)};assert len(lookup)==240
    p=[lookup[(tuple(-v for v in x),y,tuple(-v for v in z))]for x,y,z in points]
    assert all(p[p[i]]==i for i in range(240)) and p[0]==1 and p[1]==0
    return points,np.asarray(p+list(range(240,256)))
def extract_record(c):
    nodes=range(c['cell']*12,c['cell']*12+12);change={v:0 for v in nodes};eq=[];t=c['tape'];i=0
    while i<len(t):
        a,b,coin,j=t[i];delta=change[a]-change[b];first=j;oldcoin=coin
        change[a]-=j;change[b]+=j;i+=1
        while i<len(t):
            x,y,coin,j=t[i];assert(x,y)==(a,b);change[a]-=j;change[b]+=j;i+=1
            if coin!=oldcoin:eq.append((a,b,2*first+j-delta));break
        else:raise AssertionError('incomplete current contrast')
    initial={(0,14)[c['cell']]:0}
    while len(initial)<12:
        old=len(initial)
        for a,b,d in eq:
            if a in initial:initial[b]=initial[a]-d
            if b in initial:initial[a]=initial[b]+d
        assert len(initial)>old
    assert [initial[v]for v in nodes]==c['decoded']['initial']
    assert [initial[v]+change[v]for v in nodes]==c['decoded']['terminal']
    return [initial[v]for v in nodes]
def native_record(initial,event):
    q=initial[:];a,b,coin,j=event['record'];d=q[a]-q[b]
    if event['law']=='native':
        expected=(d-(2*coin-1)*(abs(d)%2))//2
    else:
        assert event['law']=='swap_hold';expected=d if coin else 0
    assert j==expected;q[a]-=j;q[b]+=j
    return q
def calibration(I,result):
    routes=read(I/'calibration/ROUTE_RECORDS.json')
    raw=(I/'calibration/inputs/bridge_RECORDS.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==read(I/'calibration/inputs/bridge_RECEIPT.json')['evidence_sha256']
    source={t['index']:t for t in json.loads(gzip.decompress(raw))['trials']if t['selection_law']=='native'}
    prefixes=[];events=0
    def digest(v):return hashlib.sha256((json.dumps(v,sort_keys=True,separators=(',',':'))+'\n').encode()).hexdigest()
    for route in routes['trials']:
        trial=source[route['bridge_trial']];q=trial['initial_audit_only'][:]
        assert sum(q)==1 and not any(q[12:])
        records=[row for e in trial['selection']for row in e['reference']+[e['body']]]+trial['census']
        for row in records:q=native_record(q,{'law':'native','record':row});events+=1
        h=trial['basis_handoff'];assert q==h['observed_source']==h['actual_source_audit_only']==route['input']
        assert sum(abs(x)for x in q)==1 and q[h['occupied']]==1
        assert digest({'selection':trial['selection'],'census':trial['census']})==h['prefix_sha256']
        for leg in route['legs']:
            for row in leg['records']:q=native_record(q,{'law':'native','record':row});events+=1
        assert q==[1]+[0]*239==route['observed_final_source']==route['actual_final_source_audit_only']
        assert route['frame']==[0,1]and route['frame_native_event'][3]==1
        prefix=digest({'bridge_prefix':h['prefix_sha256'],'route':route['legs']});assert prefix==route['prefix_sha256'];prefixes.append(prefix)
    q=[1]+[0]*239
    for leg in routes['forced_timeout']['legs']:
        for row in leg['records']:q=native_record(q,{'law':'native','record':row});events+=1
    assert q==routes['forced_timeout']['final_source']and routes['forced_timeout']['status']=='route_timeout'
    assert len(prefixes)==36 and prefixes==result['actual_source_prefixes']
    assert F(result['route_failure_upper'])==1-(1-F(1,2**32))**5<1-(1-F(1,2**32))**8
    assert result['independent_unit_chart_not_signed_pair_spectator']
    return events
def pointer_record(initial,final,matrices,effects):
    n=initial.shape[-1];prep=matrices[0];a=prep[:2,0]
    branches=np.stack((a[0]*initial.reshape(-1),a[1]*final.reshape(-1)),axis=1)/math.sqrt(n)
    rho=branches.T@branches.conj();out={}
    for setting,order in [('X',[matrices[1]]),('Y',[matrices[2],matrices[1]])]:
        state=rho.copy();failure=0
        for M in order:
            failure+=float(np.trace(M[2:,:2]@state@M[2:,:2].conj().T).real)
            state=M[:2,:2]@state@M[:2,:2].conj().T
        law=[float(np.trace(E@state).real)for E in effects]
        out[setting]=(law,failure,float(np.trace(rho).real-failure-sum(law)))
    return out,float(np.vdot(prep[2:,0],prep[2:,0]).real)
def approximate_decode(short,long,epsilon):
    # Floating independent comparison only; producer rational endpoints carry the certificate.
    r=2*math.sqrt(2)*float(epsilon)/(1-float(epsilon));delta=r/(1-r)
    theta=-math.atan2(short[1],short[0]);lo=(theta-delta)*10;hi=(theta+delta)*10
    theta=-math.atan2(long[1],long[0])
    low=math.ceil((4*lo-theta-delta)/(2*math.pi));high=math.floor((4*hi-theta+delta)/(2*math.pi))
    assert low==high
    return [(theta+2*math.pi*low-delta)/4,(theta+2*math.pi*low+delta)/4],low
def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--record',action='store_true');args=ap.parse_args()
    start=time.perf_counter();P=args.source.resolve();I=P/'inputs';r=read(P/'RECEIPT.json')
    assert r['all_gates_pass'] and all(r['gates'].values())
    frozen=read(HERE/'PREREG.json');assert sha(P/'PREREG.json')==frozen['producer_prereg_sha256']
    for name,h in read(P/'PREREG.json')['files'].items():assert sha(P/name)==h,name
    for name,item in read(P/'BRIDGE_ORIGINS.json')['inputs'].items():assert sha(P/name)==item['sha256']
    saved=np.load(P/'COMMON_MAPS.npz');origins=read(P/'BRIDGE_ORIGINS.json')
    for component,item in origins['components'].items():assert sha(I/component/'MANIFEST.sha256')==item['manifest_sha256']
    records=json.loads(gzip.decompress((I/'extraction/RECORDS.json.gz').read_bytes()))
    handoffs=read(I/'extraction/HANDOFFS.json')['handoffs'];prefixes=[];events=0;native=0;control=0;zero=0;reverse=0;censuses=0
    edges={frozenset(x)for x in read(I/'extraction/graph.json')['seams']}
    for trial in records['trials']:
        initial=[0]*240;tape=[]
        for c in trial['censuses']:initial[c['cell']*12:c['cell']*12+12]=extract_record(c);tape.extend(c['tape']);censuses+=1
        assert tape==[e['record']for e in trial['events'][:len(tape)]]
        assert all(e['law']=='native'for e in trial['events'][:len(tape)])
        # Only after reconstruction, compare separately labelled harness data.
        assert initial==trial['initial_audit_only'];K=sum(initial[:12])
        assert abs(K)<=12 and sum(initial[12:24])==-K and initial[0]==initial[14]==0 and max(map(abs,initial))<=64
        q=initial
        for e in trial['events']:
            assert frozenset(e['record'][:2])in edges;q=native_record(q,e);events+=1
        assert q==trial['observed_source']==trial['source_audit_only']
        if trial['law']=='native'and K:
            assert q==[1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]+[0]*225
            digest=hashlib.sha256(json.dumps(trial['events'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
            matches=[x for x in handoffs if x['complete_prefix_sha256']==digest];assert len(matches)==1
            h=matches[0];assert h['actual_source']==h['observed_source']==q and h['positive_frame']==[0,1]and h['pair_code']==14
            assert not h['numerical_total_supplied_to_observer'] and not h['source_reset']
            assert [e['record'][:2]for e in trial['events']if e['kind']=='frame_move'and e['record'][3]]==[[0,1],[1,0]]
            if K<0:assert any(e['kind']=='orientation_move'and e['record'][3]for e in trial['events']);reverse+=1
            prefixes.append(digest);native+=1
        elif trial['law']=='native':assert not any(q) and trial['outcome']=='no_unit';zero+=1
        else:assert trial['outcome']=='selection_rejected';control+=1
    assert(native,reverse,zero,control)==(24,12,1,25)and prefixes==r['source_join']['all24_actual_prefixes']
    calibration_events=calibration(I,r['single_calibration_join'])
    eta=F(22,2**63)+F(24*64**2//2*11+22+11*237+3*239+2,2**64);assert eta==F(8501,2**58)==F(r['source_join']['inconclusive_upper'])
    first=np.load(I/'common_v1/inputs/energy/first_cell_certified/maps.npz');second=np.load(I/'common_v1/inputs/energy/second_cell/maps.npz')
    maps=np.load(I/'common_v1/INTERFACE_MAPS.npz');ec=read(I/'energy/inputs/certificate.json');basis=np.asarray(ec['integer_modes'],int)
    assert np.array_equal(first['dark_integer_basis'],basis[:2].T)and np.array_equal(second['dark_integer_basis'],basis[2:].T)
    a=np.column_stack((basis[0]/2,(2*basis[1]-basis[0])/math.sqrt(12)));b=np.column_stack((basis[2]/2,(2*basis[3]-basis[2])/math.sqrt(12)))
    pair=np.stack([np.outer(a[:,i],b[:,j])for i in range(2)for j in range(2)],axis=-1)
    rho0=pump(first['certified_Q_proposal'],0);rho1=pump(second['certified_Q_proposal'],2);deph=pump(first['certified_Q_proposal'],0,True)
    density_error=float(np.max(abs(np.kron(rho0,rho1)-saved['same_prepared_full144_density'])));assert density_error<1e-12
    target=pair[:12,12:24].reshape(144,4);target=target@target.conj().T/4
    assert np.max(abs(target-saved['same_selected_full144_density']))<1e-14
    chart,perm=exact_chart();assert np.array_equal(perm,saved['global_Ry_address_permutation'])
    projectors=[exact_R(basis[:2,:12].T),exact_R(basis[2:,12:24].T)];p0=perm[:12];p1=perm[12:24]-12
    for i in range(12):
      for j in range(12):
       for k in range(12):
        for l in range(12):assert projectors[0][i][k]*projectors[1][j][l]==projectors[0][p0[i]][p0[k]]*projectors[1][p1[j]][p1[l]]
    pointer=np.load(I/'common_v1/inputs/flow/attempts/stage2_002/evidence/reference_preparation.npz')['pointer_effects_C']
    law_error=0
    for name,states in [('native',(rho0,rho1)),('dephasing',(deph,rho1))]:
        positions=[];flags=[];dirty=[];scratch=[]
        for i,state in enumerate(states):
            outputs=[flag_state(maps[f'loader{i}_clean_columns'],state,E)for E in pointer]
            flags.append([float(np.trace(x).real)for x in outputs]);positions.append(addresses(outputs[1],pointer))
            dirty.append(float(np.trace(maps[f'loader{i}_dirty_gram'][:12,:12]@state).real))
            scratch.append([float(np.trace(g[:12,:12]@state).real)for g in maps[f'loader{i}_scratch_word_grams']])
        law=np.multiply.outer(*positions);law_error=max(law_error,float(np.max(abs(law-saved[name+'_all65536_position_centres']))));assert law_error<1e-12
        row=r['full_joint_spatial_records'][name];rejects=[flags[0][x]*flags[1][y]for x,y in [(0,0),(0,1),(1,0)]]
        residual=1-law.sum()-sum(rejects)-dirty[0]-sum(flags[0])*dirty[1]
        assert abs(residual-row['unresolved_processor_plus_rounding_residual'])<1e-12
        assert np.max(abs(np.asarray(scratch)-row['scratch_words_comparison_centres']))<1e-12
        assert np.max(abs(np.asarray(rejects)-row['region_reject_centres']))<1e-12
    spatial=read(I/'spatial/RECEIPT.json')
    for oldkey,newkey in [('full_joint_complete_error_upper','spatial_complete_error_upper'),('complete_control_error_upper','complete_control_error_upper')]:
        old=F(spatial[oldkey]);assert ceil15(gain(old/(1+2*old)+eta))==F(r[newkey])
    assert 2*F(r['spatial_complete_error_upper'])==F(r['full_joint_halfturn_TV_upper'])
    lower=F(spatial['active_dephasing_full_joint_TV_lower'])-2*(F(r['complete_control_error_upper'])-F(spatial['complete_control_error_upper']))
    assert -ceil15(-lower)==F(r['active_dephasing_full_joint_TV_lower'])
    flow=read(I/'energy/attempts/v2/RECEIPT.json');observer=read(I/'energy/observer/RECEIPT.json');clock=read(I/'energy/clock/RECEIPT.json')
    assert all(x['all_gates_pass']and all(x['gates'].values())for x in (flow,observer,clock))
    assert observer['matrix_keys']==['Y/'+(-math.pi/4).hex(),'Y/'+(math.pi/4).hex(),'Z/'+(math.pi/4).hex()]
    mats=np.load(I/'energy/observer/NATIVE_TERMINAL_MATRICES.npz');native=[mats[str(i)]for i in range(3)]
    effects=np.load(I/'energy/inputs/native_pointer.npz')['effects'];outputs=np.load(I/'energy/attempts/v2/SOURCE_OUTPUTS.npz')
    pointer_error=0;decoded={}
    for name,initial in [('single',a),('pair',pair)]:
        parts={k:F(v)for k,v in flow['complete_budgets'][name].items()if k!='complete'}
        assert gain(sum(parts.values(),F()))==F(flow['complete_budgets'][name]['complete'])
        eps=ceil15(gain(sum(parts.values(),F())+(eta if name=='pair'else F())));assert eps==F(r['energy_complete_errors_upper'][name])
        observed={}
        for t in ('0.1','4'):
            record,prepfail=pointer_record(initial,outputs[name+'_'+t.replace('.','p')],native,effects);fringe=[]
            for setting in ('X','Y'):
                law,fail,residual=record[setting];row=r['energy_records'][t][name][setting]
                pointer_error=max(pointer_error,max(abs(x-y)for x,y in zip(law,row['pointer_output_centres'])))
                assert max(abs(fail-row['analysis_processor_failure_centre']),abs(prepfail-row['reference_preparation_failure_centre']),abs(residual-row['arithmetic_mass_residual']))<1e-12
                fringe.append((law[0]-law[1])/sum(law))
            observed[t]=fringe
        interval,winding=approximate_decode(observed['0.1'],observed['4'],eps);decoded[name]=interval
        assert winding==r['energy_decoded'][name]['winding_from_short']
        assert max(abs(x-y)for x,y in zip(interval,r['energy_decoded'][name]['energy_interval']))<1e-11
    assert pointer_error<1e-12
    C=math.lcm(1260,1280)*math.ceil(F(10*max(x['total_pulses']for x in flow['words'].values()),math.lcm(1260,1280)))
    assert C==r['common_clock']['slots_per_time_unit']==clock['clock_slots_per_time_unit']
    for row in clock['programmes']:
        N=C*F(row['duration']);assert N==row['total_clock_slots']<2**40
        assert N-row['represented_active_processor_pulses']==row['represented_idle_slots']>=0
    result=dict(producer_prereg_sha256=sha(P/'PREREG.json'),producer_receipt_sha256=sha(P/'RECEIPT.json'),actual_handoffs=24,physical_sign_reversals=12,main_source_events_reconstructed=events,record_only_censuses=censuses,separate_positive_only_calibrations=36,calibration_events_reconstructed=calibration_events,exact_joint_density_entries=20736,full_joint_bins_per_assay=65536,independent_density_error=density_error,independent_full_instrument_error=law_error,independent_native_pointer_error=pointer_error,extraction_inconclusive_upper=str(eta),joint_spatial_TV_upper=r['full_joint_halfturn_TV_upper'],active_control_lower=r['active_dephasing_full_joint_TV_lower'],pair_keV_certified=r['pair_keV_at511_anchor'],floating_decoder_compared_not_new_interval_certificate=True,clock_rate=C,all_complete_error_compositions_checked=True,producer_imports=False,literal_new_native_calls=0,other_group_G6=False,seconds=time.perf_counter()-start)
    if args.record:
        with(HERE/'RESULT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        old=read(HERE/'RESULT.json');assert {k:v for k,v in old.items()if k!='seconds'}=={k:v for k,v in result.items()if k!='seconds'}
    print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
