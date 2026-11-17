"""Frozen common-source interface check; no producer imports or new native calls."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import argparse,gzip,hashlib,json,math,time
import numpy as np
HERE=Path(__file__).resolve().parent;I=HERE/'inputs'
def read(p):return json.loads(p.read_text())
def sci(p):return read(p)['scientific']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def step(q,a,b,c):
    x=list(q);s=x[a]+x[b];before=x[a];x[a]=s//2+c*(s%2);x[b]=s-x[a]
    return x,before-x[a]
def pair_step(p,n,a,b,c):
    x={p:1,n:-1};left=x.get(a,0);s=left+x.get(b,0)
    x[a]=s//2+c*(s%2);x[b]=s-x[a];j=left-x[a];x={i:v for i,v in x.items()if v}
    out=0 if not x else 256*next(i for i,v in x.items()if v==1)+next(i for i,v in x.items()if v==-1)
    return out,j
def provenance():
    p=read(HERE/'PROVENANCE.json');manifests={}
    for k,c in p['components'].items():
        path=HERE/c['copied_manifest'];assert sha(path)==c['manifest_sha256']
        manifests[k]={n:h for h,n in(line.split('  ',1)for line in path.read_text().splitlines())}
    for n,m in p['inputs'].items():assert sha(HERE/n)==m['sha256']==manifests[m['component']][m['source_relative_path']]
    return p
def source_join(edges):
    evidence=json.loads(gzip.decompress((I/'loader/RECORDS.json.gz').read_bytes()))
    original={};body_records=0
    for trial in evidence['trials']:
        q=trial['initial_audit_only'][:];read_source=[0]*240
        for st in trial['stages']:
            base=12*st['cell']
            for block in st['blocks']:
                if block['source']is None:continue
                a,b,c,j=block['source']
                if trial['law']=='native':q,k=step(q,a,b,c)
                else:k=c*(q[a]-q[b]);q[a]-=k;q[b]+=k
                assert j==k;body_records+=1
            assert q[base:base+12]==st['decoded'][0]
            for a,b,c,j in st['census']:q,k=step(q,a,b,c);assert k==j;body_records+=1
            read_source[base:base+12]=st['decoded'][1]
        assert q==read_source==trial['observed_source']==trial['source_audit_only']
        assert trial['selected']==(trial['law']=='native')
        if trial['selected']:original[trial['index']]=trial
    assert len(original)==12
    inverse={};vacuum=outside=0
    for code,e,c,j,out in evidence['source_join']['loaded_branches']:
        p,n=divmod(code,256);actual,current=pair_step(p,n,*edges[e],c)
        assert(actual,current)==(out,j)
        key=(e,c,j,out);assert key not in inverse or inverse[key]==code;inverse[key]=code
        vacuum+=out==0
        if out:pp,nn=divmod(out,256);outside+=not(pp<12 and 12<=nn<24)
    assert len(inverse)==181440 and vacuum==2 and outside==70
    flags=evidence['source_join']['product_word_flags'];assert len(flags)==65536 and flags[0]=='collision'
    assert dict(Counter(flags))=={'collision':240,'outside_cells':57216,'accepted':144,'outside_range':7936}
    routes=json.loads(gzip.decompress((I/'route/RECORDS.json.gz').read_bytes()));final={};route_records=0
    for trial in routes['trials']:
        q=original[trial['index']]['observed_source'][:]
        for move in trial['moves']:
            for a,b,c,j in move['records']:
                if trial['law']=='native':q,k=step(q,a,b,c)
                else:k=0
                assert k==j;route_records+=1
        assert q==trial['observed_source']==trial['source_audit_only']
        if trial['law']=='native':
            wanted=[0]*240;wanted[0]=1;wanted[14]=-1
            assert q==wanted and trial['success']and trial['positive_frame']==[0,1]
            assert[(m['source'],m['target'])for m in trial['moves'][-2:]]==[(0,1),(1,0)]
            final[trial['index']]=trial
        else:assert not trial['success']
    assert len(final)==12
    return final,{'replayed_source_body_current_records':body_records,'replayed_route_records':route_records,
        'fine_action_columns':len(inverse),'contact_outputs':vacuum,'outside_cell_outputs':outside,'chart_counts':dict(Counter(flags))}
def operator_bound(n):
    x=F(2*n,10**11);e=x/(1-x);f=n*F(2,10**11)**2/(1-x)**2+e*e
    return(e*(2+e)+f)/2
def region_ops(cell):
    ops=[]
    def controlled(ctrl,zeros):
        for z in zeros:ops.append(('flip',(z,)))
        scratch=list(range(9,13));chain=[(ctrl[0],ctrl[1],scratch[0])]
        for j in range(2,len(ctrl)-1):chain.append((scratch[j-2],ctrl[j],scratch[j-1]))
        ops.extend(('ccx',t)for t in chain);ops.append(('ccx',(scratch[len(ctrl)-3],ctrl[-1],8)))
        ops.extend(('ccx',t)for t in reversed(chain));ops.extend(('flip',(z,))for z in reversed(zeros))
    if cell==0:
        controlled([3,4,5,6,7],[3,4,5,6,7]);controlled([2,3,4,5,6,7],[2,4,5,6,7])
    else:
        controlled([2,3,4,5,6,7],[4,5,6,7]);controlled([3,4,5,6,7],[3,5,6,7])
    return ops
def broad_monitor(cell,routines,pointer):
    ops=region_ops(cell);truth=np.arange(512);initial=truth.copy()
    for kind,t in ops:
        if kind=='flip':truth^=1<<t[0]
        else:truth^=(((truth>>t[0])&1)&((truth>>t[1])&1))<<t[2]
    expected=initial^((((initial&255)>=cell)&((initial&255)<cell+12)).astype(int)<<8)
    assert np.array_equal(truth,expected)
    labels=list(range(cell,cell+12))+([12,14,240,255]if cell==0 else[0,11,240,255])
    state=np.zeros((8192,len(labels)),complex)
    for j,p in enumerate(labels):state[p,j]=1
    for kind,t in ops:
        mask=sum(1<<k for k in t);bases=np.array([i for i in range(8192)if i&mask==0])
        ix=np.array([bases+sum(((k>>j)&1)<<bit for j,bit in enumerate(t))for k in range(1<<len(t))])
        block=state[ix];state[ix]=(routines[kind]@block.reshape(len(ix),-1)).reshape(block.shape)
    ideal=np.zeros_like(state)
    for j,p in enumerate(labels):ideal[p+(256 if cell<=p<cell+12 else 0),j]=1
    residual=float(np.linalg.norm(state-ideal,'fro'));assert residual<1e-8
    blocks=state[:512].reshape(2,256,len(labels))
    scratch=np.asarray([state[512*k:512*(k+1)].conj().T@state[512*k:512*(k+1)]for k in range(16)])
    dirty=scratch[1:].sum(axis=0)
    effects=np.asarray([np.einsum('iac,ij,jad->cd',blocks.conj(),e,blocks,optimize=True)for e in pointer])
    assert np.linalg.norm(effects[1][:12,:12]-np.eye(12),2)<1e-8
    assert np.max(abs(np.diag(effects[1])[12:]))<1e-8
    normres=np.eye(len(labels))-effects.sum(axis=0)-dirty
    n=sum(40 if kind=='ccx'else 2 for kind,t in ops);assert n==(680 if cell==0 else 672)
    return blocks,dirty,effects,{'cell':cell,'labels':labels,'operations':[[k,list(t)]for k,t in ops],
        'routine_applications':len(ops),'represented_native_pulses':n,'all512_clean_truth_columns':True,
        'raw_column_residual':residual,'unresolved_processor_plus_rounding_operator_residual':float(np.linalg.norm(normres,2)),
        'processor_first_failure_probability_upper':str(n*F(2,10**11)**2/(1-F(2*n,10**11))**2),
        'blank_reset_bits_charged':10,'reset_pointer_reads_charged':10,'flag_and_scratch_pointer_reads_charged':5},scratch
def conditional_source(blocks,rho,effect):
    a=blocks[:,:,:12]
    return np.einsum('ij,jac,cd,ibd->ab',effect,a,rho,a.conj(),optimize=True)
def position_read(rho,pointer):
    args=[rho.reshape((2,)*16),list(range(16))]
    for bit in range(8):args.extend([pointer,[16+bit,8+bit,bit]])
    args.append(list(range(16,24)))
    return np.einsum(*args,optimize='greedy').real.reshape(256)
def propagate(Q,initial,dephase=False):
    rho=np.zeros((12,12),complex);rho[initial,initial]=1;deficit=np.eye(12)-Q.conj().T@Q;history=[rho.copy()]
    for _ in range(512):
        rho=Q@rho@Q.conj().T+np.trace(deficit@rho)*np.eye(12)/12
        if dephase:rho=np.diag(np.diag(rho))
        history.append(rho.copy())
    return np.asarray(history)
def phase_interval(records,error):
    eps=float(error);radius=2*math.sqrt(2)*eps/(1-2*eps);angle=radius/(1-radius)
    z=records['1/10']['complex_accepted_fringe'];e=-10*math.atan2(z[1],z[0]);short=[e-10*angle,e+10*angle]
    z=records['4']['complex_accepted_fringe'];base=-math.atan2(z[1],z[0])/4
    possible=[base+k*math.pi/2 for k in range(32)if short[0]-angle/4<=base+k*math.pi/2<=short[1]+angle/4]
    assert len(possible)==1
    return{'short':short,'long':[possible[0]-angle/4,possible[0]+angle/4],'phase_error':angle}
def run():
    origins=provenance();edges=read(I/'loader/graph.json')['seams'];routes,sourceinfo=source_join(edges)
    source=np.load(I/'flow/source_columns.npz');assert np.array_equal(source['seams'],edges)
    L=np.zeros((240,240),int)
    for a,b in edges:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    assert np.array_equal(L,1260*np.eye(240,dtype=int)-source['counts'])
    a=np.load(I/'energy/first_cell_certified/maps.npz');b=np.load(I/'energy/second_cell/maps.npz')
    f=np.load(I/'spatial/inputs/first_cell_maps.npz');field=np.load(I/'spatial/FIELD_A_RECORDS.npz')
    assert np.array_equal(a['certified_Q_proposal'],f['certified_Q_proposal'])
    assert np.array_equal(a['dark_integer_basis'],f['dark_integer_basis'])
    D0=a['dark_integer_basis'];D1=b['dark_integer_basis']
    assert np.array_equal(L@D0,6*D0)and np.array_equal(L@D1,6*D1)
    assert not any(D0[0])and not any(D1[14])and not np.any(D0[12:])and not np.any(D1[:12])
    rho0=propagate(a['certified_Q_proposal'],0);rho1=propagate(b['certified_Q_proposal'],2)
    deph=propagate(a['certified_Q_proposal'],0,True)
    difference=float(np.max(abs(rho0-field['native_density_history'])))
    assert difference<1e-12 and np.max(abs(deph-field['dephasing_density_history']))<1e-12
    joined=sci(I/'energy/routed_join/RECEIPT.json');jdata=np.load(I/'energy/routed_join/maps.npz')
    for row in joined['handoffs']:
        n=row['trial_index'];assert row['actual_source']==routes[n]['observed_source']and row['native_frame']==[0,1]
        assert np.max(abs(jdata['trial_%02d_first'%n]-rho0[-1]))<1e-12
        assert np.max(abs(jdata['trial_%02d_second'%n]-rho1[-1]))<1e-12
        assert row['route_record_sha256']==hashlib.sha256(json.dumps(routes[n],sort_keys=True).encode()).hexdigest()
    # Independently recover the separated-pair source generator, not2*input6.
    paircolumns=0
    for u in D0.T:
        for v in D1.T:
            out=defaultdict(int);target={}
            for p in np.flatnonzero(u):
                for n in np.flatnonzero(v):
                    coeff=int(u[p]*v[n]);target[256*int(p)+int(n)]=1248*coeff
                    for aa,bb in edges:
                        for coin in(0,1):out[pair_step(int(p),int(n),aa,bb,coin)[0]]+=coeff;paircolumns+=1
            assert{k:v for k,v in out.items()if v}==target
    clock=read(I/'clock/RECEIPT.json');assert set(clock['native_pair_generator_values'])=={12}
    assert set(clock['native_single_generator_values'])=={6}
    reference=np.load(I/'flow/attempts/stage2_002/evidence/reference_preparation.npz');pointer=reference['pointer_effects_C']
    assert np.linalg.norm(pointer.sum(axis=0)-np.eye(2),2)<1e-10
    app=np.load(I/'energy/apparatus/maps.npz');routines={'ccx':app['native_routine_ccx'],'flip':app['native_routine_flip']}
    mon=[broad_monitor(c,routines,pointer)for c in(0,12)]
    reset=F(35,2**31);ptr=F(200000000001,20000000000000000000000)
    broad=sum((operator_bound(m[3]['represented_native_pulses'])+10*(reset+ptr)+5*ptr for m in mon),F())
    firstrec=sci(I/'energy/first_cell_certified/RECEIPT.json');secondrec=sci(I/'energy/second_cell/RECEIPT.json')
    sharp=read(I/'energy/SHARP_FLOW_CERTIFICATE.json')
    assert sharp==read(I/'spatial/inputs/SHARP_FLOW_CERTIFICATE.json')
    transfer=read(I/'spatial/TRANSFER_RESULT.json');flow=F(sharp['512_complete_U4_halfdiamond_upper'])
    assert F(transfer['flow_complete_halfdiamond_upper'])*512==flow
    assert F(transfer['apparatus512_upper'])==F(firstrec['512epoch_apparatus_error_upper'])
    raw=[flow+F(r['512epoch_apparatus_error_upper'])for r in(firstrec,secondrec)]
    gain=[x/(1-2*x)for x in raw];bad=[F(r['512epoch_ideal_bad_weight_upper'])for r in(firstrec,secondrec)]
    route_error=F(read(I/'route/RECEIPT.json')['per_trial_timeout_upper']);terminal=8*(reset+ptr)
    S=sum(gain)+sum(bad)+route_error+broad+terminal;field_error=S/(1-2*S);field_upper=2*field_error
    fa=read(I/'spatial/FIELD_A_RESULT.json');comparison=F(fa['comparison_512_error_upper'])
    delta=F(secondrec['certified_subblock_error_upper']);kap=4*delta+2*delta*delta;u=F(1,2**53);gamma=128*u/(1-128*u)
    secondcomp=512*kap/(1-512*kap)+512*12*4096*gamma/(1-512*kap)+F(1,10**12)
    controlraw=raw[0]+512*terminal;controlgain=controlraw/(1-2*controlraw)
    C=controlgain+gain[1]+route_error+broad+terminal;control_error=C/(1-2*C)+comparison+secondcomp
    control_lower=F(fa['dephasing_coordinate_halfturn_TV'])-2*control_error
    ry=fa['exact_rotations']['Ry'];assert ry[0]==1 and ry[1]==0
    spatiallaws={};arrays={}
    for name,first in [('native',rho0[-1]),('dephasing',deph[-1])]:
        probabilities=[];flags=[];dirty=[];scratch_words=[]
        for k,density in enumerate((first,rho1[-1])):
            blocks,dirt,effects,meta,scratch=mon[k]
            states=[conditional_source(blocks,density,e)for e in pointer]
            flags.append([float(np.trace(x).real)for x in states]);dirty.append(float(np.trace(dirt[:12,:12]@density).real))
            scratch_words.append([float(np.trace(g[:12,:12]@density).real)for g in scratch])
            if k==0:probabilities=position_read(states[1],pointer)
        positions=probabilities*flags[1][1]
        rejects=[flags[0][aa]*flags[1][bb]for bb in(0,1)for aa in(0,1)if not(aa==bb==1)]
        dirtyfirst=dirty[0];dirtysecond=sum(flags[0])*dirty[1]
        record={'all256_position_centres':positions.tolist(),'region_reject_flag_centres':rejects,
            'dirty_scratch_first':dirtyfirst,'dirty_scratch_second_after_first_clean':dirtysecond,
            'scratch_word_comparison_centres':scratch_words,
            'unresolved_processor_plus_rounding_residual':float(1-sum(positions)-sum(rejects)-dirtyfirst-dirtysecond)}
        assert abs(record['unresolved_processor_plus_rounding_residual'])<1e-7
        tv=float(np.sum(abs(positions[:12]-positions[ry]))/2);record['Ry_TV_centre']=tv
        if name=='native':assert tv<float(field_upper)
        else:assert tv>.01
        spatiallaws[name]=record
    # The ideal regional loader is the identity on both native/control outputs.
    # Its finite error is added without conditioning away rejected outcomes.
    oldlinear=F(joined['linear_physical_TP_composition_bound']);energyS=oldlinear+broad
    energy_error=energyS/(1-2*energyS);coupled=sci(I/'energy/coupled_energy/RECEIPT.json')
    singleS=F(coupled['complete_source_preparation_error_upper'][0])+F(coupled['complete_coupled_sector_reader_error_upper'])+route_error+broad
    single_error=singleS/(1-2*singleS)
    records=joined['handoffs'][0]['complete_energy_output_centres'];bounds={}
    for name,error in [('single',single_error),('pair',energy_error)]:bounds[name]=phase_interval({t:records[t][name]for t in records},error)
    lo,hi=bounds['single']['long'];pl,ph=bounds['pair']['long'];ratio=[pl/hi,ph/lo]
    oldratio=math.log1p(-12/1260)/math.log1p(-6/1260)
    gates={'all12_same_source_handoffs':len(routes)==12,'loaded_current_action_and_partial_injection':sourceinfo['fine_action_columns']==181440,
        'same_graph_modes_seed_and_density':difference<1e-12,'both_full_region_predicates':all(m[3]['all512_clean_truth_columns']for m in mon),
        'complete_energy_error_below_001':energy_error<F(1,1000),'energy_interval_and_oldclock_separation':ratio[0]<2<ratio[1]<oldratio,
        'full_spatial_halfturn_below_004':field_upper<F(1,250),'active_dephasing_separation_above_01':control_lower>F(1,100),
        'shape_retained':F(1,2)-field_error>F(1,3)and F(1,3)-2*field_error>F(1,5)}
    for k,m in enumerate(mon):arrays['loader%d_clean_columns'%k]=m[0];arrays['loader%d_dirty_gram'%k]=m[1];arrays['loader%d_effects'%k]=m[2];arrays['loader%d_scratch_word_grams'%k]=m[4]
    arrays.update(first_final_density=rho0[-1],second_final_density=rho1[-1],dephasing_final_density=deph[-1])
    result={'candidate':'A-old240-selected-compact-source-v1','all_gates_pass':all(gates.values()),'gates':gates,
        'component_manifest_sha256':{k:v['manifest_sha256']for k,v in origins['components'].items()},
        'source_interface':sourceinfo,'independent_pair_mode_columns':paircolumns,'field_energy_density_max_difference':difference,
        'literal_new_native_calls':0,'fresh_fullflow_executions':0,'compressed_native_routine_applications':sum(m[3]['routine_applications']for m in mon),
        'new_broad_loader_programs':[m[3]for m in mon],
        'processor_failure_semantics':'Positive first-failure complements inherited from native word/prefix theorem; not evaluated as scalar probabilities from success-map numerical deficits.',
        'scratch_read_semantics':'All16 four-bit scratch comparison outcomes retained; four native bit-reader/reset allowances per monitor bound their physical readout.','complete_broad_loader_error_upper':str(broad),
        'spatial_complete_error_upper':str(field_error),'spatial_halfturn_TV_upper':str(field_upper),
        'active_dephasing_TV_lower':str(control_lower),'isotropic_TV_lower':str(F(1,2)-field_error),'Rx_TV_lower':str(F(1,3)-2*field_error),
        'spatial_full_outcome_centres':spatiallaws,'energy_complete_error_upper':str(energy_error),
        'energy_intervals':bounds,'pair_to_single_interval':ratio,'pair_keV_at_single511_anchor_interval':[511*x for x in ratio],
        'old_fixed_count_ratio_retained':oldratio,'old_changed_clock_long_fringe_retained':coupled['fixed_count_phase_nonadditivity'],
        'charge_current_Gauss_transfer':'Exact theorem on admitted finite pair domain; no long trajectory or TV transfer of unbounded reward claimed.',
        'broad_loader_before_assay_choice':True,'spatial_assay_before_dark_sector_diagnosis':True,
        'fresh_copies_same_preparation':True,'finite_exact_halfturn_claimed':False,'other_group_G6':False,'upstream_award_or_pin_changed':False}
    return result,arrays
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args();started=time.monotonic()
    reg=read(HERE/'PREREG.json')
    for n,h in reg['sha256'].items():assert sha(HERE/n)==h,n
    result,arrays=run();result['prereg_sha256']=sha(HERE/'PREREG.json')
    if args.record:
        with(HERE/'RECEIPT.json').open('xb')as f:f.write(encode(result))
        np.savez_compressed(HERE/'INTERFACE_MAPS.npz',**arrays)
    else:
        assert read(HERE/'RECEIPT.json')==result
        old=np.load(HERE/'INTERFACE_MAPS.npz')
        for k,v in arrays.items():assert np.array_equal(old[k],v),k
    print(json.dumps({'all_gates_pass':result['all_gates_pass'],'broad_loader_error':float(F(result['complete_broad_loader_error_upper'])),
        'field_halfturn_upper':float(F(result['spatial_halfturn_TV_upper'])),'dephasing_lower':float(F(result['active_dephasing_TV_lower'])),
        'energy_complete_error':float(F(result['energy_complete_error_upper'])),'pair_keV_interval':result['pair_keV_at_single511_anchor_interval'],
        'same_density_difference':result['field_energy_density_max_difference'],'seconds':time.monotonic()-started},indent=2))
    assert result['all_gates_pass']
if __name__=='__main__':main()
