"""Independent assembly arithmetic: standard library + NumPy, no producer imports."""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
import gzip,hashlib,json,math,time,argparse
import numpy as np

HERE=Path(__file__).resolve().parent
ROW=HERE.parents[1];WORK=ROW/'work'
ENERGY=WORK/'2026-11-15_codex_source_repair_v2'
LOADER=WORK/'2026-11-15_codex_signed_pair_loading'
CLOCK=WORK/'2026-11-15_codex_poisson_source_clock'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def scientific(p):return read(p)['scientific']
def verify_manifest(p):
    for line in (p/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);assert sha(p/n)==h,(p,n)
    return sha(p/'MANIFEST.sha256')
def pair_step(p,n,a,b,c):
    x={p:1,n:-1};s=x.get(a,0)+x.get(b,0)
    x[a]=s//2+c*(s%2);x[b]=s-x[a];x={i:v for i,v in x.items() if v}
    if not x:return 0
    assert sorted(x.values())==[-1,1]
    return 256*next(i for i,v in x.items() if v==1)+next(i for i,v in x.items() if v==-1)

def exact_source():
    edges=read(CLOCK/'graph.json')['seams'];assert len(edges)==630
    L=np.zeros((240,240),dtype=np.int64)
    for a,b in edges:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    vec=[]
    for places,values in [([4,5,6,7],[1,-1,-1,1]),([5,6,8,11],[-1,-1,1,1]),
                          ([12,15,17,18],[-1,-1,1,1]),([12,15,16,19],[-1,-1,1,1])]:
        x=np.zeros(240,dtype=np.int64);x[places]=values;vec.append(x)
        assert np.array_equal(L@x,6*x)
    columns=0
    for u in vec[:2]:
        for v in vec[2:]:
            result=defaultdict(int);target={}
            for p in np.flatnonzero(u):
                for n in np.flatnonzero(v):
                    coeff=int(u[p]*v[n]);target[256*int(p)+int(n)]=1248*coeff
                    for a,b in edges:
                        for coin in (0,1):
                            result[pair_step(int(p),int(n),a,b,coin)]+=coeff;columns+=1
            assert {i:x for i,x in result.items() if x}==target
            assert result[0]==0
    clock=read(CLOCK/'RECEIPT.json')
    assert set(map(F,clock['native_single_generator_values']))=={F(6)}
    assert set(map(F,clock['native_pair_generator_values']))=={F(12)}
    assert F(clock['clock_tick_duration'])==F(1,1260)
    sampler=F(68,2**128)+F(68,math.factorial(160))+F(34,33*math.factorial(33))
    address=F(197,512)**64
    assert sampler==F(clock['per128bit_clock_sampler_failure_upper'])
    assert address==F(clock['per_source_address_failure_upper'])
    assert F(1254,1260)**2-F(1248,1260)==F(1,44100)
    return L,vec,5040*sampler+161280*address,columns

def complete_bound(pulses,target=F()):
    x=F(2*pulses,10**11);e=(x+target)/(1-x)
    failure=F(pulses)*(F(2,10**11)**2)/(1-x)**2
    return (e*(2+e)+failure)/2

def certified_cells(L,vec):
    eigen,V=np.linalg.eigh(L.astype(float))
    U=(V*np.exp(4j*1260*np.log1p(-eigen/1260)))@V.T
    reset=F(35,2147483648);pointer=F(200000000001,20000000000000000000000)
    cells=[]
    for cell,name,basis,threshold in [(0,'first_cell_certified',vec[:2],F(1,4)),(12,'second_cell',vec[2:],F(1,16))]:
        receipt=scientific(ENERGY/name/'RECEIPT.json');a=np.load(ENERGY/name/'maps.npz')
        Q=a['certified_Q_proposal'];D=np.stack(basis,axis=1)[cell:cell+12]
        assert np.array_equal(a['dark_integer_basis'],np.stack(basis,axis=1))
        P=D@np.linalg.inv(D.T@D)@D.T
        delta=F(receipt['certified_subblock_error_upper'])
        assert np.linalg.norm(Q-U[cell:cell+12,cell:cell+12],2)<float(delta)
        defect=np.eye(12)-Q.conj().T@Q;tail=np.eye(12)-P
        for _ in range(64):tail=Q.conj().T@tail@Q+np.trace(tail)/12*defect
        assert np.linalg.norm(tail-a['tail64'],2)<1e-12
        kappa=4*delta+2*delta*delta;u=F(1,2**53);gamma=128*u/(1-128*u)
        rounding=64*12*4096*gamma/(1-64*kappa)+F(1,10**12)
        bound=F(math.ceil(float(np.linalg.norm(tail,'fro'))*1e12),10**12)+rounding+64*kappa/(1-64*kappa)
        assert bound==F(receipt['64epoch_tail_operator_upper']) and bound<threshold
        assert F(receipt['512epoch_ideal_bad_weight_upper'])==threshold**8
        app=21*(reset+pointer)+complete_bound(receipt['region_word_pulses'])+pointer+8*F(2,10**11)+F(1,4)**16
        assert 512*app==F(receipt['512epoch_apparatus_error_upper'])
        cells.append({'Q':Q,'P':P,'D':D,'defect':defect,'bad_bound':threshold**8,'app_bound':512*app})
    return cells

def propagated(cell,address):
    rho=np.zeros((12,12),complex);rho[address,address]=1
    Q=cell['Q'];defect=cell['defect']
    for _ in range(512):
        loss=np.trace(defect@rho)
        rho=Q@rho@Q.conj().T+loss*np.eye(12)/12
    assert abs(np.trace(rho)-1)<1e-8
    assert np.min(np.linalg.eigvalsh((rho+rho.conj().T)/2))>-1e-10
    assert np.trace((np.eye(12)-cell['P'])@rho).real<float(cell['bad_bound'])+1e-7
    return rho

def energy_interval(record,error,t):
    z=record['complex_accepted_fringe'];assert z is not None
    s=min(record[k]['energy_success_probability_centre'] for k in ('X','Y'))
    # TV(p,q)<=error, q accepts surely: conditional TV<=error/s.
    # A bounded X/Y mean changes by at most twice that conditional TV.
    radius=math.nextafter(2*math.sqrt(2)*float(error)/s,math.inf)
    assert radius<.5
    angle=math.nextafter(math.asin(radius),math.inf)
    decoded=-math.atan2(z[1],z[0])/t
    width=angle/t+1e-8
    return [decoded-width,decoded+width],radius

def verify_law(record,expected_length,upstream=None):
    for setting in ('X','Y'):
        r=record[setting];p=np.asarray(r['diagnostic_probability_centres'])
        assert len(p)==expected_length and min(p)>-1e-10
        missing=r.get('upstream_sector_processor_failure_centre',upstream)
        assert missing is not None,'upstream sector failure must be explicit'
        total=sum(p)+r['downstream_processor_failure_centre']+missing
        assert abs(total-1)<1e-8,total
        assert abs(1-total-r['arithmetic_normalization_residual'])<1e-12
        assert r['downstream_processor_failure_centre']>=-1e-10 and missing>=-1e-10
    return True

def main():
    started=time.monotonic();bindings={}
    for directory in (CLOCK,LOADER,WORK/'2026-11-16_codex_signed_pair_routing'):
        bindings[str(directory.relative_to(ROW))+'/MANIFEST.sha256']=verify_manifest(directory)
    reg=read(ENERGY/'PREREG_COUPLED_ENERGY.json')
    for group,directory in [('files',ENERGY),('clock_files',CLOCK),('loading_files',LOADER)]:
        for name,h in reg[group].items():assert sha(directory/name)==h;bindings[str((directory/name).relative_to(ROW))]=h
    receipt=read(ENERGY/'coupled_energy/RECEIPT.json');r=receipt['scientific']
    assert receipt['prereg_sha256']==sha(ENERGY/'PREREG_COUPLED_ENERGY.json')
    for n,h in receipt['artifacts'].items():assert sha(ENERGY/'coupled_energy'/n)==h
    L,vec,clock_error,source_columns=exact_source();cells=certified_cells(L,vec)
    sector=np.load(ENERGY/'poisson_energy/source_sector_maps.npz')
    for name,cell in zip(('first','second'),cells):
        E=sector[name+'_effects'];failure=sector[name+'_processor_failure']
        assert np.linalg.norm(E.sum(axis=0)+failure-np.eye(12),2)<1e-8
        assert np.linalg.norm(E[1]-cell['P'],2)<1e-7
        assert all(np.linalg.eigvalsh((e+e.conj().T)/2).min()>-1e-9 for e in E)
    reset=F(35,2147483648);pointer=F(200000000001,20000000000000000000000)
    lower=scientific(ENERGY/'poisson_energy/RECEIPT.json')
    pulses=max(x['word_pulses']for x in lower['coarse_sector_instruments'])
    base=2*complete_bound(pulses,F(1,10**8))+20*(reset+pointer)+150*F(2,10**11)
    assert base==F(lower['complete_sector_and_phase_halfdiamond_upper'])
    costs=[]
    for cost in r['additional_coupled_costs']:
        v=complete_bound(cost['pulses'])+cost['reset_bits']*reset+cost['pointer_reads']*pointer
        assert v==F(cost['complete_bound']);costs.append(v)
    reader=base+max(costs)+8*(reset+pointer)+clock_error
    assert reader==F(r['complete_coupled_sector_reader_error_upper'])
    assert clock_error==F(r['complete_time4_countclock_error_upper'])
    flow=read(ENERGY/'SHARP_FLOW_CERTIFICATE.json');oneflow=F(flow['complete_U4_halfdiamond_upper'])
    assert 512*oneflow==F(flow['512_complete_U4_halfdiamond_upper'])
    errors=[]
    for i,c in enumerate(cells):
        x=512*oneflow+c['app_bound'];gain=x/(1-2*x);e=c['bad_bound']+gain
        assert x==F(r['linear512_preparation_word_and_apparatus_bounds'][i])
        assert gain==F(r['gain_corrected512_preparation_word_and_apparatus_bounds'][i])
        assert e==F(r['complete_source_preparation_error_upper'][i]);errors.append(e)
    total=sum(errors)+reader;assert total==F(r['complete_joint_preparation_and_energy_read_error_upper'])and total<F(1,1000)
    rho=[propagated(cells[0],0),propagated(cells[1],2)]
    flag=[];failure=[]
    for name,density in zip(('first','second'),rho):
        flag.append(np.array([np.trace(e@density).real for e in sector[name+'_effects']]))
        failure.append(float(np.trace(sector[name+'_processor_failure']@density).real))
    pairflags=np.kron(flag[1],flag[0]);fpair=failure[0]+sum(flag[0])*failure[1]
    assert np.max(abs(flag[0]-r['source_flag_laws']['first']))<1e-12
    assert np.max(abs(flag[1]-r['source_flag_laws']['second']))<1e-12
    assert np.max(abs(pairflags-r['source_flag_laws']['pair']))<1e-12
    tablechecks=0
    for timestr,group in r['native_flag_controlled_words'].items():
        t=float(F(timestr))
        for name,meta in group.items():
            ispair=name=='pair';fl=pairflags if ispair else flag[0];up=fpair if ispair else failure[0]
            rec=r['complete_source_coupled_records'][timestr][name];verify_law(rec,16 if ispair else 4)
            nf=2 if ispair else 1;ref=2 if ispair else 1;cols=1<<nf;dim=1<<meta['bits'];energy=12 if ispair else 6
            for setting in ('X','Y'):
                raw=np.asarray(meta['complete_per_input_flag_output_centres'][setting]);f=np.asarray(meta['per_input_flag_processor_failure_centres'][setting])
                assert raw.shape==(cols,dim) and max(abs(raw.sum(axis=1)+f-1))<1e-8
                ideal=np.zeros_like(raw)
                for j in range(cols):
                    angle=-energy*t if j==cols-1 else 0
                    expectation=math.cos(angle)if setting=='X'else math.sin(angle)
                    ideal[j,j]=(1+expectation)/2;ideal[j,j+(1<<ref)]=(1-expectation)/2
                assert np.max(abs(raw-ideal))<1e-7
                assert np.max(abs(fl@raw-rec[setting]['diagnostic_probability_centres']))<1e-12
                assert abs(float(fl@f)-rec[setting]['downstream_processor_failure_centre'])<1e-12
                assert abs(up-rec[setting]['upstream_sector_processor_failure_centre'])<1e-12;tablechecks+=cols
    original=json.loads(gzip.decompress((LOADER/'RECORDS.json.gz').read_bytes()))
    source_trials={t['index']:t for t in original['trials']if t['law']=='native'}
    assert len(r['all12_actual_signed_source_handoffs'])==12
    meta=r['native_flag_controlled_words']['1/10']['pair']
    for handoff in r['all12_actual_signed_source_handoffs']:
        original=source_trials[handoff['trial_index']];p,n=divmod(original['pair_code'],256)
        assert handoff['pair_code']==original['pair_code']and handoff['positive_address']==p and handoff['negative_address']==n
        densities=[propagated(cells[0],p),propagated(cells[1],n-12)]
        fl=[];f=[]
        for name,density in zip(('first','second'),densities):
            fl.append(np.array([np.trace(e@density).real for e in sector[name+'_effects']]))
            f.append(float(np.trace(sector[name+'_processor_failure']@density).real))
        joint=np.kron(fl[1],fl[0]);assert np.max(abs(joint-handoff['pair_flag_centres']))<1e-12
        record=handoff['pair_energy_record'];verify_law(record,16)
        for setting in ('X','Y'):
            computed=joint@np.asarray(meta['complete_per_input_flag_output_centres'][setting])
            assert np.max(abs(computed-record[setting]['diagnostic_probability_centres']))<1e-12
        assert abs(handoff['pair_energy_centre']-12)<1e-7
    route=WORK/'2026-11-16_codex_signed_pair_routing';routes=json.loads(gzip.decompress((route/'RECORDS.json.gz').read_bytes()))
    wanted=[0]*240;wanted[0]=1;wanted[14]=-1
    routed=[x for x in routes['trials']if x['law']=='native'];assert len(routed)==12
    for x in routed:assert x['success']and x['observed_source']==wanted and x['positive_frame']==[0,1]
    route_error=F(read(route/'RECEIPT.json')['per_trial_timeout_upper'])
    # Route failures are terminal. Native instrument error below also covers
    # conditional numerical-centre normalization and raw-map gain.
    joined_reg=read(ENERGY/'PREREG_ROUTED_JOIN.json')
    for group,directory in [('files',ENERGY),('route_files',route)]:
        for name,h in joined_reg[group].items():assert sha(directory/name)==h;bindings[str((directory/name).relative_to(ROW))]=h
    joined_receipt=read(ENERGY/'routed_join/RECEIPT.json');joined=joined_receipt['scientific']
    assert joined_receipt['prereg_sha256']==sha(ENERGY/'PREREG_ROUTED_JOIN.json')
    assert joined_receipt['maps_sha256']==sha(ENERGY/'routed_join/maps.npz')
    assert joined['route_manifest_sha256']==sha(route/'MANIFEST.sha256')
    assert joined['coupled_receipt_sha256']==sha(ENERGY/'coupled_energy/RECEIPT.json')
    assert F(joined['route_timeout_upper'])==route_error and joined['literal_new_native_calls']==0
    linear=total+route_error;final_error=linear/(1-2*linear)
    assert F(joined['linear_physical_TP_composition_bound'])==linear
    assert F(joined['complete_join_error_upper'])==final_error
    source_byindex={x['index']:x for x in routed};arrays=np.load(ENERGY/'routed_join/maps.npz')
    assert len(joined['handoffs'])==12
    for row in joined['handoffs']:
        index=row['trial_index'];prior=source_byindex[index]
        assert row['input_pair_code']==prior['input_pair_code']and row['actual_routed_pair_code']==14
        assert row['actual_source']==wanted and row['native_frame']==[0,1]
        assert row['route_record_sha256']==hashlib.sha256(json.dumps(prior,sort_keys=True).encode()).hexdigest()
        for k,name in enumerate(('first','second')):
            assert np.linalg.norm(arrays['trial_%02d_%s'%(index,name)]-rho[k],2)<1e-12
            assert np.linalg.norm(cells[k]['P'][:,0 if k==0 else 2])==0
        for timestr,group in row['complete_energy_output_centres'].items():
            for name,record in group.items():
                verify_law(record,4 if name=='single'else 16)
                canonical=r['complete_source_coupled_records'][timestr][name]
                for setting in('X','Y'):
                    assert np.max(abs(np.asarray(record[setting]['diagnostic_probability_centres'])-canonical[setting]['diagnostic_probability_centres']))<1e-12
    interval={};longinterval={}
    single_linear=errors[0]+reader+route_error
    single_final=single_linear/(1-2*single_linear)
    for name,err in [('single',single_final),('pair',final_error)]:
        rec=r['complete_source_coupled_records']['1/10'][name]
        interval[name],radius=energy_interval(rec,err,.1)
        longrec=r['complete_source_coupled_records']['4'][name]
        raw=-math.atan2(*reversed(longrec['complex_accepted_fringe']))/4
        width=math.asin(radius)/4+1e-8
        candidates=[raw+k*math.pi/2 for k in range(32)if interval[name][0]-width<=raw+k*math.pi/2<=interval[name][1]+width]
        assert len(candidates)==1
        longinterval[name]=[candidates[0]-width,candidates[0]+width]
        energy=6 if name=='single'else 12;assert longinterval[name][0]<=energy<=longinterval[name][1]
    ratio=[longinterval['pair'][0]/longinterval['single'][1],longinterval['pair'][1]/longinterval['single'][0]]
    oldratio=math.log1p(-12/1260)/math.log1p(-6/1260)
    assert ratio[0]<2<ratio[1]<oldratio
    for name,rec in r['negative_source_records'].items():
        verify_law(rec,4);assert rec['complex_accepted_fringe']is None and rec['energy_success_probability']<1e-12
    assert len(r['vacuum_native_records'])==1260
    for a,b,c,j in r['vacuum_native_records']:assert j==0 and c in(0,1)
    assert r['vacuum_generator']=='0'and r['vacuum_transfer_eigenvalue']=='1'
    assert r['fixed_count_phase_nonadditivity']>.05
    for n in ('PREREG_COUPLED_ENERGY.json','routed_join/RECEIPT.json','routed_join/maps.npz','PREREG_ROUTED_JOIN.json','coupled_energy/RECEIPT.json','coupled_energy/maps.npz','poisson_energy/RECEIPT.json','poisson_energy/source_sector_maps.npz','SHARP_FLOW_CERTIFICATE.json','first_cell_certified/RECEIPT.json','first_cell_certified/maps.npz','second_cell/RECEIPT.json','second_cell/maps.npz'):
        bindings[str((ENERGY/n).relative_to(ROW))]=sha(ENERGY/n)
    result={'independent_assembly_checks_pass':True,'producer_modules_imported':False,'exact_signed_source_mode_columns':source_columns,
        'per_flag_native_pointer_rows_checked':tablechecks,'unrouted_actual_source_handoffs':12,'same_source_routed_handoffs':12,
        'complete_join_error_with_route_timeout':float(final_error),'independent_short_energy_intervals':interval,
        'long_intervals_unwrapped_only_by_short_records':longinterval,'pair_to_single_ratio_interval':ratio,
        'old_fixed_count_ratio':oldratio,'pair_keV_at_single511_anchor_interval':[511*x for x in ratio],
        'native_route_consumer_binding_in_final_producer':'Verified separately frozen routed_join; all12 source prefixes/densities and complete records agree.',
        'scientific_scope':'Complete energy-output arithmetic and exact source values; physical interpretation/other-group G6 not certified here.',
        'bindings':bindings,'seconds':time.monotonic()-started}
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    if args.record:
        with(HERE/'FINAL_RESULTS.json').open('x')as f:f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else:
        expected=read(HERE/'FINAL_RESULTS.json');expected.pop('seconds')
        actual=dict(result);actual.pop('seconds');assert actual==expected
    print(json.dumps({k:v for k,v in result.items()if k!='bindings'},indent=2))
if __name__=='__main__':main()
