#!/usr/bin/env python3
"""Common-source integration. No producer imports; complete instruments retained."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,gzip,hashlib,json,math,sys,time
from functools import lru_cache
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
I=HERE/'inputs'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def gain(x):return x/(1-2*x)
def extend(x,eta):return gain(x/(1+2*x)+eta)
def up(x,scale=10**15):return F(-((-x.numerator*scale)//x.denominator),scale)
def tv(x,y):return float(np.sum(abs(x-y))/2)

def census(c):
    vertices=list(range(12*c['cell'],12*c['cell']+12));delta={v:0 for v in vertices}
    equations=[];tape=c['tape'];cursor=0
    while cursor<len(tape):
        a,b,coin,j=tape[cursor];before=delta[a]-delta[b];j0=j;c0=coin
        delta[a]-=j;delta[b]+=j;cursor+=1
        while cursor<len(tape):
            x,y,coin,j=tape[cursor];assert(a,b)==(x,y)
            delta[a]-=j;delta[b]+=j;cursor+=1
            if coin!=c0:equations.append((a,b,2*j0+j-before));break
        else:raise AssertionError('Unresolved census')
    known={(0,14)[c['cell']]:0}
    while len(known)<12:
        n=len(known)
        for a,b,d in equations:
            if a in known:known[b]=known[a]-d
            if b in known:known[a]=known[b]+d
        assert len(known)>n
    initial=[known[v]for v in vertices]
    assert initial==c['decoded']['initial']
    assert [known[v]+delta[v]for v in vertices]==c['decoded']['terminal']
    return initial

def source_join():
    receipt=read(I/'extraction/RECEIPT.json');raw=(I/'extraction/RECORDS.json.gz').read_bytes()
    assert sha(I/'extraction/RECORDS.json.gz')==receipt['records_sha256']
    data=json.loads(gzip.decompress(raw));expected=[0]*240;expected[0]=1;expected[14]=-1
    supplied=read(I/'extraction/HANDOFFS.json');handoffs={x['harness_K_audit_only']:x for x in supplied['handoffs']}
    edges={frozenset(x)for x in read(I/'extraction/graph.json')['seams']}
    counts=Counter();prefixes=[];zero=0;control=0;reversals=0
    for t in data['trials']:
        initial=[0]*240
        for c in t['censuses']:initial[12*c['cell']:12*c['cell']+12]=census(c)
        census_tape=[r for c in t['censuses']for r in c['tape']]
        assert [e['record']for e in t['events'][:len(census_tape)]]==census_tape
        assert all(e['law']=='native'for e in t['events'][:len(census_tape)])
        # The reader reconstructs before comparing any separately labelled audit.
        assert initial==t['initial_audit_only'] and initial[0]==initial[14]==0
        k=sum(initial[:12]);assert -12<=k<=12 and sum(initial[12:24])==-k
        assert max(map(abs,initial))<=64 and k==t['harness_K_audit_only']
        q=initial[:]
        for event in t['events']:
            a,b,coin,j=event['record'];assert frozenset((a,b))in edges
            before=q[a];old=q[a]**2+q[b]**2
            if event['law']=='native':
                total=q[a]+q[b];q[a]=total//2+coin*(total%2);q[b]=total-q[a]
                assert q[a]**2+q[b]**2<=old
            else:
                assert event['law']=='swap_hold'
                if coin:q[a],q[b]=q[b],q[a]
                assert q[a]**2+q[b]**2==old
            assert before-q[a]==j;counts[event['law']]+=1
        assert q==t['observed_source']==t['source_audit_only']
        if t['law']=='native'and k:
            h=handoffs[k];assert q==expected==h['actual_source']==h['observed_source']
            assert h['pair_code']==14 and h['positive_frame']==[0,1]
            assert not h['numerical_total_supplied_to_observer']and not h['source_reset']
            digest=hashlib.sha256(json.dumps(t['events'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
            assert digest==h['complete_prefix_sha256']
            frame=[e['record'][:2]for e in t['events']if e['kind']=='frame_move'and e['record'][3]]
            assert frame==[[0,1],[1,0]]
            if k<0:assert any(e['kind']=='orientation_move'and e['record'][3]for e in t['events']);reversals+=1
            prefixes.append(digest)
        elif t['law']=='native':assert t['outcome']=='no_unit'and not any(q);zero+=1
        else:assert t['outcome']=='selection_rejected';control+=1
    assert len(prefixes)==24 and zero==1 and control==25 and reversals==12
    eta=F(22,2**63)+F(24*64**2//2*11+22+11*237+3*239+2,2**64)
    assert eta==F(receipt['preparation_inconclusive_upper'])==F(8501,2**58)
    return {'actual_handoffs':24,'physical_sign_reversals':reversals,'no_unit':zero,'changed_law_rejections':control,
            'reconstructed_main_events':dict(counts),'all24_actual_prefixes':prefixes,'inconclusive_upper':str(eta)},eta

def calibration_join():
    old=read(I/'calibration/ROUTE_RECORDS.json')
    raw=(I/'calibration/inputs/bridge_RECORDS.json.gz').read_bytes()
    bridge_receipt=read(I/'calibration/inputs/bridge_RECEIPT.json')
    assert hashlib.sha256(raw).hexdigest()==bridge_receipt['evidence_sha256']
    trials={t['index']:t for t in json.loads(gzip.decompress(raw))['trials']if t['selection_law']=='native'}
    edges={frozenset(e)for e in read(I/'calibration/inputs/graph.json')['seams']}
    reconstructed=0;route_events=0;prefixes=[]
    def event(q,row):
        a,b,c,j=row;assert frozenset((a,b))in edges
        before=q[a];total=q[a]+q[b];q[a]=(total+c)//2;q[b]=total-q[a]
        assert before-q[a]==j
    for route in old['trials']:
        source=trials[route['bridge_trial']];q=source['initial_audit_only'][:]
        assert sum(q)==1 and not any(q[12:])
        for selection in source['selection']:
            for row in selection['reference']+[selection['body']]:event(q,row);reconstructed+=1
        for row in source['census']:event(q,row);reconstructed+=1
        handoff=source['basis_handoff']
        assert q==handoff['observed_source']==handoff['actual_source_audit_only']==route['input']
        assert sum(abs(v)for v in q)==1 and q[handoff['occupied']]==1
        prefix={'selection':source['selection'],'census':source['census']}
        assert hashlib.sha256(canonical(prefix)).hexdigest()==handoff['prefix_sha256']
        for leg in route['legs']:
            for row in leg['records']:event(q,row);route_events+=1
        assert q==[1]+[0]*239==route['observed_final_source']==route['actual_final_source_audit_only']
        assert route['frame']==[0,1]and route['frame_native_event'][3]==1
        prefix={'bridge_prefix':handoff['prefix_sha256'],'route':route['legs']}
        assert hashlib.sha256(canonical(prefix)).hexdigest()==route['prefix_sha256']
        prefixes.append(route['prefix_sha256'])
    assert len(prefixes)==36
    q=[1]+[0]*239
    for leg in old['forced_timeout']['legs']:
        for row in leg['records']:event(q,row);route_events+=1
    assert q==old['forced_timeout']['final_source']and old['forced_timeout']['status']=='route_timeout'
    route_error=F(old['route_failure_upper']);assert route_error==1-(1-F(1,2**32))**5
    return {'actual_positive_only_calibrations':36,'source_prefix_events_reconstructed':reconstructed,
      'route_events_reconstructed':route_events,'actual_source_prefixes':prefixes,
      'route_failure_upper':str(route_error),'initial_boundary':'Heralded actual positive-unit checkpoint from declared total-one calibration-source ancestry; original selection/census records and every failed alternative retained.',
      'independent_unit_chart_not_signed_pair_spectator':True},route_error

def pump(q,seed,dephase=False):
    state=np.zeros((12,12),complex);state[seed,seed]=1;escape=np.eye(12)-q.conj().T@q
    for _ in range(512):
        state=q@state@q.conj().T+np.trace(escape@state)*np.eye(12)/12
        if dephase:state=np.diag(np.diag(state))
    return state

def projector(d):
    g=d.T@d;det=int(g[0,0]*g[1,1]-g[0,1]*g[1,0])
    inverse=[[F(int(g[1,1]),det),F(-int(g[0,1]),det)],[F(-int(g[1,0]),det),F(int(g[0,0]),det)]]
    return [[sum(F(int(d[a,j]))*inverse[j][k]*int(d[b,k])for j in range(2)for k in range(2))for b in range(len(d))]for a in range(len(d))]

def loaded(block,rho,effect):
    k=block[:,:,:12]
    return sum(effect[i,j]*k[j]@rho@k[i].conj().T for i in range(2)for j in range(2))
def address(rho,pointer):
    args=[rho.reshape((2,)*16),list(range(16))]
    for bit in range(8):args.extend([pointer,[16+bit,8+bit,bit]])
    return np.einsum(*args,list(range(16,24)),optimize='greedy').real.reshape(256)

def complete_joint(rhos,m,pointer):
    flags=[];positions=[];scratch=[];dirty=[]
    for k,rho in enumerate(rhos):
        states=[loaded(m['loader%d_clean_columns'%k],rho,e)for e in pointer]
        flags.append([float(np.trace(s).real)for s in states]);positions.append(address(states[1],pointer))
        scratch.append([float(np.trace(g[:12,:12]@rho).real)for g in m['loader%d_scratch_word_grams'%k]])
        dirty.append(float(np.trace(m['loader%d_dirty_gram'%k][:12,:12]@rho).real))
    law=np.outer(*positions);rejects=[flags[0][a]*flags[1][b]for a in range(2)for b in range(2)if(a,b)!=(1,1)]
    dirty1=sum(flags[0])*dirty[1];residual=float(1-law.sum()-sum(rejects)-dirty[0]-dirty1)
    return law,{'all_region_flag_centres':flags,'region_reject_centres':rejects,'scratch_words_comparison_centres':scratch,
       'dirty_first_centre':dirty[0],'dirty_second_after_first_clean_centre':dirty1,
       'unresolved_processor_plus_rounding_residual':residual,'complete_position_mass_centre':float(law.sum()),
       'collision_address_centre':float(np.trace(law[:240,:240])),
       'invalid_address_centre':float(law.sum()-law[:240,:240].sum())}

def rotate(points):
    points=[tuple(tuple(v)for v in p)for p in points];lookup={p:i for i,p in enumerate(points)}
    assert len(lookup)==240
    p=[lookup[(tuple(-v for v in x),y,tuple(-v for v in z))]for x,y,z in points]
    assert sorted(p)==list(range(240))and all(p[p[i]]==i for i in range(240))
    return np.array(p+list(range(240,256)))

def observe(z,norm,matrices,keys,effects):
    pulses={k:matrices[str(i)]for i,k in enumerate(keys)}
    prep=pulses['Y/'+(-math.pi/4).hex()];plus=prep[:2,0]
    prepfail=float(np.vdot(prep[2:,0],prep[2:,0]).real)
    rho=np.array([[abs(plus[0])**2,plus[0]*plus[1].conjugate()*z.conjugate()],
                  [plus[1]*plus[0].conjugate()*z,abs(plus[1])**2*norm]],complex)
    result={};fringe=[]
    for setting in ('X','Y'):
        state=rho.copy();failure=0.
        order=[('Z',math.pi/4),('Y',math.pi/4)]if setting=='Y'else[('Y',math.pi/4)]
        for kind,angle in order:
            matrix=pulses[kind+'/'+angle.hex()]
            failure+=float(np.trace(matrix[2:,:2]@state@matrix[2:,:2].conj().T).real)
            state=matrix[:2,:2]@state@matrix[:2,:2].conj().T
        law=[float(np.trace(e@state).real)for e in effects]
        result[setting]={'pointer_output_centres':law,'reference_preparation_failure_centre':prepfail,
            'analysis_processor_failure_centre':failure,'arithmetic_mass_residual':float(np.trace(rho).real-failure-sum(law))}
        fringe.append((law[0]-law[1])/sum(law))
    result['fringe_centre']=fringe
    return result

def atan_small(x):
    assert abs(x)<=F(1,2)
    terms=64;value=sum(((-1)**k*x**(2*k+1)/F(2*k+1)for k in range(terms)),F())
    remainder=abs(x)**(2*terms+1)/F(2*terms+1)
    return value-remainder,value+remainder
@lru_cache(None)
def pi_interval():
    a,b=atan_small(F(1,5));c,d=atan_small(F(1,239))
    return 16*a-4*d,16*b-4*c
def atan_interval(x):
    pl,ph=pi_interval()
    if x<0:
        lo,hi=atan_interval(-x);return -hi,-lo
    if x>1:
        lo,hi=atan_interval(1/x);return pl/2-hi,ph/2-lo
    if x>F(1,2):
        lo,hi=atan_small((x-1)/(x+1));return pl/4+lo,ph/4+hi
    return atan_small(x)
def phase_interval(record):
    x,y=map(F,record['fringe_centre']);pl,ph=pi_interval()
    assert x or y
    if not x:return (pl/2,ph/2)if y>0 else(-ph/2,-pl/2)
    lo,hi=atan_interval(y/x)
    if x<0:return(lo+pl,hi+ph)if y>=0 else(lo-ph,hi-pl)
    return lo,hi
def outward_interval(lo,hi):
    return [math.nextafter(float(lo),-math.inf),math.nextafter(float(hi),math.inf)]
def decode(short,long,epsilon):
    # No expected energy, winding, anchor or ratio is an argument to this decoder.
    sqrt2upper=F(14142135623730951,10**16);assert sqrt2upper**2>2
    radius=2*sqrt2upper*epsilon/(1-epsilon);assert radius<1
    angle=radius/(1-radius);pl,ph=pi_interval()
    sl,sh=phase_interval(short);ll,lh=phase_interval(long)
    coarse_low=(-sh-angle)*10;coarse_high=(-sl+angle)*10
    raw_low,raw_high=-lh,-ll
    # Bounded positive energy from the full source norm (<=24) implies
    # nonnegative winding. Determine it solely from the measured intervals.
    assert 4*coarse_low-raw_high-angle>=0
    low=math.ceil((4*coarse_low-raw_high-angle)/(2*ph))
    high=math.floor((4*coarse_high-raw_low+angle)/(2*pl))
    assert low==high and low>=0,'Short records do not resolve alias'
    interval=outward_interval((raw_low+2*pl*low-angle)/4,(raw_high+2*ph*low+angle)/4)
    return {'energy_interval':interval,'short_interval':outward_interval(coarse_low,coarse_high),'winding_from_short':low,
       'complete_TV_radius_upper':str(epsilon),'conditional_complex_radius_upper':math.nextafter(float(radius),math.inf),
       'angle_radius_upper':math.nextafter(float(angle),math.inf),'decoder_arithmetic':'Exact rational alternating arctangent enclosures and Machin pi; outward float endpoints only.'}

def run():
    origins=read(HERE/'BRIDGE_ORIGINS.json')
    for name,component in origins['components'].items():
        assert sha(I/name/'MANIFEST.sha256')==component['manifest_sha256']
    for name,entry in origins['inputs'].items():
        manifest={n:h for h,n in(line.split('  ',1)for line in(I/entry['component']/'MANIFEST.sha256').read_text().splitlines())}
        assert sha(HERE/name)==entry['sha256']
        if entry['original_relative']!='MANIFEST.sha256':assert manifest[entry['original_relative']]==entry['sha256']
    joined,eta=source_join();calibration,calibration_route=calibration_join()
    common=read(I/'common_v1/RECEIPT.json');m=np.load(I/'common_v1/INTERFACE_MAPS.npz')
    first=np.load(I/'common_v1/inputs/energy/first_cell_certified/maps.npz')
    second=np.load(I/'common_v1/inputs/energy/second_cell/maps.npz')
    er=read(I/'energy/attempts/v2/RECEIPT.json');observer=read(I/'energy/observer/RECEIPT.json')
    clock=read(I/'energy/clock/RECEIPT.json')
    for prerequisite in (er,observer,clock):
        assert prerequisite['all_gates_pass']and all(prerequisite['gates'].values())
    rate=80640*((10*max(v['total_pulses']for v in er['words'].values())+80639)//80640)
    assert rate==clock['clock_slots_per_time_unit']and rate%1260==rate%1280==0
    for schedule in clock['programmes']:
        duration=F(schedule['duration']);pulses=er['words'][schedule['object']]['total_pulses']
        assert schedule['total_clock_slots']==rate*duration<2**40
        assert schedule['inherited_count_ticks']==1260*duration
        assert schedule['represented_active_processor_pulses']==pulses
        assert schedule['represented_idle_slots']==rate*duration-pulses>=0
        assert schedule['complete_clock_record_law']['length']==schedule['total_clock_slots']
    assert clock['native_branch_records']==[{'coin':c,'input':[0,0],'output':[0,0],'record':[c,0]}for c in(0,1)]
    ec=read(I/'energy/inputs/certificate.json');em=np.load(I/'energy/inputs/common_INTERFACE_MAPS.npz')
    assert ec['common_manifest_sha256']==origins['components']['common_v1']['manifest_sha256']
    assert ec['common_maps_sha256']==sha(I/'common_v1/INTERFACE_MAPS.npz')
    assert sha(I/'energy/inputs/common_INTERFACE_MAPS.npz')==sha(I/'common_v1/INTERFACE_MAPS.npz')
    assert read(I/'energy/inputs/graph.json')['seams']==read(I/'extraction/graph.json')['seams']
    assert er['source_census']['copied_pinned_native_kernel_sha256']==sha(I/'energy/inputs/native_event.py')
    source_manifest={n:h for h,n in(line.split('  ',1)for line in(I/'extraction/MANIFEST.sha256').read_text().splitlines())}
    assert source_manifest['native_event.py']==sha(I/'energy/inputs/native_event.py')
    assert sha(I/'calibration/native_event.py')==sha(I/'energy/inputs/native_event.py')
    assert np.array_equal(np.load(I/'calibration/inputs/first_cell_maps.npz')['certified_Q_proposal'],first['certified_Q_proposal'])
    # The inherited two-body allowance includes at least this five-move route.
    assert calibration_route<1-(1-F(1,2**32))**8
    assert F(ec['common_full_preparation_and_apparatus_allowance'])==F(common['spatial_complete_error_upper'])
    calibration_transfer=read(I/'calibration/TRANSFER_RESULT.json')
    calibration_bound=gain(F(calibration_transfer['complete_native_error_upper'])+F(calibration_transfer['ideal_bad_weight_upper'])+F(common['complete_broad_loader_error_upper']))
    assert calibration_bound<F(ec['common_full_preparation_and_apparatus_allowance'])
    calibration['complete_preparation_and_load_upper']=str(up(calibration_bound))
    calibration['conservatively_covered_by_energy_inherited_allowance']=True
    basis=np.array(ec['integer_modes'],int)
    assert np.array_equal(basis[:2].T,first['dark_integer_basis'])
    assert np.array_equal(basis[2:].T,second['dark_integer_basis'])
    a=np.column_stack((basis[0]/2,(2*basis[1]-basis[0])/math.sqrt(12)))
    b=np.column_stack((basis[2]/2,(2*basis[3]-basis[2])/math.sqrt(12)))
    pair=np.stack([np.outer(a[:,i],b[:,j])for i in range(2)for j in range(2)],axis=-1)
    rho0=pump(first['certified_Q_proposal'],0);rho1=pump(second['certified_Q_proposal'],2)
    control=pump(first['certified_Q_proposal'],0,True)
    joint=np.kron(rho0,rho1);spatial=read(I/'spatial/RECEIPT.json');sm=np.load(I/'spatial/PAIR_MAPS.npz')
    density_error=max(float(np.max(abs(joint-sm['prepared_full144_density']))),float(np.max(abs(rho0-m['first_final_density']))),
      float(np.max(abs(rho1-m['second_final_density']))),float(np.max(abs(control-m['dephasing_final_density']))))
    assert density_error<1e-12
    p=rotate(spatial['coords_exact_a_plus_b_phi']);assert np.array_equal(p,sm['global_Ry_address_permutation'])
    r0=projector(basis[:2,:12].T);r1=projector(basis[2:,12:24].T)
    target=np.kron(np.array(r0,float)/2,np.array(r1,float)/2)
    exact_return=all(r0[int(p[i])][int(p[j])]==r0[i][j]for i in range(12)for j in range(12))and all(r1[int(p[i+12]-12)][int(p[j+12]-12)]==r1[i][j]for i in range(12)for j in range(12))
    energy_target=pair[:12,12:24,:].reshape(144,4);energy_target=energy_target@energy_target.conj().T/4
    target_error=float(np.max(abs(target-energy_target)))
    assert target_error<1e-14 and np.array_equal(target,sm['selected_full144_density'])
    pointer=np.load(I/'common_v1/inputs/flow/attempts/stage2_002/evidence/reference_preparation.npz')['pointer_effects_C']
    arrays={'same_prepared_full144_density':joint,'same_selected_full144_density':target,'global_Ry_address_permutation':p}
    laws={};instrument_error=0.
    for name,x,y in [('native',rho0,rho1),('dephasing',control,rho1)]:
        law,outcomes=complete_joint((x,y),m,pointer);arrays[name+'_all65536_position_centres']=law
        instrument_error=max(instrument_error,float(np.max(abs(law-sm[name+'_all65536_position_centres']))))
        outcomes['joint_halfturn_TV_centre']=tv(law,law[np.ix_(p,p)]);laws[name]=outcomes
    assert instrument_error<1e-12
    spatial_error=up(extend(F(spatial['full_joint_complete_error_upper']),eta))
    control_error=up(extend(F(spatial['complete_control_error_upper']),eta))
    control_lower=F(spatial['active_dephasing_full_joint_TV_lower'])-2*(control_error-F(spatial['complete_control_error_upper']))
    # Round a lower bound down for compact independent reviewer data.
    control_lower=-up(-control_lower)
    native=np.load(I/'energy/observer/NATIVE_TERMINAL_MATRICES.npz')
    effects=np.load(I/'energy/inputs/native_pointer.npz')['effects'];source=np.load(I/'energy/attempts/v2/SOURCE_OUTPUTS.npz')
    records={};record_error=0.;decoded={};errors={};budget_parts={}
    for t in ('0.1','4'):
        records[t]={}
        for name,initial in [('single',a),('pair',pair)]:
            final=source[name+'_'+t.replace('.','p')]
            z=np.vdot(initial,final)/initial.shape[-1];norm=float(np.vdot(final,final).real/initial.shape[-1])
            record=observe(z,norm,native,observer['matrix_keys'],effects);records[t][name]=record
            for setting in ('X','Y'):
                record_error=max(record_error,max(abs(x-y)for x,y in zip(record[setting]['pointer_output_centres'],observer['records'][t][name][setting]['pointer_output_centres'])))
    for name in ('single','pair'):
        old=er['complete_budgets'][name];parts={k:F(v)for k,v in old.items()if k!='complete'}
        assert gain(sum(parts.values(),F()))==F(old['complete'])
        if name=='pair':parts['unknown_total_native_preparation']=eta
        errors[name]=up(gain(sum(parts.values(),F())))
        budget_parts[name]={k:str(up(v))for k,v in parts.items()}
        decoded[name]=decode(records['0.1'][name],records['4'][name],errors[name])
    one=decoded['single']['energy_interval'];two=decoded['pair']['energy_interval']
    ratio=[math.nextafter(two[0]/one[1],-math.inf),math.nextafter(two[1]/one[0],math.inf)]
    kev=[math.nextafter(511*ratio[0],-math.inf),math.nextafter(511*ratio[1],math.inf)]
    records['changed_seam']=observer['records']['changed_seam']
    records['full_HOLD']=observe(1+0j,1.,native,observer['matrix_keys'],effects)
    records['source_vacuum']=observe(1+0j,1.,native,observer['matrix_keys'],effects)
    records['unprepared_coordinate']=observer['records']['unprepared_coordinate']
    native_fringe=complex(*records['4']['single']['fringe_centre'])
    control_distance=abs(native_fringe-complex(*records['changed_seam']['fringe_centre']))
    control_radius=4*math.sqrt(2)*float(errors['single'])/(1-float(errors['single']))
    gates={'all24_actual_unknown_total_sources_join':joined['actual_handoffs']==24,
      'separate_actual_one_unit_calibration_joins_same_positive_preparation':calibration['actual_positive_only_calibrations']==36,
      'zero_and_all_controls_retained':joined['no_unit']==1 and joined['changed_law_rejections']==25,
      'same_full_preparation_and_broad_loader':density_error<1e-12,
      'same_whole_selected_body_for_both_assays':target_error<1e-14,
      'single_global_physical_halfturn_and_current_frame':exact_return and p[0]==1 and p[1]==0,
      'complete_unfiltered_full_joint_reader':instrument_error<1e-12 and all(abs(x['unresolved_processor_plus_rounding_residual'])<1e-7 for x in laws.values()),
      'finite_joint_return_bound':2*spatial_error<F(11,10000),
      'same_chart_active_dephasing_separates':control_lower>F(3,100),
      'same_native_energy_pointer_records':record_error<1e-12,
      'same_clock_realized_for_both_source_objects':len(clock['programmes'])==4,
      'source_law_control_separates_after_complete_error':control_distance>control_radius+.01,
      'two_independent_record_decoders_resolve_windings':one[0]<6<one[1]and two[0]<12<two[1],
      'single_anchor_pair_prediction_excludes_prior_excess':ratio[0]<2<ratio[1]and kev[1]<1024.45085,
      'no_new_literal_native_execution_claimed':True}
    gates={k:bool(v)for k,v in gates.items()}
    result={'candidate':'A-old240-record-selected-local-generator-v2','component_manifests':{k:v['manifest_sha256']for k,v in origins['components'].items()},
      'source_join':joined,'single_calibration_join':calibration,'gates':gates,'all_gates_pass':all(gates.values()),
      'full_prepared_density_max_difference':density_error,'full_selected_energy_spatial_density_max_difference':target_error,
      'complete_joint_instrument_max_difference':instrument_error,'native_energy_pointer_max_difference':record_error,
      'spatial_complete_error_upper':str(spatial_error),'full_joint_halfturn_TV_upper':str(2*spatial_error),
      'active_dephasing_full_joint_TV_lower':str(control_lower),'complete_control_error_upper':str(control_error),
      'full_joint_spatial_records':laws,'energy_records':records,'energy_decoded':decoded,
      'energy_complete_errors_upper':{k:str(v)for k,v in errors.items()},'energy_additive_parts_outward_upper':budget_parts,
      'pair_to_single_interval':ratio,'pair_keV_at511_anchor':kev,
      'common_clock':{'slots_per_time_unit':rate,'slots_per_inherited_count_tick':rate//1260,
        'programmes':clock['programmes'],'source_effect_error':'0 exactly','clock_and_counter_supplied':True},
      'source_law_fringe_distance_centre':control_distance,'source_law_two_arm_error_radius_upper':control_radius,
      'unresolved_upstream_records':'All extraction inconclusive, pump escape, first processor failure, padding, scratch, region rejection, invalid, and terminal bit records remain typed; no conditional success is replaced by probability one.',
      'compression':'24 identical successful actual-source handoffs share one evaluated preparation/assay map; reconstructed histories are not fresh native calls. Native word outputs and complete arbitrary-reference bounds are inherited with exact seals.',
      'literal_new_native_calls':0,'new_fullflow_executions':0,'native_terminal_bit_compositions':16,
      'finite_exact_halfturn_claimed':False,'same_group_only':True,'other_group_G6':False}
    return result,arrays

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args();start=time.monotonic()
    for n,h in read(HERE/'PREREG.json')['files'].items():assert sha(HERE/n)==h,n
    result,arrays=run()
    if args.record:
        with(HERE/'RECEIPT.json').open('xb')as f:f.write(canonical(result))
        assert not(HERE/'COMMON_MAPS.npz').exists();np.savez_compressed(HERE/'COMMON_MAPS.npz',**arrays)
    else:
        assert read(HERE/'RECEIPT.json')==json.loads(canonical(result))
        old=np.load(HERE/'COMMON_MAPS.npz')
        for k,v in arrays.items():assert np.array_equal(old[k],v),k
    print(json.dumps({'all_gates_pass':result['all_gates_pass'],'gates':result['gates'],
      'joint_TV_upper':float(F(result['full_joint_halfturn_TV_upper'])),'dephasing_lower':float(F(result['active_dephasing_full_joint_TV_lower'])),
      'pair_keV':result['pair_keV_at511_anchor'],'seconds':time.monotonic()-start},indent=2))
    assert result['all_gates_pass']
if __name__=='__main__':main()
