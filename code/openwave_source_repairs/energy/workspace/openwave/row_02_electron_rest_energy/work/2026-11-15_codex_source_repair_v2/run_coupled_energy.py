#!/usr/bin/env python3
"""Complete source-flag-controlled native energy records."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time,importlib.util,gzip
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
sys.path.insert(0,str(OW/'shared/compact_logarithmic_flow'))
import native_words as nw
import reference,certify
from run_phase import controlled_phase
from model.native_interface import pointer

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
def complete_effect(outcome,bits):
    E=np.ones((1,1),complex);local=reference.pointer_effects('C')
    for bit in reversed(range(bits)):E=np.kron(E,local[(outcome>>bit)&1])
    return E

def controlled_records(energy,t,pair=False):
    bits=4 if pair else 2;ref=2 if pair else 1;flag=3 if pair else 0;nf=2 if pair else 1
    columns=1<<nf;state=np.zeros((1<<bits,columns),complex)
    for j in range(columns):state[j,j]=1
    words,raw,failed,metadata=reference.routines()
    word=[nw.c.pulse('Y',-math.pi/4,[ref])]
    if pair:word+=nw.c.remap(words['ccx'],[0,1,3])
    theta=-float(F(t)*F(energy))
    for _ in range(32):word +=[nw.c.pulse('I',theta/32,[ref]),nw.c.pulse('J',-theta/32,[flag,ref])]
    if pair:word+=nw.c.remap(words['ccx'],[0,1,3])
    state,bad=nw.apply_word(word,state,bits)
    records={};failures={}
    for name,analysis in [('X',[nw.c.pulse('Y',math.pi/4,[ref])]),
        ('Y',[nw.c.pulse('Z',math.pi/4,[ref]),nw.c.pulse('Y',math.pi/4,[ref])])]:
        terminal,failed_read=nw.apply_word(analysis,state,bits)
        probabilities=[]
        for output in range(1<<bits):
            E=complete_effect(output,bits)
            probabilities.append(np.einsum('ai,ab,bi->i',terminal.conj(),E,terminal,optimize=True).real)
        records[name]=np.asarray(probabilities).T
        failures[name]=np.real(np.diag(bad+failed_read))
        assert max(abs(records[name].sum(axis=1)+failures[name]-1))<1e-7
    accept0=(1<<nf)-1;accept1=accept0|(1<<ref)
    # Scratch is required zero; all other outcomes remain terminal records.
    return records,failures,{'bits':bits,'reference_bit':ref,'success_outputs':[accept0,accept1],
        'source_flag_columns':columns,'pulses':len(word),'longest_analysis_pulses':2,
        'native_pointer_reads_per_setting':bits,'reset_bits_per_setting':2*bits,
        'source_flag_encoding_pulses_upper':2*nf,'word_sha256':nw.digest_word(word),
        'complete_per_input_flag_output_centres':{k:v.tolist() for k,v in records.items()},
        'per_input_flag_processor_failure_centres':{k:v.tolist() for k,v in failures.items()}}

def observe(flaglaw,instrument,upstream_failure=0.):
    records,failures,meta=instrument;out={};conditional=[];accepted=[]
    for setting in ('X','Y'):
        law=np.asarray(flaglaw)@records[setting]
        p0,p1=[float(law[i]) for i in meta['success_outputs']];success=p0+p1
        downstream=float(np.asarray(flaglaw)@failures[setting])
        out[setting]={'diagnostic_probability_centres':law.tolist(),'energy_success_probability_centre':success,
                      'downstream_processor_failure_centre':downstream,
                      'upstream_sector_processor_failure_centre':float(upstream_failure),
                      'arithmetic_normalization_residual':float(1-upstream_failure-downstream-law.sum())}
        conditional.append(None if success<1e-12 else (p0-p1)/success);accepted.append(success)
    z=None if any(x is None for x in conditional) else complex(*conditional)
    out['complex_accepted_fringe']=None if z is None else [z.real,z.imag]
    out['energy_success_probability']=min(accepted)
    out['centre_semantics']='Unmodified floating/interval-comparison centres; physical law is positive and normalized within the complete certified output radius.'
    return out,z

def pump(cell,initial=None):
    path=HERE/('first_cell_certified' if cell==0 else 'second_cell')
    a=np.load(path/'maps.npz');Q=a['certified_Q_proposal'];D=a['dark_integer_basis'][cell:cell+12]
    P=D@np.linalg.inv(D.T@D)@D.T;rho=np.zeros((12,12),complex)
    position=(0 if cell==0 else 2) if initial is None else initial-cell
    rho[position,position]=1
    escape=np.eye(12)-Q.conj().T@Q
    for _ in range(512):rho=Q@rho@Q.conj().T+np.trace(escape@rho)*np.eye(12)/12
    assert abs(np.trace(rho)-1)<1e-7 and np.trace(P@rho).real>.99
    return rho,P

def unresolved_address_reader(witness):
    blocks=witness.reshape(2,256,-1);M=pointer.local_matrix();native_effects=pointer.native_effects('C')
    densities=[np.zeros((256,256),complex) for _ in range(blocks.shape[2])]
    for z in range(4):
        K=M[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]*math.sqrt(native_effects[0][z]/2)
        amp=np.einsum('ab,bsj->asj',K,blocks)
        for j in range(blocks.shape[2]):densities[j]+=amp[:,:,j].T@amp[:,:,j].conj()
    effects=np.asarray(reference.pointer_effects('C'));result=[]
    for density in densities:
        args=[density.reshape((2,)*16),list(range(16))]
        for bit in range(8):args.extend([effects,[16+bit,8+bit,bit]])
        args.append(list(range(16,24)))
        law=np.einsum(*args,optimize='greedy').real.reshape(256)
        result.append({'complete_flag0_address_centres':law.tolist(),
            'measured_class_centres':{'charged_unresolved':float(sum(law[:240])),
                'physical_vacuum':float(law[240]),'invalid_padding':float(sum(law[241:]))}})
    assert result[0]['measured_class_centres']['charged_unresolved']>1-1e-7
    assert result[1]['measured_class_centres']['physical_vacuum']>1-1e-7
    assert result[2]['measured_class_centres']['invalid_padding']>1-1e-7
    return result

def energy_intervals(records,epsilon,pfloor):
    radius=2*math.sqrt(2)*epsilon/(pfloor-epsilon)
    assert 0<radius<1/10
    phase_bound=radius/(1-radius)
    short=complex(*records['1/10']['complex_accepted_fringe']);long=complex(*records['4']['complex_accepted_fringe'])
    es=-math.atan2(short.imag,short.real)*10
    short_interval=[es-10*phase_bound,es+10*phase_bound]
    base=-math.atan2(long.imag,long.real)/4
    candidates=[(k,base+k*math.pi/2) for k in range(32)
                if short_interval[0]-phase_bound/4<=base+k*math.pi/2<=short_interval[1]+phase_bound/4]
    assert len(candidates)==1
    k,center=candidates[0]
    return {'complete_error_used':epsilon,'acceptance_floor_used':pfloor,
        'complex_conditional_fringe_radius':radius,'phase_error_upper':phase_bound,
        'short_energy_interval':short_interval,'long_winding_from_short_interval':k,
        'long_energy_interval':[center-phase_bound/4,center+phase_bound/4]}

def witness_sector_labels(labels):
    source=np.zeros((512,len(labels)),complex)
    for j,p in enumerate(labels):source[p,j]=1
    source,failed=nw.apply_word([nw.c.pulse('Y',-math.pi/4,[8])],source,9)
    for places,values in [([4,5,6,7],[1,-1,-1,1]),([4,5,6,7,8,11],[-1,-1,-1,-1,2,2])]:
        integer=[0]*256
        for p,v in zip(places,values):integer[p]=v
        prep,err,worst=nw.preparation(integer)
        for word in [nw.c.inverse(prep),controlled_phase(math.pi),prep]:
            source,bad=nw.apply_word(word,source,9);failed+=bad
    source,bad=nw.apply_word([nw.c.pulse('Y',math.pi/4,[8])],source,9);failed+=bad
    blocks=np.stack((source[:256],source[256:]));ef=reference.pointer_effects('C')
    laws=np.asarray([np.real(np.diag(np.einsum('iac,ij,jad->cd',blocks.conj(),E,blocks,optimize=True))) for E in ef]).T
    assert max(laws[:,1])<1e-12
    return laws,source,failed

def main():
    started=time.monotonic();reg=json.loads((HERE/'PREREG_COUPLED_ENERGY.json').read_text())
    for name,digest in reg['files'].items():assert sha(HERE/name)==digest
    clock=HERE.parent/'2026-11-15_codex_poisson_source_clock'
    for name,digest in reg['clock_files'].items():assert sha(clock/name)==digest
    loading=HERE.parent/'2026-11-15_codex_signed_pair_loading'
    for name,digest in reg['loading_files'].items():assert sha(loading/name)==digest
    counted=json.loads((clock/'RECEIPT.json').read_text());E=F(counted['native_single_generator_values'][0]);EP=F(counted['native_pair_generator_values'][0])
    spec=importlib.util.spec_from_file_location('clock_native',clock/'native_event.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=np.load(OW/'shared/compact_logarithmic_flow/source_columns.npz');vacuum_records=[]
    for a,b in source['seams']:
        for coin in (0,1):
            q,r=module.native((0,)*240,(int(a),int(b)),coin)
            assert q==(0,)*240 and r==(coin,0);vacuum_records.append([int(a),int(b),coin,0])
    # Complete sector flags, computed from the actual native source maps.
    sector=np.load(HERE/'poisson_energy/source_sector_maps.npz')
    rho0,P0=pump(0);rho1,P1=pump(12)
    flags0=[float(np.trace(e@rho0).real) for e in sector['first_effects']]
    flags1=[float(np.trace(e@rho1).real) for e in sector['second_effects']]
    f0=float(np.trace(sector['first_processor_failure']@rho0).real)
    f1=float(np.trace(sector['second_processor_failure']@rho1).real)
    fpair=f0+sum(flags0)*f1
    assert flags0[1]>.99 and flags1[1]>.99
    pairflags=np.array([flags0[a]*flags1[b] for b in (0,1) for a in (0,1)])
    held=np.zeros((12,12));held[0,0]=1
    holdflags=[float(np.trace(e@held).real) for e in sector['first_effects']]
    witness_laws,witness,failed_witness=witness_sector_labels([0,240,255])
    address_records=unresolved_address_reader(witness)
    records={};instruments={};runtime_instruments={};energies={};additivity={}
    for t in (F(1,10),F(4)):
        one=controlled_records(E,t);pair=controlled_records(EP,t,True)
        a,za=observe(flags0,one,f0);b,zb=observe(pairflags,pair,fpair)
        assert za is not None and zb is not None and abs(zb-za*za)<1e-7
        records[str(t)]={'single':a,'pair':b};instruments[str(t)]={'single':one[2],'pair':pair[2]};additivity[str(t)]=abs(zb-za*za)
        runtime_instruments[str(t)]=(one,pair)
        if t==F(1,10):
            er=-math.atan2(za.imag,za.real)/float(t);pr=-math.atan2(zb.imag,zb.real)/float(t)
            assert abs(er-float(E))<1e-6 and abs(pr-float(EP))<1e-6 and abs(pr/er-2)<1e-6
            energies={'single':er,'pair':pr,'ratio':pr/er,'decoder':'-atan2(Y,X)/t;0<=single<=12,0<=pair<=24,t=1/10'}
            negatives={}
            for j,(name,flaglaw) in enumerate([('HOLD_source_pump',holdflags),('off_sector',witness_laws[0]),('vacuum',witness_laws[1]),('invalid',witness_laws[2])]):
                upstream=float(np.trace(sector['first_processor_failure']@held).real) if j==0 else float(failed_witness[j-1,j-1].real)
                rr,z=observe(flaglaw,one,upstream);assert rr['energy_success_probability']<1e-12 and z is None;negatives[name]=rr
    oldE=sum(certify.energy_interval(1260-int(E),1260))/2;oldP=sum(certify.energy_interval(1260-int(EP),1260))/2
    old_one=controlled_records(oldE,F(4));old_pair=controlled_records(oldP,F(4),True)
    olda,za=observe(flags0,old_one,f0);oldb,zb=observe(pairflags,old_pair,fpair)
    olderr=abs(zb-za*za);assert olderr>.05
    # Physical zero-load reference has zero generator. Its accepted-energy
    # reader says unresolved; a separate native address read identifies240.
    pointer=reference.pointer_effects('C');vacuum_address=[]
    for label in range(256):
        p=1.
        for bit in range(8):p*=float(pointer[(label>>bit)&1][(240>>bit)&1,(240>>bit)&1].real)
        vacuum_address.append(p)
    assert vacuum_address[240]>1-1e-7
    clock_error=5040*F(counted['per128bit_clock_sampler_failure_upper'])+161280*F(counted['per_source_address_failure_upper'])
    lower=json.loads((HERE/'poisson_energy/RECEIPT.json').read_text())['scientific']
    reset=F(35,2147483648);pointer_error=F(200000000001,20000000000000000000000)
    costs=[]
    for group in list(instruments.values())+[{'single':old_one[2],'pair':old_pair[2]}]:
        for name,m in group.items():
            pulse_total=m['pulses']+m['longest_analysis_pulses']+m['source_flag_encoding_pulses_upper']
            pulse_bound=nw.complete_bound(pulse_total)['half_diamond']
            reset_bound=m['reset_bits_per_setting']*reset
            read_bound=m['native_pointer_reads_per_setting']*pointer_error
            costs.append({'source':name,'pulses':pulse_total,'pulse_bound':str(pulse_bound),
                'reset_bits':m['reset_bits_per_setting'],'reset_bound':str(reset_bound),
                'pointer_reads':m['native_pointer_reads_per_setting'],'pointer_bound':str(read_bound),
                'complete_bound':str(pulse_bound+reset_bound+read_bound)})
    additional=max(F(c['complete_bound']) for c in costs)
    address_bound=8*(reset+pointer_error)
    phase_error=F(lower['complete_sector_and_phase_halfdiamond_upper'])+additional+address_bound+clock_error
    flow=F(json.loads((HERE/'SHARP_FLOW_CERTIFICATE.json').read_text())['512_complete_U4_halfdiamond_upper'])
    apps=[F(json.loads((HERE/p/'RECEIPT.json').read_text())['scientific']['512epoch_apparatus_error_upper']) for p in ('first_cell_certified','second_cell')]
    raw_prep=[flow+a for a in apps]
    gains=[x/(1-2*x) for x in raw_prep]
    prep_error=[F(1,65536)+gains[0],F(1,4294967296)+gains[1]]
    joint_complete=sum(prep_error)+phase_error
    assert joint_complete<F(1,1000)
    coarse={}
    for name,err in [('single',prep_error[0]+phase_error),('pair',joint_complete)]:
        rr={t:records[t][name] for t in records}
        coarse[name]=energy_intervals(rr,float(err),1-float(err))
    single_interval=coarse['single']['long_energy_interval'];pair_interval=coarse['pair']['long_energy_interval']
    ratio_interval=[pair_interval[0]/single_interval[1],pair_interval[1]/single_interval[0]]
    mass_interval=[511*x for x in ratio_interval]
    loaded=json.loads(gzip.decompress((loading/'RECORDS.json.gz').read_bytes()))
    handoffs=[]
    for trial in loaded['trials']:
        if trial['law']!='native':continue
        code=trial['pair_code'];p,n=divmod(code,256)
        assert trial['observed_source'][p]==1 and trial['observed_source'][n]==-1
        a,pa=pump(0,p);b,pb=pump(12,n)
        af=[float(np.trace(e@a).real) for e in sector['first_effects']]
        bf=[float(np.trace(e@b).real) for e in sector['second_effects']]
        joint=np.asarray([af[i]*bf[j] for j in (0,1) for i in (0,1)])
        fa=float(np.trace(sector['first_processor_failure']@a).real)
        fb=float(np.trace(sector['second_processor_failure']@b).real)
        result,z=observe(joint,runtime_instruments['1/10'][1],fa+sum(af)*fb)
        assert result['energy_success_probability']>.99 and z is not None
        handoffs.append({'trial_index':trial['index'],'pair_code':code,'positive_address':p,'negative_address':n,
            'first_flag_centres':af,'second_flag_centres':bf,'pair_flag_centres':joint.tolist(),
            'pair_energy_record':result,'pair_energy_centre':-10*math.atan2(z.imag,z.real)})
    assert len(handoffs)==12
    out=HERE/'coupled_energy';out.mkdir(exist_ok=False)
    np.savez_compressed(out/'maps.npz',source_pump_first=rho0,source_pump_second=rho1,
        source_flags_first=np.asarray(flags0),source_flags_second=np.asarray(flags1),
        source_flags_pair=pairflags,negative_source_sector_columns=witness,
        negative_source_processor_failure=failed_witness)
    result={'candidate':'A-native-poisson-compact-energy-v1','source_clock_manifest_sha256':sha(clock/'MANIFEST.sha256'),
        'source_flag_laws':{'first':flags0,'second':flags1,'pair':pairflags.tolist()},
        'complete_source_coupled_records':records,'native_flag_controlled_words':instruments,
        'decoded_short_time_energies':energies,'measured_additivity_residuals':additivity,
        'certified_energy_intervals':coarse,'certified_pair_to_single_ratio_interval':ratio_interval,
        'certified_pair_keV_at_single511_anchor_interval':mass_interval,
        'upstream_sector_processor_failure_centres':{'first':f0,'second':f1,'pair_first_failure':fpair},
        'all12_actual_signed_source_handoffs':handoffs,
        'signed_source_loader_manifest_sha256':sha(loading/'MANIFEST.sha256'),
        'actual_handoff_error_scope':'Universal compact-weight bound on classical energy outputs; no full-state PD/2 assertion for arbitrary inputaddresses.',
        'negative_source_records':negatives,'fixed_count_control_records':{'single':olda,'pair':oldb},
        'negative_unresolved_native_address_records':address_records,
        'fixed_count_phase_nonadditivity':olderr,'vacuum_native_calls':len(vacuum_records),
        'vacuum_native_records':vacuum_records,'vacuum_transfer_eigenvalue':'1','vacuum_generator':'0',
        'vacuum_complete_native_address_law':vacuum_address,
        'complete_time4_countclock_error_upper':str(clock_error),
        'complete_coupled_sector_reader_error_upper':str(phase_error),
        'additional_coupled_costs':costs,'unresolved_address_read_error_upper':str(address_bound),
        'linear512_preparation_word_and_apparatus_bounds':list(map(str,raw_prep)),
        'gain_corrected512_preparation_word_and_apparatus_bounds':list(map(str,gains)),
        'gain_formula':'x/(1-2x), x=512*(singleflow+singleapparatus); covers stored nonexactCP gain aswellasexactTP target.',
        'complete_source_preparation_error_upper':list(map(str,prep_error)),
        'complete_joint_preparation_and_energy_read_error_upper':str(joint_complete),
        'native_processor_counts':nw.COUNTS,'zero_energy_assigned_to_invalid':False,
        'raw_fullpump_execution_used':False,'source_input_representation':'intervalcertifiedfullsourceprogram withcompleteimplementationbounds',
        'probability_centre_semantics':'Tiny signed numerical centres retained verbatim; exact native complete instrument is positive/normalized and within the stated totalvariation radius. No clipping or renormalization.',
        'retained_upstream_failure_symbols':['sector_processor_complement','source_reset_timeout','source_pointer_implementation_unresolved'],
        'full_row_solved':False,'other_group_review':False}
    write(out/'RECEIPT.json',{'scientific':result,'artifacts':{'maps.npz':sha(out/'maps.npz')},
        'prereg_sha256':sha(HERE/'PREREG_COUPLED_ENERGY.json'),'elapsed_seconds':time.monotonic()-started})
    print(json.dumps({'coupled_energy_pass':True,'decoded':energies,'old_clock_mismatch':olderr,
        'joint_complete_error':float(joint_complete),'seconds':time.monotonic()-started}))

if __name__=='__main__':main()
