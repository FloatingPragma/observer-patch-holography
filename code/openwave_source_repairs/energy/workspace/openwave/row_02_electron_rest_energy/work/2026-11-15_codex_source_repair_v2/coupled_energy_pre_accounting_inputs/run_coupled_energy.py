#!/usr/bin/env python3
"""Complete source-flag-controlled native energy records."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time,importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
sys.path.insert(0,str(OW/'shared/compact_logarithmic_flow'))
import native_words as nw
import reference,certify
from run_phase import controlled_phase

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
        'source_flag_columns':columns,'pulses':len(word),'word_sha256':nw.digest_word(word)}

def observe(flaglaw,instrument):
    records,failures,meta=instrument;out={};conditional=[];accepted=[]
    for setting in ('X','Y'):
        law=np.asarray(flaglaw)@records[setting]
        p0,p1=[float(law[i]) for i in meta['success_outputs']];success=p0+p1
        out[setting]={'complete_probabilities':law.tolist(),'energy_success_probability':success,
                      'processor_failure':float(np.asarray(flaglaw)@failures[setting])}
        conditional.append(None if success<1e-12 else (p0-p1)/success);accepted.append(success)
    z=None if any(x is None for x in conditional) else complex(*conditional)
    out['complex_accepted_fringe']=None if z is None else [z.real,z.imag]
    out['energy_success_probability']=min(accepted)
    return out,z

def pump(cell):
    path=HERE/('first_cell_certified' if cell==0 else 'second_cell')
    a=np.load(path/'maps.npz');Q=a['certified_Q_proposal'];D=a['dark_integer_basis'][cell:cell+12]
    P=D@np.linalg.inv(D.T@D)@D.T;rho=np.zeros((12,12),complex);rho[0 if cell==0 else 2,0 if cell==0 else 2]=1
    escape=np.eye(12)-Q.conj().T@Q
    for _ in range(512):rho=Q@rho@Q.conj().T+np.trace(escape@rho)*np.eye(12)/12
    assert abs(np.trace(rho)-1)<1e-7 and np.trace(P@rho).real>.99
    return rho,P

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
    assert flags0[1]>.99 and flags1[1]>.99
    pairflags=np.array([flags0[a]*flags1[b] for b in (0,1) for a in (0,1)])
    held=np.zeros((12,12));held[0,0]=1
    holdflags=[float(np.trace(e@held).real) for e in sector['first_effects']]
    witness_laws,witness,failed_witness=witness_sector_labels([0,240,255])
    records={};instruments={};energies={};additivity={}
    for t in (F(1,10),F(4)):
        one=controlled_records(E,t);pair=controlled_records(EP,t,True)
        a,za=observe(flags0,one);b,zb=observe(pairflags,pair)
        assert za is not None and zb is not None and abs(zb-za*za)<1e-7
        records[str(t)]={'single':a,'pair':b};instruments[str(t)]={'single':one[2],'pair':pair[2]};additivity[str(t)]=abs(zb-za*za)
        if t==F(1,10):
            er=-math.atan2(za.imag,za.real)/float(t);pr=-math.atan2(zb.imag,zb.real)/float(t)
            assert abs(er-float(E))<1e-6 and abs(pr-float(EP))<1e-6 and abs(pr/er-2)<1e-6
            energies={'single':er,'pair':pr,'ratio':pr/er,'decoder':'-atan2(Y,X)/t;0<=single<=12,0<=pair<=24,t=1/10'}
            negatives={}
            for name,flaglaw in [('HOLD_source_pump',holdflags),('off_sector',witness_laws[0]),('vacuum',witness_laws[1]),('invalid',witness_laws[2])]:
                rr,z=observe(flaglaw,one);assert rr['energy_success_probability']<1e-12 and z is None;negatives[name]=rr
    oldE=sum(certify.energy_interval(1260-int(E),1260))/2;oldP=sum(certify.energy_interval(1260-int(EP),1260))/2
    olda,za=observe(flags0,controlled_records(oldE,F(4)));oldb,zb=observe(pairflags,controlled_records(oldP,F(4),True))
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
    # Source-sector compilation includes both cell projectors. Additional
    # coupled flag logic/pointers use far fewer than2000 primitive pulses.
    phase_error=F(lower['complete_sector_and_phase_halfdiamond_upper'])+F(1,10**7)+clock_error
    flow=F(json.loads((HERE/'SHARP_FLOW_CERTIFICATE.json').read_text())['512_complete_U4_halfdiamond_upper'])
    apps=[F(json.loads((HERE/p/'RECEIPT.json').read_text())['scientific']['512epoch_apparatus_error_upper']) for p in ('first_cell_certified','second_cell')]
    prep_error=[F(1,65536)+flow+apps[0],F(1,4294967296)+flow+apps[1]]
    joint_complete=sum(prep_error)+phase_error
    assert joint_complete<F(1,1000)
    out=HERE/'coupled_energy';out.mkdir(exist_ok=False)
    np.savez_compressed(out/'maps.npz',source_pump_first=rho0,source_pump_second=rho1,
        source_flags_first=np.asarray(flags0),source_flags_second=np.asarray(flags1),
        source_flags_pair=pairflags,negative_source_sector_columns=witness,
        negative_source_processor_failure=failed_witness)
    result={'candidate':'A-native-poisson-compact-energy-v1','source_clock_manifest_sha256':sha(clock/'MANIFEST.sha256'),
        'source_flag_laws':{'first':flags0,'second':flags1,'pair':pairflags.tolist()},
        'complete_source_coupled_records':records,'native_flag_controlled_words':instruments,
        'decoded_short_time_energies':energies,'measured_additivity_residuals':additivity,
        'negative_source_records':negatives,'fixed_count_control_records':{'single':olda,'pair':oldb},
        'fixed_count_phase_nonadditivity':olderr,'vacuum_native_calls':len(vacuum_records),
        'vacuum_native_records':vacuum_records,'vacuum_transfer_eigenvalue':'1','vacuum_generator':'0',
        'vacuum_complete_native_address_law':vacuum_address,
        'complete_time4_countclock_error_upper':str(clock_error),
        'complete_coupled_sector_reader_error_upper':str(phase_error),
        'complete_source_preparation_error_upper':list(map(str,prep_error)),
        'complete_joint_preparation_and_energy_read_error_upper':str(joint_complete),
        'native_processor_counts':nw.COUNTS,'zero_energy_assigned_to_invalid':False,
        'raw_fullpump_execution_used':False,'source_input_representation':'intervalcertifiedfullsourceprogram withcompleteimplementationbounds',
        'full_row_solved':False,'other_group_review':False}
    write(out/'RECEIPT.json',{'scientific':result,'artifacts':{'maps.npz':sha(out/'maps.npz')},
        'prereg_sha256':sha(HERE/'PREREG_COUPLED_ENERGY.json'),'elapsed_seconds':time.monotonic()-started})
    print(json.dumps({'coupled_energy_pass':True,'decoded':energies,'old_clock_mismatch':olderr,
        'joint_complete_error':float(joint_complete),'seconds':time.monotonic()-started}))

if __name__=='__main__':main()
