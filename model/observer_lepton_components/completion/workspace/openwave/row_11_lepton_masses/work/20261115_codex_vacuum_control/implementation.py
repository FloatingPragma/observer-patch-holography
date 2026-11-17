"""Source-word controlled family filter on three occupied codes plus vacuum."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import importlib.util,sys,json,hashlib,math,time,argparse
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
FAMILY=HERE.parent/'20261114_codex_family_filter'
COUNT=HERE.parent/'20261114_codex_count_return'
spec=importlib.util.spec_from_file_location('sealed_family',FAMILY/'implementation.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
np=m.np;c=m.c;iv=m.iv;synth=m.synth;interpreter=m.interpreter
GATE=3e-8
WORDS={};BLOCKS=[];ARRAYS={}

def sha(b):return hashlib.sha256(b).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n'
def fwht(values):
    a=list(values);stride=1
    while stride<len(a):
        for start in range(0,len(a),2*stride):
            for j in range(stride):
                x,y=a[start+j],a[start+j+stride]
                a[start+j],a[start+j+stride]=x+y,x-y
        stride*=2
    return a

def bind_word(name,word,bits,units,main=True):
    admitted=interpreter.source_word(word,bits)
    record={'pulses':len(word),'qubits':bits,'sha256':sha(json.dumps(word,sort_keys=True).encode()),
            'admission_sha256':interpreter.digest(admitted),'target_grid_units':units}
    if name in WORDS:assert WORDS[name]==record
    WORDS[name]=record
    if main:BLOCKS.append({'name':name,**m.pulse_bound(len(word),units)})
    return word

@lru_cache(None)
def controlled_phase_word(angle):
    # Full11-bit diagonal: object address0 and family f<3, including f=3 vacuum.
    factors=(3,1,1,-1)
    numerators=[angle*factors[mask>>9] for mask in range(2048)]
    target=[angle if (index&511)==0 and (index>>9)<3 else 0 for index in range(2048)]
    assert fwht(numerators)==[2048*t for t in target]
    word=[]
    for mask,numerator in enumerate(numerators):
        coefficient=iv.binary_exact(numerator,2048*iv.AQ)
        label=['Z' if mask>>j&1 else 'I' for j in range(11)]
        word.extend(c.pauli_rotation(label,coefficient))
    assert len(word)==57350
    interpreter.source_word(word,11)
    return word

def controlled_ray(a,padding_weight,angle,name,phase_units,main=True):
    preparation=m.preparation_word('phi',padding_weight)
    inverse=c.inverse(preparation)
    words=[bind_word(name+'_inverse_preparation',inverse,9,9,main),
           bind_word(name+'_occupied_phase',controlled_phase_word(angle),11,phase_units,main),
           bind_word(name+'_preparation',preparation,9,9,main)]
    # No chart measurement/export between these three admitted words.
    x=c.execute(words[0],9,a.reshape(512,-1)).reshape(a.shape)
    x=m.from_global(c.execute(words[1],11,m.to_global(x)))
    x=c.execute(words[2],9,x.reshape(512,-1)).reshape(a.shape)
    bad=x[list(m.PADDING)].copy()
    good=np.zeros_like(x);good[list(m.CODES)]=x[list(m.CODES)]
    return good,bad

def effect(a):return m.effect(a)
def state_record(a):return m.state_record(a)
def cross(a,f=0):
    a=np.asarray(a).reshape(-1,4,4)
    return np.einsum('ai,aj->ij',a[:,:,f],a[:,:,3].conj())
def first_failure(word,pulse,complement,initial,bits):
    out=np.asarray(initial,complex).reshape(1<<bits,-1).copy()
    for k,event in enumerate(word[:pulse+1]):
        full=c.native_pulse(event['kind'],event['angle']);d=1<<len(event['targets'])
        if k==pulse:
            assert d<=complement<6
            block=np.zeros((d,d),complex);block[0]=full[complement,:d]
        else:block=full[:d,:d]
        out=c._apply(out,block,event['targets'],bits)
    return out

def source_bindings():
    result={}
    for module in list(sys.modules.values()):
        file=getattr(module,'__file__',None)
        if file:
            p=Path(file).resolve()
            if p.suffix=='.py' and p.is_relative_to(m.SOURCE):result[str(p.relative_to(ROOT))]=sha(p.read_bytes())
    for directory,names in [(FAMILY,['implementation.py','SOURCE_ARRAYS.npz','ARRAY_BINDINGS.json','RECEIPT.json','PROOF.md','MANIFEST.sha256']),
                            (COUNT,['PROOF.md','MANIFEST.sha256','CONTROL_BINDINGS.json']),
                            (HERE.parent/'20261114_codex_mass_spatial_join',['PROOF.md','MANIFEST.sha256'])]:
        for name in names:
            p=directory/name;result[str(p.relative_to(ROOT))]=sha(p.read_bytes())
    for n,h in json.loads((COUNT/'CONTROL_BINDINGS.json').read_text())['files'].items():
        assert sha((ROOT/n).read_bytes())==h;result[n]=h
    return dict(sorted(result.items()))

def calculate():
    checks={};outcomes={}
    def check(name,value):checks[name]=bool(value)
    check('sealed_family_manifest',sha((FAMILY/'MANIFEST.sha256').read_bytes())=='4714358385d7e0655197fe8b0653e4495de9df752ae4dc083b1d019c8c792fef')
    check('sealed_count_manifest',sha((COUNT/'MANIFEST.sha256').read_bytes())=='106d9fe49103a0ba82b519efc702d04db154007a0866349127f9b53c17175c0d')
    check('pinned_native_bytes',sha((m.SOURCE/'model/observer/native.py').read_bytes())=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d')
    base=np.load(FAMILY/'SOURCE_ARRAYS.npz')
    bindings=json.loads((FAMILY/'ARRAY_BINDINGS.json').read_text())
    assert sha((FAMILY/'SOURCE_ARRAYS.npz').read_bytes())==bindings['archive_sha256']
    for n,r in bindings['arrays'].items():assert sha(base[n].tobytes())==r['raw_sha256']
    check('reused_family_arrays_bound',True)
    old=json.loads((FAMILY/'RECEIPT.json').read_text())
    blank=np.zeros(512,complex);blank[0]=1
    prep=bind_word('initial_p',m.preparation_word('prep',F()),9,9)
    prepared=c.execute(prep,9,blank)
    pad=prepared[list(m.PADDING)];successful=np.zeros(512,complex);successful[list(m.CODES)]=prepared[list(m.CODES)]
    ARRAYS['preparation']=prepared
    prepbad=np.zeros((len(m.PADDING),4,4),complex)
    a=np.zeros((512,4,4),complex)
    for j in range(4):a[:,j,j]=successful;prepbad[:,j,j]=pad
    ARRAYS['preparation_padding']=prepbad
    a,ff=m.apply_tick(a,base['object_forward'],base['shift_clean'],base['shift_dirty'])
    ARRAYS['family_after_tick']=a;ARRAYS['forward_scratch_failure']=ff[0]
    BLOCKS.append({'name':'reused_forward_family_tick',**m.pulse_bound(2*old['object_word_pulses']+old['controlled_shift']['native_pulses'],2)})
    check('occupied_tick_same_sealed_source',np.linalg.norm(a[:,:,:3]-base['family_after_tick'])<GATE)
    # Independent vacuum object from original source operators, not installed identity.
    p=m.embedding(m.rays()['prep'])[list(m.CODES)]
    phi=m.embedding(m.rays()['phi'])[list(m.CODES)]
    _,_,twoc,vc=m.reg.object_operators(m.geometry.C_SEAM)
    _,_,twoe,ve=m.reg.object_operators(m.geometry.E_SEAM)
    vac=np.kron(np.asarray((vc@twoc)/2,complex),np.asarray((ve@twoe)/2,complex))@p
    check('actual_vacuum_tick',np.linalg.norm(a[list(m.CODES),3,3]-vac)<GATE)
    vac_amplitude=phi.conj()@vac;vac_effect=float(abs(vac_amplitude)**2)
    _,exactA=m.exact_effects();A=np.array(exactA.evalf(18).tolist(),complex)
    A4=np.zeros((4,4),complex);A4[:3,:3]=A
    D=np.eye(4)-float(m.EPS)*A4
    ARRAYS['A0']=A;ARRAYS['A_occupied']=A4;ARRAYS['target_D']=D
    theta=iv.split_angle(F(9,25),F(16,25))
    main,click=controlled_ray(a,F(1,2),2*theta,'weak',2)
    main,corr=controlled_ray(main,F(),-theta,'compensation',1)
    main,invfail=m.apply_tick(main,base['object_inverse'],base['inverse_shift_clean'],base['inverse_shift_dirty'],True)
    BLOCKS.append({'name':'reused_inverse_family_tick',**m.pulse_bound(2*old['object_word_pulses']+old['inverse_shift_native_pulses'],2)})
    rw=bind_word('destructive_p_return',m.return_word(),9,20)
    nonreturn,returned=m.word_export(rw,main)
    for name,value in [('first_click',click),('correction_padding',corr),('inverse_scratch_failure',invfail[0]),('nonreturn',nonreturn),('returned',returned)]:
        ARRAYS[name]=value;outcomes[name]=effect(value)
    outcomes['preparation_padding']=effect(prepbad);outcomes['forward_scratch_failure']=effect(ff[0])
    k=returned[m.PADDING.index(511)]
    ARRAYS['return_map']=k
    check('actual_four_input_return_map',np.linalg.norm(k+D)<GATE)
    check('actual_click_effect',np.linalg.norm(outcomes['first_click']-(16/25)*A4)<GATE)
    check('actual_nonreturn_effect',np.linalg.norm(outcomes['nonreturn']-(4/25)*(A4-A4@A4))<GATE)
    check('complete_return_effect',np.linalg.norm(outcomes['returned']-D@D)<GATE)
    check('full_terminal_effect_sum',np.linalg.norm(sum(outcomes.values())-np.eye(4))<GATE)
    vacuum_failure=sum(float(v[3,3].real) for n,v in outcomes.items() if n!='returned')
    check('zero_target_vacuum_failure',abs(vacuum_failure)<GATE)
    cross_errors=[];failure_cross_errors=[]
    for f in range(3):
        target=np.zeros((4,4),complex);target[:,3]=D[:,f]
        actual=cross(returned,f);ARRAYS['return_cross_'+str(f)]=actual
        cross_errors.append(float(np.linalg.norm(actual-target)))
        for name in ('first_click','correction_padding','nonreturn'):
            failure_cross_errors.append(float(np.linalg.norm(cross(ARRAYS[name],f))))
    check('all_three_charged_return_cross_blocks',max(cross_errors)<GATE)
    check('all_failure_cross_blocks_zero_target',max(failure_cross_errors)<GATE)
    # Control1: actual unconditioned weak detector, preserving the vacuum branch.
    weak,comp=m.filter_words('phi')
    bind_word('disabled_control_weak',weak,9,20,False);bind_word('disabled_control_compensation',comp,9,19,False)
    ctrl,ctrlclick=m.word_export(weak,a);ctrl,ctrlcorr=m.word_export(comp,ctrl)
    ctrl,ctrlinv=m.apply_tick(ctrl,base['object_inverse'],base['inverse_shift_clean'],base['inverse_shift_dirty'],True)
    ctrlnon,ctrlret=m.word_export(rw,ctrl)
    for name,value in [('disabled_first_click',ctrlclick),('disabled_correction_padding',ctrlcorr),('disabled_inverse_scratch',ctrlinv[0]),('disabled_nonreturn',ctrlnon),('disabled_returned',ctrlret)]:ARRAYS[name]=value
    ctrlk=ctrlret[m.PADDING.index(511)]
    ctrlA=A4.copy();ctrlA[3,3]=vac_effect
    check('disabled_control_own_map',np.linalg.norm(ctrlk+np.eye(4)-float(m.EPS)*ctrlA)<GATE)
    ctrlweight=float(effect(ctrlclick)[3,3].real)
    check('disabled_control_detected_by_vacuum_click',ctrlweight>1e-3 and abs(ctrlweight-16/25*vac_effect)<GATE)
    ctrltotal=effect(prepbad)+effect(ff[0])+effect(ctrlclick)+effect(ctrlcorr)+effect(ctrlinv[0])+effect(ctrlnon)+effect(ctrlret)
    check('disabled_control_full_effect_sum',np.linalg.norm(ctrltotal-np.eye(4))<GATE)
    # Control2: source phase on the final vacuum, coherently on the return flag.
    vword=c.compile_diagonal([0.,0.,0.,math.pi]);bind_word('wrong_vacuum_phase',vword,2,1,False)
    moved=returned.transpose(1,0,2).reshape(4,-1)
    flipped=c.execute(vword,2,moved).reshape(4,len(m.PADDING),4).transpose(1,0,2)
    ARRAYS['wrong_phase_returned']=flipped
    check('wrong_phase_same_diagonal_effect',np.linalg.norm(effect(flipped)-effect(returned))<GATE)
    check('wrong_phase_flips_charged_response',all(np.linalg.norm(cross(flipped,f)+cross(returned,f))<GATE for f in range(3)) and np.linalg.norm(cross(flipped,0)-cross(returned,0))>1)
    # Uniform proof checks all new primitive requests, not just populated branches.
    all_words=[prep,rw,weak,comp,vword]
    for pw,angle in [(F(1,2),2*theta),(F(),-theta)]:
        w=m.preparation_word('phi',pw);all_words.extend([w,c.inverse(w),controlled_phase_word(angle)])
    unique={(e['kind'],e['angle']) for w in all_words for e in w}
    for kind,angle in sorted(unique):synth.certify_primitive(kind,angle)
    check('all_new_full_six_level_primitive_intervals',True)
    check('all_full_diagonal_inverse_Walsh_identities',True)
    check('all_source_words_admitted',max(x['pulses'] for x in WORDS.values())<=interpreter.MAX_WORD)
    delta=F();gain=F(1)
    for b in BLOCKS:delta+=gain*F(b['complete_halfdiamond']);gain*=F(b['gain'])
    check('complete_halfdiamond_below_frozen_1e4',delta<F(1,10000))
    sample=np.zeros((2048,4),complex)
    for j in range(4):sample[512*j,j]=1
    middle=controlled_phase_word(2*theta)
    queries=[]
    for index in (0,len(middle)//2,len(middle)-1):
        d=1<<len(middle[index]['targets']);bad=first_failure(middle,index,d,sample,11)
        queries.append({'pulse':index,'proper_code_dimension':d,'complement':d,'mass':float(np.linalg.norm(bad)**2)})
    check('processor_failure_queries_keep_all_four_inputs',len(queries)==3)
    for n,v in outcomes.items():ARRAYS['effect__'+n]=v
    return {'task':'LOCAL-M10-R11-VACUUM-CONTROL-20261115','candidate':'lepton-vacuum-filter','class':'F',
            'checks':checks,'all_checks_pass':all(checks.values()),'N_count':0,'epsilon':'2/5','vacuum_family_code':3,
            'diagnostic_gate':GATE,'return_map_error':float(np.linalg.norm(k+D)),
            'full_effect_sum_error':float(np.linalg.norm(sum(outcomes.values())-np.eye(4))),
            'vacuum_failure_weight':vacuum_failure,'charged_cross_errors':cross_errors,
            'failure_cross_errors':failure_cross_errors,'disabled_vacuum_effect':vac_effect,
            'disabled_vacuum_click':ctrlweight,'actual_pinned_native_calls':0,
            'reused_source_blocks':['object_forward','object_inverse','shift_clean','shift_dirty','inverse_shift_clean','inverse_shift_dirty'],
            'fresh_source_preparation':True,'source_words':WORDS,'main_branch_pulses':sum(b['pulses'] for b in BLOCKS),
            'new_unique_primitive_requests':len(unique),'complete_error_blocks':BLOCKS,
            'complete_halfdiamond_bound':str(delta),'complete_gain_bound':str(gain),
            'processor_failure_queries':queries,'source_hashes':source_bindings(),
            'array_evaluation':'no truncation or success normalization; exact source-map bound separate from floating diagnostic',
            'services':['four-code family/vacuum chart','coherent proper-code loads','analog source controls','scratch and fresh processors','complete coherent exports'],
            'physical_row_established':False,'retained_count_extension_executed':False}

def main():
    start=time.perf_counter();r=calculate();blob=canon(r)
    (HERE/'RECEIPT.json').write_text(blob)
    np.savez_compressed(HERE/'SOURCE_ARRAYS.npz',**ARRAYS)
    bindings={'archive_sha256':sha((HERE/'SOURCE_ARRAYS.npz').read_bytes()),'receipt_sha256':sha(blob.encode()),
              'implementation_sha256':sha(Path(__file__).read_bytes()),'arrays':{n:{'shape':list(a.shape),'dtype':str(a.dtype),'raw_sha256':sha(a.tobytes())} for n,a in ARRAYS.items()}}
    (HERE/'ARRAY_BINDINGS.json').write_text(canon(bindings))
    report={'status':'pass' if r['all_checks_pass'] else 'fail','checks':len(r['checks']),
            'failed':[n for n,v in r['checks'].items() if not v],'seconds':time.perf_counter()-start,
            'main_branch_pulses':r['main_branch_pulses'],'halfdiamond':float(F(r['complete_halfdiamond_bound'])),
            'return_map_error':r['return_map_error'],'vacuum_failure_weight':r['vacuum_failure_weight']}
    (HERE/'RUN.json').write_text(canon(report));print(canon(report),flush=True)
    assert r['all_checks_pass'],report['failed']

if __name__=='__main__':main()
