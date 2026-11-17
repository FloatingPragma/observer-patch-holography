"""A separately specified 144-object chart, using existing native pulse services."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import sys, math, json, hashlib, argparse
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=ROOT/'reverse-engineering-reality'
sys.path.insert(0,str(SOURCE))
import numpy as np
import sympy as sy
from model.observer.source_scalar import compiler as c
from model.native_effect_synthesis import synthesis as synth, intervals as iv
from model.native_program import interpreter
from model.openwave_coherent_native_architecture import source as arch
from model.lepton_register_round import register as reg, geometry
from model.finite_family_clock import family, primitives as fp

CODES=tuple(ic+16*ie for ic in range(12) for ie in range(12))
PADDING=tuple(sorted(set(range(512))-set(CODES)))
PAD=511
EPS=F(2,5)
SURVIVAL=1-EPS
NUM_GATE=3e-8

def sha(value):
    return hashlib.sha256(value).hexdigest()

def state_record(a):
    a=np.asarray(a,complex)
    return {'shape':list(a.shape),'real':a.real.tolist(),'imag':a.imag.tolist()}

def weights_of(vector):
    out=[]
    for x in vector:
        square=sy.simplify(x*sy.conjugate(x))
        assert square.is_Rational and sy.im(x)==0 and x.is_nonnegative is True
        out.append(F(int(sy.numer(square)),int(sy.denom(square))))
    assert sum(out)==1
    return tuple(out)

def rays():
    obj=fp.state(sy.sqrt(sy.Rational(1,3)),0,sy.sqrt(sy.Rational(2,3)))
    prep=sy.kronecker_product(obj,obj)
    phi=fp.analyzer()
    flags=fp.flags()
    hold=sy.kronecker_product(flags['C_idle'],flags['E_agree'])
    return {'prep':prep,'phi':phi,'hold':hold}

def embedding(vector):
    out=np.zeros(512,complex)
    out[list(CODES)]=np.array(vector.evalf(18).tolist(),complex).ravel()
    return out

@lru_cache(None)
def preparation_word(name,padding_weight):
    vector=rays()[name]; weights=[F()]*512
    for address,q in zip(CODES,weights_of(vector),strict=True):weights[address]=(1-padding_weight)*q
    weights[PAD]=padding_weight
    # Every relative phase of these actual source rays is zero. No diagonal
    # preparation word is required; all amplitude-tree zero requests remain.
    word=synth._state_word(weights)
    interpreter.source_word(word,9)
    assert len(word)==12291
    return word

@lru_cache(None)
def ray_phase(name,padding_weight,angle):
    prep=preparation_word(name,padding_weight)
    phases=[0]*512;phases[0]=angle
    word=c.inverse(prep)+synth._phase_word(phases)+prep
    interpreter.source_word(word,9)
    assert len(word)==35852
    return word

def filter_words(name):
    theta=iv.split_angle(SURVIVAL**2,1-SURVIVAL**2)
    return (ray_phase(name,F(1,2),2*theta),ray_phase(name,F(),-theta))

def return_word():
    return ray_phase('prep',F(1,2),2*iv.split_angle(F(),F(1)))

def flip(bit):
    return [c.pulse('I',-math.pi/2,[bit]),c.pulse('X',math.pi/2,[bit])]

def hadamard(bit):
    return [c.pulse('I',-math.pi/2,[bit]),c.pulse('H',math.pi/2,[bit])]

def phase_minus(bits,at):
    phases=[0.]*(1<<len(bits));phases[at]=math.pi
    return c.remap(c.compile_diagonal(phases),bits)

def object_word():
    # Object bits: orientation0, record1, low load2, high load3.
    o,r,l,z=0,1,2,3
    u=(hadamard(l)+hadamard(r)+phase_minus([l],1)+phase_minus([r],1)
       +phase_minus([l,r],3)+c.cnot(o,l))
    q=c.inverse(u)+phase_minus([l,r,z],0)+u
    base=flip(o)+q+phase_minus([r,o],3)
    # On load code3, Q is already identity, while the address/record words
    # would give V R2. Undo those two operations only on that padding load.
    padding_inverse=(phase_minus([o,r,l,z],15)
                     +c.remap(arch.primitive_routines()['ccx'],[l,z,o]))
    return base+padding_inverse

def shift_program():
    permutation=[]
    for index in range(64):
        ic=index&15;f=index>>4
        target=(f+reg.step_of(ic))%3 if ic<12 and f<3 else f
        permutation.append(ic+16*target)
    return arch.compile_permutation(tuple(permutation)),tuple(permutation)

def source_block(word,bits):
    interpreter.source_word(word,bits)
    return c.execute(word,bits)

def shift_blocks():
    program,permutation=shift_program()
    truth=arch.execute_basis(program)
    assert tuple(map(int,truth[:64]))==permutation
    q=program['address_bits']+program['scratch_bits']
    routines=arch.primitive_routines()
    matrices={name:source_block(routines[name],1 if name=='flip' else 3)
              for name in ('flip','ccx')}
    state=np.zeros((1<<q,64),complex);state[:64]=np.eye(64)
    pulse_count=0
    for name,targets in program['gates']:
        state=c._apply(state,matrices[name],list(reversed(targets)),q)
        pulse_count+=len(routines[name])
    expected=np.zeros_like(state);expected[list(permutation),np.arange(64)]=1
    return state[:64],state[64:],{'address_bits':6,'scratch_bits':program['scratch_bits'],
          'all_dirty_addresses_checked':len(truth),'gates':len(program['gates']),
          'native_pulses':pulse_count,'full_dirty_permutation_sha256':sha(truth.tobytes()),
          'clean_source_error':float(np.linalg.norm(state-expected)),
          'clean_source_max_abs_error':float(np.max(abs(state-expected)))},routines

def to_global(a):
    # a[object512,family4,input_family3], little-endian object then family.
    return a.transpose(1,0,2).reshape(2048,-1)

def from_global(a):
    return a.reshape(4,512,-1).transpose(1,0,2)

def apply_object_block(a,block,bits):
    return from_global(c._apply(to_global(a),block,list(reversed(bits)),11))

def apply_shift(a,clean,dirty):
    g=to_global(a)
    targets=[10,9,3,2,1,0]
    axes=[11-1-j for j in targets]
    axes += [i for i in range(11) if i not in axes]
    axes += [11]
    arranged=g.reshape([2]*11+[g.shape[1]]).transpose(axes)
    columns=arranged.reshape(64,-1)
    failed=(dirty@columns).reshape(-1,g.shape[1])
    out=(clean@columns).reshape(arranged.shape).transpose(np.argsort(axes)).reshape(g.shape)
    return from_global(out),failed

def apply_tick(a,block,clean,dirty,inverse=False):
    failures=[]
    if not inverse:
        a=apply_object_block(a,block,[0,1,2,3])
        a=apply_object_block(a,block,[4,5,6,7])
        a,failed=apply_shift(a,clean,dirty);failures.append(failed)
    else:
        a,failed=apply_shift(a,clean,dirty);failures.append(failed)
        a=apply_object_block(a,block,[4,5,6,7])
        a=apply_object_block(a,block,[0,1,2,3])
    return a,failures

def word_export(word,a):
    flat=a.reshape(512,-1)
    actual=c.execute(word,9,flat).reshape(a.shape)
    padding=actual[list(PADDING)].copy()
    success=np.zeros_like(actual);success[list(CODES)]=actual[list(CODES)]
    return success,padding

def effect(a):
    flat=np.asarray(a).reshape(-1,a.shape[-1])
    return flat.conj().T@flat

def exact_effects():
    # Contract the actual native432 tick independently of the detector words.
    parts=reg.compression()
    M=sy.zeros(3)
    for a,b in parts.values():
        M+=(sy.Matrix(a.tolist())+sy.sqrt(2)*sy.Matrix(b.tolist()))/(96*sy.sqrt(3))
    M=M.applyfunc(sy.simplify)
    R=fp.shift()
    expected=sy.eye(3)/(3*sy.sqrt(3))+(R+R.T)/(3*sy.sqrt(6))
    assert (M-expected).applyfunc(sy.simplify)==sy.zeros(3)
    A=sy.simplify(M.T*M)
    return M,A

def pulse_bound(length,target_units=0):
    x=length*synth.PULSE_GATE
    error=x/(1-x)+F(target_units,iv.AQ)
    failure=length*synth.PULSE_GATE**2/(1-x)**2
    delta=(error*(2+error)+failure)/2
    gain=1/(1-x)**2+failure
    return dict(pulses=length,target_grid_error=str(F(target_units,iv.AQ)),
                operator_error=str(error),first_failure_mass=str(failure),
                complete_halfdiamond=str(delta),gain=str(gain))

def primitive_failures(word,pulse,complement,initial):
    """One actual first-failure query on the full new object chart.

    This returns all spectators/references; only the addressed proper code
    resets to zero. The whole prefix is evaluated, without success division.
    """
    a=np.array(initial,complex,copy=True).reshape(512,-1)
    assert 0<=pulse<len(word)
    for k,event in enumerate(word[:pulse+1]):
        full=c.native_pulse(event['kind'],event['angle']);d=1<<len(event['targets'])
        if k==pulse:
            assert d<=complement<6
            block=np.zeros((d,d),complex);block[0]=full[complement,:d]
        else:block=full[:d,:d]
        a=c._apply(a,block,event['targets'],9)
    return a

def calculate():
    checks={};arrays={};all_words=[]
    checks['coherent_native_360_columns']=family.native_event_check()['exact_mismatches']==0
    ow=object_word();oi=c.inverse(ow);all_words += [ow,oi]
    B=source_block(ow,4);Bi=source_block(oi,4)
    target=np.eye(16,dtype=complex)
    for seam in (geometry.C_SEAM,geometry.E_SEAM):
        _,_,twow,v=reg.object_operators(seam)
        wanted=(v@twow)/2
        target[:12,:12]=wanted
        checks['object_word_'+str(seam)]=bool(np.linalg.norm(B-target)<NUM_GATE)
    checks['object_inverse_word']=bool(np.linalg.norm(Bi-target.T)<NUM_GATE)
    arrays['object_forward']=B;arrays['object_inverse']=Bi
    shift,dirty,shift_meta,routines=shift_blocks()
    # The inverse is recompiled from the inverse permutation, not a returned
    # floating matrix adjoint substituted into the execution.
    program,perm=shift_program();inverse_perm=tuple(perm.index(i) for i in range(64))
    ip=arch.compile_permutation(inverse_perm)
    truth=arch.execute_basis(ip);assert tuple(map(int,truth[:64]))==inverse_perm
    qi=ip['address_bits']+ip['scratch_bits'];v=np.zeros((1<<qi,64),complex);v[:64]=np.eye(64)
    local={name:source_block(routines[name],1 if name=='flip' else 3) for name in ('flip','ccx')}
    for name,targets in ip['gates']:v=c._apply(v,local[name],list(reversed(targets)),qi)
    invshift,invdirty=v[:64],v[64:]
    checks['controlled_shift_full_source_word']=shift_meta['clean_source_error']<NUM_GATE
    itarget=np.zeros_like(v);itarget[list(inverse_perm),np.arange(64)]=1
    checks['controlled_inverse_full_source_word']=bool(np.linalg.norm(v-itarget)<NUM_GATE)
    arrays['shift_clean']=shift;arrays['shift_dirty']=dirty
    arrays['inverse_shift_clean']=invshift;arrays['inverse_shift_dirty']=invdirty
    all_words += [routines['flip'],routines['ccx']]
    prepword=preparation_word('prep',F());all_words.append(prepword)
    blank=np.zeros(512,complex);blank[0]=1
    prepared=c.execute(prepword,9,blank)
    idealprep=embedding(rays()['prep'])
    checks['actual_preparation_word']=bool(np.linalg.norm(prepared-idealprep)<NUM_GATE)
    prep_padding=prepared[list(PADDING)]
    prepared_success=np.zeros(512,complex);prepared_success[list(CODES)]=prepared[list(CODES)]
    a=np.zeros((512,4,3),complex)
    for j in range(3):a[:,j,j]=prepared_success
    after,forward_fail=apply_tick(a,B,shift,dirty)
    # Independent original432-state source comparison.
    p=idealprep[list(CODES)]
    native=reg.full_tick()/4
    oracle=native@np.kron(p[:,None],np.eye(3))
    actual=np.array([after[code,f,:] for code in CODES for f in range(3)])
    checks['same_full_432_tick']=bool(np.linalg.norm(actual-oracle)<NUM_GATE)
    M,A=exact_effects();An=np.array(A.evalf(18).tolist(),complex);Mn=np.array(M.evalf(18).tolist(),complex)
    arrays['family_after_tick']=after
    outputs={}
    inverse_pulses=sum(len(routines[name]) for name,targets in ip['gates'])
    for label,name,compensate in [('native_phi','phi',True),('no_compensation','phi',False),('hold_control','hold',True)]:
        words=filter_words(name);all_words.extend(words)
        current,click=word_export(words[0],after)
        if compensate:current,corrpad=word_export(words[1],current)
        else:corrpad=np.zeros((len(PADDING),4,3),complex)
        current,inverse_fail=apply_tick(current,Bi,invshift,invdirty,True)
        rw=return_word();all_words.append(rw)
        nonreturn,returned=word_export(rw,current)
        result={'first_click':click,'correction_padding':corrpad,'destructive_return':returned,
                'nonreturn':nonreturn,'inverse_scratch_failure':inverse_fail[0]}
        for key,value in result.items():arrays[label+'__'+key]=value
        effects={key:effect(value) for key,value in result.items()}
        outputs[label]={key:state_record(value) for key,value in effects.items()}
        # Complete successful processor branches retain all object/family and
        # scratch outcomes. Processor failures have the separate complete bound.
        total=sum(effects.values())+effect(forward_fail[0])+float(np.vdot(prep_padding,prep_padding).real)*np.eye(3)
        checks[label+'_full_coarse_effect_sum']=bool(np.linalg.norm(total-np.eye(3))<NUM_GATE)
        return_address=PADDING.index(PAD)
        k=returned[return_address,:3,:]
        outputs[label]['return_map']=state_record(k)
        if label=='native_phi':
            expected=-np.eye(3)+float(EPS)*An
            checks['actual_return_map']=bool(np.linalg.norm(k-expected)<NUM_GATE)
            checks['actual_click_effect']=bool(np.linalg.norm(effects['first_click']-float(1-SURVIVAL**2)*An)<NUM_GATE)
            checks['actual_nonreturn_effect']=bool(np.linalg.norm(effects['nonreturn']-float(EPS**2)*(An-An@An))<NUM_GATE)
            checks['actual_correction_padding_zero_target']=bool(np.linalg.norm(effects['correction_padding'])<NUM_GATE)
            checks['actual_extra_family_code_retained']=bool(sum(np.linalg.norm(value[:,3,:])**2 for key,value in result.items() if value.ndim==3 and value.shape[1]==4)<NUM_GATE)
        elif label=='no_compensation':
            raw=complex(-16/25,12/25)
            checks['omitted_compensation_map']=bool(np.linalg.norm(k+np.eye(3)+raw*An)<NUM_GATE)
            checks['omitted_compensation_detected']=bool(np.linalg.norm(k.imag)>1e-3)
        else:
            checks['changed_analyzer_detected']=bool(np.linalg.norm(effects['first_click']-float(1-SURVIVAL**2)*An)>1e-3)
    # Direct non-family probes of the emitted Phi weak word.
    phi=embedding(rays()['phi']);overlap=next(i for i in CODES if abs(phi[i])>0)
    orth=next(i for i in CODES if phi[i]==0)
    probes=np.stack([phi,np.eye(512)[:,overlap],np.eye(512)[:,orth],(phi+1j*np.eye(512)[:,orth])/math.sqrt(2)],axis=1)
    probe_shape=(512,1,4)
    weak,correction=filter_words('phi')
    y,abort=word_export(weak,probes.reshape(probe_shape));y,abort2=word_export(correction,y)
    targetprobe=probes-float(EPS)*phi[:,None]*(phi.conj()@probes)[None,:]
    checks['full_chart_reference_probe_filter']=bool(np.linalg.norm(y.reshape(512,4)-targetprobe)<NUM_GATE)
    checks['full_chart_probe_complete_probability']=bool(np.linalg.norm(effect(y)+effect(abort)+effect(abort2)-probes.conj().T@probes)<NUM_GATE)
    arrays['probe_survivor']=y;arrays['probe_first_click']=abort;arrays['probe_correction_padding']=abort2
    arrays['forward_scratch_failure']=forward_fail[0];arrays['preparation_padding']=prep_padding
    unique={(x['kind'],x['angle']) for w in all_words for x in w}
    for kind,angle in unique:synth.certify_primitive(kind,angle)
    checks['all_emitted_primitive_interval_checks']=True
    # Every complete source word stays within the unchanged grammar.
    admitted=[interpreter.source_word(w,9) for w in all_words]
    checks['all_words_admitted']=all(x['schema']==interpreter.WORD_SCHEMA for x in admitted)
    queried=[]
    for pulse in (0,len(weak)//2,len(weak)-1):
        d=1<<len(weak[pulse]['targets'])
        f=primitive_failures(weak,pulse,d,probes)
        queried.append({'pulse':pulse,'complement':d,'original_mass_sum':float(np.linalg.norm(f)**2)})
    checks['failure_query_retains_reference']=all(x['original_mass_sum']>=0 for x in queried)
    sources={}
    for mod in list(sys.modules.values()):
        path=getattr(mod,'__file__',None)
        if path:
            q=Path(path).resolve()
            if q.suffix=='.py' and q.is_relative_to(SOURCE):sources[str(q.relative_to(ROOT))]=sha(q.read_bytes())
    for rel in ['model/native_coexistent_effects/PROOF.md','model/native_effect_synthesis/PROOF.md','model/finite_family_clock/PROOF.md','model/lepton_register_round/PROOF.md']:
        q=SOURCE/rel;sources[str(q.relative_to(ROOT))]=sha(q.read_bytes())
    blocks=[pulse_bound(len(prepword),9),pulse_bound(2*len(ow),1),pulse_bound(shift_meta['native_pulses'],1),
            pulse_bound(len(weak),20),pulse_bound(len(correction),19),pulse_bound(inverse_pulses,1),
            pulse_bound(2*len(oi),1),pulse_bound(len(return_word()),20)]
    delta=F();gain=F(1)
    for block in blocks:delta+=gain*F(block['complete_halfdiamond']);gain*=F(block['gain'])
    checks['complete_bound_below_fixed_gate']=delta<F(1,10000)
    result={'task':'LOCAL-M10-R11-FAMILY-FILTER-20261114','candidate':'lepton-encoded-filter','class':'F',
        'checks':checks,'all_checks_pass':all(checks.values()),'diagnostic_norm_gate':NUM_GATE,
        'chart':{'codes':list(CODES),'all_work_addresses':512,'padding_addresses':list(PADDING),'selected_padding':PAD,
                 'family_codes':[0,1,2],'extra_family_code':3,'new_chart_supplied':True,'original_H66_interface_changed':False},
        'source_coherent_columns':360,'actual_pinned_classical_native_calls':0,
        'object_word_pulses':len(ow),'controlled_shift':shift_meta,'inverse_shift_native_pulses':inverse_pulses,
        'preparation_pulses':len(prepword),'rank_one_phase_word_pulses':len(weak),
        'unique_native_primitive_requests':len(unique),'source_hashes':dict(sorted(sources.items())),
        'original_family_effect':[[str(x) for x in row] for row in A.tolist()],
        'original_family_compression':[[str(x) for x in row] for row in M.tolist()],
        'coarse_outcome_effects':outputs,'queried_processor_first_failures':queried,
        'complete_error_blocks':blocks,'complete_halfdiamond_bound':str(delta),'complete_gain_bound':str(gain),
        'physical_row_established':False,'N_count_ancillas_executed':0,
        'declared_services':['coherent proper-code loads and gluing','blank processors and scratch','full coherent chart export','analog source requests and Kraus trace semantics'],
        'output_scope':'all chart padding and extra family outputs retained; processor first failures lazy and uniformly bounded; no success renormalization'}
    return result,arrays

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result,arrays=calculate()
    text=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n'
    path=HERE/'RECEIPT.json'
    if args.check:
        assert path.read_text()==text,'full deterministic receipt differs'
        saved=np.load(HERE/'SOURCE_ARRAYS.npz')
        for k,v in arrays.items():assert np.array_equal(saved[k],v),k
    else:
        path.write_text(text);np.savez_compressed(HERE/'SOURCE_ARRAYS.npz',**arrays)
    assert result['all_checks_pass'],[k for k,v in result['checks'].items() if not v]
    print(json.dumps({'status':'pass','checks':len(result['checks']),'complete_halfdiamond_bound':result['complete_halfdiamond_bound'],
                      'physical_row_established':False},sort_keys=True))

if __name__=='__main__':main()
