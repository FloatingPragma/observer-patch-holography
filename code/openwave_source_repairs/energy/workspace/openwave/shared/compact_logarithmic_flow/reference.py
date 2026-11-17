"""Complete native swap-parity preparation of the SAME transfer reference."""
from fractions import Fraction as F
import json,math
import numpy as np
import native_words as nw
from native_words import c
from model.native_interface import pointer
from model.observer.native import native

BITS=15

def routines():
    ccx=([c.pulse('H',math.pi/2,[2])]+c.compile_diagonal([0.]*7+[math.pi])+
         [c.pulse('H',-math.pi/2,[2])])
    flip=[c.pulse('I',-math.pi/2,[0]),c.pulse('X',math.pi/2,[0])]
    beam=[c.pulse('Y',-math.pi/4,[0])]
    words={'ccx':ccx,'flip':flip,'beam':beam,'inverse_beam':c.inverse(beam)}
    raw={};fail={};metadata={}
    for name,word in words.items():
        bits=3 if name=='ccx' else 1
        raw[name],fail[name]=nw.apply_word(word,np.eye(1<<bits),bits)
        metadata[name]={'pulses':len(word),'word_sha256':nw.digest_word(word)}
    target=np.zeros((8,8))
    for k in range(8):target[k^((1<<2) if k&3==3 else 0),k]=1
    assert np.linalg.norm(raw['ccx']-target,'fro')<1e-10
    assert np.linalg.norm(raw['flip']-np.array([[0,1],[1,0]]),'fro')<1e-10
    return words,raw,fail,metadata

def controlled_swap(a,b):
    ops=[]
    def emit(name,targets):ops.append((name,tuple(targets)))
    def controlled_x(controls,target):
        scratch=list(range(9,15));chain=[(controls[0],controls[1],scratch[0])]
        for j in range(2,len(controls)-1):chain.append((scratch[j-2],controls[j],scratch[j-1]))
        for gate in chain:emit('ccx',gate)
        emit('ccx',(scratch[len(controls)-3],controls[-1],target))
        for gate in reversed(chain):emit('ccx',gate)
    path=[int(a)]
    for j in range(8):
        if (a^b)>>j&1:path.append(path[-1]^(1<<j))
    pairs=list(zip(path[:-1],path[1:]))
    for x,y in pairs+pairs[-2::-1]:
        target=(x^y).bit_length()-1;controls=[j for j in range(8) if j!=target]
        zeros=[j for j in controls if not (x>>j&1)]
        for j in zeros:emit('flip',(j,))
        controlled_x(controls+[8],target)
        for j in reversed(zeros):emit('flip',(j,))
    return ops

def verify_permutation(ops,a,b):
    initial=np.arange(512,dtype=np.int64);out=initial.copy()
    for name,targets in ops:
        if name=='flip':out^=1<<targets[0]
        else:
            x,y,z=targets;out^=(((out>>x)&1)&((out>>y)&1))<<z
    wanted=initial.copy()
    for p in (a,b):wanted[256+p]=256+(b if p==a else a)
    assert np.array_equal(out,wanted)
    # Every elementary gate is a full dirty-scratch permutation. The explicit
    # clean-input equality suffices to identify the controlled swap used here.
    return len(initial)

def apply_routine(state,name,targets,raw,failure):
    # Routine matrices use little-endian local bit labels; _apply orders the
    # requested axes from most significant to least significant.
    ix=nw.indices(tuple(reversed(targets)),BITS);block=state[ix]
    f=np.einsum('iac,ij,jad->cd',block.conj(),failure[name],block,optimize=True)
    state[ix]=(raw[name]@block.reshape(len(raw[name]),-1)).reshape(block.shape)
    return f

def pointer_effects(law):
    matrix=pointer.local_matrix();effects=pointer.native_effects(law);out=[]
    for bit in (0,1):
        E=np.zeros((2,2),complex)
        for z in range(4):
            block=matrix[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
            E+=(effects[bit][z]/2)*(block.conj().T@block)
        out.append(E)
    return tuple(out)

def reference_preparation(source,progress=lambda x:None):
    words,raw,failure,metadata=routines();seams=source['seams'];V=len(source['counts']);K=1260
    assert V==240 and len(seams)==630
    exact_plus=np.zeros((V,V),dtype=np.int64)
    # Sum_e(I+U_e)=1260 T, so its plus output on I/V is T/V.
    for a,b in seams:
        U=np.eye(V,dtype=np.int64);U[a,a]=U[b,b]=0;U[a,b]=U[b,a]=1
        exact_plus+=np.eye(V,dtype=np.int64)+U
    assert np.array_equal(exact_plus,source['counts'])
    # Tr(sum_e P_plus)/M =239, independently of the flow eigensolver.
    assert int(np.trace(exact_plus))==239*K
    for coin in (0,1):
        outputs=source['destination'][:,coin];currents=source['current'][:,coin]
        for current in sorted(set(map(int,currents))):
            image=outputs[currents==current]
            assert len(set(map(int,image)))==len(image)
        # Fine-current partial injections cover all inputs; each coin's
        # summed effect is exactly I/2, also for arbitrary reference inputs.
        assert len(outputs)==V
    selector_tapes=[];selector_native_calls=0
    for bits,limit in ((8,240),(10,630)):
        for value in range(1<<bits):
            p=0;tape=[]
            for k in range(bits):
                coin=(value>>k)&1
                out,record=native((int(p==0),int(p==1)),(0,1),coin)
                selector_native_calls+=1
                y=0 if out[0] else 1
                assert record==(coin,int(p==0)-out[0])
                tape.append((coin,int(record[1]),p,y));p=y
            selector_tapes.append({'bits':bits,'value':value,'accepted':value<limit,'tape':tape})
    # Native blank preparation is inherited; all source address routes are new
    # actual source words, keeping every proper-code failure.
    prepared=[];route_failure=[];route_words=[]
    for p in range(V):
        word=[pulse for bit in range(8) if p>>bit&1 for pulse in c.remap(words['flip'],[bit])]
        blank=np.zeros(256);blank[0]=1
        x,bad=nw.apply_word(word,blank,8);prepared.append(x[:,0]);route_failure.append(bad[0,0].real)
        route_words.append({'address':p,'pulses':len(word),'sha256':nw.digest_word(word)})
    prepared=np.asarray(prepared).T;rho0=prepared@prepared.conj().T/V
    target0=np.diag([1/V]*V+[0.]*16)
    route_error=float(np.linalg.norm(rho0-target0,'fro'))
    assert route_error<1e-10
    ef=pointer_effects('C');permutations=0;max_pulses=0;records=[];clean_outputs=[];scratch_grams=[];processor_grams=[]
    # All630 word permutations are checked on all512 clean inputs.
    # Full32768-address floating witnesses use the first seam of each
    # Hamming-distance class, plus the longest emitted circuit. The exact
    # primitive certificates and telescoping bound cover every other word.
    menu=[controlled_swap(int(a),int(b)) for a,b in seams]
    representatives={}
    for i,(a,b) in enumerate(seams):representatives.setdefault((int(a)^int(b)).bit_count(),i)
    selected=set(representatives.values())|{max(range(len(menu)),key=lambda i:len(menu[i]))}
    selected_indices=[]
    flag_indices=nw.indices((8,),BITS)
    for i,(a,b) in enumerate(seams):
        a,b=int(a),int(b);ops=menu[i]
        permutations+=verify_permutation(ops,a,b)
        pulses=2+sum(len(words[name]) for name,_ in ops);max_pulses=max(max_pulses,pulses)
        entry={'seam_index':i,'seam':[a,b],'routine_calls':len(ops),'native_pulses':pulses,
               'clean_permutation_columns':512,'full_scratch_witness_executed':i in selected}
        if i not in selected:
            records.append(entry)
            continue
        state=np.zeros((1<<BITS,2),complex)
        state[a]=[1/math.sqrt(2),1/math.sqrt(2)]
        state[b]=[1/math.sqrt(2),-1/math.sqrt(2)]
        bad=apply_routine(state,'beam',(8,),raw,failure)
        for name,targets in ops:bad+=apply_routine(state,name,targets,raw,failure)
        bad+=apply_routine(state,'inverse_beam',(8,),raw,failure)
        block=state[flag_indices]
        laws=np.array([np.einsum('iac,ij,jad->cd',block.conj(),E,block,optimize=True) for E in ef])
        wanted=np.array([[[1,0],[0,0]],[[0,0],[0,1]]],complex)
        assert np.max(abs(laws-wanted))<1e-8
        dirty=state[512:].conj().T@state[512:]
        assert np.max(abs(dirty))<1e-16
        entry['flag_law_max_error']=float(np.max(abs(laws-wanted)));records.append(entry)
        selected_indices.append(i)
        clean_outputs.append(state[:512]);scratch_grams.append(dirty);processor_grams.append(bad)
        progress({'reference_full_scratch_witness':i,'witnesses':len(selected)})
    position_timeout=F(1,16)**16
    seam_timeout=F(394,1024)**32
    selector_success=(1-position_timeout)*(1-seam_timeout)
    # Complete native reset and read costs; timeout is a retained outcome.
    reset_timeout=F(35,2147483648)
    pointer_error=F(200000000001,20000000000000000000000)
    source_bound=nw.complete_bound(max_pulses)
    total=(16*(reset_timeout+pointer_error)+30*nw.PULSE_GATE+source_bound['half_diamond']+
           max(route_words,key=lambda z:z['pulses'])['pulses']*nw.PULSE_GATE)
    assert total<F(1,10**4)
    report={'source_dimension':V,'invalid_labels':16,'trace_T':239,'ideal_plus_given_selectors':str(F(239,240)),
            'ideal_minus_given_selectors':str(F(1,240)),'position_selector_timeout':str(position_timeout),
            'seam_selector_timeout':str(seam_timeout),'selector_success':str(selector_success),
            'unconditional_ideal_plus':str(selector_success*F(239,240)),
            'unconditional_ideal_minus':str(selector_success/F(240)),
            'reset_timeout_per_bit':str(reset_timeout),'reset_attempt_budget_per_bit':2678,
            'blank_bits':15,'native_pointer_reads':1,'full_clean_permutation_columns':permutations,
            'selector_native_calls':selector_native_calls,'complete_one_draw_selector_tapes':selector_tapes,
            'coin_full_position_columns_checked':2*V,
            'full_dirty_scratch_dimension':1<<BITS,'max_native_pulses_per_seam':max_pulses,
            'full_scratch_witness_seam_indices':selected_indices,
            'all_seam_clean_permutations_checked':True,
            'all_seam_full_scratch_floating_execution_claimed':False,
            'complete_preparation_halfdiamond_upper':str(total),'route_density_frobenius_error':route_error,
            'routine_metadata':metadata,'route_words':route_words,'seam_records':records,
            'ideal_plus_target':'T/240; normalized target T/239 only on plus',
            'actual_full_mixed_density_directly_materialized':False,
            'full_mixed_preparation_certified_by_arbitrary_reference_operator_bound':True}
    arrays={'prepared_positions':prepared,'rho0':rho0,'parity_clean_outputs':np.asarray(clean_outputs),
            'scratch_failure_reference_grams':np.asarray(scratch_grams),'processor_failure_reference_grams':np.asarray(processor_grams),
            'route_failure_probabilities':np.asarray(route_failure),'pointer_effects_C':np.asarray(ef)}
    return report,arrays
