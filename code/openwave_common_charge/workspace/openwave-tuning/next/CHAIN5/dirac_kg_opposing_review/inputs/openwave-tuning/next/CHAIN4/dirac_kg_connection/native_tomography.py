"""Failure-inclusive native tomography of the actual free Fourier factor."""
from connection import *
from collections import defaultdict
from itertools import combinations

def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,default=str)+'\n');print(name,'PASS',flush=True)

def prep_word(a,b=None,phase=0):
    out=[]
    if b is not None:
        different=a^b;pivot=(different&-different).bit_length()-1
        out+=list(words.remap(native.primitive_templates()['hadamard'],(pivot,)))
        out.append(('Z',(pivot,),native.angle(F(-phase,4))))
        for bit in range(2):
            if bit!=pivot and different>>bit&1:out+=list(words.remap(native.cnot(0,1),(pivot,bit)))
    for bit in range(2):
        if a>>bit&1:out+=list(words.remap(native.primitive_templates()['flip'],(bit,)))
    return out

def pointer():
    from model.observer import native as event_source
    from model.observer.source_discrimination import circuit
    from model.native_interface import pointer as ptr
    tape=events(words.pointer_word(2,0,1));M,failed,gram,bound=execute(tape,3,np.eye(8))
    chosen=M[:,[0,4]];result={}
    for law in ('C','E'):
        groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)));metadata={}
        for code in range(4):
            K=chosen[[code,code+4]];source=circuit.decode(code)
            for coin in (0,1):
                out,(_,current)=event_source.native(source,ptr.SEAM,coin)
                key=coin,current if law=='C' else source[0]-source[1]
                groups[key][out]+=K/math.sqrt(2)
                metadata[key]={'coin':coin,'label':key[1],'current':current,'bit':int(2*current+2*coin-1<0)}
        result[law]=[(metadata[k],sum((K.conj().T@K for K in groups[k].values()),np.zeros((2,2),complex))) for k in sorted(groups)]
        total=sum((E for _,E in result[law]),np.zeros((2,2),complex))
        assert np.linalg.norm(total+failed[np.ix_([0,4],[0,4])]-np.eye(2))<2e-12
    return result,tape,gram,bound

def run():
    pins=verify();branches,ptape,pgram,pbound=pointer()
    fine_serial={law:[{'record':meta,'effect_real':E.real.tolist(),'effect_imag':E.imag.tolist()}for meta,E in entries] for law,entries in branches.items()}
    preparations=[];tapes={'pointer':ptape};maxprep=F()
    for a in range(4):
        preparations.append((f'basis{a}',a,None,0))
    for a,b in combinations(range(4),2):
        for q in range(4):preparations.append((f'pair{a}{b}_q{q}',a,b,q))
    prepared={};idealprep={};prep_bounds={}
    for name,a,b,q in preparations:
        tape=events(prep_word(a,b,q));M,failed,gram,bound=execute(tape,2,np.eye(4))
        ideal=np.zeros(4,complex);ideal[a]=1
        if b is not None:ideal[b]=1j**q;ideal/=math.sqrt(2)
        v=M[:,0];fidelity=abs(np.vdot(ideal,v));assert abs(fidelity-1)<2e-12
        prepared[name]=v;idealprep[name]=ideal;prep_bounds[name]=bound;maxprep=max(maxprep,bound);tapes['prep_'+name]=tape
    plus_tape=events(words.remap(native.primitive_templates()['hadamard'],(0,)))
    plusM,_,_,plusbound=execute(plus_tape,1,np.eye(2));plus=plusM[:,0]
    analyses={};analysis_bounds={}
    for setting in ('X','Y'):
        word=[] if setting=='X' else [('Z',(0,),native.angle(F(1,4)))]
        word+=list(words.remap(native.primitive_templates()['hadamard'],(0,)))
        tape=events(word);A,_,gram,bound=execute(tape,1,np.eye(2));analyses[setting]=A;analysis_bounds[setting]=bound;tapes['analysis_'+setting]=tape
    flavors={};flavor_bounds={}
    for flavor in range(4):
        word=[]
        for bit in range(2):
            if flavor>>bit&1:word+=list(words.remap(native.primitive_templates()['flip'],(bit,)))
        tape=events(word);M,_,_,bound=execute(tape,2,np.eye(4));flavors[flavor]=M[:,0];flavor_bounds[flavor]=bound;tapes['flavor_'+str(flavor)]=tape
    reset=native.reset_repetition_plan(10,F(1,10**8));timeout=F(reset['complete_timeout_union_bound_exact'])
    records=[];sourcecases=[];maxreconstruction=0.;maxgram=0.;fine_count=0;finehash=hashlib.sha256();zvalues={}
    cases=[(4,(0,0,0)),(4,(1,0,0)),(4,(1,1,0)),(8,(1,0,1))]
    for size,mode in cases:
        h=F(1,size);k=tuple(F(2*n,size) for n in mode);name='L'+str(size)+'_'+''.join(map(str,mode))
        source=list(reduced_matter_word(k,h));controlled_source=list(controlled(source,4));tape=events(controlled_source)
        M,failed,gram,bound=execute(tape,5,np.eye(32));maxgram=max(maxgram,gram);tapes[name]=tape
        ideal=np.eye(32,dtype=complex);targets={}
        for flavor,mass in enumerate((*MASSES,F())):
            U=bloch(np.pi*np.array(list(map(float,k))),h,mass);targets[flavor]=U
            ideal[16+4*flavor:20+4*flavor,16+4*flavor:20+4*flavor]=U
        defect=float(np.linalg.norm(M-ideal,2));assert defect<2e-10
        sourcecases.append({'name':name,'L':size,'mode_integer':mode,'h':str(h),'native_pulses':len(tape),
            'native_all_controlled_columns':32,'source_operator_residual':defect,'native_word_error_exact':str(bound),
            'complete_Gram_residual':gram,'first_export_failure_Gram_trace':float(np.trace(failed).real)})
        for n in (1,2,4):
            Mn=np.linalg.matrix_power(M,n);Fn=np.zeros_like(failed);power=np.eye(32,dtype=complex)
            for _ in range(n):Fn+=power.conj().T@failed@power;power=M@power
            assert np.linalg.norm(Mn.conj().T@Mn+Fn-np.eye(32),2)<5e-10
            ctape=events([('Z',(0,),native.angle(0,0,-h/16))]*n)
            C,Cf,cg,cb=execute(ctape,1,np.eye(2));cv=C@plus
            for flavor in range(4):
                target=np.linalg.matrix_power(targets[flavor],n)
                for prep,_,_,_ in preparations:
                    spin=prepared[prep];initial=np.kron(plus,np.kron(flavors[flavor],spin));v=Mn@initial
                    shape=v.reshape(2,16);rho=shape@shape.conj().T
                    expected=np.vdot(idealprep[prep],target@idealprep[prep]);measured={}
                    for setting in ('X','Y'):
                        after=analyses[setting]@rho@analyses[setting].conj().T
                        for clock_setting in ('X','Y'):
                            clockv=analyses[clock_setting]@cv;clockrho=np.outer(clockv,clockv.conj())
                            allowance=n*bound+cb+prep_bounds[prep]+flavor_bounds[flavor]+2*plusbound+analysis_bounds[setting]+analysis_bounds[clock_setting]+2*pbound+timeout
                            law_records={}
                            for law,entries in branches.items():
                                p=np.array([float(np.trace(E@after).real) for meta,E in entries]);assert min(p)>-1e-13
                                source_contrast=sum(x*(-1 if meta['bit'] else 1) for x,(meta,E) in zip(p,entries))
                                expected_part=expected.real if setting=='X' else expected.imag
                                assert abs(source_contrast-expected_part)<2*float(allowance)+3e-10
                                if law=='C' and clock_setting=='X':measured[setting]=source_contrast
                                for clocklaw,clockentries in branches.items():
                                    cp=np.array([float(np.trace(E@clockrho).real) for meta,E in clockentries])
                                    fine=np.outer(p,cp);fine_count+=fine.size;finehash.update(fine.astype('<f8').tobytes())
                                    clockcontrast=sum(x*(-1 if meta['bit'] else 1) for x,(meta,E) in zip(cp,clockentries))
                                    clocktarget=(math.cos if clock_setting=='X'else math.sin)(float(n*h)/8)
                                    assert abs(clockcontrast-clocktarget)<2*float(allowance)+3e-10
                                    coarse=np.zeros((2,2))
                                    for i,(meta,E) in enumerate(entries):
                                        for j,(cm,CE) in enumerate(clockentries):coarse[meta['bit'],cm['bit']]+=fine[i,j]
                                    # Component law before actual reset-success weights.
                                    assert abs(fine.sum()-np.trace(after).real*np.trace(clockrho).real)<2e-11
                                    law_records[law+clocklaw]=coarse.tolist()
                            records.append({'source':name,'macros':n,'flavor':flavor,'preparation':prep,'setting':setting,'clock_setting':clock_setting,
                                'joint_component_laws':law_records,'complete_native_TV_bound_exact':str(min(F(1),allowance)),
                                'source_first_export_probability':float(np.vdot(initial,Fn@initial).real)})
                    zvalues[name,n,flavor,prep]=complex(measured['X'],measured['Y'])
                reconstructed=np.zeros((4,4),complex)
                for a in range(4):reconstructed[a,a]=zvalues[name,n,flavor,f'basis{a}']
                for a,b in combinations(range(4),2):
                    S=zvalues[name,n,flavor,f'pair{a}{b}_q0']-zvalues[name,n,flavor,f'pair{a}{b}_q2']
                    D=zvalues[name,n,flavor,f'pair{a}{b}_q1']-zvalues[name,n,flavor,f'pair{a}{b}_q3']
                    reconstructed[a,b]=(S-1j*D)/2;reconstructed[b,a]=(S+1j*D)/2
                err=float(np.linalg.norm(reconstructed-target,2));maxreconstruction=max(maxreconstruction,err);assert err<4e-10
        print(name,'complete native tomography passed',flush=True)
    # Scalar source phases are retained by the controlled apparatus.
    scalar=events(controlled([('I',(0,),native.angle(0,0,F(3,8)))],4))
    M,_,_,b=execute(scalar,5,np.eye(32));expected=np.diag([1]*16+[np.exp(3j/8)]*16)
    scalar_error=float(np.linalg.norm(M-expected,2));assert scalar_error<2e-11
    assert abs(np.exp(3j/8)-1)>.3
    # Actual mass and clock mutations, holding all other source terms fixed.
    wrongmass=list(reduced_matter_word((0,0,0),F(1,4),(F(),F(),F())))
    mutated,_,_,_=execute(events(controlled(wrongmass,4)),5,np.eye(32))
    right=bloch((0,0,0),F(1,4),MASSES[0]);mass_separation=float(np.linalg.norm(mutated[16:20,16:20]-right,2));assert mass_separation>.05
    rightclock=events([('Z',(0,),native.angle(0,0,F(-1,64)))])
    wrongclock=events([('Z',(0,),native.angle(0,0,F(-1,32)))])
    R,_,_,_=execute(rightclock,1,np.eye(2));W,_,_,_=execute(wrongclock,1,np.eye(2));clock_separation=float(np.linalg.norm(R-W,2));assert clock_separation>.015
    save('NATIVE_PROGRAMS.json',tapes)
    save('NATIVE_RECORDS.json',{'records':records,'fine_effects':fine_serial,'fine_probabilities_evaluated':fine_count,
        'fine_probability_diagnostic_sha256':finehash.hexdigest(),'all_fine_laws_formula':'Tr(E_r rho_probe) Tr(E_s rho_clock), on unnormalized success states; compose actual reset success weights and retain first timeouts, first source and pointer export outcomes.',
        'reset_contract':reset,'source_clock_product':'Exact disjoint factor on this free terminal apparatus; arbitrary unobserved field/reference traced before the controlled reduced-factor intervention.'})
    save('NATIVE_RESULT.json',{'all_checks_pass':True,'source_cases':sourcecases,'spin_preparations':len(preparations),'physical_flavors':3,
        'retained_flavor_padding':1,'native_controlled_columns_per_case':32,'tomography_component_records':len(records),'fine_probabilities_evaluated':fine_count,
        'maximum_native_Gram_residual':maxgram,'maximum_tomographic_operator_error':maxreconstruction,
        'pointer_all8column_Gram_residual':pgram,'scalar_control_error':scalar_error,'wrong_mass_source_separation':mass_separation,
        'wrong_clock_source_separation':clock_separation,'all_first_failures_retained':True,'success_normalized':False,
        'source_instrument':'New complete controlled reduced-factor apparatus; no per-pulse controlled-full-field failure tape is transferred.',
        'maximum_complete_native_TV_bound':max(float(F(row['complete_native_TV_bound_exact'])) for row in records)})
    assert verify()==pins
if __name__=='__main__':run()
