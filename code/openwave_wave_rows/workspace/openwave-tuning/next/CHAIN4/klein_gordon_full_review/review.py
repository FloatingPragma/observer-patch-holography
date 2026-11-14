"""Independent full ROW30 review on a copied, pinned source closure.

The candidate is never run with a write-capable main entry point. This review
rebuilds the Bloch source, slow projector, derivatives, response assembly and
pointer laws independently, while reusing the declared native response API.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib, json, math, sys
import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
ROOT=HERE/'inputs'
CANDIDATE=ROOT/'openwave-tuning/next/CHAIN3/klein_gordon_completion'
sys.dont_write_bytecode=True
sys.path.insert(0,str(CANDIDATE))
import native_spectroscopy as ns
a,b=ns.a,ns.b

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,value):(HERE/name).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def peak(x):return float(np.max(abs(x)))
def complex_pair(z):return [float(z.real),float(z.imag)]

def pins():
    expected=json.loads((HERE/'INPUTS.json').read_text())['files_sha256']
    for relative,h in expected.items():assert sha(ROOT/relative)==h,relative
    assert sha(CANDIDATE/'MANIFEST.json')=='de7ae79937061bdd1686c779d7660158416d73b5505c819dc3c2ef9d4371a05e'
    assert sha(CANDIDATE/'PROOF.md')=='e9a2a6d7688b86f6526567eefc4ef5fb4ee72512dc879b8871faae26fb9d872d'
    assert sha(CANDIDATE/'RESULT.json')=='5eac1b27bab13229a56b4436f75d2379870edb74f4dfdfbee1c9eb564c776b31'
    return expected

def source(k,mu=0):
    """Build from literal seam endpoints; no producer spectral matrices."""
    S=np.eye(144,dtype=complex)
    for direction,matching in enumerate(a.common.INTERNAL_MATCHINGS):
        for u,v in matching:
            x,y=12*u+2*direction+1,12*v+2*direction+1
            S[x,x]=S[y,y]=0;S[x,y]=S[y,x]=1
    for axis,pairs in enumerate(a.spatial.links()):
        for u,v in pairs:
            x,y=12*u+11,12*v+11
            S[x,x]=S[y,y]=0
            S[x,y]=np.exp(1j*k[axis]);S[y,x]=np.exp(-1j*k[axis])
    R=np.kron(np.eye(12),np.ones((12,12))/6-np.eye(12))
    return S@(np.cos(mu)*R-1j*np.sin(mu)*np.eye(144)),S,R

def rebuild_algebra():
    W,S,R=source((0,0,0));I=np.eye(144,dtype=np.int64)
    Rn=np.rint(6*R).astype(np.int64);S0=np.rint(S.real).astype(np.int64)
    Wn=S0@Rn;C=Wn+Wn.T
    Pn=I.copy();den=4608
    for shift in (-12,8,6,4):
        factor=C-shift*I
        assert 144*int(abs(Pn).max())*int(abs(factor).max())<2**63
        Pn=Pn@factor
    gcd=math.gcd(den,int(np.gcd.reduce(Pn.ravel())));Pn//=gcd;den//=gcd
    assert np.array_equal(Pn,Pn.T) and np.array_equal(Pn@Pn,den*Pn)
    assert np.array_equal(Wn@Pn,6*Pn) and np.trace(Pn)==26*den
    assert np.count_nonzero(np.linalg.svd(W-np.eye(144),compute_uv=False)<1e-12)==26
    assert np.array_equal(Rn@Pn,Pn@Rn)
    assert np.array_equal(Rn@Pn,den//12*np.ones((144,144),np.int64)-6*Pn)
    one=np.ones(144,np.int64);assert np.array_equal(Pn@one,den*one)
    H=[];numerators=[];derivatives=[]
    for pairs in a.spatial.links():
        D=np.zeros((144,144),np.int64)
        for u,v in pairs:D[12*u+11,12*v+11]=1;D[12*v+11,12*u+11]=-1
        hn=Pn@D@Rn@Pn;hd=6*den*den
        g=math.gcd(hd,int(np.gcd.reduce(hn.ravel())));hn//=g;hd//=g
        assert 1152%hd==0;hn*=1152//hd
        assert np.array_equal(hn,hn.T) and not np.any(Rn@hn+hn@Rn)
        gn=hn@one
        # Stronger than ranks/Gram alone: exact scalar-gradient form.
        assert np.array_equal(144*hn,gn[:,None]+gn[None,:])
        assert np.array_equal(Pn@hn,den*hn)
        H.append(hn/1152);numerators.append(gn);derivatives.append(D)
    for i in range(3):
        for j in range(3):
            assert sum(int(x)*int(y)for x,y in zip(numerators[i],numerators[j]))==(3*1152**2 if i==j else 0)
    e=one/12;P=Pn/den
    # Compare only after independently deriving every matrix.
    published=b.continuum.rational_bloch_data()
    assert peak(P-published['P'])==0
    assert all(peak(x-y)==0 for x,y in zip(H,published['H']))
    save('EXACT_ALGEBRA.json',{'rank':26,'projector_denominator':den,
        'projector_numerator':Pn.tolist(),'gradient_denominator':1152,
        'gradient_vector_numerators':[v.tolist()for v in numerators],
        'exact_source_derivative_reconstruction':True,'exact_scalar_gradient_identity':True,
        'full_nullspace_not_selected_eigenvectors':True,'stationary_slow_modes':24})
    return P,R,H,e

@lru_cache(None)
def row_order(targets,bits):
    """Index by local state first, untouched environment second."""
    rest=[j for j in range(bits)if j not in targets]
    result=[]
    for local in range(1<<len(targets)):
        base=sum(((local>>(len(targets)-1-i))&1)<<t for i,t in enumerate(targets))
        for env in range(1<<len(rest)):
            result.append(base+sum(((env>>i)&1)<<t for i,t in enumerate(rest)))
    return np.array(result)

def acted(state,matrix,targets,bits,restore=True):
    order=row_order(tuple(targets),bits);n=state.shape[1]
    arranged=state[order].reshape(1<<len(targets),-1)
    result=(matrix@arranged).reshape(-1,n)
    return result[np.argsort(order)] if restore else result

@lru_cache(None)
def responses(kind,angle,precision=96):
    d=4 if kind=='J'else 2
    event=a.certified_pulse(kind,tuple(range(2 if kind=='J'else 1)),angle,precision)
    a.verify_pulse(event)
    maps=a.execute_full_response(event)
    return maps[0][:d],tuple(K[d:]for K in maps[1:])

def native(word,bits,initial):
    state=np.asarray(initial,complex).copy()
    if state.ndim==1:state=state[:,None]
    failure=[];complements=[]
    for kind,targets,angle in word:
        success,others=responses(kind,tuple(angle))
        gram=np.zeros((state.shape[1],state.shape[1]),complex)
        for K in others:
            failed=acted(state,K,targets,bits,False);gram+=failed.conj().T@failed
        failure.append(gram);complements.append(len(others))
        state=acted(state,success,targets,bits)
    return state,failure,complements

def native_controls():
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1,-1]);I=np.eye(2)
    generators={'X':X,'Y':Y,'Z':Z,'I':I,'H':(X+Z)/np.sqrt(2),'J':np.diag([1,1,1,0])}
    rows=[]
    for angle in (F(1,7),F(-2,9)):
        for kind,generator in generators.items():
            width=2 if kind=='J'else 1;d=1<<width
            word=ns.controlled_event(kind,tuple(range(width)),a.angle(0,0,angle),width)
            success,failures,counts=native(word,width+1,np.eye(2*d))
            target=np.eye(2*d,dtype=complex);target[d:,d:]=expm(1j*float(angle)*generator)
            error=peak(success-target);gram=peak(success.conj().T@success+sum(failures)-np.eye(2*d))
            assert error<3e-12 and gram<3e-12 and all(counts)
            rows.append({'kind':kind,'angle_exact':str(angle),'columns':2*d,'pulses':len(word),'matrix_residual':error,'complete_gram_residual':gram,'fine_first_export_complements':sum(counts)})
    u=np.zeros(16);u[:12]=1/np.sqrt(12);reflection=2*np.outer(u,u)-np.eye(16)
    for mass in (F(-3,4),F(0),F(1)):
        h=F(1,4);mu=float(h*mass)
        word=list(ns.controlled_word(ns.native_word.coin_template(h,mass),4))
        success,failures,counts=native(word,5,np.eye(32))
        target=np.eye(32,dtype=complex)
        target[16:,16:]=np.cos(mu)*reflection-1j*np.sin(mu)*np.eye(16)
        error=peak(success-target);gram=peak(success.conj().T@success+sum(failures)-np.eye(32))
        assert error<3e-11 and gram<3e-11 and all(counts)
        rows.append({'kind':'full_massive_coin','mass_exact':str(mass),'columns':32,'pulses':len(word),'matrix_residual':error,'complete_gram_residual':gram,'fine_first_export_complements':sum(counts)})
    print('Independent native controlled matrices and all complements passed',flush=True)
    return rows

def independent_pointer():
    initial=np.eye(8)[:,[2,6]]
    success,failures,counts=native(b.detector_template(),3,initial)
    laws={}
    for law in ('C','E'):
        groups=defaultdict(lambda:np.zeros((2,2),complex));metadata={}
        for code in range(4):
            native_source=b.circuit.decode(code);K=success[[code,code+4],:]
            for coin in (0,1):
                out,(_,current)=b.native.native(native_source,b.pointer.SEAM,coin)
                label=current if law=='C'else native_source[0]-native_source[1]
                key=(coin,label,out);groups[key]+=K/np.sqrt(2)
                metadata[key]={'coin':coin,'current':current,'label':label,'bit':int(2*current+2*coin-1<0)}
        maps=[(metadata[key],K)for key,K in groups.items()]
        assert peak(sum(K.conj().T@K for _,K in maps)+sum(failures)-np.eye(2))<1e-13
        laws[law]=(maps,failures)
    return laws

def readout(rho,quadrature,law,pointers,clock_steps=0,size=2):
    word=[('Z',(0,),a.angle(0,0,-F(1,16*size)))]*clock_steps
    if quadrature=='Y':word += [('Z',(0,),a.angle(F(1,4)))]
    word+=b.H(0);success,failures,_=native(word,1,np.eye(2))
    alive=success@rho@success.conj().T;maps,pointer_failed=pointers[law]
    probabilities=[float(np.trace(K@alive@K.conj().T).real)for _,K in maps]
    failed=[float(np.trace(g@rho).real)for g in failures]+[float(np.trace(g@alive).real)for g in pointer_failed]
    assert abs(sum(probabilities)+sum(failed)-np.trace(rho).real)<2e-12
    assert min(probabilities+failed)>=-1e-13
    return sum((1-2*meta['bit'])*p for (meta,_),p in zip(maps,probabilities)),sum(probabilities)+sum(failed)

def finite_native(pointers):
    rows=[]
    # New sign, mixed-direction, zero and high integer-label cases.
    for size,mode in ((2,(0,0,0)),(2,(1,1,1)),(4,(3,1,2))):
        bits=b.member(size)[1];initial=np.zeros(1<<bits);initial[0]=1
        word=ns.preparation_template(size,mode)
        prepared,failures,counts=native(word,bits,initial)
        psi=prepared[:,0];failed=sum(float(g[0,0].real)for g in failures)
        expected=np.zeros((size**3,16,16),complex)
        for cell in range(size**3):
            coords=(cell//size**2,(cell//size)%size,cell%size)
            expected[cell,:12,:12]=np.exp(2j*np.pi*np.dot(mode,coords)/size)/(12*size**1.5)
        err=float(np.linalg.norm(psi-expected.ravel()));assert err<2e-11
        assert abs(np.vdot(psi,psi).real+failed-1)<2e-11
        for mass in (F(-3,4),F(0),F(1)):
            walk=ns.walk.Walk(size,h=F(1,size),mass=mass)
            state=psi.reshape(walk.geometry.shape)
            W,_,_=source(2*np.pi*np.asarray(mode)/size,float(mass)/size)
            e=np.ones(144)/12;fibre=e.astype(complex)
            for step in (1,2,3):
                state=walk.step(state);fibre=W@fibre
                c=np.vdot(psi,state.ravel());assert abs(c-np.vdot(e,fibre))<2e-11
                norm=np.vdot(psi,psi).real;rho=np.array([[norm,c.conjugate()],[c,norm]])/2
                for q in ('X','Y'):
                    for law in ('C','E'):
                        observed,total=readout(rho,q,law,pointers)
                        tick,clock_total=readout(np.ones((2,2))/2,q,law,pointers,step,size)
                        target=c.real if q=='X'else c.imag
                        expected_tick=math.cos(step/(8*size))if q=='X'else math.sin(step/(8*size))
                        assert abs(observed-target)<2e-11 and abs(tick-expected_tick)<2e-11
                        assert abs(total*clock_total+failed-1)<3e-11
                        rows.append({'size':size,'mode':list(mode),'mass_exact':str(mass),'step':step,'quadrature':q,'law':law,'contrast':observed,'clock_contrast':tick,'complete_probability':total*clock_total+failed,'preparation_error':err})
    print('Independent native Fourier/source/two-pointer cases passed',flush=True)
    return rows

def audit_frozen_records():
    settings=0;maxnormalization=0;maxcontrast=0
    for case in json.loads((CANDIDATE/'COMPLETE_RECORDS.json').read_text()):
        W,_,_=source(2*np.pi*np.asarray(case['mode'])/case['size'],float(F(case['mass_exact']))/case['size'])
        e=np.ones(144)/12;c=np.vdot(e,np.linalg.matrix_power(W,case['steps'])@e)
        for name,record in case['laws'].items():
            matter,clock=record['matter'],record['clock']
            assert all(r['probability']>=-1e-14 for r in matter+clock)
            assert {r['outcome']for r in matter}=={'analysis_export_failure','pointer_export_failure','native_pointer'}
            for pair in record['joint']:
                expected=matter[pair['matter']]['probability']*clock[pair['clock']]['probability']
                assert abs(pair['probability']-expected)<1e-15
            total=sum(x['probability']for x in record['joint'])+sum(record['preparation_export_failures'])
            maxnormalization=max(maxnormalization,abs(total-1));assert abs(total-1)<3e-11
            contrast=sum((1-2*r['record']['bit'])*r['probability']for r in matter if r['outcome']=='native_pointer')
            expected=c.real if name[0]=='X'else c.imag
            maxcontrast=max(maxcontrast,abs(contrast-expected));assert abs(contrast-expected)<3e-11
            settings+=1
    rest_settings=0
    for case in json.loads((CANDIDATE/'REST_CONTROL.json').read_text())['rows']:
        z=np.exp(-1.5j*float(F(case['mass_exact'])))
        assert abs(complex(*case['return'])-z)<2e-11
        for record in case['cases']:
            for side in ('matter','clock'):
                total=sum(r['probability']for r in record[side])
                if side=='matter':total+=sum(record['preparation_export_failures'])
                assert abs(total-1)<3e-11
            rest_settings+=1
    assert settings==288 and rest_settings==12
    return {'complete_joint_settings':settings,'rest_settings':rest_settings,'maximum_normalization_residual':maxnormalization,'maximum_independent_return_residual':maxcontrast,
            'record_scope':'Certified compressed source centre; reference-H/source failures priced by complete bound, not represented as literally sampled large-stream outcomes.'}

def rational_bounds(h,m,K,n):
    x=5*K/4;L=K+abs(m);q=K/4+abs(m)
    cm=x*(L+q)+F(13,16)*L*L+q*q/2
    vector=F(5,2)*h*K+n*h*h*cm
    frequency=h*h*(F(3,4)*(m*m+F(5,64)*K*K)**2+m**4/12+F(5,128)*m*m*K*K+K**4/16)
    return vector,frequency

def full_limit(P,R,H,e):
    rows=[];eigen=np.linalg.eigh(P);basis=eigen[1][:,eigen[0]>.5]
    assert basis.shape==(144,26)
    modes=((0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,-2,1))
    for mode in modes:
        p=2*np.pi*np.asarray(mode);K=F(13 if max(map(abs,mode))==2 else 7)
        for mass in (F(-3,4),F(0),F(1)):
            generator=1j*sum(p[i]*H[i]for i in range(3))-1j*float(mass)*R@P
            restricted=basis.conj().T@generator@basis
            limiting=basis@expm(restricted)
            energy=math.sqrt(float(mass)**2+float(p@p)/48)
            expected=1+0j if energy==0 else math.cos(energy)-1j*float(mass)/energy*math.sin(energy)
            assert abs(np.vdot(e,expm(generator)@e)-expected)<2e-12
            for size in (32,128,1024):
                h=F(1,size);W,S,_=source(float(h)*p,float(h*mass))
                actual=np.linalg.matrix_power(W,size)
                err=float(np.linalg.norm(actual@basis-limiting,2))
                vb,fb=rational_bounds(h,mass,K,size)
                assert err<=float(vb)+2e-10
                A=np.kron(np.eye(12),np.ones((12,1))/np.sqrt(12))
                values,vectors=np.linalg.eigh(A.conj().T@S@A)
                transfer=float(values[-1]);v=A@vectors[:,-1]
                polynomial=float(np.linalg.norm(W@W@v-2*transfer*math.cos(float(h*mass))*W@v+v))
                assert polynomial<2e-12
                omega=math.acos(min(1,max(-1,transfer*math.cos(float(h*mass)))))
                deviation=abs((omega/float(h))**2-energy**2)
                assert deviation<float(fb)+2e-8
                # Full rank26 source, not just the admitted scalar column.
                rows.append({'mode':list(mode),'mass_exact':str(mass),'size':size,'full_26_column_operator_error':err,'vector_bound_exact':str(vb),'squared_frequency_error':deviation,'squared_frequency_bound_exact':str(fb),'spectral_polynomial_residual':polynomial})
    # No assertion that the other 24 branches propagate with this dispersion.
    nonzero=1j*sum((2*np.pi)*H[i]for i in range(3))-1j*R@P
    energies=np.linalg.eigvalsh(1j*basis.conj().T@nonzero@basis)
    assert np.count_nonzero(abs(energies+1)<1e-12)==24
    print('Independent full 26-column massive continuum and dispersion checks passed',flush=True)
    return rows

def family_and_program():
    result=json.loads((CANDIDATE/'RESULT.json').read_text());rows=[]
    for supplied in result['formal_complete_native_family']:
        size=supplied['size'];m=F(supplied['mass_exact']);steps=supplied['steps'];h=F(1,size)
        actual=ns.certificate(size,steps,m,tuple(supplied['mode']),bandwidth=F(7))
        assert actual==supplied
        dp=F(actual['delta_preparation_exact']);dw=F(actual['delta_controlled_step_exact']);dr=F(actual['delta_clock_and_readout_exact'])
        reset=actual['reset'];Q=reset['blank_carriers']+reset['ready_carriers'];q=F(reset['base_failure_exact']);r=reset['blocks_per_register']
        assert Q*q**r==F(reset['timeout_bound_exact'])<=h*h and max(dp,dw,dr)<=h*h
        bound,_=rational_bounds(h,m,F(7),steps)
        assert bound==F(actual['finite_vs_continuum_vector_upper_exact'])
        assert bound+dp+steps*dw+dr+Q*q**r==F(actual['native_vs_continuum_uncapped_exact'])
        assert not actual['large_word_executed']
        rows.append({'size':size,'native_plus_continuum_exact':actual['native_vs_continuum_uncapped_exact'],'pulses_upper':actual['whole_experiment_pulses_upper'],'memory_bits':actual['total_memory_bits'],'all_failure_allowances_included':True})
    # A zero-macro word cheaply checks actual preparation/readout layouts.
    # Source routing is covered by the pinned full-member compiler theorem;
    # its enormous all-pulse stream is not expanded for this review.
    stream=list(ns.native_program(2,0,F(1),(1,0,1),'Y','Y',96))
    reset=stream[0];assert reset['operation']=='native_configuration_reset_batch'
    memory=a.layout(2)['memory_bits']
    assert reset['blank_wires']==list(range(memory+2))
    assert reset['ready_pairs']==[[memory+2,memory+3],[memory+4,memory+5]]
    pointers=[x for x in stream if x['operation']=='native_C_or_E'];assert len(pointers)==2
    assert {x['read_bit']for x in pointers}=={memory,memory+1}
    for event in stream:
        if event['operation']=='M6':assert a.verify_pulse(event)
    # Symbolic source envelope covers every actual reversible primitive,
    # without inheriting the earlier plus_i-only selector envelope bug.
    assert set(ns.controlled_primitive_templates())==set(a.primitive_templates())
    return rows,{'zero_macro_emitted_operations':len(stream),'native_pointers':2,'all_rational_controls_verified':True,'large_words_executed':False}

def controls(P,R,H,e,pointers):
    W,_,_=source(np.array([2*np.pi,0,0])/1024,1/1024)
    c=np.vdot(e,np.linalg.matrix_power(W,1024)@e)
    E=math.sqrt(1+4*np.pi*np.pi/48)
    wrong_time=math.cos(2*E)-1j*math.sin(2*E)/E
    wrong_mass=math.cos(2*np.pi/math.sqrt(48))
    changes={'double_time':abs(c-wrong_time),'mass_omitted':abs(c-wrong_mass),'mass_sign':abs(c-c.conjugate())}
    assert min(changes.values())>.1
    z=np.exp(.37j);rho=np.array([[1,z.conjugate()],[z,1]])/2
    phase=readout(rho,'Y','C',pointers)[0]-readout(np.ones((2,2))/2,'Y','C',pointers)[0]
    assert abs(phase)>.3;changes['controlled_global_phase_omitted']=abs(phase)
    event=dict(a.certified_pulse('Z',(0,),a.angle(F(1,4)),96));event['controls_exact']=list(event['controls_exact']);event['controls_exact'][0]='1'
    try:a.verify_pulse(event)
    except ArithmeticError:changes['mutated_native_controls_rejected']=True
    else:raise AssertionError('altered controls accepted')
    return changes

def main():
    before=pins();P,R,H,e=rebuild_algebra();print('Exact source-derived algebra passed',flush=True)
    native_rows=native_controls();pointer=independent_pointer()
    finite=finite_native(pointer);frozen=audit_frozen_records();limit=full_limit(P,R,H,e)
    family,program=family_and_program();negative=controls(P,R,H,e,pointer)
    assert pins()==before
    result={'schema':'oph.row30.independent_full_review.v1','row':30,
      'criterion':'emergent ω² = k² + m² for matter waves',
      'candidate_manifest_sha256':sha(CANDIDATE/'MANIFEST.json'),
      'decision':'PASS_FULL_LITERAL_CRITERION_ON_DECLARED_NATIVE_ARCHITECTURE',
      'reviewer_independent_of_row30_producer':True,'review_type':'full_criterion',
      'native_controlled_matrices':native_rows,'native_new_joint_measurements':finite,
      'all_frozen_fine_records':frozen,'full_26_mode_limit':limit,
      'full_native_family':family,'program_layout':program,'negative_controls':negative,
      'architecture_inputs':['coherent native M6 proper-code load/export','LOCAL12 cubic routing','coherent event reuse','rational mass parameter','product reset apparatus','common source macro clock'],
      'empirical_mass_prediction':False,'interacting_field_limit_claimed':False,
      'large_native_source_stream_literally_executed':False,
      'executed_scope':'Native local responses/preparation/pointers; complete finite source arrays and full Bloch evolution; unexpanded source compiler with uniform complete-instrument error theorem.',
      'all_candidate_and_provider_bytes_unchanged':True,'frozen_input_count':len(before),
      'official_rows_changed':[],'second_independent_full_review_still_required':True}
    save('REVIEW.json',result)
    imports={}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path:
            p=Path(path).resolve()
            if p.suffix=='.pyc'and p.with_suffix('.py').exists():p=p.with_suffix('.py')
            if p.is_file()and p.is_relative_to(ROOT):
                relative=str(p.relative_to(ROOT));assert relative in before,relative;imports[relative]=sha(p)
    save('IMPORTS.json',{'files_sha256':imports})
    print(json.dumps({'decision':result['decision'],'native_columns':sum(x['columns']for x in native_rows),'new_joint_settings':len(finite),'full_26_mode_cases':len(limit),'input_files':len(before)}),flush=True)

if __name__=='__main__':main()
