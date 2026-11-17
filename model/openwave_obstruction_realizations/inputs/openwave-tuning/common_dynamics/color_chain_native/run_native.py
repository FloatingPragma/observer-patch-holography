"""Actual six-level M6 all-column blocks for the finite colour-chain compiler.

No giant accepted filter word is executed. Its exact controls are defined by
the lazy iterator and its uniform certificate; these are small block tests.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json,math
import numpy as np
from scipy.linalg import expm
import native_path as p
import source

HERE=Path(__file__).resolve().parent
PRECISION=220

def matrix(z):return dict(real=z.real.tolist(),imag=z.imag.tolist())

def embed(A,targets,bits):
    dim=1<<bits;out=np.zeros((dim,dim),complex)
    for col in range(dim):
        sub=sum(((col>>q)&1)<<(len(targets)-1-j) for j,q in enumerate(targets))
        base=col&~sum(1<<q for q in targets)
        for row in range(len(A)):
            dest=base+sum(((row>>(len(targets)-1-j))&1)<<q for j,q in enumerate(targets))
            out[dest,col]=A[row,sub]
    return out

def run(name,word,bits,target):
    generators=p.car.field_replay().source_generators()
    U=np.eye(1<<bits,dtype=complex);G=np.zeros_like(U);first=[];codes=[];unique=[];index={};cache={};bound=F()
    for i,(kind,targets,angle) in enumerate(word):
        event=p.native.certified_pulse(kind,targets,angle,PRECISION);p.native.verify_pulse(event)
        key=(kind,tuple(event['controls_exact']))
        if key not in cache:
            W=expm(sum(float(F(c))*a for c,a in zip(event['controls_exact'],generators)))
            cache[key]=W;index[key]=len(unique)
            unique.append({k:v for k,v in event.items() if k!='targets'})
        W=cache[key];codes.append(dict(control=index[key],targets=list(targets)))
        d=1<<len(targets);mask=sum(1<<q for q in targets)
        bases=np.arange(1<<bits)[(np.arange(1<<bits)&mask)==0]
        indices=np.array([bases+sum(((z>>(len(targets)-1-j))&1)<<q for j,q in enumerate(targets)) for z in range(d)])
        arranged=U[indices];complements=[]
        for q in range(d,6):
            rejected=np.einsum('b,bcr->cr',W[q,:d],arranged)
            loss=rejected.conj().T@rejected;G+=loss
            complements.append(np.diag(loss).real.tolist())
        first.append(dict(pulse=i,complement_probability_by_input=complements))
        U[indices]=np.einsum('ab,bcr->acr',W[:d,:d],arranged)
        bound+=F(event['full_response_operator_error_exact'])
    error=float(np.linalg.norm(U-target,2));gram=float(np.linalg.norm(U.conj().T@U+G-np.eye(len(U)),2))
    assert error<2e-10,(name,error)
    assert gram<2e-10,(name,gram)
    result=dict(name=name,bits=bits,columns=len(U),pulses=len(codes),unique_controls=unique,program=codes,
                final_success_matrix=matrix(U),summed_first_failure_effect=matrix(G),every_first_export=first,
                exact_full_response_operator_error_upper=str(bound),diagnostic_operator_error=error,
                diagnostic_complete_gram_error=gram,interpretation='floating raw-six-level execution; exact controls/bounds independently supplied by adapter')
    print(name,'columns',len(U),'pulses',len(codes),'error',error,'Gram',gram,flush=True)
    return result,U,G

def main():
    results=[]
    H=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
    C=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],complex)
    r,U,G=run('Bell purification',list(p.mixing_word((0,),(1,))),2,embed(C,(0,1),2)@embed(H,(0,),2));results.append(r)
    m=3;dim=1<<m;Fmat=np.exp(2j*np.pi*np.outer(np.arange(dim),np.arange(dim))/dim)/np.sqrt(dim)
    r,U,G=run('inverse QFT3',p.qft_word((0,1,2),True),3,Fmat.conj().T);results.append(r)
    # Actual physical source: every L2 colour-chain matrix entry is included.
    d=source.data(2);rounded=source.rounding(d,PRECISION);dt=F(1,1024)
    prod=np.diag(np.exp(-1j*float(dt)*np.array(rounded['diagonal'],float)))
    for i,j,c in rounded['pairs']:
        a=float(dt*c);R=np.eye(16,dtype=complex);R[i,i]=R[j,j]=np.cos(a);R[i,j]=R[j,i]=-1j*np.sin(a);prod=R@prod
    target=np.eye(32,dtype=complex);target[16:,16:]=prod
    word=p.evolution_word(d,(F(),F(),dt),1,PRECISION,control=4)
    r,U,G=run('complete L2 controlled Lie slice',word,5,target)
    exact=np.eye(32,dtype=complex);exact[16:,16:]=expm(-1j*float(dt)*source.dense(d))
    r['one_slice_vs_exact_source_diagnostic']=float(np.linalg.norm(U-exact,2))
    r['one_slice_trotter_upper_exact']=str(dt*dt*rounded['normsum_upper']**2/2)
    r['missing_control_mutation_distance']=float(np.linalg.norm(target-np.kron(np.eye(2),prod),2))
    assert r['missing_control_mutation_distance']>.01
    results.append(r)
    # QPE complete test on a finite 2-level Hamiltonian; it tests orientation,
    # scalar phases and the actual inverse Fourier circuit, not an extra row.
    toy=dict(bits=1,code_dimension=2,dimension=2,diagonal=[F(7,2),F(4)],pairs=[],dipole=[F(),F()])
    phase=(1,2,3);M=8;controlled=np.zeros((16,16),complex)
    for x in range(16):controlled[x,x]=np.exp(2j*np.pi*(x>>1)*float(toy['diagonal'][x&1])/32)
    prep=np.eye(16,dtype=complex)
    for bit in phase:prep=embed(H,(bit,),4)@prep
    inverse=np.zeros((16,16),complex)
    for a in range(16):
        for b in range(16):
            if (a&1)==(b&1):inverse[a,b]=Fmat.conj().T[a>>1,b>>1]
    r,U,G=run('complete small QPE orientation',p.qpe_unitary(toy,phase,[1,1,1],PRECISION),4,inverse@controlled@prep)
    r['reversed_source_sign_mutation_distance']=float(np.linalg.norm(inverse@controlled@prep-inverse@controlled.conj()@prep,2));assert r['reversed_source_sign_mutation_distance']>.5
    results.append(r)
    # Complete pointer unitary on all8 columns; then literal C/E fine channels.
    pointer_target=embed(C,(0,1),3)@embed(C,(2,0),3)@embed(np.array([[0,1],[1,0]],complex),(1,),3)
    r,W,G=run('native C/E pointer',p.pointer_word(2,(0,1)),3,pointer_target)
    from model.observer import native as event
    from model.native_interface import pointer
    from model.observer.source_discrimination import circuit
    selected=W[:,[0,4]];laws={}
    for law in ['C','E']:
        groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)));metadata={}
        for code in range(4):
            state=circuit.decode(code);K=selected[[code,code+4]]
            for coin in [0,1]:
                out,(_,current)=event.native(state,pointer.SEAM,coin);label=current if law=='C' else state[0]-state[1]
                key=(coin,label);groups[key][out]+=K/np.sqrt(2)
                metadata[key]=dict(coin=coin,label=label,current=current,bit=int(2*current+2*coin-1<0))
        rows=[];channels={0:np.zeros((4,4),complex),1:np.zeros((4,4),complex)}
        for key,ks in sorted(groups.items()):
            effect=sum((K.conj().T@K for K in ks.values()),np.zeros((2,2),complex))
            channel=sum((np.kron(K.conj(),K) for K in ks.values()),np.zeros((4,4),complex))
            channels[metadata[key]['bit']]+=channel
            rows.append(dict(record=metadata[key],kraus=[matrix(K) for K in ks.values()],effect=matrix(effect)))
        err=0
        for bit in [0,1]:
            P=np.diag([int(bit==0),int(bit==1)]);err=max(err,float(np.max(np.abs(channels[bit]-np.kron(P,P)))))
        assert err<2e-12
        laws[law]=dict(all_fine_outcomes=rows,decoded_projective_channel_error=err)
    r['complete_laws']=laws;results.append(r)
    for axis in ['X','Y']:
        time=F(1,8);word=list(p.had(0))+[('Z',(0,),p.native.angle(0,0,-time/16))]
        target=expm(-1j*float(time/16)*np.diag([1,-1]))@H
        if axis=='Y':word.append(('Z',(0,),p.native.angle(F(1,4))));target=expm(1j*np.pi/4*np.diag([1,-1]))@target
        word+=list(p.had(0));target=H@target
        r,U,G=run('source-time clock '+axis,word,1,target);expected=(math.cos if axis=='X' else math.sin)(float(time/8))
        r['source_time_exact']=str(time);r['expected_sign_contrast']=expected;r['native_sign_contrast']=float(abs(U[0,0])**2-abs(U[1,0])**2)
        assert abs(r['native_sign_contrast']-expected)<1e-12;results.append(r)
    out=dict(precision_bits=PRECISION,executed_blocks=results,total_executed_pulses=sum(r['pulses'] for r in results),
             full_filter_executed=False,all_first_export_complements_retained=True,
             large_word_scope='Only the finite exact compiler/error theorem and lazy iterator specify the accepted high-precision filter/hold experiment.')
    (HERE/'NATIVE_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
