"""Independent fixed-frame exact arithmetic and selected native size3 replay."""
from review import *

def up(x):
    q=F(x);r=math.isqrt(q.numerator*q.denominator)
    return F(r+(r*r<q.numerator*q.denominator),q.denominator)
def down(x):
    q=F(x);return F(math.isqrt(q.numerator*q.denominator),q.denominator)

def exact_dot(L,R):
    assert L.shape[1]*int(abs(L).max())*int(abs(R).max())<2**53
    a=L.astype(float)@R.astype(float)
    assert np.array_equal(a,np.rint(a))
    return a.astype(np.int64).astype(object)

def fixed_frame():
    d=load(SIZE/'FRAME.json');geo=load(SIZE/'GEOMETRY.json');cert=load(SIZE/'CERTIFICATE.json')['certificate']
    assert sha(SIZE/'FRAME.json')==load(SIZE/'CERTIFICATE.json')['frame_sha256']
    Q=np.array(d['frame_integers'],np.int64);lam=np.array(d['curvature_numerators'],np.int64);den=d['denominator']
    B=np.zeros((324,972),np.int64);A=np.zeros((972,702),np.int64);period=np.zeros((3,972),np.int64)
    for j,(u,v) in enumerate(geo['edges']):B[u,j]=1;B[v,j]=-1
    for i,j,v in geo['face_entries']:A[i,j]=v
    for i,j,v in geo['period_times3_entries']:period[i,j]=v
    winding=np.array(geo['winding_columns'],np.int64).T
    assert not np.any(B@A) and not np.any(period@A) and not np.any(B@winding)
    assert np.array_equal(period@winding,3*np.eye(3,dtype=int))
    K=A@A.T
    # Four centered 15-bit limbs, independent of producer's three unsigned
    # 18-bit decomposition and summation order. Integer GEMM remains exact.
    base=1<<15;tail=Q.copy();limbs=[]
    for _ in range(3):
        high=(tail+base//2)//base;limbs.append(tail-high*base);tail=high
    limbs.append(tail)
    assert np.array_equal(sum((v.astype(object)*base**i for i,v in enumerate(limbs)),np.zeros(Q.shape,dtype=object)),Q)
    gram=np.zeros((646,646),dtype=object)
    for i,left in enumerate(limbs):
        for j in range(i,len(limbs)):
            piece=exact_dot(left.T,limbs[j]);gram+=(piece if i==j else piece+piece.T)*base**(i+j)
        print('frame limb',i,'done',flush=True)
    gram-=den**2*np.eye(646,dtype=object)
    KQ=sum((exact_dot(K,x)*base**i for i,x in enumerate(limbs)),np.zeros(Q.shape,dtype=object))
    resid=den*KQ-Q.astype(object)*lam.astype(object)[None,:]
    gram2=F(sum(int(x)**2 for x in gram.flat),den**4)
    resid2=F(sum(int(x)**2 for x in resid.flat),den**4)
    assert F(cert['Gram_squared_Frobenius_exact'])==gram2
    assert F(cert['residual_squared_Frobenius_exact'])==resid2
    eps=up(gram2);r=up(resid2);minimum=F(int(lam[0]),den);maximum=F(int(lam[-1]),den)
    delta=r/minimum;gamma=eps+2*(1+eps)*delta+delta**2
    eta=r+(F(1504,17)+maximum)*(delta+gamma)
    freq=F(5,4)*eta/down(minimum);cap=324*freq
    assert F(cert['cap1_H_error_including_vacuum'])==cap and gamma<F(1,2) and eta<minimum/2
    return stage('FRAME',{'fixed_frame_sha256':sha(SIZE/'FRAME.json'),'modes':646,'physical_states':647,'padding':377,
        'independent_centered_integer_limb_bits':15,'integer_limbs':4,'Gram_squared_Frobenius_exact':str(gram2),
        'residual_squared_Frobenius_exact':str(resid2),'cap1_H_error_including_vacuum':str(cap),
        'cap1_H_error_float':float(cap),'full_integer_BA_and_period_winding':True,
        'runtime_eigensolver_used':False,'fixed_rational_eigenbasis_retained':True})

def selected_native():
    sys.path.insert(0,str(SIZE));import native_member,walsh
    program=load(SIZE/'PROGRAM.json');data=load(SIZE/'NATIVE.json')
    word=program['diagonal_word'];phases,frequencies,E0=native_member.targets()
    # Consume every scalar and variable coefficient against all phase targets.
    error=walsh.certify(word,phases)
    assert error==F(data['all_phase_operator_error_exact'])
    midpoint=np.array([float(F(x.lo+x.hi,2*walsh.nf.iv.SCALE)) for x in phases])
    chosen=[0,1,646,647,1023]
    M,fail,gram=execute(word,10,np.eye(1024,dtype=complex)[:,chosen])
    target=np.eye(1024,dtype=complex)[:,chosen]*np.exp(1j*midpoint[chosen])
    residual=float(np.linalg.norm(M-target,2));assert residual<float(error)+2e-10 and gram<3e-10
    # Replay full prepared native field/clock state once. All four settings
    # share this emitted diagonal and differ only by audited local analysis.
    fieldword=program['experiments']['XY'];state=np.zeros(2048,complex);state[0]=1
    out,failed,gg=execute(fieldword,11,state);out=out[:,0]
    stored=next(x for x in data['experiments'] if x['field_quadrature']=='X' and x['clock_quadrature']=='Y')
    branches,pg,pl=native_pointer();maxprob=0.
    for law,values in branches.items():
        E=[sum((E for meta,E in values if meta['bit']==b),np.zeros((2,2),complex)) for b in (0,1)]
        assert max(abs(e[0,1]) for e in E)<1e-13
        probs=abs(out)**2;transition=np.array([[E[b][x,x].real for x in (0,1)] for b in (0,1)])
        for bit in range(11):
            indices=slots(11,(bit,));probs[indices]=transition@probs[indices]
        frozen=stored['laws'][law]['all2048_joint_probabilities']
        delta=float(np.max(abs(probs-np.array(frozen))));assert delta<3e-10;maxprob=max(maxprob,delta)
        assert abs(sum(probs)+failed[0,0].real-1)<3e-10
    visible_scalar=abs(1-np.exp(1j*float(F(E0.lo+E0.hi,2*walsh.nf.iv.SCALE))))
    assert visible_scalar>.1
    return stage('SIZE3_NATIVE',{'all_phase_targets':1024,'diagonal_pulses':len(word),'selected_native_columns':chosen,
        'selected_phase_residual':residual,'selected_complete_Gram_residual':gram,
        'full_prepared_field_clock_word_pulses':len(fieldword),'full_record_outcomes_per_law':2048,
        'prepared_complete_Gram_residual':gg,'max_record_difference':maxprob,
        'visible_scalar_deletion_defect':visible_scalar,'word_certificate_exact':str(error),
        'all_first_export_complements_accumulated':True,'success_normalized':False})

if __name__=='__main__':
    pins=pincheck();fixed_frame();selected_native();assert pincheck()==pins
    stage('SIZE3_CHECKS',{'all_checks_pass':True,'producer_pins_unchanged':True})
