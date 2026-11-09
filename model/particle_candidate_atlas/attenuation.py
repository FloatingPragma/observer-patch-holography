"""Complete coherent attenuation of the existing native Koide analyzer."""
from fractions import Fraction as F
from functools import lru_cache
import json
import math
import sympy as s

from model.standard_model_source_map import koide_source as source
from model.observer.source_scalar import compiler


def clean(matrix):
    return s.ImmutableMatrix(s.Matrix(matrix).applyfunc(lambda x:s.factor(s.radsimp(s.expand(x)))))


def rational(value):
    if type(value) not in (int,F,str):
        raise ValueError('an exact rational input is required')
    return F(value.numerator,value.denominator) if type(value)is F else F(value)


def attenuation_parameter(value):
    r=rational(value)
    if not 0<=r<=1:raise ValueError('pass amplitude must lie in [0,1]')
    return s.Rational(r.numerator,r.denominator)


def chi_from_declared_u(value):
    """Evaluate the stipulated MCPR rate; this does not select its architecture."""
    u=rational(value)
    if not 0<=u<77:raise ValueError('declared u must lie in [0,77)')
    return u*u/(512*(1-u/77))


@lru_cache(maxsize=16)
def _factors(r,z):
    prep,fourier,blocks=source.analyzer_factors(z)
    loss=s.eye(16);ell=s.sqrt(1-r*r)
    for flag in (1,2):
        loss[flag,flag]=r;loss[flag,8+flag]=-ell
        loss[8+flag,flag]=ell;loss[8+flag,8+flag]=r
    loss=clean(loss)
    if clean(loss.H*loss)!=s.eye(16):raise ArithmeticError('complete loss rotation failed')
    return prep,fourier,tuple(blocks),loss


@lru_cache(maxsize=128)
def _block(out_flag,in_flag,r,phase):
    prep,fourier,blocks,loss=_factors(r,phase)
    lo,fo=divmod(out_flag,8);li,fi=divmod(in_flag,8)
    return clean(sum((fourier[fo,j]*loss[8*lo+j,8*li+j]*prep[j,fi]*blocks[j]
                     for j in range(8)),s.zeros(4)))


def analyzer_block(out_flag,in_flag,r='3/5',phase=1):
    """4-family-coordinate block of the complete 64-column analyzer."""
    if any(type(j)is not int or not 0<=j<16 for j in (out_flag,in_flag)):
        raise ValueError('all sixteen flag indices are explicit')
    return _block(out_flag,in_flag,attenuation_parameter(r),source._phase(phase))


def kraus(r='3/5',phase=1):
    r=attenuation_parameter(r);phase=source._phase(phase)
    a,b=s.sqrt(source.WEIGHTS[0]),s.sqrt(source.WEIGHTS[1])
    rot=source.P;omega=-s.Rational(1,2)+s.I*s.sqrt(3)/2;ell=s.sqrt(1-r*r)
    maps=[]
    for loss in (0,1):
        for flag in range(8):
            if flag<3:
                oriented=b*(omega**(-flag)*phase*rot+omega**(-2*flag)*s.conjugate(phase)*rot**2)
                k=(a*s.eye(3)+r*oriented)/s.sqrt(3) if loss==0 else ell*oriented/s.sqrt(3)
            elif loss==0 and flag in (3,4):k=s.sqrt(source.WEIGHTS[flag])*s.eye(3)
            else:k=s.zeros(3)
            maps.append(clean(k))
    return tuple(maps)


def family(r='3/5',phase=1):
    r=attenuation_parameter(r);phase=source._phase(phase)
    a,b=s.sqrt(source.WEIGHTS[0]),s.sqrt(source.WEIGHTS[1]);rot=source.P
    c=clean(a*s.eye(3)+r*b*(phase*rot+s.conjugate(phase)*rot**2))
    y=clean(c*c/3)
    ks=kraus(str(r),phase)
    if clean(sum((k.H*k for k in ks),s.zeros(3)))!=s.eye(3):
        raise ArithmeticError('attenuation lost a complementary record')
    if clean(ks[0].H*ks[0]-y)!=s.zeros(3):raise ArithmeticError('action compression mismatch')
    positive=all(s.simplify(c[:n,:n].det()).is_positive is True for n in (1,2,3))
    return {'C':c,'Y':y,'r':r,'phase':phase,'positive_roots':positive,
            'signed_Q':s.factor(s.trace(c*c)/s.trace(c)**2),
            'original_tracial_records':tuple(s.simplify(s.trace(k.H*k)/3) for k in ks)}


def attenuator_word(r='3/5'):
    """Native four-qubit loss rotation, including its unchanged complementary columns."""
    r=attenuation_parameter(r);theta=math.acos(float(r));word=[]
    # Pi_{flag=1 or2}=(I+Z2-Z0Z1-Z0Z1Z2)/4.
    for bits,sign in (((),1),((2,),1),((0,1),-1),((0,1,2),-1)):
        label=['I']*4;label[3]='Y'
        for bit in bits:label[bit]='Z'
        word+=compiler.pauli_rotation(label,-sign*theta/4)
    return word


def weak_kraus(r='3/5',phase=1,tangent='1/64'):
    t=rational(tangent)
    if t<=0:raise ValueError('positive rational half-angle tangent required')
    t=s.Rational(t.numerator,t.denominator);co=(1-t*t)/(1+t*t);si=2*t/(1+t*t)
    x=s.Matrix([[0,1],[1,0]]);e=(co-1)*s.eye(2)-s.I*si*x
    first=analyzer_block(0,0,r,phase)[:3,:3]
    return tuple(clean((s.eye(6)if flag==0 else s.zeros(6))+
                 s.kronecker_product(e,analyzer_block(0,flag,r,phase)[:3,:3].H*first))
                 for flag in range(16))


def budget(steps,time,*,higgs_norm=1,round_isometry_error=0):
    if type(steps)is not int or steps<=0:raise ValueError('positive integer rounds required')
    t,v,e=map(rational,(time,higgs_norm,round_isometry_error))
    if min(t,v,e)<0:raise ValueError('nonnegative bounds required')
    return {'logical_complete_unhalved_error':str(min(F(2),9*t*t/(4*steps)+2*steps*e)),
            'logical_nonzero_archive_probability':str(min(F(1),t*t/(4*steps)+steps*e)),
            'full_CAR_pairing_norm':str(6*v),
            'full_CAR_complete_unhalved_error':str(min(F(2),108*v*v*t*t/steps+2*steps*e)),
            'full_CAR_nonvacuum_probability':str(min(F(1),36*v*v*t*t/steps+steps*e)),
            'zero_outcome_postselection':False,'fixed_native_service_clock_equals_weak_clock':False}


def certificate():
    native=source.native_partition()
    if tuple(source.WEIGHTS)!=(s.Rational(10,72),s.Rational(5,72),s.Rational(5,72),s.Rational(2,72),s.Rational(50,72)):
        raise ArithmeticError('native source partition changed')
    rows=[]
    for r in ('1','3/5','0'):
        f=family(r)
        rows.append({'pass_amplitude':r,'Q':str(f['signed_Q']),
                     'positive_root_chamber':f['positive_roots'],
                     'original_probabilities':list(map(str,f['original_tracial_records'])),
                     'trace_Y':str(s.trace(f['Y']))})
    return {'native_partition_replayed':native,'rows':rows,'retained_flags':16,
            'old_flag_count':8,'added_neutral_loss_flag_bits':1,
            'same_original_C_E_label_law':True,'complete_instrument_proved':True,
            'full_analyzer_columns':64,'native_full_analyzer_execution_claimed':True,
            'native_execution':native_execution('3/5',s.Rational(99,101)+s.I*s.Rational(20,101)),
            'native_weak_round':native_round_execution('3/5',1,'1/64'),
            'native_loss_word_pulses':len(attenuator_word()),
            'full_CAR_extension_uses_existing_all_occupation_theorem':True,
            'full_CAR_physical_modes':90,'full_CAR_with_fresh_flags_modes':1440,
            'gauge_neutral_flag_coordinate_not_gauge_neutral_charged_modes':True,
            'MCPR_r_and_phase_are_supplied_to_controller':True,
            'chi_selected_from_primitive_source':False,'physical_matching_proved':False,
            'shared_clock_scope':'t=N h is accumulated supplied weak-kick angle; service overhead retained',
            'budget_example':budget(1024,'1/16')}


def _fourier_word():
    import numpy as np
    c=math.sqrt((1+1/math.sqrt(3))/2);d=math.sqrt((1-1/math.sqrt(3))/2)
    basis=np.eye(8)
    basis[:3,:3]=[[c,-d,0],[d/math.sqrt(2),c/math.sqrt(2),1/math.sqrt(2)],
                 [d/math.sqrt(2),c/math.sqrt(2),-1/math.sqrt(2)]]
    tape=compiler.compile_real_orthogonal(basis)
    return compiler.remap(compiler.inverse(tape)+compiler.compile_diagonal(
        [0,math.pi,-math.pi/2,0,0,0,0,0])+tape,[2,3,4])


def analyzer_word(r='3/5',phase=1):
    """Insert loss before the existing Fourier suffix; six logical qubits."""
    old=source.analyzer_word(phase);suffix=_fourier_word()
    if old[-len(suffix):]!=suffix:raise ArithmeticError('native analyzer Fourier schedule changed')
    return old[:-len(suffix)]+compiler.remap(attenuator_word(r),[2,3,4,5])+suffix


@lru_cache(maxsize=4)
def _native_execution_json(r,phase):
    """All-column native replay; immutable canonical cached result."""
    import numpy as np
    rr=attenuation_parameter(r);z=source._phase(phase)
    prep,f,blocks,loss=_factors(rr,z)
    prep,f,loss=map(lambda m:np.asarray(m,complex),(prep,f,loss))
    control=np.zeros((64,64),complex)
    for j,block in enumerate(blocks+blocks):control[4*j:4*j+4,4*j:4*j+4]=np.asarray(block,complex)
    exact=np.kron(np.kron(np.eye(2),f),np.eye(4))@np.kron(loss,np.eye(4))@control@np.kron(np.kron(np.eye(2),prep),np.eye(4))
    word=analyzer_word(r,phase);actual=compiler.execute(word,6)
    tolerance=2e-9
    if np.linalg.norm(actual-exact,2)>tolerance:raise ArithmeticError('complete native analyzer differs')
    if np.linalg.norm(actual.conj().T@actual-np.eye(64),2)>tolerance:raise ArithmeticError('native full-column Gram failed')
    probabilities=[];expected_kraus=kraus(str(rr),z)
    for flag in range(16):
        value=float(np.sum(np.abs(actual[4*flag:4*flag+3,:3])**2)/3)
        expected=s.simplify(s.trace(expected_kraus[flag].H*expected_kraus[flag])/3)
        if abs(value-float(expected))>tolerance:raise ArithmeticError('native original flag probability changed')
        probabilities.append({'flag':flag,'exact_probability':str(expected)})
    return json.dumps({'pass_amplitude':str(rr),'phase':str(z),'full_input_columns':64,
            'pulses':len(word),'loss_pulses':len(attenuator_word(r)),
            'indexed_failure_slots':sum(6-(4 if len(event['targets'])==2 else 2)for event in word),
            'all_failure_slots_retained':True,'global_phase_fitted':False,
            'numerical_admission_upper':'1/500000000','certified_round_dilation_error':None,
            'native_full_analyzer_execution_passed':True,
            'original_tracial_records':probabilities,
            'physical_full_CAR_execution_claimed':False},sort_keys=True)


def native_execution(r='3/5',phase=1):
    return json.loads(_native_execution_json(str(attenuation_parameter(r)),source._phase(phase)))


def weak_word(r='3/5',phase=1,tangent='1/64'):
    """Chronological complete native analyzer, controlled kick and inverse."""
    t=rational(tangent)
    if t<=0:raise ValueError('positive rational half-angle tangent required')
    h=2*math.atan(float(t));analyzer=analyzer_word(r,phase);kick=[]
    for mask in range(16):
        label=['I']*7;label[6]='X'
        for bit in range(4):
            if mask>>bit&1:label[2+bit]='Z'
        kick+=compiler.pauli_rotation(label,-h/16)
    return analyzer+kick+compiler.inverse(analyzer)


@lru_cache(maxsize=2)
def _native_round_json(r,phase,tangent):
    import numpy as np
    rr=attenuation_parameter(r);z=source._phase(phase)
    prep,f,blocks,loss=_factors(rr,z)
    control=np.zeros((64,64),complex)
    for j,block in enumerate(blocks+blocks):control[4*j:4*j+4,4*j:4*j+4]=np.asarray(block,complex)
    analyzer=np.kron(np.kron(np.eye(2),np.asarray(f,complex)),np.eye(4))@np.kron(np.asarray(loss,complex),np.eye(4))@control@np.kron(np.kron(np.eye(2),np.asarray(prep,complex)),np.eye(4))
    projection=analyzer[:4].conj().T@analyzer[:4]
    t=rational(tangent);co=float((1-t*t)/(1+t*t));si=float(2*t/(1+t*t))
    ideal=np.eye(128)+(co-1)*np.kron(np.eye(2),projection)-1j*si*np.kron(np.array([[0,1],[1,0]]),projection)
    tape=weak_word(r,phase,tangent);actual=compiler.execute(tape,7);tol=4e-9
    if np.linalg.norm(actual-ideal,2)>tol:raise ArithmeticError('complete native weak projection changed')
    if np.linalg.norm(actual.conj().T@actual-np.eye(128),2)>tol:raise ArithmeticError('native weak Gram changed')
    inputs=[0,1,2,64,65,66];expected=weak_kraus(r,phase,tangent);records=[]
    for flag in range(16):
        outputs=[64*n+4*flag+f for n in (0,1) for f in range(3)]
        block=actual[np.ix_(outputs,inputs)]
        if np.linalg.norm(block-np.asarray(expected[flag],complex),2)>tol:raise ArithmeticError('a native reset record changed')
        probability=s.simplify(s.trace(expected[flag].H*expected[flag])/6)
        records.append({'flag':flag,'exact_tracial_probability':str(probability)})
    return json.dumps({'full_input_columns':128,'retained_return_flags':16,
        'half_angle_tangent':tangent,'pulses':len(tape),
        'indexed_failure_slots':sum(6-(4 if len(e['targets'])==2 else 2)for e in tape),
        'original_return_records':records,'all_native_failures_retained':True,
        'global_phase_fitted':False,'numerical_admission_upper':'1/250000000',
        'certified_round_dilation_error':None,'physical_full_CAR_execution_claimed':False},sort_keys=True)


def native_round_execution(r='3/5',phase=1,tangent='1/64'):
    return json.loads(_native_round_json(str(attenuation_parameter(r)),source._phase(phase),str(rational(tangent))))
