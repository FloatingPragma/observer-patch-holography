"""Stream a complete native source-clock experiment on the full CAR dimer."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import isqrt
import math

import numpy as np
import sympy as s

from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert, interval as iv
from model.native_spectral_detector import circuits
from model.native_spin_observer import execute
from model.protected_particle_readout import execute as local
from .local_spectral import dimer_source

COUNT=2048
STEP=F(1,256)
RATIO=F(256,257)
BITS=4


def root_interval():
    grid=2**90
    lower=F(isqrt(17*grid*grid),grid)
    return lower,lower+F(1,grid)


def tensor(label):
    matrices={'I':s.eye(2),'X':s.Matrix([[0,1],[1,0]]),
              'Y':s.Matrix([[0,-s.I],[s.I,0]]),'Z':s.diag(1,-1)}
    result=s.ones(1,1)
    for letter in reversed(label):result=s.kronecker_product(result,matrices[letter])
    return result


@lru_cache(maxsize=1)
def _source():
    h,rho,a,projections=dimer_source()
    coefficients=[];rebuilt=s.zeros(16)
    for letters in product('IXYZ',repeat=4):
        label=''.join(letters);p=tensor(label);value=s.simplify(s.trace(p*h)/16)
        if value:
            if value.is_Rational is not True:raise ArithmeticError('source coefficient is not rational')
            c=F(int(value.p),int(value.q));coefficients.append((label,c));rebuilt+=value*p
    if rebuilt!=h:raise ArithmeticError('full-source Pauli reconstruction failed')
    lower,upper=root_interval();reference=F(1,2)-(lower+upper)/16
    return h,rho,a,tuple(coefficients),reference,(upper-lower)/16


@lru_cache(maxsize=1)
def _preparation():
    h,rho,_,_,_,_=_source();lower,upper=root_interval();grid=2**40
    exact=s.zeros(16,1)
    exact[3]=exact[12]=1;exact[6]=4+s.sqrt(17);exact[9]=-exact[6]
    if ((h*exact-(4-s.sqrt(17))*exact/8).applyfunc(s.simplify)!=s.zeros(16,1)
            or (exact*exact.H/(exact.H*exact)[0]-rho).applyfunc(s.simplify)!=s.zeros(16)):
        raise ArithmeticError('prepared direction is not the actual full-source ground')
    direction=[0]*16;direction[3]=direction[12]=grid
    middle=4+(lower+upper)/2
    direction[6]=(middle*grid+F(1,2)).numerator//(middle*grid+F(1,2)).denominator
    direction[9]=-direction[6]
    levels=tree.design(direction,4)
    word=[e for bit,_,co in levels for e in tree.level_word(bit,co,4)]
    raw_error=sum((sum(tree.level_error(direction,bit,raw,co,4),F())for bit,raw,co in levels),F())
    fixed=14*(2*iv.pi_angle_error(F(math.pi/2),F(1,2))+iv.pi_angle_error(F(-math.pi),F(-1)))
    port=tree.port_bound(Counter((e['kind'],e['angle'])for e in word))
    # Its four coefficients differ from (1,4+sqrt(17),-4-sqrt(17),1)
    # by less than 2^-40 each after division by grid. Normalization is 2-Lipschitz
    # relative to either vector norm, which is greater than one.
    vector_error=F(1,2**38)
    return tuple(word),raw_error+fixed+port+vector_error,tuple(direction)


def bound_template(template):
    word=circuits.emit(template);end,angle,variables=cert.consume(word,0,template)
    if end!=len(word)or variables:raise ArithmeticError('unbound source command')
    native=tree.port_bound(Counter((e['kind'],e['angle'])for e in word))
    return word,angle+native


@lru_cache(maxsize=1)
def _controls():
    _,_,_,terms,reference,reference_error=_source()
    shifted=[(label,c-reference if label=='IIII'else c)for label,c in terms]
    controlled=[(label+'I',c/2)for label,c in shifted]+[(label+'Z',-c/2)for label,c in shifted]
    template=[]
    for label,c in controlled+list(reversed(controlled)):
        template+=cert.pauli(label,-STEP*c/2)
    clock,clock_native=bound_template(template)
    norm=sum(abs(c)for _,c in shifted)
    if norm>=F(37,32):raise ArithmeticError('source norm bound changed')
    splitting=STEP**3*norm**3/3
    charge,charge_error=bound_template(cert.cnot(0,4)+[('X',(0,),cert._pi(1,2)),('I',(0,),cert._pi(-1,2))])
    hadamard,hadamard_error=bound_template(circuits.hadamard_template(4))
    # A fresh supplied blank selects stop:continue amplitudes1:16.
    # A declared binary64 control, checked against the exact 1:16 direction.
    # Its source-word digest does not depend on platform libm rounding.
    angle=float.fromhex('-0x1.82250768ac529p+0')
    coin_word=(circuits.compiler.pulse('Y',angle,[0]),)
    coin_error=cert.two_vector_error(F(angle),1,16,True)+tree.port_bound(Counter({('Y',angle):1}))
    coin=local.compile_block(local.key(coin_word),1)
    return {'clock':clock,'clock_native_error':clock_native,'clock_step_error':clock_native+splitting+STEP*reference_error,
            'splitting':splitting,'reference_error':reference_error,'norm_sum':norm,
            'charge':charge,'charge_error':charge_error,'hadamard':hadamard,'hadamard_error':hadamard_error,
            'coin':coin,'coin_error':coin_error,'coin_word':coin_word}


def _physical_slots(word):
    return tuple(2 if e['kind']=='J'else 4 for e in word)


def _failure_statistics(run,slots):
    if run.failure_weights.shape!=(len(slots),4):raise ArithmeticError('failure domain changed')
    for row,count in zip(run.failure_weights,slots,strict=True):
        if np.any(row[count:]!=0):raise ArithmeticError('padding is not a physical export')
    return float(run.failure_weights.sum()),float(run.failure_weights.max(initial=0))


@lru_cache(maxsize=1)
def _execute():
    """Execute every native pulse; retain small prefixes that reconstruct failures."""
    word,prep_error,_=_preparation();c=_controls()
    blank=np.zeros((16,1),complex);blank[0]=1
    prep=execute.complete(word,4,blank)
    tagged=np.zeros((32,1),complex);tagged[:16]=prep.success
    charge=execute.complete(c['charge'],5,tagged)
    branches=charge.success.reshape(2,16).T.copy()
    before_h=np.zeros((32,2),complex);before_h[:16]=branches
    initialized=[execute.complete(c['hadamard'],5,before_h[:,q:q+1])for q in range(2)]
    running=np.column_stack([x.success[:,0]for x in initialized])
    coin=c['coin'];stop=coin.success[0,0];continuation=coin.success[1,0]
    coin_fail=np.array([float(np.linalg.norm(x[:,0])**2)for x in coin.failures[0]])
    proper=np.zeros((COUNT,2,2));prefixes=np.zeros((COUNT+1,32,2),complex)
    failures={'preparation':prep.failure_probability,'charge':charge.failure_probability,
              'initialize':sum(x.failure_probability for x in initialized),
              'selector':0.,'readout':0.,'clock':0.}
    maximum=max(float(prep.failure_weights.max(initial=0)),float(charge.failure_weights.max(initial=0)),
                *(float(x.failure_weights.max(initial=0))for x in initialized))
    slots=_physical_slots(c['clock']);read_slots=_physical_slots(c['hadamard'])
    prefix_residual=0.
    for j in range(COUNT):
        prefixes[j]=running
        norms=np.sum(abs(running)**2,axis=0)
        selector_weights=coin_fail[:,None]*norms[None,:]
        failures['selector']+=float(selector_weights.sum());maximum=max(maximum,float(selector_weights.max(initial=0)))
        updated=[]
        for q in range(2):
            read=execute.complete(c['hadamard'],5,stop*running[:,q:q+1])
            proper[j,q]=np.sum(abs(read.success.reshape(2,16))**2,axis=1)
            total,largest=_failure_statistics(read,read_slots);failures['readout']+=total;maximum=max(maximum,largest)
            step=execute.complete(c['clock'],5,continuation*running[:,q:q+1])
            total,largest=_failure_statistics(step,slots);failures['clock']+=total;maximum=max(maximum,largest)
            prefix_residual=max(prefix_residual,read.normalization_residual,step.normalization_residual)
            updated.append(step.success[:,0])
        running=np.column_stack(updated)
    prefixes[COUNT]=running;dummy=np.sum(abs(running)**2,axis=0)
    total=float(proper.sum()+dummy.sum()+sum(failures.values()))
    if abs(total-1)>2e-8 or maximum>1e-20 or sum(failures.values())>1e-20:
        raise ArithmeticError('complete native stopped experiment lost probability or failure admission')
    for x in (proper,prefixes,dummy,before_h):x.flags.writeable=False
    return {'proper':proper,'dummy':dummy,'prefixes':prefixes,'preparation':prep,'charge':charge,
            'initialize':tuple(initialized),'before_h':before_h,'failures':failures,
            'largest_failure_weight':maximum,'normalization_residual':abs(total-1),
            'largest_stage_normalization_residual':prefix_residual,
            'signed_score':float(np.sum(proper[:,:,0]-proper[:,:,1]))}


def complete_error():
    _,prep,_=_preparation();c=_controls();tail=RATIO**COUNT
    selectors=(1-tail)/(1-RATIO);clocks=RATIO*selectors
    total=prep+c['charge_error']+c['hadamard_error']+selectors*c['coin_error']+clocks*c['clock_step_error']+(1-tail)*c['hadamard_error']
    # Failure events use only native port/angle error, not the unitary splitting.
    failure=prep**2+c['charge_error']**2+c['hadamard_error']**2+COUNT*(c['coin_error']**2+c['clock_native_error']**2+c['hadamard_error']**2)
    if total>=F(1,100000)or failure>=F(1,10**20):raise ArithmeticError('complete native error budgets failed')
    return {'complete_output_halftrace_upper':F(1,100000),
            'all_native_failure_probability_upper':F(1,10**20),
            'weighted_clock_steps_upper':F(256),'weighted_selector_uses_upper':F(257),
            'actual_rational_error':total,'actual_failure_error':failure}


def _outward(interval,grid=10**12):
    lo,hi=interval
    return F((lo*grid).__floor__(),grid),F((hi*grid).__ceil__(),grid)


def _product(x,y):
    products=[a*b for a in x for b in y]
    return min(products),max(products)


@lru_cache(maxsize=1)
def _ideal_records():
    from .local_spectral import _phase_approximation
    lo,hi=root_interval();weights=((F(1,4)+lo/68,F(1,4)+hi/68),
                                 (F(1,4)-hi/68,F(1,4)-lo/68))
    energies=((lo/8-F(1,16),hi/8-F(1,16)),(lo/8+F(1,16),hi/8+F(1,16)))
    scale=2**60;unit=[]
    for low,high in energies:
        re,im,error=_phase_approximation((low+high)*STEP/2,(high-low)*STEP/2)
        unit.append((int(re*scale),int(im*scale),error))
    phases=[(scale,0),(scale,0)];weight=1-RATIO;score=[F(),F()];rows=[]
    for j in range(COUNT):
        char=[F(),F()]
        for (x,y),(_,_,error),w in zip(phases,unit,weights,strict=True):
            if x*x+y*y>4*scale*scale:raise ArithmeticError('phase recurrence norm exceeds two')
            e=j*(2*error+F(1,scale));term=_product(w,(F(x,scale)-e,F(x,scale)+e))
            char[0]+=term[0];char[1]+=term[1]
        plus=(weight*(F(1,4)+char[0]/2),weight*(F(1,4)+char[1]/2))
        minus=(weight*(F(1,4)-char[1]/2),weight*(F(1,4)-char[0]/2))
        rows.append((plus,minus));score[0]+=2*weight*char[0];score[1]+=2*weight*char[1]
        updated=[]
        for (x,y),(re,im,_)in zip(phases,unit,strict=True):
            updated.append(((x*re-y*im+scale//2)//scale,(x*im+y*re+scale//2)//scale))
        phases=updated;weight*=RATIO
    return tuple(rows),_outward(score),RATIO**COUNT/2


def native_failure_record(stage,time_index,charge,event,slot):
    """Reconstruct one failure output vector, including exact zero prefixes."""
    run=_execute();c=_controls()
    if type(event)is not int or type(slot)is not int:raise ValueError('integer source record required')
    if stage in('preparation','charge'):
        if time_index is not None or charge is not None:raise ValueError('this record precedes charge readout')
        current=run[stage];word=_preparation()[0]if stage=='preparation'else c['charge']
    elif stage=='initialize':
        if time_index is not None or type(charge)is not int or not 0<=charge<6:raise ValueError('raw charge initialization record required')
        word=c['hadamard']
        current=run['initialize'][charge]if charge<2 else None
    else:
        if type(time_index)is not int or not 0<=time_index<COUNT or type(charge)is not int or not 0<=charge<6:raise ValueError('retained time and raw charge required')
        state=run['prefixes'][time_index,:,charge:charge+1]if charge<2 else np.zeros((32,1),complex)
        if stage=='selector':
            if event!=0 or not 0<=slot<4:raise ValueError('physical selector complement required')
            vector=c['coin'].failures[0][slot][:,0,None]*state.T
            return np.array(vector,copy=True).reshape(-1,1)
        if stage not in('clock','readout'):raise ValueError('actual source stage required')
        word=c['clock']if stage=='clock'else c['hadamard']
        amplitude=c['coin'].success[1 if stage=='clock'else 0,0]
        current=execute.complete(word,5,amplitude*state)if charge<2 else None
    if current is None:
        if not 0<=event<len(word)or not 0<=slot<_physical_slots(word)[event]:
            raise ValueError('actual source event and physical complement required')
        return np.zeros((2**(5-len(word[event]['targets'])),1),complex)
    return np.array(execute.failure(word,current,event,slot),copy=True)


def terminal_record_interval(time_index,charge_raw,clock_raw):
    """Original native probability enclosure for the complete raw record domain."""
    if type(charge_raw)is not int or not 0<=charge_raw<6:
        raise ValueError('one of six raw charge labels required')
    rows,_,dummy=_ideal_records();cap=complete_error()['complete_output_halftrace_upper']
    if time_index is None:
        if clock_raw is not None:raise ValueError('dummy clock remains unmeasured')
        interval=(dummy,dummy)
    else:
        if type(time_index)is not int or not 0<=time_index<COUNT:
            raise ValueError('actual retained time required')
        if type(clock_raw)is not int or not 0<=clock_raw<6:
            raise ValueError('one of six raw clock labels required')
        if clock_raw>=2:return F(),F()
        interval=rows[time_index][clock_raw]
    if charge_raw>=2:return F(),F()
    return _outward((max(F(),interval[0]-cap),min(F(1),interval[1]+cap)))


def _certificate():
    run=_execute();c=_controls();errors=complete_error();ideal,score,dummy=_ideal_records();cap=errors['complete_output_halftrace_upper']
    records=[]
    for j,row in enumerate(ideal):
        for q in range(2):
            for b,interval in enumerate(row):
                lo,hi=_outward((max(F(),interval[0]-cap),min(F(1),interval[1]+cap)))
                actual=run['proper'][j,q,b]
                if not float(lo)-2e-11<=actual<=float(hi)+2e-11:raise ArithmeticError('native proper record missed the complete source bound')
                records.append({'time_index':j,'charge':1 if q==0 else-1,'clock_raw':b,'signed_score':1 if b==0 else-1,'original_probability_interval':(lo,hi)})
    for q in range(2):
        box=_outward((max(F(),dummy-cap),min(F(1),dummy+cap)))
        if not float(box[0])-2e-11<=run['dummy'][q]<=float(box[1])+2e-11:raise ArithmeticError('native dummy record missed source bound')
        records.append({'time_index':None,'charge':1 if q==0 else-1,'clock_raw':None,'signed_score':0,'original_probability_interval':box})
    native_score=_outward((score[0]-2*cap,score[1]+2*cap))
    if not float(native_score[0])-2e-11<=run['signed_score']<=float(native_score[1])+2e-11:raise ArithmeticError('complete native score changed')
    time_error=STEP*F(9,16)/2+RATIO**COUNT
    total_error=2*cap+time_error
    original=_outward((max(F(),native_score[0]-total_error),min(F(1),native_score[1]+total_error)))
    window=max(F(),2*original[0]-1)
    # Independently reconstruct the exact full-source Lorentzian score.
    _,rho,a,projections=dimer_source()
    exact=s.simplify(sum(s.trace(rho*op.H*p*op)/(1+energy**2)for op in(a.H,a)for energy,p in projections))
    from .spectral_window_access import _radical_interval
    exact_box=_radical_interval(exact)
    if not original[0]<=exact_box[0]<=exact_box[1]<=original[1]:raise ArithmeticError('native score misses exact same-source spectral law')
    if window<=F(1,2):raise ArithmeticError('native original-window floor failed')
    programs={'preparation':_preparation()[0],'charge':c['charge'],'initialize':c['hadamard'],'readout':c['hadamard'],'clock':c['clock'],'selector':c['coin_word']}
    program_rows={name:{'pulse_count':len(word),'word_sha256':cert.word_digest(word),'physical_failure_slots':sum(_physical_slots(word)),
                        'event_failure_slot_counts':_physical_slots(word)}for name,word in programs.items()}
    before_charge=program_rows['preparation']['physical_failure_slots']+program_rows['charge']['physical_failure_slots']
    per_charge=program_rows['initialize']['physical_failure_slots']+COUNT*sum(program_rows[k]['physical_failure_slots']for k in('selector','readout','clock'))
    total_failure_records=before_charge+6*per_charge
    return {'source_sites':2,'full_CAR_rows':16,'source_c':F(1,8),'source_U':F(1),
            'same_centered_source_excitation_clock':True,'time_step':STEP,'time_count':COUNT,'ideal_geometric_ratio':RATIO,
            'all_proper_terminal_records':records,'proper_terminal_record_count':len(records),
            'raw_pointer_terminal_record_count':36*COUNT+6,
            'exact_zero_raw_pointer_records':{'charge_raw':(2,3,4,5),'clock_raw':(2,3,4,5),'all_indices_explicit_in_raw_product':True},
            'native_programs':program_rows,'physical_first_failure_records':total_failure_records,
            'first_failure_records_on_two_nonzero_charge_prefixes':before_charge+2*per_charge,
            'zero_charge_prefixes_follow_same_declared_controller':True,
            'all_first_failure_record_interval':(F(),F(1,10**20)),
            'all_first_failures_total_probability_interval':(F(),F(1,10**20)),
            'failure_record_reconstruction':'native_failure_record(stage,time_index,charge,event,slot)',
            'failure_prefix_rows_stored':(COUNT+1,32,2),'working_source_clock_rows':32,
            'complete_output_halftrace_upper':cap,'weighted_clock_steps_upper':F(256),'weighted_selector_uses_upper':F(257),
            'floating_complete_normalization_residual_admission':F(1,50_000_000),
            'trace_preservation_by_complete_instrument_composition':True,
            'native_signed_score_interval':native_score,'exact_ground_Lorentzian_interval':exact_box,
            'score_enclosure_method':'exact rational ideal spectral recurrence plus complete native channel bound',
            'floating_contraction_diagnostics_used_as_interval_premise':False,
            'original_Lorentzian_interval_from_native_records':original,'original_total_weight_in_0_to_1_lower':window,
            'all_native_pulses_on_potentially_nonzero_prefixes_executed':True,
            'all_source_rows_and_charge_columns_retained':True,
            'postselection_or_output_renormalization_used':False,'selector_blanks_are_supplied':True,
            'dummy_clock_environment_retained':True,'large_infinite_ground_collar_executed':False,
            'hardware_shots_sampled':False,'physical_species_or_mass_assigned':False}


def _detach(value):
    if type(value) is F:
        return F(value.numerator,value.denominator)
    if type(value) is dict:
        return {key:_detach(item)for key,item in value.items()}
    if type(value) is tuple:
        return tuple(_detach(item)for item in value)
    if type(value) is list:
        return [_detach(item)for item in value]
    return value


def certificate():
    return _detach(_certificate())


if __name__=='__main__':
    import time
    start=time.monotonic();value=certificate()
    print('native proper records',value['proper_terminal_record_count'],'failure records',value['physical_first_failure_records'])
    print('score',value['native_signed_score_interval'],'Lorentzian',value['original_Lorentzian_interval_from_native_records'],'window lower',value['original_total_weight_in_0_to_1_lower'])
    print('execution',_execute()['signed_score'],_execute()['normalization_residual'],_execute()['failures'],'seconds',time.monotonic()-start)
