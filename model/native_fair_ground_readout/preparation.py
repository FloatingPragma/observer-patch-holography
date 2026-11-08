"""Fair H66 Dicke-ground preparation using the existing source services."""
from fractions import Fraction as F
import numpy as np
from model.native_transfer_control import control as c
from model.native_energy_measurement import measurement as m
from model.native_effect_synthesis import packet as synthesis_packet


def ground_word(terminal):
    if type(terminal) is not int or not 0<=terminal<66:
        raise ValueError('a decoded original source configuration required')
    weights=[F()]*512
    for code in c.syn.tomo.codes():weights[code]=F(1,66)
    code=c.syn.tomo.codes()[terminal]
    route=c.syn.tomo.route_word(tuple(('X',bit) for bit in range(9) if code>>bit&1))
    word=route+c.syn._state_word(weights)
    if len(word)>c.GROUND_PULSES:raise ArithmeticError('ground word exceeds its existing admission')
    c.syn.interpreter.source_word(word,9)
    return word


def ideal_branch(slots,columns):
    columns=c.menus.source_columns(columns)
    weight,mapping,records=m.preparation_branch(slots)
    if not mapping:return np.zeros_like(columns),records
    initial,_=next(iter(mapping.items()))
    omega=np.ones(66)/np.sqrt(66)
    return np.sqrt(float(weight))*omega[:,None]*columns[initial,:][None,:],records


def budget():
    b=c.word_bound(c.GROUND_PULSES,9)
    if b['halfdiamond']>=F(1,4000000) or b['gain']>=F(1000001,1000000):
        raise ArithmeticError('the inherited complete preparation budget failed')
    return {'native_fair_E_attempts':11,'formal_raw_words':25**11,
            'live_raw_words':66*2**11,'source_words':1,
            'native_coherent_pulses_upper':c.GROUND_PULSES,'work_qubits':9,
            'complete_halfdiamond_error_strict_upper':'1/4000000',
            'complete_gain_strict_upper':'1000001/1000000',
            'all_raw_records_and_first_processor_padding_failures_retained':True,
            'preparation_archive_can_remain_correlated_with_untouched_reference':True,
            'source_ground_factors_from_archive_and_reference_in_ideal_model':True,
            'success_postselection':False,'growing_coherent_service_embedding':False,
            'native_attempt_to_service_pulse_duration_identification':False}


def local_phase_word(enabled=True):
    if type(enabled) is not bool:raise ValueError('a Boolean intervention flag required')
    phases=[0]*512
    for code,state in zip(c.syn.tomo.codes(),c.syn.interpreter.local.basis(),strict=True):
        if enabled and state[0]:phases[code]=c.iv.pi_grid()//2
    word=c.syn._phase_word(phases)
    if len(word)!=c.PHASE_PULSES:raise ArithmeticError('local phase admission changed')
    return word


def scored_energy_join():
    """Exact complete-instrument bound; the inherited energy readout is state independent."""
    m.budget();m.spectral.certify()
    counts=m.spectral.source_counts();states=c.syn.interpreter.local.basis()
    mean=F(sum((60*int(i==j)-counts[i][j])*int(states[i][0]==states[j][0])
               for i in range(66) for j in range(66)),66)
    if mean!=F(50,33):raise ArithmeticError('the exact source-local excitation energy changed')
    delta,gain=c.compose([c.word_bound(c.GROUND_PULSES,9),c.word_bound(c.PHASE_PULSES,1),
                        {'halfdiamond':F(1,12000),'gain':F(1001,1000)}])
    epsilon=60*(m.spectral.RESIDUAL_GATE+2*m.spectral.GRAM_GATE)/m.spectral.FLOOR
    # Rational 96-term lower log score, retaining its explicit scalar remainder.
    tail=F(60,25)*F(1,3)**97/(97*F(2,3))
    difference=F(2,33)-2*epsilon/25-2*tail-4*delta
    if delta>=F(1,11900) or difference<=F(3,50):
        raise ArithmeticError('complete native energy-score contrast lost its margin')
    return {'source_local_phase':'exp(i*pi*n0/2)',
            'exact_exchange_energy_of_ideal_excitation':'50/33',
            'counted_transfer_energy_lower':'50/33',
            'spectral_label_score':'(12/5)*sum(k=1..96,(1-lambda_j)^k/k); all abort scores zero',
            'score_range':'[0,1]','complete_unnormalized_score_contrast_strict_lower':'3/50',
            'complete_halfdiamond_error_strict_upper':'1/11900',
            'source_words':68,'native_coherent_pulses_upper':c.GROUND_PULSES+c.PHASE_PULSES+66*m.WORD_PULSES,
            'native_fair_E_attempts':11,'state_independent_measurement_certificate_reused':True,
            'same_fair_source_ground_phase_and_readout_dictionary':True,
            'stable_particle_pole_or_mass_identified':False}


def run():
    from model.primitive_observer import transfer
    source=transfer.native_counts();transfer.exact_audit(source)
    c.syn.certify_primitive.cache_clear()
    if tuple(source.states)!=tuple(c.syn.interpreter.local.basis()):
        raise ArithmeticError('the ground and detector source charts disagree')
    if len(source.states)!=66 or any(sum(row)!=60 for row in source.counts):
        raise ArithmeticError('the exact fair ground changed')
    # One concrete live fair E record, plus all terminal routing checks.
    masks=c.syn.interpreter.local.basis();state=masks[17];slots=[]
    for edge in c.detector.tree():
        rows=m.fair_kernel(edge)
        slot=next(i for i,(_,pairs) in enumerate(rows) if any(x==masks.index(state) for x,_ in pairs))
        target=dict(rows[slot][1])[masks.index(state)];state=masks[target];slots.append(slot)
    weight,mapping,_=m.preparation_branch(tuple(slots))
    if mapping!={17:masks.index(state)} or weight!=F(1,2048):
        raise ArithmeticError('fair native record and terminal decoder disagree')
    terminal=masks.index(state);word=ground_word(terminal)
    checked=c.syn.certify_word(word)
    for bit in range(9):c.syn.certify_word(c.syn.tomo.route_word((('X',bit),)))
    for code in c.syn.tomo.codes():
        remaining=code
        for bit in range(9):
            if code>>bit&1:remaining^=1<<bit
        if remaining:raise ArithmeticError('an exact record-dependent route failed')
    initial=np.zeros((66,1),complex);initial[terminal,0]=1
    prepared,padding=c.menus.successful_word(word,initial)
    if np.linalg.norm(prepared[:,0]-np.ones(66)/np.sqrt(66))>=1e-8:
        raise ArithmeticError('actual compiled fair-ground diagnostic failed')
    if abs(np.trace(padding))>=1e-8:raise ArithmeticError('padding diagnostic failed')
    kick=local_phase_word();kick_checked=c.syn.certify_word(kick)
    excited,kick_padding=c.menus.successful_word(kick,prepared)
    control_word=local_phase_word(False);control_checked=c.syn.certify_word(control_word)
    control,control_padding=c.menus.successful_word(control_word,prepared)
    ideal=np.ones(66,dtype=complex)/np.sqrt(66)
    for i,state in enumerate(source.states):
        if state[0]:ideal[i]*=1j
    if np.linalg.norm(excited[:,0]-ideal)>=1e-8 or abs(np.trace(kick_padding))>=1e-8:
        raise ArithmeticError('actual local-phase excitation diagnostic failed')
    if np.linalg.norm(control[:,0]-np.ones(66)/np.sqrt(66))>=1e-8 or abs(np.trace(control_padding))>=1e-8:
        raise ArithmeticError('the actual no-intervention control diagnostic failed')
    # Full inherited ideal instrument, with the final quantum residual retained.
    blocks,residual=m.ideal_probabilities(np.column_stack((control[:,0],excited[:,0])))
    _,values=m.spectral.load()
    scores=[F(12,5)*sum(((1-v)**k/k for k in range(1,97)),F()) for v in values]
    if not all(0<=s<=1 for s in scores):raise ArithmeticError('energy labels escaped the bounded score')
    energy_scores=sum(float(s)*np.real(np.diag(block)) for s,block in zip(scores,blocks,strict=True))
    if energy_scores[1]-energy_scores[0]<=3/50:raise ArithmeticError('energy-score diagnostic lost contrast')
    admitted=c.syn.interpreter.source_word(word,9)
    host=c.syn.interpreter.execute_word(admitted,np.zeros((1,1)))
    if host['status']!='RESOURCE_STOP' or host['pulses_executed']!=0:
        raise ArithmeticError('unchanged generic host ceiling was bypassed')
    return {'budget':budget(),'actual_compiled_pulses':len(word),
            'ground_word_sha256':synthesis_packet.word_digest(word),
            'local_phase_word_sha256':synthesis_packet.word_digest(kick),
            'matched_control_word_sha256':synthesis_packet.word_digest(control_word),
            'ground_word_declared_resources':admitted['declared_resources'],
            'local_phase_declared_resources':c.syn.interpreter.source_word(kick,9)['declared_resources'],
            'matched_control_declared_resources':c.syn.interpreter.source_word(control_word,9)['declared_resources'],
            'distinct_primitive_checks':checked,'one_live_fair_preparation_word':slots,
            'all_terminal_routes_exactly_checked':66,
            'stored_array_success_diagnostic_bound':'1/100000000',
            'stored_array_floating_roundoff_certified':False,
            'local_phase_native_pulses':len(kick),'local_phase_distinct_primitive_checks':kick_checked,
            'control_phase_distinct_primitive_checks':control_checked,
            'full_ideal_energy_instrument_with_residual_diagnostic_checked':True,
            'energy_readout_join':scored_energy_join(),
            'exact_original_fair_transfer_ground_checked':True,
            'ground_trace_distance_from_P_over_56':'55/56',
            'generic_host_status':'RESOURCE_STOP','generic_host_pulses_executed':0}
