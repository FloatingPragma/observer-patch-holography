"""Complete controlled source preparation, transfer evolution and record scoring."""
from __future__ import annotations
from fractions import Fraction as F
import numpy as np
from model.native_coexistent_effects import menus
from model.native_effect_synthesis import synthesis as syn, intervals as iv
from model.native_effect_coverage import detector
from model.native_record_response import response, transfer
from . import spectral

TIME=F(1,10)
TARGET=65
GROUND_PULSES=12309
PHASE_PULSES=11270
EVOLUTION_WORDS=66


def word_bound(pulses,grid_units):
    if type(pulses) is not int or not 0<pulses<=syn.interpreter.MAX_WORD or type(grid_units) is not int or not 0<=grid_units<=24:
        raise ValueError('an admitted bounded source word and grid allowance are required')
    x=pulses*syn.PULSE_GATE;e=x/(1-x)+F(grid_units,iv.AQ)
    fail=pulses*syn.PULSE_GATE**2/(1-x)**2
    return {'halfdiamond':(e*(2+e)+fail)/2,'gain':1/(1-x)**2+fail,
            'operator_error':e,'failure':fail}


def compose(bounds):
    delta=F();gain=F(1)
    for b in bounds:
        if b['halfdiamond']<0 or b['gain']<1:raise ValueError('complete error and gain allowances are required')
        delta+=gain*b['halfdiamond'];gain*=b['gain']
    return delta,gain


def evolution_words(time=TIME):
    time=spectral.time_parameter(time);columns,values=spectral.load()
    spectral.certify(columns,values)
    for column,value in zip(columns,values,strict=True):
        vector=tuple((x,F()) for x in column)
        yield menus.phase_word(vector,F(),spectral.phase_grid(value,time))


def ground_word(terminal):
    if type(terminal) is not int or not 0<=terminal<66:
        raise ValueError('a decoded terminal source configuration is required')
    _,pi=spectral.source_intervals()
    weights=[F()]*512
    for code,p in zip(syn.tomo.codes(),pi,strict=True):weights[code]=p
    # The literal work address zero can be padding. Export occurs after the whole word.
    code=syn.tomo.codes()[terminal]
    route=syn.tomo.route_word(tuple(('X',bit) for bit in range(9) if code>>bit&1))
    word=route+syn._state_word(weights)
    if len(word)>GROUND_PULSES:raise ArithmeticError('the source ground preparation resource changed')
    syn.interpreter.source_word(word,9)
    return word


def intervention_word(enabled=True):
    if type(enabled) is not bool:raise ValueError('a Boolean phase intervention is required')
    phases=[0]*512
    if enabled:phases[syn.tomo.codes()[TARGET]]=iv.pi_grid()//2
    word=syn._phase_word(phases)
    if len(word)!=PHASE_PULSES:raise ArithmeticError('the diagonal intervention resource changed')
    return word


def ideal_ground_branch(slots,columns):
    """One lazy preparation branch on an arbitrary untouched reference."""
    columns=menus.source_columns(columns);weight,mapping,records=detector.branch(slots)
    if not mapping:return np.zeros_like(columns),records
    initial,terminal=next(iter(mapping.items()));_,pi=spectral.source_intervals()
    omega=np.sqrt(np.array(pi,dtype=float))
    return np.sqrt(float(weight))*omega[:,None]*columns[initial,:][None,:],records


def readout_score(slots):
    weight,mapping,records=detector.branch(slots)
    if not mapping:return 0
    initial,_=next(iter(mapping.items()))
    return int(initial==TARGET)


def complete_word_failure(prefix,word,pulse,complement,columns):
    """First processor failure after all successful complete-word prefixes."""
    current=menus.source_columns(columns)
    for previous in prefix:current,_=menus.successful_word(previous,current)
    return syn.processor_failure(word,pulse,complement,current)


def complete_padding_row(prefix,word,address,columns):
    if type(address) is not int or not 0<=address<512 or address in syn.tomo.codes():
        raise ValueError('an original padding Kraus address is required')
    current=menus.source_columns(columns)
    for previous in prefix:current,_=menus.successful_word(previous,current)
    return syn.work_row(word,address)[list(syn.tomo.codes())]@current


def protocol_budget(time=TIME):
    time=spectral.time_parameter(time)
    prep=word_bound(GROUND_PULSES,9);kick=word_bound(PHASE_PULSES,1)
    evolved=[word_bound(menus.WORD_PULSES,23)]*EVOLUTION_WORDS
    delta,gain=compose([prep,kick]+evolved)
    # A unitary operator-norm error bounds half the induced diamond error.
    delta+=gain*spectral.ideal_basis_error(time)
    if delta>=F(1,11000) or gain>=F(1001,1000):
        raise ArithmeticError('the complete finite protocol exceeds its fixed envelope')
    return {'complete_halfdiamond_strict_upper':'1/11000','complete_gain_strict_upper':'1001/1000',
            'source_words':68,'source_words_max_pulses':menus.WORD_PULSES,
            'native_pulses_upper':GROUND_PULSES+PHASE_PULSES+EVOLUTION_WORDS*menus.WORD_PULSES,
            'native_fixed_E_commands':22,'work_qubits':9,
            'per_emitted_primitive_certification_required':True,
            'all_first_processor_and_padding_aborts_retained':True}


def witness():
    data=response.data();p=spectral.fraction(data['pi'][TARGET]);h0=1-spectral.fraction(data['P'][TARGET,TARGET])
    if p!=F(625,7505) or h0!=F(32,495):raise ArithmeticError('the source witness configuration changed')
    # ||H|| <=103/200 by the rational principal-log series.
    L=F(103,200)
    log=sum((transfer.DEFICIT**n/n for n in range(1,11)),F())+transfer.logarithm_tail(10)
    if log>=L:raise ArithmeticError('the transfer generator norm allowance failed')
    # u=1-it h+r, |r|<=t² L h/2; use H²<=L H.
    lower=p*h0*(2*TIME-F(3,2)*(1+F(3,2)*TIME*L)*TIME**2*L)
    if lower<=F(1,1000):raise ArithmeticError('the ideal ground response margin failed')
    actual=lower-F(4,11000)
    if actual<=F(1,2000):raise ArithmeticError('the complete source implementation erases the witness')
    return {'target_source_index':TARGET,'target_occupation_ports':[0,1],
            'ground_probability':str(p),'generator_diagonal_lower':str(h0),
            'generator_norm_upper':str(L),'time':str(TIME),
            'ideal_intervention_probability_increment_lower':str(lower),
            'ideal_increment_strict_lower':'1/1000',
            'stored_map_score_increment_lower':str(actual),
            'complete_increment_strict_lower':'1/2000',
            'all_abort_scores':0,'native_detector_score_values':[0,1],
            'control_has_no_phase_intervention':True}


def outcomes():
    return {'preparation_tree':{'raw_slots_per_command':25,'commands':11,
                                'formal_raw_words':25**11,'live_raw_words':66*2**11},
            'coherent_words':68,'first_processor_failures':'word, pulse, original complement label',
            'padding':'one observed terminal abort per word, with 446 internal Kraus rows',
            'final_tree':{'raw_slots_per_command':25,'commands':11,'formal_raw_words':25**11},
            'preparation_records_retained':True,'terminal_source_and_reference_retained':True,
            'success_postselection':False,'zero_slots_removed':False}
