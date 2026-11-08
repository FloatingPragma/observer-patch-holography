"""A counted conditional phase witness with complete apparatus readout."""
from fractions import Fraction as F
import math
import numpy as np
from . import preparation as p, coupling as c


def selected_edge(terminal):
    p.integer(terminal,0,p.FULL)
    for edge in p.SEAMS:
        if ((terminal>>edge[0])^(terminal>>edge[1]))&1:return edge
    return None


def recipe_from_tree(slots,enabled):
    if type(enabled) is not bool:raise ValueError('an explicit phase-control setting required')
    weight,pairs,records=p.readout.branch(p.FAMILY,slots)
    if not weight or not pairs:
        return {'status':'ZERO_RAW_WORD','weight':weight,'records':records}
    if len(pairs)==2:
        if pairs!=((0,0),(p.FULL,p.FULL)):raise ArithmeticError('unexpected ambiguous preparation record')
        return {'status':'INCONCLUSIVE_CONSTANT','weight':weight,'records':records}
    initial,terminal=pairs[0];edge=selected_edge(terminal)
    if edge is None:raise ArithmeticError('a singleton nonconstant tree branch lacks a crossing seam')
    return {'status':'ACTIVE','weight':weight,'records':records,
            'decoded_initial':initial,'decoded_terminal':terminal,'edge':edge,
            'incoming_contrast':((terminal>>edge[0])&1)-((terminal>>edge[1])&1),
            'word':declared_word(edge,enabled)}


def ideal_prepared_branch(slots,density):
    recipe=recipe_from_tree(slots,False)
    branch=p.readout.apply_branch(p.FAMILY,slots,density)
    if recipe['status']!='ACTIVE':return recipe['status'],branch
    y=recipe['decoded_terminal'];a,b=recipe['edge'];z=y^((1<<a)|(1<<b))
    sign=1 if (y>>a)&1 else -1
    result={}
    for (i,r,j,s),value in branch.items():
        if i!=y or j!=y:raise ArithmeticError('tree preparation failed its singleton reference factorization')
        for left,right,factor in ((y,y,1),(z,z,1),(y,z,sign),(z,y,sign)):
            result[left,r,right,s]=F(value,2)*factor
    return recipe['status'],result


def ideal_detector_branches(slots,density,enabled):
    recipe=recipe_from_tree(slots,enabled)
    branch=p.readout.apply_branch(p.FAMILY,slots,density)
    if recipe['status']!='ACTIVE':return recipe['status'],((None,0,branch),)
    y=recipe['decoded_terminal'];a,b=recipe['edge'];z=y^((1<<a)|(1<<b)) if enabled else y
    before={(z,r,z,s):v for (_,r,_,s),v in branch.items()}
    rows=p.readout.cone.fixed_instrument(p.FAMILY,before,recipe['edge'],law='E')
    return 'ACTIVE',tuple((record,int(record.law=='E' and record.label==-recipe['incoming_contrast']),block)
                          for record,block in rows)


def compact_word(enabled):
    if type(enabled) is not bool:raise ValueError('an explicit phase-control setting required')
    compiler=c.occupation.compiler;beam=c.words()['beam'][1]
    phase=[compiler.pulse('I',math.pi/2 if enabled else 0.,[0]),
           compiler.pulse('Z',-math.pi/2 if enabled else 0.,[0])]
    chart={0:0,5:1,12:2,13:3,14:4,15:5}
    repair=[dict(e,targets=[chart[x] for x in e['targets']])
            for e in c.occupation.repair_word('E',1)]
    return beam+phase+compiler.inverse(beam)+repair


def declared_word(edge,enabled):
    if type(edge) is not tuple or any(type(v) is not int for v in edge) or edge not in p.SEAMS:
        raise ValueError('an original ordered internal source seam required')
    targets=[*edge,12,13,14,15]
    word=c.occupation.compiler.remap(compact_word(enabled),targets)
    return p.source.source_word(word,16)


def pointer_certificate():
    maps=c.occupation.compiler.transfer_kraus(2,coherent=False)
    if len(maps)!=6:raise ArithmeticError('all six proper-code measurement outcomes required')
    effects=[]
    for label,matrix in enumerate(maps):
        effect=matrix[:,:2].conj().T@matrix[:,:2]
        target=np.zeros((2,2),complex)
        if label<2:target[label,label]=1
        if not np.array_equal(effect,target):
            raise ArithmeticError('the existing complete computational readout changed')
        effects.append(effect)
    if not np.array_equal(sum(effects),np.eye(2)):
        raise ArithmeticError('the complete raw-bit detector loses an outcome')
    return {'raw_bit_readout_calls':4,'labels_per_bit':6,
            'complete_label_words':6**4,'logical_raw_words':16,
            'zero_complement_label_words':6**4-16,
            'zero_complement_outcomes_retained':True}


def certificate():
    c.synthesis.certify_primitive.cache_clear();rows=[];diags=[]
    for enabled in (False,True):
        word=compact_word(enabled);unique=c.synthesis.certify_word(word)
        if len(word)!=305:raise ArithmeticError('the existing complete witness pulse count changed')
        initial=np.zeros((64,4),complex);initial[:4]=np.eye(4)
        run=c.recovery.execute_complete(word,6,initial)
        bound=c.complete_bound(len(word));delta=F(bound['complete_halfdiamond_upper'])
        for initial_pair,d in ((1,1),(2,-1)):
            score=0.;mass=0.
            for coin in (0,1):
                decoded=p.source.circuit.record_values('E',coin)
                for raw in range(8):
                    block=run.success[4*(coin+2*(raw&1)+4*((raw>>1)&1)+8*((raw>>2)&1)):
                                      4*(coin+2*(raw&1)+4*((raw>>1)&1)+8*((raw>>2)&1))+4,initial_pair]
                    probability=float(np.vdot(block,block).real);mass+=probability
                    if decoded.get(raw)==-d:score+=probability
            if abs(score-int(enabled))>2*float(delta):
                raise ArithmeticError('the full stored phase witness exceeds its certified score bound')
            if abs(mass+float(run.failure_gram[initial_pair,initial_pair].real)-1)>2*float(delta):
                raise ArithmeticError('the retained apparatus output exceeds its complete gain envelope')
            diags.append(score)
        declarations=[declared_word(e,enabled) for e in p.SEAMS]
        if any(d['declared_resources']['native_pulses']!=305 for d in declarations):
            raise ArithmeticError('an addressed witness changed its counted resources')
        rows.append({'phase_enabled':enabled,'compact_complete_input_columns':4,
                     'certified_primitive_types':unique,'original_raw_slots':25,
                     'compiled_E_raw_addresses':16,'zero_C_and_hold_addresses':9,
                     'source_word_digests':[p.source.digest(d) for d in declarations],
                     'source_resources':declarations[0]['declared_resources'],'error':bound})
    delta=F(rows[0]['error']['complete_halfdiamond_upper'])
    if 4*delta>=F(1,10**7):raise ArithmeticError('the unconditional contrast envelope is insufficient')
    return {'settings':rows,'complete_compiled_pointer':pointer_certificate(),'local_terminal_configurations':4096,
            'active_terminal_configurations':sum(selected_edge(x) is not None for x in range(4096)),
            'inconclusive_constant_configurations':2,
            'ideal_score_increment_per_active_tree_mass':'1',
            'stored_map_score_increment_lower_per_active_tree_mass':'9999999/10000000',
            'native_E_tree_attempts':11,'coherent_control_pulses':46,
            'compiled_E_readout_pulses':259,'total_source_pulses':305,
            'compiled_readout_is_not_extra_unaccounted_native_attempt':True,
            'full_work_number_leakage_and_raw_aborts_retained':True,
            'rounded_work_returned_as_exact_native_state':False,
            'constant_branches_postselected':False}
