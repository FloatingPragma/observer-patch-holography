"""A complete terminal interference observable and no-signalling control."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
import sympy as s

from model.coherent_screen_observer import construction, oracle
from model.native_action import incident_refinement as clocks
from model.observer import native


def readout_law(counts, seam_count, cap, law):
    incident = sum(counts)
    pi = tuple(F(x, incident) for x in counts)
    q = F(1) if cap is None else 1-(1-F(incident,seam_count))**cap
    chi = (F(1),F(1,2),F(1,2),F(0)) if law == 'C' else (F(1),F(0),F(0),F(0))
    rows = []
    for status, word in construction.records():
        weight, coherence = F(1), F(1)
        for mask in word:
            weight *= q*pi[mask]
            coherence *= q*pi[mask]*chi[mask]
        if status == 'timeout':
            weight *= 1-q
            coherence *= 1-q
        sign = -1 if word and word[0] == 1 else 1
        # Feedback changes the raw coherence sign. Classical readout relabeling
        # corrects it while retaining the complete original record.
        raw_coherence = sign*coherence
        for bit in (-1,1):
            probability = (weight+bit*sign*raw_coherence)/2
            if probability < 0:
                raise ArithmeticError('terminal interference probability is negative')
            rows.append({'status': status, 'word': list(word), 'readout': bit,
                         'probability': str(probability)})
    if sum(F(row['probability']) for row in rows) != 1:
        raise ArithmeticError('terminal interference instrument lost mass')
    return rows


@lru_cache(maxsize=1)
def calculate():
    native_checks=0
    for bits in range(4):
        for coin in (0,1):
            loads=tuple((bits>>j)&1 for j in range(12))
            image,record=native.native(loads,(0,1),coin)
            expected=(sum(value<<j for j,value in enumerate(image)),record[1],loads[0]-loads[1])
            if oracle._image(bits,(0,1),coin)!=expected:
                raise ArithmeticError('exact oracle differs from the literal native endpoint')
            native_checks+=1
    observable = s.zeros(8)
    observable[2,5] = observable[5,2] = 1
    effects = [(s.eye(8)+sign*observable)/2 for sign in (-1,1)]
    if sum(effects,s.zeros(8)) != s.eye(8):
        raise ArithmeticError('terminal POVM is incomplete')
    # Literal square-root Kraus matrices, including the six-dimensional
    # zero-eigenvalue subspace whose outcomes are both retained.
    zero = s.eye(8)-observable**2
    kraus = [(observable**2+sign*observable)/2+zero/s.sqrt(2) for sign in (-1,1)]
    if any(k.H*k != effect for k,effect in zip(kraus,effects,strict=True)):
        raise ArithmeticError('terminal Kraus maps differ from the effects')
    edges = ((0,1),(0,2),(1,2),(2,3),(3,4))
    executions = []
    for law in ('C','E'):
        blocks = oracle.exact_history(5,edges,root=(0,1,2),law=law,cap=2)
        expected = readout_law((1,1,1,1),5,2,law)
        reference = np.full((2,2),F(0),dtype=object)
        actual = []
        for status,word in construction.records():
            block=blocks[status,word]
            sign = -1 if word and word[0] == 1 else 1
            mass = sum(block[j,j] for j in range(8))
            cross = block[2,5]+block[5,2]
            for a in range(2):
                for b in range(2):
                    reference[a,b] += sum(block[2*j+a,2*j+b] for j in range(4))
            for bit in (-1,1):
                actual.append({'status':status,'word':list(word),'readout':bit,
                               'probability':str((mass+bit*sign*cross)/2)})
        if actual != expected or reference.tolist() != [[F(1,2),F(0)],[F(0),F(1,2)]]:
            raise ArithmeticError('native reference or joint terminal observable differs')
        executions.append({'law':law,'exact_native_ports':5,'cap':2,
                           'joint_records':actual,'reference_marginal':[['1/2','0'],['0','1/2']]})
    ports,edges,incident,_,_,_ = clocks.geometry()
    counts = tuple(sum(construction.category(edge)==m for edge in incident) for m in range(4))
    if counts != (23,5,4,1) or len(ports) != 222 or len(edges) != 561:
        raise ArithmeticError('matched native screen incidence differs')
    complete = {law:readout_law(counts,len(edges),None,law) for law in ('C','E')}
    capped = {law:readout_law(counts,len(edges),384,law) for law in ('C','E')}
    plus = {law:sum(F(x['probability']) for x in rows if x['readout']==1) for law,rows in complete.items()}
    if plus != {'C':F(61,72),'E':F(809,1089)} or plus['C']-plus['E'] != F(101,968):
        raise ArithmeticError('exact completed interference contrast differs')
    return {'dashboard_id':'exec-coherent-screen',
        'recommended_name':'Joint screen-reference interference under feedback',
        'recommended_status':'executed_theorem',
        'observable':'record-corrected effects (I plus/minus O)/2, O=|10,0_R><01,1_R|+h.c.',
        'finite_terminal_Kraus_completeness_exact':True,
        'literal_native_endpoint_cases_checked':native_checks,
        'untouched_reference_marginal_unchanged':True,
        'reference_alone_distinguishes_C_and_E':False,
        'exact_small_native_executions':executions,
        'matched_screen_incident_mask_counts':list(counts),
        'completed_joint_records':complete,'capped_joint_records':capped,
        'completed_positive_probability':{law:str(value) for law,value in plus.items()},
        'completed_probability_difference':'101/968',
        'terminal_observable_depends_on_midpoint_population_table':False,
        'source_scope':'the declared Bell pair/reference and independent Bernoulli environment, C/E native maps, two incident ticks, record feedback and supplied complete terminal joint instrument',
        'physical_service_accuracy_or_physical_clock_claimed':False}
