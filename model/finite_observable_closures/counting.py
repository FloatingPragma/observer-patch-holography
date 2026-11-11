"""An exact entire four-read law on the actual 28680-state federation."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib

import numpy as np

from model.observer_time_fields import federation as native


@lru_cache(maxsize=2)
def exact_law(feedback):
    if type(feedback) is not bool:
        raise ValueError('explicit feedback switch required')
    data = native.source()
    remote, incident, repair = native.operators()
    d = data.denominator*2
    denominator = len(data.states)*d**4
    if denominator >= np.iinfo(np.int64).max:
        raise ArithmeticError('exact declared complete count exceeds int64')
    records = native.records(4)
    output = {record:np.zeros(len(data.states),dtype=np.int64) for record in records}
    occupancy = np.array([[int(port in state) for state in data.states] for port in data.read_ports],dtype=bool)
    def continue_word(initial,word,vector):
        step = len(word)
        if step == 4:
            output[initial,word,'complete'] = vector
            return
        output[initial,word,'timeout'] = (remote@vector)*(2*d**(3-step))
        stopped = incident@vector
        for bit in (0,1):
            selected = np.where(occupancy[step%2]==bit,stopped,0)
            controlled = repair@selected if feedback and step == 0 and bit else 2*selected
            continue_word(initial,word+(bit,),controlled)
    for initial in native.INITIAL_LABELS:
        vector = ((occupancy[0]==initial[0]) & (occupancy[1]==initial[1])).astype(np.int64)
        continue_word(initial,(),vector)
    rows = []
    digest = hashlib.sha256()
    total = 0
    for initial,word,status in records:
        vector = output[initial,word,status]
        if np.any(vector<0):
            raise ArithmeticError('a complete exact count overflowed')
        count = sum(map(int,vector))
        total += count
        digest.update(vector.astype('<i8',copy=False).tobytes())
        rows.append({'initial_field_label':list(initial),'field_word':list(word),'status':status,
                     'count':str(count),'probability':str(F(count,denominator))})
    if total != denominator:
        raise ArithmeticError('the full exact finite source law is not normalized')
    timeout = sum(F(row['probability']) for row in rows if row['status']=='timeout')
    if timeout != 1-F(65,630)**4:
        raise ArithmeticError('native attempt timeout law differs')
    feedback_weight=sum(F(row['probability']) for row in rows if feedback and row['field_word'] and row['field_word'][0]==1)
    if feedback_weight != (F(13,15120) if feedback else 0):
        raise ArithmeticError('the actual record-controlled repair cost differs')
    return output, {'feedback':feedback,'complete_source_states':len(data.states),
        'cap_per_requested_interval':1,'scheduled_reads':4,'record_count':len(rows),
        'common_denominator':str(denominator),'total_count':str(total),
        'all_terminal_state_count_sha256':digest.hexdigest(),'records':rows,
        'timeout_probability':str(timeout),'success_probability':str(1-timeout),
        'clock_resources':{'ordinary_attempts_maximum':4,'additional_feedback_attempts_maximum':int(feedback),
            'ordinary_attempts_expected':str(sum(F(65,630)**k for k in range(4))),
            'feedback_attempts_expected':str(feedback_weight),
            'feedback_repair_counts_as_incident_opportunity':True},
        'every_terminal_source_state_retained':True,'rounded_or_missing_probability':False,
        'native_occupied_site_total':2}


def calculate():
    source = native.source()
    certificate = native.transfer_certificate()
    readback = native.readback_audit()
    baseline, controlled = exact_law(False)[1], exact_law(True)[1]
    if baseline['records'] == controlled['records']:
        raise ArithmeticError('record-controlled native repair has no observable effect')
    return {'dashboard_id':'exec-incident-fields','recommended_status':'executed_theorem',
        'recommended_name':'Exact four-read field statistics in the incident clock',
        'source':{'ports':source.ports,'seams':len(source.seams),'incident':len(source.incident),
                  'remote':len(source.remote),'preparation':'uniform fixed-two, both initial field bits retained'},
        'complete_exact_executions':[baseline,controlled],
        'same_source_positive_transfer_certificate':certificate,
        'native_current_readback_certificate':readback,
        'finite_premises_closed':['literal native fair endpoint identity/swap law','all full-sector integer transfer columns','complete initial field labels and all timeout prefixes','feedback uses only its retained bit','exact complete mass and every terminal state'],
        'source_scope':'cap1 exact full-sector witness; the distinct cap256 and nonzero-density refinement records retain their own certificates',
        'continuum_or_physical_clock_claimed':False}
