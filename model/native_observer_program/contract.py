"""Strict declared carrier, apparatus and policy choices for one program."""
from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
import json
import numpy as np
from model.observer.source_discrimination import circuit
from model.common_source import record_code, parameter_selection
from model.openwave_coherent_native_architecture.source import construction_descriptor
from model.native_charged_propagation.source import carrier


@lru_cache(maxsize=1)
def source_addresses():
    g=carrier()
    return tuple(circuit.encode(tuple(int(v in pair)for v in range(12)))for pair in g.pairs)


def decoded_basis():
    return tuple((address,tuple(v for v,n in enumerate(circuit.decode(address))if n))for address in source_addresses())


@lru_cache(maxsize=1)
def _manifest_template():
    weights,timeout=parameter_selection.completion_weights()
    return {'schema':'oph.single_native_observer_program.v1',
      'hardware':construction_descriptor(),
      'carrier':{'ports':12,'load':2,'binary_configurations':66,'source_addresses':list(source_addresses()),
                 'source_register':512,'active_memory':6,'memory_register':8,
                 'uniform_events':60,'event_register':64,'source_event_memory_complements':'retained complete admission outcomes'},
      'voltage':'3 transported-order parity + 2 current when seam=(0,5) and spectator=1 modulo6',
      'policies':{
       'coherent_spatial':{'event':'uniform30 seams and binary recoded identity/swap selector, distinct from original native coin','step':'S(2|uniform><uniform|-I)',
          'raw_records':'coherent; inverse recoder is required for original-record commitment','terminal':'all66 source positions times6 memory characters plus complete complement','maximum_steps':8},
       'committed_event':{'event':'uniform30 seams, fair C/E mixture','raw_slots_per_seam':25,'feedback':'active voltage','terminal':'original native fine record and complete source-memory readout'},
       'fixed_compensated_clock':{'law':'C','selector_weights':{str(k):str(v)for k,v in weights.items()},'timeout':str(timeout),
          'phase':'exp(-i t H2/N); H2 is the declared two-hop cost',
          'recovery':'inverse active feedback then original record-selected512 recovery',
          'readout':'complete X/Y matching-code quadratures with outside outcomes',
          'time_per_block':'64','steps_per_block':4096,'source_repreparation_between_blocks':False},
       'tracked_unitary_conjugation':{'frame':'declared W^n; use W^-n before and W^n after the fixed clock',
          'status':'interface contract; its transported-code observations require the consuming protocol'}},
      'preparation':{'logical_word':'equal coherent matching labels0 and1','memory_character':1,
          'event_input':'blank event address0, followed by the compiled uniform preparation',
          'source_and_character_preparation':'declared coherent services; no spontaneous preparation claim'},
      'work':{'reference_observable':'H2 on binary source, zero on source padding; memory/event references have declared zero energy',
          'controlled_operations':'source expectation changes are supplied by apparatus; no autonomous energy conservation inferred',
          'physical_energy_unit':None},
      'alternate_preparations':{'idle_apparatus_mode':'matching0, memory character1, event(|0>-|2>)/sqrt2 by explicit two native pulses; not the clock input'},
      'integer_extension':{'link_window':4,'adapter':'native_electric_gauss.extension_descriptor; exact declared descriptor required',
          'status':'explicit added signed-pair/tag and integer-link registers; same native source/memory motion',
          'preparation':'neutral pair/string, distinct from matching-code clock preparation'},
      'physical_claims':{'rest_mass_identified':False,'physical_clock_identified':False,'OpenWave_awarded':False}}


def manifest():
    return deepcopy(_manifest_template())


def validate_manifest(value):
    if type(value)is not dict or value!=manifest():
        raise ValueError('undeclared carrier, record policy, apparatus or clock substitution')
    return deepcopy(value)


def validate_state(value):
    a=np.asarray(value)
    if a.ndim!=4 or a.shape[:3]!=(512,8,64)or not 1<=a.shape[3]<=8 or a.dtype.kind not in 'iufc' or not np.all(np.isfinite(a)):
        raise ValueError('full512x8x64 source-memory-event state with one through8 reference columns required')
    return a.astype(complex,copy=True)


def matching_vector(logical):
    if type(logical)is not int or not 0<=logical<6:raise ValueError('one of six matching-code labels required')
    out=np.zeros(512,complex)
    for state in record_code.SUPPORTS[logical]:out[circuit.encode(state)]=1/np.sqrt(6)
    return out


def memory_character(k):
    if type(k)is not int or not 0<=k<6:raise ValueError('one of six active memory characters required')
    return np.exp(-2j*np.pi*k*np.arange(6)/6)/np.sqrt(6)


def descriptor_sha256():
    return sha256(json.dumps(manifest(),sort_keys=True,separators=(',',':')).encode()).hexdigest()
