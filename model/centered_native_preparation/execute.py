"""Stream exact local native instruments over every full-register row."""
from dataclasses import dataclass
import math
import hashlib
import numpy as np

from ..native_spectral_detector import contraction
from ..protected_particle_readout import execute as local
from ..reference_observation.observation import frozen
from . import tree

ADMISSION=2e-8


@dataclass(frozen=True)
class Run:
    bits:int
    initial:np.ndarray
    success:np.ndarray
    failure_weights:np.ndarray
    schedule:tuple
    failure_probability:float
    normalization_error:float
    rounded_negative_count:int
    rounded_negative_mass:float
    word_sha256:str
    native_failure_slots:int


def complete(events,bits,initial,progress=None):
    if type(bits) is not int or not 1<=bits<=20:raise ValueError('1..20 complete work bits required')
    initial=np.asarray(initial)
    if initial.ndim==1:initial=initial[:,None]
    if initial.ndim!=2 or initial.shape[0]!=1<<bits or not 1<=initial.shape[1]<=4 or not np.isfinite(initial).all():
        raise ValueError('complete field rows and 1..4 columns required')
    initial=frozen(np.asarray(initial,complex));state=np.array(initial);weights=[];schedule=[];rounded=0;negative=0.
    digest=hashlib.sha256();native_slots=0
    for block,(offset,instrument,targets) in enumerate(tree.blocks(events)):
        if max(targets)>=bits:raise ValueError('native target exceeds complete register')
        for (kind,angle,local_targets),failures in zip(instrument.tape,instrument.failures,strict=True):
            digest.update((kind+'|'+angle.hex()+'|'+','.join(str(targets[t]) for t in local_targets)+';').encode())
            native_slots+=len(failures)
        density=contraction.apply_inplace(state,instrument.success,targets,bits)
        current=np.einsum('pqab,ba->pq',instrument.effects,density).real
        if not np.isfinite(current).all() or np.min(current,initial=0)<-1e-40:raise ArithmeticError('invalid native first-failure weight')
        rounded+=int(np.count_nonzero(current<0));negative+=float(-current[current<0].sum())
        weights.append(np.maximum(current,0));schedule.append((offset,len(instrument.tape),targets))
        if progress is not None and (block+1)%256==0:progress(block+1,offset+len(instrument.tape))
    table=np.concatenate(weights) if weights else np.zeros((0,4))
    failures=float(table.sum());error=abs(float(np.vdot(state,state).real)+failures-float(np.vdot(initial,initial).real))
    if not math.isfinite(error) or error>ADMISSION*max(1.,float(np.vdot(initial,initial).real)):
        raise ArithmeticError('complete source probability accounting failed')
    return Run(bits,initial,frozen(state),frozen(table),tuple(schedule),failures,error,rounded,negative,digest.hexdigest(),native_slots)


def failure(events,run,event,slot):
    """Reconstruct a full native failure vector from its unnormalized prefix."""
    if type(event) is not int or not 0<=event<len(run.failure_weights):raise ValueError('actual source event required')
    from ..certified_reference.certificate import word_digest
    events=tuple(events)
    if word_digest(events)!=run.word_sha256:raise ValueError('retained failure program changed')
    state=np.array(run.initial)
    for offset,instrument,targets in tree.blocks(events):
        if offset<=event<offset+len(instrument.tape):
            local_event=event-offset
            if type(slot) is not int or not 0<=slot<len(instrument.failures[local_event]):raise ValueError('actual native complement slot required')
            value,_,_=local.arranged(state,targets,run.bits)
            return instrument.failures[local_event][slot]@value.reshape(len(instrument.success),-1)
        contraction.apply_inplace(state,instrument.success,targets,run.bits)
    raise ArithmeticError('event did not appear in the retained native word')
