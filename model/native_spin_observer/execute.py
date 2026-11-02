"""Streaming full-row native instruments; columns are a retained reference."""
from dataclasses import dataclass
import hashlib
import math
import numpy as np
from ..centered_native_preparation import tree
from ..native_spectral_detector import contraction
from ..protected_particle_readout import execute as local
from ..reference_observation.observation import frozen
from ..certified_reference import certificate as cert

ADMISSION=2e-8


@dataclass(frozen=True)
class Run:
    bits:int
    initial:np.ndarray
    success:np.ndarray
    failure_weights:np.ndarray
    offsets:tuple
    word_sha256:str
    native_failure_slots:int
    normalization_residual:float
    rounded_negative_count:int
    rounded_negative_mass:float

    @property
    def failure_probability(self):return float(self.failure_weights.sum())


def complete(events,bits,initial):
    if type(bits) is not int or not 1<=bits<=20:raise ValueError('1..20 complete work bits required')
    if np.ma.isMaskedArray(initial):raise ValueError('masked inputs forbidden')
    initial=np.asarray(initial)
    if initial.ndim==1:initial=initial[:,None]
    if initial.ndim!=2 or initial.shape[0]!=1<<bits or not 1<=initial.shape[1]<=64 or not np.isfinite(initial).all():
        raise ValueError('complete work rows and 1..64 retained reference columns required')
    initial=frozen(np.asarray(initial,complex));state=np.array(initial)
    weights=[];offsets=[];digest=hashlib.sha256();slots=0;rounded=0;negative=0.
    for offset,instrument,targets in tree.blocks(events):
        if max(targets)>=bits:raise ValueError('native target outside complete register')
        for (kind,angle,local_targets),failures in zip(instrument.tape,instrument.failures,strict=True):
            digest.update((kind+'|'+angle.hex()+'|'+','.join(str(targets[t]) for t in local_targets)+';').encode())
            slots+=len(failures)
        density=contraction.apply_inplace(state,instrument.success,targets,bits)
        current=np.einsum('pqab,ba->pq',instrument.effects,density).real
        if not np.isfinite(current).all() or np.min(current,initial=0)<-1e-40:raise ArithmeticError('invalid failure weight')
        rounded+=int(np.count_nonzero(current<0));negative+=float(-current[current<0].sum())
        weights.append(np.maximum(current,0));offsets.append((offset,len(instrument.tape),targets))
    table=np.concatenate(weights) if weights else np.zeros((0,4))
    residual=abs(float(np.vdot(state,state).real)+float(table.sum())-float(np.vdot(initial,initial).real))
    if not math.isfinite(residual) or residual>ADMISSION*max(1.,float(np.vdot(initial,initial).real)):
        raise ArithmeticError('full native output loses probability')
    return Run(bits,initial,frozen(state),frozen(table),tuple(offsets),digest.hexdigest(),slots,residual,rounded,negative)


def failure(events,run,event,slot):
    if type(event) is not int or not 0<=event<len(run.failure_weights):raise ValueError('actual event required')
    events=tuple(events)
    if cert.word_digest(events)!=run.word_sha256:raise ValueError('failure program changed')
    state=np.array(run.initial)
    for offset,instrument,targets in tree.blocks(events):
        if offset<=event<offset+len(instrument.tape):
            j=event-offset
            if type(slot) is not int or not 0<=slot<len(instrument.failures[j]):raise ValueError('actual complement slot required')
            a,_,_=local.arranged(state,targets,run.bits)
            return (instrument.failures[j][slot]@a.reshape(len(instrument.success),-1)).reshape(-1,state.shape[1])
        contraction.apply_inplace(state,instrument.success,targets,run.bits)
    raise ArithmeticError('missing source event')
