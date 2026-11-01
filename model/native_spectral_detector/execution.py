"""Full-row source execution for a classical preparation ensemble.

Every local block is the existing complete six-level native instrument.
Only classical preparation labels are summed in failure probabilities;
all labels remain separate columns in the successful quantum state.
"""
from dataclasses import dataclass
import math
import numpy as np

from model.protected_particle_readout import execute as parent
from model.reference_observation.observation import frozen
from . import contraction


ADMISSION=2e-8


@dataclass(frozen=True)
class Run:
    bits:int
    blocks:tuple
    initial:np.ndarray
    success:np.ndarray
    failure_weights:np.ndarray
    failure_probability:float
    normalization_error:float
    rounded_negative_weight_count:int
    rounded_negative_weight_total:float

    def failure(self,block_index,event,slot):
        if type(block_index) is not int or not 0<=block_index<len(self.blocks):raise ValueError('actual block required')
        instrument,targets=self.blocks[block_index]
        if (type(event) is not int or not 0<=event<len(instrument.tape) or
                type(slot) is not int or not 0<=slot<len(instrument.failures[event])):
            raise ValueError('actual native failure required; padding is not an outcome')
        state=np.array(self.initial)
        for block,old_targets in self.blocks[:block_index]:
            state=parent.apply_success(state,block,old_targets,self.bits)
        value,_,_=parent.arranged(state,targets,self.bits)
        return (instrument.failures[event][slot]@value.reshape(len(instrument.success),-1)).reshape(-1,state.shape[1])


def complete(blocks,bits,initial,progress=None):
    if type(bits) is not int or not 1<=bits<=20:raise ValueError('one through twenty complete work bits required')
    if np.ma.isMaskedArray(initial):raise ValueError('masked states forbidden')
    state=np.asarray(initial)
    if state.ndim==1:state=state[:,None]
    if (state.ndim!=2 or state.shape[0]!=1<<bits or not 1<=state.shape[1]<=32 or
            state.dtype.kind not in 'fiuc' or not np.isfinite(state).all()):
        raise ValueError('complete rows and bounded classical preparation ensemble required')
    initial=frozen(np.asarray(state,complex));state=np.array(initial)
    initial_norm=float(np.vdot(initial,initial).real)
    blocks=tuple(blocks);weights=[];rounded_count=0;rounded_total=0.
    for index,(instrument,targets) in enumerate(blocks):
        if (type(targets) is not tuple or len(targets)!=instrument.bits or len(set(targets))!=len(targets) or
                any(type(j) is not int or not 0<=j<bits for j in targets)):
            raise ValueError('distinct physical targets required')
        density=contraction.apply_inplace(state,instrument.success,targets,bits)
        current=np.einsum('pqab,ba->pq',instrument.effects,density).real
        if not np.isfinite(current).all():raise ArithmeticError('nonfinite native failure weight')
        if np.min(current,initial=0)<-1e-40:raise ArithmeticError('negative native failure weight')
        rounded_count+=int(np.count_nonzero(current<0));rounded_total+=float(-current[current<0].sum())
        weights.append(np.maximum(current,0))
        if progress is not None and (index+1)%256==0:progress(index+1,len(blocks))
    table=np.concatenate(weights,axis=0) if weights else np.zeros((0,4))
    failures=float(table.sum());error=abs(float(np.vdot(state,state).real)+failures-initial_norm)
    if not math.isfinite(error) or not math.isfinite(failures) or error>ADMISSION*max(1.,initial_norm):
        raise ArithmeticError('complete native ensemble loses probability')
    return Run(bits,blocks,initial,frozen(state),frozen(table),failures,error,rounded_count,rounded_total)


def record(name,run):
    native=sum(sum(len(x) for x in instrument.failures) for instrument,_ in run.blocks)
    return {'stage':name,'work_bits':run.bits,'classical_columns':run.success.shape[1],
            'pulses':run.failure_weights.shape[0],'native_abort_slots':native,
            'padding_slots':run.failure_weights.size-native,
            'failure_weight_array':name+'_failure_weights',
            'first_failure_probability':run.failure_probability,
            'complete_probability_residual':run.normalization_error,
            'rounded_negative_slot_count':run.rounded_negative_weight_count,
            'rounded_negative_slot_mass':run.rounded_negative_weight_total}
