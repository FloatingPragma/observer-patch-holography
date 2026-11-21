"""Native seam gates, retained reader histories and a time-independent controller.

The finite controller is a declared extension. Its energy is not the load
Hamiltonian, and its coherent counter is not a monotone proper-time clock.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import sparse
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import expm_multiply

from model.executable_model.carrier import cell
from model.candidate_laws.contrast_memory_spherical_v1.law import SRC, seam_term
from model.candidate_laws.contrast_memory_spherical_v1.reader import APP, CONTROLS, stages
from model.candidate_laws.contrast_memory_spherical_v1.sectors import MultisetSector


def finite(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f'{name}: expected real number')
    if np.iscomplexobj(value):
        raise ValueError(f'{name}: expected real number')
    try:
        value=float(value)
    except (OverflowError,ValueError) as exc:
        raise ValueError(f'{name}: invalid number') from exc
    if not math.isfinite(value):
        raise ValueError(f'{name}: nonfinite')
    return float(value)


def positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f'{name}: expected positive integer')
    return value


def source_sector():
    return np.asarray(MultisetSector(cell(), 2).states, dtype=int)


def native_gate(kappa, angle):
    """All nine input columns of the exact blank-memory reduced event."""
    kappa, angle = finite(kappa, 'kappa'), finite(angle, 'angle')
    if not math.isfinite(angle * kappa * 4):
        raise ValueError('native contrast phase is not representable')
    swap = np.zeros((9, 9))
    for i, (a, b) in enumerate(SRC):
        swap[SRC.index((b, a)), i] = 1
    phase = np.exp(-1j * angle)
    exchange = (1 + phase) / 2 * np.eye(9) + (1 - phase) / 2 * swap
    contrast = np.array([(a-b)**2 for a,b in SRC])
    return np.exp(-1j * angle * kappa * contrast)[:, None] * exchange


def embed_local(matrix, seam, states):
    """Embed on the complete 78-state total-two basis of twelve qutrits."""
    matrix = np.asarray(matrix)
    if matrix.shape != (9,9) or not np.isfinite(matrix).all():
        raise ValueError('local operator must be finite 9 by 9')
    states = np.asarray(states)
    if states.shape != (78,12) or not np.issubdtype(states.dtype, np.integer):
        raise ValueError('source basis must contain 78 integer rows of twelve loads')
    if np.any(states < 0) or np.any(states > 2) or np.any(states.sum(axis=1) != 2):
        raise ValueError('source basis must have qutrit loads totaling two in every row')
    a,b = seam
    if a == b or not (0 <= a < 12 and 0 <= b < 12):
        raise ValueError('invalid seam')
    lookup = {tuple(row): i for i,row in enumerate(states)}
    if len(lookup) != len(states):
        raise ValueError('source basis contains duplicate rows')
    out = np.zeros((len(states),len(states)), complex)
    for j,row in enumerate(states):
        col = SRC.index((row[a],row[b]))
        for i,(x,y) in enumerate(SRC):
            if matrix[i,col] == 0:
                continue
            target = row.copy()
            target[a],target[b] = x,y
            if tuple(target) not in lookup:
                raise ValueError('operator leaves declared total-two sector')
            out[lookup[tuple(target)],j] += matrix[i,col]
    return out


def event_matrix(seam, kappa, angle):
    return embed_local(native_gate(kappa,angle),seam,source_sector())


@dataclass(frozen=True)
class Step:
    kind: str
    ports: tuple[int, ...] = ()
    angle: float = 0.0
    reader: int = -1
    stage: int = -1


def program(rounds=8, angle=1/8):
    """Two readings and symmetric products of actual native seam/reference terms."""
    positive_integer(rounds,'rounds')
    finite(angle,'angle')
    terms = [('seam',tuple(e)) for e in cell().seams] + [('reference',(v,)) for v in range(12)]
    period = [(kind,ports,angle/2) for kind,ports in terms[:-1]]
    period += [(terms[-1][0],terms[-1][1],angle)]
    period += [(kind,ports,angle/2) for kind,ports in reversed(terms[:-1])]
    result = []
    for r,seam in enumerate(((0,1),(0,5))):
        for _ in range(rounds):
            result.extend(Step(kind,ports,a) for kind,ports,a in period)
        result.extend(Step('reader',seam,reader=r,stage=j) for j in range(4))
    return tuple(result)


def reader_permutations(control):
    if control not in CONTROLS:
        raise ValueError('unknown apparatus control')
    out = [stages(flag,**CONTROLS[control]) for flag in (0,1)]
    words = np.arange(APP)
    for by_flag in out:
        for p in by_flag:
            if not np.array_equal(np.sort(p),words) or not np.array_equal(p[p],words):
                raise ValueError('reader is not a reversible involution')
    return out


def gate_action(step, blocks, kappa, states, permutations):
    """Apply the actual gate to coherent archive-word blocks (vectors or matrices).

    Blocks with the same archive word are added as amplitudes. No measurement,
    renormalization, reset, or pruning by probability occurs.
    """
    result = {}
    if step.kind == 'seam':
        a,b = step.ports
        lookup = {tuple(row):i for i,row in enumerate(states)}
        swap = np.array([lookup[tuple(row[[b if v==a else a if v==b else v for v in range(12)]])]
                         for row in states])
        phase = np.exp(-1j * step.angle)
        diag = np.exp(-1j*step.angle*kappa*(states[:,a]-states[:,b])**2)
        for words,block in blocks.items():
            shaped = diag.reshape((-1,)+(1,)*(block.ndim-1))
            result[words] = shaped*((1+phase)/2*block+(1-phase)/2*block[swap])
    elif step.kind == 'reference':
        diag = np.exp(-1j*step.angle*kappa*states[:,step.ports[0]]**2)
        for words,block in blocks.items():
            result[words] = diag.reshape((-1,)+(1,)*(block.ndim-1))*block
    elif step.kind == 'reader':
        a,b = step.ports
        flag = (states[:,a]==1)&(states[:,b]==1)
        for words,block in blocks.items():
            for f in (0,1):
                key = list(words)
                key[step.reader] = int(permutations[f][step.stage][words[step.reader]])
                key = tuple(key)
                selected = (flag==f).reshape((-1,)+(1,)*(block.ndim-1))*block
                result[key] = result.get(key,0)+selected
    else:
        raise ValueError('unknown program operation')
    return result


def count_table(tape):
    totals=np.zeros(12,dtype=int)
    for step in tape:
        if step.kind in ('seam','reference'):
            totals[list(step.ports)]+=1
    capacities=tuple(int(x)+1 for x in totals)
    counts = np.zeros((len(tape)+1,12),dtype=int)
    for k,step in enumerate(tape):
        counts[k+1] = counter_gate(tuple(int(x) for x in counts[k]),step,capacities)
    return counts


def counter_gate(word,step,capacities,inverse=False):
    """Local modular addition, including dirty counter words and its inverse.

    A controller edge carries this permutation tensor the source/reader gate.
    Along the single-clock history the counters have one common basis word per
    prefix; count_table executes this permutation in that exact representation.
    """
    if len(word)!=12 or len(capacities)!=12:
        raise ValueError('twelve counter registers required')
    if type(inverse) is not bool:
        raise ValueError('inverse must be boolean')
    for x,cap in zip(word,capacities):
        positive_integer(cap,'capacity')
        if type(x) is not int or not 0<=x<cap:
            raise ValueError('counter outside its declared register')
    result=list(word)
    if step.kind in ('seam','reference'):
        for port in step.ports:
            result[port]=(result[port]+(-1 if inverse else 1))%capacities[port]
    return tuple(result)


def run_history(tape, kappa, control='nominal', full_operator=False):
    states = source_sector()
    selected = np.zeros(12,dtype=int)
    selected[[0,1]]=1
    start = int(np.flatnonzero(np.all(states==selected,axis=1))[0])
    initial = np.eye(len(states),dtype=complex) if full_operator else np.eye(len(states),dtype=complex)[:,start]
    blocks = {(0,0):initial}
    permutations = reader_permutations(control)
    counters = count_table(tape)
    prefixes=[]
    pairmask = np.array([any(row[a]==row[b]==1 for a,b in cell().seams) for row in states])

    def retain(k):
        vecs = {w:(v[:,start] if full_operator else v) for w,v in blocks.items()}
        prefixes.append({'step':k, 'counts':counters[k].tolist(),
            'archive':[{'words':list(w),'weight':float(np.vdot(v,v).real)} for w,v in sorted(vecs.items())],
            'no_pair':float(sum(np.vdot(v[~pairmask],v[~pairmask]).real for v in vecs.values()))})
    retain(0)
    for k,step in enumerate(tape,1):
        blocks=gate_action(step,blocks,kappa,states,permutations)
        retain(k)
    return blocks,prefixes


def couplings(length, kind='engineered'):
    positive_integer(length,'length')
    if kind == 'engineered':
        k=np.arange(length,dtype=float)
        return np.sqrt((k+1)*(length-k))/length
    if kind == 'uniform':
        return np.full(length,0.5)
    if kind == 'broken':
        result=couplings(length)
        result[length//2-1]=0
        return result
    raise ValueError('unknown clock law')


def clock_amplitudes(length,time,kind='engineered'):
    finite(time,'time')
    w,v=eigh_tridiagonal(np.zeros(length+1),couplings(length,kind))
    return v @ (np.exp(-1j*time*w)*v[0].conj())


def clock_law(length,time):
    """Stable complete binomial law, including endpoints, without renormalization."""
    from scipy.stats import binom
    positive_integer(length,'length')
    time=finite(time,'time')
    p=float(np.sin(time/length)**2)
    return binom.pmf(np.arange(length+1),length,p)


def controller_matrix(gates, weights=None):
    """Literal block Hamiltonian, used independently of the history factorization."""
    if not gates:
        raise ValueError('empty gate list')
    d=gates[0].shape[0]
    if d<1:
        raise ValueError('empty data space')
    if weights is None:
        weights=couplings(len(gates))
    weights=np.asarray([finite(w,'edge weight') for w in weights])
    if weights.shape!=(len(gates),) or np.iscomplexobj(weights) or not np.isfinite(weights).all():
        raise ValueError('invalid edge weights')
    blocks=[[None]*(len(gates)+1) for _ in range(len(gates)+1)]
    for k,g in enumerate(gates):
        g=sparse.csr_matrix(g)
        if g.shape!=(d,d) or not np.isfinite(g.data).all():
            raise ValueError('invalid gate')
        residual=g.conj().T@g-sparse.identity(d,format='csr')
        # This Hermitian residual has spectral norm bounded by its maximum
        # absolute row sum. An entrywise test misses coherent rank-one errors.
        if residual.nnz and np.asarray(abs(residual).sum(axis=1)).max()>2e-11:
            raise ValueError('nonunitary gate')
        blocks[k+1][k]=weights[k]*g
        blocks[k][k+1]=weights[k]*g.conj().T
    return sparse.bmat(blocks,format='csr')


def direct_clock_check():
    """Literal full-78 block propagation, with a retained entangled reference."""
    states=source_sector()
    gates=[event_matrix((0,1),8,1/8),event_matrix((0,5),8,1/16),
           np.diag(np.exp(-1j*states[:,0]**2)),event_matrix((1,7),8,1/8)]
    d=len(states)
    H=controller_matrix(gates)
    # All input columns test every source and arbitrary untouched references.
    initial=np.zeros((H.shape[0],d),complex)
    initial[:d]=np.eye(d)
    prefixes=[np.eye(d,dtype=complex)]
    for g in gates:
        prefixes.append(g@prefixes[-1])
    errors=[]
    for fraction in (0.5,1.,2.):
        time=fraction*np.pi*len(gates)/2
        actual=expm_multiply(-1j*time*H,initial,traceA=0)
        amplitudes=clock_amplitudes(len(gates),time)
        expected=np.vstack([a*v for a,v in zip(amplitudes,prefixes)])
        errors.append(float(np.linalg.norm(actual-expected,2)))
    return {'dimension':H.shape[0],'input_columns':d,'operator_errors':errors}
