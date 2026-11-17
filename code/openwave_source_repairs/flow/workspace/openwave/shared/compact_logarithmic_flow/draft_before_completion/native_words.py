"""Full M6 source words with directed primitive and preparation certificates."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
PROVIDER=HERE.parent/'localized_rest_join'
sys.path.insert(0,str(PROVIDER/'standalone'/'reverse-engineering-reality'))
sys.path.insert(0,str(HERE/'inputs'))
from model.observer.source_scalar import compiler as c
import intervals as iv

RAW_GATE=F(1,10**11)
PULSE_GATE=F(2,10**11)
SPLIT_GATE=F(2,10**12)
ANGLE_DENOMINATOR=1<<44
GATES={}
RAW_BOUNDS={}
COUNTS={'native_matrix_executions':0,'pulse_applications':0,'identity_applications':0}

def complex_interval(real,imag):return real,imag

def primitive_target(kind,angle):
    sine,cosine=iv.sincos(F(angle));zero=iv.exact(0);one=iv.exact(1)
    out=[[(one if i==j else zero,zero) for j in range(6)] for i in range(6)]
    if kind=='J':
        for j in range(3):out[j][j]=(cosine,sine)
    elif kind=='I':
        for j in range(2):out[j][j]=(cosine,sine)
    elif kind=='Z':
        out[0][0]=(cosine,sine);out[1][1]=(cosine,iv.neg(sine))
    elif kind=='X':
        out[0][0]=out[1][1]=(cosine,zero);out[0][1]=out[1][0]=(zero,sine)
    elif kind=='Y':
        out[0][0]=out[1][1]=(cosine,zero)
        out[0][1]=(sine,zero);out[1][0]=(iv.neg(sine),zero)
    elif kind=='H':
        root=math.isqrt(iv.Q*iv.Q//2);part=iv.mul(sine,(root,root+1))
        out[0][0]=(cosine,part);out[1][1]=(cosine,iv.neg(part))
        out[0][1]=out[1][0]=(zero,part)
    else:raise ValueError(kind)
    return out

def native_gate(kind,angle):
    key=(kind,float(angle).hex())
    if key not in GATES:
        full=c.native_pulse(kind,angle);wanted=primitive_target(kind,angle)
        squared=sum((iv.scalar_error(full[i,j].real,wanted[i][j][0])**2+
                     iv.scalar_error(full[i,j].imag,wanted[i][j][1])**2
                     for i in range(6) for j in range(6)),F())
        assert squared<=RAW_GATE**2,(kind,angle,str(squared))
        GATES[key]=full;RAW_BOUNDS[key]=squared;COUNTS['native_matrix_executions']+=1
    return GATES[key]

def sqrt_ratio(numerator,denominator):
    assert numerator>=0 and denominator>0
    # A directed120-bit interval, with no floating square-root premise.
    a=math.isqrt((numerator*iv.Q*iv.Q)//denominator)
    return a,a+1

def preparation(integer_vector):
    values=[int(x) for x in integer_vector];n=len(values);bits=n.bit_length()-1
    assert n==1<<bits and any(values)
    word=[];worst=F()
    for target in reversed(range(bits)):
        controls=bits-target-1;angles=[]
        for prefix in range(1<<controls):
            block=values[prefix*(1<<(target+1)):(prefix+1)*(1<<(target+1))]
            half=1<<target;left=sum(x*x for x in block[:half]);right=sum(x*x for x in block[half:])
            if left+right==0:angle=0
            else:
                x=float(block[0]) if target==0 else math.sqrt(left)
                y=float(block[1]) if target==0 else math.sqrt(right)
                angle=round(-math.atan2(y,x)*ANGLE_DENOMINATOR)
                si,co=iv.sincos(F(angle,ANGLE_DENOMINATOR))
                tx=sqrt_ratio(left,left+right);ty=sqrt_ratio(right,left+right)
                if target==0 and block[0]<0:tx=iv.neg(tx)
                if target==0 and block[1]<0:ty=iv.neg(ty)
                sy=iv.neg(si)
                error=max(abs(F(co[i]-tx[j],iv.Q)) for i in (0,1) for j in (0,1))
                error=max(error,max(abs(F(sy[i]-ty[j],iv.Q)) for i in (0,1) for j in (0,1)))
                assert error<=SPLIT_GATE
                worst=max(worst,error)
            angles.append(angle)
        for mask in range(1<<controls):
            numerator=sum(a*(-1)**((prefix&mask).bit_count()) for prefix,a in enumerate(angles))
            exact=F(numerator,ANGLE_DENOMINATOR*(1<<controls));angle=float(exact)
            assert F(angle)==exact
            label=['I']*bits;label[target]='Y'
            for control in range(controls):
                if mask>>control&1:label[target+1+control]='Z'
            word.extend(c.pauli_rotation(label,angle))
    assert len(word)==3*((bits-1)*(1<<bits)+1)
    return word,2*bits*SPLIT_GATE,worst

def address_phase(angle,bits=8):
    # Exact Walsh coefficients for the phase at all-zero address.
    coeff=angle/(1<<bits);assert F(coeff)*(1<<bits)==F(angle)
    return [p for mask in range(1<<bits)
            for p in c.pauli_rotation(['Z' if mask>>j&1 else 'I' for j in range(bits)],coeff)]

@lru_cache(maxsize=None)
def indices(targets,bits):
    dimension=1<<bits;mask=sum(1<<t for t in targets)
    bases=np.array([i for i in range(dimension) if i&mask==0],dtype=np.int64)
    return np.array([bases+sum(((k>>(len(targets)-1-j))&1)<<t for j,t in enumerate(targets))
                     for k in range(1<<len(targets))])

def apply_word(word,states,bits=8):
    """Actual raw success products; first-failure reference Gram retained.

    Each pair of columns is a separate f/g-reference experiment. Complete
    fine failure maps remain reproducible from the emitted word/prefix.
    """
    out=np.array(states,dtype=complex,copy=True)
    if out.ndim==1:out=out[:,None]
    failure=np.zeros((out.shape[1],out.shape[1]),complex)
    for pulse in word:
        full=native_gate(pulse['kind'],pulse['angle']);targets=tuple(pulse['targets']);d=1<<len(targets)
        COUNTS['pulse_applications']+=1
        if np.array_equal(full,np.eye(6)):
            COUNTS['identity_applications']+=1
            continue
        ix=indices(targets,bits);block=out[ix];flat=block.reshape(d,-1)
        comp=full[d:,:d]
        if np.any(comp):
            bad=(comp@flat).reshape(-1,out.shape[1])
            failure+=bad.conj().T@bad
        out[ix]=(full[:d,:d]@flat).reshape(block.shape)
    return out,failure

def digest_word(word):
    h=hashlib.sha256()
    for p in word:
        h.update((p['kind']+' '+float(p['angle']).hex()+' '+','.join(map(str,p['targets']))+'\n').encode())
    return h.hexdigest()

def complete_bound(pulses,target_error=F()):
    x=pulses*PULSE_GATE;assert x<1
    e=x/(1-x)+target_error
    f=pulses*PULSE_GATE**2/(1-x)**2
    return {'operator':e,'processor_failure':f,'half_diamond':(e*(2+e)+f)/2,
            'gain':1/(1-x)**2+f}

def gate_archive():
    keys=sorted(GATES)
    return keys,np.asarray([GATES[k] for k in keys]),[str(RAW_BOUNDS[k]) for k in keys]
