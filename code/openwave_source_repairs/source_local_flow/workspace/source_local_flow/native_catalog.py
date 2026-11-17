"""Executed M6 source controls and directed per-pulse comparison bounds."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,math,sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'provider'));sys.path.insert(0,str(HERE/'inputs'))
from model.observer.source_scalar import compiler as c
import intervals as iv

def target(kind,angle):
    sine,cosine=iv.sincos(F(angle));zero=iv.exact(0);one=iv.exact(1)
    out=[[(one if i==j else zero,zero)for j in range(6)]for i in range(6)]
    if kind in ('I','J'):
        for j in range(2 if kind=='I'else 3):out[j][j]=(cosine,sine)
    elif kind=='Z':out[0][0]=(cosine,sine);out[1][1]=(cosine,iv.neg(sine))
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

def sqrt_upper(f):
    scale=10**30
    n=math.isqrt(f.numerator*scale*scale//f.denominator)
    return F(n+1,scale)

def execute(angles):
    ccx=([c.pulse('H',math.pi/2,[2])]+c.compile_diagonal([0.]*7+[math.pi])+
         [c.pulse('H',-math.pi/2,[2])])
    flip=[c.pulse('I',-math.pi/2,[0]),c.pulse('X',math.pi/2,[0])]
    words={'ccx':ccx,'flip':flip,'beam':[c.pulse('Y',-math.pi/4,[0])],
           'unbeam':[c.pulse('Y',math.pi/4,[0])]}
    pulses=[p for w in words.values()for p in w]
    for h in angles:pulses +=[c.pulse('I',-2*h,[0]),c.pulse('J',2*h,[1,0])]
    matrices={};bounds={}
    for pulse in pulses:
        key=pulse['kind']+'/'+float(pulse['angle']).hex()
        if key in matrices:continue
        raw=c.native_pulse(pulse['kind'],pulse['angle']);want=target(pulse['kind'],pulse['angle'])
        sq=sum((iv.scalar_error(raw[i,j].real,want[i][j][0])**2+
                iv.scalar_error(raw[i,j].imag,want[i][j][1])**2
                for i in range(6)for j in range(6)),F())
        matrices[key]=raw;bounds[key]=sqrt_upper(sq)
    # Native raw code products are diagnostic centres. Complete exact source
    # pulse comparison, including proper-code export, uses the bounds above.
    def logical(word,bits):
        out=np.eye(1<<bits,dtype=complex)
        for p in word:
            raw=matrices[p['kind']+'/'+float(p['angle']).hex()]
            targets=p['targets'];d=1<<len(targets);mask=sum(1<<t for t in targets)
            base=np.array([x for x in range(1<<bits)if x&mask==0])
            ix=np.array([base+sum(((k>>(len(targets)-1-j))&1)<<t for j,t in enumerate(targets))for k in range(d)])
            block=out[ix];out[ix]=(raw[:d,:d]@block.reshape(d,-1)).reshape(block.shape)
        return out
    rawccx=logical(ccx,3);ideal=np.eye(8,dtype=complex);ideal[:,[3,7]]=ideal[:,[7,3]]
    assert np.linalg.norm(rawccx-ideal,'fro')<1e-10
    ccz=[c.pulse('H',math.pi/2,[2])]+ccx+[c.pulse('H',-math.pi/2,[2])]
    rawccz=logical(ccz,3);wantccz=np.diag([1]*7+[-1])
    assert np.linalg.norm(rawccz-wantccz,'fro')<1e-10
    b=2*max(bounds.values())+F(1,10**15)
    return {'executed_native_M6_matrices':len(matrices),'raw_frobenius_bounds':{k:str(v)for k,v in bounds.items()},
            'complete_per_pulse_allowance':str(b),'ccx_pulses':len(ccx),
            'ccx_raw_residual':float(np.linalg.norm(rawccx-ideal,'fro')),
            'ccz_raw_residual':float(np.linalg.norm(rawccz-wantccz,'fro')),
            'source_control_semantics':'Each pulse executes exp of exact source-generator coefficients on full M6. Proper-code routes are composed with all processor complements retained by prefix bounds.'},matrices

def complete(pulses,b):
    x=pulses*b;assert x<F(1,10)
    e=x/(1-x);failure=pulses*b*b/(1-x)**2+e*e
    return (e*(2+e)+failure)/2
