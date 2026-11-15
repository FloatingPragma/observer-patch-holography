"""Literal complete address analyzer, with no dense 10800-square matrix."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
from scipy import linalg as la

from model.lepton_cycle import instrument
from model.native_flavor_orbit import intervals as iv

PHASE=F(2,9)
FLAGS=3600
FAMILIES=3


@lru_cache(maxsize=1)
def data():
    records,_=instrument.source_data()
    plus=np.array([r==2 for r in records]);minus=np.array([r==3 for r in records]);blank=~(plus|minus)
    masks=(blank[:,None]&(~blank)[None,:],plus[:,None]&blank[None,:],minus[:,None]&blank[None,:])
    counts=[int(m.sum()) for m in masks]
    if counts!=[500,250,250]:raise ArithmeticError('actual native address partition changed')
    phi=(sum(m/np.sqrt(n) for m,n in zip(masks,counts,strict=True))/np.sqrt(3)).reshape(FLAGS)
    s=instrument.preparation().reshape(FLAGS)
    first=np.zeros(FLAGS);first[0]=1
    qphi=(first-phi)/la.norm(first-phi)
    first60=np.zeros(60);first60[0]=1
    q60=(first60-np.ones(60)/np.sqrt(60));q60/=la.norm(q60)
    angles=float(PHASE)+2*np.pi*np.arange(3)/3
    phase=np.ones((60,60,3),complex)
    phase[plus,:,:]=np.exp(1j*angles)[None,None,:]
    phase[minus,:,:]=np.exp(-1j*angles)[None,None,:]
    return {'phi':phi,'S':s,'qphi':qphi,'q60':q60,'phase':phase,'counts':counts}


def _source_householder(v,q):
    value=v.reshape(60,60,3,-1)
    value=value-2*np.einsum('i,jkl->ijkl',q,np.einsum('i,ijkl->jkl',q,value))
    value=value-2*np.einsum('j,ikl->ijkl',q,np.einsum('j,ijkl->ikl',q,value))
    return value.reshape(FLAGS,3,-1)


def _reflection(v):
    value=v.reshape(60,60,3,-1)
    value=2*value.mean(axis=0,keepdims=True)-value
    value=2*value.mean(axis=1,keepdims=True)-value
    return value.reshape(FLAGS,3,-1)


def apply(v,adjoint=False):
    """U=H_phi D R H_S on every flag/family input column."""
    value=np.asarray(v,complex)
    if value.ndim not in (2,3) or value.shape[:2]!=(FLAGS,3) or not np.isfinite(value).all():
        raise ValueError('complete 3600-by-three flag/family input required')
    shape=value.shape;value=value.reshape(FLAGS,3,-1);d=data()
    def hp(x):return x-2*d['qphi'][:,None,None]*np.einsum('x,xkn->kn',d['qphi'],x)[None,:,:]
    phase=d['phase'].reshape(FLAGS,3,1)
    if adjoint:
        value=_source_householder(_reflection(phase.conj()*hp(value)),d['q60'])
    else:
        value=hp(phase*_reflection(_source_householder(value,d['q60'])))
    return value.reshape(shape)


def rational_effect_intervals():
    theta=iv.I(PHASE)
    sine=iv.sine(theta)
    # cos(theta)=1-2 sin(theta/2)^2 uses only the certified sine grammar.
    cosine=1-2*iv.sine(theta/2)**2
    root3=iv.sqrt_rational(3)
    cs=(cosine,-cosine/2-root3*sine/2,-cosine/2+root3*sine/2)
    a=iv.sqrt_rational(F(5,108));b=iv.sqrt_rational(F(5,216))
    lambdas=[a+2*b*c for c in cs]
    if any(x.lo<=0 for x in lambdas):raise ArithmeticError('positive-root chamber not certified')
    lambdas=[iv.rounded_out(x,38) for x in lambdas]
    effects=[iv.rounded_out(x*x,38) for x in lambdas]
    if effects[1].lo<=F(1,20000) or not effects[1].hi<effects[2].lo<effects[2].hi<effects[0].lo:
        raise ArithmeticError('three ordered positive source coefficients not separated')
    return {'lambda_intervals_k0_k1_k2':[x.strings() for x in lambdas],
            'effect_intervals_k0_k1_k2':[x.strings() for x in effects],
            'charged_mass_floor_exact':'1/20000'}


@lru_cache(maxsize=1)
def certificate():
    d=data();basis=np.zeros((FLAGS,3,3),complex);basis[0]=np.eye(3)
    columns=apply(basis);row=apply(basis,adjoint=True)
    k=columns[0];y=k.conj().T@k
    _,e,amplitudes=instrument.analyzer_data(float(PHASE))
    if la.norm(k-np.diag(amplitudes))>2e-12 or la.norm(y-np.diag(e))>2e-12:
        raise ArithmeticError('literal native analyzer compression changed')
    # W=U† Q U is tested on all complete input columns in bounded memory.
    maximum=0.;projection=0.
    total=FLAGS*3
    for start in range(0,total,180):
        stop=min(total,start+180);b=np.zeros((total,stop-start),complex)
        b[np.arange(start,stop),np.arange(stop-start)]=1
        b=b.reshape(FLAGS,3,-1)
        out=apply(b);restored=apply(out,adjoint=True)
        maximum=max(maximum,float(np.max(abs(restored-b))))
        projected=out.copy();projected[1:]=0
        direct=apply(projected,adjoint=True)
        # The batch consists of literal basis columns, so this is exactly
        # row^dagger times the batch without a dense identity contraction.
        inner=row.reshape(total,3)[start:stop].conj().T
        independent=np.einsum('xki,ib->xkb',row,inner,optimize=True)
        projection=max(projection,float(np.max(abs(direct-independent))))
    if maximum>2e-12 or projection>2e-12:raise ArithmeticError('complete factorized analyzer failed')
    return {'address_dimension':FLAGS,'family_dimension':3,'complete_unitary_dimension':total,
            'all_inverse_columns_checked':total,'all_projector_columns_checked':total,
            'maximum_inverse_column_error':maximum,'maximum_projector_column_error':projection,
            'native_mask_counts':d['counts'],'phase_exact':str(PHASE),
            'phase_role':'inherited fixed coherent phase input; no passive-clock selection proved',
            'exact_factors':'U=H_phi D_(theta+2pi k/3) R H_S; all factors unitary on every input',
            'exact_blank_block':'K00=diag(sqrt(5/108)+2 sqrt(5/216) cos(theta+2pi k/3))',
            'exact_mean':'Y=K00†K00=diag(e_k); Tr Y=5/18',
            'effect_values_k0_k1_k2':e.tolist(),'lambda_values_k0_k1_k2':amplitudes.real.tolist(),
            'actual_native_mass_instrument_compression_error':float(la.norm(y-np.diag(e))),
            'wrong_equal_leg_frame_full_projector_difference':float(la.norm(row-row.conj())),
            'entire_unitary_lifted_before_compression':True,
            **rational_effect_intervals()}
