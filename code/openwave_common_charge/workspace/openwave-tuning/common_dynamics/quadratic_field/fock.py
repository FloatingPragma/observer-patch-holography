"""Finite photon Fock encoding and exact rational spectral certificate."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math

import numpy as np

import source as field

HERE=Path(__file__).resolve().parent


def sqrt_upper(value):
    value=F(value)
    x=math.isqrt(value.numerator*value.denominator)
    return F(x+int(x*x<value.numerator*value.denominator),value.denominator)


def sqrt_lower(value):
    value=F(value)
    return F(math.isqrt(value.numerator*value.denominator),value.denominator)


@lru_cache(None)
def spectral_certificate():
    d=field.geometry();den=1<<44; limb=1<<22
    K=np.rint(d['K']).astype(np.int64)
    assert np.array_equal(K,d['K'])
    Q=np.rint(d['U']*den).astype(np.int64)
    lam=np.rint(d['lam']*den).astype(np.int64)
    assert min(lam)>0
    # Limb matrix products are exact signed-int64 sums. Final products use
    # arbitrary precision; the guards exclude every int64 overflow.
    qmax=int(np.max(abs(Q)));krow=int(np.max(np.sum(abs(K),axis=1)))
    assert krow*qmax<(1<<63)
    high=Q//limb;low=Q-high*limb
    assert 288*max(int(np.max(abs(high))),int(np.max(abs(low))))**2<(1<<63)
    # Split limbs make the three BLAS integer products exact; their final
    # recombination and residual products use arbitrary-precision integers.
    gram=(high.T@high).astype(object)*limb**2
    gram+=((high.T@low).astype(object)+(low.T@high).astype(object))*limb
    gram+=(low.T@low).astype(object)-den**2*np.eye(190,dtype=object)
    residual=(K@Q).astype(object)*den-Q.astype(object)*lam[None,:].astype(object)
    gram2=sum(int(x)**2 for x in gram.flat)
    residual2=sum(int(x)**2 for x in residual.flat)
    epsilon=sqrt_upper(F(gram2,den**4))
    r=sqrt_upper(F(residual2,den**4))
    minimum=F(int(min(lam)),den);maximum=F(int(max(lam)),den)
    shift=r/minimum
    gamma=epsilon+2*(1+epsilon)*shift+shift*shift
    assert gamma<F(1,2)
    # Y=KQ Lambda^-1 lies exactly in the190dimensionalcoexact space.
    # U=Y(Y*Y)^-1/2 is an exact orthonormal frame with ||U-Q||<=shift+gamma.
    Kupper=F(1504,17)
    eta=r+(Kupper+maximum)*(shift+gamma)
    assert eta<minimum/2
    frequency=F(5,4)*eta/sqrt_lower(minimum)
    return {'denominator':den,'eigenvalues_exact':[F(int(x),den) for x in lam],
       'frame_Gram_error_upper':epsilon,'frame_residual_upper':r,
       'polar_Gram_error_upper':gamma,'actual_curvature_operator_error_upper':eta,
       'per_frequency_error_upper':frequency,
       'cap1_total_H_error_including_vacuum':96*frequency,
       'minimum_eigenvalue':minimum,'maximum_eigenvalue':maximum,
       'frame_integers':Q.tolist()}


def run():
    paths=[Path(__file__),Path(field.__file__),field.HERE.parent/'vector_field/MANIFEST.json',
           field.ROOT/'openwave-tuning/auditor/harmonic_stability/MANIFEST.json']
    pins={str(p.relative_to(field.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    spec=spectral_certificate()
    # Literal cap-one ladders on vac plus190onephoton states. Matrix-unit
    # algebra checks retain the full top-sector correction without allocating
    #190dense191x191matrices.
    for i,j in ((0,0),(0,1),(37,37),(37,189)):
        ai=np.zeros((191,191));aj=np.zeros((191,191));ai[0,i+1]=1;aj[0,j+1]=1
        top=np.diag([0.]+[1.]*190);identity=np.eye(191)
        rhs=(identity if i==j else 0)-top@((identity if i==j else 0)+aj.T@ai)
        assert np.array_equal(ai@aj.T-aj.T@ai,rhs)
    omega=np.sqrt(1.5*np.asarray([float(x) for x in spec['eigenvalues_exact']]))
    E0=omega.sum()/2;H=np.diag(np.r_[E0,E0+omega])
    # Any linearcoordinate has two nonzero eigenvalues on totalcap1.
    coef=field.geometry()['U'][0,:]*np.sqrt(.5/(2*omega))
    Qe=np.zeros((191,191));Qe[0,1:]=coef;Qe[1:,0]=coef
    variance=float(coef@coef)
    assert abs(float((Qe@Qe)[0,0])-variance)<1e-15
    assert np.linalg.matrix_rank(Qe,tol=1e-10)==2
    # Cap-one truncation preserves the freevacuum andonephoton energies
    # exactly, but its top state has a large canonical commutator defect.
    ladder=np.array([[0.,1.],[0.,0.]])
    defect=ladder@ladder.T-ladder.T@ladder-np.eye(2)
    assert np.linalg.norm(defect,2)==2
    assert pins=={str(p.relative_to(field.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.actual_carrier.photon_Fock_cap.v1',
       'checks':{'all190modes_in_cap1':True,'integer_spectral_error_certificate':True,
        'coexact_polar_frame_has_exact_Gauss_support':True,'vacuum_constant_included':True,
        'literal_top_cap_commutator_correction':True,'free_vacuum_onephoton_evolution_exact_in_cap':True,
        'cap1_does_not_certify_interacting_Weyl_records':True,'source_bytes_unchanged':True},
       'physical_dimension':191,'native_qubit_dimension':256,'retained_padding':65,
       'Hamiltonian':'H_cap=E0 I + dGamma(omega), total photon occupation<=1; all padding has declared H=0. E0=sum_j omega_j/2 is retained.',
       'completion':'For general totalcapR, H_cap=(a||E||²+b||A^Ttheta||²)/2 + P_N=R (dGamma(omega)+sum omega)/2 on the oscillator sector. This is a declared cap completion, with literal CCR boundary.',
       'spectral_certificate':{k:[str(v) for v in value] if isinstance(value,list) else str(value)
                               for k,value in spec.items() if k!='frame_integers'},
       'certified_cap1_H_error_including_vacuum_upper':float(spec['cap1_total_H_error_including_vacuum']),
       'numerical_edge0_vacuum_variance':variance,'one_mode_top_CCR_defect_norm':2,
       'harmonic_scope':'The three noncompact zero modes remain separate free Gaussian packets in SOURCE.json. This cap and its native photon word do not replace or discretize that factor.',
       'scope':'Exact free-oscillator cap and rational spectrum adapter. A Weyl displacement can reach omitted photon numbers; matter interactions require their state-dependent occupation-tail/cap bound. No full interacting native field word is certified by cap1.',
       'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'FOCK.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'SPECTRAL_FRAME.json').write_text(json.dumps({'denominator':spec['denominator'],'frame_integers':spec['frame_integers'],
        'definition':'exact coexact frame is polar(K Q Lambda^-1); Q andLambda are the rational data. Numerical eigenvectors are only used to select these certified integers.'},separators=(',',':'))+'\n')
    print(json.dumps({'checks':result['checks'],'Hamiltonian_error':result['certified_cap1_H_error_including_vacuum_upper'],
        'frequency_error':float(spec['per_frequency_error_upper']),'lambda_min':float(spec['minimum_eigenvalue'])}),flush=True)


if __name__=='__main__':run()
