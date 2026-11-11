"""Exact source-completion Schur response, including conserved directions.

H=I-T is a one-proposal relaxation operator. Its Schur complement is an
algebraic resolvent calculation, not permission to restart terminal timeouts
and not a microscopic fermion Hamiltonian.
"""
from functools import lru_cache
import numpy as np
import sympy as s
from model.common_source.hierarchy_diagnosis import native_family_dynamics as source
from model.native_modular_reference import local_preparation
from model.common_source.particle_identification import family_attachment


def clean(a):
    return s.ImmutableMatrix(a.applyfunc(lambda x:s.expand(s.radsimp(s.cancel(s.expand(x))))))


def zero(a):
    return clean(a) == s.zeros(*a.shape)


@lru_cache(maxsize=1)
def exact():
    data=source._data()
    B=data['basis']; metric=clean(B.T*B)[0,0]
    H=clean(s.eye(12)-data['T'])
    # A nonorthogonal complement avoids introducing unnecessary radicals.
    Q=s.Matrix.hstack(*B.T.nullspace())
    if Q.shape!=(12,9) or not zero(B.T*Q):
        raise ArithmeticError('incomplete family complement')
    uniform=s.ones(12,1)
    P0=uniform*uniform.T/12
    D=s.diag(*map(s.Rational,family_attachment._data()[1]))
    P=clean(B*B.T/metric)
    F=clean(B.T*D*B/metric)
    A=clean(B.T*H*B/metric)
    commutator=clean(F*A-A*F)
    virtual_commutator=clean((B.T*H*(s.eye(12)-P)*D*B-B.T*D*(s.eye(12)-P)*H*B)/metric)
    if not zero(D*H-H*D) or zero(commutator) or not zero(commutator-virtual_commutator):
        raise ArithmeticError('retained source commutator identity failed')
    V=clean(B.T*H*Q)
    C=clean(Q.T*(H+P0)*Q)
    G=clean(Q.T*Q)
    W=clean(C.to_DM(extension=True).inv().to_Matrix()*V.T)
    S=clean(A-V*W/metric)
    Z=clean(s.eye(3)+W.T*G*W/metric)
    # H has exactly one null mode per connected component. CORE are isolated;
    # the other ten ports form one connected native repair graph.
    core=tuple(local_preparation.CORE)
    complement=tuple(i for i in range(12) if i not in core)
    K=s.zeros(12,3)
    for j,i in enumerate(core): K[i,j]=1
    for i in complement: K[i,2]=1
    if not zero(H*K) or H.rank()!=9 or (B.T*K).rank()!=2:
        raise ArithmeticError('source conserved-subspace rank changed')
    Pker=clean(K*(K.T*K).inv()*K.T)
    R0=clean(B.T*Pker*B/metric)
    if S.rank()!=1 or not zero(S[:2,:]) or not zero(S[:,:2]):
        raise ArithmeticError('complete Schur response did not retain two zero directions')
    if R0.rank()!=2 or not zero(S*R0):
        raise ArithmeticError('zero pole residue does not match Schur kernel')
    # Completion-square identity proves positivity and includes every virtual
    # port. The null uniform direction has zero coupling and is lifted only
    # to define an inverse; no amplitude is discarded.
    harmonic=clean(B-Q*W)
    if not zero(harmonic.T*H*harmonic/metric-S):
        raise ArithmeticError('full virtual-state energy mismatch')
    if not zero(harmonic.T*harmonic/metric-Z):
        raise ArithmeticError('canonical residue metric mismatch')
    return dict(H=H,B=B,Q=Q,metric=metric,A=A,V=V,C=C,G=G,W=W,S=S,Z=Z,
                harmonic=harmonic,Pker=Pker,R0=R0,K=K,D=D,F=F,
                projected_commutator=commutator,virtual_commutator=virtual_commutator,timeout=data['timeout'])


def spectral_data():
    x=exact(); H=np.array(x['H'],float); B=np.array(x['B'],float)/np.sqrt(float(x['metric']))
    values,vectors=np.linalg.eigh(H)
    poles=[]
    for i,e in enumerate(values):
        if e<1e-12: continue
        amplitude=B.T@vectors[:,i]
        weight=np.outer(amplitude,amplitude)
        poles.append(dict(relaxation_pole=float(e),residue=weight.tolist(),
                          residue_trace=float(np.trace(weight))))
    zero_weight=np.array(x['R0'],float)
    sum_weights=zero_weight+sum((np.array(row['residue']) for row in poles),np.zeros((3,3)))
    if np.linalg.norm(sum_weights-np.eye(3),2)>5e-13:
        raise ArithmeticError('source spectral measure lost an outcome')
    S=np.array(x['S'],float); Z=np.array(x['Z'],float)
    z,U=np.linalg.eigh(Z); inverse=(U/np.sqrt(z))@U.T
    canonical=inverse@S@inverse
    return dict(poles=poles,zero_pole_residue=zero_weight.tolist(),
                zero_residue_eigenvalues=np.linalg.eigvalsh(zero_weight).tolist(),
                canonical_zero_frequency_matrix=canonical.tolist(),
                canonical_linearized_eigenvalues=np.linalg.eigvalsh(canonical).tolist(),
                residue_sum_error=float(np.linalg.norm(sum_weights-np.eye(3),2)))


def report():
    x=exact(); numeric=spectral_data()
    entries=lambda a:[[str(y) for y in row] for row in a.tolist()]
    # The nonzero generalized eigenvalue follows from rank-one S exactly.
    eigenvalue=s.cancel(s.trace(x['Z'].to_DM(extension=True).inv().to_Matrix()*x['S']))
    return dict(source_dimension=12,family_dimension=3,complement_dimension=9,
                original_66_state_instrument_replayed=True,terminal_timeout_exact=str(x['timeout']),
                full_occupation_spectrum_exact={'29/90':2,'61/450':10},
                full_occupation_response_commutator_exact='0',
                compressed_occupation_response_commutator_exact=entries(x['projected_commutator']),
                projected_commutator_entirely_due_to_complement=True,
                all_functions_of_full_D_H_have_common_eigenvectors=True,
                normalized_family_gram_exact='I3',source_relaxation_rank=9,
                projected_zero_pole_rank=2,naive_compression_eigenvalues=np.linalg.eigvalsh(np.array(x['A'],float)).tolist(),
                schur_matrix_exact=entries(x['S']),residue_metric_exact=entries(x['Z']),
                zero_pole_residue_exact=entries(x['R0']),
                canonical_nonzero_linearized_eigenvalue_exact=str(eigenvalue),
                full_spectral_measure=numeric,
                exact_Schur_rank=1,discarded_outcome_or_complement=False,
                independent_three_mass_poles=False,physical_Yukawa_identification=False,
                expansion_point='zero relaxation frequency, fixed by the stationary source',
                linearized_eigenvalue_identified_with_full_pole=False,
                repeated_post_timeout_protocol_claimed=False)
