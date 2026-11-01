"""Faithful same-action reference and complete source entropy bookkeeping.

All field operators act on the full16-state CAR space of the joint source
fixture. This does not identify a geometric cap modular charge or a metric.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from scipy.linalg import expm

from ..defect_field_bridge import joint
from ..field_thermal.thermal import car_restriction, restrict
from ..modular_flow.flow import preparation_resources
from ..observer.source_scalar import compiler
from ..reference_observation.observation import ReferenceObservation
from ..reference_recovery.recovery import execute_complete
from . import controls

BETA = F(1)  # one declared field reference inverse temperature, GeV^-1
ADMISSION = 2e-8


def positive_beta(beta):
    if isinstance(beta,(bool,str)) or not isinstance(beta,(int,float,F)):
        raise ValueError('finite positive inverse temperature required')
    try: finite=math.isfinite(float(beta))
    except (ValueError,OverflowError): finite=False
    if not finite or not 0<float(beta)<=8:
        raise ValueError('bounded inverse temperature in (0,8] required')
    return F(beta)


def faithful_log(density):
    value=np.asarray(density)
    if (value.ndim!=2 or value.shape[0]!=value.shape[1] or value.shape[0]>64
            or value.dtype.kind not in 'iufc' or not np.isfinite(value).all()
            or np.linalg.norm(value-value.conj().T,2)>1e-11):
        raise ValueError('bounded finite Hermitian matrix required')
    eigenvalues,vectors=np.linalg.eigh(value)
    if eigenvalues[0]<=0:
        raise ValueError('faithfulness is required; no support is clipped')
    return (vectors*np.log(eigenvalues))@vectors.conj().T


def gibbs(hamiltonian,beta=BETA):
    beta=positive_beta(beta)
    h=np.asarray(hamiltonian,complex)
    if h.shape!=(16,16) or not np.isfinite(h).all() or np.linalg.norm(h-h.conj().T,2)>1e-11:
        raise ValueError('full16-dimensional Hermitian same-action Hamiltonian required')
    energies,vectors=np.linalg.eigh(h)
    weights=np.exp(-float(beta)*(energies-energies[0])); weights/=weights.sum()
    if weights.min()<1e-7:
        raise ValueError('reference leaves the declared faithful numerical margin')
    rho=(vectors*weights)@vectors.conj().T
    root=(vectors*np.sqrt(weights))@vectors.conj().T
    logz=math.log(float(np.sum(np.exp(-float(beta)*(energies-energies[0])))))-float(beta)*energies[0]
    K=-faithful_log(rho)
    if np.linalg.norm(K-float(beta)*h-logz*np.eye(16),2)>ADMISSION:
        raise ArithmeticError('same-action Gibbs modular identity failed')
    return {'H':h,'rho':rho,'root':root,'K':K,'beta':beta,'log_partition':float(logz),
            'minimum_eigenvalue':float(weights.min()),'energies':energies}


@lru_cache(maxsize=1)
def same_action():
    from ..defect_field_energy.energy import small_energy_oracle
    return small_energy_oracle(phase=0)


@lru_cache(maxsize=1)
def reference():
    return gibbs(same_action()['H'])


@lru_cache(maxsize=1)
def prepare_reference():
    data=reference(); target=data['root'].T.reshape(-1)
    word=controls.load()['preparation']
    budget=preparation_resources(word,8)
    initial=np.zeros((256,1),complex);initial[0,0]=1
    run=execute_complete(word,8,initial)
    error=math.sqrt(np.linalg.norm(run.success[:,0]-target)**2+run.failure_probability)
    if error>ADMISSION:raise ArithmeticError('complete Gibbs preparation exceeds admission')
    columns=run.success[:,0].reshape(16,16).T.copy()
    return {'columns':columns,'run':run,'source_word':word,
            'resources':compiler.resources(word),'preparation_certificate':budget,
            'complete_error':float(error),'purifier_dimension':16}


def observation_tower(rho):
    outer=ReferenceObservation(rho,car_restriction(4,(0,1)))
    inner=ReferenceObservation(outer.target_reference,car_restriction(2,(0,)))
    direct=ReferenceObservation(rho,car_restriction(4,(0,)))
    return outer,inner,direct


def regional_certificate():
    data=reference();rho=data['rho'];H=data['H'];outer,inner,direct=observation_tower(rho)
    regional_min=float(np.linalg.eigvalsh(outer.target_reference)[0])
    one_mode_min=float(np.linalg.eigvalsh(inner.target_reference)[0])
    if (regional_min+1e-12<4*data['minimum_eigenvalue']
            or one_mode_min+1e-12<2*regional_min):
        raise ArithmeticError('complete partial trace lost its strict positivity bound')
    # Independent superoperator equality on every complex matrix unit.
    residual=0.
    for a in range(16):
        for b in range(16):
            unit=np.zeros((16,16),complex);unit[a,b]=1
            residual=max(residual,float(np.linalg.norm(inner.observe(outer.observe(unit))-direct.observe(unit),2)))
    if residual>ADMISSION:raise ArithmeticError('nested same-reference observation failed')
    kA=-faithful_log(outer.target_reference)
    bare=restrict(same_action()['h_sites'][0],car_restriction(4,(0,1)))/4
    center=lambda x:x-np.trace(x)*np.eye(len(x))/len(x)
    mismatch=float(np.linalg.norm(center(kA)-float(BETA)*center(bare),'fro'))
    # The global Hamiltonian has a nonzero cross term in either local chart.
    HA=restrict(H,car_restriction(4,(0,1)))/4
    HB=restrict(H,car_restriction(4,(2,3)))/4
    cross=H-np.kron(np.eye(4),HA)-np.kron(HB,np.eye(4))+np.trace(H)*np.eye(16)/16
    cross_squared=float(np.vdot(cross,cross).real)
    if abs(cross_squared-4)>ADMISSION:
        raise ArithmeticError('actual action differs from the exact hopping cross-term certificate')
    embedded=np.kron(np.eye(4),np.diag([1.,-1.,1.,-1.]))
    retraction=float(np.linalg.norm(outer.observe(embedded)-np.diag([1.,-1.,1.,-1.]),2))
    s=F(1,3);globalU=expm(1j*float(s)*faithful_log(rho));localU=expm(1j*float(s)*faithful_log(outer.target_reference))
    covariance=float(np.linalg.norm(outer.observe(globalU@embedded@globalU.conj().T)-localU@outer.observe(embedded)@localU.conj().T,2))
    return {'outer':outer,'inner':inner,'direct':direct,
            'nested_all_matrix_units_error':residual,'regional_K':kA,'bare_regional_H':bare,
            'regional_minimum_eigenvalue':regional_min,'one_mode_minimum_eigenvalue':one_mode_min,
            'regional_faithfulness_bound':'rho>=lambda I implies Tr_B(rho)>=dim(B)*lambda I; no regional support removal',
            'same_beta_bare_energy_mismatch':mismatch,'nonadditive_H_HS_squared':cross_squared,
            'retraction_control_error':retraction,'modular_covariance_control_error':covariance}


@lru_cache(maxsize=1)
def source_recovery():
    """Reuse the existing bounded native Petz compiler, all four inputs."""
    observation=ReferenceObservation(reference()['rho'],car_restriction(4,(0,1)))
    target=np.vstack([k.conj().T for k in observation.kraus])
    word=controls.load()['recovery']
    initial=np.zeros((64,4),complex);initial[:4]=np.eye(4)
    run=execute_complete(word,6,initial)
    error=math.sqrt(np.linalg.norm(run.success-target,2)**2+np.linalg.norm(run.failure_gram,2))
    if error>ADMISSION:raise ArithmeticError('complete frozen Petz source exceeds its target admission')
    return {'observation':observation,'target':target,'run':run,'word':word,
            'resources':compiler.resources(word),'complete_error':float(error)}


def entropy_from_amplitude(columns):
    """Subnormalized entropy; every singular value, including tiny ones, remains."""
    values=np.linalg.svd(np.asarray(columns),compute_uv=False)**2
    if not np.isfinite(values).all():raise ArithmeticError('nonfinite source density spectrum')
    return float(-sum(value*math.log(value) for value in values if value>0))


def entropy(density):
    return float(-np.trace(density@faithful_log(density)).real)


def relative_entropy(density,reference_density):
    """Full quantum logarithms; no replacement by basis populations."""
    return float(np.trace(density@(faithful_log(density)-faithful_log(reference_density))).real)


def entropy_bregman_remainder(reference_density,direction,epsilon):
    """Correct even when the retained numerical direction has trace drift."""
    changed=reference_density+epsilon*direction
    return relative_entropy(changed,reference_density)-epsilon*float(np.trace(direction).real)


def exact_boundary_obstruction():
    """Integer CAR computation of ||H_cross||_HS^2=4 in the bare chart.

    The actual change to the QP chart is local unitary on each site, hence
    preserves this norm and the subspace of additive regional Hamiltonians.
    """
    a=[]
    for mode in range(4):
        m=np.zeros((16,16),dtype=np.int64)
        for mask in range(16):
            if mask>>mode&1:
                m[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
        a.append(m)
    # H_cross=(i/2)*numerator. All arithmetic here is exact integer arithmetic.
    numerator=a[0].T@a[2]-a[2].T@a[0]-a[1].T@a[3]+a[3].T@a[1]
    left=np.trace(numerator.reshape(4,4,4,4),axis1=0,axis2=2)
    right=np.trace(numerator.reshape(4,4,4,4),axis1=1,axis2=3)
    if np.any(left) or np.any(right):raise ArithmeticError('hopping gained a local component')
    squared=F(int(np.sum(numerator*numerator)),4)
    if squared!=4:raise ArithmeticError('exact cross-boundary CAR norm changed')
    return {'cross_term_Hilbert_Schmidt_norm_squared_exact':str(squared),
            'both_partial_traces_exactly_zero':True,
            'local_chart_preserves_norm_and_additivity':True,
            'scope':'For beta>0 the same Gibbs density is not a regional tensor product. Thus no UCP retraction onto the full site algebra preserves this reference; the non-idempotent Petz observation remains valid.'}


@lru_cache(maxsize=1)
def field_gates():
    return tuple(joint.field_gate(j)[0] for j in (-1,0,1))


def history_entropy(law,policy,*,events=2):
    if law not in ('C','E') or policy not in ('coherent','flagged') or events not in (1,2):
        raise ValueError('one/two complete C/E histories and explicit collision policy required')
    data=reference();prepared=prepare_reference();H,rho=data['H'],data['rho']
    K,_,resources=joint.repair(law,policy)
    U=field_gates()
    gates=[U[joint.current(raw)+1]@K[coin*16+raw] for coin in range(2) for raw in range(16)]
    operators=gates if events==1 else [b@a for a in gates for b in gates]
    outputs=np.stack([a@prepared['columns'] for a in operators])
    probabilities=np.sum(abs(outputs)**2,axis=(1,2))
    unobserved=np.einsum('hij,hkj->ik',outputs,outputs.conj())
    effect=sum((a.conj().T@a for a in operators),np.zeros((16,16),complex))
    reference_marginal=np.einsum('hji,hjk->ik',outputs.conj(),outputs)
    if (np.linalg.norm(effect-np.eye(16),2)>ADMISSION
            or np.linalg.norm(reference_marginal-prepared['columns'].conj().T@prepared['columns'],2)>ADMISSION):
        raise ArithmeticError('history lost complete effects or its untouched purifier')
    total=float(probabilities.sum());logrho=faithful_log(rho)
    sref=entropy(rho);sout=entropy(unobserved);cq=sum(map(entropy_from_amplitude,outputs))
    classical=float(-sum(p*math.log(p) for p in probabilities if p>0))
    weighted_relative=0.
    for columns,p in zip(outputs,probabilities,strict=True):
        if p==0:continue
        block=columns@columns.conj().T
        weighted_relative+=-entropy_from_amplitude(columns)-float(np.trace(block@logrho).real)-p*math.log(p)
    delta_energy=float(np.trace((unobserved-rho)@H).real)
    relative=relative_entropy(unobserved,rho)
    trace_correction=(total-1)*data['log_partition']
    unread_error=sout-sref-float(BETA)*delta_energy+relative-trace_correction
    cq_error=cq-sref-classical-float(BETA)*delta_energy+weighted_relative-trace_correction
    if max(abs(unread_error),abs(cq_error))>ADMISSION:
        raise ArithmeticError('quantum/CQ entropy balance failed')
    # A true normalized variation is defined by the exact ideal map. For the
    # source diagnostic its tiny mass residual is retained in the formula.
    delta=unobserved-rho;epsilon=F(1,64);changed=rho+float(epsilon)*delta
    derivative=float(-np.trace(delta@(logrho+np.eye(16))).real)
    firstlaw_error=derivative-float(BETA)*delta_energy-(data['log_partition']-1)*(total-1)
    raw_relative=relative_entropy(changed,rho)
    # For subnormalized numerical outputs the entropy Bregman remainder
    # subtracts its linear trace drift. Never normalize the output to hide it.
    remainder=entropy_bregman_remainder(rho,delta,float(epsilon))
    bound=float(epsilon)**2*float(np.vdot(delta,delta).real)/(2*(1-float(epsilon))*data['minimum_eigenvalue'])
    return {'outputs':outputs,'probabilities':probabilities,'unread_density':unobserved,
            'effect':effect,'purifier_marginal':reference_marginal,
            'law':law,'policy':policy,'events':events,'complete_records':len(outputs),
            'total_probability':total,'unread_entropy_change':sout-sref,'cq_entropy_change':cq-sref,
            'classical_record_entropy':classical,'delta_energy_GeV':delta_energy,
            'quantum_relative_entropy':relative,'weighted_conditional_quantum_relative_entropy':weighted_relative,
            'numerical_trace_correction':trace_correction,'unread_identity_error':unread_error,'cq_identity_error':cq_error,
            'noncommuting_variation_HS':float(np.linalg.norm(delta@rho-rho@delta,'fro')),
            'firstlaw_derivative':derivative,'firstlaw_error':firstlaw_error,
            'epsilon':str(epsilon),'relative_entropy_remainder':remainder,
            'raw_relative_entropy_at_epsilon':raw_relative,
            'remainder_linear_trace_correction':float(epsilon)*(total-1),
            'analytic_remainder_upper_numerical':bound,
            'reference_stationarity_error':float(np.linalg.norm(delta,2)),
            'unitality_error':float(np.linalg.norm(sum(a@a.conj().T for a in operators)-np.eye(16),2)),
            'resource_per_repair':resources}
