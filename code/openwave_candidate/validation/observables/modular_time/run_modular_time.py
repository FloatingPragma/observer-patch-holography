"""Finite observer modular flow from an explicitly conditioned native output.

The modular parameter is dimensionless and is not identified with native clock
or subjective phenomenology. No Gibbs preparation or depolarizing floor is used.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import sys
import numpy as np
from scipy.linalg import expm,logm

BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.uniform_native import UniformFiberReduction
from oph_openwave.geometry import flat_ring


def matrix_data(value):
    return {'real':np.asarray(value).real.tolist(),'imag':np.asarray(value).imag.tolist()}


def faithful_spectrum(rho,tolerance=1e-12):
    rho=np.asarray(rho,complex)
    if rho.ndim!=2 or rho.shape[0]!=rho.shape[1] or not np.all(np.isfinite(rho)) or not np.allclose(rho,rho.conj().T) or not np.isclose(np.trace(rho),1):
        raise ValueError('finite normalized Hermitian observer state required')
    values,vectors=np.linalg.eigh(rho)
    if min(values)<=tolerance:raise ValueError('state is not numerically faithful on the declared algebra; no floor or pseudoinverse inserted')
    return values,vectors


def modular_flow(rho,observable,s):
    if not math.isfinite(s):raise ValueError('finite real modular parameter required')
    rho=np.asarray(rho,complex)
    values,vectors=faithful_spectrum(rho)
    observable=np.asarray(observable,complex)
    if observable.shape!=rho.shape or not np.all(np.isfinite(observable)):raise ValueError('observable dimension/value invalid')
    unitary=(vectors*np.exp(1j*s*np.log(values)))@vectors.conj().T
    return unitary@observable@unitary.conj().T


def source_pins():
    files=[Path(__file__).resolve(),BASE/'candidate/oph_openwave/uniform_native.py',BASE/'candidate/oph_openwave/modular_native.py',
           BASE/'candidate/oph_openwave/geometry.py',BASE/'candidate/oph_openwave/construction.py']
    return {str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def run(native_state_path=None):
    before=source_pins()
    model=UniformFiberReduction(flat_ring(1,1,.23),6)
    native_duration=.1;rate=.29/(3*.37)
    native=model.evolve_per_seam(model.initial_density(),native_duration,rate)
    # Observable corner P M_87 P, with P resolving adjacent ports0 and1.
    indices=[0,1];corner=native[np.ix_(indices,indices)]
    probability=float(np.trace(corner).real)
    rho=corner/probability
    projector=np.zeros((87,87));projector[indices,indices]=1
    complement_projector=np.eye(87)-projector
    inside_branch=projector@native@projector
    outside_branch=complement_projector@native@complement_projector
    nonselective=inside_branch+outside_branch
    archive={'stored':False,'scope':'full native density and a separately computed single terminal P/Q Lüders instrument; no sequential measurement during native evolution'}
    if native_state_path is not None:
        native_state_path=Path(native_state_path)
        np.savez_compressed(native_state_path,native_density=native,
                            corner_projector=projector,complement_projector=complement_projector,
                            unnormalized_inside=inside_branch,unnormalized_complement=outside_branch,
                            nonselective_terminal_density=nonselective,conditional_corner_density=rho,
                            position_indices=np.asarray(indices))
        archive.update({'stored':True,'filename':native_state_path.name,
                        'sha256':hashlib.sha256(native_state_path.read_bytes()).hexdigest(),
                        'native_density_shape':[87,87],'inside_branch_probability':float(np.trace(inside_branch).real),
                        'complement_branch_probability':float(np.trace(outside_branch).real)})
    values,vectors=faithful_spectrum(rho)
    sx=np.array([[0.,1.],[1.,0.]])
    observable=vectors@sx@vectors.conj().T
    log_rho=(vectors*np.log(values))@vectors.conj().T
    generator=-log_rho
    # Standard finite GNS realization on Hilbert-Schmidt operators, column vec.
    delta_operator=np.kron(np.linalg.inv(rho).T,rho)
    delta_values,delta_vectors=np.linalg.eigh(delta_operator)
    delta_unitary=(delta_vectors*np.exp(.37j*np.log(delta_values)))@delta_vectors.conj().T
    left_A=np.kron(np.eye(2),observable)
    gns_error=float(np.max(abs(delta_unitary@left_A@delta_unitary.conj().T-
                               np.kron(np.eye(2),modular_flow(rho,observable,.37)))))
    sqrt_rho=(vectors*np.sqrt(values))@vectors.conj().T
    omega=sqrt_rho.reshape(-1,order='F')
    cyclic_state_error=float(np.linalg.norm(delta_operator@omega-omega))
    times=np.linspace(-4,4,161);points=[];concordance=0.;state_invariance=0.;square_error=0.
    delta=float(np.log(values[0])-np.log(values[1]))
    for s in times:
        evolved=modular_flow(rho,observable,float(s))
        direct=expm(1j*s*logm(rho))@observable@expm(-1j*s*logm(rho))
        correlation=complex(np.trace(rho@observable@evolved))
        closed=values[0]*np.exp(-1j*s*delta)+values[1]*np.exp(1j*s*delta)
        concordance=max(concordance,float(np.max(abs(evolved-direct))),abs(correlation-closed))
        state_invariance=max(state_invariance,abs(np.trace(rho@evolved)-np.trace(rho@observable)))
        square_error=max(square_error,float(np.max(abs(evolved@evolved-np.eye(2)))))
        points.append({'s':float(s),'C_real':correlation.real,'C_imag':correlation.imag,'C_abs':abs(correlation),
                       'observable_01_real':float(evolved[0,1].real),'observable_01_imag':float(evolved[0,1].imag)})
    h=1e-6
    derivative=(modular_flow(rho,observable,h)-modular_flow(rho,observable,-h))/(2*h)
    derivative_error=float(np.max(abs(derivative-1j*(log_rho@observable-observable@log_rho))))
    sample=.37;correct=complex(np.trace(rho@observable@modular_flow(rho,observable,sample)))
    wrong=complex(np.trace(rho@observable@modular_flow(rho,observable,-sample)))
    # With sigma_s=rho^(is) A rho^(-is), the finite KMS shift is -i.
    sigma_minus_i=rho@observable@np.linalg.inv(rho)
    kms=float(abs(np.trace(rho@observable@sigma_minus_i)-np.trace(rho@observable@observable)))
    sigma_plus_i=np.linalg.inv(rho)@observable@rho
    wrong_kms=float(abs(np.trace(rho@observable@sigma_plus_i)-np.trace(rho@observable@observable)))
    commuting_error=float(np.max(abs(modular_flow(rho,rho,.37)-rho)))
    mixed=np.eye(2)/2
    mixed_error=float(np.max(abs(modular_flow(mixed,sx,.37)-sx)))
    rejected=False
    try:modular_flow(np.diag([1.,0.]),sx,.37)
    except ValueError:rejected=True
    result={'schema':'oph.observer.modular_time.v1',
            'scope':'Finite observer-algebra modular correlations derived from native candidate output; no native-clock conversion or subjective-phenomenology derivation.',
            'native_source':{'model':'local-coarse modular one-load candidate','preparation':'uniform fixed-Gauss modular fiber at port0, distinct from definite records',
                             'geometry':'flat_ring(1,1)','positions':87,'modulus':6,'mass':.23,'record_energy':'none',
                             'duration':native_duration,'per_seam_rate':rate,'scheduler':'exact finite GKSL generator; numerical matrix-exponential action',
                             'full_native_trace':float(np.trace(native).real)},
            'native_state_archive':archive,
            'terminal_observer_instrument':{'definition':'rho -> P rho P or Q rho Q; P+Q=I',
                                            'inside_probability':probability,
                                            'complement_probability':float(np.trace(outside_branch).real),
                                            'nonselective_trace':float(np.trace(nonselective).real),
                                            'scope':'Computed terminal branches on the same native output; subsequent modular-flow algebra is the conditioned corner only, not the whole87state.'},
            'observer':{'algebra':'P M_87(C) P isomorphic to M_2(C), identity P',
                        'P_position_indices':indices,'selection':'adjacent physical ports0and1',
                        'corner_probability':probability,'complement_probability':float(np.trace(native).real-probability),
                        'state':'P rho_native P / Tr(P rho_native), explicitly conditioned on inside-P outcome',
                        'conditional_density':matrix_data(rho),'eigenvalues':values.tolist(),
                        'minimum_eigenvalue':float(min(values)),'numerical_faithfulness_threshold':1e-12,
                        'faithfulness_scope':'strict positive numerical spectrum on declared corner; no certified native-propagation error enclosure claimed',
                        'depolarizing_floor_added':False,'thermal_or_Gibbs_state_supplied':False,
                        'support_restriction':'none after choosing the declared two-dimensional corner; full87state faithfulness not assumed'},
            'modular_flow':{'definition':'sigma_s(A)=rho^(i*s) A rho^(-i*s)',
                            'modular_Hamiltonian':'K=-log(rho), sigma_s(A)=exp(-i*s*K) A exp(i*s*K)',
                            'K':matrix_data(generator),
                            'standard_GNS_modular_operator':'Delta = (rho^-1)^T tensor rho acting on column-vectorized Hilbert-Schmidt matrices',
                            'Delta':matrix_data(delta_operator),'Delta_eigenvalues':delta_values.tolist(),
                            'GNS_cyclic_separating_vector':'vec(sqrt(rho)); rho strictly positive on declared observer algebra',
                            'observable':matrix_data(observable),
                            'apparatus_choice':'Hermitian eigenframe sigma_x; noncommuting observer observable selected using measured state basis',
                            'dimensionless_modular_frequency':abs(delta),'native_duration_fixed':native_duration,
                            'physical_clock_conversion':'absent; if rho=exp(-beta H)/Z, given convention corresponds to physical Heisenberg t=-beta*hbar*s, conditionally',
                            'subjective_time_status':'candidate algebraic relational parameter only; no consciousness, experienced duration or phenomenology inferred'},
            'checks':{'observer_positive_faithful':bool(min(values)>1e-12),
                      'terminal_instrument_complete':bool(np.allclose(projector+complement_projector,np.eye(87)) and abs(np.trace(nonselective)-np.trace(native))<1e-12),
                      'noncommuting_observable':bool(np.linalg.norm(rho@observable-observable@rho)>1e-3),
                      'spectral_vs_general_matrix_function':concordance<1e-11,'state_invariance':state_invariance<1e-12,
                      'automorphism_square_identity':square_error<1e-12,'generator_sign_derivative':derivative_error<1e-7,
                      'finite_KMS_minus_i_sign':kms<1e-10,'GNS_operator_matches_flow':gns_error<1e-10,
                      'GNS_state_vector_fixed':cyclic_state_error<1e-10,'commuting_control_constant':commuting_error<1e-12},
            'mutations':{'reverse_modular_sign_changes_correlation':abs(correct-wrong)>.01,
                         'wrong_KMS_plus_i_rejected':wrong_kms>.01,'nonfaithful_state_rejected_without_floor':rejected,
                         'maximally_mixed_has_trivial_flow':mixed_error<1e-12},
            'errors':{'spectral_matrix_function_max':concordance,'invariant_state_max':float(state_invariance),
                      'automorphism_max':square_error,'derivative_max':derivative_error,'KMS_minus_i':kms,'wrong_KMS_plus_i':wrong_kms,'GNS_operator':gns_error,'GNS_state_vector':cyclic_state_error},
            'correlation_definition':'C(s)=Tr(rho A sigma_s(A)); observer state unchanged while observable correlations evolve',
            'plot_data':points,'full_physical_criterion_met':False,'official_rows_promoted':0,
            'assumptions':['finite one-load local-coarse candidate with supplied coherent law and mass',
                           'uniform modular fiber preparation is specified, not generated by a proven native preparation process',
                           'selected two-port observer corner and conditioning apparatus',
                           'numerically faithful observer density with no added noise floor',
                           'state-dependent logarithmic modular generator and chosen noncommuting observable',
                           'modular parameter is distinct from native clock; no phenomenological identification'],
            'source_sha256':before}
    if source_pins()!=before:raise ValueError("native or observer source changed during evaluation")
    result["checks"]={key:bool(value) for key,value in result["checks"].items()}
    result["mutations"]={key:bool(value) for key,value in result["mutations"].items()}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--plot-data',type=Path);p.add_argument('--native-state',type=Path,default=Path(__file__).with_name('native_state.npz'))
    args=p.parse_args();result=run(args.native_state)
    if not all(result['checks'].values()) or not all(result['mutations'].values()):raise SystemExit('modular-flow check failed')
    if args.output:args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if args.plot_data:
        with args.plot_data.open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(result['plot_data'][0]));writer.writeheader();writer.writerows(result['plot_data'])
    print(json.dumps({'checks':result['checks'],'mutations':result['mutations'],'observer_eigenvalues':result['observer']['eigenvalues'],
                      'modular_frequency':result['modular_flow']['dimensionless_modular_frequency']}))
