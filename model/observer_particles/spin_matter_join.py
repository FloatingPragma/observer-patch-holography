#!/usr/bin/env python3
"""Construct and execute the conditional independent Spin x matter join.

The source spin cover and existing conditional M15 are retained. Distinct
factor actions, locked diagonal transport, complete proper-code instruments
and their relative-phase readouts are never identified with each other.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PUBLIC=ROOT
sys.path.insert(0,str(PUBLIC))
sys.path.insert(0,str(PUBLIC/'code/a5_closure'))
import super_tannakian_matter_lift_certificate as matter
from model.observer.source_scalar import compiler
from model.observer.source_discrimination import circuit


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def zero(a):return all(s.expand(v)==0 for v in a)
def from_sparse(p):return s.SparseMatrix(*p['shape'],{(i,j):s.sympify(v) for i,j,v in p['entries']})
def numeric(a):return np.array(a.evalf(),dtype=complex)


def exact_join():
    packet=json.loads((HERE/'results/representation_identification_operators.json').read_text())
    basis=[tuple(x) for x in packet['basis_subsets']]
    indices=packet['even_indices'];selected=[basis[i] for i in indices]
    # Reuse the actual public exterior functor to transport all twelve
    # supplied source-current lifts, not a new gauge-field-name table.
    fock=matter.FockSpace()
    assert fock.subsets==basis
    def c5(v):
        v=s.expand(v)
        re,im=s.expand(s.re(v)),s.expand(s.im(v))
        def f(x):return matter.F5(str(x.coeff(s.sqrt(5),0)),str(x.coeff(s.sqrt(5),1)))
        return matter.C5(f(re),f(im))
    def smatrix(a):
        return s.SparseMatrix(len(a),len(a[0]),{
            (i,j):s.Rational(str(v.re.a))+s.Rational(str(v.re.b))*s.sqrt(5)
                   +s.I*(s.Rational(str(v.im.a))+s.Rational(str(v.im.b))*s.sqrt(5))
            for i,row in enumerate(a) for j,v in enumerate(row) if not v.is_zero()})
    gauge=[]
    for encoded in packet['spin_lifted_one_body_generators']:
        a=from_sparse(encoded)
        dg=fock.dgamma([[c5(a[i,j]) for j in range(5)] for i in range(5)])
        gauge.append(smatrix(dg).extract(indices,indices))
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    spin=[s.kronecker_product(p/2,s.eye(15)) for p in pauli]
    gauge30=[s.kronecker_product(s.eye(2),g) for g in gauge]
    for j in spin:
        for g in gauge30:assert zero(j*g-g*j)
    for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
        assert zero(spin[i]*spin[j]-spin[j]*spin[i]-s.I*spin[k])
    # Two Weyl index representations have the same spatial Spin action.
    # These are finite field-index boosts, not unitary Hilbert-space boosts.
    for handedness in (-1,1):
        boost=[handedness*s.I*j for j in spin]
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
            assert zero(boost[i]*boost[j]-boost[j]*boost[i]+s.I*spin[k])
            assert zero(spin[i]*boost[j]-boost[j]*spin[i]-s.I*boost[k])
    internal_centre=s.diag(*[(-1)**sum(i>=3 for i in x) for x in selected])
    weak_centre=s.kronecker_product(s.eye(2),internal_centre)
    spin_centre=-s.eye(30)
    locked_centre=spin_centre*weak_centre
    assert list(weak_centre.diagonal()).count(-1)==16
    assert list(locked_centre.diagonal()).count(-1)==14
    # Construct an actual source lift of a pi rotation about y. Quaternion
    # y has S=-i sigma_y, and its spatial R=diag(-1,1,-1).
    u=s.Matrix([[0,-1],[1,0]])
    rotation=s.diag(-1,1,-1)
    pi=matter.carrier_implementer(
        [[matter.F5(int(rotation[i,j])) for j in range(3)] for i in range(3)],
        [[c5(u[i,j]) for j in range(2)] for i in range(2)])
    transported=smatrix(fock.exterior_lift(pi)).extract(indices,indices)
    assert transported*transported==internal_centre
    independent=s.kronecker_product(u,s.eye(15))
    locked=s.kronecker_product(u,transported)
    assert independent*independent==spin_centre
    assert locked*locked==locked_centre
    # There is no one-dimensional twist of the measured double cover that
    # flips its centre: centre=-1 is the quaternion i,j commutator.
    qi,qj=-s.I*pauli[0],-s.I*pauli[1]
    commutator=qi*qj*qi.H*qj.H
    assert commutator==-s.eye(2)
    identification={k:from_sparse(v) for k,v in packet['even_identification_operators'].items()}
    for a in identification.values():
        for j in spin:assert zero(j*s.kronecker_product(s.eye(2),a)-s.kronecker_product(s.eye(2),a)*j)
    return {
        'dimension':30,'matter_dimension':15,'physical_spin_factor_dimension':2,
        'source_current_generators':len(gauge),'exact_spin_gauge_commutators':36,
        'independent_spin_centre':{'minus_one':30,'plus_one':0},
        'internal_weak_centre':{'minus_one':16,'plus_one':14},
        'locked_diagonal_centre':{'minus_one':14,'plus_one':16},
        'both_Weyl_index_extensions_satisfy_Lorentz_algebra':True,
        'source_rotation_data_alone_selects_physical_chirality':False,
        'quaternion_commutator_is_central_minus_one':True,
        'scalar_character_cannot_flip_lift_centre':True,
        'commutation_scope':'Independent factor actions on one finite module, not a source-selected spin connection or particle Hamiltonian.',
        'selected_matter_subsets':[list(x) for x in selected],
    },numeric(transported).real


def monomial_word(matrix):
    n=len(matrix);bits=n.bit_length()-1
    assert n==1<<bits and np.array_equal(matrix.T@matrix,np.eye(n))
    mapping=np.argmax(abs(matrix),axis=0)
    assert all(np.count_nonzero(matrix[:,j])==1 for j in range(n))
    seen=set();word=[]
    for root in range(n):
        if root in seen:continue
        cycle=[];i=root
        while i not in seen:
            seen.add(i);cycle.append(i);i=int(mapping[i])
        for destination in cycle[1:]:word+=circuit.basis_swap_word(cycle[0],destination,bits)
    signs=np.ones(n)
    for j in range(n):signs[mapping[j]]=matrix[mapping[j],j]
    word+=compiler.compile_diagonal(np.where(signs<0,math.pi,0.0))
    return word


def controlled_y(angle,control,target):
    return [compiler.pulse('Y',angle/2,[target])]+compiler.cnot(control,target)+[
        compiler.pulse('Y',-angle/2,[target])]+compiler.cnot(control,target)


def execute_join(transported):
    # Tensor ordering: spin is bit4, matter0..14 in low four bits.
    # Both padding addresses15,31 are retained in every complete matrix.
    matter16=np.eye(16);matter16[:15,:15]=transported
    internal_word=monomial_word(matter16)
    spin_word=[compiler.pulse('Y',-math.pi/2,[4])]
    actual_g=compiler.execute(internal_word,5)
    actual_s=compiler.execute(spin_word,5)
    expected_g=np.kron(np.eye(2),matter16)
    expected_s=np.kron(np.array([[0,-1],[1,0]]),np.eye(16))
    joint=compiler.execute(spin_word+internal_word,5)
    reverse=compiler.execute(internal_word+spin_word,5)
    compensated=compiler.execute(spin_word+internal_word+compiler.inverse(internal_word),5)
    locked_twice=compiler.execute((spin_word+internal_word)*2,5)
    expected_locked=np.linalg.matrix_power(expected_s@expected_g,2)
    errors={'spin':float(np.linalg.norm(actual_s-expected_s,2)),
            'internal_source_transport':float(np.linalg.norm(actual_g-expected_g,2)),
            'join':float(np.linalg.norm(joint-expected_s@expected_g,2)),
            'reverse_order':float(np.linalg.norm(reverse-joint,2)),
            'chosen_gauge_compensator':float(np.linalg.norm(compensated-expected_s,2)),
            'locked_loop':float(np.linalg.norm(locked_twice-expected_locked,2))}
    assert max(errors.values())<3e-11
    # A complete relative-phase readout of the actual two-pulse spin loop.
    # Pointer is bit5; its blank and full32 work inputs are retained.
    prepare=[compiler.pulse('Y',-math.pi/4,[5])]
    loop=controlled_y(-math.pi/2,5,4)*2
    word=prepare+loop+compiler.inverse(prepare)
    inputs=np.eye(64,dtype=complex)[:,:32]
    final=compiler.execute(word,6,inputs)
    # The real Y preparation is Q=exp(-i*pi*Y/4), so Q^dagger Z Q=-X.
    # Retain that exact global sign in the phase-sensitive comparison.
    expected=np.vstack((np.zeros((32,32)),-np.eye(32)))
    phase_error=float(np.linalg.norm(final-expected,2))
    assert phase_error<3e-11
    laws=abs(final)**2
    code=[i for i in range(32) if i%16<15]
    padding=[15,31]
    leakage=float(np.linalg.norm(joint[np.ix_(padding,code)],2))
    exports={}
    for d in (2,4):
        ks=compiler.transfer_kraus(d,True)
        assert np.array_equal(sum(k.T.conj()@k for k in ks),np.eye(6))
        exports[str(d)]={'Kraus_branches':len(ks),'failure_branches':6-d,
                         'failure_probability_for_ideal_valid_code':0,'completeness':True}
    return {
        'workspace_qubits':5,'physical_code_dimension':30,'padding_dimension':2,
        'Spin_pulse_count':len(spin_word),'internal_transport_pulse_count':len(internal_word),
        'full32_column_matrix_errors_operator_norm':errors,
        'physical_code_to_padding_leakage_operator_norm':leakage,
        'loop_pointer_qubits':6,'loop_pointer_source_pulses':len(word),
        'loop_pointer_complete64_outcome_laws_by32_inputs':laws.T.tolist(),
        'loop_pointer_state_error_operator_norm':phase_error,
        'loop_pointer_all_work_and_padding_inputs_retained':True,
        'loop_pointer_plus_probability_max':float(laws[:32].sum(axis=0).max()),
        'loop_pointer_minus_probability_min':float(laws[32:].sum(axis=0).min()),
        'raw_norm_error_max':float(np.max(abs(laws.sum(axis=0)-1))),
        'complete_M6_exports':exports,
        'apparatus_contract':'Existing coherent code gluing, reset, finite finite buffers, packing and complete success/failure exports. No arbitrary physical U30 primitive is introduced.',
        'execution_scope':'Native source-response pulses on the declared code. No new spacetime connection, covariant hopping, field vacuum, complete many-fermion construction or selected mass is asserted.'}


def main():
    exact,transport=exact_join()
    result={'schema':'oph.finite.spin_matter_join.v1',
        'status':'PASS_CONSTRUCTED_AND_SOURCE_EXECUTED_CONDITIONAL_SPIN_MATTER_JOIN',
        'exact_representation':exact,'native_execution':execute_join(transport),
        'positive_completion':'A single30-dimensional representation with independent commuting source-Spin and conditional matter actions, uniform physical-spin central−I, and source-executed joint transformations exists under the existing coherent-code contract.',
        'source_selection_status':{'Spin_matter_tensor_attachment_chosen':True,
            'local_Spin_connection_selected':False,'same_charged_propagator_emitted':False,
            'statistics_selected':False,'masses_selected':False},
        'alternative_route_result':'Locking the geometric and internal source actions does not make all matter fermionic. Independent factor actions construct the valid conditional route; a scalar twist cannot repair the locked centre, while an explicitly chosen inverse internal gauge transport does recover the independent Spin action and is source-executed here.',
        'producer_sha256':sha(__file__)}
    paths=[HERE/'results/representation_identification_operators.json',
           HERE/'representation_identification.py',HERE/'fermion_species_identification.py']
    for module in tuple(sys.modules.values()):
        file=getattr(module,'__file__',None)
        if file and Path(file).is_file() and Path(file).is_relative_to(PUBLIC):paths.append(Path(file))
    paths += [PUBLIC/p for p in compiler.source_summary()['source_sha256']]
    paths += [PUBLIC/p for p in ['docs/research/COHERENT_SOURCE_CLOCKS.md',
            'docs/research/FERMIONIC_SOURCE_CLOCKS.md','Lean/Screen/FermionSectorAssembly.lean',
            'code/a5_closure/manifests/spin_statistics_semantic_artifact.json']]
    result['source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))}
    target=HERE/'results/spin_matter_join.json'
    target.write_text(json.dumps(public_receipt(result),indent=2,sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'output':str(target),
                      'execution_errors':result['native_execution']['full32_column_matrix_errors_operator_norm'],
                      'phase_error':result['native_execution']['loop_pointer_state_error_operator_norm']}))


if __name__=='__main__':main()
