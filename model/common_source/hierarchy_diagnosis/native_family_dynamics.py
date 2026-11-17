"""The static core stabilizer and the actual completion-context response.

Only existing source operations enter: the protected-core reference, its
finite mixed-cap selector, fair native endpoint repair, and the native slow
band. A timeout remains a separate terminal output. The unread work map is a
one-proposal marginal, never a license to continue after that timeout.

The resulting noncommuting compressed response is not a Yukawa matrix or a
physical mass. The three-band compression is not an autonomous closed law.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from .. import slow_geometry
from ..particle_identification import family_attachment
from ...native_modular_reference import local_preparation
from ...observer import native


ROOT=Path(__file__).resolve().parents[3]
SOURCE_FILES=tuple(sorted(set((*local_preparation.SOURCE_PATHS,
    "model/common_source/slow_geometry.py",
    "model/common_source/particle_identification/family_attachment.py",
    "model/common_source/hierarchy_diagnosis/native_family_dynamics.py"))))


def _clean(matrix):
    return s.ImmutableSparseMatrix(s.SparseMatrix(matrix).applyfunc(s.expand))


def _zero(matrix):
    return all(s.expand(x)==0 or s.simplify(x)==0 for x in matrix)


def _entries(matrix):
    return [[str(s.expand(x)) for x in row] for row in matrix.tolist()]


@lru_cache(maxsize=2,typed=True)
def proposal_branches(selector_rounds=32):
    """Complete reduced accepted/timeout outcomes for every66-state input.

    Rows: input,terminal_flag,seam,coin,current,output,exact_probability.
    Fine auxiliary selector histories are already reduced by the existing
    exact source selector contract; no terminal outcome is conditioned away.
    """
    if type(selector_rounds) is not int or not 1<=selector_rounds<=32:
        raise ValueError("selector_rounds must be an integer from1 through32")
    selection=local_preparation.selectors(selector_rounds)
    weights=selection["completion_accepted_seam_weights"]
    timeout=selection["completion_timeout_probability"]
    if sum(weights.values(),F())+timeout!=1:
        raise ArithmeticError("the source selector lost a terminal outcome")
    states=native.states(binary=True); indices={x:i for i,x in enumerate(states)}
    branches=[]
    for index,state in enumerate(states):
        for seam,weight in sorted(weights.items()):
            for coin in (0,1):
                output,(_,current)=native.native(state,seam,coin)
                branches.append((index,"accepted",seam,coin,current,indices[output],weight/2))
        branches.append((index,"selector_timeout",None,None,None,index,timeout))
    return tuple(branches)


@lru_cache(maxsize=2,typed=True)
def _data(selector_rounds=32):
    c=slow_geometry.source()
    basis=s.ImmutableMatrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a,b in row]
                            for row in slow_geometry.raw_frame()])
    gram=_clean(basis.T*basis)
    left=s.ImmutableMatrix((gram.inv()*basis.T).applyfunc(s.simplify))
    projector=_clean(basis*left)
    retained=family_attachment.report()["compression"]["matrix_Qsqrt5"]
    feature=_clean(s.Matrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a,b in row] for row in retained]))
    if not feature.is_diagonal() or len(feature.eigenvals())!=3:
        raise ArithmeticError("the existing native feature must have three simple axes")
    states=native.states(binary=True)
    sigma=s.ImmutableMatrix(list(map(s.Rational,local_preparation.reference())))
    selection=local_preparation.selectors(selector_rounds)
    weights=selection["completion_accepted_seam_weights"]
    timeout=s.Rational(selection["completion_timeout_probability"])
    source=tuple(proposal_branches(selector_rounds))
    entries=Counter(); column_mass=[s.Rational(0)]*66
    for before,flag,seam,coin,current,after,probability in source:
        if (sum(states[before][p] for p in local_preparation.CORE)
                !=sum(states[after][p] for p in local_preparation.CORE)):
            raise ArithmeticError("a completion outcome changed protected core charge")
        if sigma[before]!=sigma[after]:
            raise ArithmeticError("the reference is not constant on an accepted local swap block")
        entries[after,before]+=s.Rational(probability)
        column_mass[before]+=s.Rational(probability)
    if column_mass!=[s.Rational(1)]*66:
        raise ArithmeticError("the complete one-proposal instrument is not normalized")
    transition=s.ImmutableSparseMatrix(66,66,dict(entries))
    if transition*sigma!=sigma:
        raise ArithmeticError("the actual reduced completion map does not fix its prepared reference")
    accepted=transition-timeout*s.eye(66)
    if accepted*sigma!=(1-timeout)*sigma:
        raise ArithmeticError("success and terminal timeout reference branches were conflated")
    mean12=s.eye(12)*timeout
    for (a,b),weight in weights.items():
        endpoint=s.eye(12)
        endpoint[a,a]=endpoint[a,b]=endpoint[b,a]=endpoint[b,b]=s.Rational(1,2)
        mean12+=s.Rational(weight)*endpoint
    mean12=_clean(mean12)
    loads=s.ImmutableMatrix(12,66,lambda i,j:states[j][i])
    if loads*transition!=mean12*loads:
        raise ArithmeticError("first-moment response differs from exact native branch enumeration")
    compressed=_clean(left*mean12*basis)
    leakage=_clean((s.eye(12)-projector)*mean12*basis)
    leakage_gram=s.ImmutableSparseMatrix((leakage.T*leakage/gram[0,0]).applyfunc(s.simplify))
    leakage_norm=s.simplify(s.trace(leakage_gram))
    if leakage_norm.is_positive is not True:
        raise ArithmeticError("the completion update unexpectedly preserves the slow band")
    twice=_clean(left*mean12*mean12*basis)
    if not _zero(twice-compressed*compressed-leakage_gram):
        raise ArithmeticError("the exact second-step projection memory identity failed")

    # The reference, as distinct from the whole preparation context, has one
    # nontrivial proper-rotation stabilizer. Recover it from the source group.
    core=set(local_preparation.CORE)
    stabilizers=[g for g in c.rotations() if {g[i] for i in core}==core]
    if len(stabilizers)!=2:
        raise ArithmeticError("the protected unlabelled core stabilizer changed")
    other=next(g for g in stabilizers if tuple(g)!=tuple(range(12)))
    inverse=tuple(other.index(i) for i in range(12))
    residual=_clean(left*s.Matrix([basis[inverse[i],:] for i in range(12)]))
    if residual!=s.diag(-1,1,-1):
        raise ArithmeticError("the recovered residual half-turn is not the existing y axis")
    allowed=[(i,j) for i in range(3) for j in range(3) if residual[i,i]==residual[j,j]]
    if len(allowed)!=5 or not _zero(residual*feature-feature*residual):
        raise ArithmeticError("the static reference commutant changed")
    if _zero(residual*compressed-compressed*residual) or _zero(feature*compressed-compressed*feature):
        raise ArithmeticError("the full completion context did not break static family symmetry")
    edges=[(i,j) for i in range(3) for j in range(i+1,3) if s.expand(compressed[i,j])!=0]
    reached={0}
    for _ in range(3):
        reached|={j for i,j in edges if i in reached}|{i for i,j in edges if j in reached}
    if reached!={0,1,2}:
        raise ArithmeticError("noncommuting update does not connect all three feature axes")
    # Simple F gives coordinate spectral projectors. Sandwiching R and
    # dividing each nonzero entry yields actual matrix units, not a numeric
    # rank inference or an assumed version of Burnside's theorem.
    matrix_units=[]
    for i,j in edges:
        ei,ej=s.zeros(3),s.zeros(3);ei[i,i]=ej[j,j]=1
        unit=_clean(ei*compressed*ej/compressed[i,j])
        target=s.zeros(3);target[i,j]=1
        if not _zero(unit-target):
            raise ArithmeticError("source algebra matrix-unit witness failed")
        matrix_units.append((i,j))
    if len(edges)!=3:
        raise ArithmeticError("this receipt expects all three existing nonzero response links")
    # Generic two-family commutator cubes have zero trace. The invariant
    # third line therefore forbids a physical three-family Jarlskog phase.
    a=s.Matrix(2,2,s.symbols("a0:4"));b=s.Matrix(2,2,s.symbols("b0:4"))
    if s.expand(s.trace((a*b-b*a)**3))!=0:
        raise ArithmeticError("static2plus1 Jarlskog obstruction failed")
    return {"basis":basis,"F":feature,"R":compressed,"T":mean12,"transition":transition,
            "timeout":timeout,"sigma":sigma,"residual":residual,"residual_permutation":other,
            "allowed":allowed,"edges":edges,"leakage":leakage,"leakage_gram":leakage_gram,
            "leakage_norm":leakage_norm,"source":source}


def compressed_response(selector_rounds=32):
    """The one-proposal work-marginal compression, not closed family dynamics."""
    return _data(selector_rounds)["R"]


def proposal_operator(selector_rounds=32):
    """66-state one-proposal work marginal; its timeout flag is separately retained."""
    return _data(selector_rounds)["transition"]


def report(selector_rounds=32):
    data=_data(selector_rounds)
    timeout=data["timeout"]
    residual=data["residual"]
    mapped_caps=[sorted(data["residual_permutation"][i] for i in cap) for cap in local_preparation.CAPS]
    serial=[[before,flag,list(seam) if seam is not None else None,coin,current,after,str(probability)]
            for before,flag,seam,coin,current,after,probability in data["source"]]
    bindings=dict(family_attachment.report()["source_bindings"])
    for name in SOURCE_FILES:
        bindings[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return {
        "schema":"oph.hierarchy_diagnosis.native_family_dynamics.v1",
        "source":{"reference":"existing ideal protected-core66-state reference",
                  "update":"existing finite mixed-cap completion selector and native fair endpoint repair",
                  "caps":[list(x) for x in local_preparation.CAPS],"selector_rounds":selector_rounds,
                  "time_parameter":"one completion proposal; not native elapsed attempts, modular time or physical time",
                  "finite_reference_preparation_errors_propagated":False,
                  "new_empirical_mass_or_coupling_inputs":[]},
        "complete_reduced_instrument":{"rows":len(data["source"]),"input_states":66,
                                      "accepted_seam_coin_current_and_timeout_outcomes_retained":True,
                                      "terminal_timeout_probability":str(timeout),
                                      "accepted_probability":str(1-timeout),
                                      "terminal_timeout_work_map":"tau Id with a retained terminal flag; never continued",
                                      "unread_work_marginal":"T=sum_e w_e E_e + tau I, for one proposal only",
                                      "all_columns_normalized_exact":True,
                                      "full_diagonal_reference_stationary_exact_for_C_and_E":True,
                                      "reference_constant_on_every_accepted_local_swap_block":True,
                                      "C_E_diagonal_fix_reason":"Both native instruments have the same half-identity plus half-swap population map; sigma is constant on each swap block.",
                                      "success_reference_output":"(1-tau)*sigma",
                                      "timeout_reference_output":"tau*sigma",
                                      "fine_auxiliary_selector_tapes_enumerated":False,
                                      "existing_exact_reduced_selector_contract_used":True,
                                      "branch_table_sha256":hashlib.sha256(json.dumps(serial,separators=(",",":"),sort_keys=True).encode()).hexdigest()},
        "static_reference_obstruction":{"reference_proper_stabilizer_size":2,
                                        "nonidentity_port_permutation":list(data["residual_permutation"]),
                                        "half_turn_on_family":_entries(residual),
                                        "commutant_complex_dimension":5,
                                        "commutant":"M2(C) on native(x,z) plus C on nativey",
                                        "allowed_endomorphism_entries":[list(x) for x in data["allowed"]],
                                        "any_covariant_endomorphism_of_static_reference_commutes_half_turn":True,
                                        "quark_mixing_from_only_static_reference":"at most2plus1; one common unmixed familyline",
                                        "weak_CP_Jarlskog_from_only_static_reference_vanishes":True},
        "full_preparation_context":{"caps_after_static_half_turn":mapped_caps,
                                    "same_caps_after_half_turn":set(map(tuple,mapped_caps))==set(local_preparation.CAPS),
                                    "actual_T12":_entries(data["T"]),
                                    "first_moment_response_checked_on_all66_source_states":True,
                                    "compressed_R3":_entries(data["R"]),
                                    "F_R_commutator":_entries(data["F"]*data["R"]-data["R"]*data["F"]),
                                    "half_turn_R_commutator":_entries(residual*data["R"]-data["R"]*residual),
                                    "nonzero_family_response_links":[list(x) for x in data["edges"]],
                                    "generated_complex_matrix_algebra_dimension":9,
                                    "algebra_proof":"Simple F gives all three spectral projectors; nonzero P_i R P_j divided by R_ij gives every offdiagonal matrix unit.",
                                    "algebraic_three_family_mixing_capacity":True,
                                    "all_retained_source_matrices_real":True,
                                    "physical_CP_phase_generated":False},
        "projection_boundary":{"slow_band_preserved_by_T12":False,
                               "outside_band_map_12_by3":_entries(data["leakage"]),
                               "outside_band_squared_norm_exact":str(data["leakage_norm"]),
                               "outside_band_squared_norm_float_diagnostic":float(data["leakage_norm"]),
                               "normalized_leakage_Gram":_entries(data["leakage_gram"]),
                               "second_step_identity":"left*T^2*B = R^2 + leakage^T*leakage/(B^T B)_00",
                               "second_step_identity_checked_exact":True,
                               "second_step_identity_scope":"Algebraic one-proposal marginal comparison; not a protocol that restarts after timeout.",
                               "autonomous_three_family_repair_dynamics_established":False,
                               "discarded_complement_or_failure_modes":False},
        "scope":{"existing_full_update_context_supplies_more_than_static_F":True,
                 "Yukawa_map_selected":False,"response_eigenvalues_identified_as_particle_masses":False,
                 "coherent_quartic_fermion_interaction_derived":False,
                 "Standard_Model_masses_or_CKM_derived":False},
        "source_bindings":dict(sorted(bindings.items()))}


if __name__=="__main__":
    print(json.dumps(report(),indent=2,allow_nan=False))
