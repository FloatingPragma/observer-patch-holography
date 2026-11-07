#!/usr/bin/env python3
"""Complete observer records for the joined charged-lepton control.

The data are simulated detector outcomes for a declared finite model, not
experimental particle data and not a source-selected mass prediction.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import bmat,csr_matrix
from scipy.sparse.linalg import expm_multiply

HERE=Path(__file__).resolve().parent
from .higgs_join.interacting_carrier import build_join, primary_family_matrices
from .recorded_response import (record_law,sample_diagonal_records,
                              independent_full_control,ramsey_kraus_actions)


def complex_json(z):
    return {"real":float(complex(z).real),"imag":float(complex(z).imag)}


def direct_joined_control(model,law,t=.371):
    """Independent full-carrier branch evolution, including a cross-family pair."""
    sectors=model["sectors"];h0,h1=sectors[0]["H"],sectors[1]["H"]
    whole=bmat([[h0,None],[None,h1]],format="csr")
    n0,n1=h0.shape[0],h1.shape[0]
    psi=np.r_[law["ground"],np.zeros(n1,complex)]
    probe=list(model["addition"].values())
    errors=[]
    for j,k in ((0,0),(0,2),(1,4)):
        def apply_w(b,state,sgn):
            u,v=state[:n0],state[n0:]
            return np.r_[u-b.getH()@(b@u)-1j*sgn*(b.getH()@v),
                         v-b@(b.getH()@v)-1j*sgn*(b@u)]
        evolved=expm_multiply(-1j*t*whole,psi)
        branch0=apply_w(probe[j],evolved,1)
        expected=record_law(h0,h1,[probe[j],probe[k]],[t],
                    spectral0=(law["energies0"],law["basis0"]),
                    spectral1=(law["energies1"],law["basis1"]))
        for sindex,sgn in enumerate((1,-1)):
            branch1=expm_multiply(-1j*t*whole,apply_w(probe[k],psi,sgn))
            f=np.vdot(branch0,branch1)
            reference=expected["blocked_correlator"][0,0,1]+sgn*expected["addition_correlator"][0,0,1]
            errors.append(abs(f-reference))
            errors.extend((abs(np.vdot(branch0,branch0)-1),abs(np.vdot(branch1,branch1)-1)))
    if max(errors)>2e-10:raise ArithmeticError("Independent joined branch-evolution check")
    return {"system_dimension":n0+n1,"with_control_qubit_dimension":2*(n0+n1),
            "probe_pairs_checked":[[0,0],[0,2],[1,4]],
            "maximum_complex_signal_or_norm_error":float(max(errors)),
            "all_evolution_forward_in_time":True}


def run(out=None,shots=4096):
    if out is not None:
        out=Path(out);out.mkdir(parents=True,exist_ok=True)
    y,sigma=primary_family_matrices()
    model=build_join(family0=y,family1=sigma,yukawa=y,fermion_energy=0.,scalar_energy=1.)
    h0,h1=(model["sectors"][q]["H"] for q in (0,1))
    # A short window is adequate to check a record law; it does not resolve
    # arbitrarily close response lines. No Fourier mass fitting is performed.
    times=np.linspace(0,12.8,129)
    law=record_law(h0,h1,list(model["addition"].values()),times)
    e0,u0,e1,u1=(law[k] for k in ("energies0","basis0","energies1","basis1"))
    recorded=sample_diagonal_records(law,shots=shots)
    independent=direct_joined_control(model,law)
    # Expose and check actual retained postmeasurement states, not only p(r).
    probe=list(model["addition"].values())
    initial=np.r_[law["ground"],np.zeros(h1.shape[0],complex)]
    rng=np.random.default_rng(2062)
    random_input=rng.normal(size=len(initial))+1j*rng.normal(size=len(initial))
    random_input/=np.linalg.norm(random_input)
    branch_packet={};branch_errors=[]
    for name,state in (("prepared",initial),("arbitrary",random_input)):
        for quadrature in ("X","Y"):
            outputs=ramsey_kraus_actions(h0,h1,probe[0],probe[2],state,.371,-1,quadrature)
            probability=[float(np.vdot(x,x).real) for x in outputs]
            branch_errors.append(abs(sum(probability)-1))
            for index,x in enumerate(outputs):
                branch_packet[f"{name}_{quadrature}_{index}_unnormalized_poststate"]=x
    if max(branch_errors)>2e-10:raise ArithmeticError("Complete back-action normalization")
    tiny=independent_full_control()
    eigen_errors=[]
    for h,e,u in ((h0,e0,u0),(h1,e1,u1)):
        residual=h@u-u*e
        eigen_errors.append({"maximum_column_residual":float(np.linalg.norm(residual,axis=0).max()),
            "orthonormality_Frobenius_error":float(np.linalg.norm(u.conj().T@u-np.eye(len(e)),"fro"))})
    amps=law["amplitudes"];weights=np.einsum("ni,nj->nij",amps.conj(),amps)
    sumerr=float(np.linalg.norm(weights.sum(axis=0)-law["Gram"],2))
    if sumerr>2e-10:raise ArithmeticError("Complete addition measure")
    frequencies=e1-e0[0]
    # Group true spin degeneracies with a conservative numerical threshold.
    # No family lines below this tolerance are discarded; ungrouped arrays
    # are always retained. This is an audit display, not a fitted resolution.
    clusters=[]
    for i,w in enumerate(frequencies):
        if not clusters or abs(w-clusters[-1]["anchor"])>2e-12*max(1,abs(w)):
            clusters.append({"anchor":float(w),"indices":[i]})
        else:clusters[-1]["indices"].append(i)
    first=[]
    for group in clusters[:6]:
        ids=group["indices"]
        total=weights[ids].sum(axis=0)
        first.append({"frequency_mean":float(np.mean(frequencies[ids])),
            "numerical_grouping_spread":float(np.ptp(frequencies[ids])),
            "multiplicity":len(ids),"matrix_residue_trace":float(np.trace(total).real),
            "diagonal_residues":np.diag(total).real.tolist()})
    # Exact common-charge ambiguity:
    # epsilon_F NF+epsilon_H NH = epsilon_F Q*+(2epsilon_F+epsilon_H) NH.
    # Keep the NH coefficient fixed while changing epsilon_F by delta.
    delta=.125
    shifted=build_join(family0=y,family1=sigma,yukawa=y,fermion_energy=delta,scalar_energy=1-2*delta)
    ambiguity=[]
    for q in (0,1):
        difference=(shifted["sectors"][q]["H"]-model["sectors"][q]["H"]).toarray()-delta*q*np.eye(model["sectors"][q]["dimension"])
        ambiguity.append(float(np.linalg.norm(difference,2)))
    if max(ambiguity)>1e-11:raise ArithmeticError("Conserved-charge offset identity")
    # In this restricted lepton sector q6Y(matter+scalar)=6Q*.
    # A common gauge potential shifts the compensating reference oppositely.
    # When treated consistently, that gauge change cancels from every record.
    compensated=build_join(family0=y,family1=sigma,yukawa=y,fermion_energy=delta,
                           scalar_energy=1-2*delta,reference_energy=-delta)
    compensation_errors=[float(np.linalg.norm(
        (compensated["sectors"][q]["H"]-model["sectors"][q]["H"]).toarray(),2)) for q in (0,1)]
    if max(compensation_errors)>1e-11:raise ArithmeticError("Common matter/reference gauge-offset cancellation")
    total_weight=float(np.trace(law["Gram"]).real)
    width=float(times[-1]-times[0])
    lowest=float(frequencies[0])
    results={"scope":model["scope"],"parameters":model["parameters"],
        "prepared_family_matrices": "audited native Y and ordered-history sigma(log(5/4)); placements and strengths are declared",
        "initial_state":"unique ground of the declared interacting Q*=0 Hamiltonian",
        "physical_preparation_of_that_ground_from_native_records_not_derived":True,
        "ground_energy":float(e0[0]),"ground_gap":float(e0[1]-e0[0]),
        "clock":"dimensionless parameter of the declared Hamiltonian; native/SI calibration not derived",
        "reference_energy":0.,"reference_energy_is_a_declared_detector_calibration":True,
        "probe_order":[list(key) for key in model["addition"]],
        "record_axes":["time","relative sign (+,-)","detector j","detector k","quadrature X/Y","outcome +/-"],
        "checks":law["checks"],"independent_full_joined_control":independent,
        "full_backaction_normalization_error":float(max(branch_errors)),
        "complete_Kraus_output_vectors_retained":True,
        "independent_small_ancilla_control":tiny,"eigensystem_checks":eigen_errors,
        "spectral_sum_rule_error":sumerr,
        "total_six_probe_addition_weight":total_weight,
        "all_1596_response_lines_retained":True,"signed_frequency_minimum":lowest,
        "signed_frequency_maximum":float(frequencies[-1]),
        "first_six_numerically_grouped_lines":first,
        "numerical_grouping_is_not_experimental_resolution":True,
        "records":recorded["receipt"],
        "finite_window":{"start":float(times[0]),"end":float(times[-1]),
            "time_step":float(times[1]-times[0]),
            "lowest_line_phase_accumulation":lowest*width,
            "lowest_line_one_cycle_time":float(2*np.pi/lowest) if lowest>0 else None,
            "short_record_window_does_not_identify_all_response_lines":True},
        "charge_offset_audit":{"delta":delta,"operator_identity_residuals":ambiguity,
            "all_addition_frequencies_shift_by":delta,
            "within_sector_record_laws_unchanged":True,
            "reference_energy_was_held_fixed_in_that_comparison":True,
            "compensating_reference_shift":-delta,
            "consistent_common_gauge_shift_residuals":compensation_errors,
            "intrinsic_mass_arbitrariness_is_not_implied":True},
        "fermion_parity":"probe creates the matter fermion and an odd compensating reference; compound operation is even. Reference parity assignment is part of the declared measurement protocol.",
        "explicit_probe_CAR_order":"b_j = c_j^dagger d_reference^dagger; b_j^dagger = d_reference c_j. Reference mode is ordered after all matter modes; its two signs are +1 on this carrier. Reversing creation order changes the sign.",
        "no_measured_mass_fit":True,"no_claim_of_autonomous_source_selected_vacuum":True}
    arrays=dict(
        times=times,probabilities=law["probabilities"],
        addition_correlator=law["addition_correlator"],blocked_correlator=law["blocked_correlator"],
        frequencies=frequencies,matrix_weights=weights,probe_Gram=law["Gram"],
        packed_binary_records=recorded["packed_binary_records"],
        plus_counts=recorded["plus_counts"],plus_probabilities=recorded["plus_probabilities"],
        estimated_addition_quadratures=recorded["estimated_addition_quadratures"],
        exact_addition_quadratures=recorded["exact_addition_quadratures"],**branch_packet)
    if out is not None:
        np.savez_compressed(out/"complete_interacting_records.npz",**arrays)
        (out/"interacting_record_results.json").write_text(json.dumps(results,indent=2)+"\n")
    return results,arrays


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--shots",type=int,default=4096)
    a=p.parse_args();result,_=run(a.out,a.shots);print(json.dumps(result,indent=2))
