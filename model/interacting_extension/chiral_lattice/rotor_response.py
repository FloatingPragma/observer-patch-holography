#!/usr/bin/env python3
"""Operator-valued quantum-link interaction on one closed four-link graph.

This is an exact Gauss reduction of the native U(1) subgroup coordinate algebra,
with all five native hypercharges represented.  It is not a full nonabelian,
chiral, or 3+1D continuum calculation.  Electric fields remain quantum variables.
The finite flux boundary is retained and its unitary-link defect is reported.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from numbers import Integral, Real
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.special import gammainc
from .source_data import load, ROOT

def _integer(value, name, minimum=None):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return int(value)

def _real(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite and real")
    return float(value)

def build_flux_ring(q,K,L=4,electric=1.,loop=1.,hopping=1.):
    """Return portable H and complete physical basis |position x,loop flux ell>.

    Initial/root static charge is -q.  All oriented-edge electric fields are
    E_e=ell+q*1[e<x].  This solves EVERY abelian Gauss constraint exactly.
    H_E=(electric/2) sum E_e^2; H_B=-loop*(W+W^dag)/2;
    H_hop=-hopping*(T+T^dag).
    """
    q=_integer(q,"q");K=_integer(K,"K",0);L=_integer(L,"L",3)
    electric=_real(electric,"electric");loop=_real(loop,"loop");hopping=_real(hopping,"hopping")
    if electric<=0:raise ValueError("electric must be positive")
    labels=tuple((x,ell) for x in range(L) for ell in range(-K,K+1))
    lookup={v:i for i,v in enumerate(labels)};n=len(labels)
    flux=np.array([[ell+q*(e<x) for e in range(L)] for x,ell in labels],int)
    gauss=np.array([[q*(v==x)-q*(v==0)+flux[i,v]-flux[i,(v-1)%L]
                     for v in range(L)] for i,(x,ell) in enumerate(labels)])
    if np.any(gauss):raise ArithmeticError("Exact Gauss charge failed")
    t=np.zeros((n,n));w=np.zeros((n,n))
    for i,(x,ell) in enumerate(labels):
        target=((x+1)%L,ell+q if x==L-1 else ell)
        if target in lookup:t[lookup[target],i]=1
        target=(x,ell+1)
        if target in lookup:w[lookup[target],i]=1
    e_diag=electric*np.sum(flux**2,axis=1)/2
    h=np.diag(e_diag)-hopping*(t+t.T)-loop*(w+w.T)/2
    projectors=[np.diag([int(x==y) for x,ell in labels]).astype(float) for y in range(L)]
    initial=np.zeros(n);initial[lookup[0,0]]=1
    return {"H":h,"labels":labels,"flux":flux,"Gauss":gauss,"T":t,"W":w,
            "electric_diagonal":e_diag,"position_projectors":projectors,
            "initial_state":initial,"q":q,"K":K,"L":L,
            "couplings":{"electric":electric,"loop":loop,"hopping":hopping}}

def operator_receipt(model):
    q=model["q"];t=model["T"];w=model["W"];n=len(t);L=model["L"]
    worst=0.
    for e in range(L):
        ee=model["flux"][:,e]
        # [E_e,T]=q T Q_e, including every finite-cutoff boundary vector.
        defect=(ee[:,None]-ee[None,:])*t-q*t@model["position_projectors"][e]
        worst=max(worst,float(np.max(abs(defect))))
        defect=(ee[:,None]-ee[None,:])*w-w
        worst=max(worst,float(np.max(abs(defect))))
    loop_power=np.linalg.matrix_power(w if q>=0 else w.T,abs(q))
    return {"physical_dimension":n,"Gauss_defect_max":int(abs(model["Gauss"]).max()),
       "Hermiticity_error":float(np.max(abs(model["H"]-model["H"].T))),
       "electric_shift_Ward_error":worst,
       "T_power_L_equals_W_power_q_error":float(np.max(abs(np.linalg.matrix_power(t,L)-loop_power))),
       "T_unitarity_defect_rank":int(np.linalg.matrix_rank(np.eye(n)-t.T@t)),
       "W_unitarity_defect_rank":int(np.linalg.matrix_rank(np.eye(n)-w.T@w)),
       "T_unitarity_defect_norm":float(np.max(abs(np.linalg.eigvalsh(np.eye(n)-t.T@t)))),
       "W_unitarity_defect_norm":float(np.max(abs(np.linalg.eigvalsh(np.eye(n)-w.T@w)))),
       "all_boundary_vectors_retained":True}

def complete_binary_records(h,state,projector,step_times):
    """Keep all 2^m histories, no postselection; times are interval durations."""
    ev,vec=eigh(h);rho=np.outer(state,state.conj())
    # Evolving pure branch vectors gives exact multi-round Lueders probabilities.
    branches={"":np.asarray(state,complex)}
    ps=[np.eye(len(h))-projector,projector]
    for dt in step_times:
        u=(vec*np.exp(-1j*ev*dt))@vec.conj().T
        branches={hist+str(b):ps[b]@u@v for hist,v in branches.items() for b in range(2)}
    records=[{"history":[int(x) for x in k],"probability":float(np.vdot(v,v).real)}
             for k,v in branches.items()]
    return {"step_times":list(step_times),"records":records,
            "normalization_error":float(abs(sum(r["probability"] for r in records)-1)),
            "minimum_probability":float(min(r["probability"] for r in records))}

def spectral_receipt(model):
    ev,vec=eigh(model["H"]);ground=vec[:,0]
    obs=model["position_projectors"][0]
    amplitudes=vec.conj().T@obs@ground;weights=abs(amplitudes)**2
    mean=float(np.vdot(ground,obs@ground).real);variance=mean-mean**2
    bound_mask=np.array([abs(ell)==model["K"] for x,ell in model["labels"]])
    boundary=float(np.sum(abs(ground[bound_mask])**2))
    active=[{"energy_above_ground":float(ev[j]-ev[0]),"weight":float(weights[j])}
            for j in range(1,len(ev)) if weights[j]>1e-10]
    lowest=active[0] if active else None
    # Full finite spectral measure is exported separately, including zero weights.
    record={"ground_energy":float(ev[0]),"gap":float(ev[1]-ev[0]),
       "position_zero_mean":mean,"elastic_spectral_weight":float(weights[0]),
       "inelastic_spectral_weight":float(weights[1:].sum()),
       "projector_variance":variance,
       "spectral_sum_rule_error":float(abs(weights.sum()-mean)),
       "variance_sum_rule_error":float(abs(weights[1:].sum()-variance)),
       "ground_flux_boundary_probability":boundary,
       "first_nonzero_record_pole":lowest,
       "record_pole_reporting_threshold":1e-10,
       "number_of_nonzero_inelastic_poles_above_1e_minus_10":len(active),
       "lowest_eight_energy_levels":ev[:8].tolist()}
    return record,ev,vec,weights

def dyson_bound(K,total_time,max_charge=6,initial_flux_support=0,hopping=1.,loop=1.):
    K=_integer(K,"K",0)
    initial_flux_support=_integer(initial_flux_support,"initial_flux_support",0)
    max_charge=_integer(max_charge,"max_charge")
    total_time=_real(total_time,"total_time")
    hopping=_real(hopping,"hopping");loop=_real(loop,"loop")
    if initial_flux_support>K or total_time<0:
        raise ValueError("The initial flux support must lie inside the cutoff and time must be nonnegative")
    r=int(np.ceil((K-initial_flux_support+1)/max(1,abs(max_charge))))
    b=2*abs(hopping)+abs(loop);x=b*total_time
    tail=float(np.exp(x)*gammainc(r,x))
    return {"K":K,"total_time":total_time,"initial_flux_support":initial_flux_support,
       "maximum_flux_step":max(1,abs(max_charge)),"minimum_insertions_to_leave":r,
       "bounded_interaction_norm_bound":b,"exponential_tail":tail,
       "state_vector_error_bound":min(2.,2*tail),
       "complete_record_total_variation_error_bound":min(1.,2*tail),
       "applies_to":"localized initial state and any complete position-record instrument whose total unitary time is T",
       "does_not_bound":"ground-state record errors unless initial-state approximation is included"}

def infinite_ground_bracket(model,finite_ground):
    # All discarded ell have |ell|>=K+1.  Minimize the diagonal quadratic exactly
    # over the two exterior half-lines and every position.
    q,K,L=model["q"],model["K"],model["L"]
    candidates=[]
    for x in range(L):
        optimum=-q*x/L
        for sign in (-1,1):
            near=int(np.floor(optimum))
            allowed=[ell for ell in range(near-2,near+3) if sign*ell>=K+1]
            allowed.append(sign*(K+1))
            for ell in allowed:
                candidates.append(model["couplings"]["electric"]/2*
                                  sum((ell+q*(e<x))**2 for e in range(L)))
    b=2*abs(model["couplings"]["hopping"])+abs(model["couplings"]["loop"])
    qlower=min(candidates)-b
    lower=(finite_ground+qlower-np.hypot(qlower-finite_ground,2*b))/2
    return {"infinite_ground_lower_bound":float(lower),
            "infinite_ground_upper_bound":float(finite_ground),
            "discarded_Q_block_lower_bound":float(qlower),
            "offdiagonal_coupling_norm_bound":b}

def main(repo=ROOT,out=None):
    src=load(repo);out = Path(out) if out is not None else None
    if out is not None:
        out.mkdir(parents=True,exist_ok=True)
    cutoffs=(4,8,12,24,48);rows=[];spectra={};matrices={}
    for block in src["blocks"]:
        q=block["q6Y"];d=block["dimension"];block_rows=[]
        for K in cutoffs:
            model=build_flux_ring(q,K);ops=operator_receipt(model)
            if max(ops["Gauss_defect_max"],ops["Hermiticity_error"],
                   ops["electric_shift_Ward_error"],ops["T_power_L_equals_W_power_q_error"])!=0:
                raise ArithmeticError("Finite physical carrier algebra")
            spec,ev,vec,weights=spectral_receipt(model)
            if max(spec["spectral_sum_rule_error"],spec["variance_sum_rule_error"])>1e-10:
                raise ArithmeticError("Complete spectral measure")
            block_rows.append({"K":K,"operators":ops,"response":spec,
                "infinite_ground_bracket":infinite_ground_bracket(model,ev[0])})
            spectra[f"q{q}_K{K}_levels"]=ev
            spectra[f"q{q}_K{K}_projector_weights"]=weights
            if K==24:
                records=complete_binary_records(model["H"],model["initial_state"],
                                      model["position_projectors"][0],(.125,.125,.125,.125))
                if records["normalization_error"]>1e-10:raise ArithmeticError("Record normalization")
                block_rows[-1]["complete_four_round_binary_records"]=records
                if q==1:
                    matrices={"H":model["H"],"rho":np.outer(model["initial_state"],model["initial_state"]),
                       "psi":model["initial_state"],"position":np.array([x for x,ell in model["labels"]]),
                       "loop_flux":np.array([ell for x,ell in model["labels"]]),"edge_flux":model["flux"],
                       "P_plus":model["position_projectors"][0],"P_minus":np.eye(len(ev))-model["position_projectors"][0],
                       "ground_state":vec[:,0],"ground_rho":np.outer(vec[:,0],vec[:,0]),
                       "T":model["T"],"W":model["W"]}
        reference=block_rows[-1]["response"]
        for row in block_rows:
            r=row["response"]
            row["difference_from_K48"]={
              "ground_energy":r["ground_energy"]-reference["ground_energy"],
              "inelastic_position_spectral_weight":r["inelastic_spectral_weight"]-
                                                   reference["inelastic_spectral_weight"]}
        # A genuine fluctuating-loop response: all else unchanged, compare
        # magnetic/holonomy coupling one and zero.  E is a noncommuting
        # dynamical operator in both.  This is a fixed declared control.
        off=build_flux_ring(q,24,loop=0.)
        offspec,*_=spectral_receipt(off)
        on=block_rows[3]["response"]
        rows.append({"q6Y":q,"source_internal_multiplicity":d,
           "scope":"U(1) Gauss sector; nonabelian components are copies, not a full nonabelian Gauss reduction",
           "cutoff_sequence":block_rows,
           "loop_interaction_control":{"loop_on":1.,"loop_off":0.,
               "ground_energy_shift":on["ground_energy"]-offspec["ground_energy"],
               "inelastic_position_weight_shift":on["inelastic_spectral_weight"]-offspec["inelastic_spectral_weight"],
               "off_response":offspec}})
    if out is not None:
        np.savez_compressed(out/"complete_flux_spectra.npz",**spectra)
    if out is not None:
        np.savez_compressed(out/"portable_quantum_flux_probe.npz",**matrices)
    result={"schema":"oph.interacting_extension.quantum_flux_response.v1",
      "source":src["provenance"],
      "numeric_status":"float64 diagnostics and evaluations of analytic bounds; not interval-certified arithmetic",
      "control_scope":"One four-edge closed graph, native U(1) group-coordinate restriction with its own Haar rotor, one mobile charge plus static conjugate charge",
      "restriction_caveat":"This U(1) control is not a reducing Hilbert subspace of the full G rotor; restricted nonabelian components are charge copies",
      "exact_Gauss_solution":"E_e=ell+q*1[e<x]; G_v=q delta(v,x)-q delta(v,0)+E_v-E_(v-1)=0",
      "hamiltonian":"(1/2) sum_e E_e^2 -(T+T^dag) -(W+W^dag)/2",
      "couplings":{"electric":1.,"hopping":1.,"loop":1.,
                   "selection":"declared dimensionless controls; not source-selected or empirically fit",
                   "magnetic_choice":"primitive U(1) loop W from native q=1 block; not a uniquely selected full M15 trace action"},
      "all_native_charge_blocks":rows,
      "uniform_complete_record_cutoff_bounds":[dyson_bound(K,.5) for K in cutoffs],
      "all_finite_spectral_outcomes_file":"complete_flux_spectra.npz",
      "portable_probe_file":"portable_quantum_flux_probe.npz",
      "source_consumption":{"source_charge_blocks":True,
          "operator_valued_gauge_flux_retained":True,"Gauss_constraints_exact":True,
          "all_finite_cutoff_complements_retained":True,
          "full_nonabelian_interactions_executed":False,
          "three_spatial_dimensions_executed_in_this_control":False,
          "source_selected_joint_vacuum":False,"interacting_chiral_continuum":False,
          "measured_masses_used":False}}
    if out is not None:
        (out/"quantum_flux_results.json").write_text(json.dumps(result,indent=2)+"\n")
    if out is not None:
        print(json.dumps({"responses":[{"q":r["q6Y"],
            "K48_ground":r["cutoff_sequence"][-1]["response"]["ground_energy"],
            "K48_inelastic_weight":r["cutoff_sequence"][-1]["response"]["inelastic_spectral_weight"],
            "loop_on_minus_off":r["loop_interaction_control"]["ground_energy_shift"]}
             for r in rows],"uniform_K48_record_bound":result["uniform_complete_record_cutoff_bounds"][-1]},indent=2))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--repo",type=Path,default=ROOT)
    p.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args();main(a.repo,a.out)
