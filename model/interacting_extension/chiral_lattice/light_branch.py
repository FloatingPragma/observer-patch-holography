#!/usr/bin/env python3
"""A growing THREE-spatial-dimensional source-convention overlap control.

The flat gauge fiber is explicitly a free control of the completed rotor carrier,
not a replacement for its fluctuating operator-valued links.  All spin, momentum
and native internal multiplicities are counted.  The portable 4x4 probe is a
single momentum block; complete binary record probabilities are computed.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from .source_data import load, ROOT

def spatial_block(momentum,spacing,gamma):
    if isinstance(spacing,bool) or not np.isscalar(spacing) or not np.isfinite(spacing) or spacing<=0:
        raise ValueError("spacing must be finite and positive")
    momentum=np.asarray(momentum,float)
    if momentum.shape!=(3,) or not np.isfinite(momentum).all():
        raise ValueError("momentum must have three finite real components")
    gs,g5=gamma
    beta=np.asarray(gs[3],complex)
    k=momentum*spacing
    s=np.sin(k);u=1-np.cos(k);b=float(np.sum(u)-1)
    rad=float(np.sqrt(b*b+s@s))
    v=(b*np.eye(4)+1j*sum(s[i]*gs[i] for i in range(3)))/rad
    # Rationalize 1+B/R near the light cone, where direct subtraction loses
    # significant digits. This is algebraically identical away from k=0.
    one_plus=float(s@s/(rad*(rad-b))) if b<0 else 1+b/rad
    d=(one_plus*np.eye(4)+1j*sum(s[i]*gs[i] for i in range(3))/rad)/spacing
    h=beta@d
    q5=np.asarray(g5)@(np.eye(4)-v)/2
    e=float(np.sqrt(max(0,2*one_plus))/spacing)
    return h,q5,{"B":b,"R":rad,"one_plus_B_over_R":one_plus,"energy":e,"k":k.tolist(),"V":v}

def binary_probe(h,observable,times):
    """Full two-round Lueders law in the normalized negative-energy state."""
    ev,vec=np.linalg.eigh(h)
    if np.count_nonzero(ev < -1e-10)!=2:raise ValueError("Need nonzero +/-E block")
    rho=vec[:,:2]@vec[:,:2].conj().T/2
    projectors=[(np.eye(len(h))+s*observable)/2 for s in (-1,1)]
    records=[]
    for t in times:
        unitary=expm(-1j*h*t)
        probs=[]
        for i,j in itertools.product(range(2),repeat=2):
            branch=projectors[j]@unitary@projectors[i]
            probs.append(float(np.trace(branch@rho@branch.conj().T).real))
        correlation=sum((-1 if i==0 else 1)*(-1 if j==0 else 1)*probs[2*i+j]
                        for i,j in itertools.product(range(2),repeat=2))
        records.append({"time":float(t),"histories":[[-1,-1],[-1,1],[1,-1],[1,1]],
                        "probabilities":probs,"correlation":float(correlation)})
    return rho,projectors,records

def position_response(length,spacing,gamma,time):
    """Exact complete position law from all L^3 momentum blocks, spin I/4."""
    us=np.empty((length,length,length,4,4),complex)
    for idx in itertools.product(range(length),repeat=3):
        k=2*np.pi*np.asarray(idx)/length
        h,_,_=spatial_block(k/spacing,spacing,gamma)
        us[idx]=expm(-1j*h*time)
    amplitudes=np.fft.ifftn(us,axes=(0,1,2))
    probabilities=np.sum(np.abs(amplitudes)**2,axis=(-2,-1))/4
    return {"length":length,"spacing":spacing,"time":time,
            "initial_position":[0,0,0],"initial_spin_state":"I4/4",
            "probabilities":[{"position":list(idx),"probability":float(probabilities[idx])}
                             for idx in itertools.product(range(length),repeat=3)],
            "total_probability":float(probabilities.sum()),
            "origin_probability":float(probabilities[(0,0,0)]),
            "away_from_origin_probability":float(1-probabilities[(0,0,0)]),
            "all_native_internal_components_have_this_position_law":True}

def main(repo=ROOT,out=None):
    src=load(repo);gamma=src["gamma"];gs,g5=gamma;beta=np.asarray(gs[3])
    out = Path(out) if out is not None else None
    if out is not None:
        out.mkdir(parents=True,exist_ok=True)
    sequences=[];allrecords=[]
    for length in (32,64,128,256):
        spacing=2*np.pi/length;p=np.ones(3)
        h,q5,d=spatial_block(p,spacing,gamma);e=d["energy"]
        rho,projectors,records=binary_probe(h,beta,np.linspace(0,1,9))
        bcharge=d["one_plus_B_over_R"]/(spacing*e)
        residue=1-bcharge*bcharge
        errors={
          "Hermiticity":float(np.linalg.norm(h-h.conj().T,2)),
          "modified_axial_charge_Hermiticity":float(np.linalg.norm(q5-q5.conj().T,2)),
          "modified_axial_charge_conservation":float(np.linalg.norm(h@q5-q5@h,2)),
          "modified_axial_charge_square":float(np.linalg.norm(q5@q5-(np.eye(4)-spacing**2*h@h/4),2)),
          "dispersion":float(np.max(np.abs(np.linalg.eigvalsh(h)-np.array([-e,-e,e,e])))),
          "record_normalization":max(abs(sum(x["probabilities"])-1) for x in records),
          "record_correlation_formula":max(abs(x["correlation"]-
                   (bcharge*bcharge+residue*np.cos(2*e*x["time"]))) for x in records),
          "continuum_operator_error":float(np.linalg.norm(h-sum(1j*beta@gs[i]*p[i] for i in range(3)),2)),
        }
        bound=spacing*np.dot(p,p)
        sequences.append({"L":length,"box_length":2*np.pi,"spacing":spacing,
             "spatial_dimension":3,"physical_momentum":p.tolist(),"energy":e,
             "continuum_energy":float(np.linalg.norm(p)),"relative_energy_error":float(e/np.linalg.norm(p)-1),
             "beta_inelastic_spectral_residue":float(residue),
             "projector_inelastic_spectral_residue":float(residue/4),
             "projector_elastic_weight":float((1-bcharge)**2/4),
             "measured_frequency":float(2*e),"kernel_gap_in_dimensionless_units":d["R"],
             "low_momentum_bound_applies":bool(spacing*np.linalg.norm(p)<=.25),
             "continuum_operator_error_bound_a_p_squared":float(bound),"errors":errors})
        allrecords.append({"L":length,"records":records})
        if not all(np.isfinite(v) for v in errors.values()):
            raise ArithmeticError("Nonfinite free-field diagnostic")
        if min(z for x in records for z in x["probabilities"]) < -1e-13:
            raise ArithmeticError("Negative complete detector probability")
        if max(v for k,v in errors.items() if k!="continuum_operator_error") >= 1e-10:
            raise ArithmeticError("Free-field operator or complete-record identity failed")
        if spacing*np.linalg.norm(p)<=.25 and errors["continuum_operator_error"]>bound+1e-12:
            raise ArithmeticError("Free continuum comparison exceeded its declared band bound")
    # Complete low-resolution Brillouin zone, including every zero and doubler
    # corner.  Native internal multiplicity is fifteen, not a selected charge.
    length=8;spacing=2*np.pi/length
    momenta=[];energies=[];gaps=[];charge_extrema=[];corner=[]
    for idx in itertools.product(range(length),repeat=3):
        k=2*np.pi*np.asarray(idx)/length
        h,q5,d=spatial_block(k/spacing,spacing,gamma)
        ev=np.linalg.eigvalsh(h)
        momenta.append(idx);energies.append(ev);gaps.append(d["R"])
        charge_extrema.append(np.linalg.eigvalsh(q5))
        if all(j in (0,length//2) for j in idx):
            corner.append({"momentum_indices":list(idx),"energy":d["energy"]})
    energies=np.asarray(energies)
    if out is not None:
        np.savez_compressed(out/"complete_free_brillouin_zone.npz",
        momentum_indices=np.asarray(momenta),four_spin_eigenvalues=energies,
        modified_axial_eigenvalues=np.asarray(charge_extrema),
        native_internal_charge=np.asarray(src["charges"]),spacing=spacing)
    h,q5,d=spatial_block(np.ones(3),2*np.pi/64,gamma)
    rho,projectors,_=binary_probe(h,beta,[0])
    if out is not None:
        np.savez_compressed(out/"portable_light_probe.npz",H=h,rho=rho,observable=beta,
                        P_minus=projectors[0],P_plus=projectors[1],q5=q5,
                        momentum=np.ones(3),spacing=2*np.pi/64)
    result={"schema":"oph.interacting_extension.light_branch.v1","source":src["provenance"],
      "control_scope":"Flat gauge fiber, vectorlike three-spatial-dimensional Hamiltonian; no interacting-continuum claim",
      "numeric_status":"float64 diagnostics, analytic identities separate; no interval-certified measured values",
      "formula":{"A":"sum_i(1-cos k_i)-1 + i sum_i gamma_i sin k_i",
                 "H":"gamma4 (I+A/sqrt(A^dag A))/a",
                 "free_kernel_identity":"R^2=1+2 sum_(i<j) (1-cos k_i)(1-cos k_j)",
                 "q5":"gamma5 (I-A/sqrt(A^dag A))/2",
                 "record_correlation":"b_h^2 + (1-b_h^2) cos(2 E t)"},
      "growing_3d_lattice_sequence":sequences,
      "complete_two_round_binary_records":allrecords,
      "complete_brillouin_zone":{"L":length,"momentum_points":length**3,
          "spin_dimension_per_point":4,"native_internal_dimension_per_point":15,
          "total_one_body_dimension":length**3*4*15,
          "total_zero_modes":int(np.count_nonzero(abs(energies)<1e-10)*15),
          "minimum_kernel_gap":float(min(gaps)),
          "all_eight_corners":corner,"file":"complete_free_brillouin_zone.npz",
          "all_complements_retained":True},
      "complete_position_response":position_response(4,1.,gamma,.5),
      "control_parameters":{"rho":1,"box_length":2*np.pi,"momentum":[1,1,1],
                            "units":"dimensionless simulator controls; no derived GeV scale",
                            "reference_state":"negative-energy projector/2 for the 4-spin block, explicitly chosen control"},
      "source_consumption":{"native_representation_and_multiplicities":True,
                           "source_gamma_convention":True,"source_Wilson_rho_one":True,
                           "native_60_carrier_used_without_enlargement":False,
                           "quantum_link_fluctuations_in_this_control":False,
                           "source_selected_joint_vacuum":False,
                           "interacting_chiral_continuum":False,
                           "measured_masses_used":False}}
    if out is not None:
        (out/"light_branch_results.json").write_text(json.dumps(result,indent=2)+"\n")
    if out is not None:
        print(json.dumps({"sequence":[{k:r[k] for k in ("L","energy","beta_inelastic_spectral_residue")}
                        for r in sequences],
                        "complete_BZ":result["complete_brillouin_zone"],
                        "position_leaves_origin":result["complete_position_response"]["away_from_origin_probability"]},indent=2))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--repo",type=Path,default=ROOT)
    p.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args();main(a.repo,a.out)
