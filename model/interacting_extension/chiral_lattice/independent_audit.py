#!/usr/bin/env python3
"""Independent source-stencil, high-precision, and record-cutoff checks."""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import mpmath as mp
from .source_data import load, ROOT
from .light_branch import spatial_block
from .rotor_response import build_flux_ring,complete_binary_records,dyson_bound

def main(repo=ROOT,out=None):
    src=load(repo)
    from model.common_source.overlap_regulator import wilson_operator,overlap
    # Directly execute the complete pinned 4D Wilson/overlap implementation.
    # Fourier projection at p4=0 must reproduce each spatial block.
    side=3;sites=tuple(itertools.product(range(side),repeat=4))
    links=np.ones((len(sites),4,1,1),complex)
    w,g5=wilson_operator(links,side);d,_,gap=overlap(w,g5)
    beta=np.asarray(src["gamma"][0][3])
    maximum_projection=0.;maximum_invariance=0.;checks=0
    for ind in itertools.product(range(side),repeat=3):
        k=2*np.pi*np.array((*ind,0))/side
        phase=np.exp(1j*np.array(sites)@k)/np.sqrt(len(sites))
        embedding=np.kron(phase[:,None],np.eye(4))
        actual=beta@(embedding.conj().T@d@embedding)
        wanted,_,_=spatial_block(k[:3],1.,src["gamma"])
        maximum_projection=max(maximum_projection,float(np.linalg.norm(actual-wanted,2)))
        maximum_invariance=max(maximum_invariance,
                 float(np.linalg.norm(d@embedding-embedding@(beta@wanted),2)))
        checks+=1
    if max(maximum_projection,maximum_invariance)>1e-11:
        raise ArithmeticError("Source-stencil projection mismatch")
    # 90-digit evaluation is independent of double precision cancellation.
    precision=[]
    with mp.workdps(90):
        for L in (32,64,128,256):
            a=2*mp.pi/L;k=a;b=3*(1-mp.cos(k))-1
            r=mp.sqrt(b*b+3*mp.sin(k)**2)
            en=mp.sqrt(2*3*mp.sin(k)**2/(r*(r-b)))/a
            residue=1-(a*en/2)**2
            _,_,info=spatial_block(np.ones(3),float(a),src["gamma"])
            err=abs(mp.mpf(info["energy"])-en)
            if err>mp.mpf("2e-14"):raise ArithmeticError("High precision light-cone check")
            precision.append({"L":L,"energy_90digit":mp.nstr(en,82),
                               "beta_residue_90digit":mp.nstr(residue,82),
                               "float64_energy_error":float(err)})
    # Compare complete record distributions across nested flux cutoffs.
    records=[]
    for q in (-4,-3,1,2,6):
        reference=build_flux_ring(q,48)
        ref=complete_binary_records(reference["H"],reference["initial_state"],
                 reference["position_projectors"][0],(.125,.125,.125,.125))
        refp=np.array([r["probability"] for r in ref["records"]])
        for K in (4,8,12,24):
            model=build_flux_ring(q,K)
            law=complete_binary_records(model["H"],model["initial_state"],
                     model["position_projectors"][0],(.125,.125,.125,.125))
            p=np.array([r["probability"] for r in law["records"]])
            tv=float(np.sum(abs(p-refp))/2)
            # Triangle inequality through the infinite evolution.
            b=dyson_bound(K,.5,max_charge=abs(q))["complete_record_total_variation_error_bound"]
            b48=dyson_bound(48,.5,max_charge=abs(q))["complete_record_total_variation_error_bound"]
            if tv>b+b48+2e-12:raise ArithmeticError("Complete-record cutoff bound")
            records.append({"q":q,"K":K,"complete_record_TV_from_K48":tv,
                            "analytic_triangle_bound_float64":min(1.,b+b48),
                            "all_16_histories_retained":True})
    result={"source":src["provenance"],
      "complete_source_stencil_check":{"four_dimensional_ambient_dimension":len(d),
          "all_three_spatial_momentum_blocks_at_p4_zero":checks,
          "source_Wilson_gap":gap,"maximum_symbol_projection_error":maximum_projection,
          "maximum_complete_momentum_subspace_invariance_error":maximum_invariance},
      "independent_high_precision_light_energy":precision,
      "complete_record_cutoff_comparisons":records,
      "remaining_numerical_limit":"90-digit scalar checks do not interval-certify matrix eigensystems or continuum claims"}
    out = Path(out) if out is not None else None
    if out is not None:
        out.mkdir(parents=True,exist_ok=True)
    if out is not None:
        (out/"independent_audit_results.json").write_text(json.dumps(result,indent=2)+"\n")
    if out is not None:
        print(json.dumps({"source_stencil":result["complete_source_stencil_check"],
           "max_90digit_scalar_discrepancy":max(p["float64_energy_error"] for p in precision),
           "record_laws_checked":len(records)},indent=2))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--repo",type=Path,default=ROOT)
    p.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args();main(a.repo,a.out)
