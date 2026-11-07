"""Replay the complete local chiral-electric construction and its finite parent."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import covariant_kinetic
from .._source import PIN, ROOT, verify_sources

HERE=Path(__file__).resolve().parent
OUTPUT=HERE/"results.json"


def source_bindings():
    bindings=verify_sources(ROOT,("model/common_source/overlap_regulator.py",))
    paths=list(HERE.glob("*.py"))+list(HERE.glob("*.md"))
    paths += [HERE.parent/name for name in ("__init__.py","_source.py","source_manifest.json")]
    paths += [HERE.parent/"dependencies"/name for name in ("__init__.py","chiral_patch.py")]
    for path in paths:
        bindings[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(bindings.items()))


def report():
    before=source_bindings()
    data=covariant_kinetic.main()
    d=data["differential"]
    if (d["ambient_dimension"],d["chiral_rank"])!=(64,32):
        raise ArithmeticError("Changed complete source-projector domain")
    if not .71366<d["bare_electric_leakage_norm"]<.71369:
        raise ArithmeticError("Changed local electric leakage diagnostic")
    if not .01221<d["Born_Huang_minimum_eigenvalue"]<.01223:
        raise ArithmeticError("Changed geometric kinetic diagnostic")
    if not .50931<d["Born_Huang_maximum_eigenvalue"]<.50934:
        raise ArithmeticError("Changed geometric kinetic diagnostic")
    for run in data["complete_penalty"]["runs"]:
        covariant_kinetic.small(run["eigensystem_residual"],"Full parent eigensystem",1e-8)
    if source_bindings()!=before:
        raise ArithmeticError("Source changed during chiral-electric replay")
    return {
        "schema":"oph.interacting_extension.chiral_electric_certificate.v1",
        "source_commit":PIN,"source_sha256":before,
        "source_domain":{"charge":6,"Euclidean_sites":16,"spin_site_dimension":64,"chiral_rank":32},
        "declared_controls":{"toron_per_link":[.17,-.11,.07,.03],"holonomy_halfwidth":.16,
            "kinetic_coefficient":.01,"boundary":"Dirichlet","grid_sizes":[5,9,17],
            "penalty_factors":[10,40,160]},
        "identities":{
            "connection":"C=[P,P']; P'+[C,P]=0",
            "parallel_transport_generator":"-C=[P',P]",
            "geometric_term":"Phi=V'^dagger(I-P)V'=V'^dagger V'+A^2>=0",
            "compression":"V^dagger(-d^2)V=-(d+A)^2+Phi",
            "finite_parent_bound":"0<=a_j-lambda_j<=norm(QKV)^2/(Lambda-a_k), j<=k, Lambda>a_k",
        },
        "complete_parent":{"dimension":576,"compressed_dimension":288,"eigenvalues_compared":8,
            "all_complementary_states_retained":True,"Schur_checks_passed":True},
        "checks":{"differential_identities":True,"independent_centered_derivative":True,
            "full_compression_assembly":True,"independent_frame_covariance":True,
            "bare_identical_link_frame_covariance":False,
            "raw_float_residuals_stored_in_equality_receipt":False},
        "scope":{"local_Euclidean_projector_patch":True,"physical_axial_conservation_proved":False,
            "finite_grid_geometric_difference_positive_semidefinite":False,
            "microscopic_kinetic_law_selected":False,"global_chiral_measure_or_continuum":False,
            "physical_masses_or_clock":False,"interval_arithmetic_certificate":False},
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",nargs="?",const=OUTPUT,type=Path)
    mode.add_argument("--check",nargs="?",const=OUTPUT,type=Path)
    args=parser.parse_args();result=report()
    if args.write is not None:
        args.write.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    elif json.loads(args.check.read_text())!=result:
        raise ValueError("Chiral-electric certificate differs from fresh replay")
    print("Local chiral-electric identities, full complement and finite parent verified.")


if __name__=="__main__":
    main()
