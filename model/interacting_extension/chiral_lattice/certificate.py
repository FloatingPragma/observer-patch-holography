"""Replay finite free propagation and abelian flux controls without rounding receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from . import independent_audit, light_branch, rotor_response
from .source_data import FILES, ROOT
from .._source import PIN, verify_sources

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results.json"


def source_bindings():
    bindings = verify_sources(ROOT, FILES)
    paths = list(HERE.glob("*.py")) + list(HERE.glob("*.md"))
    paths += [HERE.parent / name for name in ("__init__.py", "_source.py", "source_manifest.json")]
    for path in paths:
        bindings[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(bindings.items()))


def report():
    before = source_bindings()
    free = light_branch.main()
    flux = rotor_response.main()
    independent = independent_audit.main()
    bz = free["complete_brillouin_zone"]
    if (bz["total_one_body_dimension"], bz["total_zero_modes"]) != (30720, 60):
        raise ArithmeticError("Incomplete spin, momentum or native internal complement")
    reference_energies = (1.7448298797, 1.7354582753, 1.7329158901, 1.7322679027)
    reference_weights = (.9706569017, .9927428157, .9981910158, .9995480921)
    for row, energy, weight in zip(free["growing_3d_lattice_sequence"], reference_energies, reference_weights, strict=True):
        if abs(row["energy"]-energy)>1e-9 or abs(row["beta_inelastic_spectral_residue"]-weight)>1e-9:
            raise ArithmeticError("Free control differs from its declared numerical check window")
    if abs(free["complete_position_response"]["total_probability"]-1)>1e-11:
        raise ArithmeticError("The complete position law lost probability")
    for block in flux["all_native_charge_blocks"]:
        for row in block["cutoff_sequence"]:
            bracket=row["infinite_ground_bracket"]
            if not np.isfinite(list(bracket.values())).all():
                raise ArithmeticError("Nonfinite diagnostic ground-energy bracket")
            if bracket["infinite_ground_lower_bound"]>bracket["infinite_ground_upper_bound"]:
                raise ArithmeticError("Reversed diagnostic ground-energy bracket")
    charged=next(b for b in flux["all_native_charge_blocks"] if b["q6Y"]==1)
    if abs(charged["cutoff_sequence"][-1]["response"]["inelastic_spectral_weight"]-.2324598386713)>1e-9:
        raise ArithmeticError("Interacting flux response differs from its numerical check window")
    if abs(charged["loop_interaction_control"]["ground_energy_shift"]+.3768864726971)>1e-9:
        raise ArithmeticError("The loop interaction correction is not reproduced")
    stencil=independent["complete_source_stencil_check"]
    if stencil["all_three_spatial_momentum_blocks_at_p4_zero"] != 27:
        raise ArithmeticError("Incomplete direct source-stencil comparison")
    if source_bindings()!=before:
        raise ArithmeticError("Source changed during the finite replay")
    return {
        "schema": "oph.interacting_extension.propagation_certificate.v1",
        "source_commit": PIN, "source_sha256": before,
        "free_control": {
            "lattice_sides": [32,64,128,256], "fixed_momentum": [1,1,1],
            "spatial_dimensions": 3, "all_component_count_L8": 30720,
            "zero_momentum_spin_internal_multiplicity": 60,
            "kernel_gap_identity": "R^2=1+2 sum_(i<j)(1-cos k_i)(1-cos k_j)>=1",
            "dispersion": "E^2=2(1+B/R)/a^2",
            "complete_binary_law": "C(t)=b_h^2+(1-b_h^2)cos(2Et)",
            "numerical_energy_and_weight_windows": "absolute error <= 1e-9 against the declared controls",
            "direct_source_momentum_blocks_checked": 27,
            "direct_source_operator_tolerance": "1e-11",
            "all_free_operator_and_record_checks_passed": True,
        },
        "flux_control": {
            "native_charges": [-4,-3,1,2,6], "cutoffs": [4,8,12,24,48],
            "dimension_per_charge": "4(2K+1)", "record_histories_per_run": 16,
            "independent_cutoff_record_comparisons": 20,
            "hamiltonian": "sum_e E_e^2/2-(T+T^dagger)-(W+W^dagger)/2",
            "record_bound_K48_T_half": "2 sum_(n>=9) (3/2)^n/n!",
            "all_Gauss_and_electric_shift_identities_passed": True,
            "all_complete_records_and_cutoff_checks_passed": True,
            "quantum_loop_correction_checked": True,
            "response_diagnostic_window_tolerance": "1e-9",
            "cutoff_boundary_vectors_retained": True,
        },
        "scope": {
            "rotor_completion_and_couplings_are_declared": True,
            "free_flat_fiber_is_a_normalizable_rotor_state": False,
            "abelian_flux_and_free_spatial_controls_are_one_interacting_theory": False,
            "source_selected_state_clock_or_P_N_map": False,
            "interacting_chiral_continuum_or_physical_masses": False,
            "floating_diagnostics_are_interval_certificates": False,
        },
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
        raise ValueError("Finite propagation certificate differs from fresh replay")
    print("Finite free propagation, complete abelian records and cutoff bounds verified.")


if __name__=="__main__":
    main()
