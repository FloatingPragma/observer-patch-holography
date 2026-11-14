"""Regenerate the actual source, the limit control, and every retained result."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from dynamical_field import HERE, ROOT, execute


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def contract(result):
    return {
        "schema": "oph.c10.dynamical_field_contract.v1",
        "source": {"baseline": "C10_interacting_allnumber384_Z6",
                   "H_sha256_le_float64": result["baseline"]["H_sha256_le_float64"],
                   "family": "BoundarySource(M) + ModularField(kappa=.37, face_coefficient=.61), all 64M states",
                   "member_parameters": "integer M>=6; mass .23; no coefficient rescaling",
                   "source_code": "dynamical_field.py",
                   "ancestry": "same actual half-hop, CAR, periodic electric, magnetic and complete MIX laws"},
        "preparation": {"baseline": "recomputed M6 global interacting ground projector and all four fair-reservoir outcomes",
                        "family_map": "source tree chord3, k=centered(Fmod6[3]), retain mask/k and reservoir label",
                        "support": "all 384 baseline labels have max_h abs(integer F_h)<=5",
                        "member_ground_preparation": False,
                        "postselection": False,
                        "outside_weight_removed": 0,
                        "reference": "retain the full purification/record labels; theorem allows arbitrary extra quantum references"},
        "clock": {"contract": "../clock/CLOCK_CONTRACT.json",
                  "unitary": "exp(-i t H_M)", "T": "37/100",
                  "native_attempt_rate": "4/5", "idle_mark_rate": "4/5",
                  "law": "MIX", "every_no_current_attempt_retained": True,
                  "Poisson_tail_discarded": False,
                  "native_H_compiler_supplied": False,
                  "physical_seconds_identified": False},
        "limit": {"carrier": "integer Gauss fiber n+D F+(3-N)delta0=G, G=(111000), C64 tensor l2(Z)",
                  "operator": "Hmass + .37/2 sum F_h^2 + .61(2I-W-W*) + original integer-record CAR halfhops",
                  "self_adjoint_domain": "maximal diagonal electric-energy domain; bounded perturbation norm <=361/50",
                  "observable_domain": "complete marked histories and retained embedded quantum terminal state; all bounded terminal tests",
                  "radius": 5,
                  "bound": "min(1,2 sum_{j>=floor((M-1)/2)-5+1} (14837/5000)^j/j!)",
                  "extra_terminal_MIX": "decrease lower summation order by one",
                  "M42_bound": result["members"][-1]["complete_marked_clock_TV_theorem_bound"],
                  "proof": "PROOF.md"},
        "boundary": {"all_finite_wrap_transitions_retained": True,
                     "all_integer_B_sectors_retained": True,
                     "spatial_chart_transfer": "same residues and augmented face degrees; energy chart distinct",
                     "M6_complete_clock_bound": 1,
                     "minimal_bounded_fix": "increase only the record modulus, keeping the full lifted preparation and clock",
                     "integer_real_minimizer_identification": False},
        "rejections": result["obstructions"]["admission_mutations_rejected"],
        "claims_excluded": ["3D Maxwell", "spatial continuum", "harmonic oscillator identification",
                            "particle or force dictionary", "native coherent pulse compiler", "unique architecture"],
        "result": "RESULT.json",
        "official_rows_changed": [],
    }


def manifest():
    paths = [
        HERE / "dynamical_field.py", HERE / "verify.py", HERE / "PROOF.md", HERE / "README.md",
        HERE / "DYNAMICAL_FIELD_CONTRACT.json", HERE / "RESULT.json", HERE / "dynamical_field_states.npz",
        ROOT / "openwave-tuning/next/C10/joint_source.py",
        ROOT / "openwave-tuning/next/C10/charged_reservoir.py",
        ROOT / "openwave-tuning/next/C10/modular_field_adapter.py",
        ROOT / "openwave-tuning/next/C9/fermionic_shared_records.py",
        ROOT / "openwave-tuning/next/C9/fermionic_uniform_reference.py",
        ROOT / "openwave-tuning/next/C11/boundary_charge.py",
        ROOT / "openwave-tuning/next/C10/spectral_readout.py",
        ROOT / "openwave-tuning/next/ROW1_NATIVE/routes/quantum_spatial/quantum_spatial.py",
        ROOT / "openwave-tuning/candidate/oph_openwave/fiber_preparation.py",
        ROOT / "openwave-tuning/candidate/oph_openwave/uniform_native.py",
        ROOT / "openwave-tuning/candidate/oph_openwave/modular_native.py",
        ROOT / "openwave-tuning/candidate/oph_openwave/construction.py",
        ROOT / "openwave-tuning/validation/bridges/providers/fock.py",
        ROOT / "openwave-tuning/next/CHAIN/clock/CLOCK_CONTRACT.json",
        ROOT / "openwave-tuning/next/CHAIN/clock/common_clock_states.npz",
        ROOT / "openwave-tuning/next/CHAIN/coordination/INTEGER_FIELD_CONTROL.json",
        ROOT / "reverse-engineering-reality/model/native_wilson_maxwell_refinement/PROOF.md",
        ROOT / "reverse-engineering-reality/model/field_record_addenda/native_wilson_maxwell_spatial_family/DYNAMICS.md",
    ]
    # Bind every project Python module actually loaded by this regeneration.
    for module in tuple(sys.modules.values()):
        filename = getattr(module, "__file__", None)
        if filename:
            path = Path(filename).resolve()
            if path.is_relative_to(ROOT) and path.suffix == ".py":
                paths.append(path)
    return {"schema": "oph.c10.dynamical_field_manifest.v1",
            "sha256": {str(p.relative_to(ROOT)): digest(p) for p in sorted(set(paths))}}


def main(write=False):
    result, arrays = execute()
    assert all(result["checks"].values()), result["checks"]
    canonical_contract = contract(result)
    if write:
        (HERE / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(HERE / "dynamical_field_states.npz", **arrays)
        (HERE / "DYNAMICAL_FIELD_CONTRACT.json").write_text(json.dumps(canonical_contract, indent=2) + "\n")
        (HERE / "MANIFEST.json").write_text(json.dumps(manifest(), indent=2) + "\n")
    assert json.loads((HERE / "RESULT.json").read_text()) == result, "result regeneration mismatch"
    assert json.loads((HERE / "DYNAMICAL_FIELD_CONTRACT.json").read_text()) == canonical_contract, "source/preparation/clock contract mismatch"
    with np.load(HERE / "dynamical_field_states.npz") as stored:
        assert set(stored.files) == set(arrays)
        for key, value in arrays.items():
            assert np.allclose(stored[key], value, atol=1e-12, rtol=0), key
    assert json.loads((HERE / "MANIFEST.json").read_text()) == manifest(), "source or artifact bytes changed"
    review = {"schema": "oph.c10.dynamical_field_regeneration.v1",
              "independent_review": False,
              "checks": {"all_source_and_dynamics_checks": True, "result_regenerated": True,
                         "complete_archived_states_regenerated": True, "exact_contract_regenerated": True,
                         "manifest_verified": True, "source_preparation_clock_mutations_rejected": True},
              "manifest_sha256": digest(HERE / "MANIFEST.json"),
              "scope": "Producer regeneration and rejecting mutations; no independent-auditor claim."}
    (HERE / "REVIEW.json").write_text(json.dumps(review, indent=2) + "\n")
    print(json.dumps(review, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="regenerate local package artifacts before verification")
    main(parser.parse_args().write)
