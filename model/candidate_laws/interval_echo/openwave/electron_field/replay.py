"""Verify the preserved records and replay the electron comparison experiments."""

import argparse
from hashlib import sha256
import json
from math import isclose
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
MODEL_ROOT = HERE.parents[4]
EXPERIMENTS = (
    ("interval_echo_electron_field", "measure", "results"),
    ("interval_echo_electron_field", "waves", "waves"),
    ("interval_echo_electron_field", "thermal_witness", "thermal_witness"),
    ("interval_echo_electron_magnetism", "measure", "results"),
    ("interval_echo_electron_magnetism", "blank_loop", "blank_results"),
    ("interval_echo_electron_spin", "spin_blind", "spin_blind"),
    ("interval_echo_electron_spin", "three_charges", "three_charges"),
    ("interval_echo_collective_spin", "search", "results"),
    ("interval_echo_collective_spin", "binding", "binding"),
    ("interval_echo_scale_audit", "laboratory", "laboratory"),
    ("interval_echo_electron_rotation", "run", "results"),
    ("interval_echo_electron_rotation", "ring_turn", "ring_turn"),
    ("interval_echo_electron_rotation", "refine", "refined"),
    ("interval_echo_electron_flux_binding", "search", "results"),
    ("interval_echo_electron_flux_binding", "energy_locality", "energy_locality"),
    ("interval_echo_electron_flux_binding", "soft_cloud", "soft_cloud"),
    ("interval_echo_electron_flux_binding", "winding", "winding"),
    ("interval_echo_electron_flux_binding", "bridge", "bridge"),
    ("interval_echo_electron_flux_binding", "free_bridge", "free_bridge"),
    ("interval_echo_electron_species", "measure", "results"),
    ("interval_echo_electron_species", "twist", "twist"),
    ("interval_echo_electron_species", "cloud", "cloud"),
    ("interval_echo_electron_symmetry", "measure", "results"),
    ("interval_echo_electron_symmetry/local_end", "measure", "results"),
)


def compare(actual, expected, path):
    """Match complete JSON records while allowing numerical-library roundoff."""
    if isinstance(expected, dict) and isinstance(actual, dict) and actual.keys() == expected.keys():
        for key in expected:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list) and isinstance(actual, list) and len(actual) == len(expected):
        for i, child in enumerate(expected):
            compare(actual[i], child, f"{path}[{i}]")
    elif isinstance(expected, float) and type(actual) in (float, int):
        # Compare the forbidden exchange probability at its own much smaller scale.
        tolerance = 2e-29 if path.endswith("maximum_two_spin_Fock_probability_over_all_inputs_and_sampled_times") else 2e-12
        if not isclose(actual, expected, abs_tol=tolerance, rel_tol=2e-11):
            raise ValueError(f"{path}: expected {expected}, got {actual}")
    elif type(actual) is not type(expected) or actual != expected:
        raise ValueError(f"{path}: expected {expected}, got {actual}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--family-only", action="store_true", help="replay only the charged-family experiments")
    mode.add_argument("--frozen-flight-only", action="store_true", help="run the matched native no-ring control only")
    mode.add_argument("--verify-only", action="store_true", help="verify the archive without numerical dependencies")
    args = parser.parse_args()
    digests = {}
    for entry in (HERE / "records.sha256").read_text().splitlines():
        digest, name = entry.split("  ", maxsplit=1)
        if name in digests or Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError(f"Invalid manifest entry {name}")
        digests[name] = digest
    with ZipFile(HERE / "records.zip") as archive:
        if sorted(archive.namelist()) != sorted(digests):
            raise ValueError("Archive inventory differs from the manifest")
        for name, digest in digests.items():
            if sha256(archive.read(name)).hexdigest() != digest:
                raise ValueError(f"Digest mismatch for {name}")
        print(f"Verified {len(digests)} archived files.", flush=True)
        if args.verify_only:
            return
        with TemporaryDirectory() as work:
            root = Path(work).resolve()
            (root / "reverse-engineering-reality").symlink_to(MODEL_ROOT, target_is_directory=True)
            archive.extractall(root)
            env = os.environ | dict.fromkeys(("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"), "1")
            if not args.family_only:
                subprocess.run([sys.executable, "-B", "-W", "error", str(HERE / "frozen_flight.py"),
                                str(root / "experiments")], cwd=root, env=env, check=True)
            if args.frozen_flight_only:
                return
            studies = [(root / "experiments" / directory, stem, record)
                       for directory, stem, record in EXPERIMENTS]
            studies.extend((HERE / directory, "search", "results")
                           for directory in ("flux_string", "dressed_exchange", "endpoint_exchange"))
            studies.append((HERE / "endpoint_exchange", "perimeter", "perimeter"))
            if args.family_only:
                studies = []
            family = HERE / "lepton_family"
            studies.extend((family / folder, "measure", "results")
                           for folder in ("bands", "family_projector", "cyclic_charge", "cyclic_hops", "covariant_field",
                                          "clocked_charge", "dirac_read", "identity_audit", "arrival_spin", "transfer", "decay"))
            studies.append((family / "hierarchy", "search", "results"))
            studies.extend((family / "arrival_spin", stem, record) for stem, record in
                           (("circulation", "circulation"), ("tomography", "counts"), ("winding_turn", "winding_turn"), ("incidence_turn", "incidence_turn"), ("rest_partner", "rest_partner"),
                            ("anchor_spectrum", "anchor_spectrum"), ("cadence_audit", "cadence_audit"),
                            ("rest_impulse", "rest_impulse"), ("recoil", "recoil"),
                            ("localized_recoil", "localized_recoil")))
            if not args.family_only:
                studies.append((family / "arrival_spin", "pair_phase", "pair_phase"))
            for local, stem, record in studies:
                output = root / f"{local.name}_{stem}.json"
                subprocess.run([sys.executable, "-B", "-W", "error", str(local / f"{stem}.py"),
                                "--output", str(output)], cwd=root, env=env, check=True)
                compare(json.loads(output.read_text()), json.loads((local / f"{record}.json").read_text()), stem)
                print(f"Reproduced {local.name}/{stem}.", flush=True)
            if args.family_only:
                return
            subprocess.run([sys.executable, "-B", "-W", "error", str(HERE / "dark_contact.py")],
                           cwd=root, env=env, check=True)


if __name__ == "__main__":
    main()
