"""Adversarial source parity and alternative diagonalization for two ports.

The portable implementations are compared to the independently constructed
native observer matrices and literal repair operation in a pinned public
checkout. A separate Hermitian eigensolver checks the circulant spectrum.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


def audit(rer: Path) -> dict:
    pin = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"
    if subprocess.check_output(["git", "-C", str(rer), "rev-parse", "HEAD"], text=True).strip() != pin:
        raise ValueError("wrong RER commit")
    subprocess.run(["git", "-C", str(rer), "diff", "--exit-code", pin, "--", "model", "code"], check=True)
    sys.path.insert(0, str(rer))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from model.lepton_register_round import register as native
    from model.charge_conjugation_dynamics import symmetry
    from oph_openwave import leptons, conjugation
    actual = native.full_tick()
    portable = leptons.register_tick()
    exact_tick = np.array_equal(actual, portable)
    altered = portable.copy()
    altered[0, 0] += 1
    broken_tick = not np.array_equal(actual, altered)
    missing_orientation = not np.array_equal(actual, leptons.register_tick(oriented=False))
    mismatch = 0
    for state in conjugation.STATES:
        for coin in (0, 1):
            target, record = conjugation.step(state, coin)
            expected = symmetry.local_step(state, coin)
            mismatch += (*target, record[1], record[2]) != expected
    # A Hermitian matrix eigensolver does not use the three-character formula.
    shift = np.roll(np.eye(3, dtype=complex), 1, axis=0)
    theta = 2 / 9
    matrix = np.eye(3) + (np.exp(1j * theta) * shift + np.exp(-1j * theta) * shift.conj().T) / np.sqrt(2)
    ev = np.linalg.eigvalsh(matrix)
    ratios = (ev[1:] / ev[0]) ** 2
    target = leptons.comparison()
    references = np.array([float(target["muon_electron_ratio"]), float(target["tau_electron_ratio"])])
    error = float(np.max(np.abs(ratios / references - 1)))
    altered_theta = theta + 1e-4
    altered_matrix = np.eye(3) + (np.exp(1j * altered_theta) * shift + np.exp(-1j * altered_theta) * shift.conj().T) / np.sqrt(2)
    altered_ev = np.linalg.eigvalsh(altered_matrix)
    altered_ratios = (altered_ev[1:] / altered_ev[0]) ** 2
    broken_phase = bool(np.max(np.abs(altered_ratios / references - 1)) > 1e-5)
    # Koide alone cannot validate the phase: the changed phase retains the identity.
    changed_phase_koide = float(np.sum(altered_ev ** 2) / np.sum(altered_ev) ** 2)
    checks = {"all_432_squared_tick_entries_match_native_observer": bool(exact_tick),
              "all_2592_branches_match_literal_native_law": mismatch == 0,
              "independent_Hermitian_spectrum_matches": error < 1e-9,
              "changed_phase_preserves_Koide_while_changing_ratios": abs(changed_phase_koide - 2 / 3) < 1e-14}
    mutations = {"changed_tick_entry_detected": bool(broken_tick),
                 "missing_oriented_record_detected": bool(missing_orientation),
                 "changed_phase_detected_by_spectrum": broken_phase}
    if not all(value is True for value in (*checks.values(), *mutations.values())):
        raise ArithmeticError("source-parity audit or broken-control sensitivity failed")
    return {"schema": "oph.openwave.source_parity_audit.v1", "source_commit": pin,
            "checks": checks, "mutations": mutations, "conjugation_mismatches": mismatch,
            "alternative_eigenvalue_ratio_relative_error": error,
            "changed_phase_Koide_ratio": changed_phase_koide,
            "verdicts": {"native_tick_and_class_conjugation": "CONFIRMED",
                         "conditional_lepton_coefficients": "CONFIRMED",
                         "Koide_identity_fixes_phase": "REFUTED",
                         "physical_mass_and_annihilation_claims": "NOT_ESTABLISHED"},
            "scope": "Two declared finite sources; source-code parity and alternative diagonalization, with no empirical prediction, annihilation or physical mass identification.",
            "auditor_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.rer.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("PASS: independent native source parity and spectrum; three mutations detected")
