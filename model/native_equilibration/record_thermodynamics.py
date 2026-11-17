"""Coherent C/E thermodynamics with the complete native repair tape.

Fix the cut in cut_completion, and choose each complementary seam and fair
coin independently of the input and previous records. A labelled branch is
K_a=V_a/sqrt(m), m=2*number_of_seams. Direct source enumeration proves V_a
is a partial permutation, V_a* is another branch with the same weight, and
[V_a,sigma_bar]=0 for sigma_bar=completed_palm(). Hence the unread channel
is unital, sigma_bar-preserving and symmetric for both inner products
Tr(sigma_bar X*Y) and Tr(sigma_bar^(1/2) X* sigma_bar^(1/2) Y). These are
quantum GNS and KMS detailed balance; the parameter counts native attempts.

For any state rho_SR and untouched reference R, set tau=sigma_bar tensor
rho_R on supp(rho_R). Every finite unread composition preserves tau and
its logarithmic energy, giving the full coherent, zero-modular-heat law

  D(rho_SR||tau)-D(rho'_SR||tau) = S(rho'_SR)-S(rho_SR) >= 0.

The inequality is quantum data processing, not an extrapolation of the
population-only Lean theorem. Each branch also preserves the cap blocks,
so the bulk and central modular observables are separately conserved.

Keep J=(all seam, coin, C-current or E-disagreement records), and let
Omega_rho be the resulting classical-quantum state. For any fixed choice
word c, its record words partition the INITIAL source basis. Each word is
a partial isometry on its domain, so the spectrum of its output block is
that of the corresponding principal input block. Let D_c be this input
pinching and Q=average_c [S(D_c(rho_SR))-S(rho_SR)]. Then, for n attempts,

  0 <= Q <= log(66),
  S(Omega_rho) = n log(m)+S(rho_SR)+Q,
  D(Omega_rho||Omega_tau) = D(rho_SR||tau)-Q,
  0 <= Delta S_SR-Q <= n log(m),
  S(J|SR) = n log(m)+Q-Delta S_SR.

Pinching fixes tau. Its entropy increase is relative entropy to the
pinched state, bounded by the log of at most 66 projectors, even with an
arbitrary untouched R. The two losses in the relative-entropy chain are
Q (unread quantum record dephasing) and Delta S_SR-Q (tracing the classical
tape). For a diagonal initial source, Q=0 exactly: endpoint plus the full
tape recovers the input, and no initial classical distinguishability is
lost in the joint state. H(X_n,J)=H(X_0)+n log(m) then holds exactly.

Physical erasure of a stored tape is unnecessary for the unread marginal
or its mixing bound. Its entropy is not omitted from the joint ledger.
The marginal prediction requires subsequent operations to ignore J; it
cannot be used for a selected branch or record-dependent feedback. Actual
memory reset, physical heat, area and laboratory time remain separate.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from . import cut_completion as cut
from ..observer import native, quantum_mixing

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "record_results.json"
SOURCE_PATHS = (
    "model/native_equilibration/record_thermodynamics.py",
    "model/native_equilibration/cut_completion.py",
    "model/observer/native.py", "model/observer/quantum_mixing.py",
    "Lean/Thermodynamics/CapFirstLaw.lean",
    "Lean/Thermodynamics/FluctuationTheorems.lean",
)


@dataclass(frozen=True)
class Branch:
    label: tuple[int, int, int, int]
    pairs: tuple[tuple[int, int], ...]  # (input, output); coefficient 1/sqrt(m)


def choices(cap=cut.CAP):
    """Independent seam/coin alphabet, not selected successful movements."""
    cut.cut_sectors(cap)  # public input validation
    return tuple((seam, coin) for seam in quantum_mixing.seams()
                 if not set(seam).intersection(cap) for coin in (0, 1))


def _law(law):
    if law not in ("C", "E"):
        raise ValueError("native record law C or E required")


def branches(cap=cut.CAP, law="C"):
    """All nonempty actual native instrument branches, including no movement."""
    _law(law)
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    result = []
    alphabet = choices(cap)
    if not alphabet:
        raise ValueError("the complement must contain a native seam")
    for seam, coin in alphabet:
        groups = defaultdict(list)
        for i, x in enumerate(basis):
            y, (_, current) = native.native(x, seam, coin)
            record = current if law == "C" else x[seam[0]]-x[seam[1]]
            groups[record].append((i, positions[y]))
        result.extend(Branch((*seam, coin, record), tuple(pairs))
                      for record, pairs in sorted(groups.items()))
    return tuple(result)


def history_domains(choice_word, law="C"):
    """Complete initial-basis partition for a fixed seam/coin word.

    This includes every record outcome; its pair lists implement the
    corresponding unnormalised partial isometries. No branch is conditioned
    away. The empty word has the single identity domain.
    """
    _law(law)
    legal = set((seam, coin) for seam in quantum_mixing.seams() for coin in (0, 1))
    if type(choice_word) is not tuple or any(c not in legal for c in choice_word):
        raise ValueError("a tuple of native (seam, coin) choices is required")
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    groups = defaultdict(list)
    for i, x in enumerate(basis):
        records = []
        for seam, coin in choice_word:
            y, (_, current) = native.native(x, seam, coin)
            records.append(current if law == "C" else x[seam[0]]-x[seam[1]])
            x = y
        groups[tuple(records)].append((i, positions[x]))
    return tuple((record, tuple(pairs)) for record, pairs in sorted(groups.items()))


def exact_balance_audit(cap=cut.CAP, law="C"):
    """Exact branch reversal and both weighted matrix-unit symmetries."""
    source = branches(cap, law)
    sigma = cut.completed_palm(cap)
    m = len(choices(cap))
    basis = native.states(binary=True)
    input_count = [0]*66
    output_count = [0]*66
    reverse = {}
    rows, cols = [], []
    for b in source:
        inputs, outputs = zip(*b.pairs, strict=True)
        if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
            raise ArithmeticError("a labelled native branch is not a partial injection")
        adjoint = tuple(sorted((j, i) for i, j in b.pairs))
        if adjoint == b.pairs:
            inverse = b
        else:
            candidates = [a for a in source if a.label[:2] == b.label[:2]
                          and a.label[2] == 1-b.label[2] and a.pairs == adjoint]
            if len(candidates) != 1:
                raise ArithmeticError("missing unique equal-weight native adjoint branch")
            inverse = candidates[0]
        reverse[b.label] = inverse.label
        for i, j in b.pairs:
            input_count[i] += 1
            output_count[j] += 1
            if (sigma[i] != sigma[j]
                    or any(basis[i][port] != basis[j][port] for port in cap)):
                raise ArithmeticError("branch changes the completed Palm weight or retained cap")
        for (i, a), (j, b_out) in product(b.pairs, repeat=2):
            rows.append(66*a+b_out)
            cols.append(66*i+j)
    if any(reverse[reverse[a]] != a for a in reverse):
        raise ArithmeticError("record reversal is not an involution")
    if set(input_count) != {m} or set(output_count) != {m}:
        raise ArithmeticError("complete tape instrument is not trace preserving and unital")
    p = coo_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                   shape=(4356, 4356)).tocsr()
    original, denominator = cut.native_cut_counts(cap, law)
    if denominator != m or not cut._zero(p-original) or not cut._zero(p-p.T):
        raise ArithmeticError("labelled instrument does not match the complete native channel")
    for row, col in zip(*p.nonzero(), strict=True):
        a, b = divmod(int(row), 66)
        i, j = divmod(int(col), 66)
        # Equal left/right weights imply both sigma_j (GNS) and
        # sqrt(sigma_i*sigma_j) (KMS) equality without rounded square roots.
        if sigma[a] != sigma[i] or sigma[b] != sigma[j]:
            raise ArithmeticError("weighted GNS/KMS symmetry fails")
    digest = hashlib.sha256(json.dumps(
        [[b.label, b.pairs, reverse[b.label]] for b in source],
        separators=(",", ":")).encode()).hexdigest()
    return {"law": law, "choice_alphabet": m, "record_branches": len(source),
            "all_matrix_units": 4356, "partial_injections_and_adjoint_involution": True,
            "completed_palm_branch_commutation": True, "GNS_and_KMS_symmetry": True,
            "complete_instrument_integer_sha256": digest}


def exact_history_audit(cap=cut.CAP, law="C", attempts=2):
    """Exhaustive finite control of the all-length partial-isometry proof."""
    _law(law)
    if type(attempts) is not int or not 0 <= attempts <= 2:
        raise ValueError("retained exhaustive audit supports zero to two attempts")
    alphabet = choices(cap)
    sigma = cut.completed_palm(cap)
    branch_lookup = {b.label: b for b in branches(cap, law)}
    words = 0
    max_partition = 0
    digest = hashlib.sha256()
    for word in product(alphabet, repeat=attempts):
        domains = history_domains(word, law)
        cover = sorted(i for _, pairs in domains for i, _ in pairs)
        if cover != list(range(66)):
            raise ArithmeticError("history omitted or duplicated an input state")
        max_partition = max(max_partition, len(domains))
        for records, pairs in domains:
            images = [j for _, j in pairs]
            if len(set(images)) != len(images):
                raise ArithmeticError("endpoint plus tape failed to recover the input")
            composed = {i: i for i in range(66)}
            for (seam, coin), record in zip(word, records, strict=True):
                step = dict(branch_lookup[(*seam, coin, record)].pairs)
                composed = {i: step[x] for i, x in composed.items() if x in step}
            if tuple(composed.items()) != pairs:
                raise ArithmeticError("history partition disagrees with composed Kraus branches")
            if any(sigma[i] != sigma[j] for i, j in pairs):
                raise ArithmeticError("equilibrium word and reversed word weights differ")
            digest.update(json.dumps([word, records, pairs], separators=(",", ":")).encode())
            words += 1
    return {"attempts": attempts, "all_choice_words": len(alphabet)**attempts,
            "all_labelled_words": words, "largest_input_partition": max_partition,
            "all_input_choice_pairs": 66*len(alphabet)**attempts,
            "initial_population_recoverable_from_endpoint_and_complete_tape": True,
            "equilibrium_forward_and_adjoint_word_weights_equal": True,
            "integer_history_sha256": digest.hexdigest()}


def run():
    return {"schema": "oph.native_equilibration.record_thermodynamics.v1",
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            "cap": list(cut.CAP), "reference": "completed_palm: E_R(sigma_Palm)",
            "laws": {law: {"detailed_balance": exact_balance_audit(law=law),
                           "complete_history": exact_history_audit(law=law)} for law in ("C", "E")},
            "all_finite_history_identities": {
                "coherent_Clausius": "Delta S_SR = D(rho_SR||tau)-D(rho_n_SR||tau) >= 0; tau=sigma_bar tensor rho_R",
                "pinching_cost": "Q=average_choice_words[S(D_c(rho_SR))-S(rho_SR)]; 0<=Q<=log(66)",
                "joint_tape_entropy": "S(Omega_rho)=n*log(30)+S(rho_SR)+Q",
                "joint_relative_entropy": "D(Omega_rho||Omega_tau)=D(rho_SR||tau)-Q",
                "tracing_tape_loss": "0 <= Delta S_SR-Q <= n*log(30)",
                "conditional_tape_entropy": "S(J|SR)=n*log(30)+Q-Delta S_SR",
                "diagonal_population": "Q=0; H(X_n,J)=H(X_0)+n*log(30); joint relative entropy is conserved"},
            "scope": {"arbitrary_untouched_reference": True,
                      "all_outcomes_and_schedule_randomness_retained": True,
                      "unread_marginal_requires_physical_tape_erasure": False,
                      "subsequent_record_feedback_covered_by_unread_bound": False,
                      "full_quantum_entropy_law_machine_checked_in_Lean": False,
                      "stationarity_under_physical_area_or_metric_variation": False,
                      "physical_heat_temperature_or_clock_calibration": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, sort_keys=True)+"\n"
    if args.check:
        if args.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit("native record thermodynamics receipt drift")
    else:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    print("native coherent and complete-record thermodynamics: all checks passed")


if __name__ == "__main__":
    main()
