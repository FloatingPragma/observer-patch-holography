# BD threshold and spectrum gate

This directory carries the reproducible threshold gate for GitHub issue 368
and the transverse-rank/isolation certificate for issue 369. The
Bouchard-Donagi papers determine a visible massless-cohomology branch. They do
not determine the stabilized moduli point, nonzero-mode spectrum, string scale,
hidden completion, or a map from its supersymmetric compactification to the
non-supersymmetric OPH low-energy branch.

The builder therefore has three outputs with different logical roles. It
reproduces the Higgs/top, stop, and one-loop gauge-running proxy coordinates
from frozen inputs. It separately emits a fail-closed threshold certificate
with compatibility_evaluated set to false and a fail-closed rank-obstruction
certificate with promotion_allowed set to false.
Supplying benchmark MSSM values does not close the gate. The OPH low-energy
target used here has no supersymmetric partner sector. The interface accepts
either a conventional SUSY-breaking route with mediation and soft boundary
conditions, or a genuinely non-supersymmetric UV deformation or continuation.
To certify this BD candidate, the latter must be proved OPH-equivalent to the
BD branch and preserve its cited visible data. It is not established by a
projection label or hash: it must also prove its worldsheet or modular
consistency, anomaly cancellation, vacuum stability, spectrum, couplings, and
thresholds. An unrelated construction is a different candidate.

The source packet is:

    code/particles/data/oph_bd_threshold_spectrum_inputs.json

The generated bundle is:

    code/particles/runs/uv/oph_bd_threshold_spectrum/

Regenerate it from the repository root with:

    python3 code/particles/uv/build_oph_bd_threshold_spectrum_receipts.py

Verify the committed bundle bytes and the recorded local dependency hashes,
selectors, and precision without writing files:

    python3 code/particles/uv/build_oph_bd_threshold_spectrum_receipts.py --check

Run the focused tests with:

    python3 -m pytest -q code/particles/uv/test_oph_bd_threshold_spectrum_receipts.py

The builder requires every non-null UV receipt to have
`status=hash_and_envelope_verified`, validates its repository-relative path and
SHA-256 value, and binds its JSON envelope to issue 368, the exact BD branch,
and the named input slot. This verifies identity and bytes, not the receipt's
scientific content. It does not download the external literature sources;
their independently checked hashes are recorded for a separate fetch audit.

This bundle is deliberately an open-evidence generator, not a compatibility
evaluator: both promotion flags remain false even if all receipt slots are
filled. An evaluator may consume a complete, hash/envelope-bound packet only after
it also performs a real forward spectrum run, threshold matching in one
declared scheme, a vacuum-stability test, and the constraint-augmented
rank/isolation test. Target-side proxy agreement
alone can never alter either promotion flag.

The issue-369 certificate records the two-complex-dimensional ambient
determinantal normal space for the codimension-two one-Higgs pullback locus.
The source packet does not contain the local matrix entries and Jacobian minor
needed to certify pullback transversality, so the corresponding rank-two
pullback statement is conditional.
After known presentation redundancies are removed, the documented
pre-completion one-Higgs source has

    11 + 11 + (51 - 2) = 71 complex = 142 real

directions. The comparison registry has five real coordinates, only three with
promoted OPH-surface status, and no common threshold scheme. Its proxy uses no
BD branch values, so the pullback is constant and has rank zero. Any map on the
full published slice with those five outputs has rank at most five and
infinitesimal nullity at least 137.

Those bounds do not determine the dimension after stabilization, vacuum, and
decoupling constraints. The packet supplies neither that completed constraint
locus nor a selected point or physical Jacobian. The moduli-locking gate
therefore fails. This retracts the operator-safe selected BD candidate status
while retaining the structural BD audit row and leaving the recovered OPH core
unchanged.

The companion simulator contract in
`oph-physics-sim/docs/STRING_VACUUM_SELECTION_RECEIPT_CONTRACT.md` defines the
required quotient, stability, spectrum, threshold, augmented-Jacobian,
interval-isolation, branch-coverage, and catalogue receipts. Its verifier
recomputes interval contraction algebra and refuses producer promotion flags;
string-specific semantic verifiers remain open.
