# Wave replay amendment v1

This local software repair follows E-134 and findings WR-01 through WR-04 in
[the sealed review](../../notes/wave_review_20261115/REVIEW.md). It builds on
E-100, E-105, E-109 and E-131. Re-execution is needed to test the repairs;
it is not a new discovery experiment or a resumption of the other row studies.
The original bundle, review, failed replays and registered scientific checks
remain unchanged. No official status or G4/G6 gate is awarded here.

## Frozen implementation and comparison policy

`POLICY.json` names the full source commit and original outer manifest digest.
The adapter extracts raw Git blobs, verifies their object IDs, all 420 outer
manifest entries, 409 copied pins and all six ancestor manifests. An explicitly
provided byte-preserving bundle can be used without Git. CRLF-converted bundles
fail custody; they are never silently normalized or resealed. Scoped root
`.gitattributes` entries disable future checkout conversion of these packages.

Only scratch copies execute. The four scientific leaves retain their original
checks. The free/orbital packaging wrappers are replaced by the explicit
custody, call-count and artifact comparisons in this adapter. Both remaining
leaves run directly. The adapter fails explicitly under `-O`, `-OO` or effective
`PYTHONOPTIMIZE`. Children receive a cleaned environment and check that their
optimization flag is zero before loading any producer code. Every failure has
a nonzero process exit; a receipt records failures without converting them to
success. Fresh output directories preserve prior attempts.

All fields, sequence lengths, JSON types, integers, strings, booleans and nulls
must match exactly. Duplicate keys, nonfinite floats and missing/extra fields
fail. Corresponding finite derived floats must satisfy
`abs(actual - expected) <= 1e-12 + 1e-12 * max(abs(actual), abs(expected))`.
This is a packaging equivalence policy, frozen before amended execution, and
does not replace any scientific threshold. At unit scale it allows roughly
4,500 binary64 machine epsilons for accumulation/eigensolver variation. Near
zero its absolute budget is 100 times smaller than the original `1e-10`
source-residual and clock gates; the original tighter `1e-12` projector gates
still execute unchanged. It is not an error theorem for arbitrary platforms:
platforms outside this budget fail for investigation.

Gzip bytes and JSON formatting may differ. Each receipt's embedded raw digest
must match its own decompressed archive before structured comparison. Only
these three explicitly named archive-digest fields receive this treatment;
all other digest strings compare exactly. Original and generated artifact
hashes and byte comparisons are recorded, including their failures.

## Orbital diagnostic amendment

Every internal band uses normalized projector column **port 0** in the existing
source port numbering. The spatial ground eigenvector has positive first
coordinate. The connected finite spatial operator has a simple strictly
positive ground state, so this convention has no tie or zero coordinate.
Only these two expressions in the orbital scratch program change. The native
law, seeds, schedules, currents, clocks, caps, populations and scientific
thresholds retain their original definitions.

The old scalar readouts are *not* accepted within a larger tolerance. The
reference for the amended diagnostics is reconstructed from original sealed
occupation counts and the sign-canonicalized stored spatial profile, using
polynomial internal projectors at the exact eigenvalues
`0, 5-sqrt(5), 6, 5+sqrt(5)`. The executing leaf instead uses its numerical
eigensolver projectors. The report retains every old/new profile and scalar
readout and hashes the patched source. All other outputs compare against the
original sealed artifacts. Full band projector readouts already used by the
free and unified sources are unchanged.

This changes an orbital diagnostic decoder convention. It must be disclosed
in the author's candidate/pin proposal; it is not silently attributed to the
historical candidate's scalar readout. The free `c=0` wave claim and its source
map retain their scope. No new candidate is selected or registered here.

## Validation and remaining decisions

Regression tests must reject corrupted custody (including a wrong origin),
changed discrete records, malformed structure, nonfinite values, bad embedded
hashes, out-of-budget floats and the former unstable reference. They also
exercise optimized entry and verify that assertions remain active in children.
Run the same frozen implementation on Windows and Linux with the original
scientific assertions enabled, then compare the two output sets using this
same policy. Record the inputs' SHA-256 hashes before the first amended run
in `REGISTRATION.json`; this local chronology is not a platform preregistration.

The proposed equivalence policy requires outside adjudication against the
repository's literal byte-replay G4 rule. Full outside-group proof/source
review, the author's candidate and pin decision, prospective platform
registration and execution, and maintainer row decisions remain outstanding.
