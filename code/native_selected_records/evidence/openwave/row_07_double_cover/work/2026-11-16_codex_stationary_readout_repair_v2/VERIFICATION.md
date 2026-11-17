# Verification scope and retained defects

The frozen v5 producer first recorded all 256 cases in 172.0413 seconds. Its default `-O` replay then passed in 173.351 seconds with all 277 scientific files' hashes and mtimes unchanged. A targeted type mutant exposed that its generic comparator treated a boolean as an integer. This is a verification defect, not a changed scientific record. V5 remains byte-identical.

`verify_v6.py` is a verification-only descendant, prospectively bound in `VERIFIER_V6_BINDING.json`. It requires recursive type equality, exact integers/booleans and exact/zero-preserving analytic control fields, validates INDEX and leaf custody before numerical comparison, and provides no first-record option. Its full `-O` replay passed 256 cases in 257.5812 seconds. Fifteen actual subprocess mutants all rejected with nonzero exits; source records were unchanged.

`audit_saved.py` adds saved-law semantic checks before numerical replay: nonnegative physical flag/path probabilities, normalized laws, stopped-path total mass, parity-only continuation and exact analytic control availability. Three further actual `-O` mutants (negative physical probability, wrong parity and tiny positive control probability erased) all reject. `SEMANTIC_VERIFIER_BINDING.json` predates those tests. Its normal read-only audit passed all 256 records. Ramsey floating centers retain a documented 2e-12 roundoff enclosure; this is separate from physical probability laws and the independent complete-service certificate.

`replay.py` composes those two checked read-only functions. The composed wrapper itself has not been redundantly executed after both constituent passes. It does not replace baselines. The official reproduction command is `python3 -B -O replay.py`.

`compact_review.py` was copied alone and executed under `-O`; its complete standard-library payload has no filesystem/network/process access. The independent Chebyshev counter calculation agrees with full source visibilities within 3.04e-14 and costs within 6.36e-12 relative, including tiny positive control branches. It excludes the full workspace/flag and outward interval audit; those have separate bound reviews.
