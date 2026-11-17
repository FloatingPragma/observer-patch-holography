# Preparation cap and exact-kernel scope

Added after execution and before seal in response to the same-group reader.
The preregistered proof and its hash are retained unchanged.

The exact symmetric outer Markov transfer and the finite-well matrix in
PROOF.md use exact preparation: repeat each relocation edge until it moves.
A probe may then be capped at256 attempts; its acceptance is exactly the
escape-by-256 probability. This is a complete almost-surely terminating
native protocol with all relocation costs charged.

The actual calibration also caps each relocation edge at128 attempts. Its
preparation failure probability conditional on any prior record is at most
25·2^-128. Because path length depends on the previous inner position, the
fully capped outer marginal need not itself be exactly Markov. Each such
conditional outer transition differs in total variation from the exact-
preparation, capped-probe transition by at most25·2^-128; over k attempts,
a coupling bounds the full outer-history discrepancy by k·25·2^-128.
Failures produce a hold and remain in the archive. The observed absence of
failures does not eliminate this theoretical error.

The finite-well eigenvalues in RECEIPT.json are exact-population consequences
(up to displayed numerical arithmetic) of the exact-preparation/capped-probe
controller. They are not a claim of an exactly Markov fully capped process.
The registered martingale comparison already includes this preparation
allowance. The native controller interaction and strict well survive this
allowance; no hydrogenic ladder or production field trajectory was tested.
