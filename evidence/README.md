# Computational evidence archives

This directory contains scientific data and receipt packages cited by the
papers. Each package states its own replay and claim boundary.

- [`source_net_causal_poset/`](source_net_causal_poset/) is the expected OPH
  causal poset: the source-record family receipts at `q <= 55`, the carrier
  realization receipts and event logs, a property table with the status of
  every property (proved in Lean, computed, or declared), a standard-library
  generator that rebuilds the poset from its definitions, and a checker.
- The [archived local-domain package](https://github.com/FloatingPragma/observer-patch-holography/tree/512b47075b0dd00a570b08cc9b576bec7559acb5/evidence/local_domain) contains the complete 11-file local-domain
  family used by the spacetime, particle, and cosmology papers. It contains
  the compressed stage-one arrays, all four staged receipts, the five
  downstream finite-interface receipts, the binding manifest, and a
  standard-library archive checker.
- [`custody/`](custody/) holds byte-exact copies of the frozen-prediction
  custody directories with their OpenTimestamps proofs; the frozen-prediction
  validator verifies them in every clone.

Evidence directories are not interchangeable. A carrier count shared by two
packages does not make one run evidence for the observable tested by another.
