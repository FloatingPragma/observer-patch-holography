# Source flavour identifiability

The actual native current and Higgs tensors admit 27 independent complex
charged family coefficients and six optional complex symmetric neutral
coefficients. Exact finite field redefinitions leave **13** charged flavour
parameters, **22** with a generic Weinberg coefficient, or **20** on its
generic rank-two stratum. These counts concern the declared canonical
finite action; physical chiral-measure and theta-angle questions are separate.

The existing shared-F3 source trial proves aligned quark squared kernels.
An exact full-width native-tensor counterexample disproves alignment as a
consequence of the gauge/Higgs identities. A complete eighteen-outcome coherent
current-loop instrument measures its nonzero rephasing-invariant CP contrast,
retaining every output family pair. This supplies a concrete additional source equation
and observable for one remaining flavour direction.

The [receipt](results.json) contains exact tangent matrices and nullspaces,
native identities, the counterexample and all readout probabilities. The
[proof](PROOF.md) also explains why the chosen massless-neutrino packet has
one independent Majorana phase rather than two. Empirical Yukawas, CKM and
oscillation inputs retain their ancestry; no physical particle assignment
or source-derived measured coupling is asserted.

```sh
VECLIB_MAXIMUM_THREADS=1 python3 -m model.source_flavour_identifiability --write
VECLIB_MAXIMUM_THREADS=1 python3 -m model.source_flavour_identifiability --check
python3 -m unittest model.source_flavour_identifiability.test_flavour -v
```
