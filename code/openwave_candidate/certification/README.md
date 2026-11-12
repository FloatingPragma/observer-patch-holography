# Local scalar-symbol mathematical audit

This directory records a definition-only recomputation and a subsequent local
comparison. The label `certification` does not confer external certification.
The result concerns the declared equal-weight twelve-direction scalar cosine
operator; it supplies no physical-sector bridge, physical lattice scale,
OpenWave geometry mapping, or empirical validation. `full_criterion_met` is
false.

Inspect these four files:

1. [DEFINITION_ONLY.md](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/certification/DEFINITION_ONLY.md): the sole mathematical input supplied
   to the two initial implementations, without expected coefficients.
2. [recompute_a.json](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/certification/recompute_a.json): exact direct-port moment expansion over
   Q(sqrt(5)), generated Legendre normalization, and 40 symbol checks at 80 digits.
3. [recompute_b.json](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/certification/recompute_b.json): sign-paired binomial moments, exact
   spherical averages and vertex normalization, and 50 checks at 110 digits.
4. [comparison.json](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/certification/comparison.json): complete polynomial comparison,
   normalization and ratios, 130-digit reconstruction of all 90 numerical rows,
   a deliberately modified-coefficient rejection, and hashes of the exact inputs.

Both initial agents were context-isolated from OPH sources, target values, and
each other's output until their implementations and receipts were saved. Both
used SymPy and mpmath and both used exact spherical monomial averages; their
independent implementations provide a local computational cross-check, not an
organizationally independent review. The comparison phase was permitted to
read both results and the pinned RER reference. B's subsequent path-provenance
cleanup was checked to leave every non-provenance receipt field unchanged.

The exact comparison accounts for A including powers of `a` inside each
homogeneous polynomial and B factoring them outside. Its negative control
changes a nonheadline sextic polynomial coefficient, so matching only the four
scalar coefficients cannot pass. The numerical checks use analytic eighth-order
Taylor bounds; A additionally uses the valid alternating-series signed interval.

The RER reference is pinned to commit
`3bcacf1276ddb3ead35be2b9adcd63b6643e2d35`. Its receipt describes a conditional
physical branch and explicitly leaves physical comparison unarmed. Agreement
with its mathematical coefficients does not establish its physical premises or
meet its dataset-specific preregistration and empirical gates. Harmonic rank
six here does not denote particle spin.

From the public RER repository root, reproduce the post-seal comparison with:

```sh
python3 code/openwave_candidate/candidate/audits/compare_dispersion.py
```

In the development workspace, run
`python3 openwave-tuning/candidate/audits/compare_dispersion.py` instead. The
script recognizes both layouts; `--rer-repo PATH` can select another checkout.
It reads the source receipt from its pinned Git object, rather than depending
on the checkout's current file contents. Regenerate the comparison after any
input receipt or implementation changes so its hashes remain current.
