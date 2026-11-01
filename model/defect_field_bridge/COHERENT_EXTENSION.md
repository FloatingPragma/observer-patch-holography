# Coherent collision-hold extension

This explicitly designed full-CAR repair law improves the earlier
[`full_domain`](FULL_DOMAIN.md) collision policy while retaining that
flagged policy as a control. Both laws are identical on every original
binary endpoint configuration and on the colored hard-core two-particle
sector.

For each fair coin, append the nonbinary endpoint identity to the
existing raw-zero branch:

\[
K_{b,0}=K^{\rm binary}_{b,0}\oplus I_{\rm nonbinary}/\sqrt2.
\]

The other branches vanish on the nonbinary subspace. Binary outputs are
binary, so they cannot collide with the appended nonbinary outputs.
Every enlarged branch is a partial injection, and every input still
appears exactly once for each fair coin. Thus the full instrument is
trace preserving on all particle-number sectors. Its fair-coin channel
is also unital. Every raw slot remains in the representation; the former
collision flag is now a zero slot.

Assigning the existing raw-zero label to nonbinary holds is a **new law
specification**, including for E. It is not inferred from the original
binary record meaning. No additional observation of the collision
predicate is performed. Recording and then forgetting that predicate
would reproduce the dephasing control instead.

On the subspace where each repaired endpoint has even occupation, a
state either has two empty endpoints or at least one nonbinary endpoint.
Both coins therefore apply `I/sqrt(2)` with raw zero. All coherences in
this subspace, including their correlations with arbitrary references,
are preserved exactly. A pair drive on initially empty endpoints stays
in this subspace. Its inverse-drive echo is one for this policy, compared
with approximately `0.6317667401` for the flagged control at
`phi=pi/2`, `tau=0.002 GeV^-1`, `j=0`.

This is a policy-dependent improvement. It does not remove all repair
backaction: a superposition involving a binary state that moves can
still acquire different records and lose coherence. Both policies and
their complete outputs remain available for comparison with observations.

The package executes all sixteen input columns of the two-site,
two-flavor source dilation and all eight columns of the three-site,
one-flavor crossing witness, for C and E. It retains all 512 source
output addresses and every coin/raw slot. The general 1,080-mode
extension is an exact finite algebraic construction; a full
`2^1080` source isometry is not materialized. Floating source bounds remain
distinct from physical hardware or interval certificates.

```sh
python3 -m model.defect_field_bridge.coherent_extension --check
python3 -m pytest -q -W error model/defect_field_bridge/test_coherent_extension.py
```
