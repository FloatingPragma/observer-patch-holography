# Compact transit state in the local-coarse modular instrument

For a connected graph with incidence B, choose a signed cycle vector chi with
B chi=0 and norm one. The bare port/transit state q=(0,chi) is an exact
H_m eigenvector at energy -m. The canonical carrier face (0,11,5) gives the
three seam coefficients +1 on (0,11), -1 on (5,11), -1 on (0,5), divided by
sqrt(3). It is supported on three transit seams in one carrier of any canonical
ring federation. The state is immobile; the offset m is supplied.

Let V map each position to its normalized complete modular Gauss fiber, as in
[uniform_native.md](uniform_native.md). Every lifted half-hop translates one
fiber bijectively to the next, so H_lift V=V H_bare. The same identity holds
for each local-coarse C/E Kraus row. Every local event treats a transit load as
"elsewhere". Its zero-record Kraus row acts as one half of the identity on q;
its other record rows annihilate q. Hence Phi=V q satisfies

```text
H_lift Phi = -m Phi;
K_(seam,law,coin,0) Phi = Phi/2;
K_(other local records) Phi = 0.
```

There are four zero-record branches per drawn seam, each with probability 1/4.
The density Phi Phi* is invariant under coherent evolution, local commit,
any Bernoulli mixture of the two, and the corresponding supplied per-seam
continuous generator. Complete position readouts have probability 1/3 on each
supported transit seam and zero elsewhere. No transit indicator is emitted.
The implementation retains all local outcomes, including the zero-probability
nonzero-current/contrast branches.

This is spatial compactness with enormous record support. Each position fiber
has 6^(M-N+1) configurations. On flat(1,1) its rank is 40; on flat(2,2) it is
157. Phi has 3*6^rank nonzero joint basis amplitudes. Its uniform record
preparation is nonlocal and supplied. The canonical calculation represents
that complete state implicitly and executes its reduced position/instrument
channel; it does not enumerate or physically prepare the enormous state.
Record measurement can leave V's image.

The 36-state modular triangle is enumerated separately. Its uniform-fiber
intertwiner and selected-state eigen-equation have zero numerical residual;
a coherent/commit step at m=0.23,t=0.37,p=0.29 changes the density by trace
distance 1.3e-15. Choosing one definite Gauss-compatible record state per
transit position instead gives port leakage 0.0831732082 after t=0.37.
Complete global position/record monitoring followed by coherent evolution gives
0.2379296315 leakage. Breaking one edge of the canonical spatial cycle gives
0.1032191680 port leakage. These are rejecting preparation/instrument controls.
Position-only measurement is another channel and can preserve V's code while
destroying this cycle's spatial coherence; code invariance and state persistence
are distinct statements.

The result corrects an unrestricted inference from a leaking definite-record
preparation: a different explicit preparation in the same modular local
instrument has exact conditional persistence. It supplies no electron identity,
physical charge/energy dictionary, moving particle, internal clock, field
magnetic moment or angular-momentum theorem. Its density has no oscillating
readout; its phase exp(i m t) alone is not a physical clock.

Implementation: [compact_native.py](../oph_openwave/compact_native.py).
Tests: [test_compact_native.py](../tests/test_compact_native.py).
The exact cycle and fiber arguments are separate from the retained numerical
exponential error checks. `full_criterion_met` is false.
