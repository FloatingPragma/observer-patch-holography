# Full-Haar coupling dependence and a retained physical clock

The carrier, Gauss basis and clock are the exact shared-link two-plaquette
source in [`openwave_coupled_haar_source`](../openwave_coupled_haar_source/PROOF.md).
Write `x=g²`, `M=2-W1-W2` and `H(x)=x E/2+2M/x` at spacing one.
Every reported coupling is a declared architecture parameter. No measured
mass, running-coupling datum or beta function is fitted.

The seven rational couplings are `1,3/2,2,3,4,6,8`. Each has its own fixed
complete dyadic frame and native ground/low-state trials. The exact
full-complement Schur proof, disjoint residual windows and Temple bound
are the same as in
[`openwave_coupled_haar_spectrum`](../openwave_coupled_haar_spectrum/PROOF.md).
The implementation recomputes all errors against the exact source.
Every representation outside the proof coordinates remains in the bound.

## Exact derivative without an unbounded-observable shortcut

The operator family has a common electric domain and is analytic for x>0.
Its ground is simple, and the separately certified first three full-source
levels identify the reported first excitation as simple and isolated.
The Feynman--Hellmann identity therefore applies at these points. The
algebraic operator identity

`H'(x)=H(x)/x-4M/x²`

gives `E_k'=E_k/x-4<M>_k/x²`. This is especially useful because M is
bounded between zero and four. No trace-distance estimate is applied to
the unbounded electric energy. The ground expectation uses the certified
ground projector. The first-excited trial's full residual, divided by its
separation from both the ground and second-excited intervals, bounds its
projector error. Bounded Wilson expectation errors then enclose both
derivatives and their difference. All seven gap-derivative lower bounds
are positive. Exact endpoint subtraction also gives the finite-step slopes
and ratios; these do not assert derivative bounds between sample points.

## A necessary constant-mass tangent

In the declared clock, a corresponding physical gap is `Delta(x(a))/a`.
Holding that gap fixed requires

`dx/dlog(a)=Delta(x)/Delta'(x)`.

Both numerator and denominator have positive certified intervals at every
reported point. Their interval quotient therefore gives a necessary local
level-set tangent. It is not a renormalization-group beta function, a
uniquely selected coupling law or an accepted spatial refinement. Similarly,
the ratio of two dimensionless gaps is the spacing ratio required to make
those two particular energies equal; it does not match all source observables.

The separate spatial Haar-blocking calculation supplies an actual change
of link and face degrees. A fixed-graph coupling scan alone cannot replace
that spatial construction or its full-Hamiltonian, volume and continuum
requirements. No GeV normalization or physical glueball assignment follows.
