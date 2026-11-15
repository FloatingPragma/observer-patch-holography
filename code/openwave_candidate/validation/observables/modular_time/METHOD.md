# Modular correlations on a native observer state

The experiment starts with an output of the actual modular local-coarse
one-load candidate on `flat_ring(1,1)`, using its explicit uniform Gauss-fiber
preparation. Propagation mass is 0.23, record modulus 6, per-seam rate
0.29/(3*0.37), and native duration 0.1 in inverse hop-amplitude units. The
87-position density is computed with the matrix-free native GKSL scheduler.
No record energy, Gibbs state, thermal bath or depolarizing floor is inserted.
The numerical propagator has no certified error enclosure in this experiment.

The observer algebra is a declared corner `P M_87(C) P`, where P resolves
adjacent ports 0 and 1. The inside-P outcome has probability 0.9035519258604874;
the complementary outcome has probability 0.0964480741395125 and is retained.
Conditioning explicitly defines `rho_obs=P rho_native P/Tr(P rho_native)`.
This is an apparatus/postselection choice, not a hidden partial trace or a
claim that the full 87-position state is faithful. The retained `native_state.npz`
contains the full 87-by-87 native density, P and Q=I-P, both unnormalized
terminal Lüders branches, their nonselective sum and the conditioned corner.
These branches are computed on the same native output, without sequential
measurements during propagation. Their probabilities sum to the native trace.

The numerical observer eigenvalues are 0.0131523564468651 and 0.986847643553135.
They are strictly above the declared 1e-12 numerical faithfulness threshold.
The code rejects a singular state without adding a floor or pseudoinverse.
Restricting a pure state to its one-dimensional support changes the algebra
to M_1, with trivial modular flow; it cannot preserve the discarded off-diagonal
observables. The receipt distinguishes this algebra issue from positivity of
the displayed conditional numerical state.

For a faithful finite density rho and real dimensionless parameter s,

```text
sigma_s(A) = rho^(i s) A rho^(-i s);
K = -log(rho);
sigma_s(A) = exp(-i s K) A exp(i s K).
```

In the standard Hilbert-Schmidt representation, the cyclic separating vector
is `Omega=vec(sqrt(rho))` and the modular operator is
`Delta=(rho^-1)^T tensor rho`. The code explicitly checks
`Delta^(is) (I tensor A) Delta^(-is)=I tensor sigma_s(A)` and
`Delta Omega=Omega`. This demonstrates the finite modular construction as well
as a spectral matrix formula.

The selected observer observable is eigenframe sigma_x, explicitly a
state-dependent apparatus choice. If p0,p1 are the observer eigenvalues and
`d=log(p0)-log(p1)`, its correlation is

```text
C(s)=Tr(rho A sigma_s(A))
    =p0 exp(-i s d)+p1 exp(i s d).
```

Its dimensionless modular frequency is `abs(d)=4.31791472430015`.
The observer state is invariant while the noncommuting observable and its
correlations evolve. A commuting observable and a maximally mixed control
have trivial flow. The latter is a control calculation, not the native
preparation used for the displayed correlation.

Two algorithms agree: spectral powers and general matrix logarithm/exponential.
The correlation also agrees with its closed formula. Group/star properties,
the infinitesimal generator, state invariance and the finite KMS relation are
checked. With the definition above, the KMS shift is minus i:

```text
Tr(rho A sigma_(-i)(B)) = Tr(rho B A).
```

Reversing the real modular sign changes the complex correlation; the wrong
plus-i KMS shift fails by approximately 73.0453 in the selected control. This
sign convention also matters for physical time: if a separately supplied
Gibbs interpretation `rho=exp(-beta H)/Z` were valid, conventional physical
Heisenberg time would be `t=-beta*hbar*s`. No such Gibbs assumption or time
conversion is made here.

The modular parameter is an algebraic state/observer parameter. The declared native
evolution parameter is fixed at duration 0.1 for preparation and is not equated to
s. The experiment supplies no theorem about experienced duration, consciousness,
or subjective phenomenology. It gives a concrete observer-algebra modular flow
and its correlations, with the interpretation and physical clock identification
kept explicit.

`MODULAR_TIME.json` retains the conditional state, modular generator/operator,
observable, assumptions, checks, rejecting controls and source hashes. The
source identities are checked before and after evaluation. `correlations.csv`
contains 161 points for s in [-4,4], suitable for plotting real/imaginary correlation
and observable matrix elements. Three tests pass, including the singular-support
and automorphism controls. No reserved acceptance configuration is executed.

Reproduce with the adjacent `run_modular_time.py` and `test_modular_time.py`.

Public reproduction uses the identical bundle-relative paths under
`validation/observables/modular_time`. Native source hashes are relative to
the scientific bundle root. The export manifest authenticates byte-identical
source, state archive, receipts, CSV and plots. Source files and numerical
artifacts are not rewritten during export.
