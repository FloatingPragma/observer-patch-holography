# An interacting spin limit of the fully recorded field repair

This package closes the repair and growing-volume premises of the strong
coupling construction in [FA86](../strong_coupling_refinement/README.md).
It constructs a three-dimensional, fixed-spacing quantum spin dynamics
and a ground state as limits of centered source Hamiltonians, with every
repair outcome retained. A second finite-volume theorem gives a family
whose first spectral spacing tends to zero. Neither statement identifies
a relativistic particle or supplies a continuum limit of the lattice.

## Specified source and protected carrier

On a finite nearest-neighbor graph let there be two CAR modes per site,
real signed hopping `t_xy=+/-c/2`, and the centered Hamiltonian

```
K=T+U D-(U/2)N,    H=K+U|V|/2,    tau=c|E|,    U>=4 tau>0.
```

The overall constant has no effect on dynamics. The projector `P` has
exactly one fermion per site. Its two internal states are a spin one-half.
This is the full-domain colored field instrument implemented by
`model/defect_field_bridge/full_domain.py`, with two flavors, site-outer
ordering, fair binary coins, and all sixteen raw slots for each coin.
It is a specified field extension, not the one-bit occupation instrument
of the bare native carrier. A number-preserving onsite spin frame can be
used to attach these modes to the signed Hubbard coordinates. Such a
frame preserves `P` and the assertion below. No large-system native pulse
execution is inferred from this source-map attachment.

For either C or E, either coin, and either internal state at each endpoint,
both sites occupied once gives exactly `(raw=0, unchanged mask, +1, current=0)`.
Consequently each complete repair Kraus map obeys

```
M_(coin,raw) P = delta_(raw,0) P/sqrt(2).                 (1)
```

This proves a two-state internal memory survives every complete repair,
including its classical record. The algebra check covers all 12,288
protected columns on the cube. An independent 768-column full-Fock check
keeps all 32 coin/raw slots and checks number preservation and injectivity
within each branch. No collision branch is deleted outside `P`.

## Uniform finite recorded-history bound

FA86 constructs the canonical low-band isometry `W` and proves

```
Hspin = sum_xy (c^2/U)(S_x.S_y-1/4),
||W-P|| <= 2 tau/U,
||W* H W-Hspin|| <= eta = 3 tau^3/U^2,
||exp(-itH)P-P exp(-itHspin)|| <= 4 tau/U+|t| eta.       (2)
```

Here `P` on the right denotes the isometric spin inclusion. Take any
initial spin state, possibly entangled with an arbitrary untouched
reference. Interleave `m` free intervals with complete recorded repairs,
using any seam/law policy adapted to the preceding records. Put
`T=sum |t_j|`. Replace the free intervals one at a time. The trace norm
distance of two isometrically evolved states is at most twice their
operator norm distance; subsequent complete instruments are contractions
of the direct-sum trace norm. Equation (1) makes the ideal repairs exactly
identity on spin, with fair coins and raw zero. Thus

```
||rho_actual,cq-rho_spin,cq||_1
 <= min(2, 8 m tau/U + 6 T tau^3/U^2).                  (3)
```

This is a full recorded-output, all-input estimate, not a postselected
success-branch estimate. Initializing `W rho W*` instead of `P rho P*`
costs at most `4 tau/U` in (3). Fixed spin interventions preserving `P`
can also be interleaved by the same telescoping argument. Their own
implementation errors, if any, must be added.

`experiment.py` separately compares full 64-dimensional three-site
evolution and all repair branches on a maximally entangled spin/reference
input. It includes fixed C, fixed E, and adaptive C/E histories. The
nonzero collision-record probabilities are retained. These floating-point
checks corroborate (3); they are not its proof.

## A joint volume/coupling limit, including the repair record

For any growing graph sequence with edge count `E>=2`, choose fixed `J>0`
and the explicit source trajectory

```
c_E=J E^4,    U_E=J E^8,    tau_E=J E^5.
```

Then `c_E^2/U_E=J`, `tau_E/U_E=E^-3`, `eta_E=3J/E`, and `U_E>=4tau_E`.
Give every edge an independent Poisson clock of fixed rate `Gamma`, and
apply its complete repair at each ring. For a fixed time horizon `T`,
there are on average `Gamma E T` repairs and one more free interval.
Integrating (3) over this common classical clock law proves

```
||rho_actual,cq-rho_spin,cq||_1
 <= min(2, 8(1+Gamma E T)/E^3 + 6JT/E) -> 0.             (4)
```

The estimate is uniform over initial states and references on the whole
finite graph. It uses a joint trajectory, so no thermodynamic limit at
fixed source `c,U` or unproved local Schrieffer-Wolff theorem is needed.
For a consistent infinite record use independent clocks and fair coins
on the countable fixed lattice, with a fixed C/E choice per edge. Every
finite region has locally finite records. Restrictions of (4) converge to
these records together with the spin dynamics. Arbitrary box-dependent
global adaptive policies satisfy (3) but are not asserted to define one
common infinite policy.

## Local dynamics and ground state on the infinite cubic lattice

The spin bond `h_xy=J(S_x.S_y-1/4)` has norm `J`. On degree-six graphs,
the first edge in a commutator chain has at most `6|X|` choices and each
subsequent edge at most 12. The Duhamel commutator-chain estimate therefore
gives, for disjoint supports at graph distance `d`,

```
||[alpha_t(A_X),B_Y]||
 <= ||A|| ||B|| |X| sum_(n>=d) (24J|t|)^n/n!
 <= ||A|| ||B|| |X| 2^-d (exp(48J|t|)-1).               (5)
```

The recursive estimate and norm-limit construction used here are the
standard ones in [Nachtergaele--Ogata--Sims (2006), equations 2.9--2.12
and Theorem 2.2](https://arxiv.org/abs/math-ph/0603064). The explicit
constants follow from the above edge counts, not from an assumed finite
propagation speed of the Hubbard Hamiltonians with diverging `c_E`.

For a norm-one observable supported on `s` sites in the radius-`r` cube,
with `0<=r<R`, cutting a radius-`R` cube crosses `6(2R+1)^2` bonds. Integrating (5) gives
the valid bound

```
B_R(T)=min(2,6J s (2R+1)^2 |T| 2^(-(R-r))
                         (exp(48J|T|)-1)).              (6)
```

It tends to zero at fixed `T,r,s`. Finite spin evolutions are therefore
Cauchy on local observables, uniformly on compact times, and define a
strongly continuous automorphism group of the quasi-local spin algebra.
It is interacting: for adjacent bonds `[h_xy,h_yz]` contains the nonzero
three-spin term `iJ^2 sum_abc epsilon_abc S_x^a S_y^c S_z^b`; a local
spin's derivative depends on its neighbor. Combining (4) and (6) gives
an explicit source-to-infinite-spin local expectation bound. Finite
products and local intervention protocols follow by telescoping.

One can use an isometric spin inclusion for these comparisons. If a
unital local algebra attachment is desired, identify the four onsite
Fock states as `C^2_charge tensor C^2_spin`, mapping the first charge
block to `(|up>,|down>)` and the second to `(|empty>,|double>)`, and send
`A` to `I_charge tensor A`. This is an even, norm-preserving onsite
*-homomorphism, commuting with `P` and restricting to the stated spin
observable. It is a declared observable dictionary, not a proof of its
native detector synthesis for unbounded volumes.

There is also a source-derived ground state. All charged sectors have
energy at least `U/2-tau>0` for `H`, whereas the `2^|V|` low-band levels
in the half-filled sector are nonpositive by the min-max principle on
`P`. Thus any full source ground lies in this low band. Its pullback
`nu_E` through `W` is a ground state of `Hlow=W*HW`. For every local `A`,

```
nu_E(A* [Hspin,A]) >= -2 eta_E ||A||^2.                 (7)
```

Indeed the corresponding `Hlow` expression is nonnegative and the
commutator with `Hspin-Hlow` costs at most `2eta_E||A||^2`.
Take a weak-* convergent subsequence of these finite-box states, extended
outside the box. The spin commutator in (7) is local once the box contains
the support and its neighbors; `eta_E` tends to zero. The limit is a
ground state for the dynamics constructed above, with a positive-energy
GNS implementation by the ground-state condition. The source's dressed
and undressed states differ by at most `4/E^3`, so the same local state is
obtained from source ground readouts. This proves existence, not a unique
phase, an effective algorithm for the subsequence, or a particle pole.
All charged source excitations instead have diverging lower bound
`JE^8/2-JE^5`; finite-energy spin observables are neutral.

## A rigorous collapsing finite spectral spacing

For this statement extend the specified signed source to periodic tori
`L x W x W`, with `L` even and `W>=3` odd. Every undirected nearest-neighbor
edge occurs once; choose any real hopping signs of magnitude `c/2`.
FA86 applies to this graph without a bipartite-ground-state assumption,
and the effective spin interaction is translation invariant regardless
of the signs. This is an explicit torus extension, not the old open cube.

Let `A=W^2`, and define the first spacing as `E_1-E_0`, **counting ground
multiplicity**. If the spin ground is degenerate this spacing is zero.
Otherwise it is a spin singlet and a translation eigenvector. Apply the
two twists

```
Q_±=exp(±(2pi i/L) sum_x x sum_(y,z) S^z_(x,y,z)).
```

On total `S^z=0`, translation changes either twist by the factor
`exp(2pi i sum_slice S^z)=(-1)^A=-1`. Both twisted ground vectors are
therefore orthogonal to the ground. Averaging their energy expectations
cancels the sine/current term. There are `LA` longitudinal bonds and
`||S^x S^x+S^y S^y||=1/2`, hence

```
gap_spin <= (JLA/2)(1-cos(2pi/L)) <= pi^2 J A/L.
```

Ordered eigenvalues of `Hlow` differ from those of `Hspin` by at most
`eta`; all other full source sectors are higher. Therefore

```
gap_source <= pi^2 J A/L+6 tau^3/U^2.                   (8)
```

Taking `L=2W^4`, `E=3LW^2=6W^6` and the same joint trajectory yields the
fully rational bound `gap_source <= 5J/W^2+J/W^6 -> 0`, using `pi^2<10`.
Small complete spin diagonalizations and exact twist/translation columns
are independent controls. A collapsing first spacing may reflect ground
degeneracy or symmetry breaking. It does not prove a vanishing gap above
the entire ground space, an infinite-GNS massless particle, Lorentz
dispersion, or measured species.

## Replay and claim boundary

Run `python -m model.protected_spin_limit --write` and then `--check`.
The cold receipt binds current sources and repeats exact algebra and
rational bounds. Historical prerequisite receipt custody is reported
separately; stale inherited bindings cannot become `all_pass` merely
because the new checks succeed. The tests additionally run full-output
finite dynamics. The analytical proof imports the explicitly named
finite reduction and commutator-chain results.

This closes finite field-repair protection, a joint growing-volume
interacting spin limit, ground-state existence, and a controlled finite
spacing-collapse family. It leaves spatial continuum scaling, a shared
relativistic massive/massless theory, physical species and unit calibration,
the native unbounded-volume detector realization, and Einstein's geometric
balance/area/metric premises open. No canonical physical claim is promoted.
