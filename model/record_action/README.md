# A record-inferred kinetic parent and the missing clock/action join

CP-04 / [#1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140).
This component starts from the complete endpoint-biased C/E source and the
record-only inference in [PR #1151](https://github.com/FloatingPragma/observer-patch-holography/pull/1151),
working revision `67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c`. It constructs a
positive local kinetic operator with a nonzero density interaction, derives
its approximation from the existing native record counts, and identifies
precisely why it is not the recovered observer's physical Hamiltonian.

## Actual source and theorem attachment

Use the existing 66 binary two-load states, fixed positive seam probabilities,
repair probability `r>0`, positive port weights `theta`, and complete C/E
outputs. The reference `pi(x)=prod(theta_i^x_i)/Z` is fixed independently of
the repair. Its exact detailed balance, connected support and coin effects
are checked against all native basis branches, including zero currents and
holds. C and E agree on populations but are not identified as quantum channels.

For the **row-stochastic** population kernel `P` and `D=diag(pi)`, define

```
H = I - sqrt(D) P invsqrt(D) = sum_e r q_e Q_e.
```

On a seam `(a,b)` and each spectator `k`, the `(ak,bk)` block is

```
Q_e = [[1-p, -sqrt(p(1-p))], [-sqrt(p(1-p)), p]],
p = theta_a/(theta_a+theta_b).
```

Empty or fully occupied endpoints have zero block. There are ten disjoint
rank-one blocks per seam and 300 moving pairs. Each block is a positive
projector annihilating `(sqrt(theta_a),sqrt(theta_b))`. Therefore `0<=H<=rI`
and `sqrt(pi)` is a common zero state. Connected native support forces every
zero state to be proportional to it: zero quadratic form makes every local
projector vanish, fixing amplitude ratios along every support edge. The
ground state is unique. Its nonzero finite spectrum is the population
relaxation spectrum of `I-P`. Eigenvalues printed in the receipt are numerical
diagnostics; positivity and uniqueness follow from the preceding exact proof.

This is an application of the established stochastic-matrix/Rokhsar--Kivelson
construction, not a new general theorem. See Castelnovo, Chamon, Mudry and
Pujol, [*Generalized Rokhsar--Kivelson Hamiltonians and the Stochastic Matrix
Form decomposition*](https://arxiv.org/abs/cond-mat/0502068). The contribution here is the native branch mapping, record-derived error bound and incompatible-clock
witness on this specific observer source. The existing
[`FiniteConditionalRepair`](../../Lean/Thermodynamics/FiniteConditionalRepair.lean)
provides the conditional-resampling detailed-balance foundation; the actual
finite mapping here is Python exact arithmetic, not a new Lean formalization.

On hard-core occupations the same local term is

```
Q_ab = [theta_b n_a + theta_a n_b - (theta_a+theta_b)n_a n_b
        - sqrt(theta_a theta_b)(b_a^* b_b + b_b^* b_a)]/(theta_a+theta_b).
```

Thus the two-particle diagonal differs from the independent one-particle
lift by `-r q_ab n_a n_b`. At the canonical uniform scheduler this coefficient
is exactly `-1/30` on all 30 native seams. This is a nonzero finite hard-core
interaction with a positive total Hamiltonian, not a relativistic QFT or a
chosen particle species. No new hopping service has been compiled or executed.

## A complete modular-energy ledger

Let `K=-log(pi)` and orient positive native current from a to b. Each complete
C/E record determines

```
Delta K = J_ab log(theta_a/theta_b).
```

The exact certificate checks the exponential of this identity against the
actual before/after states for all 7,920 source/coin/seam columns across both
laws. At stationarity, forward and reverse fluxes agree on each seam, so mean
`Delta K` and population entropy production are zero, although current second
moments are positive. This does not calculate archive erasure, apparatus heat
or an energy cost of maintenance. The ledger concerns the dimensionless
modular potential and attempted native events, including holds.

## Recover the parent from the existing records

The estimator consumes only the 30 `(trials,coin-one-count)` pairs retained
by #1151. Its eleven tree edges select relative weights; nineteen held-out
edges and the shared interval test retain their existing admission/stop rule.
There is no new acquisition claimed here. The source receipt and stream hash
bind these reused counts; the private replay weights are used only by a
separate diagnostic comparing the estimate to the known simulation.

On the simultaneous confidence event, each true coin probability lies in
its admitted interval. If that interval and its inferred center lie in
`[1/8,7/8]`, then `||Q'(p)||=1/(2 sqrt(p(1-p)))<2`. Hence the exact rational
quantity `B=(1/30) sum_e 2 max(|p_hat-lo|,|p_hat-hi|)` bounds `||H_hat-H||`
for the canonical no-hold uniform-seam continuation. Weyl's inequality gives
the same absolute bound for each ordered eigenvalue. Duhamel's identity
gives `||exp(-it H_hat)-exp(-it H)|| <= |t|B`; half the diamond distance is
at most `min(1,|t|B)`, including untouched references and subsequent complete
instruments. Add the inherited simultaneous statistical failure probability
for an unconditional channel bound against an oracle using the same
record-determined stop decisions. This does not bound distance to a protocol
that always proceeds. Stops are outputs, not postselection.
Here `t` is a **new declared parent-action parameter**. It is not automatically
the modular parameter, native event count or physical duration.

## The clock/action obstruction and surviving routes

The recovered modular generator is diagonal in the occupation basis; `H`
has 300 strictly nonzero off-diagonal hops for every admitted positive law.
Thus **no affine identification `H=aK+bI` exists** on this full source algebra.
For the actual unequal weights, an exact nonzero hop between unequal
reference probabilities also proves `[H,pi]!=0`. Consequently this `pi`
cannot be a Gibbs state of `H` at any inverse temperature: such a state must
commute with `H`. The uniform-weight control has scalar `K`, and has
nonzero kinetics; it does not rescue the affine identification. Its reference
can of course be the zero-inverse-temperature Gibbs state of any `H`.

This excludes a specific attempt to promote the reversible population
kernel to the same-reference modular Hamiltonian. It does not exclude an
implemented coherent service, a new reference, additional accessible
degrees of freedom, a non-Gibbs operational clock comparison or refinement.
Each surviving route requires its own complete state/action/clock mapping.
The existing positive modular-clock result remains intact. Common physical
normalization, a localized particle, independent kinetic dispersion and
conserved stress/Ward identities are open. No integrated gate is accepted.

```sh
python -m model.record_action.action --write
python -m model.record_action.action --check
python -m pytest -q -W error model/record_action
```
