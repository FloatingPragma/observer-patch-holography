# Preparing centered grounds and enclosing their local records

## Source and finite preparation

The source is the complete two-mode-per-site CAR system, with the original
Weyl hopping and the particle-hole-centered clock

\[
 K_\Lambda=T_\Lambda+U D_\Lambda-\tfrac U2N_\Lambda,
 \qquad Q_\Lambda=K_\Lambda+\tfrac U2|\Lambda|I
 =T_\Lambda+\tfrac U2\sum_x(N_x-1)^2.                 \tag{1}
\]

The scalar in (1) changes no state, excitation energy or Heisenberg
evolution. A spatial bond has norm at most `c`, including both modes.
All Fock sectors are retained. The onsite change of frame to real signed
hopping acts on the entire local four-state Fock space.

The actual eight-site preparation is supplied by
[`centered_native_preparation`](../centered_native_preparation/README.md).
It executes the 262137-pulse signed amplitude tree and 184-pulse inverse
frame on all 65536 field rows, from supplied blanks. Every first-export
failure is an absorbing record. Its source-certified complete-state error
is less than `10^-7`; its exact dyadic target has much smaller certified
error from the unique full centered ground. The subsequent complete
charged four-bin clock is also executed by that owner. The present
calculator authenticates that existing source-bound receipt and archive;
it does not claim a second execution of that tape.

For a larger finite box with rational `c,U`, the full matrix in (1) is a
computable rational symmetric matrix in the signed frame. Its distinct
eigenvalues are the real roots of its square-free characteristic
polynomial `q`. Exact real-root isolation selects the lowest root `a`.
The polynomial

\[
 P_0={q(Q_\Lambda)/(Q_\Lambda-aI)\over q'(a)}          \tag{2}
\]

means evaluation of the polynomial quotient `q(z)/(z-a)`, including at
`z=a`. It is the complete lowest-eigenvalue projector. Degeneracy presents
no difficulty: the normalization is its integer rank `g=Tr(P_0)`.
Every entry is algebraic and has terminating rational enclosures.
The vector `vec(P_0)/sqrt(g)` purifies `P_0/g`. Certified rational
approximations to this vector feed the existing finite amplitude-tree
compiler. The environment is retained; discarding it only specifies the
source marginal. This constructs a finite preparation program for each
box and requested error without a uniform gap or an unknown-amplitude
oracle. Its cost can be exponential in the full Fock dimension. The
growing tapes are not executed here. Blank registers, coherent addressed
control, command precision and the preparation schedule are the inherited
apparatus inputs; autonomous repair does not select these states.

For `q` logical occupation bits and a nonzero complex integer direction,
the existing compiler uses at most `8*2^q` pulses and has tree operator
error `(q+2)2^-b`. With `p` source-port bits its complete-isometry error
adds at most `29*8*2^q*2^-p`. A direction error `epsilon_v` adds directly
to the prepared-state error. These are occupation-register controls;
the spin-code wrappers are not inserted on the full CAR carrier.
The resulting error bounds every subsequent complete record event, with
no division by its success probability. Arbitrary input replacement
requires an explicit retained reset/environment instrument; supplied
blank inputs suffice for the preparation above.

## Computable energy-density bounds

Take an even-sided cube of side `L`, volume `V=L^3`, and let `E_L` be the
absolute lowest eigenvalue of `Q_L`. Disjoint such cubes tile the lattice.
There are `3/L` missing nearest-neighbour bonds per site. Their norm bound
and a product of even finite ground projectors give

\[
 E_L/V-3c/L\ \le e\ \le E_L/V.                       \tag{3}
\]

For the upper bound every cross-cube hopping expectation is zero: an
even product state has zero expectation of each odd one-mode factor.
No boundary ground is substituted for an infinite-ground marginal.
Equation (3) establishes a computable energy density, with error at most
`3c/L` plus the certified finite eigenvalue error. Onsite terms are all
retained. Their norm does not enter the surface term.

Let `A=A*` be an even local observable with algebraic matrix entries and
norm at most a known rational `a`. Its coordinate
spans are `s_i=max(X_i)-min(X_i)`. Add, only as a variational diagnostic,
`h sum_x theta_x(A)`, with rational h, to the Hamiltonian. The preparation and clock in (1)
are unchanged. The finite diagnostic keeps every translate whose support
lies inside the cube. For `L>max s_i`, the density of omitted translates is

\[
 \nu_L=1-\prod_{i=1}^3(1-s_i/L).
\]

If `[ell_L(h),u_L(h)]` encloses the full finite minimum, tiling gives

\[
 {\ell_L(h)\over V}-{3c\over L}-|h|a\nu_L
 \le e(h)\le {u_L(h)\over V}+|h|a\nu_L.             \tag{4}
\]

Parity is sufficient for the product argument, even when `A` changes
particle number. Odd expectations vanish in the even reference and do
not require odd Hamiltonian perturbations. The same surface estimates
prove existence of `e(h)`. For example, compare a larger cube to full
`L`-tiles and bound its remaining strip by the finite interaction norm;
the strip fraction tends to zero at fixed `L`.

Every translation-invariant ground of (1) has energy density `e`.
A direct justification uses the local ground inequality for the Kraus
operators of a reset to a finite-box ground. Resetting leaves disjoint
even terms unchanged and changes only boundary terms besides the box
Hamiltonian. Its possible energy improvement is bounded by twice the
boundary norm. Division by volume and exhaustion give the variational
minimum. This argument requires no thermodynamic spectral gap.

For any such ground `omega`, the variational inequality is
`e(h)<=e(0)+h omega(A)`. Thus, for `h>0`,

\[
 {e(h)-e(0)\over h}\le\omega(A)
 \le {e(0)-e(-h)\over h}.                            \tag{5}
\]

The lower bounds in (4) and the upper bound `u_L(0)/V` in (3) therefore
give explicit rational outer intervals for every selected ground's local
expectation. Finite algebraic eigenvalue isolation makes each interval
computable. Error terms scale as `c/(Lh)+a nu_L`, so narrowing `h` while
holding `L` fixed is invalid. These intervals can be intersected over
successive sizes and perturbations.

For a region containing `s` sites, use the Pauli basis on its
`4^s`-dimensional Hilbert space, with `16^s` basis elements.
If each expectation lies in `[l_j,u_j]` and a specified finite reference
has Pauli expectations `v_j`, then every ground marginal obeys

\[
 \tfrac12\|\rho_X-\sigma_X\|_1
 \le \min\left(1,\tfrac12\sqrt{\sum_{j\ne I}
       \max(|v_j-l_j|,|v_j-u_j|)^2}\right).           \tag{6}
\]

This follows from Pauli orthogonality and `||M||_1<=sqrt(4^s)||M||_2`.
Odd coefficients are exactly zero for even states. Preparation error of
`sigma_X` adds to (6). Complete quantum channels contract trace distance,
so (6), the existing local-clock cutoff error and its source-service
error add before any event is selected. All complementary preparation
and detector records are retained. An event's original probability is
bounded directly; a small event is not conditionally renormalized.

Equations (4)--(6) give computable enclosures, rather than a presupposed
state oracle. The actual default source is not proved to have a unique
translation-invariant marginal. Consequently this argument does not
assert that every requested small trace-error certificate terminates.
It does not identify an effective common subsequence or a simultaneous
volume, time and physical-scaling modulus. Finite energy-density accuracy
alone is insufficient to replace this marginal calculation.

## A uniform actual one-site marginal bound

For the translation-, gauge-, hidden-spin- and particle-hole-invariant
centered reference, write

\[
 \rho_x=\operatorname{diag}(d,\tfrac12-d,\tfrac12-d,d).
                                                               \tag{7}
\]

This is the complete four-state CAR marginal, not its singly occupied
compression. A perfect matching of the cubic lattice admits a product
of the exact two-site centered singlet ground. Cross-dimer hopping
expectations vanish. At hopping amplitude `t=c/2`, its energy per site
for (1) is `-delta`, where

\[
 \delta={\sqrt{U^2+4c^2}-U\over4}>0.                \tag{8}
\]

Hence the selected ground has `e=e_T+Ud<=-delta`. On an edge let
`P=P_single,x P_single,y`. Exact source hopping satisfies `P T_e P=0`
and `||T_e||<=c`. Put `q=omega(1-P)<=4d`. Cauchy--Schwarz on its two
off-diagonal blocks and the complementary diagonal block gives

\[
 |\omega(T_e)|\le2c\sqrt{q(1-q)}+cq
                 \le3c\sqrt q\le6c\sqrt d.
\]

Three edges per site give `e_T>=-18c sqrt(d)`. Independently, the complete
one-particle hopping matrix has zero trace and
`Tr(h_L^2)<=3c^2 V`. Minimizing `dGamma(h_L)` over the full Fock space
and applying Cauchy--Schwarz yields `e_T>=-c sqrt(3/2)`. Therefore

\[
 {\delta^2\over324c^2}\le d
 \le {c\sqrt{3/2}-\delta\over U}.                  \tag{9}
\]

At the declared default `c=1/8,U=1`, the exact rational inequalities
`sqrt(17)>103/25` and `sqrt(3/2)<49/40` imply

\[
 \boxed{\quad 1/90000<d<233/1600.\quad}             \tag{10}
\]

Every eigenvalue of (7) is consequently greater than `1/90000`.
This supplies uniform one-site faithfulness and a bounded one-site
modular generator `||-log rho_x||<log(90000)` for this symmetric
thermodynamic reference. It does not fix `d`, establish faithfulness of
every larger region, select a physical temperature, or prove a continuum
particle branch. The finite source controls verify the whole-Fock bond
identity, the dimer energy and the kinetic trace-square normalization.

## Executed complete local records of the prepared cube

The finite control consumes the authenticated, unnormalized full-row
success vector of the existing native cube preparation. Its earlier
first-failure archive is retained. Three native pulses copy the first
occupation bit to a blank tag. The tag is destructively committed, keeping
all six raw values. The two proper classical columns are kept separately
and a second supplied blank tag receives the second three-pulse copy.
Measurement is not assumed to reset a quantum register. This is an actual
six-pulse execution on 131072 field/tag rows. It does not merge the two
classical histories coherently or divide by their weights.

There are 36 final raw pointer pairs, 12 first-tag export-failure slots,
and 72 second-tag failure slots labelled by the first raw record. The
earlier complete preparation archive is additional to these records.
The same source generator bounds certify preparation plus tagging error
less than `10^-6`. Each additional failure record has original probability
at most `10^-12`, since its ideal map is zero and the complete-isometry
error bounds its amplitude. Unsupported proper-code pointer values are
retained null maps. The admitted original probabilities are

| Local occupation | Original probability interval |
|---|---|
| Vacuum | `[0.015954, 0.015958]` |
| Up only | `[0.484042, 0.484046]` |
| Down only | `[0.484042, 0.484046]` |
| Pair | `[0.015954, 0.015958]` |

The intervals enclose the exact integer target probabilities with the
complete source-control error included. The six-pulse floating execution
is checked against those bounds; its rounding diagnostics do not supply
the mathematical probability certificate. The native preparation tape
is authenticated rather than executed again in this calculation.

The same certified cube energy gives a fully source-derived instance of
(3): its `Q` ground energy lies strictly between `-6203/50000` and
`-2481/20000`, and the infinite energy-density interval adds the explicit
`3c/L=3/16` lower surface allowance. Its width is too large to certify a
small infinite-ground marginal error. The finite records in the table
are therefore statements about the prepared cube. Applying (4)--(6) to
larger boxes is the specified computable test for transferring a local
experiment to the selected infinite reference.

There is a coarse uniform transfer at one site. The certified exact cube
ground has `0.015955<d_cube<0.015957`. Combining this with (10) and the
diagonal form (7) gives

\[
 \tfrac12\|\rho_{x,\mathrm{cube}}-\rho_{x,\infty}\|_1
 =2|d_{\mathrm{cube}}-d_\infty|\le0.25934.
\]

Adding the complete native preparation error gives a strict bound `13/50`
for the one-site state together with its retained success/failure flag,
relative to any of the stated symmetric infinite grounds with the ideal
success flag. This controls any subsequent complete experiment confined
to that site. It does not control several-site correlations, subsequent
nonlocal time evolution, or certify arbitrarily small marginal error.
