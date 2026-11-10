# Unscreened charge costs and dynamical screening

## Source, state, reference and clock

The source is the complete eight-cell construction of
[charged transport](../interacting_charged_transport/PROOF.md): 240 native
fermion modes, eight five-state scalars, 24 complete 119-state links,
all local coefficients and the link excitation penalty equal to one,
and `H_lambda=H0+lambda T`, `lambda=2^-34`, `||T||<=360`.
Write `v=lambda*360`. The complete 30-mode static reference is fixed at
cell zero with Hamiltonian `mu N_ref`, `mu=4`. The selected pump branch
occupies exactly one reference mode in either native e^c or L. It does
not move the reference. The same total Hamiltonian and its clock govern
the spectral and localization statements below.

The inherited physical preparation targets the actual neutral interacting
ground g_lambda, with its complete preparation failure retained. Its
P0 localization, native elementary pump and energy-bin selection are
explicit interventions. The localization theorem holds for every density
matrix supported in the selected full-H energy bin. It therefore requires
no purity assertion and no subtraction of the preparation error from a
tiny localization probability. The physical preparation error enters the
probability of reaching the bin, through the frozen transport bounds.

Let E0 be the actual local physical ground energy. Frozen full local
bounds give `h>=l0=-11.71` on the unrestricted CAR/scalar carrier and
`h>=E0` on its complete gauge-singlet matter subspace. Here local matter
includes the scalar. Let `delta=E0-l0`, using its rational upper bound
where required. The complete source representation conjugate to the
one-particle reference has the lower bound `h>=E0+delta_r`, with

`delta_ec >= .10058529424766495`,
`delta_L >= .09688849789228088`.

These are bounds on every compensated local state, including all raw
lepton charges and spins, from
[the complete charged-sector comparison](../interacting_charged_response/PROOF.md).
They are not merely the minima of the visible elementary-removal support.
The code verifies `delta+delta_r<1`. Every source state outside the
maximal quark module is covered by the inherited local bounds.

## A complete single-charge location measurement

In the joint physical source/reference sector, define Pi_y by requiring
the full local matter gauge-singlet projector at every x other than y.
Retain the full local CAR/scalar space at y, every link state and the
reference. These projectors commute with H0, the Gauss generators and
parity. They do not require the other cells to be in their local ground.

The sum of the abelian Gauss equations over the periodic source fixes
the total matter hypercharge to minus the nonzero reference hypercharge.
Thus Pi_y forces that charge at y. For distinct y and z, Pi_y Pi_z
would force every matter cell to be a gauge singlet, contradicting the
total charge. Hence the eight Pi_y are pairwise orthogonal on this
physical one-reference sector. The ninth effect
`Pi_cloud=I-sum_y Pi_y` retains configurations with several charged
matter cells. It includes screening configurations and is not assigned
zero probability.

The abelian charge in this argument is the actual native hypercharge.
No electroweak symmetry-breaking or electromagnetic particle dictionary
is assumed. At the origin, Gauss includes the retained reference. The
source-only origin Gauss operator need not vanish.

## Flux connectivity and the linear energy bound

Resolve every link into its blank and nonblank projectors. These
projectors commute with all Gauss generators and H0. Let F be the
nonblank edge set, m its cardinality, and ell the graph distance from
zero to y. The stored N2 graph has parallel oriented edges; distances
have shell multiplicities three, three and one at ell=1,2,3.

Suppose y is outside the F-connected component containing zero. Sum
hypercharge Gauss over that component. Every internal link contributes
opposite endpoint charges which cancel. Every edge leaving the component
is blank. All its matter cells are singlets under Pi_y. Only the nonzero
reference charge is left, a contradiction. Therefore a physical state
in Pi_y has an F path from zero to y, and `m>=ell`.

If `m=ell`, F is a shortest simple path and has no extra component or
cycle. At the origin the outgoing endpoint must be conjugate to the
reference representation. At each interior cell, singlet matter and
the two endpoint Gauss equations pair conjugate irreducible link
representations. Schur's lemma propagates the same flux representation
to y. The matter at y is consequently conjugate to the reference.
The local bounds give `H0>=8E0+ell+delta_r` on this support.

If `m>=ell+1`, every cell other than y is a matter singlet, while
the unrestricted bound applies at y. Consequently

`H0 >= 7E0+l0+m >= 8E0+ell+1-delta`
`   >= 8E0+ell+delta_r`.

The last inequality is checked using exact rational isolating intervals.
It covers every additional flux component, every higher charged state
and every screening possibility consistent with Pi_y. Summing these
orthogonal supports proves the full untruncated compression bound

`Pi_y H0 Pi_y >= (8E0+ell+delta_r) Pi_y`, for y different from zero.

Every native hopping term has two different matter endpoints, at least
one different from y. At such an endpoint a single nontrivial native
fermion operator has zero compression between local matter singlets.
Therefore `Pi_y T Pi_y=0` exactly. A normalized state supported in Pi_y
has the same lower bound on its H_lambda expectation. Since the actual
neutral ground energy is at most `8E0+v`, its total excitation energy
including the reference is at least `mu+ell+delta_r-v`. This is an
energy cost for a constrained, single unscreened charge. Pi_y is not
claimed to be invariant under H_lambda.

## A bound that allows all charged matter

For a nonblank set F with m links, at most 2m cells touch F. Every
other matter cell is a singlet, except the origin if it is isolated.
An isolated origin must itself compensate the reference and gains at
least delta_r above E0. Applying l0 to all incident cells gives

`H0 >= 8E0 + gamma N_nonblank + delta_r I_origin_isolated`,
`gamma=1-2delta > .04117058849`.

The same operator inequality holds after summing all link supports.
For the interacting source, dropping the nonnegative isolated-origin
term gives `H_lambda+mu-E_ground >= mu+gamma N_nonblank-2v`.
This controls the expected number of excited links in a state of bounded
mean energy. It supplies no unrestricted charge-separation inequality.

## Full-H low-bin localization at every clock time

Let P_I be the actual total-H spectral projector onto the frozen charged
energy interval I, with upper endpoint b. On Pi_y, H0+mu has spectrum
at least a_y=`8E0+ell+delta_r+mu`. The receipt uses a rational lower
bound on a_y and verifies `a_y>b` for every nonorigin cell.

Put `X=Pi_y P_I`, viewed as a map from ran P_I to ran Pi_y,
`A=(H0+mu)|ran Pi_y` and `B=H_total|ran P_I`. H0 commutes with Pi_y,
and H_total commutes with P_I, so

`A X-X B = -lambda Pi_y T P_I`.

The ordered spectral separation gives the exact finite-dimensional
Sylvester inverse

`X=integral_0^infinity exp(-s A)(-lambda Pi_y T P_I)exp(s B) ds`.

Its integrand norm is bounded by `v exp(-s(a_y-b))`. Hence

`||Pi_y P_I|| <= v/(a_y-b)`.

For every normalized pure or mixed state supported in P_I, at every
real source time, its Pi_y probability is bounded by the square of this
quantity. No isolation or rank hypothesis is placed on the full charged
energy spectrum. All higher full-H states remain in the energy-bin
complement.

The outward rounded upper bounds, conditional on the energy bin, are:

| Native field | One cell, ell=1 | One cell, ell=2 | One cell, ell=3 | All seven remote cells |
| --- | --- | --- | --- | --- |
| e^c | 2.6384390e-15 | 2.2150808e-16 | 7.5730292e-17 | 8.6555743e-15 |
| L | 4.3996076e-16 | 1.0988276e-16 | 4.8820878e-17 | 1.6983521e-15 |

The exact receipt carries rational bounds. These estimates concern
single unscreened charge. They do not bound the cloud outcome, and
they do not imply that the charge is at the origin with probability
one minus the displayed remote bound. The positive short-time arrivals
of the transport package and this low-bin localization concern distinct
settings. The former allows the high energy associated with a flux
excitation; it does not establish resonant propagation of the low bin.

## Original preparation and complete sequential records

The detector measures physical preparation, P0, the native pump, the
total-energy bin, then charge location after any real source time.
This is a sequential instrument. Its energy and location projectors
are not asserted to commute.

Write q for physical preparation success, a for conditional P0 success,
and b for the joint P0/pump/bin probability conditional on preparation.
The frozen preparation fidelity, graph and spectral-window estimates
give q and a lower bounds and a strictly positive b interval. The
conditional pump split after P0 is exactly one-half. Let w_y and w_cloud
be the complete position probabilities after bin selection. The thirteen
records have probabilities

`1-q; q(1-a); qa/2; q(a/2-b); qb w_y (eight cells); qb w_cloud`.

All are nonnegative and their sum is exactly one. Every preparation
failure, energy complement, source copy and quantum reference is
retained. The localization upper bound applies to the actual bin-selected
state even if its preparation differs from the exact ground target.
It is not reduced by a global trace-distance subtraction.

## A screened remote pair is an exact counterexample

Take the physical atomic product density with each local state P_g/2,
all links blank and the static reference vacuum. The actual local
spin purification may be retained. This is a specified physical state
used for a counterexample; it is not substituted for g_lambda.

Apply the origin native pump B_r, and the adjoint branch of stored
edge nine, which moves matter from cell three `(0,1,1)` to cell seven
`(1,1,1)`. This edge is disjoint from the origin. The branch is an actual
term of the same gauge-invariant T with coefficient magnitude 1/4.

The spin-mixed actual local ground has one-body density I/2 on each
native spin/irrep space. Contracting the actual link vacuum Gram I/d
and the native Pauli matrix gives remote-hop norm squared `d/32`.
The origin pump has probability one-half and acts on a disjoint
factor, so the combined positive norm squared is `d/64` for both
e^c (d=1) and L (d=2). No assumed empty matter vacuum enters this
normalization.

The reference has native charge q at zero, while the matter has charges
`(-q,-q,+q)` at cells `(0,3,7)`. Only edge nine is nonblank. The
origin matter locally compensates the reference; the remote pair has
its own one-link flux. Every Gauss constraint holds. There is no
nonblank path from zero to cell three, even though that charged cell
is distance two away. All Pi_y annihilate this state: it belongs to
the retained cloud outcome.

Only three local matter factors differ from the true local ground.
The full local norm 57 bounds its atomic excitation energy by
`3(57-E0)+1`. Every hopping term flips two local fermion parities,
so its expectation in this fixed local-parity pattern is zero. The
same-H total energy relative to g_lambda is therefore at most
`mu+3(57-E0)+1+v`. The construction and its atomic energy bound do
not depend on the placement of a disjoint remote nearest-neighbor
edge. This refutes the unrestricted inference from a distant charged
cell to a string from the reference. It does not prove thermodynamic
string breaking, a scattering state or small energy for the counterexample.

## Charged spectral interpretation

The native removal field label and positive energy-bin weight identify
an odd source response with an explicit compensator. They do not
identify an asymptotically free charged lepton. A finite-regulator
mobile band requires translation-compatible dressed probes and their
positive joint energy-momentum spectral measures, spatial transport
on that same band and control of the reference/boundary energy.

A continuum Coulomb sector has a separate infrared criterion. Under
the local relativistic framework, Gauss-law charge and the asymptotic
field hypotheses of Buchholz, a nonzero charged state cannot have a
sharp mass eigenvalue. A finite-regulator isolated eigenvalue is not
in conflict with that result. Dressed or inclusive charged thresholds
and transport are appropriate continuum targets when those hypotheses
hold. The precise hypotheses, primary citation and explicit architecture
alternatives are stated in [COULOMB.md](COULOMB.md).

The finite proofs here claim neither a thermodynamic charged phase nor
an electromagnetic field limit. They retain the source current, source
clock and the location and energy of the compensator throughout.
