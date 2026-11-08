# Source translation covariance and a field-generated position observable

The source and state are the ordered cubic spin interaction and the selected
zero-field ground used by `spin_canonical_field`. Its continuum has spatial
coordinates `x=p n`, time `T=J p t`, invariant vacuum `Omega`, and strongly
continuous spatial translations `U(a)`. No second Hilbert space or state is
postulated. The construction below takes a closed subspace of this same
continuum representation. The coordinate has the units of the source envelope
coordinate. A positive constant `hbar` sets an action unit; its numerical value
in laboratory units is not calculated here.

## What the existing field construction supplies

For real smooth tests, first with compact Fourier support away from zero,
write `q(F)=Q(F)Omega` and `p(G)=P(G)Omega`. Use Fourier measure
`dnu(k)=d^3k/(2pi)^3`. The source establishes

    ||q(F)||^2 <= 3 integral |Fhat(k)|^2/|k| dnu(k),
    ||p(G)||^2 <= (1/4) integral |k| |Ghat(k)|^2 dnu(k),
    exp(is Q(F)) exp(it P(G))
      = exp(-i s t m integral F G) exp(it P(G)) exp(is Q(F)),       (1)

where the actual selected order `m` obeys `m>1/1000`. The equal-time `Q`
operators strongly commute. The state has zero transverse one-point functions.
Spatial translations act on this same state by

    U(a)Omega=Omega,
    U(a)q(F)=q(F_a),           F_a(x)=F(x-a).                      (2)

The signs in (2) fix the momentum convention used below. Reversing both the
coordinate translation and the generator convention leaves the conclusions
unchanged. These are statements about the spatial envelope translations of
the declared cubic spin source. The rank-three coordinates inside one native
carrier have no proved identification with this envelope coordinate.

The bounds in (1) place `Omega` in the domains of both fields. Differentiate
the matrix coefficient of the Weyl relation once in each parameter, putting
each derivative on a separate vector. This proves the weak form

    2 Im <q(F),p(G)> = m integral F G.                             (3)

No product-domain assumption `Omega in Dom(QP)` is needed. With inner
product linear in its second argument, the Weyl relation gives exactly

    <(exp(-isQ)-I)Omega,(exp(itP)-I)Omega>
      - <(exp(-itP)-I)Omega,(exp(isQ)-I)Omega>
      = (exp(-ist c)-1)<Omega,exp(itP)exp(isQ)Omega>,
    c=m integral F G.

The terms linear in only one parameter cancel. Divide by `s t`. The two
first-order vector differences converge in norm by the separate domain
properties. The left limit is `-2i Im<q(F),p(G)>`; the right limit is `-ic`.
This proves (3) without differentiating a product of unbounded operators.
Equation (3) and Cauchy--Schwarz imply

    ||q(F)||^2 >= m^2 |integral F G|^2
                     / integral |k| |Ghat(k)|^2 dnu(k).

Choosing `Ghat(k)=Fhat(k)/|k|` is admissible on the initial Fourier-cut test
class. Hence

    m^2 integral |Fhat(k)|^2/|k| dnu(k)
      <= ||q(F)||^2
      <= 3 integral |Fhat(k)|^2/|k| dnu(k).                       (4)

This is a source spectral lower bound, not an assumption of a free or Gaussian
field. It extends by completion to the indicated weighted norm. For real
`F,G`, strong commutativity of the `Q` fields makes `<q(F),q(G)>` real: use
their joint spectral measure and Cauchy--Schwarz to integrate the product.
Thus complex linear extension preserves (4). This step is required before
using complex Fourier amplitudes.

## The source Gram operator supplies the intertwiner

Let `H_Q` be the closure of these complex field vectors in the original GNS
space. For `f(k)=Fhat(k)/sqrt(|k|)` define

    T f = q(F).

The initial test class is dense in `L2(R^3,dnu)`. Equation (4) extends `T`
to a bounded, bounded-below map onto `H_Q`:

    m^2 I <= T* T <= 3 I.                                        (5)

The range is closed, and its closure is `H_Q`, so it is onto. Equation (2)
gives `U(a)T=T M_a`, where `M_a f(k)=exp(-i k.a)f(k)`. The positive Gram
operator `T*T` therefore commutes with every `M_a`. The spectral projections
of the three multiplication coordinates generate the scalar multiplication
algebra. Its commutant on scalar `L2` is that same algebra. Consequently

    T*T = multiplication by r(k),       m^2<=r(k)<=3 a.e.          (6)

Equivalently the actual equal-time source structure factor has density
`s_Q(k)=r(k)/|k|`, bounded above and below by (4). In particular the spatial
translation spectrum on this cyclic subspace is absolutely continuous, of
multiplicity one, and has full Lebesgue support. None of these properties is
inferred from a dimension count.

The polar factor

    V = T (T*T)^(-1/2) : L2(R^3,dnu) -> H_Q                       (7)

is unitary onto this actual source subspace and satisfies `U(a)V=V M_a`.
Its positive inverse square root fixes the whitening convention. It is
uniquely determined by the chosen field-vector map and the source Gram form.
In momentum representation,

    V* q(F) = sqrt(s_Q(k)) Fhat(k).

No numerical value or smoothness of `r(k)` is assumed. Its lack of established
smoothness causes no domain problem: all operator domains below are pulled
back through the unitary `V`. Equation (7) is a representation of the source
subspace, not a replacement of the source state by a chosen wave function.

## Exact covariance and the uncertainty relation on this subspace

Define on `H_Q`

    P_j = V (hbar k_j) V*,
    X_j = V (i partial/partial k_j) V*.                           (8)

The first operator is exactly the restriction of the existing spatial
translation generator, because `U(a)|H_Q=exp(-i a.P/hbar)` by (2) and (7).
The second is the selfadjoint Fourier coordinate, with its pulled-back
Sobolev domain. Its joint spectral measure `E_Q(B)` is a normalized sharp
position measurement on `H_Q`. Directly on momentum wave functions,

    exp(i k.a) i partial_j [exp(-i k.a)f]
       = i partial_j f + a_j f,

so, including the selfadjoint operator domains,

    U(a)* X_j U(a) = X_j + a_j I,
    U(a) E_Q(B) U(a)* = E_Q(B+a).                              (9)

The dense common invariant core `D=V Schwartz(R^3)` is preserved by each
`X_j`, `P_j`, and spatial translation. On it

    [X_i,P_j] = i hbar delta_ij I.                             (10)

For normalized vectors in this core, Cauchy--Schwarz applied to centered
operators proves

    Delta X_i Delta P_i >= hbar/2.                            (11)

The usual weak-commutator extension applies to vectors with the required
finite moments and form domains. Equation (11) is a property of this
source-generated position and the same source spatial generator. It neither
identifies a measured particle species nor gives a sharp energy--momentum
shell. The physical evolution has not been proved to preserve `H_Q`; its
full time-history space can contain additional vectors. A source field
quadrature `Q(F)` is not the coordinate operator `X_j` in (8).

## The vacuum requires a no-position outcome

The subspace `H_Q` is orthogonal to `Omega`, by the zero one-point function,
and reduces all spatial translations. Let `Pi_Q` be its projection. Extend
the detector to the full original GNS space by

    E(B)=inclusion E_Q(B) projection,
    E_dark=I-Pi_Q.                                             (12)

Then `E(R^3)+E_dark=I`, the spatial effects obey (9), and the dark effect is
translation invariant. On the vacuum the dark probability is one. On any
normalized `q(F)` the spatial probabilities sum to one and are nonzero;
the construction is not a zero detector on its specified source sector.
If `X_j` is extended by zero on the dark sector, its covariance is
`U(a)* X_j U(a)=X_j+a_j Pi_Q`. Take the direct sum of `D` with the dark
vectors of compact joint spatial-momentum support. This is a dense common
invariant core, and on it `[X_i,P_j]=i hbar delta_ij Pi_Q`. A constant full-space bound `hbar/2` does
not apply to an arbitrary state containing dark weight.

A normalized position POVM on all of `R^3` and the entire vacuum-containing
space cannot have exact translation covariance. Otherwise
`mu(B)=<Omega,E(B)Omega>` would be a translation-invariant probability
measure on `R^3`. Partition space into disjoint translates of a bounded
half-open cube. Each cube has the same mass; a positive value contradicts
finite total mass, and zero mass contradicts normalization. The same
obstruction holds for a normal invariant mixed state. Thus the dark outcome
or a sector restriction is necessary, even in an exact continuum theory.

## What makes a physical detector identification stronger

Equations (7)--(12) construct a source-defined covariant position observable
in the ambient bounded-operator algebra. They do not prove it belongs to a
specified local observer measurement algebra. Whitening is generally
nonlocal, and membership in the GNS bounded-operator algebra alone is not
operational availability. The finite cube proof supplies an exact complete
coordinate detector built from the same type of Gram operation, together
with a positive test showing that its cyclic code leaks under the actual
Hamiltonian. Those finite reflection coordinates are not continuous
translation momenta.

For a family of admissible complete source instruments `E_n`, a sufficient
operational attachment is convergence of all bounded-position statistics
on a dense set of common source histories, with the same spatial action:

    <psi,E_n(f) phi> -> <psi,E(f) phi>

for every bounded continuous spatial test `f`, with the no-position branch
included. The action used to translate those instruments must converge to
the `U(a)` of (2), rather than an independently selected representation.
Identifying the reported coordinate with the carrier rank-three readout
requires this convergence for that actual carrier record map. The existing
finite carrier Gram, and the existence of the unitary `V`, do not supply it.
To infer moment variances or an uncertainty inequality from detector data,
second-moment tail bounds are also needed; bounded-statistic convergence
alone permits a small probability at an arbitrarily distant location.
A dynamically stable particle interpretation further needs a sector
preserved, or quantitatively approximately preserved, by the same time
generator and a demonstrated particle dispersion.
